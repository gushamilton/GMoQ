#!/usr/bin/env python3
"""
FungWas Stage 1 CLI

Production-ready RIF quantile GWAS with jackknife covariance estimation.

Usage:
    fungwas-stage1 --bgen data.bgen --pheno pheno.txt --out results/chr22
    
Outputs:
    {out}.stage1.tsv.gz  - Beta/SE per SNP per tau, plus d_resid (see below)
    {out}.cov.gz         - Upper triangle covariance matrices (binary)
    {out}.cov.ids.tsv.gz - SNP IDs for covariance alignment (one per row)
    {out}.cov.meta.json  - Covariance layout and jackknife degrees-of-freedom model
    {out}.Rtau.tsv       - Population tau correlation matrix (for plugin_cor)

The final column of the TSV, d_resid, is sum_i g_resid_i^2 for that SNP: the
denominator of every per-tau slope. Stage 2 uses it to scale pooled covariance
and to calculate approximate effective degrees of freedom for jackknife
inference. It is residual genotype variation, not a literal sample count.

WARNING about the per-tau standard errors in this file.  se_tau_* are jackknife
estimates, so a normal reference can miss covariance-estimation uncertainty.
The approximate effective degrees of freedom at an extreme tau can be far
below the block count.
Measured at N = 50,000 with --blocks 1000: about 835 df at tau = 0.50, but 112
at tau = 0.02, and 38 at tau = 0.02 for a MAF 0.01 variant. Do not studentise a
per-tau slope against a normal, against blocks - 1, or against the degrees of
freedom Stage 2 reports for a fitted parameter. See docs/KNOWN_ISSUES.md 2.2.
"""

import argparse
import functools
import gzip
import logging
import os
import re
import sqlite3
import struct
import sys
import time
import zlib
from collections import Counter
from pathlib import Path
from typing import Optional, List, Set

import numpy as np

from .core import (
    compute_rif,
    residualize,
    compute_block_scores,
    compute_block_scores_multi,
    compute_block_finalize,
    compute_block_scores_panel,
    compute_block_finalize_panel,
    compute_pooled_panel,
    compute_covariance_from_blocks,
    residual_sum_squares,
    resolve_covar_cols,
    DEFAULT_N_BLOCKS,
    resolve_n_blocks,
    count_nonempty_blocks,
    check_blocks_vs_samples,
    set_numpy_fallback_allowed,
    require_cpp_extension,
    DENSITY_METHODS,
    DENSITY_METHOD_SILVERMAN_GAUSSIAN,
)
from .cov_meta import build_cov_meta, write_cov_meta
from .pooled_cov import (
    SIGMA0_ALL_CHROM,
    build_sigma0_meta,
    compute_sigma0,
    sigma0_analysis_block,
    write_sigma0,
)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(message)s",
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)

#: ``--loco-mode`` values.
#:
#: ``subtract`` (the default, and the only historical behaviour) builds the RIF
#: on ``Y - c*L(chr)``.  ``covariate`` leaves ``Y`` unadjusted and appends
#: ``L(chr)`` as an extra column of the basis the RIF columns -- and the
#: genotype -- are residualised on.  See ``compute_chr_rif_residuals`` and
#: ``docs/LOCO_AS_COVARIATE.md``.
LOCO_MODE_SUBTRACT = "subtract"
LOCO_MODE_COVARIATE = "covariate"
LOCO_MODES = (LOCO_MODE_SUBTRACT, LOCO_MODE_COVARIATE)

#: ``--cov-mode`` values.
#:
#: ``full`` (the default, and the only historical behaviour) writes the per-SNP
#: jackknife covariance payload: ``<prefix>.cov.gz``, ``<prefix>.cov.ids.tsv.gz``
#: and, for int8 storage, ``<prefix>.cov.scale.gz``.  Stage 2 reads that payload
#: only under ``cov_source = "jackknife"``.
#:
#: ``pooled`` writes everything else unchanged -- the same ``.stage1.tsv.gz``,
#: byte for byte, and the same ``Sigma_0`` sidecar -- but omits the per-SNP
#: payload, which Stage 2's production path (``cov_source = "pooled"``) never
#: opens.  The jackknife itself still runs: ``se_tau_*`` in the TSV are the same
#: jackknife standard errors, computed the same way.  Use it when no downstream
#: analysis will ask for ``cov_source = "jackknife"``; the payload cannot be
#: reconstructed from what remains without re-running Stage 1.
#:
#: ``pooled-only`` goes further and skips the jackknife *computation*.  The
#: block accumulator exists only so that leave-one-block-out fits can be
#: recovered by subtraction; the point estimate ``beta_t = N_t / D``, the scale
#: ``D = sum_i g_resid_i^2`` and the pooled ``Sigma_0`` need no block structure
#: at all.  Under this mode ``se_tau_*`` are the pooled standard errors
#: ``sqrt(Sigma_0[t,t] / d_resid)`` -- the diagonal of exactly the covariance
#: Stage 2 applies under ``cov_source = "pooled"`` -- and not jackknife
#: estimates.  That is a different number in the ``se_tau_*`` columns, which is
#: why it is a separate mode rather than a change to ``pooled``: runs already
#: accepted under ``--cov-mode pooled`` have jackknife ``se_tau_*`` and must
#: keep them.
#:
#: MEASURED, chr22 testosterone with LOCO (17,550 variants, N = 399,710, T = 45,
#: ``--blocks 1000 --batch-size 100``) on one BluePebble Xeon Gold 5120 node,
#: threads 1/2/4 bound one per physical core and 8 using the SMT siblings:
#:
#:     threads   full     pooled    pooled-only    SNPs/s (full -> pooled-only)
#:     1         451.7 s  446.3 s   345.5 s        38.9 -> 50.8
#:     2         337.1 s  332.7 s   254.2 s        52.1 -> 69.0
#:     4         234.4 s  229.6 s   172.4 s        74.9 -> 101.8
#:     8         256.7 s  240.3 s   182.5 s        68.4 -> 96.2
#:
#: The kernel itself goes 190.4 s -> 143.5 s at four threads and 190.1 s ->
#: 128.0 s at eight; the write side goes 22.7 s -> 7.6 s, which is the jackknife
#: finalisation leaving.  Note that eight threads is SLOWER than four in every
#: mode: the allocation's "8 cpus" are four physical cores plus their SMT
#: siblings, and filling all eight starves the decode and write threads the
#: ``--pipeline`` overlap depends on.  Four threads is the configuration to
#: quote.
#:
#: Output, unchanged from ``pooled``: 72.8 MB -> 6.3 MB.
#:
#: WHERE THE TIME GOES NOW.  ``--no-pipeline`` serialises decode, kernel and
#: write, so the batch loop splits cleanly.  Four threads, one per physical
#: core, same chromosome and node (17,550 variants):
#:
#:                        full        pooled-only
#:     setup              19.0 s      19.0 s        phenotype, covariates, LOCO,
#:                                                  RIF, BGEN index
#:     BGEN decode        74.3 s      71.7 s        unchanged, and the floor
#:     C++ kernel        117.4 s      82.7 s        -30%
#:     TSV write          21.0 s       6.9 s        -67%, the jackknife
#:                                                  finalisation leaving
#:     Sigma_0             0.05 s      0.05 s
#:     ------------------------------------------
#:     wall              231.7 s     180.3 s
#:
#: With the default ``--pipeline`` the decode hides behind the kernel and the
#: same two runs are 238.4 s and 158.0 s (73.6 and 111.1 SNPs/s).  The kernel is
#: still the critical path; if it went to zero the floor would be setup plus
#: decode, about 91 s.
#:
#: Rebuilding the extension with ``FUNGWAS_STAGE1_ARCH="-march=haswell
#: -mtune=skylake-avx512"`` takes pooled-only to 142.2 s (123.4 SNPs/s) and full
#: to 191.4 s, with decompressed output byte-identical to the SSE2 build in both
#: modes on this chromosome.  ``setup.py`` records that raising the ISA bought
#: nothing; that measurement predates the panel kernels and no longer holds for
#: them.  The default is deliberately left alone here.
COV_MODE_FULL = "full"
COV_MODE_POOLED = "pooled"
COV_MODE_POOLED_ONLY = "pooled-only"
COV_MODES = (COV_MODE_FULL, COV_MODE_POOLED, COV_MODE_POOLED_ONLY)

#: Cov modes that do not write the per-SNP jackknife payload.
COV_MODES_NO_PAYLOAD = (COV_MODE_POOLED, COV_MODE_POOLED_ONLY)


def get_memory_mb() -> float:
    """Get current memory usage in MB."""
    try:
        import resource
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    except:
        return 0.0


def log_mem(msg: str):
    """Log message with memory usage."""
    logger.info(f"[MEM {get_memory_mb():.0f}MB] {msg}")


def read_table_whitespace(file_path: str) -> list[dict]:
    """Read whitespace-delimited tables (space or tab) into list of dicts."""
    rows = []
    with open(file_path) as f:
        header = None
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if header is None:
                header = parts
                continue
            if len(parts) < len(header):
                continue
            rows.append(dict(zip(header, parts)))
    return rows


def parse_list_arg(values: Optional[list[str]]) -> Optional[list[str]]:
    if values is None:
        return None
    parts = []
    for v in values:
        parts.extend([p.strip() for p in v.split(',') if p.strip()])
    return parts or None


def normalize_chr_label(chrom) -> str:
    """Normalize chromosome labels across BGEN/REGENIE conventions.

    Numeric autosomes are canonicalized without leading zeroes so labels like
    ``04``, ``4``, and ``chr04`` all resolve to ``4``. Non-numeric labels are
    returned uppercase after removing any leading ``chr`` prefix.
    """
    chrom_str = str(chrom).strip()
    chrom_str = re.sub(r'^(chr|CHR)', '', chrom_str)
    if chrom_str.isdigit():
        return str(int(chrom_str))
    return chrom_str.upper()


def read_keep_ids(keep_file: str) -> set[str]:
    """Load a sample keep list using the first whitespace-delimited column."""
    keep_set = set()
    with open(keep_file) as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                keep_set.add(parts[0])
    return keep_set


def build_sample_aliases(fid: str, iid: Optional[str]) -> list[str]:
    """Build plausible sample identifiers seen across phenotype/BGEN/REGENIE files."""
    aliases = []
    if fid:
        aliases.append(fid)
    if iid and iid != fid:
        aliases.append(iid)
    if fid and iid:
        aliases.extend([
            f"{fid}_{iid}",
            f"{fid} {iid}",
            f"{fid}:{iid}",
        ])
    return aliases


def resolve_sample_column(sample_lookup: dict[str, int], aliases: list[str]) -> Optional[int]:
    """Resolve the first matching sample alias in a LOCO header lookup."""
    for alias in aliases:
        if alias in sample_lookup:
            return sample_lookup[alias]
    return None


def sanitize_pheno_name(name: str) -> str:
    """Canonicalize a phenotype name the way this project names files and directories.

    Non-alphanumeric characters other than ``-`` and ``_`` become ``_``.  This is the
    same rule used for multi-phenotype output prefixes, and it is the rule REGENIE
    runners here use for per-phenotype directories, so ``testosterone_shrink_0.2``
    canonicalizes to ``testosterone_shrink_0_2``.
    """
    return ''.join(c if c.isalnum() or c in ('-', '_') else '_' for c in name)


def resolve_loco_files(loco_path: str, pheno_cols: list[str]) -> dict[str, str]:
    """
    Resolve REGENIE LOCO inputs to a phenotype->.loco path mapping.

    Supports either a direct `.loco` file or a REGENIE `_pred.list` master file.

    Resolution is **exact and deterministic**.  A REGENIE ``*_pred.list`` line is
    ``<phenotype> <path/to/*.loco>``; the phenotype is taken from the *first* field
    and compared to the requested phenotype name by equality:

    1. equality on the raw pred.list label;
    2. failing that, equality on the canonicalized label
       (:func:`sanitize_pheno_name` on both sides), which absorbs the one naming
       divergence used in this project -- ``.`` written as ``_`` in file and
       directory components, e.g. label ``testosterone_shrink_0_2`` for phenotype
       ``testosterone_shrink_0.2``.

    Substring matching is deliberately *not* used.  The manuscript phenotype names
    are proper substrings of one another (``bmi`` of ``log_bmi``/``rint_bmi``;
    ``testosterone`` of ``testosterone_shrink_0.2``/``_0.5``/``_0.7``; ``crp`` of
    ``log_crp``), so substring matching cannot distinguish them and silently
    attaches the wrong LOCO predictor.

    Anything that does not resolve to exactly one ``.loco`` path raises
    :class:`ValueError` naming the candidates -- it is never resolved by guessing.

    The result depends only on ``pheno_cols`` and the pred.list contents, never on
    ``PYTHONHASHSEED``, dict/set iteration order, or the order of pred.list entries.
    """
    path = Path(loco_path)
    if path.suffix == '.loco':
        if len(pheno_cols) != 1:
            raise ValueError(
                "--loco-preds points to a single .loco file but multiple phenotypes were requested. "
                "Provide the REGENIE *_pred.list file instead."
            )
        return {pheno_cols[0]: str(path)}

    entries: list[list[str]] = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                entries.append(parts)

    loco_files = [parts[-1] for parts in entries if parts[-1].endswith('.loco')]
    if not loco_files:
        raise ValueError(f"No .loco files found in {loco_path}")

    if len(pheno_cols) == 1 and len(loco_files) == 1:
        # Unambiguous by construction: one phenotype, one predictor file.  This also
        # preserves support for REGENIE runs whose internal phenotype label differs
        # from the requested --pheno-col (the weight_height / weight_height2 models
        # run REGENIE on a column literally named "weight").
        return {pheno_cols[0]: loco_files[0]}

    # label -> ordered, de-duplicated candidate paths.  Lists, not sets: order here
    # is pred.list order, and it is only ever used to report candidates.
    by_label: dict[str, list[str]] = {}
    by_canonical: dict[str, list[str]] = {}
    labelled = 0
    for parts in entries:
        loco_file = parts[-1]
        if not loco_file.endswith('.loco') or len(parts) < 2:
            continue
        labelled += 1
        label = parts[0]
        for key, table in ((label, by_label), (sanitize_pheno_name(label), by_canonical)):
            bucket = table.setdefault(key, [])
            if loco_file not in bucket:
                bucket.append(loco_file)

    known_labels = sorted(by_label)
    mapping: dict[str, str] = {}
    unresolved: list[str] = []
    for pheno in pheno_cols:
        matched_on = 'phenotype label'
        candidates = by_label.get(pheno)
        if not candidates:
            candidates = by_canonical.get(sanitize_pheno_name(pheno))
            matched_on = 'canonicalized phenotype label'
        if not candidates:
            unresolved.append(pheno)
            continue
        if len(candidates) > 1:
            raise ValueError(
                f"Ambiguous LOCO resolution in {loco_path}: phenotype '{pheno}' matches "
                f"{len(candidates)} distinct .loco files on {matched_on} "
                f"({', '.join(sorted(candidates))}). Refusing to guess. Provide a "
                "pred.list with exactly one entry per phenotype."
            )
        mapping[pheno] = candidates[0]

    if unresolved:
        detail = (
            f"phenotype labels present in the pred.list: {known_labels}"
            if labelled
            else f"no pred.list line carried a phenotype label before the .loco path "
                 f"({len(loco_files)} unlabelled .loco entries)"
        )
        raise ValueError(
            f"Could not map LOCO files in {loco_path} to phenotypes {unresolved}; "
            f"{detail}. resolve_loco_files matches the pred.list's first field exactly "
            "(optionally after canonicalizing '.' and other non-alphanumerics to '_'); "
            "substring matching is deliberately not used because these phenotype names "
            "are prefixes of one another. Rebuild the pred.list so each requested "
            "phenotype appears as its own first field, or pass a single .loco file for "
            "a single phenotype."
        )

    claimed: dict[str, str] = {}
    for pheno in pheno_cols:
        loco_file = mapping[pheno]
        if loco_file in claimed:
            raise ValueError(
                f"Ambiguous LOCO resolution in {loco_path}: .loco file {loco_file} "
                f"resolved for both '{claimed[loco_file]}' and '{pheno}'. Each phenotype "
                "needs its own LOCO predictor file."
            )
        claimed[loco_file] = pheno

    return mapping


def load_loco_predictions_for_file(loco_file: str,
                                   sample_aliases: list[list[str]]) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """
    Load per-chromosome REGENIE LOCO predictions aligned to the phenotype sample order.
    """
    with open(loco_file) as f:
        header = f.readline().strip().split()
        if not header or len(header) < 2:
            raise ValueError(f"Malformed LOCO header in {loco_file}")
        sample_lookup = {sample_id: idx for idx, sample_id in enumerate(header[1:])}

        matched_cols = []
        for idx, aliases in enumerate(sample_aliases):
            col_idx = resolve_sample_column(sample_lookup, aliases)
            matched_cols.append(col_idx)
        valid_mask = np.array([idx is not None for idx in matched_cols], dtype=bool)
        missing = np.where(~valid_mask)[0]
        matched_cols_arr = np.array(
            [idx for idx in matched_cols if idx is not None],
            dtype=np.int64
        )

        if matched_cols_arr.size == 0:
            raise ValueError(f"No phenotype samples were found in LOCO file {loco_file}")
        if missing.size > 0:
            preview = ', '.join(sample_aliases[i][0] if sample_aliases[i] else f"sample_{i}" for i in missing[:5])
            logger.warning(
                "%d phenotype samples are absent from LOCO file %s and will be dropped. Examples: %s",
                int(missing.size), loco_file, preview
            )

        chr_to_pred = {}
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            chrom = normalize_chr_label(parts[0])
            values = np.asarray(parts[1:], dtype=np.float64)
            if values.shape[0] != len(header) - 1:
                raise ValueError(
                    f"Chromosome {parts[0]} in {loco_file} has {values.shape[0]} predictions, "
                    f"expected {len(header) - 1}"
                )
            chr_to_pred[chrom] = values[matched_cols_arr]

    if not chr_to_pred:
        raise ValueError(f"No chromosome predictions found in {loco_file}")
    return chr_to_pred, valid_mask


def resolve_loco_scale(mode: str,
                       y: np.ndarray,
                       per_chr: dict[str, np.ndarray],
                       pheno: str) -> float:
    """
    Return the factor that puts a REGENIE `.loco` prediction on the phenotype's scale.

    REGENIE standardises quantitative phenotypes in step 1, so the predictions it writes
    to `.loco` are on the standardised scale: their sd is ~0.19-0.36 whatever the
    phenotype's own sd is. Subtracting them from an unstandardised Y (mode 'none')
    under-corrects traits with sd > 1 and over-corrects traits with sd < 1.

    The exact factor is REGENIE's `scale_Y`, which is the sd of the phenotype *after*
    regressing out REGENIE's step-1 covariates (Pheno.cpp::residualize_phenotypes), not
    the sd of the raw phenotype. Stage 1 cannot compute that directly: REGENIE's step-1
    covariate set is not an input here and need not match Stage 1's own --covar.

    Modes:
      'regression' - OLS slope of Y on the prediction, pooled over chromosomes. The
                     default. Recovers REGENIE's analysed sd to within 0.02% on every
                     trait that could be ground-truthed, because it needs no knowledge of
                     the step-1 covariate set and absorbs REGENIE's ridge shrinkage too.
      'sd'         - multiply by sd(Y). Correct only when the step-1 covariates explain
                     little of Y; it over-subtracts by 19% for `weight`, whose step-1
                     covariates (height terms) carry R^2 = 0.30.
      'none'       - subtract as-is (pre-v0.4.0 behaviour; kept to reproduce old runs).
      a float      - use that factor verbatim.
    """
    finite = np.isfinite(y)
    y_sd = float(np.std(y[finite])) if finite.any() else float('nan')

    if mode == 'none':
        scale = 1.0
    elif mode == 'sd':
        scale = y_sd
    elif mode == 'regression':
        num = 0.0
        den = 0.0
        for pred in per_chr.values():
            ok = finite & np.isfinite(pred)
            if not ok.any():
                continue
            yc = y[ok] - y[ok].mean()
            pc = pred[ok] - pred[ok].mean()
            num += float(pc @ yc)
            den += float(pc @ pc)
        if den <= 0:
            raise ValueError(
                f"--loco-scale regression failed for '{pheno}': the LOCO predictions have "
                f"zero variance. Use --loco-scale sd."
            )
        scale = num / den
    else:
        try:
            scale = float(mode)
        except ValueError:
            raise ValueError(f"--loco-scale must be 'sd', 'regression', 'none' or a number; got {mode!r}")

    pred_sds = [float(np.std(v[np.isfinite(v)])) for v in per_chr.values() if np.isfinite(v).any()]
    pred_sd = float(np.median(pred_sds)) if pred_sds else float('nan')
    logger.info(
        "LOCO scale for '%s': mode=%s factor=%.6g (sd(Y)=%.6g, median sd(LOCO)=%.6g, "
        "effective sd of subtracted term=%.6g)",
        pheno, mode, scale, y_sd, pred_sd, scale * pred_sd
    )

    # Guard on the thing that actually matters: whether the adjustment removes variance
    # or adds it. Measured directly rather than inferred from an sd ratio, so a trait with
    # strongly predictive covariates cannot trip it spuriously.
    var_before, var_after, n_used = 0.0, 0.0, 0
    for pred in per_chr.values():
        ok = finite & np.isfinite(pred)
        if not ok.any():
            continue
        var_before += float(np.var(y[ok]))
        var_after += float(np.var(y[ok] - scale * pred))
        n_used += 1
    if n_used:
        var_before /= n_used
        var_after /= n_used
        logger.info(
            "LOCO adjustment for '%s': var(Y) %.6g -> %.6g (%+.1f%%)",
            pheno, var_before, var_after,
            100.0 * (var_after - var_before) / var_before if var_before > 0 else float('nan')
        )
        if var_after > var_before:
            logger.warning(
                "LOCO adjustment for '%s' INCREASES the analysed variance (%.6g -> %.6g). "
                "Subtracting this term injects noise instead of removing it, which biases "
                "every non-location quantile effect by sd_correct/sd_analysed. This is the "
                "signature of an unscaled REGENIE prediction (--loco-scale none); use the "
                "default --loco-scale regression.",
                pheno, var_before, var_after
            )
    return scale


def load_loco_predictions(loco_path: str,
                          pheno_cols: list[str],
                          sample_aliases: list[list[str]]) -> tuple[dict[str, dict[str, np.ndarray]], np.ndarray]:
    """Load aligned LOCO predictions for one or more phenotypes and return the shared valid sample mask."""
    loco_files = resolve_loco_files(loco_path, pheno_cols)
    raw_predictions = {}
    valid_masks = {}
    combined_valid_mask = np.ones(len(sample_aliases), dtype=bool)
    for pheno in pheno_cols:
        loco_file = loco_files[pheno]
        log_mem(f"Loading LOCO predictions for {pheno} from {loco_file}")
        pheno_preds, valid_mask = load_loco_predictions_for_file(loco_file, sample_aliases)
        raw_predictions[pheno] = pheno_preds
        valid_masks[pheno] = valid_mask
        combined_valid_mask &= valid_mask

    if not np.any(combined_valid_mask):
        raise ValueError("No samples remain after intersecting phenotype rows with LOCO predictions.")

    predictions = {}
    for pheno in pheno_cols:
        pheno_valid_idx = np.where(valid_masks[pheno])[0]
        keep_positions = np.nonzero(combined_valid_mask[pheno_valid_idx])[0]
        predictions[pheno] = {
            chrom: values[keep_positions]
            for chrom, values in raw_predictions[pheno].items()
        }

    return predictions, combined_valid_mask


def load_snps_list(snp_file: str) -> List[str]:
    """Load list of SNP rsIDs from file (one per line)."""
    snps = []
    with open(snp_file) as f:
        for line in f:
            rsid = line.strip().split()[0] if line.strip() else ''
            if rsid:
                snps.append(rsid)
    return snps


def find_multiallelic_rsids_from_bgi(bgen_path: str, requested_rsids: Optional[Set[str]] = None) -> Set[str]:
    """Return rsIDs with more than one BGEN record in the .bgi index."""
    bgi_path = f"{bgen_path}.bgi"
    if not os.path.exists(bgi_path):
        raise FileNotFoundError(f"Missing BGEN index required for multiallelic filtering: {bgi_path}")

    conn = sqlite3.connect(bgi_path)
    try:
        cur = conn.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='Variant'")
        if cur.fetchone() is None:
            raise ValueError(f"BGEN index does not contain Variant table: {bgi_path}")

        if requested_rsids is None:
            query = (
                "SELECT rsid FROM Variant "
                "WHERE rsid IS NOT NULL AND rsid != '' "
                "GROUP BY rsid HAVING COUNT(*) > 1"
            )
            return {row[0] for row in cur.execute(query)}

        if not requested_rsids:
            return set()

        cur.execute("CREATE TEMP TABLE requested_rsid (rsid TEXT PRIMARY KEY)")
        cur.executemany(
            "INSERT OR IGNORE INTO requested_rsid (rsid) VALUES (?)",
            ((rsid,) for rsid in requested_rsids),
        )
        query = (
            "SELECT v.rsid FROM Variant v "
            "JOIN requested_rsid r ON v.rsid = r.rsid "
            "GROUP BY v.rsid HAVING COUNT(*) > 1"
        )
        return {row[0] for row in cur.execute(query)}
    finally:
        conn.close()


def query_bgi_records_for_rsids(bgen_path: str, requested_rsids: Set[str]) -> List[tuple]:
    """Return the .bgi records for a set of rsIDs, in BGEN file order, in one index scan.

    Returns tuples ``(file_start_position, chromosome, position, rsid)`` for every BGEN
    record whose rsID is in ``requested_rsids``. Records are ordered by
    ``file_start_position``, i.e. exactly the order a record-by-record walk of the BGEN
    visits them in, so downstream ordering is unchanged.

    The UK Biobank .bgi has ``PRIMARY KEY (chromosome, position, rsid, allele1, allele2,
    file_start_position)`` and no index on ``rsid`` alone, so ``WHERE rsid = ?`` degrades
    to a full table scan *per lookup*. A single JOIN against a TEMP table of the requested
    rsIDs pays one scan for the whole set instead - the same pattern already used by
    find_multiallelic_rsids_from_bgi.
    """
    bgi_path = f"{bgen_path}.bgi"
    if not os.path.exists(bgi_path):
        raise FileNotFoundError(f"Missing BGEN index required for indexed variant access: {bgi_path}")

    if not requested_rsids:
        return []

    conn = sqlite3.connect(bgi_path)
    try:
        cur = conn.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='Variant'")
        if cur.fetchone() is None:
            raise ValueError(f"BGEN index does not contain Variant table: {bgi_path}")

        cur.execute("CREATE TEMP TABLE requested_rsid (rsid TEXT PRIMARY KEY)")
        cur.executemany(
            "INSERT OR IGNORE INTO requested_rsid (rsid) VALUES (?)",
            ((rsid,) for rsid in requested_rsids),
        )
        query = (
            "SELECT v.file_start_position, v.chromosome, v.position, v.rsid "
            "FROM Variant v "
            "JOIN requested_rsid r ON v.rsid = r.rsid "
            "ORDER BY v.file_start_position"
        )
        return cur.execute(query).fetchall()
    finally:
        conn.close()


def load_analysis_inputs(pheno_file: str, pheno_cols: list[str], covar_file: Optional[str],
                         keep_file: Optional[str] = None,
                         covar_cols: Optional[list[str]] = None) -> tuple:
    """Load phenotypes/covariates and return aligned raw inputs for analysis.

    Returns ``(sample_ids, sample_aliases, y_values, X, covar_cols)`` where
    ``covar_cols`` are the covariate columns resolved against the file, in the
    order they occupy ``X`` after the intercept.
    """
    keep_set = None
    if keep_file:
        log_mem(f"Loading sample keep list from {keep_file}")
        keep_set = read_keep_ids(keep_file)
        log_mem(f"Keeping {len(keep_set)} samples")

    log_mem(f"Loading phenotypes from {pheno_file}")
    pheno_data = read_table_whitespace(pheno_file)

    for col in pheno_cols:
        if col not in pheno_data[0]:
            available = [c for c in pheno_data[0].keys() if c not in ['FID', 'IID']]
            raise ValueError(f"Phenotype column '{col}' not found. Available: {available}")

    if covar_file:
        log_mem(f"Loading covariates from {covar_file}")
        covar_data = read_table_whitespace(covar_file)
        covar_map = {r['FID']: r for r in covar_data}
        available_covar_cols = [c for c in covar_data[0].keys() if c not in ['FID', 'IID']]
        covar_cols = resolve_covar_cols(covar_cols, available_covar_cols)
    elif covar_cols is not None:
        raise ValueError("covar_cols provided but no covariate file (--covar) was supplied.")
    else:
        covar_map = {}
        covar_cols = []

    sample_ids = []
    sample_aliases = []
    y_values = {col: [] for col in pheno_cols}
    X_list = []

    for row in pheno_data:
        fid = row['FID']
        iid = row.get('IID', fid)
        if keep_set is not None and fid not in keep_set:
            continue

        try:
            pheno_vals = []
            for col in pheno_cols:
                val = float(row[col]) if row[col] != 'NA' else np.nan
                pheno_vals.append(val)
            if np.any(np.isnan(pheno_vals)):
                continue

            if covar_file:
                if fid not in covar_map:
                    continue
                c_row = covar_map[fid]
                x_vals = [float(c_row[c]) if c_row[c] != 'NA' else np.nan for c in covar_cols]
                if np.any(np.isnan(x_vals)):
                    continue
                X_list.append(x_vals)

            sample_ids.append(fid)
            sample_aliases.append(build_sample_aliases(fid, iid))
            for col, val in zip(pheno_cols, pheno_vals):
                y_values[col].append(val)
        except (ValueError, KeyError):
            continue

    n_samples = len(sample_ids)
    log_mem(f"Loaded {n_samples} samples with valid phenotypes")

    if covar_file:
        X = np.column_stack([np.ones(n_samples), np.array(X_list)])
    else:
        X = np.ones((n_samples, 1))

    # covar_cols is returned because the pooled covariance sidecar has to record
    # exactly which covariate columns entered X, in order: Sigma_0 is the
    # covariate-residualised influence covariance, so a reader that does not know
    # the covariate set cannot tell whether it applies to their scan.
    return sample_ids, sample_aliases, y_values, X, list(covar_cols)


def build_rif_residuals(y_values: dict[str, list[float] | np.ndarray], pheno_cols: list[str],
                        X: np.ndarray, taus: np.ndarray,
                        density_method: str = DENSITY_METHOD_SILVERMAN_GAUSSIAN
                        ) -> tuple[list[np.ndarray], np.ndarray, list[np.ndarray]]:
    """Compute residualized RIF matrices for one or more phenotypes."""
    log_mem(f"Residualizing phenotypes on {X.shape[1]} covariates...")
    Q, _ = np.linalg.qr(X)

    rif_resid_list = []
    q_tau_list = []
    for col in pheno_cols:
        log_mem(f"Computing RIF matrix for {col}...")
        t_rif = time.perf_counter()
        y = np.asarray(y_values[col], dtype=np.float64)
        rif, q_tau, _ = compute_rif(y, taus, density_method=density_method)
        log_mem(f"RIF build for {col}: {time.perf_counter() - t_rif:.2f}s")
        rif_resid_list.append(rif - Q @ (Q.T @ rif))
        q_tau_list.append(q_tau)

    return rif_resid_list, Q, q_tau_list


class VariantGenerator:
    """Generator that yields variants one at a time while keeping bfile handle alive."""
    def __init__(self, bfile, snp_set, excluded_rsids: Optional[Set[str]] = None):
        self.bfile = bfile  # Keep handle alive!
        self.snp_set = snp_set
        self.excluded_rsids = excluded_rsids or set()
        self.found = 0
        self.target = len(snp_set)
        
    def __iter__(self):
        for var in self.bfile:
            if var.rsid in self.excluded_rsids:
                continue
            if var.rsid in self.snp_set:
                self.found += 1
                yield var
                if self.found >= self.target:
                    break


class SeekingVariantGenerator:
    """Yield the requested variants in BGEN file order by seeking to indexed offsets.

    Equivalent to VariantGenerator, but reads only the selected records instead of
    walking every record in the file. ``records`` must be
    ``(file_start_position, chromosome, position, rsid)`` tuples for the target rsIDs,
    already ordered by ``file_start_position``, so the yielded order matches the
    record-by-record walk exactly.

    Seeking goes through ``BgenReader.fetch(chrom, pos, pos)``, whose .bgi query is an
    index SEARCH on the (chromosome, position) primary-key prefix. (``with_rsid`` is
    deliberately not used: the .bgi has no rsid index, so it is a full table scan per
    lookup.) Records sharing a position that are not targets are discarded by file
    offset, and targets sharing a position are re-sorted by offset so file order is
    preserved there too.
    """

    def __init__(self, bfile, records, snp_set):
        self.bfile = bfile  # Keep handle alive!
        self.records = records
        self.snp_set = snp_set
        self.found = 0
        self.target = len(snp_set)

    def __iter__(self):
        recs = self.records
        n = len(recs)
        i = 0
        while i < n:
            chrom = recs[i][1]
            pos = recs[i][2]
            wanted = set()
            j = i
            while j < n and recs[j][1] == chrom and recs[j][2] == pos:
                wanted.add(recs[j][0])
                j += 1
            i = j

            # Materialise the fetch before yielding: BgenReader.fetch is a generator over
            # a sqlite cursor shared with the rest of the index, so it must be drained
            # before control leaves this frame.
            candidates = [var for var in list(self.bfile.fetch(chrom, pos, pos))
                          if var.fileoffset in wanted]
            candidates.sort(key=lambda var: var.fileoffset)

            for var in candidates:
                if var.rsid not in self.snp_set:
                    continue
                self.found += 1
                yield var
                if self.found >= self.target:
                    return


class FilteredBgenIterator:
    """Iterator over all BGEN records except excluded rsIDs."""
    def __init__(self, bfile, excluded_rsids: Optional[Set[str]] = None):
        self.bfile = bfile
        self.excluded_rsids = excluded_rsids or set()

    def __iter__(self):
        for var in self.bfile:
            if var.rsid not in self.excluded_rsids:
                yield var


def get_variants_bgen(bgen_path: str, sample_path: Optional[str],
                      snp_list: Optional[List[str]], allow_multiallelic_rsids: bool = False,
                      seek_variants: bool = True):
    """
    Get variants from BGEN file using the bgen package.
    
    Uses delay_parsing=True for efficient streaming and pre-filters SNPs
    using the BGEN index to avoid scanning variants not in the file.
    
    Returns
    -------
    variants : iterable of variant objects
    samples : list of sample IDs
    bfile_handle : BgenReader handle (keep alive during iteration)
    """
    from bgen import BgenReader
    
    log_mem("Using bgen package with delay_parsing...")
    bfile = BgenReader(bgen_path, sample_path=sample_path or '', delay_parsing=True)
    
    if snp_list:
        log_mem(f"Fetching {len(snp_list)} specific SNPs using bgen...")

        snp_set = set(snp_list)
        bgi_path = f"{bgen_path}.bgi"

        if seek_variants and os.path.exists(bgi_path):
            # One .bgi scan gives membership, the multiallelic-rsID record counts AND the
            # file offsets, so the selected records can be seeked to directly instead of
            # walking every record in the chromosome.
            log_mem("Reading BGEN index for variant list...")
            records = query_bgi_records_for_rsids(bgen_path, snp_set)
            record_counts = Counter(rec[3] for rec in records)
            valid_snps = set(record_counts)
            invalid_count = len(snp_set) - len(valid_snps)

            if invalid_count > 0:
                log_mem(f"{invalid_count} requested SNPs not in BGEN, processing {len(valid_snps)}")

            if len(valid_snps) == 0:
                logger.warning("No valid SNPs found in BGEN file!")

            if not allow_multiallelic_rsids:
                excluded_rsids = {rsid for rsid, count in record_counts.items() if count > 1}
                if excluded_rsids:
                    valid_snps = valid_snps - excluded_rsids
                    log_mem(
                        f"Excluded {len(excluded_rsids)} requested rsIDs with multiple BGEN records; "
                        f"processing {len(valid_snps)} unambiguous SNPs"
                    )

            records = [rec for rec in records if rec[3] in valid_snps]
            log_mem(f"Seeking to {len(records)} indexed BGEN records (skipping the full-file walk)")
            gen = SeekingVariantGenerator(bfile, records, valid_snps)
            return gen, bfile.samples, bfile

        # Fallback: no usable index, or --no-seek-variants. Walk the file record-by-record.
        log_mem("Reading BGEN index for variant list...")
        bgen_rsids = set(bfile.rsids())
        valid_snps = snp_set & bgen_rsids
        invalid_count = len(snp_set) - len(valid_snps)
        
        if invalid_count > 0:
            log_mem(f"{invalid_count} requested SNPs not in BGEN, processing {len(valid_snps)}")
        
        if len(valid_snps) == 0:
            logger.warning("No valid SNPs found in BGEN file!")

        excluded_rsids = set()
        if not allow_multiallelic_rsids:
            excluded_rsids = find_multiallelic_rsids_from_bgi(bgen_path, valid_snps)
            if excluded_rsids:
                valid_snps = valid_snps - excluded_rsids
                log_mem(
                    f"Excluded {len(excluded_rsids)} requested rsIDs with multiple BGEN records; "
                    f"processing {len(valid_snps)} unambiguous SNPs"
                )
        
        # Use generator for memory-efficient streaming
        gen = VariantGenerator(bfile, valid_snps, excluded_rsids)
        return gen, bfile.samples, bfile
    else:
        # Return the whole file for iteration
        excluded_rsids = set()
        if not allow_multiallelic_rsids:
            excluded_rsids = find_multiallelic_rsids_from_bgi(bgen_path, None)
            if excluded_rsids:
                log_mem(f"Excluded {len(excluded_rsids)} rsIDs with multiple BGEN records")
        log_mem("No SNP list provided; streaming full BGEN directly.")
        return FilteredBgenIterator(bfile, excluded_rsids), bfile.samples, bfile


def match_samples(bgen_samples: list, sample_ids: list):
    """Match phenotype samples to BGEN samples and return indices."""
    duplicated_bgen = [sample for sample, count in Counter(bgen_samples).items() if count > 1]
    if duplicated_bgen:
        preview = ", ".join(map(str, duplicated_bgen[:5]))
        raise ValueError(
            f"BGEN sample IDs are not unique ({len(duplicated_bgen)} duplicates; "
            f"examples: {preview}). Refusing ambiguous sample matching."
        )
    duplicated_pheno = [sample for sample, count in Counter(sample_ids).items() if count > 1]
    if duplicated_pheno:
        preview = ", ".join(map(str, duplicated_pheno[:5]))
        raise ValueError(
            f"Phenotype sample IDs are not unique ({len(duplicated_pheno)} duplicates; "
            f"examples: {preview}). Refusing to analyse a sample more than once."
        )
    sample_to_idx = {s: i for i, s in enumerate(bgen_samples)}
    matched_indices = [sample_to_idx.get(fid) for fid in sample_ids]
    valid_mask = np.array([idx is not None for idx in matched_indices])
    bgen_indices = np.array([idx for idx in matched_indices if idx is not None])
    return valid_mask, bgen_indices


class _NeedMoreRecordBytes(Exception):
    """The record window read from the BGEN was too small to reach the bit depth."""


#: Read-only file descriptors for BGENs we take preamble windows from, keyed by
#: path.  One per process; closed when the process exits.
_BGEN_RAW_FDS: dict = {}

#: First window size tried when reading a variant record's preamble directly.
#: The bit-depth byte sits just past the per-sample ploidy array, which in a
#: layout-2 record is a constant byte repeated n_samples times and so deflates to
#: a few hundred bytes; 64 KiB covers production records with room to spare, and
#: the window doubles on demand when it does not.
_BGEN_WINDOW_BYTES = 1 << 16


def _bgen_record_window(bgen_path: str, offset: int, size: int) -> bytes:
    """Read ``size`` bytes of a BGEN starting at ``offset``.

    A short read means end-of-file, which the caller uses to tell "my window was
    too small" apart from "this record really is truncated".
    """
    fd = _BGEN_RAW_FDS.get(bgen_path)
    if fd is None:
        fd = os.open(bgen_path, os.O_RDONLY)
        _BGEN_RAW_FDS[bgen_path] = fd
    return os.pread(fd, size, offset)


def _variant_bit_depth(variant, bgen_path: Optional[str] = None) -> int:
    """Read a BGEN variant's probability bit depth.

    ``bgen`` stores the depth in each layout-2 record but does not currently expose it
    on its Python variant object, so the layout-2 preamble is parsed directly.

    The depth is read for **every** variant, not inferred once and reused: the BGEN
    spec allows it to vary from record to record, and getting it wrong silently
    selects the wrong dosage decode path.

    Where the record's file offset is known, the preamble is read straight from the
    BGEN with ``pread``.  The fallback, ``variant.copy_data()``, returns the whole
    record as a Python ``list`` of per-byte ``int`` objects -- 799,454 of them for a
    production 400k-sample record -- and materialising that list cost 8.9 ms per
    variant, 2.4x the cost of actually decoding the dosages (3.8 ms) and about 70%
    of all per-variant BGEN time.  The bytes parsed are the same either way.
    """
    advertised = getattr(variant, "bit_depth", None)
    if advertised is not None:
        return int(advertised)

    # bgen 1.9.7 wheels do not expose layout/compression despite newer stubs
    # advertising them; the production input is layout 2.
    layout = int(getattr(variant, "layout", 2))
    if layout == 1:
        return 16
    if layout != 2:
        raise ValueError(f"Unsupported BGEN layout: {layout}")

    fileoffset = getattr(variant, "fileoffset", None)
    if bgen_path is not None and fileoffset is not None:
        window = _BGEN_WINDOW_BYTES
        while True:
            raw = _bgen_record_window(bgen_path, int(fileoffset), window)
            at_eof = len(raw) < window
            try:
                return _parse_layout2_bit_depth(raw, truncated=not at_eof,
                                                variant=variant)
            except _NeedMoreRecordBytes:
                if at_eof:
                    raise ValueError(
                        "BGEN variant is truncated while reading its genotype header"
                    ) from None
                window *= 2
    return _parse_layout2_bit_depth(bytes(variant.copy_data()), truncated=False,
                                    variant=variant)


def _parse_layout2_bit_depth(raw: bytes, truncated: bool = False, variant=None) -> int:
    """Parse a layout-2 record preamble for its probability bit depth.

    ``truncated`` says ``raw`` is a prefix of the record rather than the whole of
    it, so running off the end means "read a bigger window", not "corrupt file".
    """
    offset = 0

    def take(width: int) -> bytes:
        nonlocal offset
        end = offset + width
        if end > len(raw):
            if truncated:
                raise _NeedMoreRecordBytes
            raise ValueError("BGEN variant is truncated while reading its genotype header")
        value = raw[offset:end]
        offset = end
        return value

    def take_prefixed(width: int) -> None:
        length = int.from_bytes(take(width), "little")
        take(length)

    # Variant metadata precedes the layout-2 genotype block.
    take_prefixed(2)  # varid
    take_prefixed(2)  # rsid
    take_prefixed(2)  # chromosome
    take(4)  # position
    n_alleles = int.from_bytes(take(2), "little")
    expected_alleles = getattr(variant, "alleles", None)
    if not 1 <= n_alleles <= 1024:
        raise ValueError("BGEN variant has an invalid allele count")
    if expected_alleles is not None and n_alleles != len(expected_alleles):
        raise ValueError("BGEN variant preamble has the wrong allele count")
    for _ in range(n_alleles):
        take_prefixed(4)
    block_length = int.from_bytes(take(4), "little")
    if truncated and offset + block_length > len(raw):
        # Only a prefix of the genotype block is in this window. That is enough
        # for the bounded inflate below, which never needs the whole block; if it
        # turns out not to be, the parse asks for a bigger window.
        block = raw[offset:]
        offset = len(raw)
        block_truncated = True
    else:
        block = take(block_length)
        block_truncated = False

    compression = getattr(variant, "compression", None)
    if compression in ("zlib", 1) or compression is None:
        if len(block) < 4:
            raise ValueError("BGEN compressed genotype block is missing its uncompressed length")
        # Only the preamble is needed.  Limiting inflate avoids decoding every
        # probability just to discover the per-record bit depth.  Older bgen wheels
        # expose neither layout nor compression, so fall back to a raw layout-2 block
        # when the payload is not zlib-compressed.
        inflater = zlib.decompressobj()
        try:
            compressed = block[4:]
            block = inflater.decompress(compressed, 8)
            if len(block) < 8:
                block += inflater.decompress(inflater.unconsumed_tail, 8 - len(block))
            compression = "zlib"
        except zlib.error:
            if block[4:8] == b"\x28\xb5\x2f\xfd":
                raise ValueError(
                    "BGEN genotype block is zstd-compressed but the reader exposes no compression metadata"
                )
            block = block
            compression = None
    elif compression == "zstd":
        try:
            import zstandard as zstd
        except ImportError as exc:
            raise ValueError("Cannot inspect zstd-compressed BGEN bit depth without zstandard") from exc
        if block_truncated:
            # zstd here is decompressed in one shot, so a prefix will not do.
            raise _NeedMoreRecordBytes
        if len(block) < 4:
            raise ValueError("BGEN zstd genotype block is missing its uncompressed length")
        expected_length = int.from_bytes(block[:4], "little")
        block = zstd.ZstdDecompressor().decompress(block[4:], max_output_size=expected_length)
    elif compression not in (None, ""):
        raise ValueError(f"Unsupported BGEN compression: {compression}")

    if len(block) < 8:
        if block_truncated:
            raise _NeedMoreRecordBytes
        raise ValueError("BGEN genotype block is missing its layout-2 preamble")
    n_samples = int.from_bytes(block[:4], "little")
    expected_samples = getattr(variant, "n_samples", None)
    if expected_samples is not None and n_samples != int(expected_samples):
        raise ValueError("BGEN genotype preamble has the wrong sample count")
    min_ploidy = block[6]
    max_ploidy = block[7]
    if min_ploidy > max_ploidy or max_ploidy > 63:
        raise ValueError("BGEN genotype preamble has invalid ploidy bounds")
    # Layout 2 stores one ploidy/missingness byte for every sample even when
    # min_ploidy == max_ploidy; the reader may skip materialising those values,
    # but they remain in the on-disk preamble.
    header_offset = 8 + n_samples
    if compression == "zlib" and len(block) < header_offset + 2:
        needed = header_offset + 2 - len(block)
        block += inflater.decompress(inflater.unconsumed_tail, needed)
    elif compression == "zstd" and len(block) < header_offset + 2:
        raise ValueError("Cannot bounded-read variable-ploidy zstd BGEN preamble")
    if header_offset + 2 > len(block):
        if block_truncated:
            # The compressed prefix in this window ran out before the bit-depth
            # byte; retry against a larger slice of the record.
            raise _NeedMoreRecordBytes
        raise ValueError("BGEN genotype block is missing phased/bit-depth fields")
    phased = block[header_offset]
    bit_depth = int(block[header_offset + 1])
    if phased not in (0, 1) or not 1 <= bit_depth <= 32:
        raise ValueError("BGEN genotype preamble has invalid phased/bit-depth fields")
    return bit_depth


def extract_dosage(
    variant,
    bgen_indices: Optional[np.ndarray] = None,
    probability_bit_depth: Optional[int] = None,
):
    """
    Extract dosage from variant for the first allele (alleles[0]).
    This ensures dosage matches effect_allele which is set to alleles[0].
    
    BGEN dosage is the expected count of allele[1] (alt_dosage).
    Dosage for allele[0] = 2 - alt_dosage (assuming diploid).
    
    Parameters
    ----------
    variant : variant object from bgen
    bgen_indices : indices to extract, or None for all samples
    
    Returns
    -------
    dosages : numpy array of dosages for alleles[0]
    """
    if probability_bit_depth is None:
        probability_bit_depth = _variant_bit_depth(variant)

    if probability_bit_depth <= 8:
        # alt_dosage is the expected count of allele[1].  For diploid, dosage of
        # allele[0] = 2 - alt_dosage.  bgen's native 8-bit path is vectorised.
        dosages = 2.0 - variant.alt_dosage
    else:
        # The bgen 1.9.7 high-bit alt_dosage path can decode the stored probabilities
        # correctly while returning an inconsistent dosage.  Reconstruct allele[0]
        # directly from the decoded [AA, Aa, aa] probabilities for those variants.
        if bool(variant.is_phased):
            raise ValueError("Phased BGEN dosage extraction is unsupported for high bit depths")
        probabilities = np.asarray(variant.probabilities, dtype=np.float64)
        if probabilities.ndim != 2 or probabilities.shape[1] != 3:
            raise ValueError(
                "High-bit BGEN dosage extraction requires diploid biallelic probabilities"
            )
        dosages = 2.0 * probabilities[:, 0] + probabilities[:, 1]
    if bgen_indices is not None:
        return dosages[bgen_indices]
    return dosages


@functools.lru_cache(maxsize=8)
def _triu_indices_cached(n_taus: int):
    """``np.triu_indices(n_taus)``, built once instead of once per SNP.

    The uncached call allocated and filled two length-T(T+1)/2 index arrays for
    every variant (17550 times per chromosome) to produce an index that depends
    only on the tau count.
    """
    return np.triu_indices(n_taus)


def write_stats_batch(stats_list, batch_info, args, f_tsv_list, f_cov_list, f_cov_ids_list,
                      f_cov_scale_list, all_betas_list, n_cov_elements, n_taus):
    """Write one processed batch to the stage1 outputs."""
    triu = _triu_indices_cached(n_taus)
    want_cov = args.cov_mode not in COV_MODES_NO_PAYLOAD
    for p_idx, stats in enumerate(stats_list):
        f_tsv = f_tsv_list[p_idx]
        f_cov = f_cov_list[p_idx]
        f_cov_ids = f_cov_ids_list[p_idx]
        f_cov_scale = f_cov_scale_list[p_idx]
        all_betas = all_betas_list[p_idx]

        for snp_stats, info in zip(stats, batch_info):
            if np.isnan(snp_stats[0, 0]):
                continue

            beta, se, Sigma = compute_covariance_from_blocks(snp_stats, args.blocks)
            # Exactly the denominator the jackknife divides by, summed over blocks.
            d_resid = float(np.sum(snp_stats[:, 0]))

            row = list(info)
            for b, s in zip(beta, se):
                row.extend([f'{b:.6g}', f'{s:.6g}'])
            row.append(f'{d_resid:.9g}')
            f_tsv.write('\t'.join(map(str, row)) + '\n')

            if not want_cov:
                # --cov-mode pooled: Stage 2's production path reads only the
                # pooled Sigma_0 sidecar, so the per-SNP jackknife covariance
                # payload is computed but never written.
                continue

            cov_upper = Sigma[triu]
            if args.cov_dtype == 'int8':
                scale = np.max(np.abs(cov_upper)) / 127.0
                if not np.isfinite(scale) or scale <= 0:
                    scale = 1.0
                cov_q = np.clip(np.rint(cov_upper / scale), -127, 127).astype(np.int8)
                f_cov.write(cov_q.tobytes(order='C'))
                f_cov_scale.write(struct.pack('f', float(scale)))
            else:
                # Byte-for-byte what ``struct.pack(f'{n}f', *cov_upper)`` wrote:
                # both are little-endian IEEE-754 binary32 with round-to-nearest,
                # but this one does not build an n-element argument tuple and
                # make n individual float conversions per variant.
                f_cov.write(cov_upper.astype('<f4', copy=False).tobytes())
            f_cov_ids.write('\t'.join(map(str, info)) + '\n')

            if args.output_rtau:
                all_betas.append(beta)


def write_stats_batch_fused(fused_list, batch_info, args, f_tsv_list, f_cov_list, f_cov_ids_list,
                            f_cov_scale_list, all_betas_list, n_cov_elements, n_taus,
                            d_resid_batch=None):
    """Write one batch from the fused C++ finalizer output.

    Byte-for-byte the same output layout as write_stats_batch; only the arithmetic that
    produced beta/se/Sigma differs (C++ rather than NumPy).  The skip rule is the same
    test in a different place: write_stats_batch skips a SNP whose block statistics are
    NaN, and the fused path returns NaN beta/se/cov for exactly those SNPs (the C++
    genotype filter guarantees sum(g_resid^2) >= min_var > 0, so the only way to reach a
    non-positive jackknife denominator is to have failed that filter already).

    The one column that is not bit-identical between the two paths is d_resid. The
    fused kernel never returns the per-block statistics that carry it, so it is
    recomputed in NumPy from the same genotypes and the same covariate basis. The two
    differ only in summation order, which matters nowhere: d_resid feeds the degrees of
    freedom of a t reference, where a last-bit difference is invisible.
    """
    want_cov = args.cov_mode not in COV_MODES_NO_PAYLOAD
    for p_idx, (betas, se_all, cov_all) in enumerate(fused_list):
        f_tsv = f_tsv_list[p_idx]
        f_cov = f_cov_list[p_idx]
        f_cov_ids = f_cov_ids_list[p_idx]
        f_cov_scale = f_cov_scale_list[p_idx]
        all_betas = all_betas_list[p_idx]

        for j, info in enumerate(batch_info):
            if np.isnan(betas[j, 0]):
                continue

            beta = betas[j]
            se = se_all[j]
            cov_upper = cov_all[j]

            row = list(info)
            for b, s in zip(beta, se):
                row.extend([f'{b:.6g}', f'{s:.6g}'])
            d_resid = float(d_resid_batch[j]) if d_resid_batch is not None else float('nan')
            row.append(f'{d_resid:.9g}')
            f_tsv.write('\t'.join(map(str, row)) + '\n')

            if not want_cov:
                continue

            if args.cov_dtype == 'int8':
                scale = np.max(np.abs(cov_upper)) / 127.0
                if not np.isfinite(scale) or scale <= 0:
                    scale = 1.0
                cov_q = np.clip(np.rint(cov_upper / scale), -127, 127).astype(np.int8)
                f_cov.write(cov_q.tobytes(order='C'))
                f_cov_scale.write(struct.pack('f', float(scale)))
            else:
                # Same bytes as struct.pack(f'{n}f', *cov_upper); see
                # write_stats_batch.
                f_cov.write(np.asarray(cov_upper).astype('<f4', copy=False).tobytes())
            f_cov_ids.write('\t'.join(map(str, info)) + '\n')

            if args.output_rtau:
                all_betas.append(np.asarray(beta))


def write_stats_batch_pooled(pooled_list, batch_info, args, f_tsv_list,
                             sigma0_sd_list, all_betas_list):
    """Write one batch from the pooled-only kernel.

    Same TSV layout, same rows, same row order and the same skip rule as
    :func:`write_stats_batch` -- a variant is written iff it passed the MAC and
    variance filters, which the pooled kernel signals with a finite ``d_resid``
    exactly as the block kernel signals it with non-NaN block statistics.

    Two columns carry different numbers, and both differences are intended:

    ``beta_tau_*`` and ``d_resid`` are the same estimands summed in a different
    order.  The block path forms the totals as ``sum_b (sum_{i in b} ...)``, the
    pooled path as ``sum_i ...``; same addends, different association, so they
    agree to rounding rather than to the bit.  At ``%.6g`` and ``%.9g`` that is
    invisible except for a value that happens to sit on a formatting boundary.

    ``se_tau_*`` are the **pooled** standard errors ``sqrt(Sigma_0[t,t]/D)``,
    not jackknife estimates.  There is no cheap jackknife: the delete-one-block
    standard error is a function of the per-block statistics and nothing else,
    so a path that does not accumulate blocks cannot produce it.  What it can
    produce is the diagonal of the covariance Stage 2 actually applies under
    ``cov_source = "pooled"``, which is the self-consistent choice -- and, on
    the measurement in ``pooled_cov.py``, the better-behaved one (diagonal ratio
    against an exact HC0 sandwich, 5th-95th percentile: 0.990-1.010 pooled
    against 0.725-1.307 jackknife).

    VERIFIED on chr22 testosterone with LOCO, all 17,550 variants, against a
    full-mode run whose own output is byte-identical (decompressed) to the
    accepted reference scan:

      * ``.sigma0.tsv.gz``    byte-identical
      * identity columns      identical, 87,750/87,750
      * ``d_resid``           identical, 17,550/17,550
      * ``beta_tau_*``        identical on 789,749 of 789,750; the one that
                              differs sits on a ``%.6g`` rounding boundary
                              (relative difference 2.0e-6 in the printed value)
      * ``se_tau_*``          differ on 99.99%, median relative difference 1.5%,
                              maximum 25% -- this is the jackknife against the
                              pooled estimator, not an error

    Through Stage 2 (skew, ``cov_source="pooled"``, ``tau_cov_estimator="ols"``,
    the accepted gate's own weights) every parameter, standard error, p-value
    and degrees-of-freedom column matches the gate on 17,549 of 17,550 rows and
    agrees to 4e-11 on the remaining one; the genome-wide-significant counts are
    unchanged (20 ``beta_left``, 27 ``beta_right``).  What moves is ``Q`` and
    its calibration -- ``Q``, ``Q_df_eff``, ``Q_scale`` on every row, ``Q_pval``
    on 17,251 -- because ``Q = sum_t resid_t^2 / se_tau_t^2`` is the one place
    Stage 2's pooled path reads ``se_tau_*`` (R/gwas.R, the
    ``tau_cov_estimator == "ols"`` branch).  Maximum relative change: ``Q``
    21.7%, ``Q_pval`` 60.0%, ``|dlog10 Q_pval|`` 0.40.  No ``Q_pval`` crosses
    5e-8 in either run.  Whether the pooled-diagonal ``Q`` is the better
    calibrated statistic under a pooled covariance is a real question and is
    NOT settled here; it is only noted that it is the self-consistent one.
    """
    for p_idx, (betas, d_resid) in enumerate(pooled_list):
        f_tsv = f_tsv_list[p_idx]
        sd0 = sigma0_sd_list[p_idx]
        all_betas = all_betas_list[p_idx]

        d_resid = np.asarray(d_resid, dtype=np.float64)
        finite = np.isfinite(d_resid)
        with np.errstate(invalid='ignore', divide='ignore'):
            se_all = sd0[None, :] / np.sqrt(d_resid)[:, None]

        for j, info in enumerate(batch_info):
            if not finite[j]:
                continue
            row = list(info)
            for b, s in zip(betas[j], se_all[j]):
                row.extend([f'{b:.6g}', f'{s:.6g}'])
            row.append(f'{float(d_resid[j]):.9g}')
            f_tsv.write('\t'.join(map(str, row)) + '\n')

            if args.output_rtau:
                all_betas.append(np.asarray(betas[j]))


class _SkipCount:
    """Mutable counter so a generator can report variants it dropped."""

    __slots__ = ("non_biallelic",)

    def __init__(self):
        self.non_biallelic = 0


def _decode_variants(variants, bgen_path, bgen_indices, skipped: _SkipCount):
    """Yield ``(info, variant_chr, dosages)`` for each biallelic variant, in file order.

    This is exactly the per-variant work the batch loop used to do inline; it is
    split out so it can be run ahead of the kernel on another thread.
    """
    for variant in variants:
        alleles = variant.alleles
        if len(alleles) != 2:
            skipped.non_biallelic += 1
            continue
        bit_depth = _variant_bit_depth(variant, bgen_path=bgen_path)
        dosages = extract_dosage(variant, probability_bit_depth=bit_depth)[bgen_indices]
        info = (variant.rsid, variant.chrom, variant.pos, alleles[0], alleles[1])
        yield info, normalize_chr_label(variant.chrom), dosages


class _BatchPipeline:
    """Overlap BGEN decoding with the kernel by running batches on a worker thread.

    Decoding a BGEN record is single-threaded Python that the OpenMP kernel
    cannot help with, and it was a large share of Stage 1's wall time while the
    worker threads sat idle (a full chr22 run averaged 2.4 of 8 allocated cores).
    The kernel releases the GIL, so the main thread can decode batch k+1 while
    batch k is being scored and written.

    Ordering and arithmetic are unchanged.  A single worker runs one batch at a
    time, in submission order, so the output files are written in exactly the
    order the serial loop wrote them; genotypes are double-buffered so the batch
    in flight is never the one being filled.

    The reader stays on the calling thread throughout: ``bgen``'s .bgi index is a
    sqlite connection, and sqlite objects may only be used by the thread that
    created them.
    """

    def __init__(self, enabled: bool, n_samples: int, batch_size: int):
        self.enabled = enabled
        self._pending = None
        self._pool = None
        self.buffers = [np.zeros((n_samples, batch_size), dtype=np.float64, order='F')]
        if enabled:
            from concurrent.futures import ThreadPoolExecutor
            self._pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix="stage1-batch")
            # Second genotype buffer so the worker can read one while the reader
            # fills the other.
            self.buffers.append(
                np.zeros((n_samples, batch_size), dtype=np.float64, order='F')
            )
        self._buf_idx = 0

    @property
    def batch_buffer(self):
        return self.buffers[self._buf_idx]

    def submit(self, fn, *args, **kwargs):
        """Run one batch, returning ``(compute_s, write_s)`` of a *previous* batch.

        Returns ``(0.0, 0.0)`` when nothing had completed yet.
        """
        if not self.enabled:
            return fn(*args, **kwargs)
        timings = self.drain()
        self._pending = self._pool.submit(fn, *args, **kwargs)
        self._buf_idx ^= 1
        return timings

    def drain(self):
        """Wait for the in-flight batch, if any, and return its timings."""
        if self._pending is None:
            return 0.0, 0.0
        pending, self._pending = self._pending, None
        return pending.result()

    def close(self):
        if self._pool is not None:
            self._pool.shutdown(wait=True)
            self._pool = None


def _process_batch(G_view, RIF_resid_list, Q, block_ids, batch_info, args, use_multi,
                   f_tsv_list, f_cov_list, f_cov_ids_list, f_cov_scale_list,
                   all_betas_list, n_cov_elements, T, sigma0_sd_list=None):
    """Run the kernel over one batch and write it out. Returns (compute_s, write_s)."""
    if args.cov_mode == COV_MODE_POOLED_ONLY:
        # No block_ids, no n_blocks, no (n_blocks, panel, 1+T) accumulator: the
        # kernel returns the two things the pooled analysis needs.
        t_compute = time.perf_counter()
        pooled_list = [
            compute_pooled_panel(G_view, rif, Q, n_threads=args.threads,
                                 panel_size=args.panel_size)
            for rif in RIF_resid_list
        ]
        compute_s = time.perf_counter() - t_compute
        t_write = time.perf_counter()
        write_stats_batch_pooled(pooled_list, batch_info, args, f_tsv_list,
                                 sigma0_sd_list, all_betas_list)
        return compute_s, time.perf_counter() - t_write

    t_compute = time.perf_counter()
    if use_multi:
        stats_list = compute_block_scores_multi(
            G_view, RIF_resid_list, Q, block_ids, args.blocks, n_threads=args.threads
        )
        fused_list = None
    elif args.fused_finalize:
        finalize = compute_block_finalize_panel if args.panel_kernel else compute_block_finalize
        kw = {'panel_size': args.panel_size} if args.panel_kernel else {}
        fused_list = [
            finalize(
                G_view, RIF_resid_list[0], Q, block_ids, args.blocks,
                n_threads=args.threads, **kw
            )
        ]
        stats_list = None
    else:
        scores = compute_block_scores_panel if args.panel_kernel else compute_block_scores
        kw = {'panel_size': args.panel_size} if args.panel_kernel else {}
        stats_list = [
            scores(
                G_view, RIF_resid_list[0], Q, block_ids, args.blocks,
                n_threads=args.threads, **kw
            )
        ]
        fused_list = None
    compute_s = time.perf_counter() - t_compute

    t_write = time.perf_counter()
    if fused_list is not None:
        write_stats_batch_fused(
            fused_list, batch_info, args, f_tsv_list, f_cov_list, f_cov_ids_list,
            f_cov_scale_list, all_betas_list, n_cov_elements, T,
            d_resid_batch=residual_sum_squares(G_view, Q)
        )
    else:
        write_stats_batch(
            stats_list, batch_info, args, f_tsv_list, f_cov_list, f_cov_ids_list,
            f_cov_scale_list, all_betas_list, n_cov_elements, T
        )
    write_s = time.perf_counter() - t_write
    return compute_s, write_s


def _loco_column_for_chr(loco_predictions: dict[str, dict[str, np.ndarray]],
                         pheno_cols: list[str], chrom_key: str,
                         valid_mask: np.ndarray) -> np.ndarray:
    """The chromosome's LOCO prediction on the analysed rows, as one basis column.

    Under ``--loco-mode covariate`` every phenotype in a multi-phenotype run
    shares one residualising basis ``Q``, because ``Q`` also projects the
    genotype.  REGENIE predictions are per trait, so a multi-phenotype run
    would need one basis per trait; that is refused upstream in ``main``.
    """
    col = pheno_cols[0]
    per_chr = loco_predictions[col]
    if chrom_key not in per_chr:
        raise ValueError(f"Chromosome {chrom_key} missing from LOCO predictions for phenotype '{col}'")
    return np.asarray(per_chr[chrom_key], dtype=np.float64)[valid_mask]


def compute_chr_rif_residuals(pheno_cols: list[str], y_values_full: dict[str, np.ndarray],
                              X_full: np.ndarray, taus: np.ndarray,
                              loco_predictions: Optional[dict[str, dict[str, np.ndarray]]],
                              chrom: Optional[str], valid_mask: np.ndarray,
                              loco_scales: Optional[dict[str, float]] = None,
                              density_method: str = DENSITY_METHOD_SILVERMAN_GAUSSIAN,
                              loco_mode: str = LOCO_MODE_SUBTRACT
                              ) -> tuple[list[np.ndarray], np.ndarray, list[np.ndarray]]:
    """Build residualized RIF matrices for a chromosome, handling LOCO per ``loco_mode``.

    ``subtract`` (historical default) forms ``Y - c*L(chr)`` and builds the RIF on
    that.  ``covariate`` leaves ``Y`` alone and appends ``L(chr)`` to the
    residualising basis instead, so the RIF is built on the raw unconditional
    distribution of ``Y`` and each influence-function column -- and, through the
    same ``Q``, the genotype -- is residualised on ``[X, L(chr)]``.

    The distinction matters because a RIF/UQR estimand is defined against the
    *unconditional* distribution of ``Y``: ``RIF_i(tau) = q_tau + (tau -
    1{Y_i <= q_tau}) / f(q_tau)`` depends on ``Y`` only through each sample's
    position relative to the population's tau-th quantile.  Changing ``Y`` by a
    genotype-independent term therefore moves ``q_tau``, ``f(q_tau)`` and every
    indicator, i.e. it redefines the target rather than removing a nuisance.
    That is unlike a mean GWAS, where the same subtraction leaves the slope
    untouched and only shrinks residual variance.
    """
    if loco_predictions is None:
        y_values_chr = {col: y_values_full[col][valid_mask] for col in pheno_cols}
        X_chr = X_full[valid_mask]
    elif loco_mode == LOCO_MODE_COVARIATE:
        chrom_key = normalize_chr_label(chrom)
        y_values_chr = {col: y_values_full[col][valid_mask] for col in pheno_cols}
        # The scale factor is deliberately not applied: a covariate's column span
        # -- and hence the projection -- is invariant to rescaling it, so
        # --loco-scale has no effect in this mode.
        loco_col = _loco_column_for_chr(loco_predictions, pheno_cols, chrom_key, valid_mask)
        X_chr = np.column_stack([X_full[valid_mask], loco_col])
    else:
        chrom_key = normalize_chr_label(chrom)
        y_values_chr = {}
        for col in pheno_cols:
            per_chr = loco_predictions[col]
            if chrom_key not in per_chr:
                raise ValueError(f"Chromosome {chrom_key} missing from LOCO predictions for phenotype '{col}'")
            scale = 1.0 if loco_scales is None else loco_scales.get(col, 1.0)
            y_values_chr[col] = y_values_full[col] - scale * per_chr[chrom_key]
        y_values_chr = {col: y_values_chr[col][valid_mask] for col in pheno_cols}
        X_chr = X_full[valid_mask]

    rif_resid_list, q_chr, q_tau_list = build_rif_residuals(
        y_values_chr, pheno_cols, X_chr, taus,
        density_method=density_method
    )
    return rif_resid_list, q_chr, q_tau_list


def main():
    parser = argparse.ArgumentParser(
        description="FungWas Stage 1: Fast RIF quantile GWAS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Required arguments
    parser.add_argument('--bgen', required=True,
                        help='Path to BGEN file')
    parser.add_argument('--pheno', required=True,
                        help='Path to phenotype file (space-delimited, with FID IID and phenotype column)')
    parser.add_argument('--pheno-col', dest='pheno_col', required=False,
                        help='Name of the phenotype column to analyze')
    parser.add_argument('--pheno-cols', dest='pheno_cols', default=None,
                        help='Comma-separated phenotype columns to analyze in one pass '
                             '(uses only samples with non-missing values across all phenotypes/covariates)')
    parser.add_argument('--out', required=True,
                        help='Output prefix (creates {out}.stage1.tsv.gz, {out}.cov.gz)')
    
    # Optional arguments
    parser.add_argument('--sample', default=None,
                        help='Path to sample file (if not embedded in BGEN)')
    parser.add_argument('--covar', default=None,
                        help='Path to covariate file')
    parser.add_argument('--covar-col', dest='covar_cols', nargs='+', default=None,
                        help='Covariate column(s) to use (space or comma-separated). '
                             'If omitted, all non-ID columns in the covariate file are used.')
    parser.add_argument('--snps', default=None,
                        help='File with SNP IDs to process (one per line)')
    parser.add_argument('--keep', default=None,
                        help='File with sample IDs to keep (one per line)')
    parser.add_argument('--min-sample-match', dest='min_sample_match',
                        type=float, default=0.01,
                        help="Refuse to run if fewer than this fraction of phenotype "
                             "samples match the BGEN (default 0.01). Catches an "
                             "ID-convention mismatch, which produces a near-total "
                             "failure to match and would otherwise run silently on "
                             "whatever happened to collide. A partial overlap is "
                             "legitimate -- a phenotype file routinely holds samples "
                             "absent from a given BGEN -- so anything between this "
                             "threshold and 50%% warns rather than refuses. Set to 0 "
                             "to disable entirely.")
    parser.add_argument('--loco-scale', dest='loco_scale', default='regression',
                        help="How to put REGENIE .loco predictions on the phenotype's scale. "
                             "'regression' (default) uses the OLS slope of Y on the prediction, "
                             "pooled over chromosomes; 'sd' multiplies by sd(Y), which is only "
                             "exact when REGENIE's step-1 covariates explain little of Y; "
                             "'none' subtracts as-is (pre-v0.4.0 behaviour); or give a number "
                             "to use it verbatim.")
    parser.add_argument('--loco-mode', dest='loco_mode', default=None,
                        choices=list(LOCO_MODES),
                        help="Where the REGENIE .loco prediction enters the model. "
                             "'covariate' (DEFAULT since 0.3.0) leaves Y unadjusted and "
                             "appends L(chr) to the residualising basis, so each RIF "
                             "column and the genotype are residualised on [X, L(chr)]. "
                             "'subtract' is the historical behaviour, building the RIF on "
                             "Y - c*L(chr); it is retained for reproducing earlier runs. "
                             "RIF/UQR estimands are defined against the UNCONDITIONAL "
                             "distribution of Y, so subtracting L changes q_tau, f(q_tau) "
                             "and every indicator: it redefines the target rather than "
                             "removing a nuisance. Measured on testosterone chr22, "
                             "subtracting inflated SEs by 9.5-14.7%% and moved the median "
                             "effect by a third of a standard error. --loco-scale is "
                             "inert under 'covariate': a projection is invariant to "
                             "rescaling the column.")
    parser.add_argument('--loco-preds', dest='loco_preds', default=None,
                        help='REGENIE step1 LOCO predictions (for future use)')
    parser.add_argument('--allow-numpy-fallback', action='store_true',
                        help='Allow the uncompiled NumPy fallback. This can be ~1000x slower and is for debugging only.')
    parser.add_argument('--no-panel-kernel', dest='panel_kernel', action='store_false',
                        help='Disable the panel block-score kernel and use the per-SNP form. '
                             'Output is bit-identical either way; the panel form is about 4x '
                             'faster. Ignored for multi-phenotype runs, which always use the '
                             'per-SNP form.')
    parser.add_argument('--panel-kernel', dest='panel_kernel', action='store_true',
                        help='Use the panel form of the block-score kernel: carry --panel-size '
                             'SNPs through one sweep over the samples so the residualised RIF '
                             'matrix is streamed once per panel instead of once per SNP. '
                             'Bit-identical output; single phenotype only.')
    parser.set_defaults(panel_kernel=True)
    parser.add_argument('--panel-size', dest='panel_size', type=int, default=16,
                        help='Number of SNPs per panel for --panel-kernel (default: 16).')
    parser.add_argument('--fused-finalize', dest='fused_finalize', action='store_true',
                        help='Compute beta/SE/covariance inside the C++ OpenMP loop instead of '
                             'returning per-block statistics and redoing the work per SNP in NumPy. '
                             'Single-phenotype only. Results differ from the NumPy path in the last '
                             'bits (different summation order).')
    parser.add_argument('--no-seek-variants', dest='seek_variants', action='store_false',
                        help='Walk the BGEN record-by-record instead of seeking to .bgi file '
                             'offsets (slower; kept for validation and for BGENs with no index).')
    parser.set_defaults(seek_variants=True)
    parser.add_argument('--allow-multiallelic-rsids', action='store_true',
                        help='Allow rsIDs that map to multiple BGEN records. By default these are excluded using the .bgi index.')
    
    # Algorithm parameters
    default_taus = ",".join([f"{t:.3f}" for t in np.linspace(0.02, 0.98, 45)])
    parser.add_argument('--taus', default=default_taus,
                        help='Comma-separated tau values (default: 45 taus from 0.02 to 0.98)')
    parser.add_argument('--density-method', choices=DENSITY_METHODS,
                        default=DENSITY_METHOD_SILVERMAN_GAUSSIAN,
                        help='Kernel for the Stage 1 RIF density (default: '
                             'silverman_gaussian). fourth_order_gaussian is an '
                             'optional fourth-order Gaussian kernel using the same '
                             'Silverman bandwidth.')
    parser.add_argument('--blocks', type=int, default=None,
                        help=f'Number of jackknife blocks (default: {DEFAULT_N_BLOCKS}). '
                             'Must exceed the number of taus. Fewer blocks make the fixed-RIF '
                             'covariance and Stage 2 t calibration noisier and remain approximate.')
    parser.add_argument('--batch-size', dest='batch_size', type=int, default=500,
                        help='SNPs per batch for memory efficiency')
    parser.add_argument('--threads', type=int, default=1,
                        help='OpenMP threads for C++ kernel')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for block assignment')
    
    # Output options
    parser.add_argument('--output-rtau', dest='output_rtau', action='store_true',
                        help='Also output R_tau correlation matrix for plugin_cor')
    parser.add_argument('--no-sigma0', dest='sigma0', action='store_false',
                        help='Do not write the pooled cross-quantile covariance sidecar '
                             '(<prefix>.sigma0.tsv.gz + .sigma0.meta.json). The sidecar '
                             'is one 45x45 matrix per chromosome, costs well under a '
                             'second to compute, and is what lets Stage 2 use the pooled '
                             'covariance instead of the per-SNP jackknife.')
    parser.set_defaults(sigma0=True)
    parser.add_argument('--cov-dtype', dest='cov_dtype', choices=['int8', 'float32'],
                        default='int8',
                        help='Covariance storage dtype: int8 (with .cov.scale.gz sidecar) '
                             'or float32 (legacy). Default: int8')
    parser.add_argument('--cov-mode', dest='cov_mode', choices=list(COV_MODES),
                        default=COV_MODE_FULL,
                        help="Which covariance artefacts to write. 'full' (default) "
                             "writes the per-SNP jackknife payload (.cov.gz, "
                             ".cov.ids.tsv.gz, and .cov.scale.gz for int8), which "
                             "Stage 2 reads only under cov_source='jackknife'. "
                             "'pooled' omits that payload and writes only what "
                             "Stage 2's production path (cov_source='pooled', "
                             "tau_cov_estimator='ols') actually reads: the same "
                             ".stage1.tsv.gz, byte for byte, and the same Sigma_0 "
                             "sidecar. The jackknife still runs and se_tau_* are "
                             "unchanged. The omitted payload cannot be rebuilt "
                             "without re-running Stage 1, so use 'pooled' only when "
                             "no downstream analysis will ask for the jackknife "
                             "covariance. "
                             "'pooled-only' additionally skips the jackknife "
                             "COMPUTATION: no per-block accumulator is built, so "
                             "beta_tau_* and d_resid come from single-pass totals "
                             "and se_tau_* are the pooled standard errors "
                             "sqrt(Sigma_0[t,t]/d_resid) rather than jackknife "
                             "estimates. beta_tau_* and d_resid agree with 'full' "
                             "to rounding (different summation order, same "
                             "addends); se_tau_* are a different estimator, which "
                             "changes Stage 2's Q lack-of-fit statistic even "
                             "though it leaves the parameter estimates and their "
                             "standard errors untouched. Requires the Sigma_0 "
                             "sidecar.")
    parser.add_argument('--no-jackknife-cov', dest='cov_mode', action='store_const',
                        const=COV_MODE_POOLED,
                        help='Alias for --cov-mode pooled.')
    parser.add_argument('--no-pipeline', dest='pipeline', action='store_false',
                        help='Score and write each batch inline instead of on a '
                             'background thread. Output is identical either way; '
                             'the pipeline just stops the OpenMP threads idling '
                             'while the next batch of BGEN records is decoded.')
    parser.add_argument('--pipeline', dest='pipeline', action='store_true',
                        help='Overlap BGEN decoding with the kernel by scoring and '
                             'writing each batch on a background thread while the '
                             'next batch is read (default). Batches are processed '
                             'one at a time in order, so output order is unchanged; '
                             'the cost is one extra genotype buffer '
                             '(batch_size * n_samples * 8 bytes).')
    parser.set_defaults(pipeline=True)
    
    args = parser.parse_args()
    set_numpy_fallback_allowed(args.allow_numpy_fallback)
    require_cpp_extension("fungwas-stage1 CLI")

    if not args.pheno_col and not args.pheno_cols:
        parser.error("Provide --pheno-col or --pheno-cols")
    
    # Parse taus
    taus = np.array([float(t) for t in args.taus.split(',')])
    T = len(taus)
    try:
        args.blocks = resolve_n_blocks(args.blocks, T, flag="--blocks")
    except ValueError as exc:
        parser.error(str(exc))
    
    log_mem(f"Starting FungWas Stage 1 (taus={T}, blocks={args.blocks}, "
            f"density_method={args.density_method})")
    logger.info("Note: For Stage 2, consider mean-centering non-intercept columns of your W matrix to reduce collinearity.")
    t_start = time.time()
    t_load = time.perf_counter()
    
    if args.pheno_cols:
        pheno_cols = [c.strip() for c in args.pheno_cols.split(',') if c.strip()]
    else:
        pheno_cols = [args.pheno_col]
    covar_cols = parse_list_arg(args.covar_cols)

    use_multi = len(pheno_cols) > 1
    if args.panel_kernel and len(pheno_cols) > 1:
        # The panel kernel has no multi-phenotype counterpart. Now that it is the
        # default, fall back silently rather than failing a valid --pheno-cols run.
        args.panel_kernel = False
    if False:
        raise ValueError("--panel-kernel supports a single phenotype only; "
                         "drop it or use a single --pheno-col.")
    if args.cov_mode == COV_MODE_POOLED_ONLY and not args.sigma0:
        parser.error(
            "--cov-mode pooled-only needs the Sigma_0 sidecar: it is where "
            "se_tau_* comes from, and it is the only covariance the run emits. "
            "Drop --no-sigma0, or use --cov-mode pooled."
        )
    if args.cov_mode == COV_MODE_POOLED_ONLY and args.fused_finalize:
        parser.error(
            "--fused-finalize is a jackknife finalisation; --cov-mode pooled-only "
            "does not compute the jackknife. Drop one of the two."
        )

    if args.panel_kernel and args.panel_size < 1:
        raise ValueError("--panel-size must be >= 1")
    # --loco-mode defaults to 'covariate' (0.3.0), but only bites when there is a
    # prediction to place. With no --loco-preds the mode is moot: nothing is
    # subtracted and nothing is added to the basis. Only refuse when the user
    # explicitly asked for covariate mode and gave nothing to put in it.
    if args.loco_mode is None:
        args.loco_mode = LOCO_MODE_COVARIATE if args.loco_preds else LOCO_MODE_SUBTRACT
    elif args.loco_mode == LOCO_MODE_COVARIATE and not args.loco_preds:
        raise ValueError("--loco-mode covariate needs --loco-preds; without a "
                         "prediction there is no column to add to the basis.")
    if args.fused_finalize and len(pheno_cols) > 1:
        raise ValueError("--fused-finalize supports a single phenotype only; "
                         "drop it or use one --pheno-col.")
    if not use_multi and args.pheno_cols:
        log_mem("Single phenotype provided via --pheno-cols; using single-phenotype path.")

    if use_multi:
        log_mem(f"Loading {len(pheno_cols)} phenotypes in one pass...")
    sample_ids, sample_aliases, y_values_raw, X_full, covar_cols_used = load_analysis_inputs(
        args.pheno, pheno_cols, args.covar, keep_file=args.keep, covar_cols=covar_cols
    )
    y_values_full = {col: np.asarray(values, dtype=np.float64) for col, values in y_values_raw.items()}
    loco_predictions = None
    loco_files = {col: None for col in pheno_cols}
    loco_scales = None
    if args.loco_preds:
        # Resolved here as well as inside load_loco_predictions so the pooled
        # covariance sidecar can name the predictor that produced it. A Sigma_0
        # built without the LOCO adjustment, or with another trait's, is wrong in
        # exactly the silent way the covariate set is.
        if args.loco_mode == LOCO_MODE_COVARIATE and len(pheno_cols) > 1:
            raise ValueError(
                "--loco-mode covariate supports a single phenotype per run. The "
                "residualising basis Q is shared across phenotypes because it also "
                "projects the genotype, but REGENIE LOCO predictions are per trait, "
                "so several traits would need several bases. Run them separately."
            )
        loco_files = resolve_loco_files(args.loco_preds, pheno_cols)
        loco_predictions, loco_valid_mask = load_loco_predictions(args.loco_preds, pheno_cols, sample_aliases)
        loco_scales = None  # filled in after the sample mask is applied
        if not np.all(loco_valid_mask):
            dropped = int((~loco_valid_mask).sum())
            log_mem(f"Dropping {dropped} phenotype rows not present in LOCO predictions")
            sample_ids = [sid for sid, keep in zip(sample_ids, loco_valid_mask) if keep]
            sample_aliases = [aliases for aliases, keep in zip(sample_aliases, loco_valid_mask) if keep]
            y_values_full = {col: values[loco_valid_mask] for col, values in y_values_full.items()}
            X_full = X_full[loco_valid_mask]

        # Computed on the LOCO-matched phenotype rows, which is a superset of the rows
        # that survive the later BGEN sample match. The factor is a property of the
        # phenotype's scale, not of the analysed subset, so the difference is immaterial
        # (0.2% for log_bmi: 484,460 rows vs 461,320).
        if args.loco_mode == LOCO_MODE_COVARIATE:
            # No scale is resolved: the projection onto span([X, L]) does not
            # change if L is rescaled, so resolve_loco_scale would only emit a
            # factor and a variance warning that describe an adjustment this
            # mode never makes.
            loco_scales = None
            log_mem("LOCO predictions loaded; Y is left unadjusted and L(chr) enters "
                    "the residualising basis per chromosome (--loco-mode covariate)")
        else:
            loco_scales = {
                col: resolve_loco_scale(args.loco_scale, y_values_full[col], loco_predictions[col], col)
                for col in pheno_cols
            }
            log_mem("LOCO predictions loaded; RIFs will be recomputed per chromosome on Y - LOCO(chr)")
    # NB: the RIF and the covariate basis are NOT built here. They are built after
    # the BGEN sample match below, on exactly the rows that will be analysed. See
    # the comment at that call site.
    log_mem(f"Phenotype load time: {time.perf_counter() - t_load:.2f}s")
    N = len(sample_ids)
    
    # Setup output dir
    out_dir = os.path.dirname(args.out) or '.'
    os.makedirs(out_dir, exist_ok=True)
    
    # Load SNP list if provided
    snp_list = None
    if args.snps:
        snp_list = load_snps_list(args.snps)
        log_mem(f"Loaded {len(snp_list)} SNPs from {args.snps}")
    
    # Open BGEN and get variants
    variants, bgen_samples, bfile_handle = get_variants_bgen(
        args.bgen, args.sample, snp_list, args.allow_multiallelic_rsids,
        seek_variants=args.seek_variants
    )

    # bgen does not expose bit depth on its Python variant object.  It is a per-record
    # field, so each selected variant gets a bounded preamble probe; production 8-bit
    # records still use the native vectorised dosage path after that small check.
    variants = iter(variants)
    
    # Match samples
    valid_mask, bgen_indices = match_samples(bgen_samples, sample_ids)
    N_matched = len(bgen_indices)

    log_mem(f"Matched {N_matched}/{N} samples to BGEN")

    # A partial sample match is the dangerous case, not a total one. Zero overlap
    # fails loudly downstream; a 3% overlap runs to completion on a small,
    # non-random subset and produces attenuated effects with inflated SEs that are
    # internally consistent and wrong. match_samples keys on the phenotype file's
    # FID alone, while the LOCO reader resolves five aliases (fid, iid, fid_iid,
    # "fid iid", fid:iid), so a cohort whose BGEN IDs are FID_IID matches LOCO and
    # then silently drops every sample here.
    match_frac = (N_matched / N) if N else 0.0
    if args.min_sample_match > 0 and 0.0 < match_frac < 0.50 \
            and match_frac >= args.min_sample_match:
        logger.warning(
            "Only %d of %d phenotype samples (%.1f%%) matched the BGEN. That is "
            "legitimate when the phenotype file covers a wider cohort than this "
            "genotype release, but it is also what a partial ID-convention "
            "collision looks like. Check that the analysed subset is the one you "
            "intended.", N_matched, N, 100 * match_frac)
    if args.min_sample_match > 0 and match_frac < args.min_sample_match:
        preview = ", ".join(str(s) for s in sample_ids[:3])
        bgen_preview = ", ".join(str(s) for s in bgen_samples[:3])
        raise ValueError(
            f"Only {N_matched} of {N} phenotype samples ({100 * match_frac:.1f}%) "
            f"matched the BGEN, below --min-sample-match "
            f"({100 * args.min_sample_match:.0f}%).\n"
            f"  phenotype IDs look like: {preview}\n"
            f"  BGEN sample IDs look like: {bgen_preview}\n"
            "A failure this near-total means the two files use different identifier "
            "conventions rather than covering different samples. Sample matching uses "
            "the phenotype file's first column (FID) against the BGEN sample "
            "identifiers with NO alias resolution -- unlike the LOCO reader, which "
            "resolves fid, iid, fid_iid, 'fid iid' and fid:iid, so a mismatched "
            "cohort can match LOCO and then drop every sample here. Rewrite the "
            "phenotype file's first column to match the BGEN, or pass "
            "--min-sample-match 0 to proceed deliberately."
        )

    # Build the covariate basis and the RIF on the MATCHED rows, not on the keep
    # sample. Both depend on which rows are in the analysis, and neither survives
    # being row-subset afterwards:
    #
    #   * Q is used by the kernels as a projector, g_resid = g - Q (Q' g). Columns
    #     orthonormal over the keep sample are not orthonormal over a subset of it,
    #     so the projection becomes partial and leaks the covariate-parallel part of
    #     the genotype into the numerator of every slope. On the UK Biobank BMI
    #     inputs (486,254 keep rows, 461,320 matched) the row-subset basis reached
    #     max |Q'Q - I| = 0.97 - the ancestry-PC columns kept 3-9% of their norm,
    #     because the unmatched samples are the ones the genotype release excludes
    #     on ancestry - and shifted beta_tau by a median of 0.12 SE, up to 2.3 SE,
    #     in a consistent direction at all 45 taus, while leaving se_tau unchanged.
    #   * the RIF depends on the sample through q_tau, the bandwidth and f_hat.
    #     That channel is second-order here (mean delta-z -0.0008) because q_tau is
    #     an additive constant the intercept removes and f_hat largely cancels
    #     between beta and se, but building on the matched rows costs nothing and
    #     removes it too.
    #
    # This is what the --loco-preds path has always done via compute_chr_rif_residuals,
    # which is why no published result was affected; the call below is that same
    # function with no LOCO adjustment.
    RIF_resid_list, Q, q_tau_list = compute_chr_rif_residuals(
        pheno_cols, y_values_full, X_full, taus, None, None, valid_mask,
        density_method=args.density_method
    )

    # Pooled cross-quantile covariance, Sigma_0 = Psi_r' Psi_r / n, accumulated
    # from the same residualised RIF the kernel is about to consume -- so it is
    # a by-product of work already done, and it is built on the post-sample-match
    # objects by construction (the defect fixed at 4754794 was a basis built on
    # the keep sample and applied to a subset).
    #
    # Without --loco-preds the RIF does not depend on the chromosome, so one
    # matrix under the key "ALL" covers the whole run. With it, the phenotype is
    # Y - LOCO(chr) and the matrices are collected per chromosome in the variant
    # loop below; the build here is then thrown away before any variant is
    # scored, so it must not be recorded.
    sigma0_by_chr = [dict() for _ in pheno_cols]
    sigma0_seconds = 0.0
    # sqrt(diag(Sigma_0)) for the RIF currently in force, which is what
    # --cov-mode pooled-only writes into se_tau_* (divided by sqrt(d_resid)).
    # Rebound, never mutated, whenever the LOCO adjustment changes the RIF, so
    # a batch already in flight on the pipeline thread keeps the list it was
    # handed -- the same discipline RIF_resid_list follows.
    sigma0_sd_list = None
    if args.sigma0 and loco_predictions is None:
        t_s0 = time.perf_counter()
        for p_idx, rif in enumerate(RIF_resid_list):
            sigma0_by_chr[p_idx][SIGMA0_ALL_CHROM] = compute_sigma0(rif)
        sigma0_sd_list = [np.sqrt(np.diag(sigma0_by_chr[p_idx][SIGMA0_ALL_CHROM]))
                          for p_idx in range(len(pheno_cols))]
        sigma0_seconds += time.perf_counter() - t_s0
        log_mem(f"Sigma_0 accumulated for {len(pheno_cols)} phenotype(s) in "
                f"{sigma0_seconds:.2f}s")

    t_batches = time.perf_counter()
    
    # Assign jackknife blocks
    check_blocks_vs_samples(args.blocks, N_matched, flag="--blocks")
    np.random.seed(args.seed)
    block_ids = np.random.randint(0, args.blocks, N_matched).astype(np.int32)
    n_blocks_nonempty = count_nonempty_blocks(block_ids, args.blocks, flag="--blocks")
    
    # Output setup
    _sanitize_name = sanitize_pheno_name

    out_prefixes = []
    for col in pheno_cols:
        if len(pheno_cols) == 1:
            out_prefixes.append(args.out)
        else:
            out_prefixes.append(f"{args.out}.{_sanitize_name(col)}")

    out_tsv_list = [f"{prefix}.stage1.tsv.gz" for prefix in out_prefixes]
    out_cov_list = [f"{prefix}.cov.gz" for prefix in out_prefixes]
    out_cov_scale_list = [f"{prefix}.cov.scale.gz" for prefix in out_prefixes]
    out_cov_meta_list = [f"{prefix}.cov.meta.json" for prefix in out_prefixes]
    out_cov_ids_list = [f"{prefix}.cov.ids.tsv.gz" for prefix in out_prefixes]
    
    # Column headers
    cols = ['snp_id', 'chr', 'bp', 'effect_allele', 'other_allele']
    for t in taus:
        # Use 3 decimals to avoid tau label collisions on dense grids (e.g. 100 taus).
        cols.extend([f'beta_tau_{t:.3f}', f'se_tau_{t:.3f}'])
    # Appended last so that anything parsing the beta/se block by position is
    # unaffected. This is the genotype's effective sample size, which Stage 2
    # turns into the per-SNP degrees of freedom of the covariance.
    cols.append('d_resid')
    
    # Upper triangle size
    n_cov_elements = T * (T + 1) // 2
    
    # For R_tau computation
    all_betas_list = [([] if args.output_rtau else None) for _ in pheno_cols]
    
    # Process in batches
    total_snps = 0
    skip_counter = _SkipCount()
    pipeline = _BatchPipeline(bool(args.pipeline), N_matched, args.batch_size)
    G_batch = pipeline.batch_buffer
    batch_info = []
    
    from contextlib import ExitStack

    total_compute = 0.0
    total_write = 0.0
    active_chr = None

    write_cov_payload = args.cov_mode not in COV_MODES_NO_PAYLOAD

    with ExitStack() as stack:
        f_tsv_list = [stack.enter_context(gzip.open(path, 'wt')) for path in out_tsv_list]
        # Under --cov-mode pooled the per-SNP payload is not produced at all, so
        # its files are never created rather than created empty: an empty
        # .cov.gz would look to Stage 2 like a truncated jackknife payload.
        f_cov_list = (
            [stack.enter_context(gzip.open(path, 'wb')) for path in out_cov_list]
            if write_cov_payload else [None for _ in out_cov_list]
        )
        f_cov_scale_list = (
            [stack.enter_context(gzip.open(path, 'wb')) for path in out_cov_scale_list]
            if (write_cov_payload and args.cov_dtype == 'int8')
            else [None for _ in out_cov_scale_list]
        )
        f_cov_ids_list = (
            [stack.enter_context(gzip.open(path, 'wt')) for path in out_cov_ids_list]
            if write_cov_payload else [None for _ in out_cov_ids_list]
        )
        for f_tsv in f_tsv_list:
            f_tsv.write('\t'.join(cols) + '\n')
        for f_ids in f_cov_ids_list:
            if f_ids is not None:
                f_ids.write('snp_id\tchr\tbp\teffect_allele\tother_allele\n')
        
        idx_in_batch = 0

        decoded = _decode_variants(variants, args.bgen, bgen_indices, skip_counter)

        for info, variant_chr, dosages in decoded:
            if loco_predictions is not None and variant_chr != active_chr:
                if idx_in_batch > 0:
                    compute_s, write_s = pipeline.submit(
                        _process_batch,
                        G_batch[:, :idx_in_batch], RIF_resid_list, Q, block_ids, batch_info,
                        args, use_multi, f_tsv_list, f_cov_list, f_cov_ids_list,
                        f_cov_scale_list, all_betas_list, n_cov_elements, T,
                        sigma0_sd_list
                    )
                    total_compute += compute_s
                    total_write += write_s
                    total_snps += idx_in_batch
                    idx_in_batch = 0
                    batch_info = []
                    G_batch = pipeline.batch_buffer

                # The RIF/Q the worker reads are about to be replaced, so the
                # in-flight batch has to finish before they are.
                compute_s, write_s = pipeline.drain()
                total_compute += compute_s
                total_write += write_s

                if args.loco_mode == LOCO_MODE_COVARIATE:
                    log_mem(f"Adding the LOCO prediction for chromosome {variant_chr} "
                            f"to the residualising basis")
                else:
                    log_mem(f"Applying LOCO phenotype adjustment for chromosome {variant_chr}")
                RIF_resid_list, Q, q_tau_list = compute_chr_rif_residuals(
                    pheno_cols, y_values_full, X_full, taus, loco_predictions, variant_chr, valid_mask,
                    loco_scales, density_method=args.density_method,
                    loco_mode=args.loco_mode
                )
                if args.sigma0:
                    t_s0 = time.perf_counter()
                    for p_idx, rif in enumerate(RIF_resid_list):
                        sigma0_by_chr[p_idx][variant_chr] = compute_sigma0(rif)
                    sigma0_sd_list = [
                        np.sqrt(np.diag(sigma0_by_chr[p_idx][variant_chr]))
                        for p_idx in range(len(pheno_cols))
                    ]
                    sigma0_seconds += time.perf_counter() - t_s0
                active_chr = variant_chr

            G_batch[:, idx_in_batch] = dosages
            batch_info.append(info)

            idx_in_batch += 1
            
            # Process batch
            if idx_in_batch >= args.batch_size:
                compute_s, write_s = pipeline.submit(
                    _process_batch,
                    G_batch, RIF_resid_list, Q, block_ids, batch_info,
                    args, use_multi, f_tsv_list, f_cov_list, f_cov_ids_list,
                    f_cov_scale_list, all_betas_list, n_cov_elements, T,
                    sigma0_sd_list
                )
                total_compute += compute_s
                total_write += write_s

                total_snps += args.batch_size
                # Running throughput over the whole batch loop. The per-batch
                # figure this replaces timed one call to _process_batch, which
                # under --pipeline returns as soon as the previous batch is done
                # and so reported the submit latency rather than the rate.
                elapsed = time.perf_counter() - t_batches
                log_mem(f"Processed {total_snps} SNPs ({total_snps/elapsed:.1f} SNPs/s)")

                idx_in_batch = 0
                batch_info = []
                G_batch = pipeline.batch_buffer

        # Final partial batch
        if idx_in_batch > 0:
            compute_s, write_s = pipeline.submit(
                _process_batch,
                G_batch[:, :idx_in_batch], RIF_resid_list, Q, block_ids, batch_info,
                args, use_multi, f_tsv_list, f_cov_list, f_cov_ids_list,
                f_cov_scale_list, all_betas_list, n_cov_elements, T,
                sigma0_sd_list
            )
            total_compute += compute_s
            total_write += write_s

            total_snps += idx_in_batch

        # Nothing may be in flight when the output files are closed.
        compute_s, write_s = pipeline.drain()
        total_compute += compute_s
        total_write += write_s
        pipeline.close()

    # Output R_tau if requested
    if args.output_rtau:
        for prefix, all_betas in zip(out_prefixes, all_betas_list):
            if not all_betas:
                continue
            log_mem(f"Computing R_tau correlation matrix for {prefix}...")
            beta_mat = np.array(all_betas)
            R_tau = np.corrcoef(beta_mat.T)

            rtau_file = f"{prefix}.Rtau.tsv"
            with open(rtau_file, 'w') as f:
                f.write('\t'.join([f'tau_{t:.3f}' for t in taus]) + '\n')
                for row in R_tau:
                    f.write('\t'.join([f'{v:.6f}' for v in row]) + '\n')
            logger.info(f"Wrote R_tau: {rtau_file}")
    
    elapsed = time.time() - t_start
    # Record the jackknife block count (and hence the degrees of freedom of the
    # estimated tau-covariance) so Stage 2 can read theta/SE against a t_nu
    # reference instead of a normal one. See fungwas_stage1/cov_meta.py.
    #
    # The pooled covariance sidecar, and the provenance block that lets Stage 2
    # decide whether it may be used. Both are written per phenotype: Sigma_0 is
    # the covariance of one trait's influence functions and means nothing for
    # another's.
    sigma0_meta_list = []
    for p_idx, col in enumerate(pheno_cols):
        chroms = sorted(sigma0_by_chr[p_idx].keys())
        sigma0_meta_list.append(build_sigma0_meta(
            phenotype=col,
            covariates=covar_cols_used,
            taus=taus,
            n_samples=N_matched,
            chromosomes=chroms,
            loco_predictor=loco_files.get(col),
            loco_scale_mode=(None if loco_scales is None else args.loco_scale),
            loco_scale_factor=(None if loco_scales is None
                               else loco_scales.get(col)),
            loco_mode=(LOCO_MODE_COVARIATE
                       if (loco_predictions is not None
                           and args.loco_mode == LOCO_MODE_COVARIATE)
                       else None),
            covariate_basis_rank=int(X_full.shape[1] + (
                1 if (loco_predictions is not None
                      and args.loco_mode == LOCO_MODE_COVARIATE) else 0)),
            density_method=args.density_method,
        ))
    for p_idx, prefix in enumerate(out_prefixes):
        if not args.sigma0:
            continue
        if not sigma0_by_chr[p_idx]:
            logger.warning(
                "No Sigma_0 was accumulated for %s (no variants were scored); "
                "the pooled covariance sidecar is not written.", pheno_cols[p_idx]
            )
            continue
        m_path, s_path = write_sigma0(prefix, sigma0_by_chr[p_idx], sigma0_meta_list[p_idx])
        logger.info("Wrote pooled covariance: %s, %s", m_path, s_path)
    log_mem(f"Sigma_0 compute time: {sigma0_seconds:.2f}s")
    for p_idx, meta_path in enumerate(out_cov_meta_list):
        cov_meta = build_cov_meta(
            cov_dtype=args.cov_dtype,
            n_taus=T,
            n_blocks=args.blocks,
            seed=args.seed,
            n_blocks_nonempty=n_blocks_nonempty,
            analysis=sigma0_analysis_block(sigma0_meta_list[p_idx]),
            se_source=("pooled_sigma0" if args.cov_mode == COV_MODE_POOLED_ONLY
                       else None),
        )
        write_cov_meta(meta_path, cov_meta)
    log_mem(
        f"Batch compute time: {total_compute:.1f}s, write time: {total_write:.1f}s, "
        f"batch total: {time.perf_counter() - t_batches:.1f}s"
    )
    if args.cov_mode == COV_MODE_POOLED_ONLY:
        logger.info(
            "Reminder: --cov-mode pooled-only was in force, so se_tau_* in the TSV "
            "are the POOLED standard errors sqrt(Sigma_0[t,t]/d_resid), not "
            "jackknife estimates, and no per-SNP jackknife covariance exists for "
            "this scan. Stage 2 must be run with cov_source = \"pooled\"."
        )
    else:
        logger.info(
        "Reminder: se_tau_* in the TSV are jackknife estimates. A per-tau Wald test is "
        "t with FAR fewer degrees of freedom than --blocks-1 at extreme taus "
        "(~112 at tau=0.02 vs ~835 at tau=0.50, N=50k, B=1000; 38 at tau=0.02 with "
        "MAF 0.01). See docs/KNOWN_ISSUES.md 2.2 before studentising them."
        )
    log_mem(f"Dropped {skip_counter.non_biallelic} non-biallelic variants")
    log_mem(f"Done! {total_snps} SNPs in {elapsed:.1f}s ({total_snps/elapsed:.1f} SNPs/s)")
    logger.info(f"Output: {', '.join(out_tsv_list)}")


if __name__ == "__main__":
    main()
