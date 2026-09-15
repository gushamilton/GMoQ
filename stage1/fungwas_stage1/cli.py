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
denominator of every per-tau slope, and the genotype's effective sample size.
Stage 2 needs it to work out the effective degrees of freedom of the estimated
covariance per SNP, which is what makes its t reference correct.

WARNING about the per-tau standard errors in this file.  se_tau_* are jackknife
estimates, so beta/se at a single tau is t distributed, NOT normal -- and the
effective degrees of freedom at an extreme tau are far below the block count.
Measured at N = 50,000 with --blocks 1000: about 835 df at tau = 0.50, but 112
at tau = 0.02, and 38 at tau = 0.02 for a MAF 0.01 variant. Do not studentise a
per-tau slope against a normal, against blocks - 1, or against the degrees of
freedom Stage 2 reports for a fitted parameter. See docs/KNOWN_ISSUES.md 2.2.
"""

import argparse
import gzip
import logging
import os
import re
import sqlite3
import struct
import sys
import time
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
    compute_covariance_from_blocks,
    residual_sum_squares,
    resolve_covar_cols,
    DEFAULT_N_BLOCKS,
    resolve_n_blocks,
    count_nonempty_blocks,
    check_blocks_vs_samples,
    set_numpy_fallback_allowed,
    require_cpp_extension,
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
                        X: np.ndarray, taus: np.ndarray) -> tuple[list[np.ndarray], np.ndarray, list[np.ndarray]]:
    """Compute residualized RIF matrices for one or more phenotypes."""
    log_mem(f"Residualizing phenotypes on {X.shape[1]} covariates...")
    Q, _ = np.linalg.qr(X)

    rif_resid_list = []
    q_tau_list = []
    for col in pheno_cols:
        log_mem(f"Computing RIF matrix for {col}...")
        t_rif = time.perf_counter()
        y = np.asarray(y_values[col], dtype=np.float64)
        rif, q_tau, _ = compute_rif(y, taus)
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
    sample_to_idx = {s: i for i, s in enumerate(bgen_samples)}
    matched_indices = [sample_to_idx.get(fid) for fid in sample_ids]
    valid_mask = np.array([idx is not None for idx in matched_indices])
    bgen_indices = np.array([idx for idx in matched_indices if idx is not None])
    return valid_mask, bgen_indices


def extract_dosage(variant, bgen_indices: Optional[np.ndarray] = None):
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
    # alt_dosage is the expected count of allele[1]
    # For diploid, dosage of allele[0] = 2 - alt_dosage
    alt_dosages = variant.alt_dosage
    dosages = 2.0 - alt_dosages
    if bgen_indices is not None:
        return dosages[bgen_indices]
    return dosages


def write_stats_batch(stats_list, batch_info, args, f_tsv_list, f_cov_list, f_cov_ids_list,
                      f_cov_scale_list, all_betas_list, n_cov_elements, n_taus):
    """Write one processed batch to the stage1 outputs."""
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

            cov_upper = Sigma[np.triu_indices(n_taus)]
            if args.cov_dtype == 'int8':
                scale = np.max(np.abs(cov_upper)) / 127.0
                if not np.isfinite(scale) or scale <= 0:
                    scale = 1.0
                cov_q = np.clip(np.rint(cov_upper / scale), -127, 127).astype(np.int8)
                f_cov.write(cov_q.tobytes(order='C'))
                f_cov_scale.write(struct.pack('f', float(scale)))
            else:
                f_cov.write(struct.pack(f'{n_cov_elements}f', *cov_upper))
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

            if args.cov_dtype == 'int8':
                scale = np.max(np.abs(cov_upper)) / 127.0
                if not np.isfinite(scale) or scale <= 0:
                    scale = 1.0
                cov_q = np.clip(np.rint(cov_upper / scale), -127, 127).astype(np.int8)
                f_cov.write(cov_q.tobytes(order='C'))
                f_cov_scale.write(struct.pack('f', float(scale)))
            else:
                f_cov.write(struct.pack(f'{n_cov_elements}f', *cov_upper))
            f_cov_ids.write('\t'.join(map(str, info)) + '\n')

            if args.output_rtau:
                all_betas.append(np.asarray(beta))


def _process_batch(G_view, RIF_resid_list, Q, block_ids, batch_info, args, use_multi,
                   f_tsv_list, f_cov_list, f_cov_ids_list, f_cov_scale_list,
                   all_betas_list, n_cov_elements, T):
    """Run the kernel over one batch and write it out. Returns (compute_s, write_s)."""
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


def compute_chr_rif_residuals(pheno_cols: list[str], y_values_full: dict[str, np.ndarray],
                              X_full: np.ndarray, taus: np.ndarray,
                              loco_predictions: Optional[dict[str, dict[str, np.ndarray]]],
                              chrom: Optional[str], valid_mask: np.ndarray,
                              loco_scales: Optional[dict[str, float]] = None) -> tuple[list[np.ndarray], np.ndarray, list[np.ndarray]]:
    """Build residualized RIF matrices for a chromosome, optionally using LOCO-adjusted phenotypes."""
    if loco_predictions is None:
        y_values_chr = {col: y_values_full[col][valid_mask] for col in pheno_cols}
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

    rif_resid_list, q_chr, q_tau_list = build_rif_residuals(y_values_chr, pheno_cols, X_full[valid_mask], taus)
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
    parser.add_argument('--blocks', type=int, default=None,
                        help=f'Number of jackknife blocks (default: {DEFAULT_N_BLOCKS}). '
                             'Must exceed the number of taus. Fewer blocks stay correctly '
                             'sized under the Stage 2 t reference but give a noisier SE.')
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
    
    log_mem(f"Starting FungWas Stage 1 (taus={T}, blocks={args.blocks})")
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
    if args.panel_kernel and args.panel_size < 1:
        raise ValueError("--panel-size must be >= 1")
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
        pheno_cols, y_values_full, X_full, taus, None, None, valid_mask
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
    if args.sigma0 and loco_predictions is None:
        t_s0 = time.perf_counter()
        for p_idx, rif in enumerate(RIF_resid_list):
            sigma0_by_chr[p_idx][SIGMA0_ALL_CHROM] = compute_sigma0(rif)
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
    skipped_non_biallelic = 0
    G_batch = np.zeros((N_matched, args.batch_size), dtype=np.float64, order='F')
    batch_info = []
    
    from contextlib import ExitStack

    total_compute = 0.0
    total_write = 0.0
    active_chr = None

    with ExitStack() as stack:
        f_tsv_list = [stack.enter_context(gzip.open(path, 'wt')) for path in out_tsv_list]
        f_cov_list = [stack.enter_context(gzip.open(path, 'wb')) for path in out_cov_list]
        f_cov_scale_list = (
            [stack.enter_context(gzip.open(path, 'wb')) for path in out_cov_scale_list]
            if args.cov_dtype == 'int8' else [None for _ in out_cov_scale_list]
        )
        f_cov_ids_list = [stack.enter_context(gzip.open(path, 'wt')) for path in out_cov_ids_list]
        for f_tsv in f_tsv_list:
            f_tsv.write('\t'.join(cols) + '\n')
        for f_ids in f_cov_ids_list:
            f_ids.write('snp_id\tchr\tbp\teffect_allele\tother_allele\n')
        
        idx_in_batch = 0
        
        for variant in variants:
            if len(variant.alleles) != 2:
                skipped_non_biallelic += 1
                continue

            variant_chr = normalize_chr_label(variant.chrom)
            if loco_predictions is not None and variant_chr != active_chr:
                if idx_in_batch > 0:
                    compute_s, write_s = _process_batch(
                        G_batch[:, :idx_in_batch], RIF_resid_list, Q, block_ids, batch_info,
                        args, use_multi, f_tsv_list, f_cov_list, f_cov_ids_list,
                        f_cov_scale_list, all_betas_list, n_cov_elements, T
                    )
                    total_compute += compute_s
                    total_write += write_s
                    total_snps += idx_in_batch
                    idx_in_batch = 0
                    batch_info = []

                log_mem(f"Applying LOCO phenotype adjustment for chromosome {variant_chr}")
                RIF_resid_list, Q, q_tau_list = compute_chr_rif_residuals(
                    pheno_cols, y_values_full, X_full, taus, loco_predictions, variant_chr, valid_mask,
                    loco_scales
                )
                if args.sigma0:
                    t_s0 = time.perf_counter()
                    for p_idx, rif in enumerate(RIF_resid_list):
                        sigma0_by_chr[p_idx][variant_chr] = compute_sigma0(rif)
                    sigma0_seconds += time.perf_counter() - t_s0
                active_chr = variant_chr

            # Extract dosage
            all_dosages = extract_dosage(variant)
            dosages = all_dosages[bgen_indices]
            
            G_batch[:, idx_in_batch] = dosages

            batch_info.append((
                variant.rsid, variant.chrom, variant.pos,
                variant.alleles[0], variant.alleles[1]
            ))
            
            idx_in_batch += 1
            
            # Process batch
            if idx_in_batch >= args.batch_size:
                t_batch = time.time()
                
                compute_s, write_s = _process_batch(
                    G_batch, RIF_resid_list, Q, block_ids, batch_info,
                    args, use_multi, f_tsv_list, f_cov_list, f_cov_ids_list,
                    f_cov_scale_list, all_betas_list, n_cov_elements, T
                )
                total_compute += compute_s
                total_write += write_s
                
                total_snps += args.batch_size
                elapsed = time.time() - t_batch
                log_mem(f"Processed {total_snps} SNPs ({args.batch_size/elapsed:.1f} SNPs/s)")
                
                idx_in_batch = 0
                batch_info = []
        
        # Final partial batch
        if idx_in_batch > 0:
            compute_s, write_s = _process_batch(
                G_batch[:, :idx_in_batch], RIF_resid_list, Q, block_ids, batch_info,
                args, use_multi, f_tsv_list, f_cov_list, f_cov_ids_list,
                f_cov_scale_list, all_betas_list, n_cov_elements, T
            )
            total_compute += compute_s
            total_write += write_s
            
            total_snps += idx_in_batch
    
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
            covariate_basis_rank=int(X_full.shape[1]),
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
        )
        write_cov_meta(meta_path, cov_meta)
    log_mem(
        f"Batch compute time: {total_compute:.1f}s, write time: {total_write:.1f}s, "
        f"batch total: {time.perf_counter() - t_batches:.1f}s"
    )
    logger.info(
        "Reminder: se_tau_* in the TSV are jackknife estimates. A per-tau Wald test is "
        "t with FAR fewer degrees of freedom than --blocks-1 at extreme taus "
        "(~112 at tau=0.02 vs ~835 at tau=0.50, N=50k, B=1000; 38 at tau=0.02 with "
        "MAF 0.01). See docs/KNOWN_ISSUES.md 2.2 before studentising them."
    )
    log_mem(f"Dropped {skipped_non_biallelic} non-biallelic variants")
    log_mem(f"Done! {total_snps} SNPs in {elapsed:.1f}s ({total_snps/elapsed:.1f} SNPs/s)")
    logger.info(f"Output: {', '.join(out_tsv_list)}")


if __name__ == "__main__":
    main()
