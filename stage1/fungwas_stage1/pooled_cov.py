"""Pooled cross-quantile covariance (``Sigma_0``) written alongside Stage 1 output.

WHAT THIS REPLACES.  Stage 1 has always estimated the T x T cross-quantile
covariance of ``beta_tau`` **per variant**, by a delete-one-block jackknife.
That costs 229-289 gzipped bytes per variant (~275 MB per trait genome-wide),
depends on a block count and a random seed, and -- measured against an exact
per-individual HC0 sandwich on real chr22 data -- carries a 5th-to-95th
percentile spread of 0.84 to 1.16 on the standard error of a single SNP.

It does not need to be per variant.  Write the per-tau slope as

    beta_t = sum_i g_resid,i psi_r,i,t / D ,        D = sum_i g_resid,i^2

where ``psi_r`` is the covariate-residualised recentered influence function.
Under working independence (``g_i^2`` uncorrelated with ``psi_i psi_i'`` across
individuals, which is what an i.i.d. sample gives) the sampling covariance of
``beta`` is

    Sigma_v = Sigma_0 / D ,      Sigma_0 = (1/n) Psi_r' Psi_r

and ``Sigma_0`` contains **no genotype at all**.  It is a property of the
phenotype, the covariate set and the tau grid.  The only per-variant quantity
is the scalar ``D``, which Stage 1 already emits as the ``d_resid`` column of
the TSV (added in v0.3.0 for the jackknife degrees-of-freedom law).

MEASURED, on the published chr22 Stage 1 inputs (N = 461,320) against an exact
per-individual HC0 sandwich, 500 variants per trait:

    trait     diag ratio (med / p5 / p95)      max |corr - corr(HC0)| (median)
    bmi       1.00000 / 0.99045 / 1.00981      0.0053
    weight    1.00000 / 0.98993 / 1.01055      0.0053

against the jackknife's 0.99 / 0.725 / 1.307 and 0.30-0.35.  Downstream Stage 2
standard errors agree with the exact sandwich to within 0.03%.

WHAT IS STORED.  One T x T float64 matrix **per trait per chromosome**.  Per
chromosome because the LOCO path subtracts a per-chromosome predictor from the
phenotype before the RIF is built, so the influence functions differ by
chromosome; without ``--loco-preds`` there is one matrix under the key
``"ALL"``.  chr22 at T = 45: 8,280 bytes for the upper triangle against the
4.0 MB of int8 jackknife covariance for the same 17,550 variants.

WHY THE METADATA IS NOT OPTIONAL.  ``Sigma_0`` built with a different covariate
set, a different tau grid or a different sample than the scan it is applied to
is wrong in a way nothing downstream can detect: the standard errors stay
plausible and the point estimates are untouched.  The failure mode that matters
is precisely the one measured for the per-bin genotype-sum estimator -- dropping
the covariates inflates the variance by ``1/(1 - R2_t)``, which is 1.7% for
``bmi`` and 25.6% for ``weight``.  So the sidecar records the covariate set, the
LOCO predictor, the tau grid, the sample count and a fingerprint over all of
them, and Stage 2 refuses to use a ``Sigma_0`` whose fingerprint does not match
the covariance metadata of the scan it is analysing.  This is the same
precedent as the refusal to emit a p-value when the degrees of freedom are
unknown (``cov_meta.resolve_nu``): silence is what let the normal-reference
defect survive.

ON-DISK LAYOUT.  Two files next to the Stage 1 outputs:

``<prefix>.sigma0.meta.json``
    Provenance and the compatibility fingerprint.  JSON, readable by Stage 2
    with or without ``jsonlite``.

``<prefix>.sigma0.tsv.gz``
    The matrices themselves, one row per (chromosome, tau), ``%.17g`` so the
    float64 round-trips exactly.  A text matrix rather than a binary blob
    because it is 8 KB per chromosome, it is read by R and by Python, and it
    is the object a human most wants to be able to look at.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

#: Bumped whenever the on-disk sidecar schema gains or changes a field.
SIGMA0_META_VERSION = 1

#: ``estimator`` value recorded in the sidecar.
SIGMA0_ESTIMATOR = "pooled_influence_covariance"

#: Chromosome key used when the phenotype carries no per-chromosome adjustment,
#: i.e. when Stage 1 ran without ``--loco-preds``.  One matrix then applies to
#: every variant in the run.
SIGMA0_ALL_CHROM = "ALL"

#: Stage 1 TSV column carrying the per-variant scale ``D = sum_i g_resid_i^2``.
SIGMA0_SCALE_COLUMN = "d_resid"

#: The identity Stage 2 applies, recorded on disk so a reader need not guess.
SIGMA0_FORM = "Sigma_variant = Sigma_0 / d_resid"


def compute_sigma0(rif_resid: np.ndarray) -> np.ndarray:
    """Pooled cross-quantile covariance of the residualised influence functions.

    ``Sigma_0 = (1/n) Psi_r' Psi_r`` for ``Psi_r`` the (n, T) covariate-residualised
    RIF matrix -- exactly the matrix Stage 1 already hands to the block-score
    kernel, so this is a by-product of work already done and needs no second
    pass over the phenotype.

    Symmetric and positive semi-definite by construction: it is a Gram matrix.
    That is one of the reasons it replaces the jackknife, which needs a PSD
    repair and an eigenvalue sign rule.

    Parameters
    ----------
    rif_resid : array (n, T)
        Residualised RIF matrix for one phenotype on one chromosome.

    Returns
    -------
    array (T, T), float64
    """
    Psi = np.asarray(rif_resid, dtype=np.float64)
    if Psi.ndim != 2:
        raise ValueError(f"rif_resid must be 2-D (n, T), got shape {Psi.shape}")
    n = Psi.shape[0]
    if n < 2:
        raise ValueError(f"rif_resid needs at least 2 rows to form a covariance, got {n}")
    Sigma0 = (Psi.T @ Psi) / float(n)
    # Force exact symmetry. The gemm is symmetric in exact arithmetic but not
    # bit-for-bit, and Stage 2 reads the upper triangle only.
    return 0.5 * (Sigma0 + Sigma0.T)


def _canonical_taus(taus: Sequence[float]) -> List[float]:
    """Tau grid rounded to the precision the Stage 1 TSV column names carry.

    Stage 1 labels its columns ``beta_tau_%.3f``, so Stage 2 can only ever
    recover the grid to three decimals.  The fingerprint has to agree with what
    Stage 2 can actually check, or it would refuse every real run.
    """
    return [round(float(t), 3) for t in taus]


def analysis_fingerprint(phenotype: str,
                         covariates: Sequence[str],
                         taus: Sequence[float],
                         n_samples: int,
                         loco_predictor: Optional[str],
                         loco_scale_mode: Optional[str] = None,
                         loco_scale_factor: Optional[float] = None) -> str:
    """Stable hash over everything ``Sigma_0`` depends on but the genotype.

    Two Stage 1 artefacts with the same fingerprint were produced from the same
    phenotype, the same covariate columns in the same order, the same tau grid,
    the same number of analysed samples and the same LOCO predictor.  That is
    the whole compatibility condition for ``Sigma_v = Sigma_0 / d_resid``.

    The LOCO predictor enters by path, not by content: reading a 500 MB
    prediction file to hash it would cost more than the rest of Stage 1, and the
    path is what distinguishes a run with the predictor from one without and one
    trait's predictor from another's.

    The ``--loco-scale`` mode and factor enter too.  Two runs over the same
    phenotype, covariates, tau grid and predictor but different scales analyse
    *different variables* -- on a trait like ``log_bmi`` the effect estimates
    differ by a factor of 2.2 (``docs/LOCO_SCALE.md``) -- so a ``Sigma_0`` from
    one is not valid for the other, and without this they would hash alike.
    They default to ``None`` so that a fingerprint computed for a run with no
    LOCO adjustment is unchanged.
    """
    payload = {
        "phenotype": str(phenotype),
        "covariates": [str(c) for c in covariates],
        "taus": _canonical_taus(taus),
        "n_samples": int(n_samples),
        "loco_predictor": None if loco_predictor is None else str(loco_predictor),
    }
    # Added after schema v1 shipped; omitted entirely when absent so that a run
    # without a LOCO adjustment keeps the fingerprint it had before.
    if loco_scale_mode is not None:
        payload["loco_scale_mode"] = str(loco_scale_mode)
    if loco_scale_factor is not None:
        payload["loco_scale_factor"] = round(float(loco_scale_factor), 12)
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:32]


def build_sigma0_meta(phenotype: str,
                      covariates: Sequence[str],
                      taus: Sequence[float],
                      n_samples: int,
                      chromosomes: Sequence[str],
                      loco_predictor: Optional[str] = None,
                      loco_scale_mode: Optional[str] = None,
                      loco_scale_factor: Optional[float] = None,
                      covariate_basis_rank: Optional[int] = None,
                      package_version: Optional[str] = None,
                      matrix_file: Optional[str] = None) -> dict:
    """Build the ``.sigma0.meta.json`` payload.

    Parameters
    ----------
    phenotype : str
        Phenotype column name as given to ``--pheno-col``.
    covariates : sequence of str
        Covariate columns actually used, in the order they entered ``X``,
        **excluding** the intercept Stage 1 prepends.
    taus : sequence of float
        The tau grid.
    n_samples : int
        Number of samples actually analysed -- after the keep file, the
        phenotype/covariate completeness filter, the LOCO membership filter and
        the BGEN sample match.  This is ``n`` in ``Sigma_0 = Psi_r' Psi_r / n``.
    chromosomes : sequence of str
        Normalised chromosome labels for which a matrix was written, or
        ``["ALL"]`` when the run had no per-chromosome phenotype adjustment.
    loco_predictor : str, optional
        Resolved ``.loco`` path for this phenotype, or ``None``.
    loco_scale_mode : str, optional
        The ``--loco-scale`` mode in force (``regression``, ``sd``, ``none``, or a
        literal number), or ``None`` when no LOCO adjustment was applied.
        Recorded because it is not recoverable from the outputs: a scan made
        under the pre-v0.4.0 unscaled behaviour would otherwise be
        indistinguishable from a corrected one.  See ``docs/LOCO_SCALE.md``.
    loco_scale_factor : float, optional
        The numeric factor actually multiplied into the prediction.
    covariate_basis_rank : int, optional
        ``K``, the number of columns of the orthonormal basis ``Q`` (covariates
        plus intercept).  Recorded for diagnosis only.
    """
    if package_version is None:
        try:
            from . import __version__ as package_version  # type: ignore
        except Exception:  # pragma: no cover - version is metadata only
            package_version = "unknown"
    covariates = [str(c) for c in covariates]
    chromosomes = [str(c) for c in chromosomes]
    return {
        "sigma0_meta_version": SIGMA0_META_VERSION,
        "package_version": str(package_version),
        "estimator": SIGMA0_ESTIMATOR,
        "form": SIGMA0_FORM,
        "scale_column": SIGMA0_SCALE_COLUMN,
        "phenotype": str(phenotype),
        "n_taus": int(len(taus)),
        "taus": _canonical_taus(taus),
        "n_samples": int(n_samples),
        "n_covariates": len(covariates),
        "covariates": covariates,
        "covariate_basis_rank": (None if covariate_basis_rank is None
                                 else int(covariate_basis_rank)),
        "loco_predictor": loco_predictor,
        "loco_scale_mode": None if loco_scale_mode is None else str(loco_scale_mode),
        "loco_scale_factor": (None if loco_scale_factor is None
                              else float(loco_scale_factor)),
        "chromosomes": chromosomes,
        "matrix_file": matrix_file,
        "fingerprint": analysis_fingerprint(
            phenotype, covariates, taus, n_samples, loco_predictor,
            loco_scale_mode, loco_scale_factor
        ),
    }


def sigma0_paths(prefix: str) -> Tuple[str, str]:
    """``(matrix_path, meta_path)`` for a Stage 1 output prefix."""
    return f"{prefix}.sigma0.tsv.gz", f"{prefix}.sigma0.meta.json"


def write_sigma0(prefix: str,
                 matrices: Dict[str, np.ndarray],
                 meta: dict) -> Tuple[str, str]:
    """Write the matrices and their metadata sidecar.

    The matrix file is a gzipped TSV with columns ``chr``, ``tau`` and one
    column per tau, ``%.17g`` so float64 survives the round trip.  Rows are
    ordered by the ``chromosomes`` list in ``meta`` and then by tau index.
    """
    matrix_path, meta_path = sigma0_paths(prefix)
    taus = list(meta["taus"])
    T = int(meta["n_taus"])
    chroms = list(meta["chromosomes"])
    missing = [c for c in chroms if c not in matrices]
    if missing:
        raise ValueError(f"Sigma_0 metadata lists chromosomes with no matrix: {missing}")

    header = ["chr", "tau"] + [f"s_{t:.3f}" for t in taus]
    with gzip.open(matrix_path, "wt") as handle:
        handle.write("\t".join(header) + "\n")
        for chrom in chroms:
            S = np.asarray(matrices[chrom], dtype=np.float64)
            if S.shape != (T, T):
                raise ValueError(
                    f"Sigma_0 for chromosome {chrom} has shape {S.shape}, expected ({T}, {T})"
                )
            if not np.all(np.isfinite(S)):
                raise ValueError(f"Sigma_0 for chromosome {chrom} is not finite")
            for t in range(T):
                row = [str(chrom), f"{taus[t]:.3f}"]
                row.extend(f"{v:.17g}" for v in S[t])
                handle.write("\t".join(row) + "\n")

    meta = dict(meta)
    meta["matrix_file"] = matrix_path.rsplit("/", 1)[-1]
    with open(meta_path, "w") as handle:
        json.dump(meta, handle, indent=2)
    return matrix_path, meta_path


def read_sigma0(prefix: str) -> Tuple[dict, Dict[str, np.ndarray]]:
    """Read back what :func:`write_sigma0` wrote.

    Returns ``(meta, matrices)`` with ``matrices`` keyed by the normalised
    chromosome label.
    """
    matrix_path, meta_path = sigma0_paths(prefix)
    with open(meta_path) as handle:
        meta = json.load(handle)
    T = int(meta["n_taus"])
    matrices: Dict[str, List[np.ndarray]] = {}
    with gzip.open(matrix_path, "rt") as handle:
        handle.readline()
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != T + 2:
                raise ValueError(
                    f"{matrix_path}: expected {T + 2} fields per row, got {len(fields)}"
                )
            matrices.setdefault(fields[0], []).append(
                np.array([float(v) for v in fields[2:]], dtype=np.float64)
            )
    out = {}
    for chrom, rows in matrices.items():
        if len(rows) != T:
            raise ValueError(
                f"{matrix_path}: chromosome {chrom} has {len(rows)} rows, expected {T}"
            )
        out[chrom] = np.vstack(rows)
    return meta, out


def sigma_for_variant(sigma0: np.ndarray, d_resid: float) -> np.ndarray:
    """``Sigma_0 / d_resid`` -- the per-variant covariance of ``beta_tau``."""
    d = float(d_resid)
    if not np.isfinite(d) or d <= 0:
        raise ValueError(f"d_resid must be finite and positive, got {d_resid!r}")
    return np.asarray(sigma0, dtype=np.float64) / d


def sigma0_analysis_block(meta: dict) -> dict:
    """The subset of the sidecar that ``.cov.meta.json`` also carries.

    Stage 2 compares the two blocks field by field and refuses on any
    disagreement, so the two writers must produce the same shape.  Keeping the
    construction in one function is what makes that true by construction rather
    than by inspection.
    """
    return {
        "phenotype": meta["phenotype"],
        "n_samples": meta["n_samples"],
        "n_covariates": meta["n_covariates"],
        "covariates": list(meta["covariates"]),
        "loco_predictor": meta["loco_predictor"],
        "loco_scale_mode": meta.get("loco_scale_mode"),
        "loco_scale_factor": meta.get("loco_scale_factor"),
        "fingerprint": meta["fingerprint"],
    }


def check_fingerprint(sigma0_meta: dict, cov_meta: dict) -> None:
    """Refuse a ``Sigma_0`` that does not belong to the scan being analysed.

    Raises ``ValueError`` naming the first field that differs.  Python-side
    counterpart of the Stage 2 check, used by the Stage 1 tests and by anyone
    building a pooled covariance outside R.
    """
    analysis = cov_meta.get("analysis") if isinstance(cov_meta, dict) else None
    if not isinstance(analysis, dict):
        raise ValueError(
            "The Stage 1 .cov.meta.json carries no 'analysis' block, so there is "
            "nothing to check this Sigma_0 against. It was written by a Stage 1 "
            "older than the pooled covariance. Re-run Stage 1, or use the "
            "jackknife covariance that scan did emit."
        )
    theirs = sigma0_analysis_block(sigma0_meta)
    for key in ("phenotype", "n_samples", "n_covariates", "covariates",
                "loco_predictor", "fingerprint"):
        a, b = theirs.get(key), analysis.get(key)
        if a != b:
            raise ValueError(
                f"Sigma_0 does not match the scan it is being applied to: {key} is "
                f"{a!r} in the Sigma_0 sidecar and {b!r} in the covariance metadata. "
                f"Sigma_v = Sigma_0 / d_resid is only the sampling covariance of "
                f"beta_tau when both were built from the same phenotype, covariate "
                f"set, tau grid, LOCO predictor and sample. Using a mismatched one "
                f"biases every standard error in the same direction with nothing "
                f"downstream able to detect it."
            )


def check_taus(sigma0_meta: dict, taus: Sequence[float]) -> None:
    """Refuse a ``Sigma_0`` built on a different tau grid."""
    theirs = [round(float(t), 3) for t in sigma0_meta.get("taus", [])]
    ours = _canonical_taus(taus)
    if theirs != ours:
        raise ValueError(
            f"Sigma_0 was built on a tau grid of {len(theirs)} levels that does not "
            f"match the {len(ours)} being analysed. Sigma_0 is indexed by tau; a "
            f"different grid makes every element refer to a different quantile."
        )
