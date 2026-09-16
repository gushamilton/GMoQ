"""
Core Stage 1 computation functions.

Provides:
- RIF computation (Recentered Influence Functions)
- Block score accumulation for jackknife
- Covariance estimation from block statistics
"""

import numpy as np
from typing import Tuple, Optional
import logging

logger = logging.getLogger(__name__)

#: Default number of delete-one-block jackknife blocks.
#:
#: Under a fixed-RIF, regular-block approximation, increasing B reduces
#: jackknife SE noise; ``sd(SE_hat)/SE = 1/sqrt(2*(B-1))`` is a
#: design-variance approximation, not a guarantee for the full plug-in RIF
#: estimator. Stage 2's ``t_{B-1}`` reference is also an approximation whose
#: size depends on the covariance source, effective-df model, and estimand; it
#: is not guaranteed correctly sized for every B. Measured cost of this default
#: relative to the old B = 65: about +30% Stage 1 wall time and peak memory
#: around 1 GB.
DEFAULT_N_BLOCKS = 1000

# Density estimators used to construct the Stage 1 RIF.  The Gaussian path is
# the historical/default estimator.  The fourth-order Gaussian kernel is an
# explicit optional method: it has zero second moment and therefore removes
# the leading smooth KDE bias, but its negative lobes mean that a density floor
# is still required at the RIF denominator.  Keeping the choice explicit is
# important because it changes the estimand's nuisance construction.
DENSITY_METHOD_SILVERMAN_GAUSSIAN = "silverman_gaussian"
DENSITY_METHOD_FOURTH_ORDER_GAUSSIAN = "fourth_order_gaussian"
DENSITY_METHODS = (
    DENSITY_METHOD_SILVERMAN_GAUSSIAN,
    DENSITY_METHOD_FOURTH_ORDER_GAUSSIAN,
)

ALLOW_NUMPY_FALLBACK = False
NUMPY_FALLBACK_WARNING = (
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    "WARNING: NUMPY FALLBACK ENABLED FOR STAGE 1\n"
    "WARNING: THIS PATH IS SUPER SUPER SLOW (~1000X SLOWER THAN THE C++ PATH)\n"
    "WARNING: DO NOT USE THIS FOR PRODUCTION GWAS UNLESS YOU ARE DEBUGGING\n"
    "WARNING: IT IS ONLY RUNNING BECAUSE --allow-numpy-fallback WAS EXPLICITLY SET\n"
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
)

# Try to import C++ extension
try:
    from . import _stage1_cpp as cpp_ext
    HAVE_CPP = True
    logger.info("C++ extension loaded successfully")
except ImportError:
    HAVE_CPP = False
    logger.warning("C++ extension not found; Stage 1 will fail unless NumPy fallback is explicitly enabled.")


def set_numpy_fallback_allowed(allow: bool):
    """Control whether the very slow NumPy fallback may be used."""
    global ALLOW_NUMPY_FALLBACK
    ALLOW_NUMPY_FALLBACK = allow
    if allow and not HAVE_CPP:
        logger.warning(NUMPY_FALLBACK_WARNING)


def resolve_n_blocks(n_blocks: Optional[int], n_taus: int, flag: str = "--blocks") -> int:
    """Resolve and validate the jackknife block count.

    ``None`` means "use the default". The hard constraint is ``B > T``: the
    jackknife covariance has rank ``B - 1``, so fewer blocks than taus makes it
    singular. For a fixed RIF and regular independent blocks, the covariance
    interpretation is approximate; fewer blocks make the denominator noisier.
    The Stage 2 ``t_{B-1}`` reference is conditional on its effective-df
    assumptions and does not provide universal size guarantees, especially
    when q/f/RIF nuisance estimates are held outside the block sums.
    """
    if n_blocks is None:
        n_blocks = DEFAULT_N_BLOCKS
        logger.info("Using default %s=%d jackknife blocks", flag, n_blocks)
    n_blocks = int(n_blocks)

    if n_blocks <= n_taus:
        raise ValueError(
            f"{flag} ({n_blocks}) must be greater than the number of taus ({n_taus}): "
            f"the jackknife covariance has rank {flag}-1 and would be singular. "
            f"Increase {flag} or reduce the tau grid."
        )
    if n_blocks < DEFAULT_N_BLOCKS:
        logger.warning(
            "%s=%d is below the default of %d. The fixed-RIF covariance and "
            "Stage 2's t_%d reference remain approximations; SE_hat "
            "carries %.1f%% relative noise (1/sqrt(2*(B-1))) against %.1f%% at the "
            "default, which costs power.",
            flag, n_blocks, DEFAULT_N_BLOCKS, n_blocks - 1,
            100.0 / np.sqrt(2.0 * (n_blocks - 1)),
            100.0 / np.sqrt(2.0 * (DEFAULT_N_BLOCKS - 1)),
        )
    return n_blocks


#: Tolerance for the covariate-basis orthonormality check performed at every
#: kernel entry point (see :func:`check_orthonormal_basis`).
#:
#: Chosen to sit far above float64 QR noise and far below the smallest instance
#: of the defect it exists to catch. Both bounds are measured, on the real
#: UK Biobank BMI covariate matrix (486,254 x 13 = intercept + age + sex + 10
#: ancestry PCs, of which 461,320 rows are in the chr22 genotype release):
#:
#: * **Floor - a correctly built basis.** LAPACK's Householder QR is
#:   unconditionally backward stable in the orthogonality sense: the deviation
#:   is O(K**1.5 * eps) irrespective of cond(X) - including for a rank-deficient
#:   X, whose Q still has orthonormal columns - and forming the Gram matrix adds
#:   O(sqrt(N) * eps). Measured: 2.3e-14 over all 486,254 rows, 8.4e-15 over the
#:   461,320 matched rows. 1e-8 is 5.6 orders of magnitude above that, so it
#:   cannot fire on a basis that was built on the rows it is used on.
#: * **Ceiling - the smallest realisable instance of the defect.** Dropping even
#:   a single row from the sample the basis was built on perturbs Q'Q by
#:   q_i q_i'. Measured: 2.3e-6 for an arbitrary row and 1.3e-4 for the most
#:   influential one - both at least two orders of magnitude ABOVE 1e-8. The
#:   production defect (the full 24,934-row mismatch) measured 9.7e-1, eight
#:   orders above.
#:
#: The tolerance therefore leaves no gap: no numerically sound basis fails it,
#: and no mis-built one passes it. That is why this raises rather than warns -
#: see :func:`check_orthonormal_basis`.
#:
#: Cost: one N-by-K Gram product per kernel call, measured at 14 ms for
#: (461320, 13) against a ~3.5 s production batch of 500 SNPs, i.e. ~0.4%.
ORTHONORMAL_BASIS_TOL = 1e-8


def check_orthonormal_basis(Q: np.ndarray, tol: float = ORTHONORMAL_BASIS_TOL,
                            context: str = "Stage 1 kernel") -> float:
    """Fail loudly unless ``Q`` is orthonormal over the rows it will be applied to.

    Every Stage 1 kernel residualises the genotype as ``g - Q (Q' g)``. That is
    the orthogonal projection onto the complement of ``col(X)`` **only if**
    ``Q'Q = I``. If it is not, the covariate-parallel part of ``g`` is only
    partly removed and what survives leaks into the numerator of every slope: a
    deterministic bias in ``beta_tau``, with standard errors that are left
    almost untouched and therefore give no warning that anything is wrong.

    The way this happens in practice is a basis built on one sample and applied
    to a subset of it - a QR computed on the keep file, then row-subset to the
    samples that matched the genotype file. Orthonormality over N rows says
    nothing about orthonormality over n < N of them. On the UK Biobank BMI
    inputs (486,254 keep rows, 461,320 matched) the row-subset basis reached
    ``max |Q'Q - I| = 9.7e-1``, with the ancestry-PC columns retaining 3-9% of
    their norm, and shifted ``beta_tau`` by a median of 0.12 SE and up to
    2.3 SE - in a consistent direction at all 45 taus.

    Raising is deliberate. The failure is silent by construction, so a warning
    would be swallowed by the thousands of lines a Stage 1 run already logs; no
    legitimate caller passes a non-orthonormal basis (the kernels all document
    ``Q`` as the QR factor); and the tolerance sits 5.6 orders of magnitude above
    measured QR noise, so a spurious failure is not a realistic risk.

    Parameters
    ----------
    Q : array (N, K)
        Covariate basis, as handed to the kernel.
    tol : float
        Maximum tolerated ``max |Q'Q - I|``. See ``ORTHONORMAL_BASIS_TOL``.
    context : str
        Included in the error message to identify the caller.

    Returns
    -------
    float
        The measured ``max |Q'Q - I|``.

    Raises
    ------
    ValueError
        If the deviation exceeds ``tol``, or is not finite.
    """
    Q = np.asarray(Q, dtype=np.float64)
    if Q.ndim != 2:
        raise ValueError(f"Covariate basis Q must be 2-D, got shape {Q.shape} ({context})")
    gram = Q.T @ Q
    deviation = float(np.max(np.abs(gram - np.eye(Q.shape[1]))))
    if not np.isfinite(deviation) or deviation > tol:
        raise ValueError(
            f"Covariate basis is not orthonormal over the {Q.shape[0]} rows being "
            f"analysed: max |Q'Q - I| = {deviation:.6g} > {tol:g} ({context}). "
            f"The kernels residualise as g - Q (Q' g), which is a projection only if "
            f"Q'Q = I; otherwise beta_tau is biased while se_tau is not, so nothing "
            f"downstream can detect it. The usual cause is a basis built on a larger "
            f"sample and then row-subset (e.g. QR on the keep file, subset to the "
            f"genotype-matched samples). Build it as np.linalg.qr(X[rows]) on exactly "
            f"the rows being analysed."
        )
    return deviation


def residual_sum_squares(G: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Per-SNP ``D = sum_i g_resid_i^2`` for genotypes residualised on ``Q``.

    ``D`` is the denominator of every per-tau slope, and it is also the
    genotype's effective sample size: ``N_eff = D^2 / sum_i g_resid_i^4``
    collapses to ``D`` because a genotype's fourth central moment equals its
    second (measured ``N_eff/D`` = 1.003 on real imputed dosages). Stage 2
    needs it to compute the effective degrees of freedom of the jackknife
    covariance per SNP.

    Used only on the ``--fused-finalize`` path, where the per-block statistics
    that carry ``D`` are never materialised. Everywhere else ``D`` is read off
    the block statistics directly, which is exactly the quantity the jackknife
    divides by.

    Parameters
    ----------
    G : array (N, M)
        Genotype dosages; NaN is mean-imputed per SNP, as in the kernels.
    Q : array (N, K)
        Orthonormal covariate basis.
    """
    check_orthonormal_basis(Q, context="residual_sum_squares")
    G = np.asarray(G, dtype=np.float64)
    col_means = np.nanmean(G, axis=0)
    Gi = np.where(np.isnan(G), col_means, G)
    R = Gi - Q @ (Q.T @ Gi)
    return np.einsum("ij,ij->j", R, R)


def count_nonempty_blocks(block_ids: np.ndarray, n_blocks: int, flag: str = "--blocks") -> int:
    """Number of jackknife blocks that actually received a sample.

    Deleting an empty block reproduces the full-sample estimate, so empty
    blocks contribute no replicate and the degrees of freedom are set by the
    blocks that are populated, not by the requested count. At the default of
    1000 blocks this only bites on small datasets; it could not happen at the
    old default of 65.
    """
    n_nonempty = int(np.unique(np.asarray(block_ids)).size)
    if n_nonempty < n_blocks:
        logger.warning(
            "%d of %d jackknife blocks are empty (%d samples for %s=%d). "
            "Degrees of freedom are set by the %d populated blocks. Lower %s.",
            n_blocks - n_nonempty, n_blocks, len(block_ids), flag, n_blocks,
            n_nonempty, flag,
        )
    return n_nonempty


def check_blocks_vs_samples(n_blocks: int, n_samples: int, flag: str = "--blocks") -> None:
    """Refuse a block count that cannot be populated by the sample."""
    if n_blocks > n_samples:
        raise ValueError(
            f"{flag} ({n_blocks}) exceeds the number of analysed samples ({n_samples}). "
            f"Most blocks would be empty and the recorded degrees of freedom would be "
            f"wrong. Pass an explicit {flag} no larger than {n_samples}."
        )


def require_cpp_extension(context: str = "Stage 1 computation"):
    """Fail fast unless the compiled extension is present or fallback was explicitly allowed."""
    if HAVE_CPP or ALLOW_NUMPY_FALLBACK:
        return
    raise RuntimeError(
        f"{context} requires the compiled C++ extension, but it is not available. "
        "Build it with `python setup.py build_ext --inplace` in `fungwas/stage1`, or rerun "
        "the CLI with `--allow-numpy-fallback` if you intentionally want the ~1000x slower path."
    )


def resolve_covar_cols(requested: Optional[list[str]], available: list[str]) -> list[str]:
    """
    Resolve requested covariate columns against available columns.

    Parameters
    ----------
    requested : list of str or None
        Covariate columns requested by the user.
    available : list of str
        Available covariate columns in the file (excluding FID/IID).

    Returns
    -------
    list of str
        Columns to use as covariates.
    """
    if requested is None:
        return available
    missing = [c for c in requested if c not in available]
    if missing:
        raise ValueError(f"Covariate column(s) not found: {missing}. Available: {available}")
    return requested


def compute_rif(Y: np.ndarray, taus: np.ndarray,
                density_floor: float = 1e-8,
                density_method: str = DENSITY_METHOD_SILVERMAN_GAUSSIAN
                ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute Recentered Influence Function matrix.
    
    Parameters
    ----------
    Y : array (N,)
        Phenotype values
    taus : array (T,)
        Quantile levels
    density_floor : float
        Minimum density value to avoid division issues
    density_method : {"silverman_gaussian", "fourth_order_gaussian"}
        Kernel used for the density at each empirical quantile.  The default
        is the historical Silverman Gaussian KDE.  ``fourth_order_gaussian``
        uses ``(1.5 - 0.5*u**2) * phi(u)`` with the same Silverman bandwidth.
        
    Returns
    -------
    RIF : array (N, T)
        RIF values for each sample and tau
    q_tau : array (T,)
        Empirical quantiles at each tau
    fhat : array (T,)
        Kernel density estimates at each quantile
    """
    Y = np.asarray(Y, dtype=np.float64)
    taus = np.asarray(taus, dtype=np.float64)
    if density_method not in DENSITY_METHODS:
        raise ValueError(
            f"density_method must be one of {DENSITY_METHODS}, got {density_method!r}"
        )
    if not np.isfinite(density_floor) or density_floor <= 0:
        raise ValueError(f"density_floor must be finite and positive, got {density_floor!r}")
    if Y.ndim != 1 or Y.size < 2 or not np.all(np.isfinite(Y)):
        raise ValueError("Y must be a finite one-dimensional array with at least two values")
    if taus.ndim != 1 or taus.size == 0 or np.any(~np.isfinite(taus)) or np.any((taus <= 0) | (taus >= 1)):
        raise ValueError("taus must be a non-empty one-dimensional array strictly inside (0, 1)")
    N = len(Y)
    T = len(taus)
    
    # Empirical quantiles (type 8 for unbiased)
    q_tau = np.quantile(Y, taus, method='median_unbiased')
    
    # Kernel density at each quantile
    # Silverman's rule for bandwidth
    h = 1.06 * np.std(Y) * N ** (-0.2)
    
    fhat = np.zeros(T)
    normalizer = np.sqrt(2 * np.pi)
    for t, q in enumerate(q_tau):
        u = (Y - q) / h
        kernel = np.exp(-0.5 * u ** 2)
        if density_method == DENSITY_METHOD_FOURTH_ORDER_GAUSSIAN:
            kernel *= 1.5 - 0.5 * u ** 2
        fhat[t] = np.mean(kernel) / (h * normalizer)
    fhat = np.maximum(fhat, density_floor)
    
    # RIF = q + (tau - I(Y <= q)) / f(q)
    RIF = np.zeros((N, T))
    for t in range(T):
        indicator = (Y <= q_tau[t]).astype(float)
        RIF[:, t] = q_tau[t] + (taus[t] - indicator) / fhat[t]
    
    return RIF, q_tau, fhat


def residualize(Y: np.ndarray, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Residualize Y on covariates X using QR decomposition.
    
    Parameters
    ----------
    Y : array (N,) or (N, T)
        Outcome(s) to residualize
    X : array (N, K)
        Covariate matrix (should include intercept)
        
    Returns
    -------
    Y_resid : array, same shape as Y
        Residualized outcome
    Q : array (N, K)
        Orthonormal basis for X
    """
    Q, R = np.linalg.qr(X)
    if Y.ndim == 1:
        Y_resid = Y - Q @ (Q.T @ Y)
    else:
        Y_resid = Y - Q @ (Q.T @ Y)
    return Y_resid, Q


def compute_block_scores(G: np.ndarray, Y_resid: np.ndarray, Q: np.ndarray,
                         block_ids: np.ndarray, n_blocks: int,
                         min_mac: float = 5.0, min_var: float = 1e-8,
                         n_threads: int = 1) -> np.ndarray:
    """
    Compute per-block sufficient statistics for jackknife.
    
    Parameters
    ----------
    G : array (N, M)
        Genotype dosages (may contain NaN for missing)
    Y_resid : array (N, T)
        Residualized RIF matrix
    Q : array (N, K)
        QR basis for covariates
    block_ids : array (N,)
        Block assignment for each sample (0 to n_blocks-1)
    n_blocks : int
        Number of jackknife blocks
    min_mac : float
        Minimum minor allele count filter
    min_var : float
        Minimum variance filter
    n_threads : int
        Number of OpenMP threads (if C++ available)
        
    Returns
    -------
    stats : array (M, n_blocks, 1+T)
        Per-SNP, per-block statistics
        [:, :, 0] = sum(G_resid^2) in block
        [:, :, 1:] = sum(G_resid * Y_resid) in block per tau
    """
    check_orthonormal_basis(Q, context="compute_block_scores")
    require_cpp_extension("Block-score computation")
    if HAVE_CPP:
        G_f = np.asfortranarray(G, dtype=np.float64)
        Y_f = np.asfortranarray(Y_resid, dtype=np.float64)
        Q_f = np.asfortranarray(Q, dtype=np.float64)
        block_ids_i32 = block_ids.astype(np.int32)
        if hasattr(cpp_ext, 'stage1_block_scores_raw_static'):
            return cpp_ext.stage1_block_scores_raw_static(
                G_f,
                Y_f,
                Q_f,
                block_ids_i32,
                n_blocks,
                min_mac,
                min_var,
                n_threads,
            )
        return cpp_ext.stage1_block_scores(
            G_f,
            Y_f,
            Q_f,
            block_ids_i32,
            n_blocks, min_mac, min_var, n_threads
        )
    else:
        return _compute_block_scores_numpy(G, Y_resid, Q, block_ids, n_blocks, min_mac, min_var)


def compute_block_scores_panel(G: np.ndarray, Y_resid: np.ndarray, Q: np.ndarray,
                               block_ids: np.ndarray, n_blocks: int,
                               min_mac: float = 5.0, min_var: float = 1e-8,
                               n_threads: int = 1, panel_size: int = 16) -> np.ndarray:
    """Panel form of :func:`compute_block_scores`.

    Identical output to :func:`compute_block_scores` - **bit for bit**, not merely to
    rounding.  The only difference is loop nesting: ``panel_size`` SNPs are carried
    through one sweep over the samples together, so the (N, T) residualised RIF matrix
    is streamed from DRAM once per panel instead of once per SNP.  Every floating-point
    accumulator still receives the same addends in the same order, so no reassociation
    takes place.

    ``Y_resid`` and ``Q`` are taken in C (row-major) order, which is the order NumPy
    produces them in, so the per-batch ``asfortranarray`` copies of the per-SNP path
    (166 MB + 48 MB at production N) are not made either.
    """
    check_orthonormal_basis(Q, context="compute_block_scores_panel")
    require_cpp_extension("Panel block-score computation")
    if not HAVE_CPP or not hasattr(cpp_ext, 'stage1_block_scores_panel_static'):
        raise RuntimeError(
            "The C++ extension does not export stage1_block_scores_panel_static; "
            "rebuild the extension or use compute_block_scores instead."
        )
    G_f = np.asfortranarray(G, dtype=np.float64)
    Y_c = np.ascontiguousarray(Y_resid, dtype=np.float64)
    Q_c = np.ascontiguousarray(Q, dtype=np.float64)
    block_ids_i32 = block_ids.astype(np.int32)
    return cpp_ext.stage1_block_scores_panel_static(
        G_f, Y_c, Q_c, block_ids_i32, n_blocks, min_mac, min_var, n_threads, panel_size,
    )


def compute_block_finalize_panel(G: np.ndarray, Y_resid: np.ndarray, Q: np.ndarray,
                                 block_ids: np.ndarray, n_blocks: int,
                                 min_mac: float = 5.0, min_var: float = 1e-8,
                                 n_threads: int = 1, panel_size: int = 16):
    """Panel kernel fused with the jackknife finalisation.

    The block statistics are bit-identical to :func:`compute_block_scores`; the
    finalisation is the same ``finalize_stats_flat`` used by
    :func:`compute_block_finalize`, so the beta/SE/covariance output is bit-identical to
    that function's.
    """
    check_orthonormal_basis(Q, context="compute_block_finalize_panel")
    require_cpp_extension("Fused panel block finalization")
    if not HAVE_CPP or not hasattr(cpp_ext, 'stage1_finalize_block_scores_panel_static'):
        raise RuntimeError(
            "The C++ extension does not export stage1_finalize_block_scores_panel_static; "
            "rebuild the extension or use compute_block_finalize instead."
        )
    G_f = np.asfortranarray(G, dtype=np.float64)
    Y_c = np.ascontiguousarray(Y_resid, dtype=np.float64)
    Q_c = np.ascontiguousarray(Q, dtype=np.float64)
    block_ids_i32 = block_ids.astype(np.int32)
    return cpp_ext.stage1_finalize_block_scores_panel_static(
        G_f, Y_c, Q_c, block_ids_i32, n_blocks, min_mac, min_var, n_threads, panel_size,
    )


def compute_pooled_panel(G: np.ndarray, Y_resid: np.ndarray, Q: np.ndarray,
                         min_mac: float = 5.0, min_var: float = 1e-8,
                         n_threads: int = 1, panel_size: int = 16):
    """Point estimates and the per-variant scale, with no jackknife blocks.

    Returns ``(betas, d_resid)`` with ``betas`` of shape (M, T) and ``d_resid`` of
    shape (M,), NaN for variants that fail the MAC or variance filter -- the same
    variants the block kernel leaves NaN, decided by the same two tests on the same
    inputs.

    This is the whole Stage 1 computation that ``cov_source = "pooled"`` needs:
    ``beta_t = (sum_i g_resid_i psi_r_i_t) / D`` and ``D = sum_i g_resid_i^2``, with
    ``Sigma_v = Sigma_0 / D`` supplying the covariance from a matrix that contains no
    genotype.  The ``(n_blocks, panel_size, 1+T)`` accumulator the jackknife needs is
    never allocated.

    Not bit-identical to the block path, and cannot be: the block path forms the
    totals as ``sum_b (sum_{i in b} ...)`` over blocks, this one as ``sum_i ...``
    over samples.  Same addends, different association.  Everything upstream of the
    totals -- the MAC filter, the imputed mean, ``Qtg``, and hence ``g_resid``
    itself -- is bit-for-bit what the block path computes.
    """
    check_orthonormal_basis(Q, context="compute_pooled_panel")
    require_cpp_extension("Pooled-only panel computation")
    if not HAVE_CPP or not hasattr(cpp_ext, 'stage1_pooled_panel_static'):
        raise RuntimeError(
            "The C++ extension does not export stage1_pooled_panel_static; "
            "rebuild the extension (python setup.py build_ext --inplace)."
        )
    G_f = np.asfortranarray(G, dtype=np.float64)
    Y_c = np.ascontiguousarray(Y_resid, dtype=np.float64)
    Q_c = np.ascontiguousarray(Q, dtype=np.float64)
    return cpp_ext.stage1_pooled_panel_static(
        G_f, Y_c, Q_c, min_mac, min_var, n_threads, panel_size,
    )


def compute_block_finalize(G: np.ndarray, Y_resid: np.ndarray, Q: np.ndarray,
                           block_ids: np.ndarray, n_blocks: int,
                           min_mac: float = 5.0, min_var: float = 1e-8,
                           n_threads: int = 1):
    """Fused C++ path: block statistics, jackknife beta/SE and covariance in one pass.

    Same arithmetic as ``compute_block_scores`` followed by
    ``compute_covariance_from_blocks`` per SNP, but done inside the OpenMP loop, so the
    (M, n_blocks, 1+T) intermediate is never materialised and the covariance is computed
    in parallel instead of one SNP at a time in NumPy.

    The summation order differs from the NumPy path (the C++ code accumulates the
    jackknife cross-products with explicit loops where NumPy calls BLAS), so results
    agree only to floating-point rounding, not bit-for-bit.

    Returns
    -------
    betas : array (M, T)
    se : array (M, T)
    cov_upper : array (M, T*(T+1)//2)
        Packed upper triangle of Sigma in row-major order, i.e. the same order as
        ``Sigma[np.triu_indices(T)]``.

    Rows for SNPs that fail the genotype filters (or have a non-positive total
    denominator) are NaN throughout, matching the NaN convention of
    ``compute_block_scores``.
    """
    check_orthonormal_basis(Q, context="compute_block_finalize")
    require_cpp_extension("Fused block finalization")
    if not HAVE_CPP or not hasattr(cpp_ext, 'stage1_finalize_block_scores_static'):
        raise RuntimeError(
            "The C++ extension does not export stage1_finalize_block_scores_static; "
            "rebuild the extension or use compute_block_scores instead."
        )
    G_f = np.asfortranarray(G, dtype=np.float64)
    Y_f = np.asfortranarray(Y_resid, dtype=np.float64)
    Q_f = np.asfortranarray(Q, dtype=np.float64)
    block_ids_i32 = block_ids.astype(np.int32)
    return cpp_ext.stage1_finalize_block_scores_static(
        G_f,
        Y_f,
        Q_f,
        block_ids_i32,
        n_blocks,
        min_mac,
        min_var,
        n_threads,
    )


def compute_block_scores_multi(G: np.ndarray, RIF_resid_list: list[np.ndarray], Q: np.ndarray,
                               block_ids: np.ndarray, n_blocks: int,
                               min_mac: float = 5.0, min_var: float = 1e-8,
                               n_threads: int = 1) -> list[np.ndarray]:
    """
    Compute per-block sufficient statistics for multiple phenotypes.

    Parameters
    ----------
    G : array (N, M)
        Genotype dosages (may contain NaN for missing)
    RIF_resid_list : list of arrays (N, T)
        Residualized RIF matrices
    Q : array (N, K)
        QR basis for covariates
    block_ids : array (N,)
        Block assignment for each sample (0 to n_blocks-1)
    n_blocks : int
        Number of jackknife blocks
    min_mac : float
        Minimum minor allele count filter
    min_var : float
        Minimum variance filter
    n_threads : int
        Number of OpenMP threads (if C++ available)

    Returns
    -------
    list of arrays (M, n_blocks, 1+T)
        Per-SNP, per-block statistics for each phenotype
    """
    if not isinstance(RIF_resid_list, (list, tuple)):
        raise ValueError("RIF_resid_list must be a list of arrays")

    check_orthonormal_basis(Q, context="compute_block_scores_multi")
    require_cpp_extension("Multi-phenotype block-score computation")
    if HAVE_CPP and hasattr(cpp_ext, 'stage1_block_scores_multi'):
        G_f = np.asfortranarray(G, dtype=np.float64)
        rif_f = [np.asfortranarray(rif, dtype=np.float64) for rif in RIF_resid_list]
        Q_f = np.asfortranarray(Q, dtype=np.float64)
        block_ids_i32 = block_ids.astype(np.int32)
        if hasattr(cpp_ext, 'stage1_block_scores_multi_raw_static'):
            return cpp_ext.stage1_block_scores_multi_raw_static(
                G_f,
                rif_f,
                Q_f,
                block_ids_i32,
                n_blocks,
                min_mac,
                min_var,
                n_threads,
            )
        return cpp_ext.stage1_block_scores_multi(
            G_f,
            rif_f,
            Q_f,
            block_ids_i32,
            n_blocks, min_mac, min_var, n_threads
        )

    return [
        _compute_block_scores_numpy(G, rif, Q, block_ids, n_blocks, min_mac, min_var)
        for rif in RIF_resid_list
    ]


def _compute_block_scores_numpy(G, Y_resid, Q, block_ids, n_blocks, min_mac, min_var):
    """Pure numpy fallback (slower)."""
    N, M = G.shape
    T = Y_resid.shape[1]
    
    stats = np.full((M, n_blocks, 1 + T), np.nan)
    
    for j in range(M):
        g = G[:, j].copy()
        valid = ~np.isnan(g)
        n_valid = np.sum(valid)
        
        if n_valid < 100:
            continue
            
        g_mean = np.nanmean(g)
        mac = min(np.nansum(g), 2 * n_valid - np.nansum(g))
        
        if mac < min_mac:
            continue
            
        # Mean impute missing
        g[~valid] = g_mean
        
        # Residualize on covariates
        g_resid = g - Q @ (Q.T @ g)
        
        g2 = np.sum(g_resid ** 2)
        if g2 < min_var:
            continue
        
        # Accumulate per block
        block_stats = np.zeros((n_blocks, 1 + T))
        for i in range(N):
            b = block_ids[i]
            if b >= n_blocks:
                continue
            block_stats[b, 0] += g_resid[i] ** 2
            block_stats[b, 1:] += g_resid[i] * Y_resid[i, :]
        
        stats[j] = block_stats
    
    return stats


def compute_covariance_from_blocks(stats: np.ndarray, n_blocks: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute jackknife estimates from block statistics.
    
    Parameters
    ----------
    stats : array (n_blocks, 1+T)
        Block statistics for one SNP
    n_blocks : int
        Number of jackknife blocks
        
    Returns
    -------
    beta : array (T,)
        Point estimates
    se : array (T,)  
        Standard errors
    Sigma : array (T, T)
        Full covariance matrix
    """
    T = stats.shape[1] - 1
    D_k = stats[:, 0]  # sum(G^2) per block
    N_k = stats[:, 1:]  # sum(G*Y) per block per tau
    
    D_total = np.sum(D_k)
    N_total = np.sum(N_k, axis=0)
    
    if D_total <= 0:
        return np.full(T, np.nan), np.full(T, np.nan), np.full((T, T), np.nan)
    
    # Leave-one-out betas.
    #
    # This is the delete-one-block jackknife in its O(B) subtraction form: the
    # per-block sufficient statistics are accumulated once by the C++ kernel and
    # each leave-one-out fit is recovered as (N_total - N_k) / (D_total - D_k).
    # No block-deleted fit is ever recomputed from the sample level.
    #
    # The loop this replaces ran a Python iteration per block, i.e. B * M
    # iterations per chromosome (17.5 million at B = 1000, M = 17550), each one
    # doing three NumPy calls on a length-T vector.  The vectorised form below
    # feeds every accumulator exactly the same addends: each output element is
    # the same scalar subtraction followed by the same scalar division, so it is
    # bit-for-bit identical to the loop, not merely equal to rounding.  Blocks
    # with a non-positive (or NaN) denominator are left at zero exactly as the
    # `if denom > 0` guard left them.
    denom = D_total - D_k
    usable = denom > 0
    jk_betas = np.zeros((n_blocks, T))
    np.divide(N_total - N_k, denom[:, None], out=jk_betas,
              where=usable[:, None])

    # Point estimate
    beta = N_total / D_total
    
    # Jackknife covariance
    centered = jk_betas - np.mean(jk_betas, axis=0)
    Sigma = ((n_blocks - 1) / n_blocks) * (centered.T @ centered)
    
    se = np.sqrt(np.diag(Sigma))
    
    return beta, se, Sigma


def run_stage1(G: np.ndarray, Y: np.ndarray, 
               taus: np.ndarray,
               covariates: Optional[np.ndarray] = None,
               n_blocks: Optional[int] = None,
               seed: int = 42,
               n_threads: int = 1,
               allow_numpy_fallback: bool = False) -> dict:
    """
    Run full Stage 1 quantile GWAS pipeline.
    
    Parameters
    ----------
    G : array (N, M)
        Genotype dosage matrix
    Y : array (N,)
        Phenotype
    taus : array (T,)
        Quantile levels to analyze
    covariates : array (N, K), optional
        Covariate matrix (intercept added automatically)
    n_blocks : int, optional
        Number of jackknife blocks. If None, uses DEFAULT_N_BLOCKS (1000).
    seed : int
        Random seed for block assignment
    n_threads : int
        OpenMP threads for C++ kernel
        
    Returns
    -------
    dict with keys:
        betas : array (M, T) - tau-level effect estimates
        se : array (M, T) - standard errors
        covariances : list of (T, T) arrays - full covariance per SNP
        q_tau : array (T,) - empirical quantiles
        taus : array (T,) - tau levels used
        d_resid : array (M,) - sum of squared residualised dosages per SNP
    """
    set_numpy_fallback_allowed(allow_numpy_fallback)
    require_cpp_extension("run_stage1")

    N, M = G.shape
    T = len(taus)
    n_blocks = resolve_n_blocks(n_blocks, T, flag="n_blocks")
    check_blocks_vs_samples(n_blocks, N, flag="n_blocks")
    
    # Compute RIF
    logger.info(f"Computing RIF matrix for {T} taus...")
    RIF, q_tau, fhat = compute_rif(Y, taus)
    
    # Setup covariates
    if covariates is None:
        X = np.ones((N, 1))
    else:
        X = np.column_stack([np.ones(N), covariates])
    
    # Residualize RIF
    logger.info("Residualizing RIF on covariates...")
    RIF_resid, Q = residualize(RIF, X)
    
    # Assign blocks
    np.random.seed(seed)
    block_ids = np.random.randint(0, n_blocks, N).astype(np.int32)
    n_blocks_nonempty = count_nonempty_blocks(block_ids, n_blocks, flag="n_blocks")
    
    # Compute block scores
    logger.info(f"Computing block scores for {M} SNPs...")
    all_stats = compute_block_scores(G, RIF_resid, Q, block_ids, n_blocks, n_threads=n_threads)
    
    # Aggregate results
    betas = np.zeros((M, T))
    ses = np.zeros((M, T))
    d_resid = np.full(M, np.nan)
    covariances = []
    
    for j in range(M):
        if np.isnan(all_stats[j, 0, 0]):
            betas[j] = np.nan
            ses[j] = np.nan
            covariances.append(np.full((T, T), np.nan))
        else:
            beta, se, Sigma = compute_covariance_from_blocks(all_stats[j], n_blocks)
            betas[j] = beta
            ses[j] = se
            covariances.append(Sigma)
            d_resid[j] = float(np.sum(all_stats[j][:, 0]))
    
    return {
        'betas': betas,
        'se': ses,
        'covariances': covariances,
        'q_tau': q_tau,
        'taus': taus,
        'n_blocks': n_blocks,
        'n_blocks_nonempty': n_blocks_nonempty,
        'd_resid': d_resid,
    }
