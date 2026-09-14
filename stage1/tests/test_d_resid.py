"""Tests for the per-SNP jackknife denominator D = sum_i g_resid_i^2.

Stage 2 turns D into the effective degrees of freedom of the estimated
covariance, per SNP, via

    1/nu = 1/(B - 1) + c / N_eff,      N_eff = D^2 / sum_i g_resid_i^4 ~ D

so D has to be exactly the quantity the jackknife divides by, and it has to
reach the output file.
"""

import numpy as np
import pytest

from fungwas_stage1 import cov_meta
from fungwas_stage1.core import (
    compute_covariance_from_blocks,
    residual_sum_squares,
    run_stage1,
)


def _orthonormal_covariate_basis(N, K, rng):
    X = np.column_stack([np.ones(N), rng.standard_normal((N, K))])
    Q, _ = np.linalg.qr(X)
    return Q


def test_residual_sum_squares_matches_an_explicit_projection():
    rng = np.random.default_rng(0)
    N, M = 300, 7
    G = rng.binomial(2, 0.3, size=(N, M)).astype(float)
    Q = _orthonormal_covariate_basis(N, 2, rng)

    got = residual_sum_squares(G, Q)
    for j in range(M):
        g = G[:, j]
        r = g - Q @ (Q.T @ g)
        assert got[j] == pytest.approx(float(r @ r), rel=1e-12)


def test_residual_sum_squares_mean_imputes_missing_dosages():
    rng = np.random.default_rng(1)
    N = 200
    Q = _orthonormal_covariate_basis(N, 1, rng)
    g = rng.binomial(2, 0.4, size=N).astype(float)
    g_with_nan = g.copy()
    g_with_nan[[3, 17, 42]] = np.nan

    imputed = g_with_nan.copy()
    imputed[np.isnan(imputed)] = np.nanmean(g_with_nan)
    r = imputed - Q @ (Q.T @ imputed)

    got = residual_sum_squares(g_with_nan[:, None], Q)
    assert got[0] == pytest.approx(float(r @ r), rel=1e-12)


def test_d_resid_is_the_jackknife_denominator():
    """D must be the same number compute_covariance_from_blocks divides by."""
    rng = np.random.default_rng(2)
    n_blocks, T = 40, 6
    stats = np.empty((n_blocks, 1 + T))
    stats[:, 0] = rng.uniform(0.5, 3.0, n_blocks)
    stats[:, 1:] = rng.normal(size=(n_blocks, T))

    beta, _, _ = compute_covariance_from_blocks(stats, n_blocks)
    d_total = float(np.sum(stats[:, 0]))
    # beta_t = (sum_k N_kt) / D, so D is recoverable from any tau
    n_total = np.sum(stats[:, 1:], axis=0)
    assert np.allclose(beta, n_total / d_total, rtol=1e-12)


def test_run_stage1_returns_d_resid_matching_the_genotypes():
    rng = np.random.default_rng(3)
    N, M, T = 1500, 12, 7
    G = rng.binomial(2, 0.25, size=(N, M)).astype(float)
    covar = rng.standard_normal((N, 2))
    Y = rng.standard_normal(N) + 0.1 * covar[:, 0]
    taus = np.linspace(0.2, 0.8, T)

    out = run_stage1(G, Y, taus, covariates=covar, n_blocks=120, seed=11, n_threads=1)

    assert "d_resid" in out
    assert out["d_resid"].shape == (M,)
    assert np.all(np.isfinite(out["d_resid"]))
    assert np.all(out["d_resid"] > 0)

    X = np.column_stack([np.ones(N), covar])
    Q, _ = np.linalg.qr(X)
    expected = residual_sum_squares(G, Q)
    assert np.allclose(out["d_resid"], expected, rtol=1e-8)


def test_d_resid_tracks_the_genotype_effective_sample_size():
    """N_eff = D^2/Q4 collapses to D for a genotype, which is why one float is enough."""
    rng = np.random.default_rng(4)
    N = 60000
    Q = np.ones((N, 1)) / np.sqrt(N)
    for maf in (0.01, 0.05, 0.2, 0.5):
        g = rng.binomial(2, maf, size=N).astype(float)
        r = g - Q @ (Q.T @ g)
        D = float(r @ r)
        n_eff = D ** 2 / float(np.sum(r ** 4))
        assert n_eff / D == pytest.approx(1.0, abs=0.05)
        # and D itself tracks 2 N p(1-p)
        assert D / (2 * N * maf * (1 - maf)) == pytest.approx(1.0, abs=0.05)


def test_effective_nu_falls_with_maf_through_d_resid():
    """The MAF dependence of nu is entirely carried by D."""
    N, B = 461320, 1000
    nus = []
    for maf in (0.01, 0.05, 0.20, 0.50):
        d = 2 * N * maf * (1 - maf)
        nus.append(cov_meta.effective_nu(B, d, cov_meta.C_LOCATION))
    assert nus == sorted(nus)
    # production values quoted in BP1_effective_df.md for the location parameter
    assert round(nus[0]) == 911
    assert round(nus[1]) == 979
    assert round(nus[3]) == 995


def test_about_1100_blocks_would_make_nu_1000_at_production():
    """Recorded because it is the block count that would make nu = 1000 real.

    The worst case is the rarest variant, MAF 0.01 at production N. B = 1000
    leaves nu at 911 there; the law says B = 1108 reaches a true 1000, so
    ~1100 is the block count to revisit if a round nu is wanted.
    """
    N = 461320
    d_rare = 2 * N * 0.01 * 0.99
    assert cov_meta.effective_nu(1000, d_rare, cov_meta.C_LOCATION) < 1000
    assert cov_meta.effective_nu(1108, d_rare, cov_meta.C_LOCATION) >= 1000
    # ~1100 is within 1% of the target across the MAF range
    for maf in (0.01, 0.05, 0.5):
        d = 2 * N * maf * (1 - maf)
        assert cov_meta.effective_nu(1100, d, cov_meta.C_LOCATION) > 990
