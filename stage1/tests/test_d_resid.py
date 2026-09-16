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


