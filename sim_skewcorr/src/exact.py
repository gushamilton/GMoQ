#!/usr/bin/env python
"""Exact (no Monte Carlo) behaviour of the genotype-skewness bias correction.

For each cell we compute, to quadrature accuracy:

  * theta_pop  -- the exact Stage-2 estimand (A @ b_pop) for OLS and GLS
  * theta_dgp  -- what the generating parameters say it should be
  * S_true = A Sh A'   -- the exact sampling covariance of theta_hat
  * S_rep  = A Sp A'   -- the covariance Stage 2 REPORTS (pooled)
  * the exact rejection rate of the uncorrected and the skew-corrected
    dispersion statistic at any nominal level, by 1-D quadrature over
    beta_mu_hat (the corrected statistic is a quadratic form in a bivariate
    normal, so conditioning on beta_mu_hat leaves a normal tail probability).

The correction is
    T = beta_sigma2_hat - k * (beta_mu_hat^2 - se_mu_rep^2),     k = 1 - 2p
    Var(T) = se_sigma2^2 + 4 k^2 beta_mu^2 se_mu^2 + 2 k^2 se_mu^4
             - 4 k beta_mu cov_mu_sigma2
all of which use the REPORTED (pooled) covariance, exactly as Stage 2 would.
"""
from __future__ import annotations

import sys
import os
import numpy as np
from scipy import optimize, stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from truth import Cell, TAUS
from make_truth import sigma0_pop, sigma_hc_pop

ZC05 = stats.norm.isf(0.05 / 2)
ZCGW = stats.norm.isf(5e-8 / 2)


def cell_matrices(maf, shape, delta, gamma, N, covar_sd=0.0):
    """Everything Stage 2 would see, exactly."""
    taus = TAUS
    null = Cell(maf, shape, 0.0, 0.0, covar_sd)
    cell = Cell(maf, shape, delta, gamma, covar_sd)
    q0 = null.quantile(taus)
    m0, v0 = null.moments()
    b_pop, q, f = cell.b_tau(taus)
    b_dgp = delta * np.ones(len(taus)) + gamma * (q0 - m0)
    # exactly make_weights_vqtl() for a normal baseline
    W = np.column_stack([np.ones(len(taus)), (q0 - m0) / (2.0 * v0)])

    D = N * null.var_g
    Sp = sigma0_pop(taus, f) / D                                   # pooled (reported)
    Sh = sigma_hc_pop(cell, taus, q, f) / (N * null.var_g ** 2)    # true sandwich

    out = {}
    Spi = np.linalg.inv(Sp)
    for est in ("ols", "gls"):
        if est == "ols":
            A = np.linalg.solve(W.T @ W, W.T)
        else:
            A = np.linalg.solve(W.T @ Spi @ W, W.T @ Spi)
        out[est] = dict(
            A=A,
            theta_pop=A @ b_pop,
            theta_dgp=A @ b_dgp,
            S_true=A @ Sh @ A.T,
            S_rep=A @ Sp @ A.T,
        )
    out["_"] = dict(maf=maf, shape=shape, delta=delta, gamma=gamma, N=N,
                    v0=v0, var_g=null.var_g, W=W, b_pop=b_pop, b_dgp=b_dgp)
    return out


def rate(mean, S_true, S_rep, k, zc, correct, ngrid=40001, span=16.0,
         debias=True):
    """Exact P(reject) for the (un)corrected dispersion statistic."""
    m1, m2 = mean
    s1t = np.sqrt(S_true[0, 0])
    s2t = np.sqrt(S_true[1, 1])
    rho = S_true[0, 1] / (s1t * s2t)
    s1r = np.sqrt(S_rep[0, 0])
    s2r = np.sqrt(S_rep[1, 1])
    cr = S_rep[0, 1]

    if not correct:
        return float(stats.norm.sf((zc * s2r - m2) / s2t)
                     + stats.norm.cdf((-zc * s2r - m2) / s2t))

    x = np.linspace(m1 - span * s1t, m1 + span * s1t, ngrid)
    wdens = stats.norm.pdf(x, m1, s1t)
    cmu = m2 + rho * (s2t / s1t) * (x - m1)
    csd = s2t * np.sqrt(max(1.0 - rho ** 2, 0.0))
    corr = k * (x ** 2 - (s1r ** 2 if debias else 0.0))
    varT = (s2r ** 2 + 4 * k ** 2 * x ** 2 * s1r ** 2 + 2 * k ** 2 * s1r ** 4
            - 4 * k * x * cr)
    seT = np.sqrt(np.maximum(varT, 1e-300))
    hi = corr + zc * seT
    lo = corr - zc * seT
    if csd <= 0:
        inner = ((cmu > hi) | (cmu < lo)).astype(float)
    else:
        inner = stats.norm.sf((hi - cmu) / csd) + stats.norm.cdf((lo - cmu) / csd)
    trap = getattr(np, "trapezoid", np.trapz)
    return float(trap(wdens * inner, x))


def delta_for_z(maf, shape, N, z_target, est="gls", covar_sd=0.0):
    """Find the location effect whose TRUE location Z equals z_target."""
    def fz(d):
        M = cell_matrices(maf, shape, d, 0.0, N, covar_sd)[est]
        return M["theta_pop"][0] / np.sqrt(M["S_true"][0, 0]) - z_target
    lo, hi = 1e-7, 0.02
    while fz(hi) < 0 and hi < 6.0:
        hi *= 1.6
    return optimize.brentq(fz, lo, hi, xtol=1e-12, rtol=1e-14)


def gamma_for_z(maf, shape, N, z_target, est="gls", covar_sd=0.0):
    """Find the log-scale effect whose TRUE dispersion Z equals z_target."""
    def fz(g):
        M = cell_matrices(maf, shape, 0.0, g, N, covar_sd)[est]
        return M["theta_pop"][1] / np.sqrt(M["S_true"][1, 1]) - z_target
    lo, hi = 1e-7, 0.02
    while fz(hi) < 0 and hi < 6.0:
        hi *= 1.6
    return optimize.brentq(fz, lo, hi, xtol=1e-12, rtol=1e-14)
