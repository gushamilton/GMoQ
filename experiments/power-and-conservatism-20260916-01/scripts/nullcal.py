"""Exact null calibration of the G-MoQ Stage 2 pooled/OLS parameter test.

Established by validate2.R to machine precision (max rel err 5e-15 on theta, se and p):

    theta_k = a_k' b,     se_k = sqrt(a_k' Sigma_0 a_k / D),     p = 2 * pt(-|t|, N)

with b = Psi_r' g / D and D = ||g_resid||^2, so

    t_k = (u_k' g) / (sqrt(D) * sd_u_k),   u_k = Psi_r a_k,   sd_u_k^2 = a_k' Sigma_0 a_k
                                                                       = sum_i u_ki^2 / N.

Conditioning on the realised phenotype (which is what a permutation null also does),
sd_u is EXACT, not estimated: Var(u'g) = Var(g) * sum_i u_i^2 identically for iid g.
So every calibration error at a given D comes from one place -- the departure of the
weighted genotype sum S = sum_i u_i g_i from normality. That departure is computed
here exactly by a Lugannani-Rice saddlepoint on the exact CGF, not by Edgeworth.
"""
from __future__ import annotations
import numpy as np
from scipy.stats import norm


# --------------------------------------------------------------------------
# CGF of S = sum_i u_i g_i for g_i iid on support `vals` with probs `probs`
# --------------------------------------------------------------------------
def cgf(s, u, vals, probs):
    """K(s), K'(s), K''(s). u is (N,), vals/probs are (B,)."""
    z = s * np.outer(u, vals)                 # N x B
    m = z.max(axis=1, keepdims=True)
    e = probs[None, :] * np.exp(z - m)
    den = e.sum(axis=1)
    K = (np.log(den) + m[:, 0]).sum()
    uv = u[:, None] * vals[None, :]           # N x B
    m1 = (e * uv).sum(axis=1) / den
    m2 = (e * uv * uv).sum(axis=1) / den
    return K, m1.sum(), (m2 - m1 * m1).sum()


def _solve_saddle(t, u, vals, probs, s0):
    """Newton on K'(s) = t."""
    s = s0
    for _ in range(60):
        K, k1, k2 = cgf(s, u, vals, probs)
        r = k1 - t
        if k2 <= 0 or not np.isfinite(k2):
            return np.nan, np.nan, np.nan
        step = r / k2
        s_new = s - step
        if not np.isfinite(s_new):
            return np.nan, np.nan, np.nan
        if abs(step) < 1e-13 * max(1.0, abs(s)):
            s = s_new
            break
        s = s_new
    K, k1, k2 = cgf(s, u, vals, probs)
    return s, K, k2


def spa_upper(t, u, vals, probs, var_s):
    """P(S >= t) by Lugannani-Rice. t > 0, E[S] = 0."""
    if t <= 0:
        raise ValueError("spa_upper expects t > 0")
    s, K, k2 = _solve_saddle(t, u, vals, probs, t / var_s)
    if not np.isfinite(s) or s <= 0:
        return norm.sf(t / np.sqrt(var_s)), np.nan
    arg = 2.0 * (s * t - K)
    if arg <= 0:
        return norm.sf(t / np.sqrt(var_s)), s
    w = np.sqrt(arg)
    v = s * np.sqrt(k2)
    if v <= 0:
        return norm.sf(t / np.sqrt(var_s)), s
    p = norm.sf(w) + norm.pdf(w) * (1.0 / v - 1.0 / w)
    return max(p, 0.0), s


def spa_two_sided(t, u, vals, probs, var_s):
    """P(|S| >= t). Both tails get their own saddlepoint (the law is not symmetric)."""
    pu, su = spa_upper(t, u, vals, probs, var_s)
    pl, sl = spa_upper(t, -u, vals, probs, var_s)   # P(-S >= t) = P(S <= -t)
    return pu + pl, pu, pl


# --------------------------------------------------------------------------
def hwe(f):
    return np.array([0.0, 1.0, 2.0]), np.array([(1 - f) ** 2, 2 * f * (1 - f), f ** 2])


def cumulants_of_S(u, f):
    """Exact cumulants 2,3,4 of S and the induced skew / excess kurtosis."""
    p = 2 * f * (1 - f)
    k2g = p
    k3g = p * (1 - 2 * f)
    k4g = p * (1 - 6 * f * (1 - f))
    s2 = k2g * np.sum(u ** 2)
    s3 = k3g * np.sum(u ** 3)
    s4 = k4g * np.sum(u ** 4)
    return s2, s3 / s2 ** 1.5, s4 / s2 ** 2
