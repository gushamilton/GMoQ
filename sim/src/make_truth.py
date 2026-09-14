#!/usr/bin/env python
"""Exact population quantities for one cell, written for the R side.

Also computes, with NO Monte Carlo at all, the exact population versions of the
two covariance matrices whose ratio is the whole question:

    Sigma_pooled = Sigma_0 / D,            Sigma_0[s,t] = (min(tau_s,tau_t) - tau_s tau_t)/(f_s f_t)
    Sigma_HC0    = E[gc^2 psi^c psi^c'] / (N Var(g)^2)

Sigma_0 is exact for ANY distribution of Y because psi_s is an affine function
of the indicator 1{Y <= q_s} and q_s is the unconditional tau_s-quantile, so
E[1{Y<=q_s} 1{Y<=q_t}] = min(tau_s, tau_t) whatever the shape.  That is exactly
why the pooled estimator looks unimpeachable under the null.  Sigma_HC0 is the
true sandwich and does depend on the conditional law of Y given g; the ratio of
the two IS the working-independence error, with zero Monte Carlo error.
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

from truth import Cell, TAUS


def sigma0_pop(taus, f):
    tm = np.minimum.outer(taus, taus) - np.outer(taus, taus)
    return tm / np.outer(f, f)


def sigma_hc_pop(cell, taus, q, f):
    """E[gc^2 psi^c psi^c'] -- the numerator of the exact HC0 sandwich.

    E[psi_s] = q_s exactly, so the centred influence function is simply
    psi^c_s = (tau_s - 1{Y <= q_s}) / f_s, and

        E[gc^2 psi^c_s psi^c_t]
          = sum_k pr_k gc_k^2 [ tau_s tau_t - tau_s F_t(k) - tau_t F_s(k)
                                + F_{min(s,t)}(k) ] / (f_s f_t)

    with F_s(k) = P(Y <= q_s | g = g_k).  Quantiles increase with tau, so
    min(q_s, q_t) = q_{min(s,t)}.
    """
    T = len(taus)
    imin = np.minimum.outer(np.arange(T), np.arange(T))
    out = np.zeros((T, T))
    for k in range(3):
        Fg = cell._cdf_g(q, k)                       # F(q_s | g_k), length T
        M = (np.outer(taus, taus)
             - np.outer(taus, Fg)
             - np.outer(Fg, taus)
             + Fg[imin])
        out += cell.pr[k] * cell.gc[k] ** 2 * M
    out = out / np.outer(f, f)
    return 0.5 * (out + out.T)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--maf", type=float, required=True)
    ap.add_argument("--shape", required=True)
    ap.add_argument("--delta", type=float, default=0.0)
    ap.add_argument("--gamma", type=float, default=0.0)
    ap.add_argument("--covar-sd", type=float, default=0.0)
    a = ap.parse_args()

    taus = TAUS
    null = Cell(a.maf, a.shape, 0.0, 0.0, a.covar_sd)
    cell = Cell(a.maf, a.shape, a.delta, a.gamma, a.covar_sd)

    q0 = null.quantile(taus)
    m0, v0 = null.moments()
    b_pop, q, f = cell.b_tau(taus)
    m, v = cell.moments()

    # design matrices, built from the NULL (baseline) population law -- the
    # location/log-scale design that is correctly specified for that baseline.
    # For a normal baseline (q0 - m0)/(2 v0) == qnorm(tau)/(2 sd) exactly, i.e.
    # this IS the shipped make_weights_vqtl().
    W_both = np.column_stack([np.ones(len(taus)), (q0 - m0) / (2.0 * v0)])
    W_add = np.ones((len(taus), 1))

    # first-order population profile implied by the generating parameters
    b_dgp = a.delta * np.ones(len(taus)) + a.gamma * (q0 - m0)

    S0 = sigma0_pop(taus, f)
    Shc = sigma_hc_pop(cell, taus, q, f)

    base = a.out
    for nm, arr in [("taus", taus), ("q0", q0), ("q", q), ("f", f),
                    ("bpop", b_pop), ("bdgp", b_dgp),
                    ("Wboth", W_both), ("Wadd", W_add),
                    ("S0pop", S0), ("Shcpop", Shc)]:
        np.ascontiguousarray(np.asarray(arr, dtype=np.float64)).tofile(f"{base}.{nm}.bin")
    with open(f"{base}.json", "w") as fh:
        json.dump(dict(T=len(taus), N=a.N, maf=a.maf, shape=a.shape,
                       delta=a.delta, gamma=a.gamma, covar_sd=a.covar_sd,
                       m0=m0, v0=v0, m=m, v=v, var_g=null.var_g), fh)
    print(f"[make_truth] {base}")


if __name__ == "__main__":
    main()
