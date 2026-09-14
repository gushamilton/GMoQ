#!/usr/bin/env python
"""Exact false-dispersion rate for the pooled covariance estimator.

No Monte Carlo anywhere: every number is a deterministic population quantity.

Logic
-----
Stage 2 forms theta_hat = A b_hat with A fixed (OLS) or built from the pooled
Sigma (GLS).  Asymptotically

    theta_hat ~ N( theta_pop , v_true ),      v_true = diag(A Sh A')

with Sh the TRUE HC0 sandwich, but Stage 2 REPORTS the standard error

    SE_rep = sqrt(v_rep),                     v_rep  = diag(A Sp A')

with Sp = Sigma_0 / d_resid the pooled (working-independence) covariance.
Write r = SE_rep / SE_true = sqrt(v_rep / v_true).  Then the reported statistic

    Z = theta_hat / SE_rep  ~  N( z_true / r , 1 / r^2 ),   z_true = theta_pop/sqrt(v_true)

and the two-sided rejection rate at critical value c is EXACTLY

    P(|Z| > c) = Phi(z_true - r c) + Phi(-z_true - r c).

Under a PURE LOCATION effect (gamma = 0) the dispersion parameter's true value
is zero by construction, so every rejection of beta_sigma2 is a FALSE dispersion
finding.  Two distinct mechanisms drive it and we separate them:

  * SE error only   (z_true := 0) -> 2 Phi(-r c).  N-independent.
  * estimand gap only (r := 1)    -> Phi(z_true-c)+Phi(-z_true-c).  Grows with N.
  * combined                      -> the rate Stage 2 actually pays.
"""
from __future__ import annotations

import argparse
import numpy as np
from scipy import optimize, stats

from truth import Cell, TAUS
from make_truth import sigma0_pop, sigma_hc_pop

# two-sided critical values
ALPHAS = {"5e-2": 0.05, "5e-8": 5e-8}


def _phi(x):
    """Standard normal CDF, accurate deep in the tail."""
    return stats.norm.sf(-np.asarray(x, dtype=float))


def rate_exact(z_true, r, c):
    """P(|Z| > c) with Z ~ N(z_true/r, 1/r^2).  Exact, tail-safe."""
    return float(_phi(z_true - r * c) + _phi(-z_true - r * c))


def core(maf, shape, delta, gamma, N, covar_sd=0.0):
    """All population quantities for one (maf, shape, delta, gamma, N) cell."""
    taus = TAUS
    null = Cell(maf, shape, 0.0, 0.0, covar_sd)
    cell = Cell(maf, shape, delta, gamma, covar_sd)
    q0 = null.quantile(taus)
    m0, v0 = null.moments()
    b_pop, q, f = cell.b_tau(taus)

    # design built from the NULL law = the shipped make_weights_vqtl()
    W2 = np.column_stack([np.ones(len(taus)), (q0 - m0) / (2.0 * v0)])
    W1 = np.ones((len(taus), 1))

    D = N * null.var_g
    Sp = sigma0_pop(taus, f) / D                              # pooled
    Sh = sigma_hc_pop(cell, taus, q, f) / (N * null.var_g**2)  # true sandwich

    out = dict(maf=maf, shape=shape, delta=delta, gamma=gamma, N=N)

    # ---- per-tau SE ratios (requirement 4) -------------------------------
    dr = np.sqrt(np.diag(Sp) / np.diag(Sh))       # per-tau SE ratio
    out["tau_se_ratio_min"] = float(dr.min())
    out["tau_se_ratio_p05"] = float(np.quantile(dr, 0.05))
    out["tau_se_ratio_med"] = float(np.median(dr))
    out["tau_se_ratio_p95"] = float(np.quantile(dr, 0.95))
    out["tau_se_ratio_max"] = float(dr.max())
    out["tau_se_ratio_spread"] = float(dr.max() - dr.min())

    # worst case over ALL linear contrasts of b(tau)
    ev, V = np.linalg.eigh(Sh)
    keep = ev > ev.max() * 1e-10
    H = V[:, keep] / np.sqrt(ev[keep])
    lam = np.linalg.eigvalsh(H.T @ Sp @ H)
    out["worst_se_ratio_lo"] = float(np.sqrt(max(lam.min(), 0.0)))
    out["worst_se_ratio_hi"] = float(np.sqrt(lam.max()))

    # ---- parameter-level ---------------------------------------------------
    # Three covariance SOURCES, matching the three Stage 2 can actually use:
    #   pooled -- A built from Sp (GLS) and the SE reported from Sp
    #   jack   -- the block jackknife makes NO working-independence assumption,
    #             so it converges on Sh: A built from Sh, SE reported from Sh.
    #             This is the exact asymptotic behaviour of the B=1000 jackknife
    #             and is the control that separates a COVARIANCE failure from an
    #             ESTIMAND (lack-of-fit) failure.
    # OLS does not use a covariance to build A, so its A is shared; only the
    # reported SE differs between sources.
    Spi = np.linalg.inv(Sp)
    Shi = np.linalg.inv(Sh)
    out["geno_skew"] = float((null.pr * null.gc**3).sum() / null.var_g**1.5)
    for wname, W in (("both", W2), ("add", W1)):
        A_ols = np.linalg.solve(W.T @ W, W.T)
        # GLS exactly as Stage 2 does it: weights come from the POOLED Sigma
        A_gls_p = np.linalg.solve(W.T @ Spi @ W, W.T @ Spi)
        # GLS as a jackknife-driven Stage 2 would do it: weights from the truth
        A_gls_j = np.linalg.solve(W.T @ Shi @ W, W.T @ Shi)
        arms = (("pooled", "ols", A_ols, Sp), ("pooled", "gls", A_gls_p, Sp),
                ("jack", "ols", A_ols, Sh), ("jack", "gls", A_gls_j, Sh))
        for src, est, A, Srep in arms:
            v_rep = np.diag(A @ Srep @ A.T)
            v_true = np.diag(A @ Sh @ A.T)
            th = A @ b_pop
            for i in range(W.shape[1]):
                nm = ["beta_mu", "beta_sigma2"][i]
                r = float(np.sqrt(v_rep[i] / v_true[i]))
                z_true = float(th[i] / np.sqrt(v_true[i]))
                z_rep = float(th[i] / np.sqrt(v_rep[i]))
                p = f"{wname}.{src}.{est}.{nm}"
                out[f"{p}.se_ratio"] = r
                out[f"{p}.theta_pop"] = float(th[i])
                out[f"{p}.z_true"] = z_true
                out[f"{p}.z_rep"] = z_rep
                for aname, a in ALPHAS.items():
                    c = stats.norm.isf(a / 2.0)
                    out[f"{p}.rate_{aname}"] = rate_exact(z_true, r, c)
                    out[f"{p}.rate_se_only_{aname}"] = rate_exact(0.0, r, c)
                    out[f"{p}.rate_bias_only_{aname}"] = rate_exact(z_true, 1.0, c)
    out["max_z_tau"] = float(np.max(np.abs(b_pop) / np.sqrt(np.diag(Sh))))
    return out


def solve_effect(maf, shape, N, target_z, kind, covar_sd=0.0):
    """Find the location (or dispersion) effect size giving |Z_true| = target_z
    for the corresponding parameter of the correctly-specified 'both' model."""
    def zof(x):
        d = x if kind == "loc" else 0.0
        g = 0.0 if kind == "loc" else x
        k = core(maf, shape, d, g, N, covar_sd)
        nm = "beta_mu" if kind == "loc" else "beta_sigma2"
        return abs(k[f"both.pooled.ols.{nm}.z_true"])
    lo, hi = 1e-5, 0.02
    while zof(hi) < target_z and hi < 6.0:
        hi *= 1.6
    return optimize.brentq(lambda x: zof(x) - target_z, lo, hi, xtol=1e-10, rtol=1e-12)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--N", type=int, default=461320)
    ap.add_argument("--mafs", default="0.005,0.01,0.05,0.2,0.3,0.5")
    ap.add_argument("--shapes", default="normal,skew")
    ap.add_argument("--zladder", default="3,5,8,12,16,22,30,40,50,60")
    a = ap.parse_args()

    mafs = [float(x) for x in a.mafs.split(",")]
    shapes = a.shapes.split(",")
    zl = [float(x) for x in a.zladder.split(",")]

    rows = []
    for shape in shapes:
        for maf in mafs:
            rows.append(dict(core(maf, shape, 0.0, 0.0, a.N), kind="null", target_z=0.0))
            for kind in ("loc", "disp"):
                for tz in zl:
                    try:
                        x = solve_effect(maf, shape, a.N, tz, kind)
                    except Exception as e:                       # noqa: BLE001
                        print(f"  !! {shape} maf={maf} {kind} z={tz}: {e}", flush=True)
                        continue
                    d = x if kind == "loc" else 0.0
                    g = 0.0 if kind == "loc" else x
                    rows.append(dict(core(maf, shape, d, g, a.N),
                                     kind=kind, target_z=tz))
                print(f"  {shape} maf={maf} {kind} done", flush=True)

    keys = list(rows[0].keys())
    with open(a.out, "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[k]:.12g}" if isinstance(r[k], float) else str(r[k])
                               for k in keys) + "\n")
    print(f"[fdr_exact] {len(rows)} rows -> {a.out}")


if __name__ == "__main__":
    main()
