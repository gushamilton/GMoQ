#!/usr/bin/env python
"""A1. Is the induced bias on the fitted scale coefficient exactly delta^2 (1-2p),
identically for OLS and GLS?  And where does the O(delta^3) term bite?

Prints, for a grid of MAF and delta:
   bias_ols, bias_gls  = theta_pop[1] - theta_dgp[1]
   pred                = delta^2 (1 - 2p)
   the ratio, and the induced bias on the LOCATION parameter (should be ~0).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from exact import cell_matrices

N = 461320
print("=" * 108)
print("A1. EXACT induced bias on the fitted parameters, pure location effect, normal baseline")
print("    pred = delta^2 * (1 - 2p).  d/sd is the per-allele effect in phenotype SD units.")
print("=" * 108)
hdr = (f"{'MAF':>6} {'delta':>7} {'d/sd':>7} {'bias_sig2_ols':>14} {'bias_sig2_gls':>14} "
       f"{'pred':>12} {'ols/pred':>9} {'gls/pred':>9} {'bias_mu_gls':>12}")
print(hdr)
for maf in (0.5, 0.3, 0.1, 0.05, 0.01):
    k = 1.0 - 2.0 * maf
    for d in (0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.35, 0.60, 1.00):
        M = cell_matrices(maf, "normal", d, 0.0, N)
        bo = M["ols"]["theta_pop"][1] - M["ols"]["theta_dgp"][1]
        bg = M["gls"]["theta_pop"][1] - M["gls"]["theta_dgp"][1]
        bmu = M["gls"]["theta_pop"][0] - M["gls"]["theta_dgp"][0]
        pred = d * d * k
        sd = np.sqrt(M["_"]["v0"])
        r1 = bo / pred if abs(pred) > 1e-18 else float("nan")
        r2 = bg / pred if abs(pred) > 1e-18 else float("nan")
        print(f"{maf:6.3f} {d:7.3f} {d/sd:7.3f} {bo:14.6e} {bg:14.6e} {pred:12.5e} "
              f"{r1:9.4f} {r2:9.4f} {bmu:12.4e}")
    print()

print()
print("=" * 108)
print("A1b. Does the bias vector lie EXACTLY in the column space of W?")
print("     If so the correction is identical whether applied to b(tau) or to the coefficient,")
print("     for ANY weighting.  resid_frac = ||(I-P_W) (b_pop - b_dgp)|| / ||b_pop - b_dgp||")
print("=" * 108)
print(f"{'MAF':>6} {'delta':>7} {'resid_frac':>12} {'rel_2nd_order_gap':>18}")
for maf in (0.3, 0.05, 0.01):
    for d in (0.01, 0.05, 0.20, 0.60):
        M = cell_matrices(maf, "normal", d, 0.0, N)
        W = M["_"]["W"]
        r = M["_"]["b_pop"] - M["_"]["b_dgp"]
        P = W @ np.linalg.solve(W.T @ W, W.T)
        res = r - P @ r
        # the part of r explained by the pure 2nd-order prediction
        pred_vec = d * d * (1 - 2 * maf) * W[:, 1]
        gap = np.linalg.norm(r - pred_vec) / max(np.linalg.norm(r), 1e-300)
        print(f"{maf:6.3f} {d:7.3f} {np.linalg.norm(res)/max(np.linalg.norm(r),1e-300):12.3e} {gap:18.3e}")
