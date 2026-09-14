#!/usr/bin/env python
"""A2. Exact Type-I of the dispersion statistic under a PURE LOCATION effect
(true dispersion effect exactly zero), uncorrected vs skew-corrected.

No Monte Carlo: rates are quadrature-exact, so 5e-8 is reachable.

Self-check first: reproduce the uncorrected rates published in
docs/FINDINGS_pooled_dispersion.md Table 1.
"""
import sys, os, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from exact import cell_matrices, delta_for_z, rate, ZC05, ZCGW

ap = argparse.ArgumentParser()
ap.add_argument("--shape", default="normal")
ap.add_argument("--N", type=int, default=461320)
ap.add_argument("--out", default=None)
a = ap.parse_args()
N = a.N

MAFS = (0.5, 0.3, 0.1, 0.05, 0.01)
ZS = (4, 8, 12, 16, 22, 30, 40, 60)

rows = []
print("=" * 126)
print(f"A2. Pure location effect, ZERO true dispersion.  shape={a.shape}  N={N}")
print("    Every rejection is FALSE.  unc = uncorrected, cor = skew-corrected.")
print("    infl = rate / nominal.  d/sd = per-allele location effect in phenotype SD units.")
print("=" * 126)
for est in ("ols", "gls"):
    print(f"\n--- estimator: {est.upper()} ---")
    print(f"{'MAF':>6} {'Z_mu':>5} {'d/sd':>7} {'G':>7} | "
          f"{'unc@.05':>9} {'cor@.05':>9} | {'unc@5e-8':>10} {'cor@5e-8':>10} | "
          f"{'infl_unc':>10} {'infl_cor':>10} | {'seinfl':>7}")
    for maf in MAFS:
        k = 1.0 - 2.0 * maf
        gskew = (1 - 2 * maf) / np.sqrt(2 * maf * (1 - maf))
        for z in ZS:
            try:
                d = delta_for_z(maf, a.shape, N, z, est=est)
            except Exception as e:
                print(f"{maf:6.3f} {z:5d}   ---- delta_for_z failed: {e}")
                continue
            M = cell_matrices(maf, a.shape, d, 0.0, N)[est]
            mean = M["theta_pop"]
            St, Sr = M["S_true"], M["S_rep"]
            sd = np.sqrt(cell_matrices(maf, a.shape, d, 0.0, N)["_"]["v0"])
            G = 0.6855 * gskew * z ** 2 / np.sqrt(N)
            r = {}
            for lvl, zc in (("05", ZC05), ("gw", ZCGW)):
                r[f"unc{lvl}"] = rate(mean, St, Sr, k, zc, correct=False)
                r[f"cor{lvl}"] = rate(mean, St, Sr, k, zc, correct=True)
            # SE inflation cost at this cell, evaluated at beta_mu = theta_pop[0]
            s1r = np.sqrt(Sr[0, 0]); s2r = np.sqrt(Sr[1, 1]); cr = Sr[0, 1]
            bm = mean[0]
            v = s2r ** 2 + 4 * k ** 2 * bm ** 2 * s1r ** 2 + 2 * k ** 2 * s1r ** 4 - 4 * k * bm * cr
            seinfl = np.sqrt(v) / s2r
            print(f"{maf:6.3f} {z:5d} {d/sd:7.4f} {G:7.3f} | "
                  f"{r['unc05']:9.4f} {r['cor05']:9.4f} | "
                  f"{r['uncgw']:10.3e} {r['corgw']:10.3e} | "
                  f"{r['unc05']/0.05:10.1f} {r['cor05']/0.05:10.2f} | {seinfl:7.3f}")
            rows.append(dict(shape=a.shape, N=N, est=est, maf=maf, z_mu=z, dsd=d / sd,
                             G=G, seinfl=seinfl, **r))
        print()

if a.out:
    keys = list(rows[0].keys())
    with open(a.out, "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for rr in rows:
            fh.write("\t".join(str(rr[kk]) for kk in keys) + "\n")
    print("wrote", a.out)
