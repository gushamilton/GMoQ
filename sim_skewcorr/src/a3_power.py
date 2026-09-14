#!/usr/bin/env python
"""A3. Exact cost of the correction: null calibration, power under a GENUINE
dispersion effect, and recovery under a MIXED location+dispersion effect.

No Monte Carlo.
"""
import sys, os, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from exact import cell_matrices, delta_for_z, gamma_for_z, rate, ZC05, ZCGW

ap = argparse.ArgumentParser()
ap.add_argument("--shape", default="normal")
ap.add_argument("--N", type=int, default=461320)
a = ap.parse_args()
N = a.N
MAFS = (0.5, 0.3, 0.1, 0.05, 0.01)


def summarise(maf, shape, delta, gamma, est):
    M_all = cell_matrices(maf, shape, delta, gamma, N)
    M = M_all[est]
    k = 1.0 - 2.0 * maf
    mean, St, Sr = M["theta_pop"], M["S_true"], M["S_rep"]
    tgt = M["theta_dgp"][1]
    s1r = np.sqrt(Sr[0, 0]); s2r = np.sqrt(Sr[1, 1]); cr = Sr[0, 1]
    # E[corrected] = theta_pop[1] - k (E[bmu^2] - s1r^2),  E[bmu^2]=mean0^2+St00
    e_cor = mean[1] - k * (mean[0] ** 2 + St[0, 0] - s1r ** 2)
    bm = mean[0]
    se_cor = np.sqrt(s2r ** 2 + 4 * k ** 2 * bm ** 2 * s1r ** 2
                     + 2 * k ** 2 * s1r ** 4 - 4 * k * bm * cr)
    out = dict(
        tgt=tgt,
        bias_unc=mean[1] - tgt, bias_cor=e_cor - tgt,
        zbias_unc=(mean[1] - tgt) / s2r, zbias_cor=(e_cor - tgt) / se_cor,
        se_ratio=se_cor / s2r,
        r05_unc=rate(mean, St, Sr, k, ZC05, False),
        r05_cor=rate(mean, St, Sr, k, ZC05, True),
        rgw_unc=rate(mean, St, Sr, k, ZCGW, False),
        rgw_cor=rate(mean, St, Sr, k, ZCGW, True),
        zmu=mean[0] / np.sqrt(St[0, 0]),
    )
    return out


for est in ("ols", "gls"):
    print("=" * 122)
    print(f"A3-NULL. delta = 0, gamma = 0.  shape={a.shape} N={N} est={est.upper()}")
    print("    The correction must not break a variant that has nothing wrong with it.")
    print("=" * 122)
    print(f"{'MAF':>6} | {'unc@.05':>9} {'cor@.05':>9} | {'unc@5e-8':>10} {'cor@5e-8':>10} | {'se_ratio':>9}")
    for maf in MAFS:
        r = summarise(maf, a.shape, 0.0, 0.0, est)
        print(f"{maf:6.3f} | {r['r05_unc']:9.5f} {r['r05_cor']:9.5f} | "
              f"{r['rgw_unc']:10.3e} {r['rgw_cor']:10.3e} | {r['se_ratio']:9.5f}")
    print()

    print("=" * 122)
    print(f"A3-POWER. GENUINE dispersion effect, NO location effect.  est={est.upper()}")
    print("    Z_sig is the true dispersion Z.  bias is of beta_sigma2 about its target 2*v0*gamma.")
    print("=" * 122)
    print(f"{'MAF':>6} {'Z_sig':>6} | {'pow_unc@5e-8':>13} {'pow_cor@5e-8':>13} | "
          f"{'pow_unc@.05':>12} {'pow_cor@.05':>12} | {'relbias_unc':>12} {'relbias_cor':>12} | {'se_ratio':>9}")
    for maf in MAFS:
        for zs in (3, 5.45, 8, 12, 20, 30):
            g = gamma_for_z(maf, a.shape, N, zs, est=est)
            r = summarise(maf, a.shape, 0.0, g, est)
            rb_u = r["bias_unc"] / abs(r["tgt"]) if abs(r["tgt"]) > 0 else np.nan
            rb_c = r["bias_cor"] / abs(r["tgt"]) if abs(r["tgt"]) > 0 else np.nan
            print(f"{maf:6.3f} {zs:6.2f} | {r['rgw_unc']:13.5f} {r['rgw_cor']:13.5f} | "
                  f"{r['r05_unc']:12.5f} {r['r05_cor']:12.5f} | "
                  f"{rb_u:12.5f} {rb_c:12.5f} | {r['se_ratio']:9.5f}")
        print()

    print("=" * 122)
    print(f"A3-MIXED. BOTH a location and a dispersion effect.  est={est.upper()}")
    print("    zbias is the remaining bias in beta_sigma2 measured in its own reported SEs.")
    print("=" * 122)
    print(f"{'MAF':>6} {'Z_mu':>6} {'Z_sig':>6} | {'zbias_unc':>10} {'zbias_cor':>10} | "
          f"{'pow_unc@5e-8':>13} {'pow_cor@5e-8':>13} | {'se_ratio':>9}")
    for maf in (0.3, 0.05, 0.01):
        for zm in (8, 16, 30):
            for zs in (0.0, 5.45, 12):
                d = delta_for_z(maf, a.shape, N, zm, est=est)
                g = gamma_for_z(maf, a.shape, N, zs, est=est) if zs > 0 else 0.0
                r = summarise(maf, a.shape, d, g, est)
                print(f"{maf:6.3f} {zm:6.1f} {zs:6.2f} | {r['zbias_unc']:10.3f} {r['zbias_cor']:10.4f} | "
                      f"{r['rgw_unc']:13.5f} {r['rgw_cor']:13.5f} | {r['se_ratio']:9.5f}")
        print()
