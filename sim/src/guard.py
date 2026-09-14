#!/usr/bin/env python
"""Fit and test a per-variant guard rule for the false-dispersion artifact.

The artifact is driven by an estimand gap whose leading term is

    b(tau) = delta - (delta^2/2) * gskew * f0'(q)/(f0(q) sqrt(Var g)) + ...

For a normal baseline f0'/f0 = -(q-m)/v, which is EXACTLY proportional to the
dispersion weight column (q-m)/(2v).  So the second-order part of a pure
location effect lands entirely on beta_sigma2, with a size set by the GENOTYPE
SKEWNESS gskew = (1-2p)/sqrt(2p(1-p)).  That predicts

    z_sig_false  ~=  c * gskew * Z_mu^2 / sqrt(N)

Every input is available to Stage 2 per variant: MAF -> gskew, N, and the
estimated location Z.  We fit c on the exact sweep and check its stability.
"""
import glob, csv
import numpy as np

rows = []
for f in glob.glob("/user/work/fh6520/wt_dispersion/sim/out/exact/*.tsv"):
    rows += list(csv.DictReader(open(f), delimiter="\t"))

dat = []
for r in rows:
    if r["kind"] != "loc" or r["shape"] != "normal":
        continue
    gs = float(r["geno_skew"]); N = float(r["N"]); zl = float(r["target_z"])
    zs = float(r["both.pooled.ols.beta_sigma2.z_true"])
    if gs == 0:
        continue
    dat.append((float(r["maf"]), N, zl, gs, zs, zs*np.sqrt(N)/(gs*zl**2),
                float(r["both.pooled.ols.beta_sigma2.rate_5e-8"])))

print("Stability of c = z_sig sqrt(N) / (gskew Z^2), by regime")
print("%-10s %6s %8s %8s %8s" % ("regime", "n", "median_c", "p10", "p90"))
for lab, f in (("Z<=12 (expansion valid)", lambda d: d[2] <= 12),
               ("Z<=22", lambda d: d[2] <= 22),
               ("all Z<=60", lambda d: True)):
    cs = np.array([d[5] for d in dat if f(d)])
    print("%-24s %6d %8.4f %8.4f %8.4f" % (lab, len(cs), np.median(cs),
                                           np.quantile(cs, .1), np.quantile(cs, .9)))

C = float(np.median([d[5] for d in dat if d[2] <= 12]))
print(f"\nadopted c = {C:.4f}   (1/sqrt(2) = {1/np.sqrt(2):.4f})\n")

print("Guard  G = c * gskew * Zmu^2 / sqrt(N)  = predicted false-dispersion |Z|")
print("%-7s %-9s %5s %9s %9s %11s %9s" % ("maf", "N", "Zmu", "G_pred", "z_actual",
                                          "rate@5e-8", "infl"))
for d in sorted(dat, key=lambda d: (d[0], d[1], d[2])):
    maf, N, zl, gs, zs, c, rate = d
    G = C * gs * zl**2 / np.sqrt(N)
    if zl not in (8, 12, 16, 22, 30) or N != 461320:
        continue
    print("%-7s %-9.0f %5.0f %9.4f %9.4f %11.4g %9.4g"
          % (maf, N, zl, G, zs, rate, rate/5e-8))

# where does the guard need to trip?  find the largest G with inflation <= 2
print("\nInflation of the 5e-8 rate as a function of the guard statistic G:")
print("%-9s %-12s" % ("G bin", "max inflation"))
bins = [0, .05, .1, .2, .3, .5, .75, 1.0, 1.5, 2.0, 1e9]
for lo, hi in zip(bins[:-1], bins[1:]):
    sub = [d[6]/5e-8 for d in dat if lo <= C*d[3]*d[2]**2/np.sqrt(d[1]) < hi]
    if sub:
        print("%-9s %-12.4g" % (f"[{lo},{hi})", max(sub)))
