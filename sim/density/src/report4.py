#!/usr/bin/env python
"""Tail / genome-wide behaviour under the null."""
import glob, json, os
import numpy as np
from scipy import stats as sps

SUMM = "/user/work/fh6520/kde_lumpy/summ"
ORDER = ["normal_smooth", "lognorm_smooth", "bimodal_smooth", "bimodal_uneven",
         "round_0p10", "round_0p25", "round_0p50", "round_1p00",
         "floor_05", "floor_20", "floor_40", "ceil_20",
         "count_8", "count_4", "bimodal_round"]
S = {}
for f in sorted(glob.glob(os.path.join(SUMM, "*.json"))):
    nm = os.path.basename(f)[:-5]
    if nm != "TEST":
        S[nm] = json.load(open(f))
runA = {k: v for k, v in S.items() if k.startswith("runA__")}


def key(nm):
    c = nm.split("__")[-1]
    return (ORDER.index(c) if c in ORDER else 99, nm)


def mc(p, n):
    return float(np.sqrt(max(p * (1 - p), 0.0) / max(n, 1)))


print("=" * 112)
print("STAGE 2 NULL TAIL, scale parameter. n = 200,000 null SNPs per cell.")
print("Observed rejection at alpha, as a MULTIPLE of alpha (1.00 = exact). "
      "MC error on the ratio in brackets.")
print("At 1e-4 the expectation is 20 hits, so the ratio has ~22% MC error; "
      "1e-6 (0.2 expected) is indicative only.")
print("=" * 112)
print("%-16s %-14s | %-16s %-16s %-16s" %
      ("cell", "estimator", "a=1e-3", "a=1e-4", "a=1e-6"))
for nm in sorted(runA, key=key):
    d = S[nm]
    for est in ("jackknife_ols", "jackknife_gls", "pooled_ols", "pooled_gls"):
        r = d["s2"]["full|%s" % est]["null"]
        n = r["n"]
        cells = []
        for lab, al in (("rej_1e3", 1e-3), ("rej_1e4", 1e-4), ("rej_1e6", 1e-6)):
            p = r[lab][1]
            cells.append("%6.2f (%.2f)" % (p / al, mc(p, n) / al))
        print("%-16s %-14s | %s %s %s" % (d["cell"], est, *cells))

print()
print("=" * 112)
print("STAGE 1 PER-TAU NULL TAIL: empirical |z| quantile / t_999 reference quantile. "
      "n = 200,000 x 100 null SNPs.")
print("1.00 = the reference distribution is exactly right in that tail.")
print("=" * 112)
taus = np.array(runA[list(runA)[0]]["taus"])
for lev in ("0.9999", "0.99999"):
    print("\nlevel %s" % lev)
    print("tau".ljust(18) + "".join("%7.2f" % t for t in taus))
    for nm in sorted(runA, key=key):
        d = S[nm]
        zq = d["s1_null_all"]["z_q"][lev]
        print(("%-12s" % d["cell"][:12]).ljust(18) +
              "".join("%7.3f" % x for x in zq["ratio"]))
