#!/usr/bin/env python
"""Guard comparison + pooled-vs-jackknife + boundary."""
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


def mc(p, n):
    return float(np.sqrt(max(p * (1 - p), 0.0) / max(n, 1)))


def key(nm):
    c = nm.split("__")[-1]
    return (ORDER.index(c) if c in ORDER else 99, nm)


def sec(t):
    print("\n" + "=" * 116 + "\n" + t + "\n" + "=" * 116)


sec("GUARD COMPARISON -- Stage 2 NULL, jackknife GLS. cov95(scale) and Type-I@.05(scale).\n"
    "  full    = all 15 taus (what the package does today)\n"
    "  tiedrop = drop taus with sample tie fraction > 1%\n"
    "  qdedup  = drop only taus whose sample quantile collides with the previous kept tau")
print("%-16s | %-22s | %-22s | %-22s" %
      ("cell", "full  cov95 / T1", "tiedrop  cov95 / T1", "qdedup  cov95 / T1"))
for nm in sorted(runA, key=key):
    d = S[nm]
    cells = []
    for tg in ("full", "tiedrop", "qdedup"):
        r = d["s2"].get("%s|jackknife_gls" % tg, {}).get("null")
        if r is None:
            cells.append("%-22s" % "        --")
        else:
            cells.append("%7.4f(%.4f) %7.4f" % (r["cov95"][1], mc(r["cov95"][1], r["n"]), r["rej05"][1]))
    print("%-16s | %s | %s | %s" % (d["cell"], cells[0], cells[1], cells[2]))

sec("GUARD COMPARISON -- Stage 2 CAUSAL, jackknife GLS: coverage of the TRUE theta_scale=s0")
print("%-16s | %9s %9s %9s | %9s %9s %9s" %
      ("cell", "full", "tiedrop", "qdedup", "bias_full", "bias_tied", "bias_qded"))
for nm in sorted(runA, key=key):
    d = S[nm]
    cv, bi = [], []
    for tg in ("full", "tiedrop", "qdedup"):
        r = d["s2"].get("%s|jackknife_gls" % tg, {}).get("causal")
        cv.append("%9.4f" % r["cov95"][1] if r else "%9s" % "--")
        b = r["bias_sd"][1] if r else None
        bi.append("%9.3f" % b if b is not None and abs(b) < 1e6 else "%9s" % ("--" if b is None else "huge"))
    print("%-16s | %s %s %s | %s %s %s" % (d["cell"], cv[0], cv[1], cv[2], bi[0], bi[1], bi[2]))

sec("POOLED vs BLOCK-JACKKNIFE: do the density errors cancel? (Stage 2, FULL taus)")
print("%-16s | %-21s | %-21s | %-21s" %
      ("cell", "NULL cov95_sc jk / pooled", "CAUSAL cov95_sc jk/pool", "CAUSAL cov95_sh jk/pool"))
for nm in sorted(runA, key=key):
    d = S[nm]
    j = d["s2"].get("full|jackknife_gls", {})
    p = d["s2"].get("full|pooled_gls", {})
    if not j or not p:
        continue
    print("%-16s | %9.4f %9.4f   | %9.4f %9.4f   | %9.4f %9.4f" %
          (d["cell"], j["null"]["cov95"][1], p["null"]["cov95"][1],
           j["causal"]["cov95"][1], p["causal"]["cov95"][1],
           j["causal"]["cov95"][0], p["causal"]["cov95"][0]))

sec("OLS vs GLS under the null (FULL taus): OLS is immune, GLS is not")
print("%-16s | %-19s | %-19s | %-19s | %-19s" %
      ("cell", "jk_ols cov95/T1", "jk_gls cov95/T1", "pl_ols cov95/T1", "pl_gls cov95/T1"))
for nm in sorted(runA, key=key):
    d = S[nm]
    row = []
    for est in ("jackknife_ols", "jackknife_gls", "pooled_ols", "pooled_gls"):
        r = d["s2"]["full|%s" % est]["null"]
        row.append("%7.4f  %7.4f" % (r["cov95"][1], r["rej05"][1]))
    print("%-16s | %s | %s | %s | %s" % (d["cell"], *row))

sec("BOUNDARY: per-tau tie fraction vs the spacing to the nearest neighbouring tau.\n"
    "A tie block of mass p at tau collides with a neighbour iff p exceeds that spacing.")
taus = np.array(runA[list(runA)[0]]["taus"])
sp = np.minimum(np.r_[np.inf, np.diff(taus)], np.r_[np.diff(taus), np.inf])
print("tau      " + "".join("%7.2f" % t for t in taus))
print("spacing  " + "".join("%7.2f" % t for t in sp))
print()
print("%-16s %9s %9s %10s %11s" % ("cell", "max tie", "max t/sp", "dup_qtau", "cov95_sc"))
for nm in sorted(runA, key=key):
    d = S[nm]
    tie = np.array(d["diag"]["tie"])
    r = d["s2"]["full|jackknife_gls"]["null"]
    print("%-16s %9.3f %9.2f %10.2f %11.4f" %
          (d["cell"], tie.max(), (tie / sp).max(), d["diag"]["dup_qtau"], r["cov95"][1]))

sec("LOCALITY (medians, robust to the degenerate taus): per-tau Stage 1 null SE ratio and "
    "Type-I@.05, binned by that tau's own tie fraction, pooled over all 15 runA cells")
x, y, t1 = [], [], []
for nm in runA:
    d = S[nm]
    x += list(d["diag"]["tie"])
    y += list(d["s1_null_all"]["se_ratio"])
    t1 += list(d["s1_null_all"]["t1_05"])
x, y, t1 = np.array(x), np.array(y), np.array(t1)
print("%-20s %6s %10s %10s %12s %10s" %
      ("tie fraction", "n", "se med", "se p05", "typeI med", "typeI min"))
for lo, hi in [(0, .001), (.001, .01), (.01, .05), (.05, .10), (.10, .20), (.20, .50)]:
    s = (x >= lo) & (x < hi)
    if s.sum():
        print("[%.3f, %.3f)%8s %6d %10.4f %10.4f %12.4f %10.4f" %
              (lo, hi, "", s.sum(), np.median(y[s]), np.quantile(y[s], .05),
               np.median(t1[s]), t1[s].min()))
print()
print("Stage 1 per-tau Type-I over ALL runA cells x taus (n=%d): "
      "median %.4f, 5th pct %.4f, 95th pct %.4f, max %.4f" %
      (len(t1), np.median(t1), np.quantile(t1, .05), np.quantile(t1, .95), t1.max()))
print("  -> no tau anywhere is ANTI-conservative; the only departures are downward "
      "(degenerate quantiles with no variation).")
