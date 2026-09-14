#!/usr/bin/env python
"""Population-truth density / estimand tables from the Stage 1 meta json."""
import glob, json, os
import numpy as np

BASE = "/user/work/fh6520/kde_lumpy"

def arr(x):
    return np.array([np.nan if v is None else v for v in x], dtype=float)

rows = []
for d in sorted(glob.glob(os.path.join(BASE, "run[ABCD]", "*/"))):
    mf = glob.glob(os.path.join(d, "meta_*.json"))
    if not mf:
        continue
    m = json.load(open(mf[0]))
    run = os.path.basename(os.path.dirname(os.path.dirname(d)))
    rows.append((run, os.path.basename(d.rstrip("/")), m))

taus = arr(rows[0][2]["taus"])
CW = 22

def head(title):
    print()
    print("=== " + title)
    print("cell".ljust(CW) + "".join("%7.2f" % t for t in taus))

def line(nm, v, fmt="%7.3f"):
    print(nm.ljust(CW) + "".join("    nan" if not np.isfinite(x) else fmt % x for x in v))

runA = [(nm, m) for run, nm, m in rows if run == "runA"]

head("POPULATION TIE FRACTION (atom) at q_tau  [runA, Silverman, N=20000]")
for nm, m in runA:
    line(nm, arr(m["atom"]))

head("f_ship / f_true   KDE population target vs TRUE observed density (nan = atom at q_tau)")
for nm, m in runA:
    line(nm, arr(m["f_ship"]) / arr(m["f_obs_true"]))

head("f_ship / f_latent   KDE target vs the latent density the user believes in")
for nm, m in runA:
    line(nm, arr(m["f_ship"]) / arr(m["f_latent"]))

head("b_ship / b_latent   ESTIMAND SHIFT (1.000 = shipped target equals latent truth)")
for nm, m in runA:
    bl = arr(m["b_latent"]); bs = arr(m["b_ship"])
    line(nm, np.where(np.abs(bl) > 1e-10, bs / bl, np.nan))

head("b_ship / b_obs   shipped target vs the CORRECT observed-scale estimand")
for nm, m in runA:
    bo = arr(m["b_obs"]); bs = arr(m["b_ship"])
    line(nm, np.where(np.abs(bo) > 1e-10, bs / bo, np.nan))

print()
print("=== SCALARS")
print("%-28s %-11s %7s %9s %9s %9s %8s" %
      ("cell", "bw_rule", "N", "h_real", "h_pop", "sd_obs", "maxatom"))
for run, nm, m in rows:
    at = arr(m["atom"])
    print("%-28s %-11s %7d %9.5f %9.5f %9.5f %8.3f" %
          (run + "/" + nm, m.get("bw_rule", "silverman"), m["N"],
           m["h_realised"], m["h_pop"], m["sd_obs"], np.nanmax(at)))

# bandwidth-rule comparison on the same cells (runB)
print()
print("=== BANDWIDTH RULES: f_ship/f_latent at the median tau (0.50, index 7) and max |ratio-1|")
print("%-24s %-11s %10s %10s %12s" % ("cell", "bw_rule", "h", "ratio@.50", "max|ratio-1|"))
for run, nm, m in rows:
    if run not in ("runA", "runB"):
        continue
    r = arr(m["f_ship"]) / arr(m["f_latent"])
    cell = m["cell"]
    print("%-24s %-11s %10.5f %10.3f %12.3f" %
          (cell, m.get("bw_rule", "silverman"), m["h_realised"], r[7],
           np.nanmax(np.abs(r - 1))))
