#!/usr/bin/env python
"""Validate the exact false-dispersion formula against the Monte Carlo that
drove the SHIPPED Stage 2 code.

For every MC cell we recompute, with no simulation, the predicted
  * population estimand theta_pop
  * SE ratio  reported/true
  * rejection rate at 0.05
and compare against the MC value with its Monte Carlo standard error.  A
discrepancy larger than ~3 MCSE would mean the exact engine is wrong.
"""
from __future__ import annotations
import csv, sys, math
import numpy as np
from scipy import stats
import fdr_exact as F

MCF = "/user/work/fh6520/wt_dispersion/sim/out/mc/params.tsv"

# cell name -> (maf, shape, delta, gamma, N, covar_sd)
def parse_cell(nm):
    covar = 0.0; N = 100000
    s = nm
    if s.startswith("bigN_"):
        N = 461320; s = s[5:]
    if s.startswith("cov_"):
        covar = 0.2; s = s[4:]      # placeholder; cov cells excluded below
    parts = s.split("_")
    shape = parts[0]
    maf = float(parts[1][3:])
    eff = parts[2]
    delta = gamma = 0.0
    tab = {"null": (0.0, 0.0), "loc02": (0.02, 0.0), "loc05": (0.05, 0.0),
           "loc10": (0.10, 0.0), "loc20": (0.20, 0.0), "loc35": (0.35, 0.0),
           "disp01": (0.0, 0.01), "disp02": (0.0, 0.02), "disp05": (0.0, 0.05),
           "disp10": (0.0, 0.10), "disp20": (0.0, 0.20), "mix": (0.05, 0.05)}
    delta, gamma = tab[eff]
    return maf, shape, delta, gamma, N, covar


def main():
    rows = list(csv.DictReader(open(MCF), delimiter="\t"))
    cache = {}
    out = []
    for r in rows:
        cell = r["cell"]
        if cell.startswith("cov_"):      # covariate cells: not modelled here
            continue
        if r["model"] != "both" or r["source"] not in ("pooled", "jackknife"):
            continue
        key = parse_cell(cell)
        if key not in cache:
            maf, shape, d, g, N, cv = key
            cache[key] = F.core(maf, shape, d, g, N, cv)
        k = cache[key]
        src = "pooled" if r["source"] == "pooled" else "jack"
        p = f'both.{src}.{r["estimator"]}.{r["param"]}'
        if p + ".se_ratio" not in k:
            continue
        n = int(r["n"])
        mc_rej = float(r["rej05"]); mcse = math.sqrt(mc_rej*(1-mc_rej)/n)
        ex_rej = k[p + ".rate_5e-2"]
        mc_sr = float(r["se_ratio"]); mcse_sr = mc_sr / math.sqrt(2*(n-1))
        ex_sr = k[p + ".se_ratio"]
        mc_tp = float(r["mean_theta"]); ex_tp = k[p + ".theta_pop"]
        sd_tp = float(r["emp_sd"]) / math.sqrt(n)
        out.append(dict(cell=cell, est=r["estimator"], src=r["source"],
                        param=r["param"], n=n,
                        mc_rej=mc_rej, mcse_rej=mcse, ex_rej=ex_rej,
                        z_rej=(mc_rej-ex_rej)/mcse if mcse > 0 else 0.0,
                        mc_se_ratio=mc_sr, ex_se_ratio=ex_sr,
                        z_se_ratio=(mc_sr-ex_sr)/mcse_sr,
                        mc_theta=mc_tp, ex_theta=ex_tp,
                        z_theta=(mc_tp-ex_tp)/sd_tp if sd_tp > 0 else 0.0))
    w = csv.DictWriter(sys.stdout, fieldnames=list(out[0].keys()), delimiter="\t")
    w.writeheader()
    for o in out:
        w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in o.items()})
    zr = np.array([o["z_rej"] for o in out])
    zs = np.array([o["z_se_ratio"] for o in out])
    zt = np.array([o["z_theta"] for o in out])
    sys.stderr.write(f"\nn compared = {len(out)}\n")
    for nm, z in (("rej05", zr), ("se_ratio", zs), ("theta_pop", zt)):
        sys.stderr.write(f"  {nm:10s} |z| med={np.median(np.abs(z)):6.2f} "
                         f"p95={np.quantile(np.abs(z),0.95):6.2f} max={np.abs(z).max():7.2f}\n")


if __name__ == "__main__":
    main()
