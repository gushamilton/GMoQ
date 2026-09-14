#!/usr/bin/env python
"""Independent check of the guard against the Monte Carlo that drove the
SHIPPED Stage 2 code.  Nothing here reuses the exact rate formula: the observed
column is a raw rejection count from simulation, with its Monte Carlo SE."""
import csv, math
import numpy as np
import fdr_exact as F
from validate_exact import parse_cell

rows = list(csv.DictReader(open("/user/work/fh6520/wt_dispersion/sim/out/mc/params.tsv"),
                           delimiter="\t"))
C = 0.6855
print("%-26s %-10s %6s %7s %9s %11s %11s %8s" % (
    "cell", "source", "Zmu", "G", "mc_rej05", "mcse", "exact_rej", "n"))
seen = set()
for r in rows:
    if (r["model"] != "both" or r["param"] != "beta_sigma2"
            or r["estimator"] != "ols" or r["cell"].startswith("cov_")):
        continue
    maf, shape, d, g, N, cv = parse_cell(r["cell"])
    if g != 0.0 or d == 0.0 or shape != "normal":
        continue
    k = F.core(maf, shape, d, g, N, cv)
    zmu = abs(k["both.pooled.ols.beta_mu.z_true"])
    gs = k["geno_skew"]
    G = C * gs * zmu**2 / math.sqrt(N)
    n = int(r["n"]); mc = float(r["rej05"])
    src = "pooled" if r["source"] == "pooled" else ("jack" if r["source"] == "jackknife" else None)
    if src is None:
        continue
    ex = k[f'both.{src}.ols.beta_sigma2.rate_5e-2']
    print("%-26s %-10s %6.1f %7.4f %9.5f %11.5f %11.5f %8d" % (
        r["cell"], r["source"], zmu, G, mc, math.sqrt(mc*(1-mc)/n), ex, n))
