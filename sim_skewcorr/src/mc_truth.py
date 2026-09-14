#!/usr/bin/env python
"""Exact population truth for each MC cell, and the exact predicted rates.

Writes, per cell:
  truth/truth_cellNN.txt  -- m0, v0, then the 45 null quantiles q0 (for building W)
and one combined table:
  out/mc_exact.tsv        -- the exact prediction the MC should reproduce
"""
import sys, os, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import stats
from exact import cell_matrices, rate

ap = argparse.ArgumentParser()
ap.add_argument("--cells", required=True)
ap.add_argument("--truthdir", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

os.makedirs(a.truthdir, exist_ok=True)
cells = json.load(open(a.cells))
LEVELS = [("05", 0.05), ("1e3", 1e-3)]

rows = []
for i, c in enumerate(cells):
    M = cell_matrices(c["maf"], c["shape"], c["delta"], c["gamma"], c["N"])
    info = M["_"]
    from truth import Cell
    null = Cell(c["maf"], c["shape"], 0.0, 0.0)
    q0 = null.quantile(np.linspace(0.02, 0.98, 45))
    m0, v0 = null.moments()
    with open(os.path.join(a.truthdir, f"truth_cell{i:02d}.txt"), "w") as fh:
        fh.write(f"{m0!r} {v0!r}\n")
        fh.write(" ".join(repr(x) for x in q0) + "\n")

    k = 1.0 - 2.0 * c["maf"]
    for est in ("ols", "gls"):
        E = M[est]
        mean, St, Sr = E["theta_pop"], E["S_true"], E["S_rep"]
        tgt = E["theta_dgp"]
        s1r = np.sqrt(Sr[0, 0]); s2r = np.sqrt(Sr[1, 1]); cr = Sr[0, 1]
        e_cor = mean[1] - k * (mean[0] ** 2 + St[0, 0] - s1r ** 2)
        bm = mean[0]
        se_cor = np.sqrt(s2r**2 + 4*k**2*bm**2*s1r**2 + 2*k**2*s1r**4 - 4*k*bm*cr)
        r = dict(cell=i, name=c["name"], est=est, maf=c["maf"], shape=c["shape"],
                 delta=c["delta"], gamma=c["gamma"], N=c["N"],
                 tgt_mu=tgt[0], tgt_s2=tgt[1],
                 pop_mu=mean[0], pop_s2=mean[1], e_cor=e_cor,
                 se_mu_rep=s1r, se_s2_rep=s2r, se_cor=se_cor,
                 se_mu_true=np.sqrt(St[0, 0]), se_s2_true=np.sqrt(St[1, 1]),
                 z_mu=mean[0]/np.sqrt(St[0, 0]),
                 sd_ratio_cor=se_cor/s2r)
        for lab, lv in LEVELS:
            zc = stats.norm.isf(lv/2)
            r[f"unc_{lab}"] = rate(mean, St, Sr, k, zc, False)
            r[f"cor_{lab}"] = rate(mean, St, Sr, k, zc, True)
        rows.append(r)
    print(i, c["name"], "done", flush=True)

keys = list(rows[0].keys())
with open(a.out, "w") as fh:
    fh.write("\t".join(keys) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[kk]) for kk in keys) + "\n")
print("wrote", a.out)
