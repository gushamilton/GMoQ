"""Task 2, self-consistent version.

The first pass mixed two things: the shape of the null law of S = sum_i u_i g_i, and a
variance-model mismatch (I set the threshold with the REAL d_resid, whose median is
0.5% below 2f(1-f)N because real dosages are imputed, while modelling g as a hard HWE
genotype with variance exactly 2f(1-f)). That 0.5% is a property of the genotype data,
not of the test: conditional on the realised genotype vector -- which is what a
permutation null does, and what the package's own D is -- Var(t) is exactly N/(N-1).

So here the threshold is set self-consistently, D = E[D] = (N-K)*Var(g), and the only
thing left in rate_exact/nominal is the departure of S from normality.
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
from scipy.stats import norm, t as tdist
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nullcal import spa_two_sided, hwe, cumulants_of_S

S = "/user/work/fh6520/gmoq_power_conserv"; OUT = f"{S}/out"
Psi = np.load(f"{S}/psi/psi_chr22.npy"); N, T = Psi.shape
Q = np.load(f"{S}/psi/Q_chr22.npy"); K = Q.shape[1]
A = pd.read_csv(f"{OUT}/A_contrasts.tsv", sep="\t", index_col=0)
params = list(A.index)
U = Psi @ A.values.astype(np.float64).T
sd_u = np.sqrt((U ** 2).sum(axis=0) / N); su2 = (U ** 2).sum(axis=0)
del Psi

bins = pd.read_csv(f"{OUT}/spa_bins.tsv", sep="\t")
bins = bins[(bins.param == params[0]) & (bins.alpha == 0.05)][["bin", "n", "maf", "Dratio"]]

ALPHAS = [0.05, 1e-3, 1e-5, 5e-8]
rows = []
for k, pn in enumerate(params):
    u = np.ascontiguousarray(U[:, k])
    for al in ALPHAS:
        c = float(tdist.ppf(1 - al / 2, df=N))
        for _, b in bins.iterrows():
            f = float(b.maf)
            vals, probs = hwe(f)
            var_g = 2 * f * (1 - f)
            var_s = var_g * su2[k]
            D = var_g * (N - K)                      # self-consistent E[D]
            thr = c * np.sqrt(D) * sd_u[k]
            p2, pu, pl = spa_two_sided(thr, u, vals, probs, var_s)
            pn_ = 2 * norm.sf(thr / np.sqrt(var_s))
            _, sk, ek = cumulants_of_S(u, f)
            rows.append(dict(param=pn, alpha=al, bin=int(b["bin"]), n=int(b.n), maf=f,
                             mac=2 * N * f, rate_exact=p2, rate_upper=pu, rate_lower=pl,
                             rate_normal=pn_, shape_factor=p2 / pn_,
                             total_factor=p2 / al, skew_S=sk, exkurt_S=ek))
        print(f"  {pn} alpha={al:g}", flush=True)
R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/spa_bins_selfconsistent.tsv", sep="\t", index=False)

R["w"] = R.n / R.groupby(["param", "alpha"]).n.transform("sum")
agg = (R.groupby(["param", "alpha"])
        .apply(lambda d: pd.Series(dict(
            rate_exact=float((d.rate_exact * d.w).sum()),
            rate_normal=float((d.rate_normal * d.w).sum()),
            worst_bin_maf=float(d.loc[(d.rate_exact / d.rate_normal).abs().idxmax(), "maf"]),
            worst_shape_factor=float((d.rate_exact / d.rate_normal).max()),
            best_shape_factor=float((d.rate_exact / d.rate_normal).min()))), include_groups=False)
        .reset_index())
agg["factor_vs_nominal"] = agg.rate_exact / agg.alpha
agg["expected_hits_17550"] = agg.rate_exact * 17550
agg.to_csv(f"{OUT}/spa_scan_selfconsistent.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(agg.to_string(), flush=True)
print("\n--- 5e-8, by MAF bin ---", flush=True)
print(R[R.alpha == 5e-8].pivot_table(index=["bin", "maf", "mac"], columns="param",
                                     values="shape_factor").to_string(), flush=True)
print("\n--- cumulants of S at 5e-8 ---", flush=True)
print(R[R.alpha == 5e-8].pivot_table(index=["maf"], columns="param",
                                     values=["skew_S", "exkurt_S"]).to_string(), flush=True)
print("DONE", flush=True)
