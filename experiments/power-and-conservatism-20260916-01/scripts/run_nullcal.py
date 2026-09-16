"""Task 2 driver: exact null rejection rates for the G-MoQ pooled/OLS Stage 2 test."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
from scipy.stats import norm, t as tdist
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nullcal import cgf, spa_two_sided, hwe, cumulants_of_S

S = "/user/work/fh6520/gmoq_power_conserv"
OUT = f"{S}/out"
os.makedirs(OUT, exist_ok=True)

# ---- inputs ---------------------------------------------------------------
Psi = np.load(f"{S}/psi/psi_chr22.npy")                     # N x T, residualised RIF
N, T = Psi.shape
A = pd.read_csv(f"{OUT}/A_contrasts.tsv", sep="\t", index_col=0)  # 3 x T
params = list(A.index)
Am = A.values.astype(np.float64)

# validate Sigma_0 against the shipped Stage 1 matrix
sg = pd.read_csv(f"/user/work/fh6520/gmoq_locovar/stage1/covariate/testosterone_chr22.sigma0.tsv.gz",
                 sep="\t")
Sig0 = sg[[c for c in sg.columns if c.startswith("s_")]].values.astype(np.float64)
Sig0_mine = (Psi.T @ Psi) / N
rel = np.max(np.abs(Sig0_mine - Sig0) / np.abs(Sig0))
print(f"[validate] N={N} T={T}  max rel err Sigma_0 (rebuilt vs shipped) = {rel:.3e}", flush=True)

U = Psi @ Am.T                                              # N x 3 contrast scores
sd_u = np.sqrt((U ** 2).sum(axis=0) / N)
print("[validate] sd_u =", sd_u, " (validate2.R: 1.97353 0.67139 0.435424)", flush=True)
print("[u] mean:", U.mean(axis=0), flush=True)
mom = {}
for k, p in enumerate(params):
    u = U[:, k]; m2 = (u**2).mean()
    mom[p] = dict(sd=float(np.sqrt(m2)),
                  skew=float((u**3).mean()/m2**1.5),
                  exkurt=float((u**4).mean()/m2**2 - 3.0),
                  max_abs_over_sd=float(np.abs(u).max()/np.sqrt(m2)))
print("[u moments]", json.dumps(mom, indent=1), flush=True)

# ---- real chr22 SNP spectrum ---------------------------------------------
s1 = pd.read_csv(f"/user/work/fh6520/gmoq_locovar/stage1/covariate/testosterone_chr22.stage1.tsv.gz",
                 sep="\t", usecols=["snp_id", "d_resid"])
rg = pd.read_csv("/user/work/fh6520/gmoq_fair_comparator/step2/loco/"
                 "testosterone_female_chr22_testosterone.regenie", sep=r"\s+")
rg = rg[rg.TEST == "ADD"][["ID", "A1FREQ", "INFO"]].rename(columns={"ID": "snp_id"})
snp = s1.merge(rg, on="snp_id")
snp["maf"] = np.minimum(snp.A1FREQ, 1 - snp.A1FREQ)
snp["mac"] = 2 * N * snp.maf
print(f"[snps] {len(snp)} matched; MAF min {snp.maf.min():.5f} "
      f"p1 {snp.maf.quantile(.01):.4f} med {snp.maf.median():.4f}; "
      f"MAC min {snp.mac.min():.0f}", flush=True)
snp["Dratio"] = snp.d_resid / (2 * snp.maf * (1 - snp.maf) * N)
print(f"[snps] d_resid / (2f(1-f)N): med {snp.Dratio.median():.4f} "
      f"p1 {snp.Dratio.quantile(.01):.4f} p99 {snp.Dratio.quantile(.99):.4f}", flush=True)

# ---- MAF bins -------------------------------------------------------------
NBIN = 30
qs = np.unique(np.quantile(snp.maf, np.linspace(0, 1, NBIN + 1)))
snp["bin"] = pd.cut(snp.maf, bins=qs, include_lowest=True, labels=False)
bins = (snp.groupby("bin")
           .agg(n=("maf", "size"), f=("maf", "median"), Dr=("Dratio", "median"),
                fmin=("maf", "min"), fmax=("maf", "max"))
           .reset_index())
print(bins.to_string(), flush=True)

ALPHAS = [0.05, 1e-3, 1e-5, 5e-8]
rows = []
for k, pn in enumerate(params):
    u = np.ascontiguousarray(U[:, k])
    su2 = float((u ** 2).sum())
    for al in ALPHAS:
        c = float(tdist.ppf(1 - al / 2, df=N))          # the package's own critical value
        c_norm = float(norm.ppf(1 - al / 2))
        for _, b in bins.iterrows():
            f = float(b.f); Dr = float(b.Dr)
            var_s = 2 * f * (1 - f) * su2
            D = Dr * 2 * f * (1 - f) * N
            thr = c * np.sqrt(D) * sd_u[k]              # threshold on |S|
            p2, pu, pl = spa_two_sided(thr, u, *hwe(f), var_s)
            pnorm_ = 2 * norm.sf(thr / np.sqrt(var_s))
            s2c, sk, ek = cumulants_of_S(u, f)
            rows.append(dict(param=pn, alpha=al, bin=int(b["bin"]), n=int(b.n),
                             maf=f, mac=2*N*f, Dratio=Dr,
                             t_crit=c, t_crit_normal=c_norm,
                             rate_exact=p2, rate_upper=pu, rate_lower=pl,
                             rate_normal_at_same_thr=pnorm_,
                             skew_S=sk, exkurt_S=ek))
        print(f"  done {pn} alpha={al:g}", flush=True)

R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/spa_bins.tsv", sep="\t", index=False)

# scan-level rate = SNP-count-weighted mean over bins
agg = (R.assign(w=lambda d: d.n / d.groupby(["param", "alpha"]).n.transform("sum"))
        .groupby(["param", "alpha"])
        .apply(lambda d: pd.Series(dict(
            rate_exact=float((d.rate_exact * d.w).sum()),
            rate_normal=float((d.rate_normal_at_same_thr * d.w).sum()),
            n_expected_17550=float((d.rate_exact * d.n).sum()))), include_groups=False)
        .reset_index())
agg["nominal"] = agg.alpha
agg["factor_vs_nominal"] = agg.rate_exact / agg.alpha
agg.to_csv(f"{OUT}/spa_scan_rates.tsv", sep="\t", index=False)
print(agg.to_string(), flush=True)
print("DONE", flush=True)
