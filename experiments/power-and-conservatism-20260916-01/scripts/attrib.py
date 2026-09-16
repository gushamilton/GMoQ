"""Attribution arithmetic: where could conservatism come from, and what would it cost."""
import numpy as np, pandas as pd, json
from scipy.stats import norm, t as tdist

S = "/user/work/fh6520/gmoq_power_conserv"; OUT = f"{S}/out"
N = 399710
out = {}

# 1. reference distribution: the pooled path's t_399710 against a normal
for al in (0.05, 1e-3, 1e-5, 5e-8):
    c = tdist.ppf(1 - al / 2, df=N)
    out[f"t_vs_normal_{al:g}"] = float(2 * norm.sf(c) / al)
# and what the per-SNP jackknife df law WOULD give if it applied here
for nu in (802, 911, 3000):
    c = tdist.ppf(1 - 5e-8 / 2, df=nu)
    out[f"if_nu_{nu}_at_5e-8_normal_rate"] = float(2 * norm.sf(c) / 5e-8)

# 2. contrast-score moments
Psi = np.load(f"{S}/psi/psi_chr22.npy")
A = pd.read_csv(f"{OUT}/A_contrasts.tsv", sep="\t", index_col=0)
U = Psi @ A.values.astype(np.float64).T
del Psi
mom = {}
for k, p in enumerate(A.index):
    u = U[:, k]; m2 = (u ** 2).mean()
    mom[p] = dict(mean=float(u.mean()), sd=float(np.sqrt(m2)),
                  skew=float((u ** 3).mean() / m2 ** 1.5),
                  exkurt=float((u ** 4).mean() / m2 ** 2 - 3),
                  max_over_sd=float(np.abs(u).max() / np.sqrt(m2)),
                  sum_over_sd=float(u.sum() / np.sqrt(m2)))
out["u_moments"] = mom

# 3. the critical value a conservatism of factor k would effectively impose.
#    cost.R turns these into hits lost on the real chr21+chr22 scans.
cost = []
for k in (1.0, 1.5, 2.5, 5.0):
    zc = norm.isf(5e-8 / (2 * k))
    cost.append(dict(conservatism_k=k, effective_crit=float(zc),
                     crit_at_nominal=float(norm.isf(2.5e-8)),
                     extra_z_needed=float(zc - norm.isf(2.5e-8))))
out["cost_of_conservatism"] = cost

print(json.dumps(out, indent=1))
with open(f"{OUT}/attribution.json", "w") as fh:
    json.dump(out, fh, indent=1)
