"""Aggregate the Monte Carlo null and the importance-sampling tail estimates."""
import glob, json
import numpy as np, pandas as pd
from scipy.stats import norm

S = "/user/work/fh6520/gmoq_power_conserv"
files = sorted(glob.glob(f"{S}/mc/plain_*.json"))
print(f"plain tasks: {len(files)}")
tot = None
for fn in files:
    d = json.load(open(fn))
    if tot is None:
        tot = dict(n=0, counts=np.zeros((3, 4)), s1=np.zeros(3), s2=np.zeros(3),
                   s4=np.zeros(3), params=d["params"], alphas=d["alphas"], crit=d["crit"])
    tot["n"] += d["n"]
    tot["counts"] += np.array(d["counts"])
    tot["s1"] += np.array(d["sum_t"]); tot["s2"] += np.array(d["sum_t2"])
    tot["s4"] += np.array(d["sum_t4"])

M = tot["n"]
rows = []
for k, p in enumerate(tot["params"]):
    var_t = tot["s2"][k] / M - (tot["s1"][k] / M) ** 2
    kurt = tot["s4"][k] / M / (tot["s2"][k] / M) ** 2
    for a_i, al in enumerate(tot["alphas"]):
        obs = tot["counts"][k, a_i]
        rate = obs / M
        se = np.sqrt(max(obs, 1e-12)) / M          # Poisson
        rows.append(dict(param=p, alpha=al, M=M, observed=int(obs), rate=rate,
                         rate_se=se, expected=al * M, ratio=rate / al,
                         ratio_se=se / al, var_t=var_t, mean_t=tot["s1"][k] / M,
                         kurt_t=kurt))
R = pd.DataFrame(rows)
R.to_csv(f"{S}/out/mc_plain.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(R.to_string())

isf = sorted(glob.glob(f"{S}/mc/is_*.tsv"))
print(f"\nIS tasks: {len(isf)}")
if isf:
    I = pd.concat([pd.read_csv(f, sep="\t") for f in isf], ignore_index=True)
    I["ratio_is"] = I.rate_is / I.alpha
    I["ratio_se"] = I.se_is / I.alpha
    I["ratio_spa"] = I.rate_spa / I.alpha
    I["is_over_spa"] = I.rate_is / I.rate_spa
    I = I.sort_values(["param", "alpha", "maf"])
    I.to_csv(f"{S}/out/mc_is.tsv", sep="\t", index=False)
    print(I.to_string())
