"""Compare two Stage 1 chr22 outputs per column class.

Usage: compare_stage1.py <ref_prefix> <new_prefix>
where <prefix> is the path stem, e.g. .../testosterone_chr22
"""
import gzip, sys
import numpy as np
import pandas as pd

ref_p, new_p = sys.argv[1], sys.argv[2]


def load(p):
    return pd.read_csv(p + ".stage1.tsv.gz", sep="\t", dtype={"snp_id": str})


a = load(ref_p)
b = load(new_p)
assert list(a.columns) == list(b.columns), "column sets differ"
assert len(a) == len(b), f"row counts differ: {len(a)} vs {len(b)}"
assert (a["snp_id"].values == b["snp_id"].values).all(), "SNP order differs"

classes = {
    "beta_tau_*": [c for c in a.columns if c.startswith("beta_tau_")],
    "se_tau_*": [c for c in a.columns if c.startswith("se_tau_")],
    "d_resid": ["d_resid"],
}

print(f"rows={len(a)}  cols={len(a.columns)}")
print(f"{'column class':<14} {'n':>10} {'n NaN ref':>10} {'NaN same':>9} "
      f"{'max |diff|':>12} {'max rel':>12} {'n differing text':>17}")
for name, cols in classes.items():
    x = a[cols].to_numpy(dtype=float)
    y = b[cols].to_numpy(dtype=float)
    nan_same = bool(np.array_equal(np.isnan(x), np.isnan(y)))
    m = ~np.isnan(x) & ~np.isnan(y)
    d = np.abs(x[m] - y[m])
    denom = np.maximum(np.abs(x[m]), np.finfo(float).tiny)
    rel = d / denom
    ndiff = int((d > 0).sum())
    print(f"{name:<14} {m.sum():>10} {int(np.isnan(x).sum()):>10} {str(nan_same):>9} "
          f"{(d.max() if d.size else 0):>12.3e} {(rel.max() if rel.size else 0):>12.3e} "
          f"{ndiff:>17}")

# Where the largest relative differences sit, in units of the SE (the scale the
# estimator is actually read on).
bc = classes["beta_tau_*"]
sc = classes["se_tau_*"]
xb, yb = a[bc].to_numpy(float), b[bc].to_numpy(float)
xs = a[sc].to_numpy(float)
m = ~np.isnan(xb) & ~np.isnan(yb) & ~np.isnan(xs) & (xs > 0)
z = np.abs(xb[m] - yb[m]) / xs[m]
print(f"\nmax |delta beta| / SE            : {z.max():.3e}")
print(f"99.999th pct |delta beta| / SE   : {np.quantile(z, 0.99999):.3e}")
ys = b[sc].to_numpy(float)
ms = ~np.isnan(xs) & ~np.isnan(ys) & (xs > 0)
print(f"max |delta SE| / SE              : {np.max(np.abs(xs[ms]-ys[ms])/xs[ms]):.3e}")

# Exact text comparison, field by field, to see how much survives %.6g.
with gzip.open(ref_p + ".stage1.tsv.gz", "rt") as f1, gzip.open(new_p + ".stage1.tsv.gz", "rt") as f2:
    nline = ndiffield = 0
    for l1, l2 in zip(f1, f2):
        if l1 != l2:
            nline += 1
            ndiffield += sum(1 for u, v in zip(l1.split("\t"), l2.split("\t")) if u != v)
print(f"\ntext: {nline} lines differ, {ndiffield} fields differ out of {len(a)*len(a.columns)}")
