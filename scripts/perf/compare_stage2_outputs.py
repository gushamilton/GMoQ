"""Compare two Stage 2 chr22 tables column by column, and check discovery status.

Usage: compare_stage2.py <ref.tsv.gz> <new.tsv.gz>
"""
import sys
import numpy as np
import pandas as pd

a = pd.read_csv(sys.argv[1], sep="\t", dtype={"snp_id": str})
b = pd.read_csv(sys.argv[2], sep="\t", dtype={"snp_id": str})
assert list(a.columns) == list(b.columns)
assert len(a) == len(b)
assert (a["snp_id"].values == b["snp_id"].values).all()

num = [c for c in a.columns if pd.api.types.is_numeric_dtype(a[c]) and c not in ("chr", "bp")]
print(f"rows={len(a)}  numeric columns compared={len(num)}")
print(f"{'column':<16} {'max |diff|':>12} {'max relative':>13} {'n differing':>12}")
worst = 0.0
for c in num:
    x = a[c].to_numpy(float)
    y = b[c].to_numpy(float)
    assert np.array_equal(np.isnan(x), np.isnan(y)), f"{c}: NaN pattern differs"
    m = np.isfinite(x) & np.isfinite(y)
    d = np.abs(x[m] - y[m])
    rel = d / np.maximum(np.abs(x[m]), np.finfo(float).tiny)
    worst = max(worst, float(rel.max()) if rel.size else 0.0)
    print(f"{c:<16} {d.max():>12.3e} {rel.max():>13.3e} {int((d>0).sum()):>12}")
print(f"\nworst relative difference over all numeric columns: {worst:.3e}")

print("\n=== discovery status ===")
for pcol in [c for c in a.columns if c.endswith("_p") or c == "Q_pval"]:
    x = a[pcol].to_numpy(float)
    y = b[pcol].to_numpy(float)
    for thr in (5e-8, 1e-5):
        sa, sb = x < thr, y < thr
        flips = int((sa != sb).sum())
        print(f"{pcol:<14} thr={thr:<8g} ref hits={int(sa.sum()):<6} new hits={int(sb.sum()):<6} "
              f"status changes={flips}")
        if flips:
            idx = np.where(sa != sb)[0]
            for i in idx[:10]:
                print(f"    {a['snp_id'][i]}  ref={x[i]:.6e}  new={y[i]:.6e}")
