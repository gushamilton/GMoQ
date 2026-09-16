"""Compare two Stage 1 packed covariance files (.cov.gz, float32 records).

Usage: compare_cov.py <ref.cov.gz> <new.cov.gz> <n_cov_elements>
"""
import gzip, sys
import numpy as np

a_p, b_p = sys.argv[1], sys.argv[2]
n_el = int(sys.argv[3]) if len(sys.argv) > 3 else 1035

with gzip.open(a_p, "rb") as f:
    a = np.frombuffer(f.read(), dtype=np.float32)
with gzip.open(b_p, "rb") as f:
    b = np.frombuffer(f.read(), dtype=np.float32)
assert a.size == b.size, f"sizes differ: {a.size} vs {b.size}"
n_rec = a.size // n_el
a = a.reshape(n_rec, n_el).astype(np.float64)
b = b.reshape(n_rec, n_el).astype(np.float64)

m = np.isfinite(a) & np.isfinite(b)
assert np.array_equal(np.isfinite(a), np.isfinite(b)), "NaN pattern differs"
d = np.abs(a[m] - b[m])
rel = d / np.maximum(np.abs(a[m]), np.finfo(float).tiny)
n_ulp_off = int((d > 0).sum())

# The diagonal of the packed upper triangle (entry (t, t) of a T x T matrix).
T = int(round((np.sqrt(8 * n_el + 1) - 1) / 2))
diag = []
idx = 0
for i in range(T):
    diag.append(idx)
    idx += T - i
diag = np.array(diag)
ad, bd = a[:, diag], b[:, diag]
md = np.isfinite(ad) & (ad > 0)

print(f"records={n_rec} elements/record={n_el} (T={T})")
print(f"elements compared         : {m.sum()}")
print(f"elements differing at all : {n_ulp_off} ({100.0*n_ulp_off/m.sum():.4f}%)")
print(f"max |diff|                : {d.max():.4e}")
print(f"max relative diff         : {rel.max():.4e}")
print(f"diagonal (variances) max relative diff : "
      f"{np.max(np.abs(ad[md]-bd[md])/ad[md]):.4e}")
print(f"implied max relative diff in SE = sqrt(var) : "
      f"{np.max(np.abs(np.sqrt(ad[md])-np.sqrt(bd[md]))/np.sqrt(ad[md])):.4e}")
print(f"float32 eps for reference : {np.finfo(np.float32).eps:.3e}")
