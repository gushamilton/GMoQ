"""Patch GMoQ_v04 stage1 core.compute_rif to take the bandwidth rule from an env var.

FUNGWAS_BW unset / "silverman"  -> byte-identical to the shipped rule
           "robust"             -> 0.9 * min(sd, IQR/1.349) * N^-0.2
           "scale:X"            -> X * shipped
           "isj"                -> improved Sheather-Jones (Botev et al. 2010)
"""
import re, sys, os

PATCH = '''

# ------------------------------------------------------------------ BW PATCH
import os as _os
from scipy.fft import dct as _dct


def _isj_bandwidth(Y, n_bins=2 ** 14):
    """Improved Sheather-Jones plug-in bandwidth (Botev, Grotowski & Kroese 2010)."""
    Y = np.asarray(Y, dtype=np.float64)
    N = len(Y)
    lo, hi = Y.min(), Y.max()
    rng = hi - lo
    if not np.isfinite(rng) or rng <= 0:
        return None
    lo -= rng / 10.0
    hi += rng / 10.0
    R = hi - lo
    idx = np.clip(((Y - lo) / R * n_bins).astype(np.int64), 0, n_bins - 1)
    hist = np.bincount(idx, minlength=n_bins).astype(np.float64)
    hist /= hist.sum()
    a = _dct(hist, type=2)
    I_sq = np.arange(1, n_bins, dtype=np.float64) ** 2
    a2 = (a[1:] / 2.0) ** 2

    def _f(s, t):
        y = I_sq ** s * a2 * np.exp(-I_sq * np.pi ** 2 * t)
        return 2.0 * np.pi ** (2.0 * s) * np.sum(y)

    def _fixed_point(t):
        l = 7
        f = _f(l, t)
        for s in range(l - 1, 1, -1):
            K0 = np.prod(np.arange(1, 2 * s, 2)) / np.sqrt(2 * np.pi)
            const = (1 + 0.5 ** (s + 0.5)) / 3.0
            if f <= 0:
                return np.nan
            tt = (2.0 * const * K0 / (N * f)) ** (2.0 / (3.0 + 2.0 * s))
            f = _f(s, tt)
        if f <= 0:
            return np.nan
        return t - (2.0 * N * np.sqrt(np.pi) * f) ** (-2.0 / 5.0)

    # bracket and bisect
    x0 = 1e-12
    x1 = 0.1
    f0 = _fixed_point(x0)
    for _ in range(60):
        f1 = _fixed_point(x1)
        if np.isfinite(f0) and np.isfinite(f1) and f0 * f1 < 0:
            break
        x1 *= 0.5
        if x1 < 1e-14:
            return None
    else:
        return None
    for _ in range(200):
        xm = 0.5 * (x0 + x1)
        fm = _fixed_point(xm)
        if not np.isfinite(fm):
            return None
        if f0 * fm <= 0:
            x1 = xm
        else:
            x0, f0 = xm, fm
        if x1 - x0 < 1e-16:
            break
    t_star = 0.5 * (x0 + x1)
    return float(np.sqrt(t_star) * R)


def _bandwidth(Y, N):
    rule = _os.environ.get("FUNGWAS_BW", "silverman") or "silverman"
    silver = 1.06 * np.std(Y) * N ** (-0.2)
    if rule == "silverman":
        return silver
    if rule == "robust":
        q75, q25 = np.percentile(Y, [75, 25])
        return 0.9 * min(np.std(Y), (q75 - q25) / 1.349) * N ** (-0.2)
    if rule.startswith("scale:"):
        return float(rule.split(":", 1)[1]) * silver
    if rule == "isj":
        h = _isj_bandwidth(Y)
        if h is None or not np.isfinite(h) or h <= 0:
            logger.warning("ISJ bandwidth failed; falling back to Silverman")
            return silver
        return h
    raise ValueError(f"unknown FUNGWAS_BW rule: {rule}")
# ---------------------------------------------------------------- END PATCH
'''


def main(path):
    s = open(path).read()
    assert "h = 1.06 * np.std(Y) * N ** (-0.2)" in s, "anchor not found"
    s = s.replace("    h = 1.06 * np.std(Y) * N ** (-0.2)",
                  "    h = _bandwidth(Y, N)")
    # insert helper before compute_rif
    s = s.replace("def compute_rif(Y: np.ndarray, taus: np.ndarray, ",
                  PATCH + "\n\ndef compute_rif(Y: np.ndarray, taus: np.ndarray, ")
    open(path, "w").write(s)
    print("patched", path)


if __name__ == "__main__":
    main(sys.argv[1])
