"""Empirical check that b_hat/se_hat is exactly invariant to the bandwidth,
and that b_hat and se_hat each scale exactly as 1/f_hat."""
import os, numpy as np, sys
sys.path.insert(0, "/user/work/fh6520/kde_lumpy/src")
from cells import build, TAUS
from shapes import make_pheno
rng = np.random.default_rng(99)
N, M = 20000, 40
for cell in ["normal_smooth", "bimodal_smooth", "round_0p50", "floor_20", "count_4"]:
    base, D = build(cell)
    G = np.stack([rng.binomial(2, 0.2, size=N).astype(float) for _ in range(M)], 1)
    Y, _ = make_pheno(base, D, G[:, 0], 0.0, 0.06, N, rng)
    cov = rng.normal(size=(N, 2))
    res = {}
    for rule in ["silverman", "scale:0.5", "scale:2.0"]:
        os.environ["FUNGWAS_BW"] = rule
        import importlib, fungwas_stage1.core as core
        importlib.reload(core)
        r = core.run_stage1(G, Y, TAUS, covariates=cov, n_blocks=1000, seed=5, n_threads=1)
        _, _, f = core.compute_rif(Y, TAUS)
        res[rule] = (r["betas"], r["se"], f)
    b0, s0, f0 = res["silverman"]
    out = [cell]
    for rule in ["scale:0.5", "scale:2.0"]:
        b, s, f = res[rule]
        sc = (f0 / f)[None, :]
        out.append(f"{rule}: max|b/(b0*sc)-1|={np.nanmax(np.abs(b/(b0*sc)-1)):.2e}"
                   f" max|se/(se0*sc)-1|={np.nanmax(np.abs(s/(s0*sc)-1)):.2e}"
                   f" max|z-z0|={np.nanmax(np.abs(b/s - b0/s0)):.2e}")
    print(" | ".join(out), flush=True)
