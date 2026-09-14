import resource
import sys
import time

import numpy as np

sys.path.insert(0, "/user/work/fh6520/wt_nscale/stage1")
from fungwas_stage1 import core, pooled_cov

TAUS = np.linspace(0.02, 0.98, 45)
for N in (30000, 100000, 461320):
    X = np.ones((N, 1))
    Qb, _ = np.linalg.qr(X)
    nr = 6 if N > 200000 else 15
    t0 = time.time()
    for r in range(nr):
        rng = np.random.default_rng([99, r, 1])
        g = rng.binomial(2, 0.30, size=N).astype(np.float64)
        y = rng.normal(size=N)
        res = core.run_stage1(g.reshape(-1, 1), y, TAUS, n_blocks=1000)
        rif, qt, fh = core.compute_rif(y, TAUS)
        rr = rif - Qb @ (Qb.T @ rif)
        S0 = pooled_cov.compute_sigma0(rr)
        del rif, rr
    dt = (time.time() - t0) / nr
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    print("N=%7d  %7.3f s/rep  peakRSS=%.2f GB  d_resid=%.1f  fhat[0]=%.5f"
          % (N, dt, rss, res["d_resid"][0], fh[0]), flush=True)
