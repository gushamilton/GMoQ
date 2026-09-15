"""Stage 1 replicate driver for the cov_source head-to-head.

Drives the SHIPPED Stage 1 code on simulated (g, Y) with a known population
truth and writes, per replicate, everything BOTH Stage 2 covariance sources
need:

  * betas          -- core.run_stage1()["betas"]        (shipped)
  * se_tau         -- core.run_stage1()["se"]           (shipped, jackknife)
  * d_resid        -- core.run_stage1()["d_resid"]      (shipped)
  * jackknife cov  -- core.run_stage1()["covariances"]  (shipped, B blocks)
  * Sigma_0        -- pooled_cov.compute_sigma0()       (shipped)

The jackknife covariance is stored as a row-major upper triangle in float32,
which is exactly the on-disk layout and precision of the production .cov.gz,
so Stage 2 sees the same bytes it would see in a real scan.  Sigma_0 is stored
in float64 upper triangle (the production sidecar is a text TSV, full
precision).

One replicate = one whole scan: its own phenotype, its own Sigma_0, its own
jackknife.  That is the object the default is about.
"""
import argparse, json, os, sys, time
import numpy as np

WS = "/user/work/fh6520/covsource_run"
sys.path.insert(0, WS + "/pkgsrc/stage1")
from fungwas_stage1 import core, pooled_cov

TAUS = np.linspace(0.02, 0.98, 45)
T = len(TAUS)
IU = np.triu_indices(T)            # row-major upper triangle, incl. diagonal
NU = len(IU[0])                    # 1035
N_BLOCKS = 1000

GAMMA_SHAPE = {"gamma2": 2.0, "gamma8": 8.0}


def draw_eps(rng, n, shape):
    if shape == "normal":
        return rng.normal(size=n)
    k = GAMMA_SHAPE[shape]
    x = rng.gamma(k, 1.0, size=n)
    return (x - k) / np.sqrt(k)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", default=WS + "/cells.json")
    ap.add_argument("--cell", type=int, required=True)
    ap.add_argument("--chunk", type=int, required=True)
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cells = json.load(open(args.cells))
    cell = cells[args.cell]
    N = int(cell["N"])
    rep0 = args.chunk * args.reps

    tag = f"cell{args.cell:03d}_chunk{args.chunk:04d}"
    done_marker = os.path.join(args.out, tag + ".json")
    if os.path.exists(done_marker):
        print(tag, "already done", flush=True)
        return

    B  = np.empty((args.reps, T))
    SE = np.empty((args.reps, T))
    D  = np.empty(args.reps)
    MAC = np.empty(args.reps)
    S0 = np.empty((args.reps, NU))                 # float64 upper triangle
    JK = np.empty((args.reps, NU), dtype=np.float32)

    X = np.ones((N, 1))
    Qb, _ = np.linalg.qr(X)
    t0 = time.time()
    for r in range(args.reps):
        rng = np.random.default_rng([args.cell, rep0 + r, 20260915])
        g = rng.binomial(2, cell["maf"], size=N).astype(np.float64)
        eps = draw_eps(rng, N, cell["shape"])
        y = cell["a"] * g + (1.0 + cell["c"] * g) * eps

        res = core.run_stage1(g.reshape(-1, 1), y, TAUS, n_blocks=N_BLOCKS)
        B[r]  = res["betas"][0]
        SE[r] = res["se"][0]
        D[r]  = res["d_resid"][0]
        MAC[r] = g.sum()
        JK[r] = res["covariances"][0][IU].astype(np.float32)

        rif, _, _ = core.compute_rif(y, TAUS)
        rif_resid = rif - Qb @ (Qb.T @ rif)
        S0[r] = pooled_cov.compute_sigma0(rif_resid)[IU]

    os.makedirs(args.out, exist_ok=True)
    B.astype("<f8").tofile(os.path.join(args.out, tag + ".beta.bin"))
    SE.astype("<f8").tofile(os.path.join(args.out, tag + ".setau.bin"))
    D.astype("<f8").tofile(os.path.join(args.out, tag + ".dresid.bin"))
    MAC.astype("<f8").tofile(os.path.join(args.out, tag + ".mac.bin"))
    S0.astype("<f8").tofile(os.path.join(args.out, tag + ".sigma0u.bin"))
    JK.astype("<f4").tofile(os.path.join(args.out, tag + ".jkcovu.bin"))
    with open(done_marker, "w") as fh:
        json.dump(dict(cell=cell, cell_index=args.cell, chunk=args.chunk,
                       reps=args.reps, rep0=rep0, T=T, n_upper=NU,
                       n_blocks=N_BLOCKS, secs=time.time() - t0), fh)
    print(f"{tag}: {args.reps} reps in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
