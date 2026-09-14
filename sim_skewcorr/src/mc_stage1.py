"""Stage 1 replicate driver for the skew-correction validation.

Adapted from /user/work/fh6520/gls_coverage/src/sim_stage1.py: identical shipped
calls (core.run_stage1 for betas + d_resid, pooled_cov.compute_sigma0 for
Sigma_0) and identical binary output format, so gls_coverage's Stage 2 driver
pattern reads it unchanged.  Two changes only:

  * the DGP is Y = delta*gc + exp(gamma*gc) * eps, matching the exact population
    engine in src/truth.py exactly, so MC and exact are comparable cell by cell;
  * the error shapes are truth.py's (normal, skew=standardised chi2_4,
    heavy=standardised t_6), again so the two agree.
"""
import argparse, json, os, sys, time
import numpy as np

STAGE1 = os.environ.get("STAGE1_PATH",
                        "/user/work/fh6520/gls_coverage/tree/GMoQ_v04/stage1")
sys.path.insert(0, STAGE1)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fungwas_stage1 import core, pooled_cov
from truth import BASES

TAUS = np.linspace(0.02, 0.98, 45)
T = len(TAUS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True)
    ap.add_argument("--cell", type=int, required=True)
    ap.add_argument("--chunk", type=int, required=True)
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cells = json.load(open(args.cells))
    cell = cells[args.cell]
    N = int(cell["N"])
    p = float(cell["maf"])
    delta = float(cell["delta"])
    gamma = float(cell["gamma"])
    base = BASES[cell["shape"]]
    rep0 = args.chunk * args.reps

    tag = f"cell{args.cell:02d}_chunk{args.chunk:04d}"
    if os.path.exists(os.path.join(args.out, tag + ".json")):
        print(tag, "already done", flush=True)
        return

    B = np.empty((args.reps, T)); D = np.empty(args.reps)
    S0 = np.empty((args.reps, T, T)); MAC = np.empty(args.reps)

    X = np.ones((N, 1))
    Qb, _ = np.linalg.qr(X)
    t0 = time.time()
    for r in range(args.reps):
        rng = np.random.default_rng([args.cell, rep0 + r, 20260914])
        g = rng.binomial(2, p, size=N).astype(np.float64)
        gc = g - 2.0 * p
        eps = base.rvs(rng, N)
        y = delta * gc + np.exp(gamma * gc) * eps

        res = core.run_stage1(g.reshape(-1, 1), y, TAUS, n_blocks=1000)
        B[r] = res["betas"][0]; D[r] = res["d_resid"][0]; MAC[r] = g.sum()

        rif, _, _ = core.compute_rif(y, TAUS)
        rif_resid = rif - Qb @ (Qb.T @ rif)
        S0[r] = pooled_cov.compute_sigma0(rif_resid)

    os.makedirs(args.out, exist_ok=True)
    B.astype("<f8").tofile(os.path.join(args.out, tag + ".beta.bin"))
    D.astype("<f8").tofile(os.path.join(args.out, tag + ".dresid.bin"))
    S0.astype("<f8").tofile(os.path.join(args.out, tag + ".sigma0.bin"))
    MAC.astype("<f8").tofile(os.path.join(args.out, tag + ".mac.bin"))
    with open(os.path.join(args.out, tag + ".json"), "w") as fh:
        json.dump(dict(cell=cell, cell_index=args.cell, chunk=args.chunk,
                       reps=args.reps, rep0=rep0, T=T, secs=time.time() - t0), fh)
    print(f"{tag}: {args.reps} reps in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
