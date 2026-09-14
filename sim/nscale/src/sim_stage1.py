"""Stage 1 replicate driver for the N-scaling coverage measurement.

Drives the SHIPPED Stage 1 code (core.run_stage1 for betas + d_resid,
pooled_cov.compute_sigma0 for Sigma_0) on simulated (g, Y) with a known
population truth.  Writes raw float64 binaries the R Stage 2 driver reads.

Identical in structure to gls_coverage/src/sim_stage1.py, with three changes:
  * workspace and code tree are arguments, not hard-coded;
  * n_blocks comes from the cell (so B is held at 1000 across the N ladder);
  * fhat and q_tau are also stored per replicate, so the KDE density estimate
    -- whose bandwidth is 1.06 sd(Y) N^-0.2, the one explicitly N-dependent
    piece of the RIF -- can be compared to the exact population density.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

TAUS = np.linspace(0.02, 0.98, 45)
T = len(TAUS)
GAMMA_SHAPE = {"gamma2": 2.0, "gamma8": 8.0}


def draw_eps(rng, n, shape):
    if shape == "normal":
        return rng.normal(size=n)
    k = GAMMA_SHAPE[shape]
    x = rng.gamma(k, 1.0, size=n)
    return (x - k) / np.sqrt(k)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ws", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--cell", type=int, required=True)
    ap.add_argument("--chunk", type=int, required=True)
    ap.add_argument("--reps", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sys.path.insert(0, args.tree)
    from fungwas_stage1 import core, pooled_cov

    cells = json.load(open(os.path.join(args.ws, "cells.json")))
    cell = cells[args.cell]
    N = int(cell["N"])
    nb = int(cell.get("n_blocks", 1000))
    rep0 = args.chunk * args.reps

    tag = "cell%02d_chunk%04d" % (args.cell, args.chunk)
    if os.path.exists(os.path.join(args.out, tag + ".json")):
        print(tag, "already done", flush=True)
        return

    B = np.empty((args.reps, T))
    D = np.empty(args.reps)
    S0 = np.empty((args.reps, T, T))
    MAC = np.empty(args.reps)
    FH = np.empty((args.reps, T))
    QT = np.empty((args.reps, T))

    X = np.ones((N, 1))
    Qb, _ = np.linalg.qr(X)
    t0 = time.time()
    for r in range(args.reps):
        # seed stream distinct from the gls_coverage study (different 3rd entry)
        rng = np.random.default_rng([args.cell, rep0 + r, 20260915])
        g = rng.binomial(2, cell["maf"], size=N).astype(np.float64)
        eps = draw_eps(rng, N, cell["shape"])
        y = cell["a"] * g + (1.0 + cell["c"] * g) * eps

        res = core.run_stage1(g.reshape(-1, 1), y, TAUS, n_blocks=nb)
        B[r] = res["betas"][0]
        D[r] = res["d_resid"][0]
        MAC[r] = g.sum()

        rif, qt, fh = core.compute_rif(y, TAUS)
        QT[r] = qt
        FH[r] = fh
        rif_resid = rif - Qb @ (Qb.T @ rif)
        S0[r] = pooled_cov.compute_sigma0(rif_resid)
        del rif, rif_resid

    os.makedirs(args.out, exist_ok=True)
    B.astype("<f8").tofile(os.path.join(args.out, tag + ".beta.bin"))
    D.astype("<f8").tofile(os.path.join(args.out, tag + ".dresid.bin"))
    S0.astype("<f8").tofile(os.path.join(args.out, tag + ".sigma0.bin"))
    MAC.astype("<f8").tofile(os.path.join(args.out, tag + ".mac.bin"))
    FH.astype("<f8").tofile(os.path.join(args.out, tag + ".fhat.bin"))
    QT.astype("<f8").tofile(os.path.join(args.out, tag + ".qtau.bin"))
    with open(os.path.join(args.out, tag + ".json"), "w") as fh_:
        json.dump(dict(cell=cell, cell_index=args.cell, chunk=args.chunk,
                       reps=args.reps, rep0=rep0, T=T, n_blocks=nb,
                       secs=time.time() - t0), fh_)
    print("%s: %d reps in %.1fs" % (tag, args.reps, time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
