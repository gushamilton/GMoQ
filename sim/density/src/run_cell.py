#!/usr/bin/env python
"""Drive the SHIPPED Stage 1 over replicates of one cell; dump raw results."""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shapes import Base, Dmap, Truth, make_pheno
from cells import CELLS, TAUS, build

from fungwas_stage1.core import compute_rif, residualize, run_stage1
from fungwas_stage1.pooled_cov import compute_sigma0
try:
    from fungwas_stage1.core import _bandwidth as _bw
except ImportError:
    def _bw(Y, N):
        return 1.06 * np.std(Y) * N ** (-0.2)

MAF_CAUSAL = 0.20
MAF_NULL = [0.05, 0.20, 0.50]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--s0", type=float, default=0.06)
    ap.add_argument("--N", type=int, default=20000)
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--rep0", type=int, default=0)
    ap.add_argument("--n-null", type=int, default=201)
    ap.add_argument("--blocks", type=int, default=1000)
    ap.add_argument("--n-covar", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--bw", default=None, help="FUNGWAS_BW rule for the patched tree")
    a = ap.parse_args()
    if a.bw:
        os.environ["FUNGWAS_BW"] = a.bw

    base, D = build(a.cell)
    taus = TAUS
    T = len(taus)
    nel = T * (T + 1) // 2
    iu = np.triu_indices(T)

    # -------- population truth (once)
    tr_alt = Truth(base, D, MAF_CAUSAL, 0.0, a.s0, taus, a.N)
    # b_latent depends only on (base, maf, a, s), never on D
    # W: column 1 = pure shift (flat), column 2 = exact latent scale shape / s0
    u = tr_alt.b_latent / a.s0
    W = np.column_stack([np.ones(T), u])
    theta_true = np.array([0.0, a.s0])

    n_null = a.n_null
    mafs = np.array([MAF_CAUSAL] + [MAF_NULL[j % len(MAF_NULL)] for j in range(n_null)])
    M = len(mafs)

    rng_master = np.random.default_rng(a.seed * 100003 + a.rep0)

    B = np.empty((a.reps, M, T), dtype=np.float64)
    S = np.empty((a.reps, M, T), dtype=np.float64)
    C = np.empty((a.reps, M, nel), dtype=np.float32)
    DR = np.empty((a.reps, M), dtype=np.float64)
    SIG = np.empty((a.reps, nel), dtype=np.float64)
    QT = np.empty((a.reps, T), dtype=np.float64)
    FH = np.empty((a.reps, T), dtype=np.float64)
    # diagnostics computable without truth
    TIE = np.empty((a.reps, T), dtype=np.float64)      # sample fraction exactly == q_tau
    FGAP = np.empty((a.reps, T), dtype=np.float64)     # |F_n(q_tau) - tau|
    SPC = np.empty((a.reps, T), dtype=np.float64)      # local spacing of distinct values
    NDH = np.empty((a.reps, T), dtype=np.float64)      # distinct values within +-h of q_tau
    HB = np.empty(a.reps, dtype=np.float64)
    NDIST = np.empty(a.reps, dtype=np.int64)
    MAXATOM = np.empty(a.reps, dtype=np.float64)

    t0 = time.time()
    for r in range(a.reps):
        rng = np.random.default_rng(rng_master.integers(1 << 62))
        G = np.empty((a.N, M))
        for j, p in enumerate(mafs):
            G[:, j] = rng.binomial(2, p, size=a.N).astype(float)
        Y, _ = make_pheno(base, D, G[:, 0], 0.0, a.s0, a.N, rng)
        cov = rng.normal(size=(a.N, a.n_covar)) if a.n_covar else None

        res = run_stage1(G, Y, taus, covariates=cov, n_blocks=a.blocks,
                         seed=int(rng.integers(1 << 30)), n_threads=a.threads)
        B[r] = res["betas"]; S[r] = res["se"]; DR[r] = res["d_resid"]
        QT[r] = res["q_tau"]
        for j in range(M):
            C[r, j] = np.asarray(res["covariances"][j])[iu].astype(np.float32)

        # re-derive RIF / Sigma0 with the SAME shipped functions
        RIF, q_tau, fhat = compute_rif(Y, taus)
        X = np.ones((a.N, 1)) if cov is None else np.column_stack([np.ones(a.N), cov])
        RIF_resid, _ = residualize(RIF, X)
        SIG[r] = compute_sigma0(RIF_resid)[iu]
        FH[r] = fhat
        assert np.allclose(q_tau, res["q_tau"])

        Ys = np.sort(Y)
        vals, cnts = np.unique(Y, return_counts=True)
        NDIST[r] = len(vals); MAXATOM[r] = cnts.max() / a.N
        HB[r] = _bw(Y, a.N)
        for t in range(T):
            q = q_tau[t]
            TIE[r, t] = np.count_nonzero(Y == q) / a.N
            FGAP[r, t] = np.count_nonzero(Y <= q) / a.N - taus[t]
            i = np.searchsorted(vals, q)
            lo = max(i - 1, 0); hi = min(i + 1, len(vals) - 1)
            SPC[r, t] = (vals[hi] - vals[lo]) / max(hi - lo, 1) if len(vals) > 1 else np.nan
            NDH[r, t] = np.count_nonzero((vals >= q - HB[r]) & (vals <= q + HB[r]))

    el = time.time() - t0
    # recompute the shipped estimand at the bandwidth actually realised
    h_real = float(np.mean(HB))
    tr_real = Truth(base, D, MAF_CAUSAL, 0.0, a.s0, taus, a.N, h_override=h_real)
    os.makedirs(a.out, exist_ok=True)
    np.savez_compressed(os.path.join(a.out, f"raw_{a.rep0:05d}.npz"),
                        beta=B, se=S, cov=C, d_resid=DR, sigma0=SIG, q_tau=QT,
                        fhat=FH, tie=TIE, fgap=FGAP, spacing=SPC, hbw=HB, ndist_h=NDH,
                        ndistinct=NDIST, maxatom=MAXATOM, mafs=mafs)
    meta = dict(cell=a.cell, s0=a.s0, N=a.N, reps=a.reps, rep0=a.rep0,
                n_null=n_null, blocks=a.blocks, n_covar=a.n_covar, taus=taus.tolist(),
                W=W.tolist(), theta_true=theta_true.tolist(),
                b_latent=tr_alt.b_latent.tolist(),
                b_obs=np.where(np.isfinite(tr_alt.b_obs), tr_alt.b_obs, None).tolist(),
                b_ship=tr_real.b_ship.tolist(),
                b_ship_silverman=tr_alt.b_ship.tolist(),
                h_realised=h_real,
                f_latent=tr_alt.f_latent.tolist(),
                f_obs_true=np.where(np.isfinite(tr_alt.f_obs_true), tr_alt.f_obs_true, None).tolist(),
                f_ship=tr_real.f_ship.tolist(),
                f_ship_silverman=tr_alt.f_ship.tolist(),
                q_obs=tr_alt.q_obs.tolist(), q_latent=tr_alt.q_latent.tolist(),
                atom=tr_alt.atom.tolist(), h_pop=tr_alt.h, sd_obs=tr_alt.sd_obs,
                elapsed=el, mafs=mafs.tolist(), bw_rule=(a.bw or "silverman"))
    with open(os.path.join(a.out, f"meta_{a.rep0:05d}.json"), "w") as f:
        json.dump(meta, f)
    print(f"{a.cell} rep0={a.rep0} reps={a.reps} M={M} elapsed={el:.1f}s", flush=True)


if __name__ == "__main__":
    main()
