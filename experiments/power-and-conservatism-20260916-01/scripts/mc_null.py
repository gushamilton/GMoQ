"""Monte Carlo null for the G-MoQ pooled/OLS Stage 2 test at the production N.

Two engines:
  plain  -- iid HWE genotypes, full pipeline (D from the realised genotype, residualised
            on the same basis Q Stage 1 uses). Gives the rejection rate at moderate alpha
            and the variance of t with sqrt(M) error.
  is     -- exponentially tilted importance sampling at the saddlepoint for a given
            threshold. Gives a direct, unbiased MC estimate of the rejection rate at
            5e-8 with an honest standard error, so no extrapolation is needed.

Usage: mc_null.py plain <task> <ntask> <m_per_task>
       mc_null.py is    <task> <ntask> <m_per_task>
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
from scipy.stats import norm, t as tdist
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nullcal import cgf, spa_two_sided, spa_upper, hwe, _solve_saddle

S = "/user/work/fh6520/gmoq_power_conserv"
OUT = f"{S}/out"; MCD = f"{S}/mc"; os.makedirs(MCD, exist_ok=True)

MODE, TASK, NTASK, MPER = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(20260916 + 7919 * TASK)

Psi = np.load(f"{S}/psi/psi_chr22.npy"); N, T = Psi.shape
Q = np.load(f"{S}/psi/Q_chr22.npy")
A = pd.read_csv(f"{OUT}/A_contrasts.tsv", sep="\t", index_col=0)
params = list(A.index)
U = np.asfortranarray(Psi @ A.values.astype(np.float64).T)     # N x 3
sd_u = np.sqrt((U ** 2).sum(axis=0) / N)
su2 = (U ** 2).sum(axis=0)
del Psi
Uf = np.ascontiguousarray(U.T, dtype=np.float64)             # 3 x N
Qf = np.ascontiguousarray(Q.T, dtype=np.float64)             # K x N
K = Qf.shape[0]

bins = pd.read_csv(f"{OUT}/spa_bins.tsv", sep="\t")
bins = bins[(bins.param == params[0]) & (bins.alpha == 0.05)][["bin", "n", "maf"]]
FS = bins.maf.values; WS = bins.n.values / bins.n.sum()

ALPHAS = np.array([0.05, 1e-3, 1e-5, 5e-8])
CRIT = tdist.ppf(1 - ALPHAS / 2, df=N)

CHUNK = 100


def stats_for(G):
    """t statistic (3 x m) for a genotype block G (N x m, float64)."""
    Sv = Uf @ G                                  # 3 x m
    gg = np.einsum("ij,ij->j", G, G)
    QG = Qf @ G                                  # K x m
    D = gg - np.einsum("ij,ij->j", QG, QG)
    return Sv / (np.sqrt(D)[None, :] * sd_u[:, None]), D


if MODE == "plain":
    cnt = np.zeros((3, len(ALPHAS)))
    s1 = np.zeros(3); s2 = np.zeros(3); s4 = np.zeros(3); ntot = 0
    zs = []          # keep a subsample of |t| for the qq/lambda
    per_maf = {}
    for rep in range(0, MPER, CHUNK):
        m = min(CHUNK, MPER - rep)
        f = rng.choice(FS, size=1, p=WS)[0]
        G = rng.binomial(2, f, size=(N, m)).astype(np.float64)
        tt, D = stats_for(G)
        for a_i, c in enumerate(CRIT):
            cnt[:, a_i] += (np.abs(tt) > c).sum(axis=1)
        s1 += tt.sum(axis=1); s2 += (tt ** 2).sum(axis=1); s4 += (tt ** 4).sum(axis=1)
        ntot += m
        k = round(float(f), 4)
        e = per_maf.setdefault(k, [0, np.zeros(3), np.zeros((3, len(ALPHAS)))])
        e[0] += m; e[1] += (tt ** 2).sum(axis=1)
        for a_i, c in enumerate(CRIT):
            e[2][:, a_i] += (np.abs(tt) > c).sum(axis=1)
        if len(zs) < 40:
            zs.append(tt[:, :50].copy())
    res = dict(mode="plain", task=TASK, n=ntot, params=params,
               alphas=ALPHAS.tolist(), crit=CRIT.tolist(),
               counts=cnt.tolist(), sum_t=s1.tolist(), sum_t2=s2.tolist(),
               sum_t4=s4.tolist(),
               per_maf={str(k): [v[0], v[1].tolist(), v[2].tolist()] for k, v in per_maf.items()})
    with open(f"{MCD}/plain_{TASK:03d}.json", "w") as fh:
        json.dump(res, fh)
    print("plain done", ntot, cnt.sum(axis=1), flush=True)

elif MODE == "is":
    # cells: (param, maf, alpha) -- the far tail only
    CELLS = [(pi, float(f), float(al))
             for pi in range(len(params))
             for f in (FS[0], FS[15])
             for al in (5e-8, 1e-5)]
    out = []
    for ci, (pi, f, al) in enumerate(CELLS):
        if ci % NTASK != TASK:
            continue
        c = float(tdist.ppf(1 - al / 2, df=N))
        u = np.ascontiguousarray(U[:, pi])
        vals, probs = hwe(f)
        var_s = 2 * f * (1 - f) * su2[pi]
        Dexp = 2 * f * (1 - f) * (N - K)
        thr = c * np.sqrt(Dexp) * sd_u[pi]
        tot = 0.0; totvar = 0.0
        for sgn in (+1.0, -1.0):
            uu = sgn * u
            s, Kv, k2 = _solve_saddle(thr, uu, vals, probs, thr / var_s)
            if not np.isfinite(s):
                continue
            z = s * np.outer(uu, vals)
            w = probs[None, :] * np.exp(z - z.max(axis=1, keepdims=True))
            w /= w.sum(axis=1, keepdims=True)
            cum = np.ascontiguousarray(np.cumsum(w, axis=1)[:, :2])
            acc = 0.0; acc2 = 0.0
            step = 100
            for r0 in range(0, MPER, step):
                mm = min(step, MPER - r0)
                Uu = rng.random((N, mm), dtype=np.float32)
                G = (Uu > cum[:, 0:1].astype(np.float32)).astype(np.float64)
                G += (Uu > cum[:, 1:2].astype(np.float32))
                Sv = uu @ G
                ind = np.where(Sv >= thr, np.exp(-s * Sv + Kv), 0.0)
                acc += ind.sum(); acc2 += (ind ** 2).sum()
            mean = acc / MPER
            tot += mean
            totvar += (acc2 / MPER - mean ** 2) / MPER
        p_spa, _, _ = spa_two_sided(thr, u, vals, probs, var_s)
        p_norm = 2 * norm.sf(thr / np.sqrt(var_s))
        out.append(dict(param=params[pi], maf=f, alpha=al, thr=float(thr),
                        rate_is=tot, se_is=float(np.sqrt(max(totvar, 0))),
                        rate_spa=p_spa, rate_normal=p_norm, m=MPER))
        print(out[-1], flush=True)
    pd.DataFrame(out).to_csv(f"{MCD}/is_{TASK:03d}.tsv", sep="\t", index=False)
    print("is done", flush=True)
