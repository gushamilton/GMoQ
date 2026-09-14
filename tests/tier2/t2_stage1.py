#!/usr/bin/env python3
"""One tier 2 array task: the exact truth for a cell, and Stage 1 for one chunk.

    python tests/tier2/t2_stage1.py --cells cells.json --task <i> --out <dir>

The design is identical across replicates within a cell (the MAFs and the
effect sizes are deterministic functions of the cell spec), so the exact
population estimand b(tau) is computed ONCE per cell and cached. What varies
between replicates is the genotype draw and the phenotype draw.

Output, in the flat-binary layout the existing sim drivers already use so the R
side reads it with readBin and nothing new has to be invented:

    <tag>.beta.bin    (reps, M, T)  float64
    <tag>.se.bin      (reps, M, T)  float64
    <tag>.cov.bin     (reps, M, T(T+1)/2) float32   packed upper triangles
    <tag>.dresid.bin  (reps, M)     float64
    <tag>.sigma0.bin  (reps, T, T)  float64
    <tag>.json        shapes, the cell spec, and timings
    <cell>.btrue.bin  (M, T)        float64        written once per cell
    <cell>.truth.json q_pop, f_pop, theta_latent, maf, causal indices
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import gmoq_sim as gs  # noqa: E402


def cell_taus(c):
    if c["n_taus"] == 45:
        return gs.PROD_TAUS
    return np.round(np.linspace(0.02, 0.98, c["n_taus"]), 3)


def make_design(c, taus):
    rng = np.random.default_rng(0)
    return gs.build_design(
        c["family"], c["maf"], c["n_causal"], c["z_target"], c["N"], taus, rng,
        baseline_skew=c["baseline_skew"], baseline_kurt=c["baseline_kurt"],
        n_null=c["n_null"], discretise=c["discretise"])


def ensure_truth(c, taus, design, out_dir, nbg):
    """Compute (and cache) the exact population estimand for this cell."""
    bt = os.path.join(out_dir, f"{c['name']}.btrue.bin")
    tj = os.path.join(out_dir, f"{c['name']}.truth.json")
    if os.path.exists(tj) and os.path.exists(bt):
        return
    M = len(design.maf)
    b_true = np.zeros((M, len(taus)))
    q_pop = np.full(len(taus), np.nan)
    f_pop = np.full(len(taus), np.nan)
    if c["discretise"] is None:
        tr = gs.Truth(design, nbg=nbg, seed=12345)
        q_pop, f_pop, b_c = tr.b_true(taus)
        if len(design.causal):
            b_true[design.causal] = b_c
    # A discretised phenotype changes the estimand as well as the estimator
    # (docs/BP2 section 4), and the pins that use one are all null cells, so no
    # truth is claimed for them: b_true stays zero, which IS the truth there.
    b_true.astype("<f8").tofile(bt)
    with open(tj, "w") as fh:
        json.dump(dict(
            name=c["name"], taus=[float(t) for t in taus],
            q_pop=[float(x) for x in q_pop], f_pop=[float(x) for x in f_pop],
            maf=[float(x) for x in design.maf],
            causal=[int(x) + 1 for x in design.causal],
            n_variants=int(M),
            theta_latent={k: [float(x) for x in getattr(design, k)]
                          for k in ("a", "b", "c", "d", "e", "m")},
            mix=design.mix, discretise=c["discretise"],
            truth_is_zero=bool(c["discretise"] is not None)), fh, indent=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True)
    ap.add_argument("--task", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--nbg", type=int, default=20000)
    ap.add_argument("--seed0", type=int, default=7_000_000)
    a = ap.parse_args()

    spec = json.load(open(a.cells))
    task = spec["tasks"][a.task]
    c = spec["cells"][task["cell"]]
    chunk = task["chunk"]

    os.makedirs(a.out, exist_ok=True)
    tag = f"{c['name']}__chunk{chunk:03d}"
    done = os.path.join(a.out, tag + ".json")
    if os.path.exists(done):
        print(tag, "already done", flush=True)
        return

    taus = cell_taus(c)
    T = len(taus)
    design = make_design(c, taus)
    M = len(design.maf)

    t_truth0 = time.time()
    ensure_truth(c, taus, design, a.out, a.nbg)
    t_truth = time.time() - t_truth0

    gs.add_stage1_to_path()
    from fungwas_stage1 import core, pooled_cov

    n_rep = min(c["chunk_reps"], c["reps"] - chunk * c["chunk_reps"])
    n_elem = T * (T + 1) // 2
    triu = np.triu_indices(T)

    beta = np.zeros((n_rep, M, T))
    se = np.zeros((n_rep, M, T))
    cov = np.zeros((n_rep, M, n_elem), dtype=np.float32)
    dres = np.zeros((n_rep, M))
    s0m = np.zeros((n_rep, T, T))
    ysd = np.zeros(n_rep)

    t0 = time.time()
    for r in range(n_rep):
        seed = a.seed0 + 1000 * task["cell"] + 10 * chunk + r
        rng = np.random.default_rng(seed)
        G = rng.binomial(2, np.broadcast_to(design.maf, (c["N"], M))).astype(np.float64)
        y = gs.simulate_phenotype(design, G, rng)
        C = rng.normal(size=(c["N"], c["n_covar"])) if c["n_covar"] else None

        res = core.run_stage1(G, y, taus, covariates=C, n_blocks=c["blocks"],
                              seed=seed, n_threads=a.threads)
        beta[r] = res["betas"]; se[r] = res["se"]; dres[r] = res["d_resid"]
        for j in range(M):
            cov[r, j] = np.asarray(res["covariances"][j])[triu].astype(np.float32)

        rif, _, _ = core.compute_rif(y, taus)
        X = np.ones((c["N"], 1)) if C is None else np.column_stack([np.ones(c["N"]), C])
        rif_r, _ = core.residualize(rif, X)
        s0m[r] = pooled_cov.compute_sigma0(rif_r)
        ysd[r] = float(y.std(ddof=1))
        print(f"  rep {r + 1}/{n_rep} {time.time() - t0:.1f}s", flush=True)

    p = os.path.join(a.out, tag)
    beta.astype("<f8").tofile(p + ".beta.bin")
    se.astype("<f8").tofile(p + ".se.bin")
    cov.astype("<f4").tofile(p + ".cov.bin")
    dres.astype("<f8").tofile(p + ".dresid.bin")
    s0m.astype("<f8").tofile(p + ".sigma0.bin")
    with open(done, "w") as fh:
        json.dump(dict(cell=c, chunk=chunk, task=a.task, reps=int(n_rep),
                       M=int(M), T=int(T), n_elem=int(n_elem),
                       taus=[float(t) for t in taus],
                       y_sd=[float(v) for v in ysd],
                       stage1_seconds=round(time.time() - t0, 2),
                       truth_seconds=round(t_truth, 2)), fh)
    print(f"{tag}: {n_rep} reps in {time.time() - t0:.1f}s "
          f"(truth {t_truth:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
