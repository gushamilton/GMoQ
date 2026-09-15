#!/usr/bin/env python3
"""Bounded covariate-basis probe for the shipped Stage 1 path.

This intentionally reuses the pilot's covariate DGP and seeds.  It compares
the production linear C basis with C,C^2,C^3; it is a diagnostic, not a
preprocessing change.
"""
from __future__ import annotations
import argparse, hashlib, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pilot_stage1 as p  # noqa: E402


def basis(C, kind):
    c = C[:, 0]
    if kind == "linear":
        return c[:, None]
    return np.column_stack([c, c * c, c * c * c])


def run_arm(G, y, C, kind, path, seed, blocks, threads):
    from fungwas_stage1 import core, pooled_cov
    Z = basis(C, kind)
    X = np.column_stack([np.ones(len(y)), Z])
    if path == "current":
        res = core.run_stage1(G, y, p.TAUS, covariates=Z, n_blocks=blocks,
                              seed=seed, n_threads=threads)
        rif, q, f = core.compute_rif(y, p.TAUS)
        rif_r, _ = core.residualize(rif, X)
        ymean, ysd = float(np.mean(y)), float(np.std(y, ddof=1))
    else:
        y_r, _ = core.residualize(y, X)
        rif, q, f = core.compute_rif(y_r, p.TAUS)
        np.random.seed(seed)
        bids = np.random.randint(0, blocks, len(y)).astype(np.int32)
        res = p.stage1_from_rif(G, rif, X, bids, blocks, threads)
        rif_r, _ = core.residualize(rif, X)
        ymean, ysd = float(np.mean(y_r)), float(np.std(y_r, ddof=1))
    res["q_tau"], res["fhat"] = q, f
    res["sigma0"] = pooled_cov.compute_sigma0(rif_r)
    res["y_mean"], res["y_sd"] = ymean, ysd
    return res


def projection(cell, kind, path, Nref, seed):
    from fungwas_stage1 import core
    rng = np.random.default_rng(seed)
    G, y, C = p.draw_covariate_data(cell, Nref, rng)
    Z = basis(C, kind); X = np.column_stack([np.ones(Nref), Z])
    if path == "y_first":
        y, _ = core.residualize(y, X)
    rif, _, _ = core.compute_rif(y, p.TAUS)
    rr, Q = core.residualize(rif, X)
    gr = G - Q @ (Q.T @ G)
    return np.sum(gr[:, :, None] * rr[:, None, :], axis=0) / np.sum(gr * gr, axis=0)[:, None]


def put(path, arr): np.asarray(arr, dtype="<f8").tofile(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--repo", required=True)
    ap.add_argument("--reps", type=int, default=100); ap.add_argument("--N", type=int, default=10000)
    ap.add_argument("--blocks", type=int, default=1000); ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--Nref", type=int, default=200000); ap.add_argument("--seed", type=int, default=9142026)
    ap.add_argument("--source-commit", default="")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    p.gs.add_stage1_to_path(a.repo)
    from fungwas_stage1 import core
    if not core.HAVE_CPP: raise RuntimeError("probe requires shipped Stage 1 C++")
    stage = os.path.join(a.out, "stage1"); truth = os.path.join(a.out, "truth")
    os.makedirs(stage, exist_ok=True); os.makedirs(truth, exist_ok=True)
    cells = [p.cells()[4], p.cells()[9]]
    kinds = ["linear", "rich"]
    script_hash = hashlib.sha256(open(__file__, "rb").read()).hexdigest()
    meta = dict(source_commit=a.source_commit, script_sha256=script_hash, N=a.N, reps=a.reps, blocks=a.blocks,
                taus=p.TAUS.tolist(), M=p.M, Nref=a.Nref,
                basis="linear C versus rich C,C^2,C^3", mc_unit="dataset")
    json.dump(meta, open(os.path.join(a.out, "probe_meta.json"), "w"), indent=2)
    triu = np.triu_indices(len(p.TAUS)); NE = len(triu[0])
    for ci, cell in enumerate(cells):
        for kind in kinds:
            proj = {path: projection(cell, kind, path, a.Nref, a.seed + 50000 + 100 * ci + 10 * (kind == "rich") + (path == "y_first"))
                    for path in ("current", "y_first")}
            json.dump(dict(cell=cell, kind=kind, projection={k: v.tolist() for k, v in proj.items()},
                           scientific_a=[0.0] * p.M, scientific_b=[0.0] * p.M,
                           taus=p.TAUS.tolist()), open(os.path.join(truth, f"{cell['name']}__{kind}.json"), "w"), indent=2)
            for path in ("current", "y_first"):
                vals = {k: [] for k in ("beta", "se", "cov", "dresid", "sigma0", "q_tau", "y_mean", "y_sd")}
                for r in range(a.reps):
                    seed = a.seed + 100000 * (4 + ci * 5) + r
                    rng = np.random.default_rng(seed)
                    G, y, C = p.draw_covariate_data(cell, a.N, rng)
                    res = run_arm(G, y, C, kind, path, seed, a.blocks, a.threads)
                    vals["beta"].append(res["betas"]); vals["se"].append(res["se"])
                    vals["cov"].append(np.stack([np.asarray(c)[triu] for c in res["covariances"]]))
                    vals["dresid"].append(res["d_resid"]); vals["sigma0"].append(res["sigma0"])
                    vals["q_tau"].append(res["q_tau"]); vals["y_mean"].append(res["y_mean"]); vals["y_sd"].append(res["y_sd"])
                prefix = os.path.join(stage, f"{cell['name']}__{kind}__{path}")
                for k, v in vals.items(): put(prefix + f".{k}.bin", np.asarray(v))
                json.dump(dict(reps=a.reps, N=a.N, blocks=a.blocks, M=p.M, T=len(p.TAUS), ne=NE),
                          open(prefix + ".json", "w"), indent=2)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__": main()
