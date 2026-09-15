#!/usr/bin/env python3
"""Paired density/RIF diagnosis for the Normal variance-linear pilot.

This diagnostic keeps the shipped Stage 1 block-score and covariance kernels
in the loop while replacing only the RIF inputs.  For every frozen dataset it
constructs three paired arms:

  A: empirical type-8 quantile plus the shipped Silverman KDE;
  B: the same empirical quantile plus exact mixture density at that quantile;
  C: exact population quantile and exact mixture density.

The output layout intentionally matches pilot_stage1.py, so the existing
pilot_stage2.R wrapper can fit the installed Stage 2 package on every arm.
This file is a diagnostic harness; it does not change production code.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time

import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import norm


TAUS = np.round(np.linspace(0.02, 0.98, 45), 3)
Z_TAU = np.asarray(norm.ppf(TAUS), dtype=float)
M = 1
ARMS = ("A_empirical_kde", "B_empirical_exact_f", "C_population_exact")


def git_value(repo: str, *args: str) -> str:
    try:
        return subprocess.check_output(["git", "-C", repo, *args], text=True).strip()
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{exc}"


def write_array(path: str, arr: np.ndarray, dtype: str = "<f8") -> None:
    np.asarray(arr, dtype=dtype).tofile(path)


def cells() -> list[dict]:
    return [
        {"name": "density_v0", "family": "normal_variance", "maf": 0.30, "a": 0.0, "v": 0.0},
        {"name": "density_v06", "family": "normal_variance", "maf": 0.30, "a": 0.0, "v": 0.06},
        {"name": "density_v24", "family": "normal_variance", "maf": 0.30, "a": 0.0, "v": 0.24},
    ]


def mixture_components(cell: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    p = float(cell["maf"])
    a = float(cell["a"])
    v = float(cell["v"])
    g = np.arange(3.0)
    x = g - 2.0 * p
    pg = np.array([(1.0 - p) ** 2, 2.0 * p * (1.0 - p), p ** 2])
    sd = np.sqrt(1.0 + v * x)
    return x, pg, sd


def exact_qf(cell: dict, taus: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Exact q_tau and f(q_tau) for the three-component Normal mixture."""
    p = float(cell["maf"])
    a = float(cell["a"])
    x, pg, sd = mixture_components(cell)

    def cdf(y: float) -> float:
        return float(np.sum(pg * ndtr((y - a * x) / sd)))

    def pdf(y: float) -> float:
        z = (y - a * x) / sd
        return float(np.sum(pg * np.exp(-0.5 * z * z) /
                           (np.sqrt(2.0 * np.pi) * sd)))

    q = np.array([brentq(lambda yy: cdf(yy) - float(t), -12.0, 12.0)
                  for t in np.asarray(taus, dtype=float)])
    f = np.maximum(np.array([pdf(yy) for yy in q]), 1e-300)
    return q, f


def exact_density(cell: dict, y: np.ndarray) -> np.ndarray:
    """Exact marginal mixture density evaluated at arbitrary y values."""
    p = float(cell["maf"])
    a = float(cell["a"])
    x, pg, sd = mixture_components(cell)
    yy = np.asarray(y, dtype=float)[:, None]
    z = (yy - a * x[None, :]) / sd[None, :]
    return np.sum(pg[None, :] * np.exp(-0.5 * z * z) /
                  (np.sqrt(2.0 * np.pi) * sd[None, :]), axis=1)


def exact_projection(cell: dict, taus: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return exact q, f, and the tau-level population RIF slope."""
    p = float(cell["maf"])
    a = float(cell["a"])
    x, pg, sd = mixture_components(cell)
    q, f = exact_qf(cell, taus)
    cond_cdf = np.asarray([ndtr((yy - a * x) / sd) for yy in q])
    vg = 2.0 * p * (1.0 - p)
    cov = np.sum(pg[None, :] * x[None, :] * cond_cdf, axis=1)
    slope = -cov / (f * vg)
    return q, f, slope


def simulate(cell: dict, G: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    p = float(cell["maf"])
    x = G[:, 0] - 2.0 * p
    variance = 1.0 + float(cell["v"]) * x
    if np.min(variance) <= 0:
        raise ValueError("non-positive conditional variance")
    return float(cell["a"]) * x + np.sqrt(variance) * rng.normal(size=len(x))


def stage1_from_rif(core, G: np.ndarray, rif: np.ndarray, block_ids: np.ndarray,
                    n_blocks: int, threads: int) -> dict:
    """Run the shipped residualisation, C++ block scores, and JK covariance."""
    X = np.ones((G.shape[0], 1), dtype=float)
    rif_resid, Q = core.residualize(rif, X)
    stats = core.compute_block_scores(
        G, rif_resid, Q, block_ids, n_blocks, n_threads=threads,
    )
    beta, se, cov = core.compute_covariance_from_blocks(stats[0], n_blocks)
    d_resid = float(np.sum(stats[0][:, 0]))
    from fungwas_stage1 import pooled_cov
    sigma0 = pooled_cov.compute_sigma0(rif_resid)
    return dict(beta=beta, se=se, cov=cov, dresid=d_resid,
                sigma0=sigma0, rif_resid=rif_resid)


def make_rif(core, arm: str, y: np.ndarray, cell: dict,
             taus: np.ndarray, q_pop: np.ndarray, f_pop: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Construct one diagnostic RIF arm and return RIF, q, f_used, f_ref."""
    if arm == "A_empirical_kde":
        rif, q, fhat = core.compute_rif(y, taus)
        return rif, q, fhat, fhat / exact_density(cell, q)
    _, q_emp, fhat = core.compute_rif(y, taus)
    if arm == "B_empirical_exact_f":
        f_used = exact_density(cell, q_emp)
        rif = q_emp[None, :] + (taus[None, :] - (y[:, None] <= q_emp[None, :])) / f_used[None, :]
        return rif, q_emp, f_used, fhat / f_used
    if arm == "C_population_exact":
        rif = q_pop[None, :] + (taus[None, :] - (y[:, None] <= q_pop[None, :])) / f_pop[None, :]
        return rif, q_pop, f_pop, fhat / f_pop
    raise ValueError(f"unknown arm: {arm}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--N", type=int, default=10000)
    ap.add_argument("--blocks", type=int, default=1000)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--seed", type=int, default=9142026)
    ap.add_argument("--source-commit", default="")
    ap.add_argument("--source-status", default="")
    ap.add_argument("--run-id", default="")
    ap.add_argument("--cell", choices=("v0", "v06", "v24"), default="")
    args = ap.parse_args()

    repo = os.path.abspath(args.repo)
    sys.path.insert(0, os.path.join(repo, "stage1"))
    from fungwas_stage1 import core
    core.require_cpp_extension("density probe")

    out = os.path.abspath(args.out)
    run_id = args.run_id or os.path.basename(out.rstrip(os.sep))
    if os.path.exists(out) and os.listdir(out):
        raise SystemExit(f"refusing non-empty output directory: {out}")
    stage_dir = os.path.join(out, "stage1")
    truth_dir = os.path.join(out, "truth")
    os.makedirs(stage_dir, exist_ok=True)
    os.makedirs(truth_dir, exist_ok=True)
    source_commit = args.source_commit or git_value(repo, "rev-parse", "HEAD")
    source_status = args.source_status or git_value(repo, "status", "--short", "--branch")
    script_sha256 = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()
    meta = dict(
        design="paired_density_rif_diagnosis",
        source_commit=source_commit,
        source_status=source_status,
        script_sha256=script_sha256,
        repo_head=git_value(repo, "rev-parse", "HEAD"),
        python=sys.version,
        numpy=np.__version__,
        scipy=__import__("scipy").__version__,
        host=platform.node(),
        N=args.N,
        reps=args.reps,
        blocks=args.blocks,
        threads=args.threads,
        taus=TAUS.tolist(),
        M=M,
        arms=list(ARMS),
        seed=args.seed,
        run_id=run_id,
        selected_cell=args.cell or "all",
        mc_unit="dataset (one shared phenotype/genotype draw; arms are paired)",
        target_definition=("projection is exact population RIF slope projected onto W; "
                           "scientific target is (a,v) for Var(Y|g)=1+v(g-2p)"),
    )
    with open(os.path.join(out, "density_probe_meta.json"), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)

    run_cells = cells()
    if args.cell:
        run_cells = [cell for cell in run_cells if cell["name"] == f"density_{args.cell}"]
    W = np.column_stack([np.ones(len(TAUS)), Z_TAU / 2.0])
    A = np.linalg.solve(W.T @ W, W.T)
    aggregate = []
    for ci, cell in enumerate(run_cells):
        q_pop, f_pop, b_tau = exact_projection(cell, TAUS)
        target_proj = A @ b_tau
        target_sci = np.array([float(cell["a"]), float(cell["v"])])
        for arm in ARMS:
            name = f"{cell['name']}_{arm}"
            arm_cell = dict(cell, name=name, arm=arm)
            truth = dict(
                cell=arm_cell,
                maf=[float(cell["maf"])],
                scientific_parameterisation="normal_variance_linear",
                scientific_a=[float(cell["a"])],
                scientific_logsd_b=[0.0],
                scientific_variance_coefficient=[float(cell["v"])],
                scientific_variance_ratio=[1.0 + float(cell["v"]) * (1.0 - 2.0 * float(cell["maf"]))],
                projection={"current": [b_tau.tolist()]},
                projection_truth_kind="normal_linear_exact",
                projection_reference_N=None,
                taus=TAUS.tolist(),
                diagnostic_arm=arm,
                target_projected_parameters=target_proj.tolist(),
            )
            with open(os.path.join(truth_dir, name + ".json"), "w", encoding="utf-8") as fh:
                json.dump(truth, fh, indent=2)

            beta_rows = []
            se_rows = []
            cov_rows = []
            d_rows = []
            sigma_rows = []
            q_rows = []
            f_rows = []
            ymean_rows = []
            ysd_rows = []
            f_ratio_rows = []
            q_delta_rows = []
            t0 = time.time()
            for r in range(args.reps):
                seed = args.seed + 100000 * ci + r
                rng = np.random.default_rng(seed)
                G = rng.binomial(2, float(cell["maf"]), size=(args.N, 1)).astype(float)
                y = simulate(cell, G, rng)
                block_rng = np.random.RandomState(seed)
                block_ids = block_rng.randint(0, args.blocks, args.N).astype(np.int32)
                rif, q_used, f_used, f_ratio = make_rif(core, arm, y, cell, TAUS, q_pop, f_pop)
                result = stage1_from_rif(core, G, rif, block_ids, args.blocks, args.threads)
                beta_rows.append(result["beta"][None, :])
                se_rows.append(result["se"][None, :])
                cov_rows.append(result["cov"][np.triu_indices(len(TAUS))][None, :])
                d_rows.append([result["dresid"]])
                sigma_rows.append(result["sigma0"])
                q_rows.append(q_used)
                f_rows.append(f_used)
                ymean_rows.append(float(np.mean(y)))
                ysd_rows.append(float(np.std(y, ddof=1)))
                f_ratio_rows.append(f_ratio)
                q_delta_rows.append(q_used - q_pop)
                theta = A @ result["beta"]
                aggregate.append(dict(cell=cell["name"], arm=arm, rep=r, N=args.N,
                                      theta_mu=float(theta[0]), theta_sigma2=float(theta[1]),
                                      target_projection_mu=float(target_proj[0]),
                                      target_projection_sigma2=float(target_proj[1]),
                                      target_scientific_mu=float(target_sci[0]),
                                      target_scientific_sigma2=float(target_sci[1]),
                                      f_ratio_mean=float(np.mean(f_ratio)),
                                      q_delta_rms=float(np.sqrt(np.mean((q_used - q_pop) ** 2)))))
                if (r + 1) % max(1, min(10, args.reps)) == 0:
                    print(f"{name}: {r + 1}/{args.reps} ({time.time() - t0:.1f}s)", flush=True)

            pfx = os.path.join(stage_dir, name + "__current")
            write_array(pfx + ".beta.bin", np.asarray(beta_rows))
            write_array(pfx + ".se.bin", np.asarray(se_rows))
            write_array(pfx + ".cov.bin", np.asarray(cov_rows))
            write_array(pfx + ".dresid.bin", np.asarray(d_rows))
            write_array(pfx + ".sigma0.bin", np.asarray(sigma_rows))
            write_array(pfx + ".q_tau.bin", np.asarray(q_rows))
            write_array(pfx + ".fhat.bin", np.asarray(f_rows))
            write_array(pfx + ".y_mean.bin", np.asarray(ymean_rows))
            write_array(pfx + ".y_sd.bin", np.asarray(ysd_rows))
            write_array(pfx + ".maf0.bin", np.full(args.reps, float(cell["maf"])))
            json.dump(dict(cell=arm_cell, path="current", reps=args.reps, N=args.N,
                           blocks=args.blocks, M=M, T=len(TAUS),
                           ne=len(TAUS) * (len(TAUS) + 1) // 2,
                           source_commit=source_commit, diagnostic_arm=arm),
                      open(pfx + ".json", "w", encoding="utf-8"), indent=2)
            write_array(os.path.join(out, name + ".f_ratio.bin"), np.asarray(f_ratio_rows))
            write_array(os.path.join(out, name + ".q_delta.bin"), np.asarray(q_delta_rows))

    with open(os.path.join(out, "density_probe_rows.jsonl"), "w", encoding="utf-8") as fh:
        for row in aggregate:
            fh.write(json.dumps(row) + "\n")
    aux_fields = ["cell", "arm", "rep", "N", "theta_mu", "theta_sigma2",
                  "target_projection_mu", "target_projection_sigma2",
                  "target_scientific_mu", "target_scientific_sigma2",
                  "f_ratio_mean", "q_delta_rms"]
    import csv
    with open(os.path.join(out, "density_aux.csv"), "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=aux_fields)
        writer.writeheader()
        writer.writerows(aggregate)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
