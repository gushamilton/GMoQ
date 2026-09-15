#!/usr/bin/env python3
"""Diagnostic fourth-order KDE repair prototype for the Normal variance arm.

This is deliberately a test-only harness.  It reuses the density diagnosis
helpers and leaves production ``compute_rif`` unchanged.  Each frozen dataset
is analysed twice with the shipped Stage 1 residualisation, C++ block-score,
and covariance kernels:

  A: production empirical type-8 quantile plus Gaussian Silverman KDE;
  D: the same empirical quantile and bandwidth with
     K(u) = (1.5 - 0.5*u^2) * phi(u).

The fourth-order estimate is floored at 1e-8 only for the RIF denominator.
Raw non-positive and floored tau densities are recorded explicitly.  The
existing pilot_stage2.R wrapper can consume the output layout and fit the
installed pooled+OLS and jackknife+OLS paths.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
import time

import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from density_probe import (  # noqa: E402
    TAUS,
    exact_projection,
    git_value,
    simulate,
    stage1_from_rif,
    write_array,
)


M = 1
ARMS = ("A_production_kde", "D_fourth_order_kde")
DENSITY_FLOOR = 1e-8


def cells() -> list[dict]:
    return [
        {"name": "repair_v0", "family": "normal_variance", "maf": 0.30, "a": 0.0, "v": 0.0},
        {"name": "repair_v06", "family": "normal_variance", "maf": 0.30, "a": 0.0, "v": 0.06},
        {"name": "repair_v24", "family": "normal_variance", "maf": 0.30, "a": 0.0, "v": 0.24},
    ]


def production_density_raw(y: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Raw production Gaussian KDE, before its shipped density floor."""
    n = len(y)
    h = 1.06 * np.std(y) * n ** (-0.2)
    u = (y[:, None] - q[None, :]) / h
    return np.mean(np.exp(-0.5 * u * u), axis=0) / (h * np.sqrt(2.0 * np.pi))


def fourth_order_rif(y: np.ndarray, taus: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return fourth-order RIF plus q, raw f, used f, and raw u moments."""
    n = len(y)
    q = np.quantile(y, taus, method="median_unbiased")
    h = 1.06 * np.std(y) * n ** (-0.2)
    u = (y[:, None] - q[None, :]) / h
    phi = np.exp(-0.5 * u * u) / np.sqrt(2.0 * np.pi)
    raw = np.mean((1.5 - 0.5 * u * u) * phi, axis=0) / h
    used = np.maximum(raw, DENSITY_FLOOR)
    indicator = (y[:, None] <= q[None, :]).astype(float)
    rif = q[None, :] + (taus[None, :] - indicator) / used[None, :]
    return rif, q, raw, used, np.full(len(taus), h)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=os.path.abspath(os.path.join(HERE, "..", "..")))
    ap.add_argument("--out", required=True)
    ap.add_argument("--run-id", default="")
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--N", type=int, default=10000)
    ap.add_argument("--blocks", type=int, default=1000)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--seed", type=int, default=20260915)
    ap.add_argument("--source-commit", default="")
    ap.add_argument("--source-status", default="")
    args = ap.parse_args()

    repo = os.path.abspath(args.repo)
    sys.path.insert(0, os.path.join(repo, "stage1"))
    from fungwas_stage1 import core

    core.require_cpp_extension("density repair probe")
    out = os.path.abspath(args.out)
    run_id = args.run_id or os.path.basename(out.rstrip(os.sep))
    if os.path.exists(out):
        raise SystemExit(f"refusing existing output directory: {out}")
    stage_dir = os.path.join(out, "stage1")
    truth_dir = os.path.join(out, "truth")
    os.makedirs(stage_dir, exist_ok=True)
    os.makedirs(truth_dir, exist_ok=True)

    source_commit = args.source_commit or git_value(repo, "rev-parse", "HEAD")
    source_status = args.source_status or git_value(repo, "status", "--short", "--branch")
    script_sha256 = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()
    meta = dict(
        design="paired_density_repair_diagnosis",
        run_id=run_id,
        source_commit=source_commit,
        source_status=source_status,
        script_sha256=script_sha256,
        repo_head=git_value(repo, "rev-parse", "HEAD"),
        host=platform.node(),
        python=sys.version,
        numpy=np.__version__,
        N=args.N,
        reps=args.reps,
        blocks=args.blocks,
        threads=args.threads,
        taus=TAUS.tolist(),
        M=M,
        arms=list(ARMS),
        seed=args.seed,
        density_floor=DENSITY_FLOOR,
        kernel_repair="(1.5 - 0.5*u^2) * phi(u), same h=1.06*sd(Y)*N^(-1/5)",
        mc_unit="dataset (one shared phenotype/genotype draw; arms are paired)",
        target_definition=("projection is the exact population RIF slope projected onto "
                           "the Normal mean/variance W; scientific target is (a,v)"),
        claim_status="diagnostic only; no production repair is implied",
    )
    with open(os.path.join(out, "density_repair_meta.json"), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)

    runtime_rows = []
    aux_rows = []
    for ci, cell in enumerate(cells()):
        q_pop, f_pop, b_tau = exact_projection(cell, TAUS)
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
                density_floor=DENSITY_FLOOR,
                kernel=("production Gaussian Silverman KDE" if arm == "A_production_kde"
                        else "fourth-order (1.5 - 0.5*u^2) * phi(u) KDE"),
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
            raw_nonpositive_rows = []
            floor_rows = []
            raw_min_rows = []
            h_rows = []
            t0 = time.time()
            for r in range(args.reps):
                seed = args.seed + 100000 * ci + r
                rng = np.random.default_rng(seed)
                G = rng.binomial(2, float(cell["maf"]), size=(args.N, 1)).astype(float)
                y = simulate(cell, G, rng)
                block_rng = np.random.RandomState(seed)
                block_ids = block_rng.randint(0, args.blocks, args.N).astype(np.int32)

                if arm == "A_production_kde":
                    rif, q_used, f_used = core.compute_rif(y, TAUS)
                    raw = production_density_raw(y, q_used)
                    h = np.full(len(TAUS), 1.06 * np.std(y) * args.N ** (-0.2))
                else:
                    rif, q_used, raw, f_used, h = fourth_order_rif(y, TAUS)
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
                f_ratio_rows.append(raw / np.maximum(exact_density_at_q(cell, q_used), 1e-300))
                raw_nonpositive_rows.append(int(np.count_nonzero(raw <= 0.0)))
                floor_rows.append(int(np.count_nonzero(raw <= DENSITY_FLOOR)))
                raw_min_rows.append(float(np.min(raw)))
                h_rows.append(float(h[0]))
                aux_rows.append(dict(cell=cell["name"], arm=arm, rep=r, N=args.N,
                                     raw_nonpositive_tau=int(np.count_nonzero(raw <= 0.0)),
                                     floor_tau=int(np.count_nonzero(raw <= DENSITY_FLOOR)),
                                     raw_density_min=float(np.min(raw)),
                                     raw_density_mean=float(np.mean(raw)),
                                     raw_over_exact_mean=float(np.mean(f_ratio_rows[-1])),
                                     q_delta_rms=float(np.sqrt(np.mean((q_used - q_pop) ** 2))),
                                     bandwidth=float(h[0])))
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
                           source_commit=source_commit, diagnostic_arm=arm,
                           density_floor=DENSITY_FLOOR,
                           kernel=("production Gaussian Silverman KDE" if arm == "A_production_kde"
                                   else "fourth-order Gaussian kernel")),
                      open(pfx + ".json", "w", encoding="utf-8"), indent=2)
            write_array(os.path.join(out, name + ".raw_over_exact.bin"), np.asarray(f_ratio_rows))
            write_array(os.path.join(out, name + ".raw_nonpositive.bin"), np.asarray(raw_nonpositive_rows))
            write_array(os.path.join(out, name + ".floor_failures.bin"), np.asarray(floor_rows))
            write_array(os.path.join(out, name + ".raw_min.bin"), np.asarray(raw_min_rows))
            write_array(os.path.join(out, name + ".bandwidth.bin"), np.asarray(h_rows))
            runtime_rows.append(dict(cell=cell["name"], arm=arm, N=args.N, reps=args.reps,
                                     elapsed_seconds=time.time() - t0,
                                     source_commit=source_commit))

    with open(os.path.join(out, "density_repair_aux.csv"), "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(aux_rows[0]))
        writer.writeheader()
        writer.writerows(aux_rows)
    with open(os.path.join(out, "density_repair_runtime.csv"), "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(runtime_rows[0]))
        writer.writeheader()
        writer.writerows(runtime_rows)


def exact_density_at_q(cell: dict, q: np.ndarray) -> np.ndarray:
    """Exact marginal density at q for the Normal-mixture reference."""
    p = float(cell["maf"])
    a = float(cell["a"])
    g = np.arange(3.0)
    x = g - 2.0 * p
    pg = np.array([(1.0 - p) ** 2, 2.0 * p * (1.0 - p), p ** 2])
    sd = np.sqrt(1.0 + float(cell["v"]) * x)
    z = (np.asarray(q)[:, None] - a * x[None, :]) / sd[None, :]
    return np.sum(pg[None, :] * np.exp(-0.5 * z * z) /
                  (np.sqrt(2.0 * np.pi) * sd[None, :]), axis=1)


if __name__ == "__main__":
    main()
