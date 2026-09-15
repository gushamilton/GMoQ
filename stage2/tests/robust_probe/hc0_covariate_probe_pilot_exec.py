#!/usr/bin/env python3
"""Independent genotype-weighted HC0 check for the rich covariate probe.

This is a diagnostic script, not a Stage 1 or Stage 2 implementation.  It
recreates the covariate-basis probe's two conditional-null cells and seeds,
uses the shipped Stage 1 kernels for the point estimates and delete-one-block
jackknife, and forms the HC0 sandwich here from individual residualised
scores.  The OLS weight matrix is fixed at the two-parameter Normal design
W = [1, qnorm(tau) / 2].

The output has one row per dataset/SNP/method/parameter and a dataset-level
summary.  The first SNP is the conditional-null SNP correlated with C; SNPs
2--6 are independent marginal-null controls.  A dataset is the Monte Carlo
unit, so the summary's MCSE is based on 100 (or more) datasets rather than on
the six SNP rows.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from typing import Dict, Iterable, Tuple

import numpy as np
from scipy.special import ndtri


def _load_pilot(validity_repo: str):
    """Load the exact pilot_stage1 module without mutating its source."""
    here = os.path.abspath(os.path.join(validity_repo, "tests", "validity_pilot"))
    if here not in sys.path:
        sys.path.insert(0, here)
    import pilot_stage1 as p  # noqa: E402
    return p


def _basis(C: np.ndarray) -> np.ndarray:
    c = np.asarray(C)[:, 0]
    return np.column_stack((c, c * c, c * c * c))


def _fit_from_rif(p, G, rif, X, block_ids, blocks, threads):
    """Use the shipped block kernel for the y-first comparator."""
    return p.stage1_from_rif(G, rif, X, block_ids, blocks, threads)


def _arm(p, G, y, C, path: str, seed: int, blocks: int, threads: int):
    """Return Stage 1 output and the individual arrays needed by HC0."""
    from fungwas_stage1 import core, pooled_cov

    Z = _basis(C)
    X = np.column_stack((np.ones(len(y)), Z))
    t0 = time.perf_counter()
    if path == "current":
        res = core.run_stage1(
            G, y, p.TAUS, covariates=Z, n_blocks=blocks,
            seed=seed, n_threads=threads,
        )
        rif, q_tau, fhat = core.compute_rif(y, p.TAUS)
    elif path == "y_first":
        y_resid, _ = core.residualize(y, X)
        rif, q_tau, fhat = core.compute_rif(y_resid, p.TAUS)
        np.random.seed(seed)
        block_ids = np.random.randint(0, blocks, len(y)).astype(np.int32)
        res = _fit_from_rif(p, G, rif, X, block_ids, blocks, threads)
    else:
        raise ValueError(path)
    stage_seconds = time.perf_counter() - t0
    rif_resid, Q = core.residualize(rif, X)
    g_resid, _ = core.residualize(G, X)
    sigma0 = pooled_cov.compute_sigma0(rif_resid)
    res["sigma0_hc0_probe"] = sigma0
    res["q_tau_hc0_probe"] = q_tau
    res["fhat_hc0_probe"] = fhat
    res["rif_resid_hc0_probe"] = rif_resid
    res["g_resid_hc0_probe"] = g_resid
    res["stage_seconds_hc0_probe"] = stage_seconds
    return res


def _theta_and_se(beta: np.ndarray, cov_tau: np.ndarray, A: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    theta = A @ beta
    cov_theta = A @ cov_tau @ A.T
    diag = np.diag(cov_theta)
    if np.any(diag < -1e-9):
        raise FloatingPointError(f"negative parameter covariance diagonal: {diag}")
    return theta, np.sqrt(np.maximum(diag, 0.0))


def _covariances(res: dict) -> Tuple[Dict[str, np.ndarray], float]:
    """Form pooled and independently computed genotype-weighted HC0 covariances."""
    rr = res["rif_resid_hc0_probe"]
    gr = res["g_resid_hc0_probe"]
    betas = np.asarray(res["betas"])
    sigma0 = np.asarray(res["sigma0_hc0_probe"])
    dresid = np.asarray(res["d_resid"], dtype=float)
    pooled = []
    hc0 = []
    t0 = time.perf_counter()
    # The score is g_r * (RIF_r - g_r beta), exactly the fixed-genotype OLS
    # residual score requested for this probe.  Keep this in a variant loop so
    # peak memory is O(N*T), not O(N*M*T).
    for j in range(gr.shape[1]):
        gj = gr[:, j]
        d = dresid[j]
        pooled.append(sigma0 / d)
        e = rr - gj[:, None] * betas[j][None, :]
        score = gj[:, None] * e
        hc0.append((score.T @ score) / (d * d))
    seconds = time.perf_counter() - t0
    return {"pooled": pooled, "jackknife": list(res["covariances"]), "hc0": hc0}, seconds


def _write_rows(path: str, rows: Iterable[dict]):
    import csv
    rows = list(rows)
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def _summary(rows: list[dict], out_path: str, reps: int):
    """Summarise with dataset as the MC unit, retaining SNP1 separately."""
    import csv
    from collections import defaultdict

    groups = defaultdict(list)
    for row in rows:
        key = tuple(row[k] for k in ("cell", "basis", "path", "snp_group", "method", "parameter"))
        groups[key].append(row)
    out = []
    for key, vals in sorted(groups.items()):
        # Rows contain one value per SNP within a dataset.  Aggregate SNPs
        # 2--6 within each dataset before taking MC means/MCSEs.
        by_rep = defaultdict(list)
        for v in vals:
            by_rep[int(v["rep"])].append(v)
        metric = {"estimate": [], "se": [], "z": [], "reject": [], "cover": []}
        for rep in sorted(by_rep):
            vv = by_rep[rep]
            for k in metric:
                metric[k].append(float(np.mean([float(x[k]) for x in vv])))
        mean_est = float(np.mean(metric["estimate"]))
        emp_sd = float(np.std(metric["estimate"], ddof=1))
        mean_se = float(np.mean(metric["se"]))
        mean_z = float(np.mean(metric["z"]))
        out.append(dict(
            cell=key[0], basis=key[1], path=key[2], snp_group=key[3],
            method=key[4], parameter=key[5], datasets=len(metric["estimate"]),
            mean_estimate=mean_est, empirical_sd_dataset_mean=emp_sd,
            mean_se=mean_se, mean_se_over_empirical_sd=(mean_se / emp_sd if emp_sd else np.nan),
            mean_z=mean_z, mcse_z=float(np.std(metric["z"], ddof=1) / np.sqrt(len(metric["z"]))),
            rejection=float(np.mean(metric["reject"])),
            rejection_mcse=float(np.std(metric["reject"], ddof=1) / np.sqrt(len(metric["reject"]))),
            coverage=float(np.mean(metric["cover"])),
            coverage_mcse=float(np.std(metric["cover"], ddof=1) / np.sqrt(len(metric["cover"]))),
        ))
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--validity-repo", required=True)
    ap.add_argument("--repo", required=True, help="Stage 1 repository containing the current Python package")
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--N", type=int, default=10000)
    ap.add_argument("--blocks", type=int, default=1000)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--seed", type=int, default=9142026)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    p = _load_pilot(args.validity_repo)
    p.gs.add_stage1_to_path(args.repo)
    from fungwas_stage1 import core
    if not core.HAVE_CPP:
        raise RuntimeError("HC0 probe requires the shipped Stage 1 C++ extension")

    # Fixed two-parameter Normal vQTL OLS design, matching qnorm(tau)/(2 sd)
    # with sd=1 and the probe's oracle W arm.
    W = np.column_stack((np.ones(len(p.TAUS)), ndtri(p.TAUS) / 2.0))
    A = np.linalg.solve(W.T @ W, W.T)
    parameter_names = ("beta_mu", "beta_sigma2")
    cells = [p.cells()[4], p.cells()[9]]
    basis_name = "rich_C_C2_C3"
    source_script = os.path.join(args.validity_repo, "tests", "validity_pilot", "covariate_basis_probe.py")
    source_hash = hashlib.sha256(open(source_script, "rb").read()).hexdigest()
    meta = dict(
        design="independent_HC0_fixed_OLS_two_parameter_rich_covariate_null",
        source_covariate_probe=source_script, source_covariate_probe_sha256=source_hash,
        validity_repo=args.validity_repo, stage1_repo=args.repo,
        N=args.N, reps=args.reps, blocks=args.blocks, threads=args.threads,
        seed=args.seed, taus=[float(x) for x in p.TAUS], W=W.tolist(),
        basis=basis_name, cells=[c["name"] for c in cells],
        mc_unit="dataset; SNP1 is C-correlated and SNP2-SNP6 are independent marginal nulls",
        covariance_form="HC0: score_i = g_resid_i*(RIF_resid_i - g_resid_i*beta); score' score / d_resid^2",
        target="conditional-null scientific target zero for every SNP after rich C,C^2,C^3 residualisation",
        host=platform.node(), generated_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    )
    with open(os.path.join(args.out, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)

    rows = []
    timing = []
    for ci, cell in enumerate(cells):
        for rep in range(args.reps):
            seed = args.seed + 100000 * (4 + ci * 5) + rep
            rng = np.random.default_rng(seed)
            G, y, C = p.draw_covariate_data(cell, args.N, rng)
            for path in ("current", "y_first"):
                res = _arm(p, G, y, C, path, seed, args.blocks, args.threads)
                covs, hc0_seconds = _covariances(res)
                timing.append(dict(cell=cell["name"], path=path, rep=rep,
                                   stage_seconds=float(res["stage_seconds_hc0_probe"]),
                                   hc0_seconds=float(hc0_seconds)))
                for j in range(G.shape[1]):
                    group = "snp1_correlated" if j == 0 else "snp2_6_independent"
                    for method, matrices in covs.items():
                        theta, se = _theta_and_se(np.asarray(res["betas"])[j], matrices[j], A)
                        for k, pname in enumerate(parameter_names):
                            est = float(theta[k]); se0 = float(se[k])
                            z = est / se0 if se0 > 0 else np.nan
                            reject = float(abs(z) > 1.959963984540054) if np.isfinite(z) else np.nan
                            cover = float(abs(est) <= 1.959963984540054 * se0) if np.isfinite(z) else np.nan
                            rows.append(dict(cell=cell["name"], basis=basis_name, path=path,
                                             rep=rep, snp=j + 1, snp_group=group,
                                             method=method, parameter=pname, estimate=est,
                                             se=se0, z=z, reject=reject, cover=cover))
    _write_rows(os.path.join(args.out, "rows.csv"), rows)
    _summary(rows, os.path.join(args.out, "summary.csv"), args.reps)
    _write_rows(os.path.join(args.out, "timing.csv"), timing)
    print(json.dumps(dict(meta=meta, rows=len(rows), timing_rows=len(timing)), indent=2))


if __name__ == "__main__":
    main()
