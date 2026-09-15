#!/usr/bin/env python3
"""Small paired known-truth pilot for the shipped G-MoQ engines.

The pilot is deliberately separate from tier 2.  It is cheap enough to run
before a confirmatory grid and adds the covariate scenario that tier 2 does
not contain: a conditionally null SNP correlated with a covariate that shifts
both the mean and the scale of the phenotype.

All point estimates and covariance matrices come from the shipped Stage 1
functions.  The ``y_first`` arm is an explicitly labelled comparator built
from the same shipped RIF, residualisation and block-score kernels; it is not
a replacement for the production path.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import time

import numpy as np
from scipy.special import expit, logit

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.abspath(os.path.join(HERE, "..", "lib"))
sys.path.insert(0, LIB)
import gmoq_sim as gs  # noqa: E402


# Production grid.  A deliberately smaller grid is not used for the main
# pilot because shape leakage and covariance-source differences depend on the
# tau range and rank.
TAUS = np.round(np.linspace(0.02, 0.98, 45), 3)
M = 6


def git_value(repo: str, *args: str) -> str:
    try:
        return subprocess.check_output(["git", "-C", repo, *args], text=True).strip()
    except Exception as exc:  # pragma: no cover - provenance must not stop a pilot
        return f"unavailable:{exc}"


def cells() -> list[dict]:
    out = []
    for maf in (0.30, 0.01):
        suffix = "common" if maf > 0.1 else "rare"
        out.extend([
            dict(name=f"null_{suffix}", family="null", maf=maf, a=0.0, b=0.0),
            dict(name=f"mean_{suffix}", family="location", maf=maf, a=0.12, b=0.0),
            dict(name=f"scale_{suffix}", family="scale", maf=maf, a=0.0, b=0.12),
            dict(name=f"combined_{suffix}", family="locscale", maf=maf, a=0.12, b=0.12),
            dict(name=f"cov_both_{suffix}", family="cov_both", maf=maf, a=0.0, b=0.0),
        ])
    # The .12 direct-effect cells are deliberately a strong-effect stress
    # contrast.  Retain one routine-strength location+scale cell so a failure
    # there is distinguishable from nonlinearity at the stress setting.
    out.append(dict(name="combined_small_common", family="locscale", maf=0.30,
                    a=0.03, b=0.03))
    return out


def make_design(cell: dict) -> gs.Design:
    maf = float(cell["maf"])
    causal = np.array([0], dtype=int) if cell["family"] != "cov_both" else np.array([], dtype=int)
    d = gs.Design(
        family="null" if cell["family"] == "cov_both" else cell["family"],
        maf=np.full(M, maf), causal=causal, mu0=0.0, s0=1.0,
        c0=0.0, d0=0.0, mix=None,
    )
    if cell["family"] != "cov_both":
        d.a[0] = float(cell["a"])
        d.b[0] = float(cell["b"])
    return d


def draw_covariate_data(cell: dict, N: int, rng: np.random.Generator):
    """Draw C first, then G|C, so the null SNP has no direct Y effect."""
    maf = float(cell["maf"])
    C = rng.normal(size=N)
    # A logistic genotype-probability model gives a clean C -> G association:
    # conditional on C, Y is independent of the SNP.  The marginal MAF is
    # close to the requested value, and the realised value is recorded in the
    # output rather than silently treated as exact.
    alpha = 1.0 if maf > 0.1 else 1.5
    p = expit(logit(maf) + alpha * C)
    G = rng.binomial(2, p[:, None], size=(N, 1)).astype(float)
    if M > 1:
        G = np.column_stack([G, rng.binomial(2, maf, size=(N, M - 1))])
    # Mean and log-SD effects are both real nuisance effects.  The SNP has no
    # direct effect conditional on C; its scientific target is exactly zero.
    Y = 0.50 * C + np.exp(0.45 * C) * rng.normal(size=G.shape[0])
    return G, Y, C[:, None]


def stage1_from_rif(G, rif, X, block_ids, n_blocks, threads):
    """Run the shipped block-score and covariance kernels on a supplied RIF."""
    from fungwas_stage1 import core

    rif_resid, Q = core.residualize(rif, X)
    stats = core.compute_block_scores(
        G, rif_resid, Q, block_ids, n_blocks, n_threads=threads,
    )
    betas = np.full((G.shape[1], rif.shape[1]), np.nan)
    ses = np.full_like(betas, np.nan)
    covs = []
    d_resid = np.full(G.shape[1], np.nan)
    for j in range(G.shape[1]):
        if np.isnan(stats[j, 0, 0]):
            covs.append(np.full((rif.shape[1], rif.shape[1]), np.nan))
            continue
        beta, se, cov = core.compute_covariance_from_blocks(stats[j], n_blocks)
        betas[j], ses[j] = beta, se
        covs.append(cov)
        d_resid[j] = float(np.sum(stats[j][:, 0]))
    return dict(
        betas=betas, se=ses, covariances=covs, d_resid=d_resid,
        sigma0=None, q_tau=None, fhat=None,
    )


def run_current_and_yfirst(G, y, C, taus, seed, n_blocks, threads, include_yfirst=True):
    """Paired production and residualise-outcome-first Stage 1 arms."""
    from fungwas_stage1 import core, pooled_cov

    # Production arm: this is exactly core.run_stage1, including its own
    # block assignment and RIF-before-covariate-residualisation order.
    current = core.run_stage1(
        G, y, taus, covariates=C, n_blocks=n_blocks, seed=seed,
        n_threads=threads,
    )
    X = np.ones((len(y), 1)) if C is None else np.column_stack([np.ones(len(y)), C])

    # Pooled Sigma_0 is formed from the same residualised RIF as the current
    # regression, reproducing the shipped --sigma0 path.
    rif, q, f = core.compute_rif(y, taus)
    rif0, _ = core.residualize(rif, X)
    current["q_tau"], current["fhat"] = q, f
    current["sigma0"] = pooled_cov.compute_sigma0(rif0)
    current["_y_mean"] = float(np.mean(y))
    current["_y_sd"] = float(np.std(y, ddof=1))
    if not include_yfirst:
        return current, None

    # Comparator uses the same shipped kernels and the same block IDs as the
    # production call, but residualises y before calling core.compute_rif.
    y_resid, _ = core.residualize(y, X)
    rif_y, q_y, f_y = core.compute_rif(y_resid, taus)
    np.random.seed(seed)
    block_ids = np.random.randint(0, n_blocks, len(y)).astype(np.int32)
    yfirst = stage1_from_rif(G, rif_y, X, block_ids, n_blocks, threads)
    yfirst["q_tau"] = q_y
    yfirst["fhat"] = f_y
    yfirst["sigma0"] = pooled_cov.compute_sigma0(core.residualize(rif_y, X)[0])
    yfirst["_y_mean"] = float(np.mean(y_resid))
    yfirst["_y_sd"] = float(np.std(y_resid, ddof=1))
    return current, yfirst


def projection_reference(cell, path, Nref, seed, threads):
    """Independent reference projection for the covariate cells.

    This is not a fitted Stage 1 replicate.  It evaluates the population
    RIF-regression target using a fresh, much larger draw and the shipped
    ``compute_rif`` and residualisation functions.  It is retained as an
    approximation with its reference N in metadata, never labelled exact.
    """
    from fungwas_stage1 import core

    rng = np.random.default_rng(seed)
    maf = float(cell["maf"])
    G, y, C = draw_covariate_data(cell, Nref, rng)
    X = np.column_stack([np.ones(Nref), C])
    if path == "y_first":
        y, _ = core.residualize(y, X)
    rif, _, _ = core.compute_rif(y, TAUS)
    rif_r, Q = core.residualize(rif, X)
    G_r = G - Q @ (Q.T @ G)
    den = np.sum(G_r * G_r, axis=0)
    return np.sum(G_r[:, :, None] * rif_r[:, None, :], axis=0) / den[:, None]


def write_array(path, arr, dtype="<f8"):
    np.asarray(arr, dtype=dtype).tofile(path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--repo", default=os.path.abspath(os.path.join(HERE, "..", "..")))
    ap.add_argument("--source-commit", default="")
    ap.add_argument("--source-status", default="")
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--N", type=int, default=5000)
    ap.add_argument("--blocks", type=int, default=1000)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--nbg", type=int, default=100000)
    ap.add_argument("--Nref", type=int, default=200000)
    ap.add_argument("--seed", type=int, default=9142026)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    stage_dir = os.path.join(args.out, "stage1")
    os.makedirs(stage_dir, exist_ok=True)

    gs.add_stage1_to_path(args.repo)
    from fungwas_stage1 import core
    if not core.HAVE_CPP:
        raise RuntimeError("pilot requires the shipped Stage 1 C++ extension")

    meta = dict(
        design="paired_known_truth_location_scale_covariate",
        source_commit=(args.source_commit or git_value(args.repo, "rev-parse", "HEAD")),
        source_status=(args.source_status or git_value(args.repo, "status", "--short")),
        python=sys.version,
        numpy=np.__version__, scipy=__import__("scipy").__version__,
        host=platform.node(), platform=platform.platform(),
        generated_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        N=args.N, reps=args.reps, blocks=args.blocks, threads=args.threads,
        taus=[float(x) for x in TAUS], M=M, Nref=args.Nref, nbg=args.nbg,
        mc_unit="dataset (one shared phenotype/genotype draw; SNP rows are not independent MC replicates)",
        targets=dict(
            scientific="generative conditional mean coefficient a and log-SD coefficient b in s(g)=s0*exp(b*g/s0), with s0=1; the shipped beta_sigma2 basis is the local variance coefficient 2*b and the one-dosage variance ratio is exp(2*b); cov_both SNP target is exactly zero conditional on C",
            projection="population slope of the shipped Stage-1 RIF after the named preprocessing, then projection onto the fitted W; exact for regular cells up to Truth background MC and reference-MC for cov_both",
        ),
    )
    json.dump(meta, open(os.path.join(args.out, "pilot_meta.json"), "w"), indent=2)

    truth_dir = os.path.join(args.out, "truth")
    os.makedirs(truth_dir, exist_ok=True)
    for cell_i, cell in enumerate(cells()):
        name = cell["name"]
        d = make_design(cell)
        if cell["family"] == "cov_both":
            proj = {
                "current": projection_reference(cell, "current", args.Nref, args.seed + 10000 + cell_i, args.threads),
                "y_first": projection_reference(cell, "y_first", args.Nref, args.seed + 20000 + cell_i, args.threads),
            }
            truth_kind = "reference_mc"
        else:
            tr = gs.Truth(d, nbg=args.nbg, seed=args.seed + 30000 + cell_i)
            _, _, bc = tr.b_true(TAUS)
            proj0 = np.zeros((M, len(TAUS)))
            if len(d.causal):
                proj0[d.causal] = bc
            proj = {"current": proj0, "y_first": proj0}
            truth_kind = "conditional_cdf_background_mc"
        sci_a = np.zeros(M); sci_logsd_b = np.zeros(M); sci_var_b = np.zeros(M)
        if cell["family"] != "cov_both":
            sci_a[0] = cell["a"]
            sci_logsd_b[0] = cell["b"]
            # The shipped Normal vQTL basis is parameterised by a variance
            # effect: W_sigma2 = z/(2*sd).  Convert the DGM's log-SD b into
            # the scientific local variance coefficient 2*b at s0=1.
            sci_var_b[0] = 2.0 * cell["b"]
        json.dump(dict(cell=cell, maf=[cell["maf"]] * M,
                       scientific_a=sci_a.tolist(),
                       scientific_logsd_b=sci_logsd_b.tolist(),
                       scientific_variance_coefficient=sci_var_b.tolist(),
                       scientific_variance_ratio=[float(np.exp(x)) for x in sci_var_b],
                       projection={k: np.asarray(v).tolist() for k, v in proj.items()},
                       projection_truth_kind=truth_kind,
                       projection_reference_N=(args.Nref if truth_kind == "reference_mc" else args.nbg),
                       taus=TAUS.tolist()),
                  open(os.path.join(truth_dir, name + ".json"), "w"), indent=2)

    for cell_i, cell in enumerate(cells()):
        name = cell["name"]
        d = make_design(cell)
        paths = ["current", "y_first"] if cell["family"] == "cov_both" else ["current"]
        arrays = {}
        for path in paths:
            arrays[path] = {k: [] for k in ("beta", "se", "cov", "dresid", "sigma0", "q_tau", "fhat", "y_mean", "y_sd", "maf0")}
        t0 = time.time()
        for r in range(args.reps):
            seed = args.seed + 100000 * cell_i + r
            rng = np.random.default_rng(seed)
            G = rng.binomial(2, np.broadcast_to(d.maf, (args.N, M))).astype(float)
            if cell["family"] == "cov_both":
                G, y, C = draw_covariate_data(cell, args.N, rng)
                current, yfirst = run_current_and_yfirst(G, y, C, TAUS, seed, args.blocks, args.threads)
                chosen = {"current": current, "y_first": yfirst}
            else:
                y = gs.simulate_phenotype(d, G, rng)
                current, _ = run_current_and_yfirst(G, y, None, TAUS, seed, args.blocks, args.threads,
                                                     include_yfirst=False)
                chosen = {"current": current}
            for path in paths:
                res = chosen[path]
                arrays[path]["beta"].append(res["betas"])
                arrays[path]["se"].append(res["se"])
                triu = np.triu_indices(len(TAUS))
                arrays[path]["cov"].append(np.stack([np.asarray(c)[triu] for c in res["covariances"]]))
                arrays[path]["dresid"].append(res["d_resid"])
                arrays[path]["sigma0"].append(res["sigma0"])
                arrays[path]["q_tau"].append(res["q_tau"])
                arrays[path]["fhat"].append(res["fhat"])
                arrays[path]["y_mean"].append(float(res.get("_y_mean", np.mean(y))))
                arrays[path]["y_sd"].append(float(res.get("_y_sd", np.std(y, ddof=1))))
                arrays[path]["maf0"].append(float(np.mean(G[:, 0]) / 2.0))
            if (r + 1) % max(1, min(10, args.reps)) == 0:
                print(f"{name}: {r + 1}/{args.reps} ({time.time() - t0:.1f}s)", flush=True)
        for path in paths:
            pfx = os.path.join(stage_dir, f"{name}__{path}")
            arr = arrays[path]
            for key in ("beta", "se", "cov", "dresid", "sigma0"):
                write_array(pfx + f".{key}.bin", np.asarray(arr[key]))
            for key in ("q_tau", "fhat", "y_mean", "y_sd", "maf0"):
                write_array(pfx + f".{key}.bin", np.asarray(arr[key]))
            json.dump(dict(cell=cell, path=path, reps=args.reps, N=args.N,
                           blocks=args.blocks, M=M, T=len(TAUS),
                           ne=len(TAUS) * (len(TAUS) + 1) // 2,
                           source_commit=meta["source_commit"]),
                      open(pfx + ".json", "w"), indent=2)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
