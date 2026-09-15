#!/usr/bin/env python3
"""Tier 1 fixtures: one clean Stage 1 artefact set, plus one per refusal path.

Runs the SHIPPED Stage 1 (fungwas_stage1.core.run_stage1) on a small simulated
dataset and writes the exact file set the Stage 1 CLI writes. The truth is
written alongside as truth.npz so the R side can check recovery against a known
theta rather than against itself.

    python tests/tier1/t1_make_fixtures.py --out <dir>

Everything is deterministic given --seed.
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


# The refusal fixtures. Each name is the `corrupt` argument to
# write_stage1_files and the directory the artefacts land in; the R side looks
# them up by the same name.
CORRUPTIONS = [
    "cov_ids_shuffled",       # covariance rows do not line up with the TSV
    "no_cov_ids",             # the alignment sidecar is missing entirely
    "no_d_resid",             # no per-variant scale for the pooled covariance
    "no_cov_meta",            # nothing to check the Sigma_0 sidecar against
    "sigma0_wrong_taus",      # Sigma_0 built on a different tau grid
    "sigma0_wrong_covariates",# Sigma_0 built on a different covariate set
    "sigma0_wrong_n",         # Sigma_0 built on a different sample
    "sigma0_wrong_loco_scale",# Sigma_0 built under a different LOCO scaling
    "future_schema",          # a sidecar schema this Stage 2 cannot read
    "duplicate_snp_ids",      # multi-allelic / duplicated rsIDs
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--N", type=int, default=20000)
    ap.add_argument("--n-null", type=int, default=200)
    ap.add_argument("--n-causal", type=int, default=12)
    ap.add_argument("--n-taus", type=int, default=15)
    ap.add_argument("--blocks", type=int, default=200)
    ap.add_argument("--maf", type=float, default=0.30)
    ap.add_argument("--z", type=float, default=14.0)
    ap.add_argument("--n-covar", type=int, default=2)
    ap.add_argument("--nbg", type=int, default=8000)
    ap.add_argument("--seed", type=int, default=20260914)
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    taus = np.round(np.linspace(0.05, 0.95, a.n_taus), 3)
    t0 = time.time()

    rng = np.random.default_rng(a.seed)
    design = gs.build_design("locscale", a.maf, a.n_causal, a.z, a.N, taus, rng,
                             n_null=a.n_null)
    res = gs.run_stage1_on(design, a.N, taus, a.seed, n_blocks=a.blocks,
                           n_covar=a.n_covar, n_threads=a.threads)
    t_stage1 = time.time() - t0

    # exact population truth for the causal variants
    tr = gs.Truth(design, nbg=a.nbg, seed=a.seed + 991)
    q_pop, f_pop, b_c = tr.b_true(taus)
    b_true = np.zeros_like(res["betas"])
    b_true[design.causal] = b_c

    clean = os.path.join(a.out, "clean")
    os.makedirs(clean, exist_ok=True)
    paths = gs.write_stage1_files(os.path.join(clean, "sim"), taus, res,
                                  phenotype="simpheno", cov_dtype="float32")

    # the same scan written with the int8 covariance, to exercise that reader
    i8 = os.path.join(a.out, "int8")
    os.makedirs(i8, exist_ok=True)
    gs.write_stage1_files(os.path.join(i8, "sim"), taus, res,
                          phenotype="simpheno", cov_dtype="int8")

    for c in CORRUPTIONS:
        d = os.path.join(a.out, c)
        os.makedirs(d, exist_ok=True)
        gs.write_stage1_files(os.path.join(d, "sim"), taus, res,
                              phenotype="simpheno", cov_dtype="float32", corrupt=c)

    # A degenerate tau grid: two taus whose sample quantiles collide, which is
    # what a tied or censored phenotype produces and what makes Sigma singular.
    dg = os.path.join(a.out, "degenerate_taus")
    os.makedirs(dg, exist_ok=True)
    rng2 = np.random.default_rng(a.seed + 5)
    d2 = gs.build_design("null", a.maf, 0, 0.0, a.N, taus, rng2, n_null=40)
    G2 = rng2.binomial(2, np.broadcast_to(d2.maf, (a.N, len(d2.maf)))).astype(float)
    y2 = gs.simulate_phenotype(d2, G2, rng2)
    y2 = np.round(y2 * 2.0) / 2.0            # heavy rounding -> tied quantiles
    gs.add_stage1_to_path()
    from fungwas_stage1 import core, pooled_cov
    res2 = core.run_stage1(G2, y2, taus, covariates=None, n_blocks=a.blocks,
                           seed=a.seed, n_threads=a.threads)
    rif2, q2, f2 = core.compute_rif(y2, taus)
    X2 = np.ones((a.N, 1))
    rif2r, _ = core.residualize(rif2, X2)
    res2["sigma0"] = pooled_cov.compute_sigma0(rif2r)
    res2["n_samples"] = a.N
    res2["covariate_names"] = []
    gs.write_stage1_files(os.path.join(dg, "sim"), taus, res2,
                          phenotype="tiedpheno", cov_dtype="float32")

    # A covariance that is not positive semi-definite, which the int8
    # quantisation can produce and which the GLS solve must not silently accept.
    npsd = os.path.join(a.out, "non_psd")
    os.makedirs(npsd, exist_ok=True)
    res3 = dict(res)
    covs = [np.array(c, dtype=float, copy=True) for c in res["covariances"]]
    for j in range(len(covs)):
        w, V = np.linalg.eigh(covs[j])
        w[0] = -abs(w[-1]) * 0.5          # force one strongly negative eigenvalue
        covs[j] = V @ np.diag(w) @ V.T
    res3["covariances"] = covs
    gs.write_stage1_files(os.path.join(npsd, "sim"), taus, res3,
                          phenotype="simpheno", cov_dtype="float32")

    # Plain text as well as npz: the R side reads the text, and nothing in the
    # suite should need a Python round trip to check a number.
    np.savetxt(os.path.join(a.out, "b_true.tsv"), b_true, delimiter="\t")
    np.savetxt(os.path.join(a.out, "theta_latent.tsv"),
               np.column_stack([design.a, design.b, design.c, design.d]),
               delimiter="\t", header="a\tb\tc\td", comments="")
    np.savez(os.path.join(a.out, "truth.npz"),
             taus=taus, b_true=b_true, q_pop=q_pop, f_pop=f_pop,
             causal=design.causal, maf=design.maf,
             a=design.a, b=design.b, c=design.c, d=design.d,
             y_sd=res["y_sd"], y_mean=res["y_mean"],
             d_resid=res["d_resid"], betas=res["betas"], se=res["se"],
             q_tau_hat=res["q_tau"], fhat=res["fhat"],
             sigma0=res["sigma0"])

    meta = dict(N=a.N, n_null=a.n_null, n_causal=a.n_causal, n_taus=a.n_taus,
                blocks=a.blocks, maf=a.maf, z_target=a.z, n_covar=a.n_covar,
                seed=a.seed, taus=[float(t) for t in taus],
                phenotype="simpheno", covariates=res["covariate_names"],
                stage1_seconds=round(t_stage1, 3),
                corruptions=CORRUPTIONS + ["degenerate_taus", "non_psd", "int8"],
                paths={k: os.path.relpath(v, a.out) for k, v in paths.items()})
    with open(os.path.join(a.out, "fixtures.json"), "w") as fh:
        json.dump(meta, fh, indent=2)

    print(f"tier1 fixtures in {a.out}: Stage 1 took {t_stage1:.2f}s, "
          f"{res['betas'].shape[0]} variants x {len(taus)} taus, "
          f"{len(CORRUPTIONS) + 3} refusal fixtures")


if __name__ == "__main__":
    main()
