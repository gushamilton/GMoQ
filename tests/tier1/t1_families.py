#!/usr/bin/env python3
"""Tier 1 fixtures for the eight model families.

One Stage 1 run per family, written as a full Stage 1 artefact set plus the
exact population estimand b(tau) for every variant. The R side then fits the
whole candidate model set to each and checks both halves of the tournament
claim: the true model is selected, and the models that cannot represent the
truth are rejected.

    python tests/tier1/t1_families.py --out <dir>
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--N", type=int, default=50000)
    ap.add_argument("--n-null", type=int, default=100)
    ap.add_argument("--n-causal", type=int, default=30)
    ap.add_argument("--n-taus", type=int, default=15)
    ap.add_argument("--blocks", type=int, default=1000)
    ap.add_argument("--maf", type=float, default=0.30)
    ap.add_argument("--z", type=float, default=9.0)
    ap.add_argument("--z-shape", type=float, default=6.0)
    ap.add_argument("--nbg", type=int, default=8000)
    ap.add_argument("--seed", type=int, default=606)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--families", default=",".join(gs.FAMILIES))
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    taus = np.round(np.linspace(0.05, 0.95, a.n_taus), 3)
    fams = [f.strip() for f in a.families.split(",") if f.strip()]

    index = {}
    for i, fam in enumerate(fams):
        t0 = time.time()
        rng = np.random.default_rng(a.seed + 37 * i)
        # The shape families get a smaller target: their weight columns are
        # weakly identified, so a |Z| = 12 design asks for a skewness or tail
        # weight at which the linearisation the whole method rests on visibly
        # breaks, and the TRUE model then fails its own lack-of-fit test.
        z_fam = a.z_shape if fam in ("skew", "kurtosis") else a.z
        design = gs.build_design(fam, a.maf, a.n_causal, z_fam, a.N, taus, rng,
                                 n_null=a.n_null)
        res = gs.run_stage1_on(design, a.N, taus, a.seed + 37 * i,
                               n_blocks=a.blocks, n_covar=0, n_threads=a.threads)
        tr = gs.Truth(design, nbg=a.nbg, seed=a.seed + 37 * i + 991)
        q_pop, f_pop, b_c = tr.b_true(taus)
        b_true = np.zeros_like(res["betas"])
        if len(design.causal):
            b_true[design.causal] = b_c

        d = os.path.join(a.out, f"fam_{fam}")
        os.makedirs(d, exist_ok=True)
        gs.write_stage1_files(os.path.join(d, "sim"), taus, res,
                              phenotype=f"sim_{fam}", cov_dtype="float32")
        np.savetxt(os.path.join(d, "b_true.tsv"), b_true, delimiter="\t")
        np.savetxt(os.path.join(d, "q_pop.tsv"), q_pop, delimiter="\t")

        index[fam] = dict(
            dir=f"fam_{fam}", y_sd=res["y_sd"], y_mean=res["y_mean"],
            q_pop=[float(x) for x in q_pop],
            q_tau_hat=[float(x) for x in res["q_tau"]],
            causal=[int(x) + 1 for x in design.causal],   # 1-based for R
            z_target=float(z_fam), shrink=float(design.shrink),
            n_variants=int(res["betas"].shape[0]),
            mix=design.mix,
            seconds=round(time.time() - t0, 2))
        print(f"  {fam}: {time.time() - t0:.1f}s", flush=True)

    with open(os.path.join(a.out, "families.json"), "w") as fh:
        json.dump(dict(taus=[float(t) for t in taus], N=a.N, blocks=a.blocks,
                       maf=a.maf, z_target=a.z, n_null=a.n_null,
                       n_causal=a.n_causal, families=index), fh, indent=2)
    print(f"tier1 family fixtures in {a.out}: {len(fams)} families")


if __name__ == "__main__":
    main()
