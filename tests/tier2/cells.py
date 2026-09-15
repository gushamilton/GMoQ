#!/usr/bin/env python3
"""The tier 2 cell grid. Written down once, here, and frozen.

    python tests/tier2/cells.py tests/tier2/cells.json

WHY A FILE. Tier 2's numbers are what decide whether a release is sound, so the
design has to be fixed before the numbers are seen. The grid and the pass
criteria (tests/criteria.json) are both committed, both read by the collator,
and neither is computed from a run.

THE GRID.

Main grid, 28 cells, at N = 50,000 on the production 45-tau grid:

    {null, location, scale, locscale, skew, kurtosis} x {MAF 0.30, 0.01}
                                                      x {normal, skewed}   = 24
    {mix_weight, mix_mean}                            x {MAF 0.30, 0.01}   =  4

Each cell carries 190 null variants and 10 causal ones per replicate, so one
cell at 300 replicates gives 57,000 null observations (Monte Carlo error on a
0.05 rate: 0.0009) and 3,000 causal ones (error on a 0.95 coverage: 0.004).

Ten causal variants, not thirty, and |Z| = 8, not 12. Thirty variants at |Z| = 12
in one phenotype make the causal effects a large part of the trait's own
variance, which is not the regime a GWAS is ever in and which pushes the
estimand far enough from linear that the departure, rather than the estimator,
dominates every coverage number.

Every cell is analysed under four estimator arms at no extra Stage 1 cost:
{jackknife, pooled} x {ols, gls}.

Pin cells. These reproduce established failure modes; the suite fails if any
of them stops being reproducible, because that would mean the thing the guards
guard against has changed shape:

    pin_vqtl_maf0.01 / 0.30 / 0.50   a PURE LOCATION effect, large, at three
        genotype skewnesses. The false-dispersion residue is quadratic in the
        mean effect and proportional to genotype skewness, so it must be
        violent at MAF 0.01, present at 0.30, and EXACTLY ABSENT at 0.50 --
        and the shipped G guard must separate the three.
        (docs/FINDINGS_pooled_dispersion.md)

    pin_ties      a phenotype rounded hard enough that two taus share a sample
        quantile. GLS must still collapse and OLS must still be unaffected.
        (docs/BP2_lumpy_phenotypes.md section 2)

    pin_ceiling   a ceiling-censored phenotype, the ONE anti-conservative
        failure found in BP2: GLS manufactures far-tail false positives.
        Run at 1,000 replicates because the statistic of interest is a
        rejection rate at 1e-4.
"""
from __future__ import annotations

import json
import sys

PROD_TAUS_N = 45
N_MAIN = 50_000
REPS_MAIN = 300
CHUNK_MAIN = 30

CF_FAMILIES = ["null", "location", "scale", "locscale", "skew", "kurtosis"]
MIX_FAMILIES = ["mix_weight", "mix_mean"]
MAFS = [0.30, 0.01]
BASELINES = {"normal": 0.0, "skewed": 0.8}      # baseline sinh-arcsinh skewness


def cell(name, **kw):
    base = dict(
        name=name, family="null", maf=0.30, N=N_MAIN, n_taus=PROD_TAUS_N,
        blocks=1000, n_null=190, n_causal=10, z_target=8.0,
        baseline_skew=0.0, baseline_kurt=0.0, discretise=None,
        reps=REPS_MAIN, chunk_reps=CHUNK_MAIN, n_covar=2,
        fit_models=["LS"], kind="main")
    base.update(kw)
    return base


def build():
    cells = []
    for fam in CF_FAMILIES:
        for maf in MAFS:
            for bname, bskew in BASELINES.items():
                cells.append(cell(
                    f"{fam}_maf{maf:g}_{bname}", family=fam, maf=maf,
                    baseline_skew=bskew,
                    n_causal=0 if fam == "null" else 10,
                    n_null=200 if fam == "null" else 190,
                    fit_models=["LS", "TRUE"]))
    for fam in MIX_FAMILIES:
        for maf in MAFS:
            cells.append(cell(f"{fam}_maf{maf:g}", family=fam, maf=maf,
                              fit_models=["LS", "TRUE"]))

    # --- pins -------------------------------------------------------------
    # A pure location effect at three genotype skewnesses. z_target is chosen
    # so the guard statistic G = 0.6855 gskew Z^2 / sqrt(N) lands well above
    # its 0.2 threshold at MAF 0.01 and 0.30, and is exactly zero at 0.50.
    for maf, z in ((0.01, 8.0), (0.30, 30.0), (0.50, 60.0)):
        cells.append(cell(f"pin_vqtl_maf{maf:g}", kind="pin_vqtl",
                          family="location", maf=maf, z_target=z,
                          n_null=170, n_causal=30, fit_models=["LS"]))

    # Ties and ceiling censoring, on the 15-tau grid and N of
    # docs/BP2_lumpy_phenotypes.md so the published numbers are comparable.
    cells.append(cell("pin_ties", kind="pin_ties", family="null", maf=0.30,
                      N=20_000, n_taus=15, blocks=1000, n_null=200, n_causal=0,
                      discretise={"kind": "round", "h": 0.5},
                      reps=400, chunk_reps=50, fit_models=["LS"]))
    cells.append(cell("pin_ceiling", kind="pin_ceiling", family="null", maf=0.30,
                      N=20_000, n_taus=15, blocks=1000, n_null=200, n_causal=0,
                      discretise={"kind": "ceiling", "tau": 0.80},
                      reps=1000, chunk_reps=50, fit_models=["LS"]))
    return cells


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "cells.json"
    cells = build()
    tasks = []
    for i, c in enumerate(cells):
        n_chunks = (c["reps"] + c["chunk_reps"] - 1) // c["chunk_reps"]
        for k in range(n_chunks):
            tasks.append(dict(cell=i, chunk=k))
    with open(out, "w") as fh:
        json.dump(dict(cells=cells, tasks=tasks,
                       n_cells=len(cells), n_tasks=len(tasks)), fh, indent=2)
    print(f"{len(cells)} cells, {len(tasks)} array tasks -> {out}")


if __name__ == "__main__":
    main()
