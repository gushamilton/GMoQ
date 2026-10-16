# tests/ -- tier 0

Tier 0 is the fast cross-package unit suite: the `tests/tier0/*.R` checks, the
Stage 2 assertion scripts and the testthat suite, against a freshly installed
Stage 2. About a minute; R only.

```bash
export GMOQ_WORK=/some/scratch          # never inside the repo
export GMOQ_DEP_LIB=/path/to/cran/library   # optional: where the CRAN dependencies live
bash tests/run_tier0.sh
Rscript --vanilla tests/report.R $GMOQ_WORK/results
```

`tests/report.R` writes `RELEASE_REPORT.md` and `RELEASE_REPORT.json`. Their
verdict describes the regression checks that ran; it does not certify
scientific validity. For recovery of known effects on real genotypes, run
`docs/examples/1kg_chr22/`.

Every script is standalone (`Rscript --vanilla tests/tier0/xxx.R`) and derives
its own path from `--file=`; sourcing one via `-e source(...)` is not supported.
`tests/criteria.json` holds the preregistered pass criteria.

The tier 1-4 simulation suites (integration, statistical validity on a Slurm
grid, scale on UK Biobank chr22) and the known-truth validity pilot are at the
tag `archive/research-records-20261016`.
