# tests/ -- the G-MoQ validation suite

The historical tier scripts are regression checks. The compact paired
known-truth pilot is documented in [`validity_pilot/README.md`](validity_pilot/README.md)
and is the relevant diagnostic for scientific validity. The historical suite
design is retained in [`docs/TEST_SUITE.md`](../docs/TEST_SUITE.md); do not treat
the legacy tier bands below as a scientific validity certificate.

Four tiers, each runnable on its own:

| tier | entry point | cost | needs |
|---|---|---|---|
| 0 unit | `bash tests/run_tier0.sh` | ~70 s | R only |
| 1 integration | `bash tests/run_tier1.sh` | ~1 min | R + the Stage 1 conda env |
| 2 statistical validity | `bash tests/run_tier2.sh submit` | Slurm array | cluster |
| 3 scale and timing | `bash tests/run_tier3.sh` | one large job | cluster + plink2 |

```bash
export GMOQ_WORK=/some/scratch                 # never inside the repo
export GMOQ_DEP_LIB=/path/to/cran/library
export GMOQ_CONDA_SH=$HOME/miniforge3/etc/profile.d/conda.sh
bash tests/run_tier0.sh
Rscript --vanilla tests/report.R $GMOQ_WORK/results
```

`tests/report.R` writes the single artefact a human reads:
`RELEASE_REPORT.md` and `RELEASE_REPORT.json`, with a verdict of
**REGRESSION CHECKS PASSED**, **REGRESSION CHECKS PARTIAL** or
**REGRESSION CHECKS FAILED**, every check's statistic and Monte Carlo error,
and an explicit list of what was not established. These labels describe the
historical checks that ran; they do not certify scientific validity.

Two things worth knowing before editing anything here:

- Every script is standalone (`Rscript --vanilla tests/tierN/xxx.R`) and derives
  its own path from `--file=`, which is the convention the pre-existing
  `stage2/tests/*.R` scripts already use. Sourcing one via `-e source(...)` is
  not supported.
- `tests/criteria.json` holds the preregistered tier 2 pass criteria and
  `tests/tier2/cells.py` the preregistered grid. Both are committed before the
  run they judge. If a band has to move, move it in a commit of its own, with
  the reason.
