# tests/ -- the G-MoQ validation suite

Full documentation: [`../docs/TEST_SUITE.md`](../docs/TEST_SUITE.md). That file
says what each tier establishes, **what it does not establish**, what it costs,
and the pass criteria.

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
`RELEASE_REPORT.md` and `RELEASE_REPORT.json`, with a verdict of **SOUND**,
**SOUND (PARTIAL)** or **NOT SOUND**, every check's statistic and Monte Carlo
error, and an explicit list of what was not established.

Two things worth knowing before editing anything here:

- Every script is standalone (`Rscript --vanilla tests/tierN/xxx.R`) and derives
  its own path from `--file=`, which is the convention the pre-existing
  `stage2/tests/*.R` scripts already use. Sourcing one via `-e source(...)` is
  not supported.
- `tests/criteria.json` holds the preregistered tier 2 pass criteria and
  `tests/tier2/cells.py` the preregistered grid. Both are committed before the
  run they judge. If a band has to move, move it in a commit of its own, with
  the reason.
