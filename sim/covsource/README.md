# `sim/covsource` — the pooled-vs-jackknife head-to-head

Measures all four Stage 2 inference configurations

    {cov_source = pooled, jackknife}  x  {tau_cov_estimator = ols, gls}

through the **shipped** Stage 1 and Stage 2, at N = 461,320 and the shipped
default `--blocks 1000`, in order to decide what `cov_source` should default to.
The comparison the default actually turns on — pooled+OLS against jackknife+OLS
— had never been run: `docs/FINDINGS_pooled_dispersion.md` has columns for
jackknife+GLS, jackknife+OLS and pooled+GLS and none for pooled+OLS.

Findings: `docs/FINDINGS_covsource_default.md`. Committed output tables:
`sim/covsource/tables/` (the directory is named `tables/` and not `results/`
because the repository `.gitignore` excludes `results/`).

## What runs where

Scripts assume a run directory `WS = /user/work/fh6520/covsource_run` holding
`lib/` (the shipped stage2 installed with `R CMD INSTALL`), `pkgsrc/` (copies of
`stage1/` and `stage2/` from this tree), `cells.json`, `raw/`, `out/`, `logs/`.
Nothing here writes into the git tree.

| file | what it is |
|---|---|
| `setup_cs.sh` | installs the shipped stage2 into `WS/lib`, times Stage 1 at N = 461,320 |
| `truth_cs.R` | EXACT population truth: `beta_tau`, `Sigma_0` (what pooled targets) and `Sigma_true` (what the jackknife targets), plus population theta and SE for every configuration. No Monte Carlo anywhere. |
| `make_cells.R` | builds `cells.json`. Effect sizes are CALIBRATED by root-finding on the exact true Z, not guessed. |
| `exact_cs.R` | the exact pooled surface: SE ratio, coverage, Type-I at 0.05 and at 5e-8. The jackknife has no closed form and deliberately does not appear here. |
| `sim_stage1_cs.py` | one replicate = one whole scan. Drives shipped `core.run_stage1` (betas, se, d_resid and the **jackknife covariance**) and shipped `pooled_cov.compute_sigma0`. The jackknife covariance is stored as a float32 row-major upper triangle: the exact on-disk layout and precision of a production `.cov.gz`. |
| `stage2_driver_cs.R` | runs the shipped `param_gwas()` over the full 2x2 on the SAME replicates, so every pooled-vs-jackknife comparison is PAIRED. Pooled `cov_data` is assembled by the package's own `.pooled_cov_data()`; no estimator is reimplemented. |
| `summarise_cs.R` | counts coverage, Type-I, SE ratio, power at matched ACTUAL Type-I, and paired differences, all with Monte Carlo error. |
| `tables_cs.R` | renders the markdown tables in the findings document. |
| `run_stage1_cs.sh`, `run_stage2_cs.sh` | Slurm arrays. |

## Reproduce

    sbatch setup_cs.sh                                  # install + time
    Rscript --vanilla make_cells.R                      # writes cells.json
    Rscript --vanilla exact_cs.R                        # exact table, no cluster needed
    sbatch --array=0-1039%640 run_stage1_cs.sh          # 52 cells x 20 chunks x 750 reps
    sbatch --dependency=afterok:<id> --array=0-51%22 run_stage2_cs.sh
    Rscript --vanilla summarise_cs.R
    Rscript --vanilla tables_cs.R > tables.md

## Two traps this code exists to avoid

* `~/.Rprofile` carries a bare `.libPaths(...)` that silently overrides
  `R_LIBS_USER`. Every R script here runs under `--vanilla` **and** asserts
  `.libPaths()[1]`.
* Degrees of freedom are not interchangeable between the two covariance sources.
  The pooled path is run with `nu = n_samples`, `nu_source =
  "pooled_sigma0_n_samples"`; the jackknife path with `nu = NULL`, `n_blocks =
  1000`, `n_eff = d_resid`, which is rung 2 of `.param_reference()` (the per-SNP
  Satterthwaite law). Running the jackknife against the pooled `nu` would
  flatter it.
