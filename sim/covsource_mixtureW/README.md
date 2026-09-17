# `sim/covsource_mixtureW` — does the pooled-covariance defect transfer to the W the paper ships?

> `sim/covsource`, referred to throughout below, was archived off `main` on
> 2026-09-17 and is frozen at tag `archive/sim-campaigns-20260917`. Its
> aggregates also live in `Gmoq_paper/outputs/simulations/biobank_covariance/full_audit/`.
> The `patch_*.py` scripts here read their anchors from the BluePebble run root
> `/user/work/fh6520/covsource_run/src/`, not from the git tree, so nothing here
> stopped working.

`sim/covsource` decided `cov_source` using the **2-column** `make_weights_vqtl`
W (`beta_mu`, `beta_sigma2`). The paper ships a **3-column** mixture/skew W from
`build_skew_mixture()` (`gamma`, `beta_left`, `beta_right`). `docs/EFFICIENCY_BOUND.md`
(register `projection-efficiency-20260916-01`) showed that the membership column
reshapes how `W` meets `Sigma`: at testosterone separation `col(W)` is
`Sigma`-invariant. That is a structural difference, so the published
pooled-vs-jackknife figures could not be assumed to transfer in either direction.

This campaign measures both W families on a common DGP.

Findings: `docs/FINDINGS_covsource_mixtureW.md`. Tables: `tables/`.

## What runs where

Run directory `WS = /user/work/fh6520/gmoq_covmixW`, holding `lib/` (the shipped
Stage 2, copied from the `sim/covsource` run root), `pkgsrc/`, `cells_mix.json`,
`raw/`, `out/`, `logs/`. Nothing here writes into the git tree.

| file | what it is |
|---|---|
| `truth_mix.R` | Exact population truth, extending `sim/covsource/truth_cs.R` to a bimodal `eps` (`"mix<S>"`, `"mix<S>p<p1>"`) and to the 3-column W. Proves and uses the identity `Sigma_true - Sigma_pooled = diag(1/f) E diag(1/f)`. |
| `verify_W.R` | The seven controls that must pass before any number is drawn (see below). |
| `exact_mix.R` | The exact pooled surface — SE ratio, coverage, Type-I at 0.05 and 5e-8 — for `vqtl2` / `mix2` / `mix3` over 158 cells. No Monte Carlo error anywhere. |
| `mechanism.R` | Why. `diag(E) = Delta(tau)(1-2tau)`; the eigenvalue agreement; the `p1 != 0.5` attempt to break it. |
| `mixdgp.R` | The DGPs the 3-column W exists for: component-specific and **membership** effects, where the null parameter is a component mean rather than `gamma`. |
| `robust.R` | Does the protection survive a W built from mis-fitted mixture parameters (+/- 20% on each)? |
| `onset.R` | The `|Z|` at which pooled first falls 1/5/10% short, and the chr22-testosterone operating point. |
| `make_cells_mix.R` | The Monte Carlo grid. Effect sizes are CALIBRATED by root-finding on the exact true Z of the W that SHIPS for that phenotype shape, not guessed. |
| `patch_eps.py` | Teaches `sim/covsource/sim_stage1_cs.py`'s `draw_eps` the standardised 2-normal mixture law. Nothing else in Stage 1 changes. |
| `patch_driver.py` | Derives the Stage 2 driver from `sim/covsource/stage2_driver_cs.R` by five localised substitutions. The shipped `param_gwas()`, the package's own `.pooled_cov_data()` and both df conventions are untouched. |
| `patch_summ.py` | Derives the summariser from `sim/covsource/summarise_cs.R`. Only paths and the power block change; every rate, its MC error, the paired differences and the 5e-8 tail extrapolation are the published code verbatim. |
| `run_stage1_mix.sh`, `run_stage2_mix.sh` | Slurm arrays. |
| `tables_mix.R` | Renders the markdown tables in the findings document. |

The three `patch_*.py` scripts exist instead of forked copies so that the
derivation from the published campaign is auditable line by line: each asserts
its anchor appears exactly once before substituting.

## Controls (`verify_W.R`, all must PASS)

1. The new `pop_truth()` reproduces the published one bit-for-bit at `shape="normal"` (rel err 1e-16).
2. The bimodal `eps` law is self-consistent: mean 0, var 1, `pdf == d/de cdf`.
3. **`make_weights_normal_mixture(include_membership=TRUE)` equals `build_skew_mixture()`'s W at `alpha = 0`** (max rel diff 3.9e-10), with column names `gamma, beta_left, beta_right`. This is what licenses using the package function in place of the production builder.
4. Column counts are asserted: `vqtl2` K=2, `mix2` K=2, `mix3` K=3.
5. Under the null `Delta(tau) == 0` and `Sigma_0 == Sigma_true` exactly — the reason null simulations cannot see this defect.
6. `B_w - B = Delta(min) - t_s Delta_t - t_t Delta_s` holds to 4e-16. The decomposition is **exact**, not first-order.
7. The `projeff` geometry is reproduced: `angle(W_gamma, v1)` = 0.019 deg at S=5.89, leading eigenvalue share 0.980.

## Reproduce

    Rscript --vanilla src/verify_W.R                     # controls first
    Rscript --vanilla src/exact_mix.R                    # exact surface, no cluster
    Rscript --vanilla src/mechanism.R
    Rscript --vanilla src/mixdgp.R
    Rscript --vanilla src/make_cells_mix.R               # writes cells_mix.json
    sbatch --array=0-189%120 src/run_stage1_mix.sh       # 19 cells x 10 chunks x 500 reps
    sbatch --dependency=afterok:<id> --array=0-18 src/run_stage2_mix.sh
    Rscript --vanilla src/summarise_mix.R
    Rscript --vanilla src/tables_mix.R
