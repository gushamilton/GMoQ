# `loco-as-covariate-20260916-01`

Does the REGENIE LOCO prediction belong in the covariate basis of the RIF
columns rather than subtracted from the phenotype? Full statement of the
question, the result and the reasoning: [`docs/LOCO_AS_COVARIATE.md`](../../docs/LOCO_AS_COVARIATE.md).

## What ran

BluePebble compute nodes only, allocation 18979029 on bp1-compute019.
Exclusive scratch root `/user/work/fh6520/gmoq_locovar`; nothing was written into
any existing `gmoq_*` or `G-MoQ`/`Gmoq_paper` tree.

| File | Slurm | What |
|---|---|---|
| `stage1.sh` | 18981251 | Two Stage 1 chr22 scans from this branch, `--loco-mode subtract` and `--loco-mode covariate`, otherwise identical to the paper configuration. Also asserts the subtract arm reproduces the paper reference. |
| `stage2.sh` | 18981301 | Matched weights and Stage 2 (skew, pooled covariance, OLS) for both arms. The covariate arm's Stage 2 baseline is raw `Y`, so its weights are built with the existing no-LOCO weights builder; the subtract arm passes `--loco`. |
| `fourway.R` | interactive | The four-way comparison against sex-stratified REGENIE, plus the SE-swap diagnostic. |
| `diag.R` | interactive | Between-run agreement, and lambda with the signal regions masked. |

Stage 1 ran 229 s (subtract) and 227 s (covariate) at 8 threads, N = 399,710,
17,550 HapMap3 chr22 variants, 45 taus, 1,000 jackknife blocks, seed 42.

## Reproduction check

The subtract arm, built from this branch, is **byte-identical** to
`paper-core-genome-20260915-02` on `stage1.tsv.gz`, `cov.gz` and
`cov.ids.tsv.gz`. Its `sigma0.tsv.gz` differs, which is a pre-existing difference
between the repository head and the paper package and not caused by this change:
every Stage 2 summary in `fourway.tsv` agrees with production run (a) to four
decimals.

The covariate arm's matched weights have baseline sd 1.28347632785109, identical
to the no-LOCO weights, confirming that leaving `Y` unadjusted puts Stage 2 on
the raw-`Y` baseline as intended.

## Results

- `fourway.tsv` — median SE relative to REGENIE, correlation over all variants
  and at REGENIE `|z| >= 4`, and genomic inflation, for `beta_left`,
  `beta_right` and `gamma` in each configuration.
- `se_swap.tsv` — each run's point estimates given the subtraction's standard
  errors, lambda recomputed. This is the diagnostic that identifies the
  subtract-LOCO lambda advantage as SE arithmetic.
- `between_run_and_lambda_masked.tsv` — how close each pair of runs is, lambda
  with 1 Mb around every REGENIE `|z| >= 5` hit removed, and the `gamma`
  shift under covariate mode broken into its estimate and SE parts.

All committed output is aggregate over SNPs. No participant-level data left
BluePebble.

## Load-bearing caveat

All three REGENIE references (`testosterone_{female,male}_testosterone_all` and
`sex_all`) were run with `--ignore-pred` — they carry no LOCO of their own. The
SE and correlation columns therefore measure efficiency and agreement against a
*matching unadjusted* estimator. They are the right yardstick for precision and
for estimand distortion; they cannot certify structure control.
