# REGENIE LOCO predictions are standardised; Stage 1 subtracted them raw

Status: confirmed against REGENIE v4.1.2 source and on real data, 2026-09-14.
Fixed in Stage 1 by `--loco-scale` (default `regression`). Pre-v0.4.0 behaviour is
`--loco-scale none`.

## What was wrong

REGENIE residualises quantitative phenotypes on the step-1 covariates and then divides
by the residual sd (`src/Pheno.cpp::residualize_phenotypes`, the step that logs
`-residualizing and scaling phenotypes...done`):

```cpp
pheno_data->phenotypes -= ((pheno_data->new_cov * beta.transpose()).array() * masked.array()).matrix();
pheno_data->scale_Y = pheno_data->phenotypes.colwise().norm().array() / sqrt(Neff - ncov);
pheno_data->phenotypes.array().rowwise() /= pheno_data->scale_Y.array();
```

`src/Data.cpp::write_predictions` writes `.loco` straight from `predictions[0]` and never
multiplies back by `scale_Y`, so the file is on the **standardised** scale. Step 2 is
consistent with this: it residualises and scales again before `res = phenotypes - blups`,
then converts betas back to original units at `Step2_Models.cpp:173`.

Stage 1 computed

    Y_adj = Y - LOCO(chr)

with no rescaling. That is only correct when the standardising constant is 1.

## Evidence

**sd(LOCO) is invariant to the phenotype's units.** chr22, each trait on its own keep set:

| trait | n | sd(Y) | sd(Y \| step-1 covars) | cov R² | sd(LOCO) | cor |
|---|---|---|---|---|---|---|
| bmi | 484,460 | 4.78617 | 4.74888 | 0.016 | 0.34987 | 0.343 |
| log_bmi | 484,460 | 0.16655 | 0.16483 | 0.021 | 0.35850 | 0.349 |
| rint_bmi | 484,460 | 0.99983 | 0.98729 | 0.025 | 0.35785 | 0.348 |
| weight | 484,460 | 15.91019 | 13.32596 | 0.298 | 0.23332 | 0.199 |
| vitamin_d | 442,632 | 21.11313 | 20.63443 | 0.045 | 0.24368 | 0.234 |
| testosterone | 419,890 | 1.28313 | 1.28213 | 0.002 | 0.18718 | 0.159 |

`sd(Y)` spans 0.167 to 21.1 — a 127-fold range — while `sd(LOCO)` stays within 0.187 to
0.359. The predictions carry no memory of the units.

**Ground truth from REGENIE step 2's own standard errors.** Inverting
`SE = sd(analysed residual) / (sd(G) sqrt(N-p))` with `var(G) = INFO * 2f(1-f)` over
chr22 SNPs with INFO > 0.9 and MAF > 0.05, on REGENIE's step-2 sample (n = 461,320):

| trait | REGENIE implied sd | c=0 | c=1 (the bug) | c=sd(Y) | c=sd(Y\|cov) | c=regression |
|---|---|---|---|---|---|---|
| log_bmi | 0.154269 | 1.0684 | **2.2065** | 0.9998 | 0.9998 | 0.9997 |
| bmi | 4.451434 | 1.0660 | 1.0412 | — | 0.9998 | — |
| weight | 13.057705 | 1.0184 | 1.0156 | 1.0017 | 1.0003 | 1.0002 |

(entries are ratio to REGENIE's implied sd). The back-out reproduces REGENIE to 0.02%.

## Consequence, by trait

- **rint_bmi** — right by accident: `sd(Y) = 0.99983`. Strictly the correct factor is the
  residual sd 0.98729, so it was 1.3% under-corrected — a 0.013% variance difference.
- **bmi, weight** — under-corrected. Of the variance reduction that correct scaling would
  achieve, `c=1` delivers **38%** for bmi (22.519 → 21.480 against a target of 19.807) and
  **15.6%** for weight (176.84 → 175.87 against 170.59). Shape effects attenuate by only
  4% and 1.5%.
- **log_bmi** — over-corrected 6.2x, and destructive. Variance goes 0.027168 → 0.115871
  (sd 0.16483 → 0.34040); the LOCO-derived share of the analysed variance is **76.6%**.

Note on "under-corrected is harmless": LOCO's job in UK Biobank is controlling relatedness
and population structure, not power. Delivering 16-38% of the intended correction means the
test statistics are correspondingly under-protected against that inflation. That is a
validity question, not only a power one.

## The bug distorts effect estimates, not just their standard errors

The mean RIF is `Y` itself, so the tau-averaged effect is unbiased:
`beta(Y - cL) = beta(Y) - c*beta(L~G) = beta(Y)`, since LOCO excludes the test chromosome.
**That argument does not extend to the quantile effects.** Adding a term independent of G
convolves the conditional distribution, and any scale or shape effect attenuates by
`sd_correct / sd_analysed`. Simulation at N = 300,000 x 20 replicates with injected noise
of 1.8x:

```
pure location shift (b = 0.20)
  tau      0.10   0.25   0.50   0.75   0.90
  ratio   1.019  1.012  0.991  0.988  0.972   mean 0.995   <- unbiased
pure variance effect (log-sd slope 0.20)
  ratio   0.522  0.538     -    0.536  0.506   theory 1/sqrt(1+1.8^2) = 0.486
```

For log_bmi the attenuation is `0.154234 / 0.340398 = 0.453`: dispersion effects are more
than halved *and* their SEs inflate 2.21x, so z on a shape effect drops by about 4.9x.

Confirmed end-to-end on chr22 against REGENIE step 2 (17,550 variants matched):

| mode | cor(beta, REGENIE) | median \|beta\| / \|REGENIE beta\| |
|------|--------------------|-----------------------------------|
| `none` (old) | 0.500 | 1.746 |
| `sd` | 0.993 | 1.007 |
| `regression` (default) | 0.993 | 1.009 |

Between modes: `cor(mean beta)` none vs regression is 0.4947, and 0.3272 at the median tau.
The distortion is worst at the median, where the density moves most.

## The fix

`--loco-scale` on `fungwas_stage1.cli`:

- `regression` (**default**) — OLS slope of `Y` on the prediction, pooled over chromosomes.
- `sd` — multiply by `sd(Y)`.
- `none` — subtract as-is; reproduces pre-v0.4.0 runs bit-for-bit.
- a number — use it verbatim.

**Why `regression` and not `sd`.** The exact factor is REGENIE's `scale_Y`, the sd of the
trait *after* its step-1 covariates. Stage 1 cannot compute that: REGENIE's step-1 covariate
set is not an input and need not match Stage 1's `--covar` — for `weight` it does not
(sex, age, height, height_z2, log_height, height_z3, R² = 0.298), and there `sd(Y)`
over-subtracts by 19.4%. `regression` needs no knowledge of that set, is the
variance-minimising factor by construction, and matched REGENIE's analysed sd to within
0.02% on every trait that could be ground-truthed.

`regression` is not equal to `scale_Y`: the slope is `scale_Y * rho / sd(pred)`, which
coincides only for a perfectly calibrated prediction. REGENIE's are mildly over-dispersed,
so the slope sits below `scale_Y` — by 11% for weight, 3% for the BMI traits. Near the
variance minimum this costs ~0.1%, and it also absorbs REGENIE's ridge shrinkage. Earlier
notes here attributed the whole `weight` anomaly to shrinkage; most of it is the covariate
R², with shrinkage tracking the predictor's strength (cor 0.35 → 0.97; cor 0.16-0.20 → 0.85-0.89).

Stage 1 logs the factor and the variance before and after, and warns when the adjustment
*increases* the analysed variance — measured directly rather than inferred from an sd
ratio, so a trait with predictive covariates cannot trip it spuriously.

## Affected outputs

Published `log_bmi` Stage 1 and Stage 2 results, and anything derived from them (LDSC,
lead SNPs), were produced under `none`: effect estimates distorted, shape effects
attenuated 2.2x, SEs inflated 2.2x. They require regeneration. `bmi` and `weight` received
38% and 16% of the intended correction — under-powered and under-protected against
structure, not wrong in their estimates. `rint_bmi` is unaffected.

The paper pipeline invokes the affected path at
`src/paper/06_run_fungwas_stage1_loco_chr.sh:89`.
