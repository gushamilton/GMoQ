# REGENIE LOCO predictions are standardised; Stage 1 subtracted them raw

Status: confirmed on real data 2026-09-14. Fixed in Stage 1 by `--loco-scale` (default `sd`).
Pre-v0.4.0 behaviour is `--loco-scale none`.

## What was wrong

REGENIE standardises quantitative phenotypes to unit variance in step 1, so the
per-chromosome predictions it writes to `.loco` are on the **standardised** scale.
Stage 1 computed

    Y_adj = Y - LOCO(chr)

with no rescaling. That is only correct when `sd(Y) == 1`.

## Evidence

Regression of each phenotype on its own chr1 LOCO prediction, UK Biobank, n = 484,460:

| trait      | sd(Y)   | sd(LOCO) | slope(Y ~ LOCO) | slope / sd(Y) | cor   |
|------------|---------|----------|-----------------|---------------|-------|
| bmi        | 4.7862  | 0.3356   | 4.6798          | 0.978         | 0.328 |
| log_bmi    | 0.1665  | 0.3438   | 0.1613          | 0.969         | 0.333 |
| rint_bmi   | 0.9998  | 0.3431   | 0.9663          | 0.967         | 0.332 |
| weight     | 15.9102 | 0.2110   | 12.1951         | 0.766         | 0.162 |

Two things fall out:

1. `sd(LOCO)` sits in 0.21-0.34 for every trait while `sd(Y)` ranges over two orders
   of magnitude (0.17 to 15.9). The predictions carry no memory of the phenotype's
   units - the signature of a standardised predictor.
2. The regression slope tracks `sd(Y)`, at 0.97 of it for the three BMI traits. The
   3% shortfall is REGENIE's own ridge shrinkage. So the factor that puts a `.loco`
   prediction back on the phenotype's scale is `sd(Y)`.

## Consequence, by trait

The sign of the error depends on whether `sd(Y)` is above or below 1.

- **rint_bmi** — `sd(Y) = 1.0000`, so `Y - LOCO` was already right. Correct by accident.
- **bmi, weight** — under-corrected by 4.7x and 12.2x. The subtraction removes ~21% and
  ~8% of the genetic term it should. Effectively these ran *without* LOCO: valid, because
  LOCO is uncorrelated with the test chromosome's genotypes, but less powerful than intended.
- **log_bmi** — over-corrected by 6.2x, and this one is destructive. Subtracting a term of
  sd 0.344 from a phenotype of sd 0.167 does not remove variance, it injects it:

      Var(Y - L) = Var(Y) + Var(L) - 2 Cov(Y, L)
                 = 0.02772 + 0.11820 - 0.03812 = 0.1078,  sd = 0.328

  The analysed variable has ~2x the sd of log(BMI) and is roughly 90% injected LOCO noise.

## Why this matters more for G-MoQ than for a mean GWAS

Subtracting a variable uncorrelated with the test genotype leaves the mean effect
unbiased, so it costs power rather than validity. But G-MoQ's object of study *is* the
shape of the phenotype distribution. Inflating log(BMI)'s variance two-fold compresses
the density at every quantile, which rescales every RIF effect and changes the quantile
grid itself. For `log_bmi` the distribution being analysed is not the log-BMI
distribution.

## The fix

`--loco-scale` on `fungwas_stage1.cli`:

- `sd` (**default**) — multiply by `sd(Y)` on the analysis sample; the exact inverse of
  REGENIE's standardisation, and what REGENIE step 2 does internally.
- `regression` — use the OLS slope of `Y` on the prediction, pooled over chromosomes;
  also absorbs the ~3% ridge shrinkage.
- `none` — subtract as-is; reproduces pre-v0.4.0 runs.
- a number — use it verbatim.

Stage 1 now logs the factor, `sd(Y)`, the median `sd(LOCO)` and the sd of the term it
actually subtracts, and warns when that term exceeds 0.6 `sd(Y)` — the signature of the
bug above.

## Affected outputs

Published `log_bmi` Stage 1 and Stage 2 results, and anything derived from them
(LDSC, lead SNPs), were produced under `none` and should be regenerated.
`bmi` and `weight` are under-powered rather than wrong. `rint_bmi` is unaffected.
LOCOEOF
ssh bp 'bash -lc "cd /user/work/fh6520/GMoQ_loco && git add -A && git -c user.name=\"Fergus Hamilton\" -c user.email=fh6520@bristol.ac.uk commit -q -F - <<EOF
stage1: rescale REGENIE LOCO predictions onto the phenotype scale

REGENIE standardises quantitative phenotypes in step 1, so the predictions
it writes to .loco are on the standardised scale. Stage 1 computed
Y - LOCO(chr) with no rescaling, which is only correct when sd(Y) == 1.

Measured on UK Biobank (n = 484,460), sd(LOCO) is 0.21-0.34 for every trait
while sd(Y) spans 0.17 to 15.9, and the slope of Y on LOCO is 0.97 sd(Y) for
the three BMI traits. The factor is sd(Y).

The sign of the error depends on sd(Y): rint_bmi (sd 1.0) was right by
accident; bmi and weight were under-corrected 4.7x and 12.2x, so LOCO was
close to a no-op; log_bmi (sd 0.167) was over-corrected 6.2x, which raises
the analysed variance from 0.0277 to 0.1078 - the adjusted phenotype is
roughly 90% injected LOCO noise.

Adds --loco-scale {sd,regression,none,<float>}, default sd. none reproduces
the old behaviour. Stage 1 now logs the factor and warns when the subtracted
term exceeds 0.6 sd(Y).

Evidence and per-trait consequences in docs/LOCO_SCALE.md.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
git log --oneline -1"' 2>&1 | tail -5
## Validation against REGENIE step 2

chr22, `log_bmi`, 17,550 variants matched to REGENIE's 17,637 ADD results. REGENIE
step 2 is the ground truth: it consumes its own LOCO predictions correctly. The
G-MoQ quantity compared is the mean of `beta_tau` over the 45-point grid, which
estimates the effect on `E[Y]` and so is directly comparable to REGENIE's `BETA`.

| mode | cor(beta, REGENIE) | median \|beta\| / \|REGENIE beta\| | median SE ratio |
|------|--------------------|-----------------------------------|-----------------|
| `none` (old) | 0.500 | 1.746 | 3.343 |
| `sd` (new default) | **0.993** | **1.007** | 1.534 |
| `regression` | 0.993 | 1.009 | 1.534 |

Under `sd` the tau-averaged effect recovers REGENIE's additive BETA to 0.7%. Under
`none` it is 75% too large and only half-correlated with the truth.

`sd` and `regression` are indistinguishable here: cor(beta) between them is 0.99967
and the median SE ratio 1.00028. `sd` is kept as the default because it is the exact
inverse of REGENIE's standardisation rather than a quantity fitted from the data.

### Correction: the bug distorts the effect estimates, not just their standard errors

An earlier reading of this defect held that point estimates stay unbiased because LOCO
excludes the test chromosome and is therefore uncorrelated with the test genotype, so
only the standard errors inflate. **That is wrong for G-MoQ**, and the chr22 run shows it:

```
cor(mean beta)        none vs sd = 0.4947   median ratio 0.9523
cor(median-tau beta)  none vs sd = 0.3272   median ratio 0.7570
median SE ratio       none / sd  = 2.1784
median-tau SE ratio   none / sd  = 2.3260
```

The unbiasedness argument is a property of linear regression and does not survive the
RIF. `RIF_i(tau) = q_tau + (tau - 1{Y_i <= q_tau}) / f(q_tau)` depends on Y through the
quantile `q_tau` and the density `f(q_tau)`, both of which move when the analysed
variable's variance quadruples. The distortion is worst at the median (cor 0.327), where
the density shifts most.

So the pre-v0.4.0 `log_bmi` results are wrong in their effect estimates as well as
under-powered, and the SEs inflate by 2.2x on top. Regeneration is required, not optional.
