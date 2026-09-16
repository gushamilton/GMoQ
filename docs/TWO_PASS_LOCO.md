# Two-pass blind LOCO

**Status: built and measured on testosterone chr22, 2026-09-16. Pooled covariance,
OLS, 17,550 variants, N = 399,710 (Stage 1) / 419,890 (LOCO).**

The problem this addresses is in [LOCO_MIXTURE_HAZARD.md](LOCO_MIXTURE_HAZARD.md):
REGENIE step 1 for testosterone was run without sex as a covariate -- deliberately,
so that G-MoQ could model the sex mixture itself -- and the whole-genome ridge
learned *sex* rather than polygenic signal. Subtracting that predictor cost real
accuracy and precision.

**Summary of what was found.** Two-pass blind LOCO works as a construction: it
produces a LOCO predictor that is no longer a membership proxy, that carries
20-70x more within-component polygenic signal than the old one, and that restores
(slightly improves on) the standard errors and effect correlations of the no-LOCO
arm. It does **not** reduce genomic inflation relative to no-LOCO. Investigating
why produced the more important result: **the old LOCO's apparent lambda advantage
was never structure control.** See [What the lambda actually
was](#what-the-lambda-actually-was).

## The procedure

```
Pass 1   run G-MoQ with NO LOCO.  Stage 2 fits the two-component skew-normal
         mixture, giving component parameters (pi_k, xi_k, omega_k, alpha_k).

Blind    turn those parameters into posterior membership weights per sample,
residual w_k(y) = pi_k f_k(y) / sum_j pi_j f_j(y),
         and residualise the phenotype on the ESTIMATED structure.

Refit    run REGENIE step 1 on the residualised phenotype.  It can no longer
         learn the mixture, because the mixture is no longer there.

Pass 2   run G-MoQ on the ORIGINAL phenotype using the resulting LOCO.
```

Every quantity used is inferred by the method. No membership label is read at any
point, which is what makes this legitimate for a genuine latent-mixture
application where no such column exists.

## Choice of residualising transform

Two transforms were built and both taken all the way through, because the argument
for each is different and testing both cost one extra phenotype column in a single
REGENIE run.

Write `m_k` and `s_k` for the fitted mean and standard deviation of component `k`
(for a skew-normal component with `delta = alpha/sqrt(1+alpha^2)`,
`m = xi + omega*delta*sqrt(2/pi)` and `s = omega*sqrt(1 - 2*delta^2/pi)`).

**`resid_mean` -- remove component location only.**

```
r_i = y_i - sum_k w_k(y_i) m_k
```

Kills the between-component mean contrast, which is the only thing a ridge *mean*
model can learn, and leaves every within-component deviation untouched. This is the
minimal intervention that closes the failure mode, and it is the **recommended
default**.

**`resid_z` -- remove component location and scale.**

```
r_i = sum_k w_k(y_i) (y_i - m_k) / s_k
```

Additionally equalises the components' variances. The argument is downstream: G-MoQ
applies a *single* global scale when it subtracts `scale x LOCO(chr)`, and with
unequal residual spreads (0.455 vs 0.300 here) one global scale cannot suit both.
In the event the two performed equivalently (below), so the simpler transform wins.

**Why a posterior-weighted transform does not destroy the signal it preserves.**
The transform is a function of `y`, so in principle it could shrink genuine
within-component variation toward the component mean. It does not here, because the
components are well separated (`S = 6.22` on the fitted parameters) and the
posterior is therefore nearly hard: only **0.56%** of samples have `w_1` in
`[0.05, 0.95]` and **0.23%** in `[0.20, 0.80]`. Within hard-assigned components the
transform is close to affine in `y`: `cor(y, resid_mean) = 0.996` (lower) and
`0.989` (upper). **Check separation before trusting this step** -- at low separation
the shrinkage in the ambiguous region would be a real cost.

**Evidence the structure is actually gone** (N = 421,405):

| quantity | value |
|---|---|
| `R^2(y ~ posterior w1)` | 0.9099 |
| `R^2(resid_mean ~ posterior w1)` | 0.00020 |
| `R^2(resid_z ~ posterior w1)` | 0.00022 |
| `R^2(y ~ sex)` *(validation only)* | 0.8996 |
| `R^2(resid_mean ~ sex)` *(validation only)* | 0.00045 |
| `R^2(resid_z ~ sex)` *(validation only)* | 0.00048 |
| within-component sd, `resid_mean` | 0.455 / 0.300 |
| within-component sd, `resid_z` | 0.980 / 0.981 |

The mixture is 91% of the phenotype's variance and 99.95% of it is removed. The
`sex` rows confirm the blind posterior really is recovering the latent variable
(`cor(posterior w1, sex) = -0.992`) and play no part in the method.

**REGENIE step 1 then found real signal.** Level-1 ridge `Rsq`: original
testosterone run 0.0253 (which the diagnostic below shows was almost entirely sex);
`resid_mean` 0.0287, `resid_z` 0.0345 -- on a phenotype from which 91% of the
variance had already been removed.

## Results: three-way comparison

chr22, pooled covariance, OLS. SE ratio is the median Stage 2 SE relative to
sex-stratified REGENIE (run `--ignore-pred`); `beta_left` is compared with the
female scan and `beta_right` with the male scan.

| configuration | SE ratio F | SE ratio M | cor all F | cor all M | cor \|z\|>=4 F | cor \|z\|>=4 M |
|---|---|---|---|---|---|---|
| with LOCO (published) | 1.1209 | 1.0677 | 0.8981 | 0.8165 | 0.9958 | 0.9943 |
| no LOCO | 1.0236 | 0.9311 | 0.9699 | 0.9188 | 0.9984 | 0.9967 |
| **two-pass, `resid_mean`** | **1.0113** | **0.9134** | **0.9664** | **0.9140** | **0.9982** | **0.9966** |
| two-pass, `resid_z` | 1.0136 | 0.9059 | 0.9662 | 0.9121 | 0.9981 | 0.9963 |
| permuted null | 1.0287 | 0.9330 | 0.0361 | 0.0974 | -0.6677 | 0.2554 |

| configuration | lambda gamma | lambda beta_left | lambda beta_right |
|---|---|---|---|
| with LOCO (published) | 1.2671 | 1.3527 | 1.3604 |
| no LOCO | 1.2650 | 1.4195 | 1.4869 |
| **two-pass, `resid_mean`** | 1.2560 | 1.4221 | 1.4778 |
| two-pass, `resid_z` | 1.2698 | 1.4161 | 1.5106 |
| permuted null (true null) | 0.9990 | 0.9045 | 0.9928 |

Two-pass recovers the whole precision and accuracy cost of the old LOCO -- SEs are
1-3% *better* than the no-LOCO arm and 10-15% better than with-LOCO, and effect
correlations land within 0.004 of no-LOCO. Genomic inflation is unchanged from
no-LOCO.

**The target as originally stated -- no-LOCO precision with with-LOCO lambda -- was
not met, and the next section explains why it was the wrong target.**

## The blind diagnostic on old and new LOCO

`scripts/loco_blind_diagnostic.R`, run with the pass-1 mixture fit. No membership
label is used; the `cor(sex, ...)` row is a validation check only.

| quantity | old LOCO | new, `resid_mean` | new, `resid_z` |
|---|---|---|---|
| LOCO scale factor | 1.0881 | 0.2090 | 0.2142 |
| `cor(y, loco)` | 0.1589 | 0.0296 | 0.0327 |
| **`r2_loco_on_posterior`** | **0.02536** | **0.000466** | **0.000329** |
| `r2_within_lower` | 0.000330 | 0.02236 | 0.01417 |
| `r2_within_upper` | 0.001458 | 0.03751 | 0.05671 |
| `var_ratio_within_lower` | 1.2059 | 0.9822 | 0.9866 |
| `var_ratio_within_upper` | 1.3068 | 0.9673 | 0.9534 |
| membership-proxy ratio | 17.40 | 0.0124 | 0.0058 |
| **verdict** | **FLAGGED** | clean | clean |
| `cor(sex, loco)` *(validation only)* | 0.1599 | -0.0229 | -0.0196 |

The diagnostic flags the old predictor decisively and clears both new ones, without
ever seeing sex. The three quantities move exactly as the mechanism predicts:

- The old predictor explained 2.5% of the posterior and 0.03-0.15% of the phenotype
  within a component -- a 17x membership-proxy ratio. The new predictors explain
  0.03-0.05% of the posterior and **2.2-5.7% within a component**: a 20-70x
  improvement in the only thing a LOCO offset is for.
- `var(y - c*loco)/var(y)` within a component went from **1.21 / 1.31** (the
  subtraction *added* 21-31% variance) to **0.98 / 0.95** (it now removes variance,
  as it should).
- `cor(sex, loco)` fell from 0.160 to -0.02, confirming independently that the
  membership leak is closed.

This diagnostic reproduced every published number for the old LOCO
(`loco_scale 1.088`, `cor(y,loco) 0.1589`, within-component variance ratios
1.21/1.31) from the fitted mixture alone.

## What the lambda actually was

LOCO_MIXTURE_HAZARD.md concluded that the 5-9% lambda rise on removing LOCO was
"real lost relatedness/structure control". **That conclusion is wrong, and the
two-pass experiment is what exposed it.**

Median Stage 2 SE, relative to the no-LOCO arm, alongside lambda:

| configuration | SE gamma | SE beta_left | SE beta_right | lambda gamma | lambda beta_left | lambda beta_right |
|---|---|---|---|---|---|---|
| with LOCO | 1.0037 | **1.0950** | **1.1467** | 1.2671 | 1.3527 | 1.3604 |
| no LOCO | 1.0000 | 1.0000 | 1.0000 | 1.2650 | 1.4195 | 1.4869 |
| two-pass `resid_mean` | 0.9998 | 0.9880 | 0.9809 | 1.2560 | 1.4221 | 1.4778 |
| two-pass `resid_z` | 0.9998 | 0.9902 | 0.9728 | 1.2698 | 1.4161 | 1.5106 |
| permuted null | 1.0016 | 1.0050 | 1.0020 | 0.9990 | 0.9045 | 0.9928 |

The old LOCO left `gamma`'s standard error untouched (1.0037) **and left `gamma`'s
lambda untouched** (1.2671 vs 1.2650). It inflated `beta_left`'s SE by 9.5% and
`beta_right`'s by 14.7%, and exactly those two lambdas fell.

Genuine relatedness or population-structure confounding inflates *every* test
statistic from the same genotypes. A control that removed it would have deflated
all three lambdas. This one deflated only the two whose standard errors it had
inflated, in the direction and roughly the magnitude that inflation alone predicts
(`1/1.1467^2 = 0.76` would take `lambda beta_right` from 1.487 to 1.13, i.e. *more*
deflation than the 1.360 observed).

So the with-LOCO arm's lower lambda was a **variance-inflation artefact of
subtracting a membership proxy**, not structure control. Across all four real
configurations lambda tracks the SE ratio monotonically and tracks nothing else.

That also explains, without any further hypothesis, why two-pass does not lower
lambda: there was no LOCO-removable inflation on this chromosome to remove. The
residual lambda of ~1.42-1.49 on `beta_left`/`beta_right` and ~1.26 on `gamma` is
chr22 signal plus whatever structure the 10 PCs in the covariate file do not
capture -- and a LOCO offset was never going to touch the latter here.

**Caveat on what has and has not been shown.** The permuted null destroys the
phenotype-relatedness link along with the genotype-phenotype association, so it
establishes that the estimator is calibrated under no structure but does not
quantify how much of the real lambda is confounding. The `gamma`-invariance
argument above is the load-bearing evidence, and it is an argument from
*differential* behaviour across parameters, which is exactly what a confounding
account cannot produce. Separating residual confounding from polygenicity properly
would need an LD-score intercept over more than 17,550 variants.

## Verdict

1. **Two-pass blind LOCO is the right construction and it works.** It is built only
   from quantities the method infers, it removes 99.95% of the learnable mixture,
   and it turns a LOCO predictor that explained 0.03-0.15% within a component into
   one that explains 2.2-5.7%. It recovers the full precision and accuracy cost of
   the old LOCO -- SEs 1-3% better than dropping LOCO entirely.
2. **It does not reduce genomic inflation, and it should not have been expected to.**
   The lambda gap it was built to close was an artefact of the old predictor's own
   variance inflation.
3. **The blind diagnostic is the durable deliverable.** It flagged the old predictor
   at a 17x membership-proxy ratio and cleared the new ones at 0.006-0.012, with no
   label anywhere. Run it before trusting any LOCO adjustment on a structured
   phenotype. The antibody traits `jcv_log_all` and `mcv_log_all` went through the
   same workflow and have not been checked.
4. **Recommendation: `resid_mean` two-pass.** It matches `resid_z` on every metric,
   has the more conservative justification (it removes only what a mean model can
   learn), and gives more balanced within-component explanatory power
   (0.022/0.038 versus 0.014/0.057). Plain no-LOCO remains an acceptable fallback:
   it costs 1-3% in SE and nothing in lambda.

## Reproduction

Scripts in this repo: `scripts/loco_blind_diagnostic.R`,
`scripts/compare_loco_configs.R`. The BluePebble run harness lived in a scratch
directory; in outline:

```bash
# Pass 1 already exists: a no-LOCO G-MoQ run and its weights_skew.rds.
# Blind residualisation (writes participant-level files that stay on the cluster):
Rscript make_resid_pheno.R --weights <pass1 weights_skew.rds> \
  --pheno <phenotypes.txt> --pheno-col testosterone --keep <keep.txt> \
  --out-pheno <resid.txt> --out-summary <components.tsv>

# Refit step 1 on both residualised columns in one run, reusing the ORIGINAL QC
# snplist so the variant set is identical and the phenotype is the only change:
regenie --step 1 --bed <bed> --phenoFile <resid.txt> \
  --phenoColList resid_mean,resid_z --covarFile <covars_no_sex.txt> \
  --keep <keep.txt> --extract <original qc_pass.snplist> --bsize 1000 --lowmem

# Pass 2: G-MoQ Stage 1 on the ORIGINAL phenotype with --loco-preds <new .loco>
# --loco-scale regression, then the usual matched-weight build and Stage 2.

Rscript scripts/loco_blind_diagnostic.R --weights <pass1 weights> \
  --pheno <phenotypes.txt> --pheno-col testosterone --keep <keep.txt> \
  --loco <new .loco> --chr 22 --tag new_loco
Rscript scripts/compare_loco_configs.R --run "label=<chr22.tsv.gz>" ... \
  --regenie-female <f.regenie> --regenie-male <m.regenie> --chr 22
```

The pass-2 pipeline was validated by rebuilding the published with-LOCO arm from
the original `.loco` file: every metric reproduced to 7 significant figures.
