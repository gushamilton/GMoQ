# The REGENIE LOCO workflow is hazardous for the phenotypes G-MoQ targets

**Status: measured on testosterone chr22, 2026-09-16. The mechanism is general;
the size of the effect is phenotype-specific.**

## The finding

REGENIE step 1 fits a whole-genome ridge to predict the phenotype. G-MoQ then
subtracts `scale x LOCO(chr)` before building RIFs. That is the standard
workflow and it is fine for an ordinary mean GWAS.

It is not fine when the phenotype's variance is dominated by latent structure
that the step-1 model cannot see. The ridge will learn **that structure**,
because it is the largest signal available to it, and the subtraction then
erodes exactly what the distributional model exists to estimate.

This is not a corner case for G-MoQ. It is the central use case: G-MoQ is for
phenotypes with strong distributional structure, and REGENIE step 1 is built to
capture polygenic *mean* signal. The two assumptions are in tension by
construction.

## What it did to testosterone

Step 1 was run with `--covarFile covars_no_sex_paper.txt` -- age and 10 PCs, no
sex. Sex was deliberately excluded so that G-MoQ could model the sex mixture
itself. The consequence was that the ridge learned sex instead:

| quantity | value |
|---|---|
| `cor(y, loco)` | 0.159 |
| `cor(sex, loco)` | 0.160 |
| within-sex R^2 of the predictor | 0.00027 (F), 0.0012 (M) |
| `var(y - 1.09*loco) / var(y)` **within sex** | 1.203 (F), 1.264 (M) |
| `var(y - 1.09*loco) / var(y)` overall | 0.976 |

It explains essentially nothing within a sex, and adds 20-26% to within-sex
variance, while appearing to help on the overall variance. Sex separation for
this phenotype is `S = 5.89`, so sex is most of the variance and the ridge had
every incentive to spend its capacity there.

## Measured cost

> **Comparator caveat (added 2026-09-16).** The table below is measured against
> a REGENIE that was run `--ignore-pred`, i.e. with **no polygenic adjustment of
> its own**. That is a handicapped reference: anyone holding the sex labels
> would fit a within-stratum step 1. The numbers here remain correct for the
> question this section asks (*what did G-MoQ's broken LOCO cost?*), because
> both arms share the same reference. They are **not** the right numbers for
> *how does G-MoQ compare to a sex-stratified GWAS at its best*. For that, and
> for the full register of every comparison in this project, see
> [COMPARATOR_LEDGER.md](COMPARATOR_LEDGER.md).

chr22 testosterone, pooled + OLS, median SE relative to sex-stratified REGENIE
(which was itself run `--ignore-pred`, so no-LOCO is the like-for-like
comparison):

| component | with LOCO | **no LOCO** |
|---|---|---|
| female (`beta_left`) SE ratio | 1.1209 | **1.0236** |
| male (`beta_right`) SE ratio | 1.0677 | **0.9311** |
| female cor with REGENIE, all 17,550 | 0.8981 | **0.9699** |
| male cor with REGENIE, all | 0.8165 | **0.9188** |
| female cor, \|z\|>=4 | 0.9958 | **0.9984** |
| male cor, \|z\|>=4 | 0.9943 | **0.9967** |

Removing it improved precision *and* accuracy, which rules out a variance
shrinkage artefact. The male ratio below 1.0 is expected rather than anomalous:
at high separation G-MoQ uses all 399,710 samples for both components while
stratified analysis discards the other sex, and joint estimation of the shared
structure recovers more than membership uncertainty costs.

## Genomic inflation: measured, and it revises the recommendation

Two lambda measures on real chr22 are each confounded -- over all variants by
genuine signal, over REGENIE-null variants (|z|<1) by a selection artefact,
since that selects variants where G-MoQ also sees nothing. Neither level is
interpretable. A permuted null settles it: the testosterone column was shuffled
across 421,405 people, preserving the phenotype distribution, the mixture
structure and the RIF construction exactly while destroying genotype
association.

| | lambda gamma | lambda beta_left | lambda beta_right |
|---|---|---|---|
| **permuted, true null** | **0.999** | **0.904** | **0.993** |
| real, no LOCO | 1.265 | 1.420 | 1.487 |
| real, with LOCO | 1.267 | 1.353 | 1.360 |

| | p<0.05 | p<1e-3 | min p over 17,550 |
|---|---|---|---|
| permuted null | 0.042 | 0.0004 | 2.9e-5 |
| real, no LOCO | 0.105 | 0.0148 | 1.7e-13 |
| real, with LOCO | 0.101 | 0.0116 | 2.1e-11 |

**The estimator is calibrated and slightly conservative** -- 4.2% against 5%
nominal, 0.04% against 0.1%. So the real-data lambda of 1.35-1.49 is genuine
chr22 signal, not miscalibration, and it also disposes of the concern that the
no-LOCO standard-error improvement might be an artefact of SEs that are simply
too small.

> **SUPERSEDED, 2026-09-16.** The paragraph below concluded that the 5-9% lambda
> rise on removing LOCO was real lost relatedness/structure control. **That is
> wrong.** Building and running two-pass blind LOCO showed the lambda gap was a
> variance-inflation artefact of subtracting a membership proxy: the old LOCO left
> `gamma`'s SE *and* `gamma`'s lambda untouched while inflating `beta_left`/
> `beta_right` SEs by 9.5%/14.7%, and only those two lambdas moved. Genuine
> confounding inflates every statistic from the same genotypes; this deflated only
> the two it had inflated. See [TWO_PASS_LOCO.md](TWO_PASS_LOCO.md#what-the-lambda-actually-was).
> The original paragraph is kept below for the record.

~~**RETRACTED 2026-09-16 (later the same day).** This section previously argued
that the 5-9% lambda rise on removing LOCO was real lost relatedness control.
That is wrong, and lambda could not have shown it either way. Masking the signal
-- dropping 1 Mb around every REGENIE |z|>=5 hit, 5,114 variants including SHBG
-- leaves lambda unchanged or higher in every arm, so at N=400,000 this lambda is
polygenicity, not confounding. Two further lines agree: lambda moved only for the
two parameters whose SEs the subtraction inflated, and gamma, whose SE ratio
between arms was 1.0037, had lambda 1.2650 against 1.2671. The rise was the
arithmetic consequence of inflated standard errors. The machinery is calibrated and the signal is
held constant, so the change is attributable to LOCO. The offset was failing at
polygenic prediction (within-sex R^2 = 0.0003) while genuinely doing structure
control -- two different jobs, one of them working.~~

**Revised recommendation: two-pass blind LOCO (option 1) is the best available
option and plain no-LOCO is an acceptable fallback.** Two-pass recovers the whole
precision and accuracy cost -- SEs 1-3% better than no-LOCO, 10-15% better than
with-LOCO -- while a blind diagnostic clears its predictor of membership proxying.
It does not lower lambda, because there was no LOCO-removable inflation to lower.

Residual caveat: conservatism of ~2.5x at p<1e-3 costs power, and its behaviour
at 5e-8 cannot be measured from 17,550 variants. It may partly explain the
recall gaps against sex-stratified REGENIE (0.41-0.47 on the female component).

## What was ruled out

- **chrX leakage is NOT the cause.** Step 1 used `chr1-22_merged`; the bim's
  highest chromosome is 22. The 23rd row of the `.loco` file is REGENIE emitting
  a prediction for a chromosome with no variants to leave out. An earlier draft
  of this analysis asserted chrX leakage; it is wrong.
- **Adding sex to step 1 is not an acceptable fix.** It works mechanically, but
  sex *is* the mixture membership, so using it contaminates the claim that G-MoQ
  recovers component effects without knowing membership -- and in a genuine
  latent-mixture application (disease subtypes, antibody classes) no such
  column exists.

## The blind diagnostic

Membership is not needed to detect this. Stage 2 already fits the mixture, so
posterior membership weights `w(y)` are available per sample. Then:

1. **Regress the LOCO predictor on the posterior weights.** High R^2 means the
   predictor is a membership proxy and subtracting it will erode signal.
2. **Compute the predictor's R^2 within posterior-weighted components.** Near
   zero means it is contributing no polygenic adjustment to pay for the damage.

Both use only what the method infers. Neither needs a label. This check would
have caught the problem immediately and should run before any LOCO adjustment
is trusted on a structured phenotype.

**Implemented as `scripts/loco_blind_diagnostic.R`.** On the old testosterone
LOCO it returns a membership-proxy ratio of 17.4 and flags it; on the two-pass
LOCO it returns 0.006-0.012 and clears it. Measured values in
[TWO_PASS_LOCO.md](TWO_PASS_LOCO.md#the-blind-diagnostic-on-old-and-new-loco).

## Options, in order of preference

1. **Two-pass blind LOCO** -- built, tested and recommended. Run G-MoQ with no
   LOCO, take the estimated mixture posteriors, residualise the phenotype on the
   *estimated* structure, refit step 1 on that residual, rerun. Uses only
   inferred quantities, so it stays legitimate. Measured: SEs 1-3% better than
   no-LOCO, effect correlations within 0.004 of it, and a LOCO predictor that
   explains 2.2-5.7% within a component instead of 0.03-0.15%. Genomic inflation
   is unchanged from no-LOCO. Full results and the recommended transform in
   [TWO_PASS_LOCO.md](TWO_PASS_LOCO.md).
2. **No LOCO.** Acceptable fallback: it costs 1-3% in SE against two-pass and
   nothing in lambda. The "+5% female, +9% male" inflation penalty recorded in
   an earlier draft of this document was the with-LOCO arm's SE inflation, not a
   no-LOCO penalty.
3. **Step 1 with membership as a covariate.** Correct only when membership is
   genuinely known AND the analysis does not claim to work without it. Not
   available for real latent-mixture applications.

## Other traits at risk

Any trait whose latent structure dominates its variance and whose step-1 run
omitted that structure. In this project that includes the `mixture`-family
antibody traits `jcv_log_all` and `mcv_log_all`, which are bimodal
seropositive/seronegative and went through the same workflow. Their LOCO
predictors have not been checked. Run the blind diagnostic above before using
those results.
