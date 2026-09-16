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

Genomic inflation with and without LOCO has NOT yet been checked. Until it has,
dropping LOCO trades a known bias for an unquantified structure/relatedness
risk.

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

## Options, in order of preference

1. **No LOCO**, once genomic inflation is confirmed acceptable. Simplest, and
   matches the `--ignore-pred` comparator.
2. **Two-pass blind LOCO.** Run G-MoQ with no LOCO, take the estimated mixture
   posteriors, residualise the phenotype on the *estimated* structure, refit
   step 1 on that residual, rerun. Uses only inferred quantities, so it stays
   legitimate, and recovers genuine polygenic and relatedness control.
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
