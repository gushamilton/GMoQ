# LOCO belongs in the covariate basis, not in the phenotype

Status: measured on testosterone chr22, 2026-09-16. Implemented as
`--loco-mode covariate`; the default stays `subtract`, and no published result
changes.

## The semantic argument

A RIF/UQR estimand is defined against the **unconditional** distribution of `Y`:

```
RIF_i(tau) = q_tau + (tau - 1{Y_i <= q_tau}) / f(q_tau)
```

`Y` enters only through where each sample sits relative to the *population's*
tau-th quantile. Change `Y` and you move `q_tau`, you move `f(q_tau)`, and you
flip indicators — you have redefined the target, not removed a nuisance. This is
exactly where a mean GWAS's intuition fails to carry over: there, subtracting any
`G`-independent term leaves the slope alone and only cuts residual variance, so
the adjustment really is nuisance.

Stage 1 already treated its two adjustments asymmetrically, and its own comment in
`scripts/stage2/prepare_matched_weights.R` said so:

- **covariates** are residualised out of the RIF *columns*, after the RIF is
  built. The unconditional distribution is untouched. Correct for UQR.
- **LOCO** was subtracted from `Y` *itself*, before the RIF was built. The
  distribution moves. The estimand moves with it.

`--loco-mode covariate` makes the second look like the first: per chromosome, the
LOCO prediction `L(chr)` becomes an extra column of the basis the RIF columns are
residualised on, and `Y` is left alone. Because Stage 1's `Q` also projects the
genotype in the scoring kernel, `L` is partialled out of both sides — this is
Frisch-Waugh, i.e. the coefficient on `G` in `RIF(tau) ~ G + X + L(chr)`.

Two consequences fall out of the mechanism:

- `--loco-scale` is inert under `covariate`. A projection is invariant to
  rescaling the column, so the entire `docs/LOCO_SCALE.md` hazard — an unscaled
  REGENIE prediction, or a scale factor that over- or under-corrects — cannot
  arise. Stage 1 resolves no factor in this mode.
- the basis gains one column (`covariate_basis_rank` 12 -> 13) and the
  `Sigma_0` fingerprint changes, so a covariate-mode sidecar cannot be silently
  consumed by a subtract-mode scan.

## Why testosterone is the case that exposes it

For this trait the REGENIE step-1 ridge learned **sex**, not polygenic signal:
`cor(sex, loco) = 0.160` against `cor(y, loco) = 0.159`, and the within-sex
`R^2` is 0.0003. Subtracting `1.088 * loco` therefore removes part of the sex
separation and *adds* 20-26% to within-component variance. Separation is what
governs how much information exists about component-specific effects, so the
subtraction attacks the estimand at its most informative point.

## Measurement

testosterone chr22, skew family, pooled covariance + OLS, 17,550 HapMap3
variants, N = 399,710. Compared against sex-stratified REGENIE
(`beta_left`<->female, `beta_right`<->male, `gamma`<->sex with `gamma` negated,
since REGENIE codes Male = 1); alleles harmonised against `ALLELE1`/`ALLELE0`
with the sign flipped when swapped.

Configurations: **(a)** subtract-from-`Y` (production);
**(a')** the same, rebuilt from this branch — Stage 1 `stage1.tsv.gz` and
`cov.gz` are byte-identical to the paper reference, and every summary below
agrees with (a) to four decimals; **(b)** no LOCO at all; **(c)** two-pass blind
LOCO; **(d)** `--loco-mode covariate`.

`SE` is the median ratio to REGENIE's SE; `cor` is over all variants; `cor_hiz`
restricts to REGENIE `|z| >= 4`; `lambda` is genomic inflation.

| | `beta_left` SE / cor / cor_hiz / lambda | `beta_right` SE / cor / cor_hiz / lambda | `gamma` SE / cor / cor_hiz / lambda |
|---|---|---|---|
| (a) subtract from Y | 1.1209 / 0.8981 / 0.9958 / 1.3527 | 1.0677 / 0.8165 / 0.9943 / 1.3604 | 1.0754 / 0.9036 / 0.9999 / 1.2671 |
| (b) no LOCO | 1.0236 / 0.9699 / 0.9984 / 1.4195 | 0.9311 / 0.9188 / 0.9967 / 1.4869 | 1.0714 / 0.9085 / 0.9999 / 1.2650 |
| (c) two-pass blind | 1.0113 / 0.9664 / 0.9982 / 1.4221 | 0.9134 / 0.9140 / 0.9966 / 1.4778 | 1.0712 / 0.9062 / 0.9999 / 1.2560 |
| **(d) LOCO as covariate** | **1.0235 / 0.9699 / 0.9984 / 1.4200** | **0.9309 / 0.9187 / 0.9968 / 1.4858** | **1.0574 / 0.8990 / 0.9999 / 1.2239** |

A permuted true null gives lambda 0.9045 / 0.9928 / 0.9990, so the estimator
itself is calibrated.

**(d) recovers no-LOCO precision and accuracy exactly.** Against (b), the point
estimates correlate 0.99993 (`beta_left`) and 0.99978 (`beta_right`), and the
median SE ratio is 1.0001 and 1.0003. The 9.5% and 14.7% SE inflation that
subtraction imposed is gone, and agreement with REGENIE rises from 0.898/0.817
to 0.970/0.919.

## Was the subtract-LOCO lambda advantage ever structure control?

No. Three lines of evidence, all pointing the same way.

**It moved only where it inflated SEs.** From (a) to (b) the median SE ratio
changes by 1.0950 (`beta_left`) and 1.1467 (`beta_right`) but 1.0037 (`gamma`).
lambda tracks that exactly: 1.42 -> 1.35 and 1.49 -> 1.36 for the two inflated
parameters, 1.2650 -> 1.2671 for the one that was not. Confounding from
relatedness does not politely confine itself to two of three parameters; SE
arithmetic does.

**Handing (b) the subtraction's SEs over-explains the gap.** Give each run's own
point estimates run (a')'s standard errors and recompute lambda: (b) lands at
1.1839 and 1.1308, *below* (a')'s own 1.3527 and 1.3604. The SE inflation
accounts for the whole lambda difference and then some — per unit of SE, the
subtracted estimates carry *more* median chi-square, not less.

**Masking the signal does not lower lambda.** Dropping every variant within 1 Mb
of a REGENIE `|z| >= 5` hit (5,114 of 17,550, including SHBG) leaves lambda
unchanged or slightly higher in every configuration (e.g. (d) `beta_left`
1.4200 -> 1.4724). At N = 400,000 on a polygenic trait, a chromosome-wide lambda
near 1.4 is what polygenicity looks like. **lambda cannot separate confounding
from signal at this sample size, in any of these configurations** — which is the
main reason it should not have been read as evidence for the subtraction.

## Does the covariate form actually control relatedness, or only appear to?

Partly, in a way that is real but strictly bounded by what `L` contains.

**The mechanism is sound.** Structure confounding acts on the genotype-phenotype
association: it is `cov(G, RIF_resid) != 0` driven by shared ancestry. Because
`Q` residualises *both* the RIF columns and `G`, the estimate is
`cov(RIF ⊥ L, G ⊥ L) / var(G ⊥ L)`. Any confounding that lies in the linear span
of `L` is removed from both sides. That is genuine adjustment, and it is strictly
more general than the subtraction, which is the same adjustment with the
coefficient pinned to `c` and `G` left unprojected. `L` is also a legitimate
covariate here rather than a collider or mediator: LOCO excludes the test
chromosome, so nothing on chr22 causes `L`.

**And it does something the subtraction never did.** `gamma` is the only place
in this analysis where structure control was ever visible, and (d) is the only
configuration that moves it: lambda 1.2650 -> 1.2239, with the estimates
shrinking (median `|beta|` ratio to (b) 0.985, regression slope 0.971) *while the
SEs also shrink* (median ratio 0.987). Estimates pulled toward the null with
precision improving is the signature of deconfounding; it is the opposite of the
SE-inflation artefact diagnosed above, which it cannot be confused with. At the
genuine sex loci (REGENIE `|z| >= 4`) the high-z slope against REGENIE moves
1.124 -> 1.090, i.e. toward 1, so the shrinkage is not indiscriminate
attenuation of real signal.

**But its reach is only as long as `L`.** Residualising on `L` removes the part
of each RIF column that is *linear in `L`*. Structure that acts on the shape of
the conditional distribution rather than on its mean — heteroskedasticity by
ancestry, tail behaviour, differential component membership beyond the linear
projection — is untouched. The subtraction does not reach it either, and does so
while corrupting the estimand, so this is not an argument for the old default.
It is an argument against claiming more than was measured.

For testosterone specifically the bound bites hard: because the ridge learned
essentially only sex, `L` carries almost no within-component information
(within-sex `R^2` = 0.0003), and the covariate correction is consequently a
near-no-op for `beta_left` and `beta_right` (r = 0.9999 against no LOCO). The
correction lands where the predictor actually has content — the between-component
axis, `gamma` — and nowhere else. The right reading is that `--loco-mode
covariate` routes whatever the predictor knows to the parameter it is about,
instead of pushing it through the phenotype where it damages every parameter. A
trait whose step-1 ridge learned real polygenic signal should see more from it;
this one had little to give.

## Limits of this evidence

- One chromosome, one trait, one family. No replication.
- All three REGENIE references were run with `--ignore-pred` — they carry no
  LOCO themselves. So the SE and correlation columns measure efficiency and
  agreement against a *matching unadjusted* estimator. They are the right
  yardstick for precision and for estimand distortion; they cannot certify
  structure control, and (d)'s slightly lower `gamma` correlation (0.8990 vs
  0.9085) is exactly what moving away from an unadjusted reference looks like,
  whether the move is deconfounding or over-adjustment. The two are not separated
  here.
- The `gamma` effect is ~3% on the estimates. It is consistent across three
  statistics but it is small.
- The pooled-covariance sidecar still lists only `--covar` columns in
  `covariates`; the LOCO column shows up in `covariate_basis_rank` and
  `loco_mode`, not in that list.

## Reproducing

Stage 1: add `--loco-mode covariate` to the usual `--loco-preds` invocation.
Stage 2: with `Y` unadjusted, the matched baseline distribution for `W` is raw
`Y`, so the weights are built exactly as for a no-LOCO scan (the Stage 2 weights
builder is invoked without `--loco`). The weights are byte-for-byte the no-LOCO
weights: baseline sd 1.28347632785109 in both.
