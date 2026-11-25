# How G-MoQ works

G-MoQ (Genomic Method of Quantiles) runs in two stages. Stage 1 estimates, for every variant,
its effect on a grid of quantiles of the phenotype and the covariance of those effects. Stage 2
maps the quantile effects onto the parameters of a distributional model that you choose. This
page describes both stages, the choices behind the defaults, and the known limitations.

What the reported effects mean is covered separately in [ESTIMAND.md](ESTIMAND.md).

## Stage 1: quantile effects by RIF regression

For a phenotype `y` measured on `N` people and a quantile level `tau`, the recentered influence
function (RIF) of the `tau`-th quantile is

    RIF_i(tau) = q_tau + (tau - 1{y_i <= q_tau}) / f(q_tau)

where `q_tau` is the sample quantile (type 8, median-unbiased) and `f(q_tau)` is the density of
`y` at that quantile. Stage 1 estimates `f` with a Gaussian kernel and Silverman's bandwidth,
`1.06 * sd(y) * N^(-1/5)`. Regressing the RIF on genotype gives the effect of the genotype on the
unconditional `tau`-th quantile (Firpo, Fortin and Lemieux, 2009).

For each variant Stage 1:

1. builds the RIF at every quantile on the grid (by default 45 levels from 0.02 to 0.98);
2. projects the covariates (plus an intercept) out of every RIF column and out of the genotype;
3. regresses each residualised RIF column on the residualised genotype, giving one slope
   `beta_tau` per quantile.

Step 2 means covariates are handled by partialling out (Frisch–Waugh–Lovell): `beta_tau` is the
genotype coefficient in `RIF(tau) ~ G + covariates`. The phenotype itself is never adjusted, so
the quantiles stay those of the phenotype you supplied.

The genotype is the expected count of the first BGEN allele (see
[ALLELE_CONVENTION.md](ALLELE_CONVENTION.md)). Missing genotype calls are replaced by the
variant's mean dosage. Variants with a minor allele count below 5, with more than two alleles,
or whose rsID appears on several BGEN records are skipped.

Stage 1 also writes `d_resid = sum_i g_resid,i^2`, the residual genotype sum of squares. It is
the denominator of every slope and scales the covariance below.

## The covariance of the quantile effects

Stage 2 needs the covariance of the `T` slopes of each variant across quantiles. Stage 1
provides two estimates.

**Pooled covariance (the default).** Write `Psi` for the `N x T` matrix of residualised RIFs.
Stage 1 writes one matrix per phenotype (per chromosome when a LOCO prediction is used),

    Sigma_0 = Psi' Psi / N,

and the covariance of a variant's slopes is `Sigma_0 / d_resid`. This is the heteroscedasticity-
robust (HC0) sandwich under working independence: it assumes the genotype does not change the
spread of the residual RIF. In simulation it reproduces the exact per-individual sandwich with a
median ratio of 1.00. It is positive semi-definite by construction, needs no blocks or random
seed, and costs a few kilobytes per chromosome.

**Block jackknife covariance.** Stage 1 also assigns samples at random to `--blocks` blocks
(default 1,000) and computes a delete-one-block jackknife covariance for every variant. It does
not assume working independence, but it is noisier, and its noise matters in the tail (below).

The two sources give identical OLS estimates; only standard errors, p-values and the Q test
differ. Use the pooled covariance for the primary analysis. Check strong signals with
`cov_source = "jackknife"` when the variant also affects the variance (a dispersion effect), and
particularly when it is rare: there the working-independence assumption fails and the pooled
standard error can be too small. Component-specific mixture parameters behave similarly.

Stage 1 records which phenotype, covariates, sample, quantile grid and LOCO predictor each
covariance belongs to. Stage 2 refuses a pooled covariance whose record does not match the scan
it is applied to, because a mismatch would leave the estimates unchanged and the standard errors
plausible but wrong.

## Stage 2: from quantile effects to model parameters

A model for the phenotype's distribution with parameters `theta` implies a quantile function
`q_tau(theta)`. Its derivatives form the weight matrix

    W[t, k] = d q_tau / d theta_k = -(dF/dtheta_k)(q_tau) / f(q_tau),

one row per quantile and one column per parameter. For a small genetic effect the quantile
effects are approximately `b = W delta`, where `delta` is the change in the parameters per copy of
the effect allele. Stage 2 estimates `delta` by ordinary least squares,

    delta_hat = A b,    A = (W'W)^(-1) W',

with covariance `A Sigma A'`, where `Sigma` is the Stage 1 covariance of `b`. Builders for common
models are in the R package: `make_weights_vqtl()` (Normal mean and variance),
`make_weights_normal_mixture()` (two-component Normal mixture, optionally with a membership
parameter), `make_weights_normal_mixture_vqtl()` and `make_weights_generic()` (any family for
which you can write the CDF, density and parameter gradients). Build `W` from the quantiles of
the same phenotype values Stage 1 analysed.

For mixtures with a membership parameter, the density that divides each column matters; see
[MIXTURE_WEIGHT_DENSITY.md](MIXTURE_WEIGHT_DENSITY.md).

### OLS or GLS

`tau_cov_estimator = "gls"` weights the fit by the inverse covariance instead. When the model is
correctly specified GLS is more efficient (in simulation its estimates were roughly 20% less
variable, with the same coverage). When the model does not fit the quantile profile, GLS is
pulled further from the truth than OLS and its standard errors look misleadingly small; a sandwich
standard error does not repair that, because the problem is the point estimate. OLS is the
default. Use GLS only for a model whose fit you have checked with Q.

### P-values and degrees of freedom

Parameter p-values are two-sided t tests of `delta_hat / SE`. The degrees of freedom reflect how
precisely the covariance was estimated:

- **Pooled covariance.** `Sigma_0` is estimated from all `N` people, so the reference is
  effectively normal.
- **Jackknife covariance.** With `B` blocks the effective degrees of freedom are well below
  `B - 1` and vary by variant and by parameter. Stage 2 uses
  `1/nu = 1/(B - 1) + c / d_resid`, with `c = 0.88` for location-type parameters and `c = 2.62`
  for scale-type parameters, constants fitted to simulations. Stage 1 records `B` in the
  `.cov.meta.json` sidecar. A normal reference with a noisy jackknife inflates the false-positive
  rate in the far tail, and more so with fewer blocks.

If Stage 2 cannot establish the degrees of freedom it reports no p-value rather than falling
back to a normal reference (`unknown_nu` controls this).

Do not test the per-quantile slopes in the Stage 1 file (`beta_tau / se_tau`) against a normal
distribution or against `B - 1` degrees of freedom: at extreme quantiles their effective degrees
of freedom can be far smaller.

### The lack-of-fit statistic Q

`Q` measures how much of a variant's quantile profile the model leaves unexplained,
`Q = sum_t r_t^2 / se_t^2` with `r = b - W delta_hat`. Under OLS its reference distribution is a
weighted sum of chi-squares; `Q_pval` uses a saddlepoint approximation to it. Under GLS, `Q` is
Hotelling's `T^2` and `Q_pval` uses the exact F reference. Both treat the covariance as known,
so read `Q_pval` as a diagnostic for model adequacy rather than as an exact test, and compare
models with `Q / (Q_df_eff * Q_scale)` (null mean 1) rather than with the raw `Q / (T - K)`.

### Functions of the parameters

For an effect on a derived quantity (for example the overall mean of a mixture, or a tail
probability), define the function on the baseline parameters and use `functional_grad()` and
`functional_estimate()`: they apply the delta method to the fitted parameters and their
covariance. Do not average the quantile effects over the grid to estimate such a quantity.

## LOCO predictions

`--loco-preds` reads REGENIE step 1 predictions (a `.loco` file or a `_pred.list`). By default
(`--loco-mode covariate`) the per-chromosome prediction is added to the covariates and the
phenotype is left unchanged, so the estimand stays defined on the phenotype's own distribution.
`--loco-mode subtract` subtracts the prediction from the phenotype before the RIF is built; that
changes the quantiles, the density and every indicator, and so changes what is estimated. If you
subtract, keep `--loco-scale regression` (the default), which puts REGENIE's standardised
prediction on the phenotype's scale.

A LOCO predictor can only remove structure it represents, and it can do harm. If the phenotype is
a mixture of groups (for example a trait whose distribution differs sharply by sex), a step 1
model fitted without the grouping variable can learn group membership rather than polygenic
signal. Adjusting for it then removes the very contrast a mixture model estimates. For such
phenotypes, fit step 1 with the grouping variable as a covariate when it is known, or run without
LOCO.

## Phenotype distribution and dispersion artefacts

Two mechanisms can make a variant with a large mean effect look like a variance effect:

- **Mean-variance coupling.** On a scale where the variance grows with the mean (most raw,
  right-skewed measurements), any mean effect drags the variance with it. Transform the phenotype
  (for example log or rank-based inverse normal) or use `estimate_rav()` and
  `residualise_dispersion()` to remove the coupling genome-wide.
- **Genotype skewness.** At low allele frequency the genotype distribution is skewed, and a pure
  mean effect leaves a second-order residue in the variance parameter. `dispersion_skew_guard()`
  flags variants at risk and `correct_dispersion_skew()` corrects the Normal mean/variance design
  (check it first with `check_vqtl_design()`). `dispersion_report()` runs these together.

A strongly skewed phenotype can also inflate dispersion tests at common variants, where neither
tool applies; transform such phenotypes.

## Known limitations

- **Pooled covariance under dispersion effects.** See above: confirm strong dispersion or mixture
  component signals, especially at low frequency, with the jackknife.
- **Uncertainty from the quantiles and the density.** Both covariance estimates hold `q_tau` and
  `f(q_tau)` fixed. That is correct under the null but omits some uncertainty for non-zero
  effects, so confidence intervals for large effects can be slightly too narrow. Null calibration
  does not establish calibration away from the null.
- **Covariate adjustment is linear.** A separate linear coefficient at each quantile does not
  control nonlinear confounding of the mean or the variance. Include the covariate terms that the
  phenotype's distribution actually depends on.
- **Density estimation.** Keep `--density-method silverman_gaussian`. The alternative
  `fourth_order_gaussian` is unreliable for multimodal, rounded or count phenotypes. Heavily tied
  phenotypes (few distinct values) are outside what has been validated.
- **Relatedness and other model families.** Validation used unrelated samples and smooth,
  small-effect location, variance and mixture settings. Strong effects, relatedness and other
  families need their own checks.

## Computation

Stage 1 is multithreaded with OpenMP (`--threads`, set to the number of physical cores). With a
`.bgi` index and `--snps`, the selected BGEN records are read and decoded in parallel, and
dosages are the exact expected allele count. The quantile cross-products are formed from each
sample's quantile bin, so their cost does not grow with the number of quantiles. The fast
kernels are algebraically identical to the reference kernel but sum in a different order, so
results differ in the last bits; `--exact-kernel` reproduces the reference kernel bit for bit.
`--cov-mode pooled` skips writing the per-variant jackknife covariance, and `--cov-mode
pooled-only` also skips computing it; then `se_tau_*` are the pooled standard errors. Several
phenotypes can share one pass over the genotypes (`--jobs`).

## Reference

Firpo S, Fortin NM, Lemieux T (2009). Unconditional quantile regressions. *Econometrica*
77:953–973.
