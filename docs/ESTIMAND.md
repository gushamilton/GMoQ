# What a G-MoQ effect estimates

## Stage 1: unconditional quantile effects

The Stage 1 coefficient at quantile `tau` is an *unconditional* quantile effect: how the
`tau`-th quantile of the phenotype in the population would move under a small shift in the
genotype distribution, with the covariate basis held fixed. It is not the conditional quantile
effect of classical quantile regression, which describes quantiles of the phenotype within strata
of the covariates.

The estimand is defined on the distribution of the phenotype you analysed. Two consequences:

- Anything done to the phenotype before Stage 1 — a transformation, subtracting a LOCO
  prediction, subtracting covariate effects — changes the estimand rather than removing a
  nuisance. Covariates belong in `--covar`; Stage 1 projects them out of the RIF and the
  genotype, leaving the phenotype's quantiles untouched.
- Record the analysed phenotype (its transformation and sample) with the results, and build the
  Stage 2 weights from its quantiles.

Effects are per copy of the first BGEN allele; see [ALLELE_CONVENTION.md](ALLELE_CONVENTION.md).

## Stage 2: projection onto a model

Stage 2 reports the change in model parameters that best explains a variant's quantile profile,
`(W'W)^(-1) W' b`. When the model fits (Q is unremarkable) this is, to first order, the change in
the parameters per copy of the effect allele. When it does not, it is the least-squares projection
of the quantile profile onto the model's directions, which is still well defined but should not be
read as a parameter change.

For the Normal mean/variance model (`make_weights_vqtl()`):

- `beta_mu` is on the phenotype scale, per allele;
- `beta_sigma2` is on the squared phenotype scale (the variance), per allele.

These are local, linear interpretations. For example, if the standard deviation is generated as
`SD(Y | g) = exp(b * g)`, the reported variance coefficient is not `b`: the derivative of the
variance at `g = 0` is `2b` times the baseline variance, and the one-allele variance ratio is
`exp(2b)`. Strong effects can depart from the linear approximation.

## Uncertainty

The primary analysis uses OLS with the pooled covariance. With the same quantile effects and the
same `W`, the point estimate is identical under the pooled and the jackknife covariance; only the
uncertainty differs. In a paired simulation of a Normal phenotype with a variance effect at minor
allele frequency 0.05, nominal 95% intervals for the variance parameter covered 92.9% of the time
with the pooled covariance and 94.5% with the jackknife. Changing the covariance does not repair
confounding or change what is estimated. See [METHOD.md](METHOD.md) for when to use which, and for
the other known limitations.

Null calibration, agreement with another program, or a non-significant Q do not by themselves
show that the parameter used to generate a phenotype is recovered without bias.
