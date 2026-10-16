# Background

Short notes on what the method estimates and on the choices behind the defaults. Each links to
the detailed document.

## What a G-MoQ effect means

Stage 1 is a recentered influence function (RIF) regression, so the coefficient at quantile `tau`
is an *unconditional* quantile effect: how the `tau`-th quantile of Y in the population would move
under a marginal shift in genotype, holding the covariate basis fixed. It is not the conditional
quantile effect of classical quantile regression. Stage 2 reports the parameters of your chosen
model that best explain the whole quantile profile.

Two consequences:
- Anything subtracted from Y before Stage 1 (a LOCO prediction, a transformation) changes the
  estimand rather than removing a nuisance.
- The estimand is defined on the analysed phenotype, so record that phenotype with the result,
  and build Stage 2's `W` from its quantiles.

Detail: [ESTIMAND.md](ESTIMAND.md).

## Standard errors: pooled or jackknife

The primary analysis is OLS with pooled covariance (`cov_source = "pooled"`, the default). Under
OLS the two covariance sources give identical estimates and differ only in the standard error.
Stage 1 writes both.

The pooled SE can be too small for rare variants with strong dispersion effects and for
component-specific mixture parameters. For those, use `cov_source = "jackknife"` as a check.
Detail: [FINDINGS_covsource_default.md](FINDINGS_covsource_default.md),
[OLS_VS_GLS.md](OLS_VS_GLS.md).

## LOCO predictions

With `--loco-preds`, the per-chromosome prediction goes into the covariate basis
(`--loco-mode covariate`, the default) rather than being subtracted from Y, because subtracting
changes the estimand. Check a LOCO predictor with `scripts/loco_blind_diagnostic.R` before
trusting it: a step-1 model fitted without the mixture variable can learn the mixture.
Detail: [LOCO_AS_COVARIATE.md](LOCO_AS_COVARIATE.md), [LOCO_SCALE.md](LOCO_SCALE.md),
[LOCO_MIXTURE_HAZARD.md](LOCO_MIXTURE_HAZARD.md).

## Mixture weights

For a mixture with a membership parameter, the default `density = "hybrid"` divides the
membership column of `W` by the Stage 1 kernel density, which is the density the Stage 1
coefficients carry. It needs `stage1_density = stage1_kde_density(y, q_tau)`.
`density = "model"` divides every column by the fitted model density instead.
Detail: [MIXTURE_WEIGHT_DENSITY.md](MIXTURE_WEIGHT_DENSITY.md).

## Stage 1 density

Leave `--density-method` at its default (`silverman_gaussian`). The alternative,
`fourth_order_gaussian`, is unreliable for multimodal, rounded and count phenotypes.

Known limits: [KNOWN_ISSUES.md](KNOWN_ISSUES.md).
