# Estimated density and quantile nuisance in RIF slopes

This note records the inference boundary for the shipped RIF construction. It
is a design note for a future validation or implementation pass; it does not
change the OLS point projection, pooled covariance, or block jackknife.

The relevant production path is `stage1/fungwas_stage1/core.py`:
`compute_rif()` builds the sample RIF, `residualize()` removes the covariate
basis, `compute_block_scores()` forms genotype-weighted block scores, and
`compute_covariance_from_blocks()` estimates the covariance while holding the
RIF fixed. `pooled_cov.compute_sigma0()` likewise consumes the already formed
residualised RIF. The density prototype that motivated this note is
`tests/validity_pilot/density_probe.py` in the validity worktree; its A/B/C
arms separate KDE-density, empirical-quantile, and population-quantile
effects.

## What the fitted slope actually uses

Let `x_i` be one genotype column after projection off the intercept and
covariates, and let

\[
 d_n = n^{-1}\sum_i x_i^2,
 \qquad
 s_n(q) = n^{-1}\sum_i x_i\{\tau-1(Y_i\leq q)\}.
\]

Because the intercept is in the projection, `sum_i x_i = 0` up to numerical
roundoff. Therefore the shipped per-tau slope can be written as

\[
 \hat\beta_\tau
 = \frac{s_n(\hat q_\tau)}{d_n\,\hat f_\tau},
\]

where `compute_rif()` uses the type-8 sample quantile and the Silverman KDE at
that quantile. This identity is useful because it exposes the nuisance terms
that are hidden when the RIF is treated as a fixed response.

Conditioning on the realised genotype and covariate rows, a first-order
expansion around the population values `q_tau` and `f_tau` is

\[
 \hat\beta_\tau-\beta_\tau \approx
 \frac{s_n(q_\tau)-s(q_\tau)}{d f}
 + \frac{\dot s(q_\tau)}{d f}(\hat q_\tau-q_\tau)
 - \frac{s(q_\tau)}{d f^2}(\hat f_\tau-f_\tau)
 + \text{higher-order and smoothing bias},
\]

with

\[
 \dot s(q)=-f(q)\,E[x\mid Y=q],
 \qquad s(q)=d f(q)\,\beta_\tau,
 \qquad d=E[x^2].
\]

The last term is consequently `-beta_tau (fhat_tau - f_tau) / f_tau`. The
quantile estimator has the
usual influence term

\[
 \hat q_\tau-q_\tau
 = n^{-1}\sum_i\{\tau-1(Y_i\leq q_\tau)\}/f_\tau
 + \text{remainder},
\]

and evaluating a density estimate at `qhat_tau` adds
`f'_tau (qhat_tau - q_tau)` to its expansion. Thus estimating the density
and quantile changes both the point estimate and its sampling covariance.

Under a conditional null for which `E[x | Y=q_tau] = 0` and
`beta_tau = 0`, the two displayed nuisance terms vanish at first order. This
explains why the density often cancels from a single-quantile null Wald
statistic and why null diagnostics can look sound even when nonnull inference
is incomplete. Under a nonzero effect, the density term is present, and the
quantile boundary term can be present as well.

For a second-order KDE with bandwidth `h`, the leading density error has
deterministic bias `O(h^2)` and stochastic size `O_p((nh)^(-1/2))`. With
the shipped `h proportional to n^(-1/5)`, both are of order `n^(-2/5)`, slower
than a root-n slope term. A valid fourth-order bias correction can reduce the
smooth deterministic term to `O(h^4)`, but it does not remove the KDE's
stochastic error, the quantile-estimation term, or density-evaluation-at-
\(\hat q\) uncertainty. For a nonzero slope, these terms can therefore remain
material after bias correction.

This is the same three-part boundary described for RIF-OLS by Firpo, Fortin,
and Lemieux: density uncertainty, population-quantile uncertainty, and the
average-derivative term all enter the asymptotic expansion; kernel density's
slower rate can dominate under the alternative, whereas testing a zero UQPE
can avoid estimating the density. See [Firpo, Fortin, and Lemieux (2009),
pp. 961–962](https://economics.ubc.ca/wp-content/uploads/sites/38/2013/05/pdf_paper_thomas-lemieux-unconditional-quantile-regressions.pdf).

## What the existing covariance estimators cover

The current block jackknife and the pooled covariance estimate the variability
of the score conditional on the realised RIF:

\[
 \psi_{i\tau}=\widehat{RIF}_{i\tau}
       -Q_iQ'\widehat{RIF}_{\cdot\tau},
 \qquad
 \widehat\beta_\tau=d_n^{-1}n^{-1}\sum_i x_i\psi_{i\tau}.
\]

`compute_covariance_from_blocks()` deletes blocks from the sums
(\sum x_i^2) and (\sum x_i\psi_i), but does not recompute `compute_rif()`
inside each deletion. `compute_sigma0()` forms the Gram matrix of this same
fixed (psi). Consequently:

* the existing JK captures genotype-score covariance and block dependence
  conditional on ((\hat q,\hat f));
* the HC0 diagnostic in `stage2/tests/robust_probe/` adds genotype-weighted
  score heteroscedasticity but also holds the RIF nuisance estimates fixed;
* neither estimator propagates the sampling variation of (hat q), (hat f),
  or the KDE bandwidth choice into a nonnull slope covariance; and
* replacing the KDE by a less biased density changes the fitted RIF and can
  reduce point bias, but does not make the fixed-RIF covariance a full
  nuisance-aware covariance.

A delete-one-block procedure that recomputes the quantile, bandwidth, density,
RIF, residualisation, and scores would target the full plug-in estimator more
closely. It is a different estimator from the shipped JK and is much more
expensive: a naive implementation repeats the (N\times T) RIF/KDE work for
each block. A fast leave-one-out KDE update would still need a validation of
quantile ties, bandwidth changes, and the resulting covariance rank.

## Minimal evidence before any correction

The next evidence pass should keep the OLS projection and estimand fixed and
measure nuisance effects directly:

1. Use the paired density prototype at each relevant sample size and LOCO or
   covariate stratum. Compare A (empirical (q), shipped KDE), B (empirical
   (q), high-accuracy or known density at (hat q)), and C (population
   (q,f) in simulation). Record (\hat f/f_{ref}), (\hat q-q_{ref}), the
   projected point shift, and empirical SD against the fixed-RIF JK/pooled SE.
   A–B isolates density substitution; B–C isolates quantile and finite-sample
   nuisance. The existing `density_bias_n10k_pair.csv` already provides this
   paired structure for the Normal variance-linear cells.
2. On a small set of causal and null variants, run a nuisance-aware bootstrap
   or delete-one-block pilot that rebuilds (q,f,) and the RIF inside every
   resample. Compare its covariance with the fixed-RIF JK and genotype-weighted
   HC0. Use dataset-level Monte Carlo units and report projection and
   scientific targets separately.
3. If the only proposed production change is fourth-order KDE bias correction,
   repeat step 1 after correction and label the result as a point-bias and
   Q-shape diagnostic. Do not claim that it supplies the missing nonnull
   covariance unless step 2 supports that claim.

This gives a bounded decision rule: correct the density bias only when its
paired point shift is material for the stated estimand; use a nuisance-aware
bootstrap or an explicitly derived augmented influence function when nonnull
SEs are needed. There is no justified scalar inflation that recovers the
missing density and quantile terms.

## Alternatives and their estimands

Estimating (f) on exactly the matched phenotype rows and within the same LOCO
stratum is necessary for compatibility. It prevents a sample or residual
distribution mismatch, but it does not remove the sampling uncertainty of the
matched density estimate. A held-out or cross-fitted density can reduce
same-sample coupling, at the cost of a new finite-sample estimator and a
variance-combination problem; it needs the paired check above.

Dropping (1/\hat f_\tau) and regressing the indicator term produces a
probability-scale slope (f_\tau\beta_\tau), not the unconditional quantile
partial effect represented by the RIF. Transforming the phenotype to an
empirical CDF or rank-normal score similarly changes the target to a rank or
probability effect and changes the appropriate weight matrix. Such a route may
be useful as a separately validated analysis, but it is not a drop-in repair
for the current `beta_mu`/`beta_sigma2` OLS projection.
