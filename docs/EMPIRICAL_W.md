# Should W be built from the empirical baseline? No.

Status: settled 2026-09-15 by exact quadrature, harness validated against the Normal
control first. **Do not adopt the empirical weight column.** Keep
`make_weights_vqtl`'s `qnorm(tau)/(2 sd)` and keep `check_vqtl_design()` gating on it.

## The question

`make_weights_vqtl(taus, q_tau, mu, sd)` validates `q_tau`'s length and then never uses
`q_tau` or `mu`, hardcoding a Normal baseline. Since Stage 1 already computes the
empirical quantiles and density, would building the dispersion column from them --
`(q_tau - mu) / (2 var)` -- fix the dispersion-skew correction on skewed phenotypes?

## Answer

No. It helps in exactly one cell and hurts nearly everywhere else.

Harness validation first: Normal, MAF 0.01, `Z_mu` 30, N = 461,320, corrected at 5e-8
gives **1.438e-09 (OLS)** and **1.537e-09 (GLS)** against uncorrected 0.9124 / 0.9161 --
reproducing the established 0.912 -> 1.5e-09. The shipped and empirical columns agree to
8.75e-14 on a Normal baseline, and the reconstructed shipped column matches the real R
`make_weights_vqtl()` to 4.44e-16 over all 45 taus.

Corrected false-dispersion rate at 5e-8, MAF 0.01, N = 461,320, `Z_mu` 30, pure location
effect with zero true dispersion, shipped `W` against empirical `W`:

| shape | OLS | GLS |
|---|---|---|
| Normal | 1.438e-09 -> 1.438e-09 | 1.537e-09 -> 1.537e-09 |
| t(6) | 1.415e-07 -> 5.850e-07 | 3.111e-06 -> 1.376e-05 |
| chi2(20) | 8.249e-08 -> 1.524e-07 | 3.753e-10 -> **2.379e-07** |
| chi2(4) | **4.511e-04 -> 5.424e-06** | 1.000 -> 1.000 |

The only gain is chi2(4) under OLS, 83x better -- and still 108x nominal, so not a fix.
The worst loss is chi2(20) under GLS, **630x worse**. Power at matched *actual* Type-I
falls in every non-Normal cell (chi2(20)/GLS 0.820 -> 0.544; t(6)/OLS 0.556 -> 0.466),
and is unchanged on Normal.

**The chi2(4)/OLS gain is a cancellation, not a repair.** Bias decomposition in
corrected-SE units at `Z_mu` 30: shipped S = +0.325, H = +1.326, terms adding to 1.650;
empirical S = -0.743, H = +1.064, terms cancelling to 0.321. The correction efficiency is
*further* from 1 for empirical (0.873 against 1.048).

Adopting the empirical column would also require disabling `check_vqtl_design()`, which
refuses elementwise anything that is not `qnorm(tau)/(2 sd)`.

## Two corrections to the existing record

**1. The prior non-Normal numbers were already empirical-W.**
`sim_skewcorr/src/exact.py` computes `W = [ones, (q0-m0)/(2 v0)]`, which *is* the
empirical column. Every non-Normal figure in `DISPERSION_SKEW_CORRECTION.md` section 5
was therefore already an empirical-W result, while the surrounding prose quotes
shipped-W figures. That section mixes two engines.

**2. The t(6) "62x residual" is not a second mechanism -- it is that mixing.**
Measured for t(6), MAF 0.01, `Z_mu` 30, corrected at 5e-8, as a multiple of nominal:

| estimator | shipped W | empirical W |
|---|---|---|
| OLS | 2.83x | 11.70x |
| GLS | **62.21x** | 275.18x |

The doc quotes 62x (shipped, GLS) in prose while the table above it gives 1.4e-05, the
empirical-W value (275x), for the same cell. Same cell, two engines, quoted together.

## The collinearity diagnosis survives -- but Pearson rho is the wrong statistic

`rho(f'/f, W_sigma2)` is invariant to the *scale* of `f'/f`; the correction is not, since
a wrongly scaled column over-subtracts proportionally. The right quantity is the
projection coefficient in the estimator's own metric, `eff = A_2 . u` with
`u(tau) = -f0'(q0) / (2 f0(q0))`. The residual O(delta^2) bias is `k delta^2 (eff - 1)`,
and `eff = 1` exactly on a Normal baseline for both estimators and both columns.

| shape | rho | sd ratio | eff_ols | eff_gls |
|---|---|---|---|---|
| Normal | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| t(6) | 0.9730 | **1.2086** | **1.1759** | 1.2959 |
| chi2(20) | 0.9335 | 1.0641 | **0.9933** | 1.2501 |
| chi2(4) | 0.6694 | 1.3035 | 0.8726 | 3.2020 |

This resolves the ordering that `rho` inverted. t(6)'s `f'/f` has nearly the right shape
but is 21% too steep; chi2(20)'s shape deficit (-6.7%) and scale excess (+6.4%) almost
cancel, so chi2(20) under OLS is accidentally near-perfect. **`eff` orders the measured
rates correctly; `rho` does not.** Documentation using `rho` as the diagnostic should use
`eff` instead.

## Is there a heavy-tail-specific mechanism?

Under OLS, no. The t(6) residual is pure O(delta^2) (higher-order term H = -0.014 against
S = 0.885), a variance-only counterfactual puts t(6) *below* nominal, and at MAF 0.50 --
where the correction is identically zero -- t(6) is perfectly calibrated (5.10e-08 at
`Z_mu` 30, 5.43e-08 at 60). The separate mechanism that hits chi2(4) at MAF 0.50
(4.478e-06 at `Z_mu` 60) does not touch t(6), so it is driven by phenotype asymmetry, not
by tails.

Under GLS there is a genuine second contributor, but it is not tail-specific: `H` becomes
large and cancels `S` for chi2(20) (1.777 and -0.857) while not cancelling for t(6)
(1.140 and +0.026), and the true-to-reported SE ratio differs on the test column (t(6)
1.0320, chi2(20) 0.9660) -- at a z = 5.45 threshold, that 6.8% gap alone is worth about 8x.

## Not run

No Monte Carlo through the shipped Stage 2 R pipeline, so the exact shipped-W numbers are
not reconciled against that path -- which `DISPERSION_SKEW_CORRECTION.md` section 9
already lists as unreconciled. MAF values other than 0.01 / 0.30 / 0.50, `Z_mu` outside
{4, 8, 16, 30, 40, 60}, and shapes beyond the four above were not tested.

Working files: `/user/work/fh6520/empw2/` (`src/`, `logs/`, `out/`).
