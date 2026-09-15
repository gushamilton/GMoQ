# Why GLS coverage gets worse as the sample grows, and OLS does not

Status: derived and confirmed numerically 2026-09-15. Companion to `N_SCALING.md`,
which measured the effect, and `OLS_VS_GLS.md`, which recorded the decision.

## The observation

Coverage of a nominal 95% interval on the dispersion parameter, under a true scale
effect, MAF 0.30, 720,000 replicates (Monte Carlo error 0.0011-0.0018):

| | N = 30,000 | N = 100,000 | N = 461,320 |
|---|---|---|---|
| OLS | 0.9440 | 0.9434 | 0.9413 |
| GLS | 0.9316 | 0.9064 | **0.8437** |

Under the null both hold at every N. `N_SCALING.md` established why misfit matters
more at large N: the absolute lack of fit is a property of the data-generating
process and is the same at every N, while standard errors shrink, so misfit measured
in SE units grows as `sqrt(N_eff)`. That does not explain why GLS is hurt so much
more than OLS. Both face the same misfit. This document explains the asymmetry.

## The mechanism

Write `b = W theta + r`, with `r` the systematic residual that does not lie in the
column space of `W`. Both estimators are biased by it:

    bias(OLS) = (W'W)^-1 W' r
    bias(GLS) = (W'Sigma^-1 W)^-1 W' Sigma^-1 r

and both report a standard error:

    SE(OLS) = sqrt(diag[ (W'W)^-1 W' Sigma W (W'W)^-1 ])     sandwich
    SE(GLS) = sqrt(diag[ (W'Sigma^-1 W)^-1 ])                model-based

Two things differ, and they compound.

**The bias is amplified.** `Sigma` is the cross-quantile covariance of the quantile
process, so it has Brownian-bridge structure: a decaying spectrum with no gap.
`Sigma^-1` therefore weights the *small*-eigenvalue directions most heavily. The
misspecification residual carries **9-15x disproportionate mass in exactly those
directions**, and `Sigma^-1` amplifies them **15.7x**. OLS, projecting in the
identity metric, gives every tau uniform weight and dilutes those directions instead.

**The standard error does not know.** `SE(GLS)` is the true sampling variance only
when the model is correct. `SE(OLS)` is a sandwich built from the observed
variability of `b` and is blind to the *direction* of the misfit, so it stays valid
under misspecification (White 1982).

Together these give a bias-to-SE ratio **3.6x larger for GLS, constant in N**:

| N | OLS bias/SE | GLS bias/SE | ratio |
|---|---|---|---|
| 30,000 | 0.073 | 0.265 | 3.63 |
| 461,320 | 0.287 | 1.038 | 3.62 |

Both ratios grow as `sqrt(N)` because the bias is fixed and the SE shrinks. GLS
simply starts 3.6x closer to the point where coverage fails, and reaches it first.

## A sandwich SE on the GLS path does not fix it

The obvious repair is to give GLS a sandwich standard error, trading its efficiency
claim for honest coverage. It is a no-op:

    (W'S^-1 W)^-1 W'S^-1 . S . S^-1 W (W'S^-1 W)^-1  =  (W'S^-1 W)^-1

The GLS sandwich collapses algebraically to the model-based form already reported.
GLS is *already* sandwich-like in structure.

**The problem is bias, not variance, and no standard error repairs a biased
estimator.** Note also that OLS's robustness does not come from its sandwich form
as such; it comes from the identity metric spreading weight uniformly across taus,
so the residual is averaged down rather than amplified.

## What this is, generically

Standard behaviour whenever an estimator is biased under misspecification and its
reported SE assumes the model is correct. GLS is efficient when the model *is*
correct; when it is not, and the residual aligns with the directions `Sigma^-1`
amplifies, the efficiency gain becomes a liability. What is specific here is the
quantile covariance's eigenstructure, which puts `Sigma^-1`'s largest weights
exactly where the residual lives.

## The practical rule

`Q` is the lack-of-fit statistic Stage 2 already computes:

    Q = (b - W theta_hat)' Sigma^-1 (b - W theta_hat)

**The raw ratio `Q/(T-K)` is the wrong quantity on the GLS path, and the
thresholds below do not transfer to it unrestated.** They were measured on the
pooled covariance path, where `nu = n_samples` is so large that `Sigma` is
effectively known. On the jackknife path `Sigma` is estimated on a finite `nu`,
and `Q_gls` is then Hotelling's `T^2`, whose null mean is `nu(T-K)/(nu-(T-K)-1)`
and not `T - K`. Measured null mean of `Q/(T-K)` under a CORRECTLY SPECIFIED
model, T = 45, K = 2 (400,000 replicates per cell):

| `nu` | `E[Q]/(T-K)` | null fraction with `Q/(T-K) > 1.1` | `> 1.5` |
|---|---|---|---|
| infinite | 1.001 | 0.304 | 0.019 |
| 1000 | 1.045 | 0.378 | 0.034 |
| 260 | 1.203 | 0.611 | 0.148 |
| 100 | 1.787 | 0.937 | 0.680 |
| 60 | 3.753 | 0.999 | 0.984 |

At `nu = 260` -- the production scale-contrast value at MAC 922 -- the raw rule
would send 61% of perfectly well specified variants to the `> 1.1` bucket.

**Use the calibrated ratio instead.** Stage 2 now returns `Q_df_eff` and
`Q_scale` on the GLS path, chosen so that `Q_df_eff * Q_scale` is the null mean
of `Q`. The quantity that carries the same meaning at every `nu` is

    R = Q / (Q_df_eff * Q_scale)        null mean 1 on BOTH estimator paths

whose null exceedance is 0.30 (`R > 1.1`) and 0.02-0.03 (`R > 1.5`) across the
whole range above. On `R` the original thresholds keep their intent:

| `R` | guidance |
|---|---|
| < 1.1 | OLS and GLS both safe |
| 1.1 - 1.5 | OLS safer; GLS acceptable only at large effect size |
| > 1.5 | OLS only, particularly for dispersion and shape parameters |

At production N the margin is thin: even `R = 1.05` is marginal for GLS,
because the bias-to-SE threshold for coverage failure falls as N rises.

Note what those exceedance rates mean: `R > 1.1` fires on 30% of correctly
specified variants even with `Sigma` known. These are screening cuts, not
tests. **`Q_pval` is the calibrated quantity**; prefer it when a rate matters.

**The rule is now enforceable on the GLS path.**
`param_gwas(tau_cov_estimator = "gls")` returns `Q_pval`, `Q_df_eff` and
`Q_scale`, referenced against `F_{T-K, nu-(T-K)+1}` rather than against the OLS
eigenvalue spectrum. GLS's `Q` reads 30-49 against OLS's 8-18 on identical
inputs because its null mean is `T - K = 43` while OLS's is
`tr(D^-1/2 M Sigma M' D^-1/2)`, far below `T - K` once the taus are correlated;
they are different statistics, and only the GLS one is pivotal. The derivation
is at the foot of `stage2/R/tail-reference.R`; `N_SCALING.md` section 7 has the
measured power.

## Consequence

GLS is not repairable here by changing the standard error. It is conditionally
inapplicable, and the condition is measurable per variant. OLS remains the default.
