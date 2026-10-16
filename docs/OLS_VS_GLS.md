# OLS or GLS: measured against known truth

Status: measured 2026-09-14 on 28 simulation cells x 20,000 replicates, N = 30,000.
Conclusion: **keep `ols` as the default; expose `gls` as the efficiency option; guard on lack of fit.**

## Setup

Stage 2 fits `b ~ W theta`, where `b` holds 45 per-quantile effects and `Sigma` is their
45x45 covariance. `ols` ignores `Sigma` in the point estimate and uses it only for the
sandwich SE; `gls` uses `Sigma^-1` in both. All cells use the pooled covariance.

Design: `{null, location, scale, mixed} x {normal, gamma(2), gamma(8)} x {MAF 0.30, 0.01}`,
each fitted under four models (`ADD`, `LS`, `LSK`, `TRU`) so that some cells are correctly
specified and others are not. Truth is known, so coverage is COUNTED, not inferred from
agreement between two estimates. Monte Carlo error on a 0.95 rate at 20,000 reps is 0.0015.

## Coverage of 95% intervals, each estimator against its own estimand

| truth | GLS | OLS |
|---|---|---|
| null | 0.9499 | 0.9499 |
| location | 0.9504 | 0.9484 |
| mixed | 0.9409 | 0.9406 |
| scale | 0.9373 | 0.9409 |

Both are exact under the null. Both sag to ~0.94 under a true scale effect — that is the
pooled covariance's working-independence assumption (`Var(psi_i | g_i)` free of `g_i`),
which is false at a genuine dispersion signal. **It costs both estimators equally**; GLS is
not more exposed to it.

Null Type-I error extrapolated to the genome-wide threshold, from 480,000 z-values per cell:

| cell | sd(z) | implied p at z = 5.4514 | ratio to 5e-8 |
|---|---|---|---|
| MAF 0.30, gls | 0.9998 | 4.96e-8 | 0.992 |
| MAF 0.01, gls | 1.0009 | 5.13e-8 | 1.027 |
| MAF 0.01, ols | 1.0003 | 5.05e-8 | 1.010 |

## Efficiency: empirical sd(GLS) / sd(OLS)

This is the ratio of ACTUAL variability of the estimates, not of reported standard errors.

| parameter | ratio |
|---|---|
| `beta_add` | 0.76 - 0.83 |
| `beta_sigma2` | 0.78 - 0.81 |
| `beta_kurt` | 0.80 - 0.82 |
| `beta_scale` | 0.88 - 0.90 |
| `beta_mu` | 0.97 - 1.00 |

GLS estimates are **17-24% less variable** on the parameters that carry the science. In only
3 of 224 cells is GLS less efficient, all at ratio 1.002 — a tie, not a loss.

## Where OLS wins: when the fitted model is wrong

Binned by lack of fit, measured as `max|residual| / median SE per tau`:

| lack of fit | cells | sd(GLS)/sd(OLS) | OLS cov95 | GLS cov95 | OLS \|bias\| | GLS \|bias\| |
|---|---|---|---|---|---|---|
| fits well (< 0.5) | 184 | 0.917 | 0.9456 | 0.9459 | 0.021 | 0.040 |
| slight (0.5-1) | 14 | 0.880 | 0.9406 | 0.9414 | 0.051 | 0.089 |
| moderate (1-3) | 10 | 0.890 | 0.9408 | 0.9415 | 0.020 | 0.048 |
| **badly (> 3)** | 16 | **0.750** | **0.9448** | **0.9359** | 0.033 | **0.128** |

When the model fits, the two are indistinguishable on coverage and GLS is ~8% tighter.
When it fits badly, GLS's coverage falls below OLS's and its bias is **4x worse** — 0.128
against 0.033 empirical standard deviations. GLS looks *most* efficient there (0.750), which
is the trap: the estimates cluster tightly around the wrong answer.

This is the textbook behaviour. GLS is efficient *given the model*; under misspecification
the `Sigma^-1` weighting concentrates on the directions that amplify the misfit rather than
averaging it away. OLS's equal weighting is more forgiving.

Worst individual cells, all GLS, all rare-MAF or skewed with a scale effect:

| cell | fit | param | cov95 | SE ratio | bias (sd) | lack of fit |
|---|---|---|---|---|---|---|
| `scale_maf0.01_gamma2` | ADD | `beta_add` | 0.879 | 0.794 | -0.078 | 3.03 |
| `scale_maf0.3_gamma2` | ADD | `beta_add` | 0.916 | 0.931 | -0.351 | 13.38 |
| `scale_maf0.3_gamma2` | LS | `beta_sigma2` | 0.919 | 0.961 | 0.415 | 4.08 |

That is precisely the regime where a vQTL claim would live, which is the main argument
against making GLS silent and automatic.

## Recommendation

- **`ols` stays the default.** Barely worse than GLS when the model is right (0.9456 vs
  0.9459 coverage), meaningfully safer when it is not.
- **`gls` is a supported option** for analyses where the model fit has been checked, worth
  17-24% on the parameters that matter.
- **A lack-of-fit guard is needed either way.** `Q` already measures exactly the quantity that
  separates the safe cells from the dangerous ones, so the diagnostic is computable per
  variant without knowing the truth.
- Independently: the GLS path has **no guard against a numerically singular `Sigma`** — it
  emits a z-score of 498 with every diagnostic counter clean. That must be fixed before GLS
  is recommended for production use.

## Reproducing

The campaign's sources were `src/cells.py` (design), `src/sim_stage1.py`,
`src/stage2_driver.R`, `src/truth.R`, `src/summarise.R`, `src/tail.R`.
Per-cell results in `out/stage2_cell*.rds`; the assembled table in `out/summary.tsv`.
