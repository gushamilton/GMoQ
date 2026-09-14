# G-MoQ: validated behaviour and known limitations

Status as of 2026-09-13, against tag `v0.2.0` (`cdbc395`), plus the
tail-calibration fix described in 2.1 and 2.2, which post-dates that tag and
rests on `fungwas_paper/docs/BP1_effective_df.md` (2026-09-14).

This records what has been tested, what was found, and what is not yet
established. Full evidence lives in the analysis tree under
`fungwas_paper/docs/`; this file is the summary that travels with the package.

---

## 1. What has been validated

| Question | Answer | Evidence |
|---|---|---|
| Is the block-jackknife covariance biased? | **No** — ratio to its true target 1.0007 / 1.0016 / 1.0006 | 8 configurations, R=100, N=50k |
| Does the jackknife recover the correlation structure? | **Yes** — max deviation 0.0057, rms 0.0023 | incl. cases far from the analytic form |
| Are Stage 2 parameter estimates biased? | **No** — -0.0004 +/- 0.0022 (null), +0.0030 +/- 0.0030 (0.02 SD/allele) | measured against the RIF projection estimand |
| Is Type-I error correct at conventional thresholds? | **Yes** — 1.085x at p=0.05 | |

Software correctness is separately covered by the test suite: 43 Stage 1
tests, 18 Stage 2 tests.

## 2. The one substantive defect: tail calibration

`SE_hat` is unbiased but **9% noisy** at `--blocks 65`. A noisy denominator is
punished severely in the tail, because the ratio of nominal to actual error
rate grows with the square of the threshold's z-score.

Type-I error inflation, nominal vs actual:

| nominal p | `--blocks 65` | `--blocks 200` | `--blocks 500` |
|---|---|---|---|
| 0.05 | 1.085 | 1.026 | 1.009 |
| 1e-4 | 2.413 | 1.367 | 1.139 |
| 1e-6 | 7.141 | 2.089 | 1.374 |
| **5e-8** | **17.44** | **3.02** | **1.62** |

**At production settings a nominal 5e-8 threshold delivers approximately
8.7e-7.** Consistent across three data-generating mechanisms (17.4 / 18.1 /
16.7). The point estimates and the covariance are both fine; only the
studentised tail is affected.

**The 5e-8 row is an extrapolation**, computed via the factorisation
`z = t/r` with `Var(theta|G)` known analytically. The normality of `t` is
verified to 4.5 sigma; genome-wide significance is 5.45 sigma. The estimator
reproduces direct counts to 0.4% at 0.05, 0.2% at 1e-3 and 2.4% at 1e-4.

**More blocks helps but does not fully fix it.** `sd(r)` falls as
`1/sqrt(B-1)`. From `ln(ratio) ~ 435 * sd(r)^2`, reaching 1.1x needs
B ~ 2,300 and 1.05x needs B ~ 4,500. Cost beyond B=500 is not measured;
the linear cost fit in the analysis tree is not verified that far out.

### 2.1 Fixed for the parameter Wald test

The table above is what happens when `theta/SE` is read against a **normal**
reference. It should not be. `Sigma` is estimated, so for a fixed contrast `a`,

    a' Sigma_jk a / a' Sigma a  ~  chi2_nu / nu

and `theta/SE` is exactly `t_nu`. Raising `B` was never the fix; it only shrank
the gap between the wrong reference and the right one.

Stage 2 now emits a parameter p-value itself, as
`2 * pt(-|theta/SE|, df = nu)`, instead of leaving downstream code to compute
`2 * pnorm(-|beta/se|)`.

**`nu` is computed per SNP and per parameter, not read as a constant.** A
constant was tried and rejected on measurement (`BP1_effective_df.md`,
2026-09-14): `nu` is not `B - 1` (801.8 against 999 at `B = 1000`, `N = 50k`),
it varies with the genotype (530 at MAF 0.01 to 979 at MAF 0.50 at `B = 1000`),
and it varies with the contrast (566 for a scale-type parameter against 802 for
a location-type one, a gap that *widens* with `B`). The stable object is the
law, fitted over 28 cells to a median error of 0.8%:

    1/nu = 1/(B - 1) + c / N_eff        c = 0.88 location, 2.62 scale

The second term is a Satterthwaite penalty for unequal block contributions.
`N_eff = D^2 / sum_i g_ri^4` collapses to `D = sum_i g_ri^2` because a
genotype's fourth central moment equals its second (`N_eff/D` = 1.003 measured
on real dosages), and `D` is the denominator Stage 1 already forms. Stage 1
emits it as the **`d_resid`** column; `B` comes from `<prefix>.cov.meta.json`
under `jackknife`, together with the coefficients in `jackknife.df_model`.

**`c` is measured, not derived.** The derivation exists — the per-individual
score is a step function of the phenotype rank, so `c = (1 - 1/B) + gamma_e/2`
with `gamma_e` its excess kurtosis — but on the production grid it gives 0.75
and 1.79 against measured 0.88 and 2.62: right mechanism, right ordering, wrong
magnitude by 15% and 32%, and wrong **anti-conservatively**. It also needs
`f_tau`, which Stage 1 does not emit. The working is in
`stage2/R/tail-reference.R`. `c` is therefore selected by the flatness of the
contrast, `|sum(a)| / (sqrt(T) ||a||)`, and **anything that is not clearly flat
gets the scale value, which is the conservative one**. The constants were
calibrated for `W = [1, q_tau - mean(q_tau)]` on a 45-tau grid; **a mixture or
GxE design gets that conservative default rather than a calibrated `c`**, and
`df_coef` exists to override it.

Every result row records `ref_dist`, `ref_df_source` and a `<param>_df` column,
so which reference produced which p-value is visible in the output file.

**What this is worth.** Actual Type-I error over nominal at `5e-8`, at
production `N = 461,320`, `B = 1000`:

| reference used | `mu`, MAF 0.01 | `mu`, MAF 0.05 | `mu`, MAF 0.50 | `scale`, MAF 0.01 |
|---|---|---|---|---|
| normal (the old behaviour) | 1.288 | 1.266 | 1.261 | 1.344 |
| `t_{B-1}` (naive) | 1.023 | 1.005 | 1.001 | **1.069** |
| per-SNP `nu` (what ships) | 1.000 | 1.000 | 1.000 | 1.000 |

The residual case that justifies the per-SNP form over a naive `t_999` is the
**scale parameter at low MAF**. At `N = 50,000` the naive form is far worse
throughout, because `N_eff` scales with `N`.

**The `t` identity also reproduces the earlier simulation exactly.**
`2*pt(-5.4513, B-1) / 5e-8` predicts 17.14 / 2.94 / 1.57 / 1.26 / 1.11 at
`B` = 65 / 200 / 500 / 1000 / 2300 against the simulated 17.44 / 3.02 / 1.62 /
1.30 / 1.15 in the table above. The small consistent excess of measured over
predicted is the `nu_eff` deficit, which is what the law now corrects.

The Stage 1 default block count is **1000**, not 65. At production that
delivers `nu` between 911 and 995 for location contrasts; **`B = 1108` would
deliver a true `nu` of 1000 at MAF 0.01**, and there is no case for going much
beyond that, since `nu` is bounded above by `N_eff/c` however many blocks are
used.

**Covariance files written before the `jackknife` field existed carry no block
count.** Stage 2 refuses to emit a parameter p-value for them and says so,
rather than falling back to a normal reference; `theta`, `SE`, `Q` and `Q_pval`
are unaffected. Stage 1 files written before `d_resid` existed fall back to a
flat `nu` computed at `N_eff = 9134` — the production MAF 0.01 value, so
conservative at commoner variants — with a warning and
`ref_df_source = "flat_fallback_n_eff"`. That costs at most 2.4% of nominal
Type-I for a location parameter and 6.5% for a scale one.

**Still outstanding: `Q_pval`.** The Lugannani-Rice saddlepoint reference for
`Q` conditions on `Sigma` as known and inherits the same defect. It is not a
`pt()` substitution — `Q` is a quadratic form in a correlated normal vector,
and the estimation of `Sigma` would have to be propagated through its
eigenvalue spectrum — and it is untouched for now. Note also that
`.q_pvalue_calibrated` already computes a `df_eff` of its own; that is the df
of the quadratic form, not of the variance estimator, and the interaction
between the two has not been examined.

**No joint multi-parameter Wald test exists in the package**, so nothing
currently needs an `F_{k, nu-k+1}` (Hotelling) reference. `Q` is a lack-of-fit
statistic on `T - k` residual degrees of freedom, not a joint test of `k`
parameters. If one is added it must use `F`, not chi-square, and with a
per-SNP `nu`.

### 2.2 Do not studentise the per-tau SEs against the block count

`se_tau_*` in the Stage 1 TSV are jackknife estimates, so a per-tau Wald test
is `t` — and its effective degrees of freedom are **far** below the block
count, because the RIF residual at an extreme tau is a scaled centred
Bernoulli whose excess kurtosis is `(1 - 6t(1-t)) / (t(1-t))`: **-2 at the
median but +45 at `tau = 0.02`**.

Measured at `N = 50,000` (`BP1_effective_df.md` 3):

| `B` | `tau=0.02` | `tau=0.107` | `tau=0.50` | `tau=0.98` |
|---|---|---|---|---|
| 65 | 42.5 | 59.0 | 62.9 | 42.8 |
| 500 | 100.9 | 306.5 | 454.1 | 96.6 |
| **1000** | **111.7** | **445.4** | **834.8** | **106.8** |
| 2000 | 118.8 | 572.5 | 1444.8 | 112.8 |

and at `tau = 0.02`, `B = 1000`, by MAF: **37.8** (MAF 0.01) rising to 312.9
(MAF 0.50).

**Per-tau `nu` at the extreme quantiles barely improves with `B` at all** —
100.9 at `B = 500` to 118.8 at `B = 2000`, an 18% gain for four times the
blocks. Stage 2 does not emit per-tau p-values, so this is not a live defect;
it is a trap for anyone who studentises the Stage 1 output directly. If you do,
use neither a normal, nor `B - 1`, nor the `<param>_df` Stage 2 reports for a
fitted parameter.

These are **measurements at `N = 50,000`, not a law**: the second-order
kurtosis expansion under-predicts the magnitude at `tau = 0.02` by about 2.4x,
and per-tau `nu` at production `N` is unmeasured.


## 3. The analytic covariance

The per-SNP cross-tau *correlation* matrix is, under stated conditions, the
analytic Brownian-bridge form:

    R[s,t] = (min(s,t) - s*t) / sqrt(s*(1-s)*t*(1-t))

Where this holds, `Sigma` can be reconstructed as
`diag(se_tau) %*% R %*% diag(se_tau)` — no stored covariance, positive
definite by construction, no block count, no seed, and `r = 1` identically,
which gives exact tail calibration.

**Validity is governed by `sd_i(F_i(q_tau))`**, computable from a phenotype
and its covariates with no GWAS. Predicted vs measured on real traits:

| trait | predicted | measured |
|---|---|---|
| `bmi` | 0.0199 | 0.0201 |
| `log_bmi` | 0.0209 | 0.0143 |
| `weight` | 0.1459 | 0.1387 |
| `weight_height` | 0.1501 | 0.1519 |
| `testosterone_shrink_0.5` | 0.0033 | 0.0145 |

**The rule is counterintuitive and easy to get backwards.** Heterogeneity that
is *adjusted for* creates the departure; heterogeneity from an *omitted*
variable is absorbed into the between-genotype term and restores the bridge.
The same bimodal phenotype gives deviation 0.61 with sex as a covariate and
0.028 without.

Consequently:

- Safe: `bmi`, `log_bmi`, and testosterone **as currently analysed**.
- Unsafe: `weight`, `weight_height` (height covariates drive `sd_i(F_i)` from
  0.045 to 0.197).
- **Contingent on the covariate set.** Testosterone is safe only because
  production does not adjust it for sex. Adding sex would break it. Re-check
  this diagnostic whenever covariates change.

## 4. Model selection

`src/89_tau_cov_model_compare.R` selects via `AIC = Q + 2k` using the raw `Q`
statistic, not its p-value. Therefore `q_method` (the null-distribution
correction added in v0.2.0) does **not** affect model selection.

The criterion itself is mis-scaled. Under a known additive truth, recovery of
the correct model:

    AIC = Q + 2k (as used)      44.0%
    naive chi-sq p-rank         40.4%
    BIC                         59.4%
    saddlepoint p-rank          61.0%
    scale-corrected Q/a + 2k    87.8%

A separately-tested hypothesis — that skewness of the phenotype manufactures
this result — is **not** substantiated: the effect is 4 points and in the
wrong direction (a symmetric marginal recovers the additive model *less*).
The `additive` weight matrix is W = 1.0 at all 45 taus by construction, so a
raw-scale additive effect produces a constant quantile shift regardless of
skew.

## 5. Not established

- Behaviour at biobank N. All validation is at N = 50,000.
- Any comparison against direct model fitting (accuracy or speed). No
  comparator has been run.
- Any speed claim against an alternative method.
- Robustness to LOCO residualisation, which production uses and the
  simulations deliberately omit.
- Calibration at 5e-8 by direct measurement rather than extrapolation.
- Per-tau effective degrees of freedom at production `N`. The numbers in 2.2
  are measured at `N = 50,000`, where a `B = 1000` block holds ~1 individual
  below `q_0.02`; at production it holds ~9, which is a different regime.
- `nu_eff` under the alternative, and for weight matrices other than
  `[1, q_tau - mean(q_tau)]`. Both `c` values were calibrated on that design.
- Interaction between the covariance degrees of freedom and the `df_eff` that
  `.q_pvalue_calibrated` computes for the quadratic form `Q`.

## 6. Traps for anyone extending this work

1. **Conditional vs unconditional variance.** The jackknife deletes
   individuals, so it targets the *unconditional* variance. A design that
   fixes genotypes measures the conditional one. Comparing the two produced
   three separate false "defects" during validation, each reversed on
   re-checking. Always check against the unconditional target first.
2. **`~/.Rprofile` may override `R_LIBS_USER`.** A bare `.libPaths(...)` call
   replaces the search path rather than prepending, so setting `R_LIBS_USER`
   silently loads the wrong package. Assert the loaded package path; do not
   trust the environment variable.
