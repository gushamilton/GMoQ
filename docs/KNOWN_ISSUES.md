# G-MoQ: validated behaviour and known limitations

Status as of 2026-09-13, against tag `v0.2.0` (`cdbc395`), plus the
tail-calibration fix described in 2.1, which post-dates that tag.

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

and `theta/SE` is exactly `t_nu`. Raising B was never the fix; it only shrank
the gap between the wrong reference and the right one.

Stage 2 now emits a parameter p-value itself, as
`2 * pt(-|theta/SE|, df = nu)`, instead of leaving downstream code to compute
`2 * pnorm(-|beta/se|)`. `nu` travels from Stage 1 in
`<prefix>.cov.meta.json` under `jackknife`, as `nu` (currently
`n_blocks - 1`, tagged `nu_source: "blocks_minus_one"`) alongside the block
count itself. Every result row records `ref_dist` and `ref_df`, so which
reference produced a p-value is visible in the output file.

**The t identity reproduces the measured inflation from a completely
different direction.** `2*pt(-5.4513, B-1) / 5e-8`, computed from nothing but
the degrees of freedom, against the simulated values in the table above:

| B | predicted by `t_{B-1}` | measured |
|---|---|---|
| 65 | 17.14 | 17.44 |
| 200 | 2.94 | 3.02 |
| 500 | 1.57 | 1.62 |
| 1000 | 1.26 | 1.30 |
| 2300 | 1.11 | 1.15 |

Agreement to a few percent across a 35-fold range of B. The small, consistent
excess of the measured values over the prediction is what you would see if the
*effective* degrees of freedom were slightly below `B - 1`; that is being
measured separately, and `jackknife.nu` is a separate field from
`jackknife.n_blocks` precisely so a measured `nu_eff` can replace `B - 1`
without a format change.

The Stage 1 default block count is now **1000**, not 65.

**Covariance files written before the `jackknife` field existed carry no `nu`.**
Stage 2 refuses to emit a parameter p-value for them and says so, rather than
falling back to a normal reference; `theta`, `SE`, `Q` and `Q_pval` are
unaffected. Pass `nu = n_blocks - 1` to `param_gwas_from_file()` to use such a
file, or `unknown_nu = "normal"` to reproduce a pre-fix run, which is then
recorded as `ref_dist = "normal"` in the output.

**Still outstanding: `Q_pval`.** The Lugannani-Rice saddlepoint reference for
`Q` conditions on `Sigma` as known and inherits the same defect. It is not a
`pt()` substitution -- `Q` is a quadratic form in a correlated normal vector,
and the estimation of `Sigma` would have to be propagated through its
eigenvalue spectrum -- and it is untouched for now.

**No joint multi-parameter Wald test exists in the package**, so nothing
currently needs an `F_{k, nu-k+1}` (Hotelling) reference. `Q` is a lack-of-fit
statistic on `T - k` residual degrees of freedom, not a joint test of `k`
parameters. If one is added it must use `F`, not chi-square.

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
