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

*(This section is about the explicitly selected block-jackknife covariance.
The pooled covariance is the Stage 2 default; it does not have this defect:
its degrees of freedom are 5e04 to 9e05 rather than 62, so its reference is
normal.)*

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

**Still outstanding on the OLS path: `Q_pval`.** The Lugannani-Rice
saddlepoint reference for the OLS `Q` conditions on `Sigma` as known and
inherits the same defect. It is not a `pt()` substitution — the OLS `Q` is a
quadratic form in a correlated normal vector weighted by `D^-1`, and the
estimation of `Sigma` would have to be propagated through its per-SNP
eigenvalue spectrum — and it is untouched. Note also that
`.q_pvalue_calibrated` already computes a `df_eff` of its own; that is the df
of the quadratic form, not of the variance estimator, and the interaction
between the two has not been examined.

**Resolved on the GLS path (2026-09-15).** The GLS `Q` weights by `Sigma^-1`
rather than `D^-1`, which collapses that eigenvalue spectrum to a flat one and
makes the statistic pivotal, so the estimation of `Sigma` CAN be propagated
exactly. With `Sigma` known `Q_gls ~ chi2_{T-K}` exactly; with `Sigma`
estimated on `nu` degrees of freedom it is Hotelling's `T^2` on the residual
subspace,

    Q_gls ~ (nu m / (nu - m + 1)) F_{m, nu - m + 1},    m = T - K

Verified uniform through `param_gwas()` over 5,000,000 correctly specified
replicates down to `p = 1e-6`. The chi-square reference that a naive port of
the OLS path would give is anti-conservative at every finite `nu`: 0.0825 at
nominal 0.05 when `nu = 1000`, 0.250 when `nu = 260`, 0.991 when `nu = 60`.
Derivation and measurements at the foot of `stage2/R/tail-reference.R`.

**A joint multi-parameter Wald test still does not exist in the package**, but
the Hotelling reference is no longer hypothetical: the GLS lack-of-fit
statistic IS one, on `m = T - k` residual dimensions rather than `k` parameter
dimensions. `nu` must exceed `m`, or the residual block of the estimated
`Sigma` is singular and `Q` is degenerate rather than merely mis-referenced
— with `T = 45` that rules out `B = 65` outright, and Stage 2 refuses there
rather than returning a number. If a joint Wald test is added it must use `F`,
not chi-square, and with a per-SNP `nu`.

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
3. **A projector is only a projector on the rows it was built from.** Anything
   used as `M - Q (Q' M)` must have `Q` built by QR of the covariate matrix
   restricted to exactly the rows being swept. Row-subsetting a basis built on a
   larger sample destroys orthonormality and biases the point estimate while
   leaving the standard error intact, so the output looks entirely normal. This
   happened -- see section 9 -- and `core.check_orthonormal_basis()` now enforces
   it at every kernel entry point.

---

## 7. Residual calibration gap for scale parameters (as of 2026-09-14)

The `t`-reference fix removes most, but not all, of the tail miscalibration.
Measured on a single simulation containing both null and causal variants
(200 replicates, 2.24M null tests), actual/nominal Type-I error at B=1000:

| parameter | p=0.05 | p=1e-3 | p=1e-4 |
|---|---|---|---|
| location (`mu`) | 1.001 | 1.042 | 1.125 |
| scale | 1.006 | **1.074** | **1.357** |

Location parameters are within 1.3 standard errors of nominal everywhere.
Scale parameters sit 2.6-3.1 standard errors above it; pooling two independent
arms gives 1.066 +/- 0.023, i.e. 2.9 SE. That is real, not Monte Carlo noise.

For context the fix removes 93% of the gap: at p=5e-8 the sequence is
normal 1.566, `t_{B-1}` 1.251, per-SNP `nu` 1.037.

The cause is not the degrees-of-freedom law, which was checked independently
on the same panel and is correct to within 1% (implied `c` 0.86 against 0.88
shipped for location, 2.72 against 2.62 for scale). The residual is that `r =
SE_hat/SE_true` is slightly heavier-tailed than `sqrt(chi2_nu/nu)` for a scale
contrast, which the `t` reference assumes exactly.

**Consequence**: scale and shape parameters are mildly anti-conservative in the
tail. Treat a scale-parameter p-value near a threshold with more caution than a
location-parameter one. This was not measured at production sample size, where
the margin is expected to be smaller.

## 8. Kernel density bandwidth (as of 2026-09-14)

**Scope update, 2026-09-15:** the Q and model-selection comparisons below do
not establish unbiased non-null parameter estimates or calibrated confidence
intervals. The later [paired density diagnosis](../tests/validity_pilot/evidence/DENSITY_DIAGNOSIS.md)
isolated smoothing bias in the Normal variance projection at both N=10,000
and N=100,000. Cancellation in a single-quantile zero-effect Wald statistic
does not remove that non-null estimation problem. No density repair has yet
been validated for production; consult the [experiment register](../experiments/README.md)
before repeating or extending these investigations.

`stage1/fungwas_stage1/core.py` selects the bandwidth as
`h = 1.06 * np.std(Y) * N**(-0.2)` -- Silverman's rule without a robustness
guard. For bimodal phenotypes the standard deviation is inflated by the
separation between modes, making `h` 3.2-3.4x too wide.

Measured distortion in `f_hat` on real phenotypes: `bmi` 0.41%, `log_bmi`
0.17%, `weight` 0.64%, `jcv_log_all` 2.78%, **`testosterone` 15.85%**,
**`mcv_log_all` 22.07%**.

Downstream consequence is small. Because `b_tau` and `se_tau` are both
proportional to `1/f_hat`, the density cancels from any single-quantile Wald
test. It does not cancel from `Q`. Recomputing `Q` under a Sheather-Jones
bandwidth on 8,000 real SNPs per trait moves the median by +0.08 (bmi/additive)
and +0.003 (testosterone/skew), but +0.82, or +6.7%, for `mcv_log_all`. Model
selection over the four BMI models flips 34 of 7,992 SNPs (0.43%), every one a
tie.

The direction is **conservative**: correcting the bandwidth makes `Q` larger,
so the current rule does not manufacture apparent lack of fit.

**Do not** replace it with `0.9 * min(sd, IQR/1.34)`. For the bimodal traits the
mode separation inflates the IQR more than the standard deviation, so the
`min()` never binds and the rule degenerates to a flat rescaling -- inert
exactly where the problem is. A Sheather-Jones plug-in is the recommended
replacement at the next re-scan; it removes 92% of the bimodal bias and reduces
to approximately the current bandwidth for near-normal traits.

`mcv_log_all` (N=9,365) is the one trait where a shape conclusion could move.

## 9. Covariate basis vs sample match (defect, fixed 2026-09-14)

**Affects every version up to and including `v0.3.0`, on the non-LOCO Stage 1
path only. No published result is affected.**

Up to `v0.3.0`, `stage1/fungwas_stage1/cli.py` built the QR of the covariate
matrix and the RIF on the full keep sample, then row-subset both to the samples
that matched the BGEN. The kernels use that basis as a projector,
`g_resid = g - Q (Q' g)`, which is a projection only if `Q'Q = I`. Columns
orthonormal over `N` rows are not orthonormal over `n < N` of them, so the
residualisation was partial and the covariate-parallel part of the genotype that
survived it leaked into the numerator of every per-tau slope.

Measured on the UK Biobank BMI inputs (486,254 keep rows, 461,320 in the chr22
genotype release):

| | `max \|Q'Q - I\|` |
|---|---|
| basis row-subset from the keep sample | **9.69e-01** |
| basis built on the matched rows | 8.44e-15 |

This is far worse than the ~5% that `n/N = 0.949` suggests: `diag(Q_sub' Q_sub)`
falls to 0.031, i.e. the PC1-PC3 columns retain 3-9% of their norm. The genotype
release is ancestry-restricted while the keep file is built on phenotype
availability alone, so the 24,934 unmatched samples carry most of the ancestry
signal (PC1 differs by 3.57 sd between the two groups).

**The damage is a bias, not extra variance.** Non-LOCO Stage 1 over 17,550 real
chr22 HapMap3 SNPs x 45 taus, before the fix against a correctly rebuilt
reference: median `|b - b_ref|/se_ref` **0.139**, max **2.98**, mean
`(z - z_ref)` **+0.173**, and the per-tau signed mean is **positive at all 45
taus** (+0.033 to +0.277). Standard errors barely move (median
`|se/se_ref - 1|` 4.9e-3), so the point estimate shifts while the uncertainty
does not admit it -- nothing downstream can detect it.

**Why no published result is affected.** Every paper Stage 1 run passes
`--loco-preds`, whose branch rebuilds both the RIF and the QR on the matched
subset before the first variant is processed. That was confirmed across all 247
Stage 1 logs, and the fix was verified to leave that path byte-identical.

**The fix** reorders `cli.py` to match first and then build both objects on the
matched rows -- i.e. the non-LOCO path now calls the same
`compute_chr_rif_residuals` the LOCO path has always used. It also closes a
second, much smaller channel for free: `q_tau`, the Silverman bandwidth and
`f_hat` are now computed on the analysed sample rather than the keep sample
(separately measured at mean `dz` = -0.0008, so it was never the problem).

**The guard.** `core.check_orthonormal_basis()` now measures `max |Q'Q - I|` at
every kernel entry point and raises above `core.ORTHONORMAL_BASIS_TOL = 1e-8`.
The tolerance is bracketed by measurement, not chosen by taste: float64 QR on
the real matrix gives 2.3e-14, while dropping even a single row gives 2.3e-6.
1e-8 therefore sits 5.6 orders above the noise and 2.4 orders below the smallest
realisable instance of the defect, leaving no basis that is numerically sound
but fails and none that is mis-built but passes. It raises rather than warns
because the failure is silent by construction and no legitimate caller passes a
non-orthonormal basis. Cost is about 0.4% of kernel time.

Full analysis and verification: `fungwas_paper/docs/BP1_sample_match.md`.

---

## 10. The pooled cross-quantile covariance (new path, 2026-09-14)

**The pooled covariance is the Stage 2 default, selectable as
`cov_source = "pooled"`. The block jackknife is unchanged and remains
available explicitly as `cov_source = "jackknife"`; every published result
used it and this is not a one-way door.**

Stage 1 now also writes, per trait per chromosome,

```
Sigma_0 = Psi_r' Psi_r / n          Psi_r = the covariate-residualised RIF
```

into `<prefix>.sigma0.tsv.gz` plus a `<prefix>.sigma0.meta.json` sidecar, and
Stage 2 forms the per-variant covariance as `Sigma_v = Sigma_0 / d_resid`.
`Sigma_0` contains no genotype: the cross-quantile correlation is a property of
the phenotype, the covariates and the tau grid, and the only per-variant
quantity is the scalar `d_resid = sum_i g_resid_i^2`, which the TSV has carried
since `v0.3.0`. `--no-sigma0` suppresses the sidecar.

### What it costs and what it buys

| | jackknife (`B = 65`, int8) | pooled |
|---|---|---|
| storage, chr22 `bmi` | 227.5 gzipped **B/variant** | 13.3 **KB/chromosome** (0.76 B/variant over 17,550) |
| Stage 1 compute | (the block-score kernel) | **0.06 s** per chromosome, 0.17% of a 35.6 s run |
| diagonal ratio to an exact per-individual HC0 sandwich | 0.993 (p5 0.727, p95 1.313) | **1.00001** (p5 0.979, p95 1.022) |
| max cross-tau correlation deviation from that sandwich | 0.306 | **0.008** |
| Stage 2 SE, median ratio to exact | 0.992-1.002 (p5-p95 0.84-1.15) | **0.9997-1.0004** (p5-p95 ±1-4%) |
| `df_eff` ratio to exact | 0.911-0.958 | **0.9998-1.0002** |
| positive semi-definite | in exact arithmetic; **not after int8 quantisation** | by construction (a Gram matrix) |
| GLS solved, 1,986 chr22 variants | **0** | **1,986** |
| seed / block count | required | none |

Measured on real chr22 `bmi` and `weight`, N = 461,320, 1,000-1,002 variants
per trait; full evidence in
`fungwas_paper/docs/BP1_pooled_implementation.md`.

### Degrees of freedom

Section 2's `t_nu` machinery exists because the jackknife covariance is
estimated from `B` replicates. `Sigma_0` is estimated from `n` individuals and
is shared across variants; its Satterthwaite degrees of freedom, computed
exactly from the influence functions for the paper's own Stage 2 contrasts,
are **50,439 to 945,964**, against 62 for the jackknife at `B = 65`. At
`p = 5e-8` the `t` p-value is 1.0002-1.0047 times the normal one, so the
reference degenerates to normal. Stage 2 records `nu = n_samples` with
`ref_df_source = "pooled_sigma0_n_samples"` and the law of section 2 is left
untouched for the jackknife path, which still needs it.

### Where it does not hold

`Sigma_v = Sigma_0 / d_resid` assumes `Var(psi_i | X_i, g_i)` is free of
`g_i`. That is exact under the null and false at a variant that genuinely
changes the **spread** of the phenotype. The drift is conservative (the pooled
standard error is the larger one) and only becomes material well past the
strongest reported vQTL effects, but it is a bias, not noise. A lead locus with
a strong dispersion signal should have its covariance recomputed - the
jackknife path is still there.

At low minor allele count the estimator is fine and the *reference* is not: on
both traits and over five MAC decades the pooled estimator's distance from the
exact sandwich is 0.87-1.06 times the sandwich's own sampling error, and the
two obey the same `MAC^-0.48` law. Below about MAC 60 no per-variant estimator
has enough information to be a reference at all.

### The refusals

`Sigma_v = Sigma_0 / d_resid` is only the covariance of `beta_tau` if `Sigma_0`
came from the same phenotype, covariate set, tau grid, LOCO predictor and
sample as the scan. A mismatch leaves the point estimates untouched and the
standard errors plausible, so nothing downstream can detect it; the measured
size of the equivalent failure for a covariance blind to the covariates is
+1.7% on the `bmi` variance and **+25.6% on `weight`**. Stage 1 therefore
writes an `analysis` block into `.cov.meta.json` (schema version 4) carrying a
fingerprint over all five, and **Stage 2 stops rather than proceeds** on any
disagreement, on a missing `d_resid` column, on a chromosome with no matrix, or
on a sidecar schema it does not understand. `require_covariates` and
`require_n_samples` let a caller assert independently. Same precedent as
refusing to emit a p-value when `nu` is unknown.

### Not done

Stage 2 still expands `Sigma_0` into the same per-variant `cov_vec` the
jackknife path builds, so genome-wide Stage 2 memory is unchanged (~10 GB at
1.2M variants). The pooled identity makes all of OLS, GLS, `tr(B)` and `Q`
factor into one `T x T` operation per chromosome plus one scalar per variant;
exploiting that would mean a second implementation of every downstream
quantity, and was deliberately left out of the first change.

---

## 11. The genotype-skewness dispersion correction: what it costs and what it refuses (2026-09-14)

`correct_dispersion_skew()` subtracts the exact bias `delta^2 * k`,
`k = 1 - 2*eaf`, that a pure location effect leaves in `beta_sigma2` at a skewed
genotype. On a Normal baseline it removes the artefact completely — the worst
cell in the sweep falls from 0.916 to 1.5e-09 at 5e-8. The full measurement is
`docs/DISPERSION_SKEW_CORRECTION.md`; the numbers that constrain *use* are here.

### What it costs

Nothing measurable on null variants, and nothing measurable on a genuine
dispersion effect with no location effect (power at 5e-8 unchanged in the fifth
decimal). The one real cost is at a rare variant carrying a colossal location
effect, where the correction has genuine variance of its own:

| MAF | Z_mu | power unc @5e-8 | power cor @5e-8 | honest power | SE ratio |
|---|---|---|---|---|---|
| 0.01 | 16 | 0.976 | 0.562 | 0.50 | 1.025 |
| 0.01 | 30 | 1.000 | **0.372** | 0.50 | 1.094 |

**About 26% of honest power is lost at MAF 0.01 / Z_mu 30.** That is the worst
cell in the whole sweep, it is the same cell where the *uncorrected* rate is
0.916, and `dispersion_skew_guard()` flags it with `G` = 6.3. `G` has therefore
inverted in meaning: it no longer bounds false positives, it marks where honest
power has been paid away.

### The skewness envelope is coupled to the estimator

Corrected rates at 5e-8, MAF 0.01, Z_mu 30, **GLS**: skewness 0.63 -> 3.8e-10;
0.71 -> 3.6e-09; 0.82 -> 5.0e-07; 0.89 -> 4.8e-05; 0.95 -> **1.00**. The usable
envelope under GLS is about **|skewness| 0.65**; under OLS it runs out to about
**1.0**. Choose the estimator and the limit together.

The failure is **not** summarised by skewness alone, nor by kurtosis alone. At
matched skewness 0.80 under GLS, skew-normal (excess kurtosis 0.65) leaves
2.0e-04 while chi2(12) (excess kurtosis 1.00) leaves 5.0e-07 — 390x apart, with
the *lower*-kurtosis distribution worse. Symmetric t(6), skewness exactly 0,
still leaves a 62x residual. What governs it is how far `f0'/f0` departs from
proportionality to `qnorm(tau)`, which no single moment captures. Transform the
phenotype (`docs/phenotype_preparation.md`) and none of this is live.

### Three ways to get a wrong answer, all now refused

Two of the three are **worse than not correcting at all**, which is why they
refuse rather than warn.

1. **A folded MAF in `eaf`.** `k = 1 - 2*eaf` is signed, so a folded frequency at
   a variant whose effect allele is the major one makes the correction *add* the
   bias: |Z| 3.27 -> 6.20, p 1.1e-03 -> 5.6e-10; at MAF 0.01 / Z_mu 30 the rate
   at 5e-8 goes to **1.000**, against 0.916 uncorrected. Stage 1 emits no
   frequency column, so a folded MAF is the likeliest input. Refused: an `eaf`
   below 0.5 at every variant, `eaf` and `maf` being the same values in
   `dispersion_report()`, and any call that has not set `eaf_orientation`.
2. **A `W` merely *named* `beta_sigma2`.** The proof is about the numbers.
   `check_vqtl_design()` verifies `W[, mu_name] == 1` and
   `W[, sigma2_name] == qnorm(taus)/(2*sd)` elementwise, which is why
   `dispersion_report()` requires `W_sd`. A shape check is not enough: the same
   Normal model parameterised in the SD rather than the variance is still
   proportional to `qnorm(taus)` and over-subtracts by `2*sd` — at sd = 1.3, from
   +1.51e-02 (truth 0) to -2.41e-02, 1.6x the original bias with the sign
   flipped. The 5-parameter mixture under `mu_name = "beta_mu1"` over-subtracts
   by 2.4x.
3. **A `w_sd_ratio` that does not match the `sd` that built `W`.** It scales the
   whole correction linearly, with no symptom. It is 1 for the usual
   `make_weights_vqtl(taus, q_tau, mu = mean(y), sd = sd(y))`.

### Two things nothing in the package sees

* **Dispersion inflation from the phenotype's shape at MAF 0.50.** On a chi2(4)
  baseline with Z_mu = 60 it is 0.163 at nominal 0.05 and 4.5e-06 at 5e-8. The
  correction is identically zero there because `k = 0`, and
  `dispersion_skew_guard()` is identically zero there because `gskew = 0`. Both
  tools read zero exactly where this mechanism is worst.
* **An asymmetric tau grid.** The `O(delta^3)` term is orthogonal to the
  dispersion column only on a grid symmetric about 0.5. On
  `linspace(0.02, 0.70, 45)` the order of the residual drops from 4.00 to 3.01
  and its relative size grows about 70x, sign depending on which tail is
  truncated. `check_vqtl_design()` warns; it does not refuse.

### Not done

The correction is applied *after* the projection, which is exact only for
`make_weights_vqtl`'s design. Correcting `b(tau)` before the projection is valid
for every `W` and costs the same O(1) per variant, but it makes the phenotype SD
a required input (it stops cancelling) and needs a declared location contrast for
a general `W` (`beta_mu` is exactly `delta` only because `W[, "beta_mu"] = 1`).
Both are prerequisites Stage 1 and Stage 2 do not currently supply, so the class
is closed by refusing unverified designs rather than by generalising. Reasoning
in `docs/DISPERSION_SKEW_CORRECTION.md` section 2.1.

## Version skew is undetectable below 0.2.0 (as of 2026-09-16)

Every `fungwasStage2` build before 2026-09-16 reports `Version: 0.1.0`,
including builds whose behaviour differs. Two installs five hours and four
merges apart were indistinguishable by version string, and the only symptom
was a refusal at run time.

The concrete failure: Stage 1 now writes `sigma0_meta_version: 2`, and a
`fungwasStage2` build predating that change caps at 1 and refuses the input
with "which this version of fungwasStage2 does not support". Current source
sets `SIGMA0_META_VERSION_MAX <- 2L` and reads both, so the fix is always to
reinstall -- but nothing in the version string tells you that you must.

This bites hardest when a long-lived analysis tree keeps its own R library
(for example a paper repository's `.local/Rlib`) and Stage 1 is later
regenerated from a newer package: Stage 2 then refuses its own Stage 1 output.

The version is now 0.4.0, matching the tag series; DESCRIPTION had been left at 0.1.0 while tags reached v0.3.0, so the string was no guide to behaviour at all. Bump it whenever the on-disk formats or the
estimator behaviour change, so a stale install is visible without running it.

## The LOCO workflow can destroy the structure being estimated (2026-09-16)

Measured on testosterone: a REGENIE step-1 ridge fit without the mixture
variable as a covariate learns the mixture instead of polygenic signal, and
subtracting it costs 8-10% in SE and 0.07-0.10 in effect correlation. See
[LOCO_MIXTURE_HAZARD.md](LOCO_MIXTURE_HAZARD.md) for the measurement, the blind
diagnostic that detects it without knowing membership, and which other traits in
this project are at risk.

Related: [EFFICIENCY_BOUND.md](EFFICIENCY_BOUND.md) shows the mixture
information cost itself is only ~2% at testosterone's separation, so the LOCO
offset was the larger term by far.

Both of those measurements use a REGENIE comparator run `--ignore-pred`, with no
LOCO of its own. [COMPARATOR_LEDGER.md](COMPARATOR_LEDGER.md) rebuilds the
comparator with a proper within-stratum step 1 and registers every
G-MoQ-vs-REGENIE estimate in this project, marking which are like-for-like.
