# What should `cov_source` default to?

**Keep the shipped default. `cov_source = "jackknife"` with `tau_cov_estimator =
"ols"` weakly dominates `pooled` + `"ols"` on calibration: it is
indistinguishable from pooled wherever pooled is right, and it is the only one
of the four configurations that holds nominal coverage when the dispersion
effect is real. Its one measured cost is 0.6-7% relative power at matched actual
Type-I. The answer is not regime-dependent.**

`docs/FINDINGS_pooled_dispersion.md` §8 recommends making pooled the scan
default. That recommendation is not supported by the measurement below, and the
evidence it rests on has a hole: its tables carry `jackknife GLS | jackknife OLS
| pooled GLS` and **no pooled-OLS column at all**, so the comparison that
actually decides the default was never run. It is run here.

---

## 0. What was measured, and how

All four configurations — `{pooled, jackknife} x {ols, gls}` — through the
**shipped** Stage 1 (`core.run_stage1`, `pooled_cov.compute_sigma0`) and the
**shipped** Stage 2 (`param_gwas`), at `N = 461,320` and the shipped default
`--blocks 1000` (`DEFAULT_N_BLOCKS`), over 52 cells x 15,000 replicates =
780,000 simulated scans. One replicate is one entire scan: its own phenotype,
its own `Sigma_0`, its own jackknife.

Everything below is marked **MEASURED** (counted from shipped-code output, with
Monte Carlo error) or **EXACT** (closed form / root-finder, no Monte Carlo error
at all) or **NOT RUN**. No estimator is reimplemented here: the pooled
`cov_data` is built by the package's own `.pooled_cov_data()`, and the jackknife
`cov_data` is the float32 row-major upper triangle Stage 1 emits, handed over in
exactly the layout and precision `.read_tau_cov_data()` would produce from a
production `.cov.gz`.

Degrees of freedom follow each source's own shipped convention — pooled with
`nu = n_samples, nu_source = "pooled_sigma0_n_samples"`, jackknife with
`nu = NULL, n_blocks = 1000` and per-SNP `n_eff = d_resid`, which is rung 2 of
`.param_reference()` (the per-SNP Satterthwaite law). Running the jackknife
against the pooled `nu` would have flattered it.

The model fitted is the package's own vQTL weight matrix (`make_weights_vqtl`),
parameters `beta_mu` and `beta_sigma2`. `W` is built at the population
quantiles, so the comparison isolates the covariance question from noise in `W`;
this matches the `gls_coverage` study's convention. Coverage is against each
configuration's **own** population estimand (OLS and GLS target different points
when `W` is misspecified), and uses the package's own `t_nu` reference.

**Pooled and jackknife share their replicates**, so every comparison between
them is paired and the Monte Carlo error on a *difference* is far below the
error on either level.

### Engine validation (EXACT vs the prior document)

The exact engine here reproduces `FINDINGS_pooled_dispersion.md`'s published
pooled coverage to four decimals:

| cell | that document | this engine (EXACT) |
|---|---|---|
| MAF 0.30, Z 16, dispersion, `beta_mu` | 0.9477 | 0.94777 |
| MAF 0.30, Z 16, dispersion, `beta_sigma2` | 0.9464 | 0.94647 |
| MAF 0.005, Z 60, dispersion, `beta_mu` | 0.6350 | 0.63523 |
| MAF 0.005, Z 60, dispersion, `beta_sigma2` | 0.4748 | 0.47502 |

The `OLS_VS_GLS.md` / `N_SCALING.md` control (pooled+OLS under a true scale
effect, MAF 0.30, `c = 0.05`) comes out **EXACT** at 0.9435 for `beta_sigma2`
against their Monte Carlo 0.9409 (N = 30,000) and 0.9413 (N = 461,320) — within
1.5 and 1.3 MC SE. The exact ratio is *independent of N*, which is why those two
numbers agree with each other: the pooled error is a working-independence bias,
not a finite-sample one.

---

## 1. The structural fact that settles most of it

With `tau_cov_estimator = "ols"` the point estimate is
`theta = (W'W)^-1 W' beta_hat`. **It does not touch the covariance at all.**
So pooled+OLS and jackknife+OLS produce *the same numbers* for every parameter;
they differ only in the standard error reported next to them.

**MEASURED**: the ratio of empirical sd of `theta_hat`, jackknife over pooled,
is `1.00000` in **all 44 main cells**, both parameters. There is no efficiency
trade-off in the point estimate, no estimand difference, and therefore no room
for a "scan with one, follow up with the other" answer on the OLS path. The only
question is which SE is closer to the truth.

---

## 2. Where the two are identical

**MEASURED**, 15,000 reps, MC SE in brackets. Under the null, and under a pure
location effect (the false-vQTL scenario), the two covariance sources are
indistinguishable:

| cell | param | pooled+OLS | jackknife+OLS | paired difference |
|---|---|---|---|---|
| null, MAF 0.30, coverage | `beta_sigma2` | 0.9502 (0.0018) | 0.9501 (0.0018) | -0.0001 |
| null, MAF 0.05, coverage | `beta_sigma2` | 0.9484 (0.0018) | 0.9484 (0.0018) | 0.0000 |
| null, MAF 0.01, coverage | `beta_sigma2` | 0.9459 (0.0018) | 0.9459 (0.0018) | -0.0001 |
| null, MAF 0.005, coverage | `beta_sigma2` | 0.9479 (0.0018) | 0.9481 (0.0018) | +0.0003 |
| null, MAF 0.30, Type-I @0.05 | `beta_sigma2` | 0.0498 (0.0018) | 0.0499 (0.0018) | +0.0001 |
| null, MAF 0.005, Type-I @0.05 | `beta_sigma2` | 0.0521 (0.0018) | 0.0519 (0.0018) | -0.0003 |
| location Z 16, MAF 0.30, coverage | `beta_sigma2` | 0.9522 (0.0017) | 0.9529 (0.0017) | +0.0007 |
| location Z 16, MAF 0.05, coverage | `beta_sigma2` | 0.9473 (0.0018) | 0.9485 (0.0018) | +0.0013 |

SE ratio (mean reported SE / empirical sd) at the null is 0.9916-1.0025 for
pooled+OLS and 0.9914-1.0023 for jackknife+OLS — the same to three decimals,
over MAF 0.005-0.30 and both phenotype shapes.

**So there is no calibration cost to the jackknife in the overwhelming majority
of a scan.**

---

## 3. Where they differ: a true dispersion effect

### 3.1 SE ratio and its direction — **MEASURED**

`mean reported SE / empirical sd`. Below 1 is **anti-conservative**.

| MAF | Z (dispersion) | param | pooled+OLS | jackknife+OLS |
|---|---|---|---|---|
| 0.30 | 5 | `beta_sigma2` | 0.9966 (0.0056) | 1.0007 (0.0058) |
| 0.30 | 16 | `beta_sigma2` | 0.9900 (0.0060) | 1.0040 (0.0061) |
| 0.30 | 60 | `beta_sigma2` | 0.9468 (0.0054) | 0.9964 (0.0058) |
| 0.05 | 16 | `beta_sigma2` | 0.9206 (0.0057) | 0.9860 (0.0055) |
| 0.05 | 60 | `beta_sigma2` | 0.7648 (0.0042) | 0.9851 (0.0058) |
| 0.01 | 16 | `beta_sigma2` | 0.8442 (0.0052) | 0.9953 (0.0060) |
| 0.01 | 30 | `beta_sigma2` | 0.7339 (0.0044) | 0.9961 (0.0060) |
| 0.005 | 10 | `beta_sigma2` | 0.8707 (0.0051) | 1.0087 (0.0059) |
| 0.005 | 16 | `beta_sigma2` | 0.7862 (0.0047) | 0.9940 (0.0057) |
| 0.005 | 30 | `beta_sigma2` | 0.6519 (0.0038) | 0.9955 (0.0056) |

Direction: **pooled+OLS is anti-conservative in every dispersion cell and never
conservative**, exactly as the prior document's concavity argument predicts.
**jackknife+OLS sits on 1.00** across the whole realistic range; over all 44 main
cells its SE ratio stays within 0.986-1.012 except in the two unphysical cells
of §6.

### 3.2 Coverage of 95% intervals — **MEASURED**

| MAF | Z | param | pooled+OLS | jackknife+OLS | pooled+GLS | jackknife+GLS |
|---|---|---|---|---|---|---|
| 0.30 | 16 | `beta_mu` | 0.9489 (0.0018) | 0.9515 (0.0018) | 0.9473 (0.0018) | 0.9413 (0.0019) |
| 0.30 | 16 | `beta_sigma2` | 0.9467 (0.0018) | 0.9507 (0.0018) | 0.9451 (0.0019) | 0.9379 (0.0020) |
| 0.30 | 60 | `beta_sigma2` | 0.9294 (0.0021) | 0.9421 (0.0019) | 0.9225 (0.0022) | 0.9321 (0.0021) |
| 0.05 | 16 | `beta_sigma2` | 0.9285 (0.0021) | 0.9447 (0.0019) | 0.9265 (0.0021) | 0.9347 (0.0020) |
| 0.05 | 60 | `beta_sigma2` | 0.8521 (0.0029) | 0.9367 (0.0020) | 0.8574 (0.0029) | 0.9058 (0.0024) |
| 0.01 | 16 | `beta_sigma2` | 0.9027 (0.0024) | 0.9491 (0.0018) | 0.8980 (0.0025) | 0.9380 (0.0020) |
| 0.01 | 30 | `beta_sigma2` | 0.8475 (0.0029) | 0.9452 (0.0019) | 0.8393 (0.0030) | 0.9169 (0.0023) |
| 0.005 | 10 | `beta_sigma2` | 0.9115 (0.0023) | 0.9511 (0.0018) | 0.9057 (0.0024) | 0.9399 (0.0019) |
| 0.005 | 16 | `beta_sigma2` | 0.8742 (0.0027) | 0.9488 (0.0018) | 0.8661 (0.0028) | 0.9360 (0.0020) |
| 0.005 | 30 | `beta_sigma2` | 0.7919 (0.0033) | 0.9437 (0.0019) | 0.7856 (0.0034) | 0.6222 (0.0040) |

The paired jackknife-minus-pooled coverage difference is **positive in every
dispersion and every mixed cell**, from +0.0012 (MAF 0.30, Z 5) to +0.2104
(MAF 0.005, Z 60).

### 3.3 Type-I error on the location parameter under a true dispersion effect

`beta_mu` is truly zero here, so every rejection is false. This is the exact
quantity `FINDINGS_pooled_dispersion.md` §4 makes its jackknife claim about.

**MEASURED at nominal 0.05:**

| MAF | Z | pooled+OLS | jackknife+OLS | pooled+GLS | jackknife+GLS |
|---|---|---|---|---|---|
| 0.30 | 16 | 0.0511 (0.0018) | 0.0485 (0.0018) | 0.0527 (0.0018) | 0.0587 (0.0019) |
| 0.05 | 16 | 0.0617 (0.0020) | 0.0489 (0.0018) | 0.0611 (0.0020) | 0.0602 (0.0019) |
| 0.05 | 60 | 0.1022 (0.0025) | 0.0490 (0.0018) | 0.1015 (0.0025) | 0.1003 (0.0025) |
| 0.01 | 16 | 0.0769 (0.0022) | 0.0467 (0.0017) | 0.0775 (0.0022) | 0.0646 (0.0020) |
| 0.01 | 30 | 0.1183 (0.0026) | 0.0498 (0.0018) | 0.1219 (0.0027) | 0.0875 (0.0023) |
| 0.005 | 16 | 0.0959 (0.0024) | 0.0495 (0.0018) | 0.0992 (0.0024) | 0.0731 (0.0021) |
| 0.005 | 30 | 0.1632 (0.0030) | 0.0497 (0.0018) | 0.1677 (0.0031) | 0.1297 (0.0027) |
| 0.005 | 60 | 0.3920 (0.0040) | 0.0656 (0.0020) | 0.3323 (0.0038) | 0.6958 (0.0038) |

**Extrapolated to 5e-8**, as inflation over nominal. The extrapolation scales the
Wald statistic against the package's own `t_nu` reference; the scale is estimated
from `sd(T)`. (A second estimate from the 99.9th percentile of `|T|` agrees to
within a factor of ~1.5 in the moderate cells; see `out/tables.md` §9.)

| MAF | Z | pooled+OLS | **jackknife+OLS** | pooled+GLS | jackknife+GLS |
|---|---|---|---|---|---|
| 0.30 | 16 | 1.26x | **0.94x** | 1.29x | 3.45x |
| 0.30 | 60 | 3.61x | **1.22x** | 3.65x | 4.95x |
| 0.05 | 16 | 4.38x | **1.15x** | 4.72x | 4.81x |
| 0.05 | 60 | 98x | **0.71x** | 96x | 85x |
| 0.01 | 16 | 18.2x | **0.74x** | 20.9x | 5.18x |
| 0.01 | 30 | 295x | **1.08x** | 351x | 40x |
| 0.005 | 16 | 82x | **1.04x** | 101x | 12.7x |
| 0.005 | 30 | 2086x | **0.99x** | 2525x | 457x |

---

## 4. The prior document's jackknife claim is right — but only for OLS

> "the location parameter's false-positive rate at 5e-8 ... is **5.00e-08 —
> exactly nominal — for the jackknife in every cell** ... So the jackknife is a
> genuine and complete fix for the dispersion-effect covariance failure"
> — `FINDINGS_pooled_dispersion.md` §4

**That sentence never says which estimator, and it is only true of one of them.**

* With **OLS** it holds, and it now has measured backing rather than an exact
  column: jackknife+OLS lands at 0.7x-1.3x nominal at 5e-8 in every cell of the
  realistic range, and its nominal-0.05 rate is 0.0467-0.0505 against a target
  of 0.05.
* With **GLS** it is false. jackknife+GLS is inflated **3.45x-4.95x at 5e-8 even
  at MAF 0.30**, 12.7x at MAF 0.005/Z 16, 457x at MAF 0.005/Z 30, and at
  MAF 0.005/Z 60 it is **worse than pooled** (0.6958 vs 0.3920 at nominal 0.05).

The cause is the one that document identified but did not carry through to its
headline claim: a 45x45 covariance estimated from 1000 blocks is noisy enough
for GLS to overfit it. **MEASURED** directly here, the relative noise in the
reported SE:

| source | measured `sd(SE)/mean(SE)` | theory `1/sqrt(2(B-1))` |
|---|---|---|
| jackknife, B = 1000 | 0.0224 - 0.0311 | 0.0224 |
| pooled | 0.0024 - 0.0077 | n/a (one matrix, N individuals) |

Its own §4 table is reproduced and completed. **MEASURED** `mean_se/emp_sd` at
the null, over MAF 0.005-0.30 and both phenotype shapes:

| configuration | range over null cells | null Type-I @0.05 |
|---|---|---|
| pooled+OLS | 0.9916 - 1.0025 | 0.0460 - 0.0541 |
| **jackknife+OLS** | **0.9914 - 1.0023** | **0.0460 - 0.0541** |
| pooled+GLS | 0.9919 - 1.0043 | 0.0487 - 0.0541 |
| jackknife+GLS | 0.9374 - 0.9581 | 0.0605 - 0.0680 |

That document's `jackknife GLS 0.898-0.954` column is confirmed. **The
pooled-OLS column it was missing reads 0.99-1.00 at the null and 0.65-0.99 under
dispersion; the jackknife-OLS column beside it reads 0.99-1.01 throughout.**

---

## 5. What the jackknife does *not* fix, confirmed

`FINDINGS_pooled_dispersion.md` §4's other half — the jackknife is "no help
whatsoever for the false-vQTL problem" — is **confirmed exactly**, and for the
reason given there: the artefact is in the estimand, and (§1) OLS's point
estimate is identical under the two covariance sources.

**MEASURED**, `beta_sigma2` rejection rate at nominal 0.05 under a **pure
location** effect with zero true dispersion:

| MAF | Z_mu | pooled+OLS | jackknife+OLS |
|---|---|---|---|
| 0.30 | 16 | 0.0526 (0.0018) | 0.0520 (0.0018) |
| 0.30 | 30 | 0.0838 (0.0023) | 0.0838 (0.0023) |
| 0.30 | 60 | 0.6163 (0.0040) | 0.6112 (0.0040) |
| 0.05 | 16 | 0.1195 (0.0026) | 0.1177 (0.0026) |
| 0.01 | 16 | 0.4507 (0.0041) | 0.4363 (0.0040) |
| 0.005 | 16 | 0.7523 (0.0035) | 0.7300 (0.0036) |
| 0.005 | 30 | 1.0000 (0.0000) | 1.0000 (0.0000) |

The jackknife buys at most 1-3% relative. The `G` guard of that document's §5
remains the right instrument for this failure, and changing `cov_source` is not
an alternative to it.

---

## 6. Power: the jackknife's one real cost

Power to detect a genuine dispersion effect, at **matched ACTUAL Type-I** — the
threshold is the empirical alpha-quantile of `|theta/SE|` from the matching null
cell, so each configuration is held to the same true error rate rather than the
same nominal one. **MEASURED**; Z above about 4 saturates at power 1, so this is
the only informative range.

Actual Type-I = 0.05:

| MAF | Z | pooled+OLS | jackknife+OLS | jackknife loss (relative) |
|---|---|---|---|---|
| 0.30 | 2 | 0.5182 (0.0041) | 0.5151 (0.0041) | 0.6% |
| 0.30 | 3 | 0.8516 (0.0029) | 0.8502 (0.0029) | 0.2% |
| 0.05 | 2 | 0.5237 (0.0041) | 0.5175 (0.0041) | 1.2% |
| 0.01 | 2 | 0.5231 (0.0041) | 0.5082 (0.0041) | 2.8% |
| 0.005 | 2 | 0.5305 (0.0041) | 0.5051 (0.0041) | 4.8% |
| 0.005 | 3 | 0.8642 (0.0028) | 0.8509 (0.0029) | 1.5% |

Actual Type-I = 0.01:

| MAF | Z | pooled+OLS | jackknife+OLS | jackknife loss (relative) |
|---|---|---|---|---|
| 0.30 | 2 | 0.2900 (0.0037) | 0.2840 (0.0037) | 2.1% |
| 0.30 | 3 | 0.6733 (0.0038) | 0.6641 (0.0039) | 1.4% |
| 0.01 | 3 | 0.6902 (0.0038) | 0.6637 (0.0039) | 3.8% |
| 0.005 | 2 | 0.3077 (0.0038) | 0.2635 (0.0036) | 14.4% |
| 0.005 | 3 | 0.7083 (0.0037) | 0.6594 (0.0039) | 6.9% |

This is the whole price of the jackknife, and it is exactly the price the block
count buys down: a noisier denominator, `sd(SE)/SE` of 2.2-3.1% against 0.2-0.8%,
costs 0.2-7% relative power at common-to-rare MAF, rising to 14% in the
hardest corner (MAF 0.005, Z 2, alpha 0.01). Set against pooled's 82x-2086x
inflation of the 5e-8 false-positive rate and 0.79 coverage in the same MAF
band when the dispersion effect is real, it is not a close trade. Note also that
matched-actual-Type-I *flatters* pooled: in a real scan you cannot match actual
error, you only get nominal, and at nominal thresholds pooled's extra rejections
are mostly error rather than power.

---

## 7. Where each configuration stops being trustworthy

**Realism first.** At N = 461,320 all four MAF values used here are inside the
production range — the chr22 release floor is MAC ~927, and these cells carry
MAC 276,792 (MAF 0.30), 46,132 (0.05), 9,226 (0.01) and 4,613 (0.005). **Every
MAF cell here could occur.** What cannot occur is the top of the effect ladder.
The SD ratio between the two homozygotes required to reach each Z (**EXACT**):

| MAF | Z = 5 | Z = 16 | Z = 30 | Z = 60 |
|---|---|---|---|---|
| 0.30 | 1.02x | 1.06x | 1.11x | 1.23x |
| 0.05 | 1.04x | 1.12x | 1.24x | 1.59x |
| 0.01 | 1.08x | 1.29x | 1.67x | **3.50x** |
| 0.005 | 1.12x | 1.45x | **2.17x** | **12.0x** |

**The `MAF 0.005 / Z 60` cell that the prior document leads with — the
0.6350 / 0.4748 collapse — requires the minor homozygote to carry twelve times
the phenotypic standard deviation of the major homozygote.** It is a fair stress
test of the estimator and a real property of it, but it is not a cell a human
GWAS will produce, and no recommendation should be anchored on it. The cells
that matter are MAF 0.01-0.05 at Z 10-30, SD ratios 1.1x-1.7x — and there pooled
is already 18x-295x inflated at 5e-8 while jackknife+OLS is at 0.7x-1.2x.

| configuration | trustworthy up to | fails beyond |
|---|---|---|
| **jackknife+OLS** | everything in the realistic grid: all four MAF, Z <= 30, both phenotypes, all four truths | SD ratio >= 3.5x (MAF 0.01 Z 60: SE ratio 0.917; MAF 0.005 Z 60: 0.829) |
| pooled+OLS | MAF 0.30 at any Z (5e-8 inflation <= 3.6x); MAF 0.05 to about Z 16 (4.4x) | MAF <= 0.01 with any real dispersion effect: 18x at 5e-8 at Z 16, 295x at Z 30 |
| pooled+GLS | tracks pooled+OLS, slightly worse everywhere | same, plus the GLS skew failure below |
| jackknife+GLS | nothing. 3.45x at 5e-8 even at MAF 0.30 under a mild dispersion effect, and 21-36% inflated at the null | everywhere |

### The skewed phenotype is a GLS problem, not a `cov_source` problem

**MEASURED**, standardised gamma(2) phenotype (skew 1.41), coverage of
`beta_sigma2`:

| cell | pooled+OLS | jackknife+OLS | pooled+GLS | jackknife+GLS |
|---|---|---|---|---|
| null, MAF 0.30 | 0.9527 (0.0017) | 0.9526 (0.0017) | 0.9502 (0.0018) | 0.9395 (0.0019) |
| dispersion Z 16, MAF 0.30 | 0.9450 (0.0019) | 0.9479 (0.0018) | **0.8189 (0.0031)** | **0.8066 (0.0032)** |
| dispersion Z 16, MAF 0.01 | 0.8988 (0.0025) | 0.9472 (0.0018) | **0.6859 (0.0038)** | **0.8068 (0.0032)** |
| location Z 16, MAF 0.01 | 0.9305 (0.0021) | 0.9471 (0.0018) | 0.9537 (0.0017) | **0.8061 (0.0032)** |

and the false-vQTL rate at nominal 0.05 for `beta_sigma2` under a pure location
effect at MAF 0.01 on that phenotype is **1.0000 for both GLS configurations**
against 0.3688 (pooled+OLS) and 0.3065 (jackknife+OLS). On a skewed outcome GLS
rejects a false dispersion hypothesis every single time, with either covariance.

---

## 8. Does the answer depend on the estimator?

Yes, but not in the shape the question is usually asked. `cov_source` and
`tau_cov_estimator` are not two independent knobs to tune per regime:

* **On the OLS path**, `cov_source` changes only the reported SE. Jackknife is
  right where pooled is right and right where pooled is wrong. Pooled's only
  compensating advantage is a steadier denominator, worth 0.2-7% relative power.
* **On the GLS path**, `cov_source` changes the point estimate too, and both
  choices are bad: pooled+GLS inherits pooled's working-independence error
  essentially undiminished (its 5e-8 inflation tracks pooled+OLS to within 15%
  in every cell), and jackknife+GLS adds an overfitting error large enough to be
  anti-conservative at the null, before any dispersion effect exists.

So the regime-dependent answer is the wrong shape. **One configuration
dominates: `jackknife` + `ols`.** The shipped defaults already are that pair,
and they should not be changed.

---

## 9. Recommendations

1. **Leave `cov_source = c("jackknife", "pooled")` as it stands.** Do not
   promote pooled to the default. There is no cell in this grid, realistic or
   otherwise, where pooled+OLS beats jackknife+OLS on calibration.
2. **Withdraw recommendation 1 of `FINDINGS_pooled_dispersion.md` §8** ("Keep
   pooled as the scan default"). It was reached without the pooled-OLS column
   and it is contradicted here. What survives from that document is its
   diagnosis: the false-vQTL artefact is an estimand gap, the `G` guard is the
   right instrument for it, and no covariance estimator fixes it.
3. **Amend §4 of that document to say *jackknife with OLS*.** As written, "the
   jackknife ... is 5.00e-08 — exactly nominal — in every cell" reads as a
   property of the covariance source. It is a property of the pair, and it is
   false for the other pairing by up to 457x.
4. **Leave `tau_cov_estimator = c("ols", "gls")` as it stands too.** The GLS
   path is anti-conservative at the null with the jackknife and collapses on a
   skewed phenotype with either covariance. It deserves a warning in the
   documentation, not a promotion.
5. **If the pooled path is kept for cost reasons** — it is what lets a scan skip
   the 229-289 bytes/SNP covariance payload — document that its SEs are only
   trustworthy where there is no true dispersion effect, which is precisely the
   hypothesis a vQTL scan is testing. Any variant the scan flags as a dispersion
   signal should have its inference recomputed against the jackknife before an
   interval or a p-value is reported. That is a defensible two-stage design; it
   is the opposite of the one §8 proposed.

---

## 10. What this does not cover

* **NOT RUN: the production `B = 65` jackknife**, or the int8-quantised
  covariance. Only `B = 1000` (the shipped `DEFAULT_N_BLOCKS`) and float32 (the
  shipped on-disk precision). `B = 65` quadruples the SE noise (8.8% against
  2.2%), so the power cost in §6 would grow by roughly that factor and the GLS
  path would get worse; the ordering under a true dispersion effect at rare MAF
  would not change, because pooled's error there is a bias and not noise.
* **NOT RUN: real genotypes, LD, covariates with real effects, relatedness.**
  Everything uses a 3-point HWE genotype law and an intercept-only covariate
  basis.
* **NOT RUN: `N` other than 461,320.** The pooled SE ratio is exactly
  N-invariant (§0), so that surface carries over; the absolute Type-I rates do
  not, because `Z` at fixed effect size scales with `sqrt(N)`.
* **NOT RUN: baselines other than normal and standardised gamma(2).** Heavier
  tails were not tried, and the GLS skew failure in §7 suggests that is where to
  look next.
* **The 5e-8 numbers are extrapolations, not counts.** 15,000 replicates reach
  1e-3, where the measured null rates are 0.00080-0.00140 against a nominal
  0.001 for every configuration except jackknife+GLS (0.00127-0.00233). Below
  that, rates come from scaling the Wald statistic against the package's own
  `t_nu` reference, which assumes the statistic stays in that family into the
  far tail. Two independent scale estimates are reported in `out/tables.md` §9;
  where they disagree, trust the ordering rather than the magnitude.
* **Two cells where the exact engine and the shipped pipeline disagree.** At
  MAF 0.005 with Z 30 and Z 60 the measured pooled+OLS coverage is *higher* than
  the asymptotic exact value (0.7919 vs 0.7498, and 0.5486 vs 0.4750). The
  asymptotic calculation is pessimistic there because `beta_hat` at MAC 4,613
  with a 2.2x-12x SD ratio has not reached its asymptotic law. Elsewhere the two
  agree to within 2.6 MC SE. Worth recording, because
  `FINDINGS_pooled_dispersion.md` quotes the exact 0.4748 without that caveat.

---

## 11. Reproducing

Scripts in `sim/covsource/` (see its `README.md`); committed output in
`sim/covsource/tables/`. The eighteen full tables — coverage, SE ratio, Type-I,
5e-8 extrapolation, power, paired differences, MC-vs-exact, deep 99.9% coverage
— are in `sim/covsource/tables/tables.md`, rendered from `summary_cs.tsv`
(624 rows), `power_cs.tsv`, `paired_cs.tsv` and `exact_table.tsv` beside it.
`cells.json` is the calibrated grid.
