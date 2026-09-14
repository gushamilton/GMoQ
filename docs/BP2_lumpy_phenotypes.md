# BP2 — Does G-MoQ hold up on phenotypes that are not smooth and continuous?

**Question.** Stage 1 builds `RIF_i(tau) = q_tau + (tau - 1{Y_i <= q_tau}) / f_hat(q_tau)` with a Gaussian
KDE at Silverman bandwidth `h = 1.06 sd(Y) N^(-1/5)` (`stage1/fungwas_stage1/core.py:308`). Every effect
and SE divides by `f_hat(q_tau)`. `BP1_density_estimate.md` showed the bandwidth is 3.2–3.4x too wide for
bimodal traits and distorts `f_hat` by up to 32%, but nobody translated that into coverage or Type-I
error. This brief does, across five families of non-smooth phenotype: ties/rounding, floor and ceiling
censoring, bimodality, count-like traits, and heavy skew.

Written 2026-09-14. All compute under `/user/work/fh6520/kde_lumpy/`. Stage 1 and Stage 2 are the
**shipped** code, unmodified (`tree/stage1` is byte-identical to `stage1/fungwas_stage1/core.py` on this
branch; `fungwasStage2` installed into a private lib with `useDynLib` verified intact). Simulated
genotypes at MAF 0.05/0.20/0.50; no real data read.

---

## Bottom line

**The failure is not where the density doc pointed.** The KDE error is real and large — up to 2.0x at a
bimodal median trough, up to 38x below a censoring limit — but it **cancels exactly** in every
single-quantile test statistic, and it never inflates Type-I error at Stage 1. What breaks is Stage 2
**GLS**, and it breaks for a different reason: when a block of tied values spans two requested taus,
those taus share the same empirical quantile, their RIF columns become identical, and the tau-covariance
becomes numerically singular (condition number 1e2 → 1e13). Ridge repair then fires on 100% of SNPs and
inflates the SEs by 10^2–10^3x.

| shape | verdict |
|---|---|
| Smooth, symmetric (normal) | **Safe.** The reference case. |
| Bimodal, no ties | **Safe for calibration, mild bias.** Null coverage 0.9438; causal coverage of the true scale parameter 0.920 (jackknife GLS). |
| Rounded, ties below the tau spacing (≤ ~4% here) | **Safe.** Null coverage 0.9433, Type-I 0.0567 — indistinguishable from the smooth control. |
| Rounded, ties above the tau spacing | **Not safe in GLS.** Coverage → 0.997–1.000, Type-I → 0.000. Total power loss. |
| Floor / ceiling censored (≥5% at the limit) | **Not safe.** Broken at 5% censoring. Ceiling censoring is additionally **anti-conservative in the far tail** (see below). |
| Count-like (≤8 distinct values) | **Not safe.** |
| Heavily skewed but smooth (lognormal) | **Calibrated but biased.** Null coverage 0.9455, but causal coverage of the latent truth falls to 0.277 at tau = 0.10. No tie guard can see this. |

**Where the boundary sits.** Not at a tie fraction. A tau is safe iff its tie fraction is **smaller than
the gap to its nearest neighbouring tau**:

| cell | max tie | max tie/spacing | mean duplicated q_tau | GLS null cov95 (scale) |
|---|---|---|---|---|
| round_0p10 | 0.039 | **0.53** | 0.00 | 0.9433 (0.0005) |
| round_0p25 | 0.097 | **1.46** | 0.19 | 0.9524 (0.0005) |
| floor_05 | 0.054 | **1.80** | 1.00 | 0.9997 (0.0000) |
| round_0p50 | 0.193 | **2.42** | 5.99 | 0.9973 (0.0001) |
| floor_20 | 0.205 | **6.85** | 3.98 | 1.0000 (0.0000) |

`round_0p10` has ties at *every* tau and is fine. `floor_05` has ties at only two taus, 5.4% of the
sample, and is broken. The discriminator is collision, not tie mass. With the 15-tau grid used here
(minimum spacing 0.03) the practical threshold is **~3% ties at a closely spaced tau, ~10% at a widely
spaced one**.

---

## 1. The density error, and why it does not do what you would expect

Population truth, computed by quadrature (`src/shapes.py`), `f_ship` = the Gaussian-KDE convolution the
shipped estimator actually converges to, against the true density of the observed trait:

**`f_ship / f_true`** (`nan` = an atom sits at `q_tau`, so no density exists there)

| cell | .02 | .05 | .10 | .20 | .50 | .80 | .95 | .98 |
|---|---|---|---|---|---|---|---|---|
| normal_smooth | 1.034 | 1.018 | 1.007 | 0.997 | 0.989 | 0.997 | 1.018 | 1.034 |
| **bimodal_smooth** | 1.189 | 1.051 | 0.968 | 0.914 | **2.013** | 0.914 | 1.051 | 1.189 |
| **lognorm_smooth** | 1.018 | 0.789 | **0.754** | 0.824 | 1.013 | 1.030 | 1.015 | 1.010 |
| floor_20 | nan | nan | nan | nan | 0.993 | 0.998 | 1.012 | 1.024 |

The bimodal median trough is estimated at **2.0x** the truth — six times worse than the 32% the BP1 doc
measured on real testosterone, because this synthetic bimodal has a deeper trough. The lognormal shoulder
is **25% low** at tau = 0.10 with no ties involved at all.

Against the *latent* density a user believes they are working with, the errors are far larger:
`f_ship/f_latent` reaches **38.3** at tau = 0.02 for `floor_40` and **19.9** for `count_4`.

**But this cancels.** Perturbing the bandwidth by 0.5x and 2.0x and re-running the shipped Stage 1
(`src/invariance.py`):

```
normal_smooth  scale:0.5  max|b/(b0*sc)-1|=3.28e-13  max|se/(se0*sc)-1|=2.58e-14  max|z-z0|=5.11e-14
round_0p50     scale:0.5  max|b/(b0*sc)-1|=6.76e-13  max|se/(se0*sc)-1|=9.60e-14  max|z-z0|=2.32e-13
floor_20       scale:2.0  max|b/(b0*sc)-1|=1.56e-13  max|se/(se0*sc)-1|=7.93e-14  max|z-z0|=1.45e-13
```

`b` and `se` each scale exactly as `1/f_hat`, so `z = b/se` is invariant to **5e-14**. Whatever the
density estimate does, the single-quantile Wald test does not notice.

And indeed it does not. Across all 15 runA cells x 15 taus (225 combinations, 200,000 null SNPs each),
**Stage 1 Type-I at 0.05 has median 0.0501, 95th percentile 0.0521, maximum 0.0575**. Not one tau
anywhere is anti-conservative. The only departures are *downward*, at quantiles that are degenerate
because the trait has no variation there (`ceil_20` above the ceiling: 0.000; `count_4` at tau 0.95/0.98:
0.000).

The per-tau jackknife SE ratio is likewise flat at 0.936–0.942 for every shape — that ~6% conservatism
is the block jackknife's own property at B = 1000, present identically in the smooth control.

---

## 2. Where it does break: Stage 2 GLS and the singular tau-covariance

Stage 2, null SNPs, n = 200,000 per cell, `tau_cov_estimator` both ways. Coverage at 95% and Type-I at
0.05 for the scale parameter (MC error in brackets):

| cell | jk_ols cov95 / T1 | **jk_gls cov95 / T1** | pl_ols cov95 / T1 | **pl_gls cov95 / T1** |
|---|---|---|---|---|
| normal_smooth | 0.9498 / 0.0502 | 0.9437 / 0.0563 | 0.9504 / 0.0496 | 0.9503 / 0.0497 |
| lognorm_smooth | 0.9486 / 0.0514 | 0.9455 / 0.0545 | 0.9499 / 0.0501 | 0.9507 / 0.0493 |
| bimodal_smooth | 0.9504 / 0.0496 | 0.9438 / 0.0562 | 0.9502 / 0.0498 | 0.9504 / 0.0496 |
| round_0p10 | 0.9495 / 0.0505 | 0.9433 / 0.0567 | 0.9498 / 0.0502 | 0.9500 / 0.0500 |
| round_0p25 | 0.9498 / 0.0502 | 0.9524 / 0.0460 | 0.9499 / 0.0501 | 0.9586 / 0.0415 |
| round_0p50 | 0.9502 / 0.0498 | **0.9973 / 0.0000** | 0.9503 / 0.0497 | **1.0000 / 0.0000** |
| round_1p00 | 0.9495 / 0.0505 | **1.0000 / 0.0000** | 0.9496 / 0.0504 | **1.0000 / 0.0000** |
| floor_05 | 0.9500 / 0.0500 | **0.9997 / 0.0003** | 0.9505 / 0.0495 | **0.9999 / 0.0001** |
| floor_20 | 0.9501 / 0.0499 | **1.0000 / 0.0000** | 0.9508 / 0.0492 | **1.0000 / 0.0000** |
| floor_40 | 0.9496 / 0.0504 | **1.0000 / 0.0000** | 0.9506 / 0.0494 | **1.0000 / 0.0000** |
| ceil_20 | 0.9495 / 0.0505 | **0.9971 / 0.0029** | 0.9499 / 0.0501 | **0.9624 / 0.0376** |
| count_8 | 0.9500 / 0.0500 | **1.0000 / 0.0000** | 0.9499 / 0.0501 | **1.0000 / 0.0000** |
| count_4 | 0.9499 / 0.0501 | **1.0000 / 0.0000** | 0.9496 / 0.0504 | **1.0000 / 0.0000** |
| bimodal_round | 0.9491 / 0.0509 | **0.9994 / 0.0000** | 0.9493 / 0.0507 | **1.0000 / 0.0000** |

MC error is 0.0005 on every coverage and 0.0005 on every Type-I. **OLS is immune everywhere** — 0.9491 to
0.9508 across every shape tested, including a 4-valued count trait.

The mechanism is visible directly in the covariance conditioning (median over 200 SNPs):

| cell | cond(Sigma), full taus | duplicate correlation pairs | cond after tie-drop |
|---|---|---|---|
| normal_smooth | 97.5 | 0.0 | — |
| round_0p10 | 105 | 0.0 | 6.25 |
| round_0p50 | **8.8e+12** | 6.0 | — |
| floor_05 | **1.2e+13** | 1.0 | **90.3** |
| floor_20 | **1.1e+13** | 10.0 | **86.6** |
| bimodal_round | **6.7e+13** | 14.0 | — |

When two taus share a sample quantile their RIF columns are literally the same vector, the tau-covariance
is rank-deficient, and the GLS solve falls back to ridge repair on 100% of SNPs. The Stage 2 log says so
explicitly: `ridge repair attempted 2020 (100.0000%), rescued 2020`.

**Why OLS escapes, and why that matters.** The OLS path is a *fixed-contrast projection*
(`stage2/R/tail-reference.R:433`): the contrast `(W'W)^-1 W'` does not depend on the estimated
covariance, and `Sigma` enters only through the sandwich `W' Sigma W`. It is never inverted, so a
rank-deficient `Sigma` costs nothing. GLS uses `Sigma^-1` as the weight matrix, which is precisely the
operation that fails. **`tau_cov_estimator = "gls"` is the package default**
(`stage2/R/gwas.R:223`), so the broken path is the one a user gets without asking.

**This is a power failure, not a false-positive failure — with one exception.** Rejection in the far tail
as a multiple of nominal (scale parameter, n = 200,000):

| cell | estimator | a = 1e-3 | a = 1e-4 | a = 1e-6 |
|---|---|---|---|---|
| normal_smooth | jackknife_gls | 1.53 (0.09) | 1.70 (0.29) | 5.00 (5.00) |
| round_0p50 | jackknife_gls | 0.00 | 0.00 | 0.00 |
| floor_20 | jackknife_gls | 0.00 | 0.00 | 0.00 |
| **ceil_20** | **jackknife_gls** | 0.80 (0.06) | **4.90 (0.49)** | **270 (37)** |
| **ceil_20** | **pooled_gls** | 1.42 (0.08) | **4.60 (0.48)** | **220 (33)** |

`ceil_20` is conservative in the body (Type-I 0.0029 at 0.05) but **270x nominal at 1e-6**. Ridge repair
gives a handful of SNPs catastrophically small SEs. **Ceiling-censored traits can manufacture
genome-wide-significant false positives in GLS.** This is the only anti-conservative failure found and it
is the one that actually threatens a published result.

---

## 3. Pooled covariance vs the block jackknife — the errors do not cancel

The question was whether a bad density, which scales both `b` and `Sigma_0`, partly cancels in the
z-statistic. It does not, and pooled is the worse of the two.

| cell | NULL cov95 jk / pooled | CAUSAL cov95 (scale) jk / pooled | CAUSAL cov95 (shift) jk / pooled |
|---|---|---|---|
| normal_smooth | 0.9437 / 0.9503 | 0.9445 / 0.9355 | 0.9475 / 0.9435 |
| lognorm_smooth | 0.9455 / 0.9507 | **0.8400 / 0.7450** | **0.7940 / 0.7100** |
| bimodal_smooth | 0.9438 / 0.9504 | **0.9200 / 0.8570** | 0.9440 / 0.9340 |
| bimodal_uneven | 0.9439 / 0.9494 | **0.9285 / 0.8820** | 0.9380 / 0.9210 |
| bimodal_round | 0.9994 / 1.0000 | 0.0000 / 0.0000 | **0.1380 / 0.0190** |

Under the null pooled is slightly *better* (it lacks the jackknife's 6% SE conservatism). Under the
alternative pooled is uniformly **worse**, by 4–9 coverage points on the bimodal and skewed traits. The
reason is structural: the pooled `Sigma_0` is an exact Gram matrix of the residualised RIF, so duplicated
quantiles make it *exactly* rank-deficient, whereas the jackknife covariance is singular only up to Monte
Carlo noise and retains a little usable signal. **No mitigating cancellation exists.**

---

## 4. The estimand shift — the problem no guard can see

Stage 1 is *consistent for what it estimates*. Coverage of `b_ship`, the population KDE-convolution
target, is 0.94–0.96 at essentially every cell and tau. Coverage of `b_latent`, the latent quantile effect
the user believes they are estimating, is not:

**Stage 1 per-tau coverage of the latent truth (95% nominal)**

| cell | .02 | .05 | .10 | .20 | .50 | .80 | .95 | .98 |
|---|---|---|---|---|---|---|---|---|
| normal_smooth | 0.942 | 0.945 | 0.942 | 0.951 | 0.952 | 0.954 | 0.951 | 0.950 |
| **lognorm_smooth** | 0.949 | **0.288** | **0.277** | 0.751 | 0.947 | 0.951 | 0.947 | 0.950 |
| **bimodal_smooth** | **0.771** | 0.925 | 0.941 | 0.905 | 0.954 | 0.903 | 0.927 | **0.760** |
| floor_20 | **0.000** | **0.000** | **0.000** | **0.000** | 0.952 | 0.955 | 0.950 | 0.950 |
| ceil_20 | 0.943 | 0.946 | 0.942 | 0.952 | 0.952 | **0.000** | **0.000** | **0.000** |
| round_1p00 | 0.216 | 0.718 | 0.000 | 0.067 | 0.625 | 0.000 | 0.000 | 0.000 |

**The failure is exactly local**, as predicted: `floor_20` is destroyed at taus 0.02–0.20, which are the
taus inside the censored block, and perfect (0.947–0.955) everywhere above it. `ceil_20` is the mirror
image. The corresponding estimand ratios `b_ship/b_latent` are 0.129/0.166/0.216/0.270/0.335 at
taus 0.02–0.20 for `floor_20` — a 3–8x attenuation — and exactly **0.000** for `ceil_20` above the
ceiling, where `F_k(q) = 1` for every genotype and there is no signal left to find.

**`lognorm_smooth` is the uncomfortable result.** It is perfectly smooth, has zero ties, passes every
diagnostic, is well calibrated under the null (0.9455), and has latent coverage of **0.277** at tau = 0.10
and Stage 2 causal coverage of **0.840 (scale) / 0.794 (shift)**, bias −0.91 and −1.10 empirical sd. This
is pure KDE bandwidth bias at a sharp density shoulder. It extends the BP1 doc's finding: the bandwidth
problem is not confined to bimodality, and it is not detectable from any tie-based guard.

---

## 5. Diagnostics and fixes

### 5.1 Q does not warn you

`Q` lack-of-fit rejection at 0.05 is **0.0000 in every one of the 15 cells**, including those where the
model is catastrophically misspecified. Q is useless as a warning here.

### 5.2 Bandwidth rules do not fix it

Re-running Stage 1 with three alternative rules (`src/bwpatch.py`), 500 reps each:

| cell | bw | cov95 (scale) | Type-I |
|---|---|---|---|
| round_0p50 | Silverman | 0.9973 | 0.0000 |
| round_0p50 | robust `0.9 min(sd, IQR/1.349)` | 1.0000 | 0.0000 |
| round_0p50 | ISJ (Sheather–Jones) | 0.9970 | 0.0000 |
| round_0p50 | 0.3x Silverman | 1.0000 | 0.0000 |
| floor_20 | robust | 1.0000 | 0.0000 |
| floor_20 | ISJ | **0.9999** (and the *shift* parameter also breaks: 0.9998) | 0.0001 |
| count_8 | ISJ | 0.9910 | 0.0000 |

None of them helps, because the bandwidth is not the problem. **ISJ actively misbehaves on tied data**: it
detects the discreteness and selects `h = 2e-5`, inflating the density estimate by up to 10^4x
(`f_ship/f_latent = 10195` for `count_8`, `3527` for `round_0p50`), and it breaks a parameter Silverman
left intact. Do not ship ISJ as a fix.

### 5.3 Sample size does not fix it

At N = 80,000, `round_0p50` and `floor_20` are still at coverage 1.0000 / Type-I 0.0000. The failure is
structural, not asymptotic. Effect size is irrelevant too (s0 = 0.02 and 0.18 both broken).

### 5.4 Two candidate guards, both computable without knowing the truth

**Guard A — tie-drop.** Drop any tau whose sample tie fraction at `q_tau` exceeds 1%. (`src/keepset.py`)

**Guard B — quantile de-duplication.** Drop only a tau whose sample `q_tau` collides with the previous
kept tau's. Strictly weaker; removes exactly the rank deficiency. (`src/keepset_q.py`)

| cell | full cov95 / T1 | **tie-drop cov95 / T1** | qdedup cov95 / T1 |
|---|---|---|---|
| round_0p10 | 0.9433 / 0.0567 | 0.9465 / 0.0535 | (not triggered) |
| round_0p50 | 0.9973 / 0.0000 | (n/a) | 0.9523 / 0.0460 |
| round_1p00 | 1.0000 / 0.0000 | (n/a) | 0.9479 / 0.0521 |
| floor_05 | 0.9997 / 0.0003 | **0.9460 / 0.0540** | 0.9450 / 0.0550 |
| floor_20 | 1.0000 / 0.0000 | **0.9470 / 0.0530** | 0.9456 / 0.0544 |
| floor_40 | 1.0000 / 0.0000 | **0.9478 / 0.0522** | 0.9458 / 0.0542 |
| ceil_20 | 0.9971 / 0.0029 | **0.9465 / 0.0535** | **0.3317 / 0.0016** |
| count_4 | 1.0000 / 0.0000 | (n/a) | **0.3366 / 0.0000** |
| count_8 | 1.0000 / 0.0000 | (n/a) | **0.7578 / 0.0000** |

**Guard A works and Guard B is dangerous.** Tie-drop restores null coverage to 0.9460–0.9478 (control:
0.9437) *and* restores causal coverage of the true scale parameter to 0.9405–0.9525 with bias −0.07 to
−0.20 sd. Quantile de-duplication restores the *null* calibration more widely — but because it keeps taus
whose quantile still sits inside a tie block, it leaves the estimand shifted (causal coverage stays
0.000–0.08), and on `ceil_20` and `count_4` it retains degenerate taus and turns a conservative failure
into an **anti-conservative** one (coverage 0.33). The weaker guard is the worse guard.

Tie-drop's limitation is that it needs at least 3 surviving taus, so it cannot rescue `round_0p25` and
worse, `count_8`, `count_4`, or `bimodal_round` — for those there is no usable tau set and the answer must
be refusal, not repair.

### 5.5 Recommended Stage 1 diagnostic

Stage 1 already computes everything needed. Emit, per tau:

1. `tie_frac(tau)` = fraction of the sample exactly equal to `q_tau`.
2. `n_distinct(Y)` and whether any two requested taus share a `q_tau`.
3. `h / local_spacing` — the bandwidth over the gap between neighbouring distinct values. This is a clean
   continuous discreteness detector: it is 10^2–10^4 for every continuous trait tested and **≤ 1.5** for
   every tied one (`round_0p10` 1.50, `round_0p25` 0.60, `round_0p50` 0.30, `count_4` 0.13).

Proposed policy:
- **any two taus sharing a `q_tau`** → refuse GLS, or drop to tie-drop, and say so;
- **`tie_frac(tau) > ` the gap to the neighbouring tau** → warn and drop that tau;
- **`h / spacing < 2`** → warn that the trait is discrete and the density estimate is not meaningful;
- OLS (`tau_cov_estimator = "ols"`) should be the default fallback: it was correctly calibrated in every
  single cell tested.

---

## 6. What I could not test

- **Real phenotypes.** Everything here is simulated. The antibody titres this was commissioned for have
  not been checked against these thresholds — that is the obvious next step and it is cheap (tie fraction
  per tau is a one-line computation on the existing phenotype files).
- **LD and relatedness.** Genotypes are independent binomial draws. No LD block structure, no cryptic
  relatedness, no LOCO interaction.
- **Other tau grids.** The boundary result is expressed *relative to the tau spacing of the grid used*
  (15 taus, minimum spacing 0.03). A denser grid collides sooner; the rule generalises but the numeric
  threshold does not.
- **A tie-aware density estimator.** Only bandwidth rules were varied. A genuinely discrete-aware
  estimator (e.g. a smoothed empirical mass function, or the mid-distribution transform) was not built or
  tested, and is the one remaining candidate fix.
- **N > 80,000, MAF < 0.05, and W with more than two columns.**
- **Far-tail rates below 1e-6.** With 200,000 null SNPs, 1e-6 has 0.2 expected hits; the `ceil_20` 270x
  figure is real but its magnitude has ~14% MC error and the extrapolation to 5e-8 is not measured.

---

## Appendix — provenance

52 Stage 1 cell runs were inherited complete from a previous agent and **re-used, not re-run**: runA (15
shapes x 2000 reps), runB (8 shapes x 3 bandwidth rules x 500 reps), runC (5 shapes at N = 80,000), runD
(4 shapes x 2 effect sizes x 1000 reps) — about 5.3M Stage 1 SNP fits, plus its script suite, its private
`fungwasStage2` install and its bandwidth-patch harness. Stage 2 had never been run: this brief adds the
52-job Stage 2 array (collate → tau guard → `param_gwas` x {jackknife, pooled} x {ols, gls} x {full,
tie-drop} → summarise), a 20-job array for the quantile-dedup guard, one bug fix in `analyse.py`, and the
analysis.
