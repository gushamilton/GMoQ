# Which covariance source should Stage 2 default to?

**Campaign:** `sim/covsource` (BP `/user/work/fh6520/covsource_run`), 52 cells x
15,000 replicates = 780,000 simulated scans, N = 461,320, `--blocks 1000`,
45 taus, paired replicates across all four `(cov_source, tau_cov_estimator)`
combinations. Aggregates: `sim/covsource/tables/`.

**Answer: `cov_source = "jackknife"` with `tau_cov_estimator = "ols"`.**

This *reverses* the recommendation in `FINDINGS_pooled_dispersion.md` section 8
item 1, which is withdrawn below.

---

## 1. Why this is not a trade-off

Under OLS the point estimate never touches the covariance: `b_hat = (W'W)^-1 W'b`
has no `Sigma` in it. The simulation confirms it numerically -- the ratio of
empirical sd(estimate) between the jackknife and pooled arms is **1.00000 in all
44 main cells**.

So the two arms produce *identical estimates* and differ only in the SE they
report. There is no efficiency to trade away. One is simply better calibrated.

## 2. The calibration gap, pure dispersion, normal phenotype

`beta_sigma2`, 95% CI coverage (table 1) and SE ratio = mean reported SE /
empirical sd, where < 1 is anti-conservative (table 5):

| MAF | Z | cov95 pooled+OLS | cov95 **jk+OLS** | SE ratio pooled | SE ratio **jk** |
|---|---|---|---|---|---|
| 0.30 | 16 | 0.9467 | **0.9507** | 0.9900 | **1.0040** |
| 0.05 | 16 | 0.9285 | **0.9447** | 0.9206 | **0.9860** |
| 0.01 | 16 | 0.9027 | **0.9491** | 0.8442 | **0.9953** |
| 0.01 | 30 | 0.8475 | **0.9452** | 0.7339 | **0.9961** |
| 0.005 | 16 | 0.8742 | **0.9488** | 0.7862 | **0.9940** |
| 0.005 | 30 | 0.7919 | **0.9437** | 0.6519 | **0.9955** |

The jackknife holds 0.944-0.951 everywhere below the unphysical `Z60` cell.
Pooled decays monotonically with rarity and effect size.

## 3. Type-I inflation at genome-wide significance

Extrapolated to 5e-8, `beta_mu` under a true dispersion effect (table 9):

| MAF | Z | pooled+OLS | **jk+OLS** |
|---|---|---|---|
| 0.30 | 16 | 1.26x | **0.94x** |
| 0.01 | 30 | 295x | **1.08x** |
| 0.005 | 30 | 2086x | **0.99x** |

## 4. The null is not where this shows up

At the null both OLS arms are exactly nominal and indistinguishable -- 0.047-0.054
at nominal 0.05 across every MAF and both phenotype shapes (table 6). **Null
simulations cannot detect this defect.** It only appears once a real dispersion
effect is present, which is why it survived earlier validation.

Note also that **GLS + jackknife is inflated at the null**: 0.060-0.068. That arm
is broken before any effect exists.

## 5. Skewed phenotype

Standardised gamma(2), true dispersion, `beta_sigma2` coverage (table 4):

| MAF | pooled+OLS | **jk+OLS** | pooled+GLS | jk+GLS |
|---|---|---|---|---|
| 0.01 | 0.8988 | **0.9472** | 0.6859 | 0.8068 |
| 0.30 | 0.9450 | **0.9479** | 0.8189 | 0.8066 |

GLS collapses on skew under *either* covariance. Only OLS + jackknife holds.

## 6. False vQTLs under a pure location effect

Extrapolated to 5e-8 (table 16). Two distinct mechanisms are superimposed here
and they must not be conflated:

| cell | pooled+OLS | jk+OLS |
|---|---|---|
| normal MAF 0.30 Z16 | 0.80 | 0.78 |
| normal MAF 0.05 Z16 | 1.35 | 1.11 |
| normal MAF 0.01 Z16 | 2.29 | 0.79 |
| normal MAF 0.005 Z30 | 35.24 | 0.05 |
| gamma2 MAF 0.01 Z16 | 6.34 | 0.30 |

- The **estimand artefact** (genotype skewness, `DISPERSION_SKEW_CORRECTION.md`)
  is what remains where both arms agree and exceed 1 -- roughly 1.1-1.35x at
  MAF 0.05. The jackknife does **not** fix this, and the `G` guard remains the
  tool for it.
- The **covariance failure** is the rest, and it is the larger term at rare MAF.
  That part the jackknife does fix.

At MAF 0.005 / Z30 the jackknife swings to conservative (0.05x, ~20x
under-rejection). Both arms fail at `Z60`, which is not a physical variant --
see section 8.

## 7. What the jackknife costs

1. **Reported-SE noise.** `sd(SE)/mean(SE)` is 0.023-0.029 for the jackknife
   versus 0.005-0.008 for pooled (table 17) -- roughly 4x noisier, on a small
   base.
2. **Power.** At matched *actual* Type-I, 0.2-7% relative power (tables 11-12).
   That is the honest comparison; an unmatched power comparison would be
   rewarding the anti-conservative arm for its own miscalibration.
3. **Compute.** 3.73-3.86 s and ~345 MiB versus 2.26-2.69 s and ~168 MiB on full
   chr22 (paper benchmark `stage2-benchmark-20260915-05`). Stage 1 already
   computes the jackknife either way, so this is Stage 2 read cost only.

None of these is close to the calibration gap in sections 2-3.

## 8. On the `Z60` cells

`FINDINGS_pooled_dispersion.md` leads with `MAF 0.005 / Z 60`. That cell requires
the minor homozygote to carry about **12x the phenotypic SD** and is not a real
variant; both arms fail there and it should not drive the recommendation. The
cells that *are* plausible for a real scan -- MAF 0.01-0.05, Z 10-30 -- are
exactly where pooled is already 2-295x inflated and the jackknife is nominal.

## 9. Corrections to `FINDINGS_pooled_dispersion.md`

1. **Section 8 item 1 ("Keep pooled as the scan default") is withdrawn.** It is
   defensible only for `beta_mu` at common MAF; it does not hold for
   `beta_sigma2`, which is the parameter the method exists to estimate.
2. Its claim that 5e-8 is "exactly nominal for the jackknife" is
   **estimator-specific and does not say so**. True for jackknife+**OLS**
   (0.94-1.08x). False for jackknife+**GLS**: 3.45-4.95x at MAF 0.30 and 457x at
   MAF 0.005/Z30, and inflated at the null before any effect exists.
3. What survives unchanged: the diagnosis of the false-vQTL artefact as an
   estimand gap, and the `G` guard as its remedy (section 6 above).

## 10. Status of the default

The evidence above is unambiguous, but flipping `cov_source` changes the
**primary analysis of the paper**, which currently pins pooled. The default has
therefore **not** been changed in code by this commit; the decision is recorded
in `BACKLOG.md` item B1 as reopened pending that call.
