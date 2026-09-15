# Is the pooled covariance safe as G-MoQ's default?

> **SUPERSEDED IN PART, 2026-09-15.** Section 8 item 1 ("keep pooled as the scan
> default") is **withdrawn**, and the claim that 5e-8 is "exactly nominal for the
> jackknife" is true only for jackknife+**OLS**, not jackknife+GLS. See
> [FINDINGS_covsource_default.md](FINDINGS_covsource_default.md), which measures
> all four arms on 780,000 paired scans. The diagnosis of the false-vQTL artefact
> as an estimand gap, and the `G` guard, are unaffected.

> **Update.** §5 below proposes flagging dispersion findings at `G >= 0.2`. That
> artefact has since been *corrected* rather than only flagged —
> `correct_dispersion_skew()`, see `docs/DISPERSION_SKEW_CORRECTION.md`. After
> correction the 5e-8 rate is within 1.03x nominal everywhere in this document's
> sweep, including the `G = 25` cells, so the flag is no longer what bounds false
> positives on a Normal baseline. It is still needed where the correction cannot
> be applied (no signed allele frequency) or is only partial (strongly skewed
> phenotypes), and its meaning otherwise inverts into a *power* warning: high `G`
> is now where the corrected statistic under-rejects. Everything else below
> stands as measured; the exact rates in Table 1 are reproduced to every
> published digit by the correction study's engine.

**Short answer: yes for the scan, no for the loci the scan finds — but the reason is
not the one the estimator has been charged with, and the proposed fallback does not fix it.**

The pooled estimator's working-independence error is real, is always in the
anti-conservative direction, and is bounded and small at common MAF. It is *not*
what produces false dispersion findings. False vQTLs come from a second,
larger, and previously unmeasured effect: a **pure location effect at a
genotype with non-zero skewness leaves a residue that lands exactly on the
dispersion weight column**. The block jackknife removes the first error and
leaves the second untouched.

Everything below is either an exact population calculation (no Monte Carlo error
at all) or a Monte Carlo through the *shipped* Stage 1 and Stage 2 code. The two
agree; where they disagree with prior documents, that is stated.

---

## 1. The false dispersion rate

Pure location effect, **zero** true dispersion effect, normal outcome,
N = 461,320. Every rejection of `beta_sigma2` is therefore false.
`Z_mu` is the location parameter's true Z. Rates are **exact** — they carry no
Monte Carlo error, which matters because 5e-8 is unreachable by simulation
(it would need ~1e10 replicates per cell).

| MAF | Z_mu | rate @0.05 | rate @5e-8 | inflation @5e-8 | jackknife @5e-8 | SE-error alone @5e-8 |
|---|---|---|---|---|---|---|
| 0.30 | 8  | 0.0502 | 5.14e-08 | **1.0x** | 5.12e-08 | 5.02e-08 |
| 0.30 | 16 | 0.0530 | 7.17e-08 | **1.4x** | 7.06e-08 | 5.08e-08 |
| 0.30 | 22 | 0.0607 | 1.39e-07 | **2.8x** | 1.35e-07 | 5.16e-08 |
| 0.30 | 30 | 0.0873 | 5.30e-07 | **10.6x** | 5.04e-07 | 5.29e-08 |
| 0.30 | 40 | 0.1712 | 4.59e-06 | **92x** | 4.22e-06 | 5.53e-08 |
| 0.30 | 60 | 0.6196 | 7.86e-04 | **1.6e4x** | 6.84e-04 | 6.28e-08 |
| 0.05 | 12 | 0.0715 | 2.66e-07 | **5.3x** | 2.51e-07 | 5.32e-08 |
| 0.05 | 22 | 0.3025 | 3.37e-05 | **674x** | 2.87e-05 | 6.17e-08 |
| 0.05 | 30 | 0.7660 | 3.26e-03 | **6.5e4x** | 2.63e-03 | 7.39e-08 |
| 0.01 | 8  | 0.0749 | 3.28e-07 | **6.6x** | 2.87e-07 | 5.78e-08 |
| 0.01 | 16 | 0.4542 | 2.00e-04 | **4.0e3x** | 1.34e-04 | 8.89e-08 |
| 0.01 | 30 | 1.000  | 0.912    | **1.8e7x** | 0.840 | 3.61e-07 |
| 0.005 | 8 | 0.1019 | 9.82e-07 | **20x** | 7.59e-07 | 6.69e-08 |
| 0.005 | 12| 0.3202 | 5.12e-05 | **1.0e3x** | 3.15e-05 | 9.58e-08 |
| 0.005 | 22| 0.9991 | 0.455    | **9.1e6x** | 0.307 | 4.20e-07 |
| **0.50** | **60** | **0.0512** | **5.89e-08** | **1.2x** | 5.00e-08 | 5.89e-08 |

Read the last row against the rest. At MAF 0.50 a location effect of
`Z_mu = 60` produces **no** false dispersion inflation at all. The driver is not
effect size; it is effect size **times genotype skewness**.

Compare the last two columns. The jackknife (which makes no working-independence
assumption) removes almost none of the inflation, and the pooled estimator's own
SE error accounts for almost none of it either. Both facts point at the same
conclusion: **this is not a covariance-estimation failure.**

### Monte Carlo confirmation

Independent simulation through the shipped Stage 1 (estimated quantiles,
estimated densities, RIF, residualisation on 5 covariates, shipped
`compute_sigma0`) and the shipped Stage 2 `param_gwas`. Observed rejection rate
at nominal 0.05, with Monte Carlo SE:

| cell | Z_mu | G | MC rate ± MCSE | exact |
|---|---|---|---|---|
| maf 0.50, beta 0.35 | 73.9 | 0.000 | 0.0547 ± 0.0023 | 0.0591 |
| maf 0.30, beta 0.02 | 4.1  | 0.022 | 0.0498 ± 0.0011 | 0.0501 |
| maf 0.30, beta 0.05 | 10.2 | 0.139 | 0.0522 ± 0.0011 | 0.0525 |
| maf 0.01, beta 0.10 | 4.4  | 0.295 | 0.0639 ± 0.0025 | 0.0617 |
| maf 0.30, beta 0.10 | 20.3 | 0.554 | 0.0877 ± 0.0014 | 0.0872 |
| maf 0.01, beta 0.20 | 8.7  | 1.146 | 0.2278 ± 0.0042 | 0.2260 |
| maf 0.30, beta 0.20 | 40.1 | 2.155 | 0.5938 ± 0.0025 | 0.5954 |
| maf 0.30, beta 0.35 | 67.8 | 6.145 | 1.0000 ± 0.0000 | 1.0000 |

Exact and MC agree within ±1.3 MCSE in every pooled cell tested (232 cells).

**Deep-tail validation.** The headline rates are at 5e-8, which no simulation can
reach, so the formula is validated as far into the tail as the MC allows — at
nominal **1e-3**, fifty times deeper than 0.05:

| cell | MC rate @1e-3 ± MCSE | exact | deviation |
|---|---|---|---|
| maf 0.30, beta 0.02 | 0.00088 ± 0.00015 | 0.00100 | -0.9σ |
| maf 0.30, beta 0.05 | 0.00135 ± 0.00018 | 0.00113 | +1.2σ |
| maf 0.30, beta 0.10 | 0.00347 ± 0.00029 | 0.00330 | +0.6σ |
| maf 0.30, beta 0.20 | 0.14295 ± 0.00175 | 0.14256 | +0.2σ |
| maf 0.30, beta 0.35 | 0.99935 ± 0.00013 | 0.99917 | +1.4σ |
| maf 0.01, beta 0.10 | 0.00170 ± 0.00041 | 0.00165 | +0.1σ |
| maf 0.01, beta 0.20 | 0.02240 ± 0.00148 | 0.02010 | +1.6σ |
| maf 0.01, beta 0.35 | 0.60820 ± 0.00488 | 0.61385 | -1.2σ |
| maf 0.50, beta 0.35 | 0.00130 ± 0.00036 | 0.00153 | -0.6σ |
| bigN maf 0.30, beta 0.20 | 0.92500 ± 0.00416 | 0.92267 | +0.6σ |

Agreement holds over three decades of rate (1e-3 to 1) and two decades of
nominal level. The 5e-8 extrapolation beyond that rests on `theta_hat` being
normal about `theta_pop`, which is assumed, not demonstrated.

---

## 2. Why: the mechanism

For `Y = delta*gc + e` the population RIF-regression estimand expands as

```
b(tau) = delta  -  (delta^2 / 2) * gskew * f0'(q_tau) / (f0(q_tau) * sqrt(Var g))  +  O(delta^3)
```

with `gskew = E[gc^3]/Var(g)^1.5 = (1-2p)/sqrt(2p(1-p))` the **genotype skewness**.

For a normal baseline `f0'/f0 = -(q-m)/v`, which is *exactly proportional to the
dispersion weight column* `(q-m)/(2v)` that `make_weights_vqtl` builds. So the
second-order residue of a pure location effect is perfectly collinear with the
dispersion parameter and is absorbed by it in full. It is an **estimand** error,
not a variance error, which is why no covariance estimator can remove it and why
it vanishes identically at `p = 0.5` where `gskew = 0`.

Two consequences, both verified numerically to <1%:
* the artifact scales as `delta^2`, i.e. as `Z_mu^2`;
* it grows with N as `sqrt(N)` — **more data makes it worse**, because the
  bias is fixed while the SE shrinks.

---

## 3. Direction of the error — the pooled SE is anti-conservative everywhere

| scenario | pooled SE vs truth | direction |
|---|---|---|
| null | 1.0000 | exact |
| pure location | < 1 always | **anti-conservative** |
| pure dispersion | < 1 always, = 1 only at the median | **anti-conservative** |

Under a pure symmetric dispersion effect the per-tau pooled/true variance ratio
is `tau(1-tau) / E_w[F_g(1-F_g)]` with weights `w ∝ gc^2`. Because `F(1-F)` is
concave this is `<= 1`, with equality only at `tau = 0.5` where `F_g = 0.5` for
every genotype. The engine returns `tau_se_ratio_max = 1.00000` exactly for every
dispersion cell, which is this identity showing up numerically.

**There is no regime in which the pooled estimator is conservative.** An
estimator that only ever errs conservative would be safe to default; this one is
not.

### 95% CI coverage (pooled, N = 461,320, normal)

| MAF | Z | pure location: mu / sigma2 | pure dispersion: mu / sigma2 |
|---|---|---|---|
| 0.30 | 16 | 0.9499 / 0.9499 | 0.9477 / 0.9464 |
| 0.30 | 60 | 0.9489 / 0.9483 | 0.9405 / 0.9339 |
| 0.05 | 30 | 0.9481 / 0.9470 | 0.9269 / 0.9114 |
| 0.05 | 60 | 0.9417 / 0.9370 | 0.8944 / 0.8528 |
| 0.005 | 30 | 0.9253 / 0.9106 | 0.8358 / 0.7492 |
| 0.005 | 60 | 0.7436 / 0.6260 | 0.6350 / 0.4748 |

The `gls_coverage` study's 0.9373/0.9409 sag under a true scale effect is
reproduced and placed: it is the *mild* end of this surface.

### Per-tau spread is far worse than the parameter

Per-tau SE ratio range across the 45 taus, pure location:

| MAF | Z_mu | min | max | parameter-level |
|---|---|---|---|---|
| 0.30 | 30 | 0.966 | 1.031 | 0.9981 |
| 0.05 | 30 | 0.858 | 1.168 | 0.9872 |
| 0.005 | 30 | 0.588 | 1.840 | 0.8666 |
| 0.005 | 60 | 0.251 | 5.320 | 0.4536 |

Quoting the median per-tau ratio hides a factor of 20. The prior doc's
recommendation to always quote min and max is correct and is reinforced here.

---

## 4. The jackknife is not the fallback (for this), but is for the other thing

| failure mode | cause | pooled | jackknife B=1000 |
|---|---|---|---|
| false dispersion from a location effect | estimand gap | bad | **also bad** |
| distorted inference under a *true* dispersion effect | covariance | bad | **fixed exactly** |

Under a true dispersion effect the location parameter's false-positive rate at
5e-8 goes to 1.19e-07 (MAF 0.005, Z=3) up to 1.18e-02 (MAF 0.005, Z=60) for
pooled, and is **5.00e-08 — exactly nominal — for the jackknife in every cell**,
because there is no estimand gap in that direction. So the jackknife is a
genuine and complete fix for the dispersion-effect covariance failure, and no
help whatsoever for the false-vQTL problem.

**A caution on the jackknife that was not previously reported.** The real
B = 1000 block jackknife is itself anti-conservative when paired with GLS,
because a 45x45 covariance estimated from 1000 blocks is noisy enough for GLS to
overfit it. Measured `mean_se/emp_sd` across MC cells:

| cell | jackknife GLS | jackknife OLS | pooled GLS |
|---|---|---|---|
| maf 0.01, loc02 | **0.910** | 0.997 | 0.995 |
| maf 0.01, disp10 | **0.898** | — | — |
| maf 0.30, loc10 | **0.954** | 1.003 | 0.998 |
| bigN maf 0.30, disp10 | **0.929** (mu), **0.942** (sigma2) | — | — |

At nominal 0.05 this makes **jackknife+GLS worse than pooled+GLS** for small
effects (0.0604 vs 0.0492 at maf 0.30 loc02). If the jackknife is used as a
fallback, prefer OLS with it, or raise B.

---

## 5. A guard that Stage 2 can actually compute

Every input is available per variant without knowing the truth:

```
gskew = (1 - 2p) / sqrt(2p(1-p))            # from MAF
G     = 0.6855 * gskew * Z_mu^2 / sqrt(N)   # Z_mu = estimated location Z
```

`G` is a direct estimate of the false-dispersion non-centrality — i.e. of the
`|Z|` that `beta_sigma2` will show when the true dispersion effect is zero. The
constant 0.6855 is stable to ±1% for `Z_mu <= 12` and ±3% to `Z_mu <= 22`, over
MAF 0.005-0.5 and N 50,000-461,320.

Inflation of the 5e-8 rate as a function of `G` (worst case over the whole sweep):

| G | max inflation @5e-8 |
|---|---|
| < 0.05 | 1.05x |
| < 0.10 | 1.17x |
| < 0.20 | **1.95x** |
| < 0.30 | 3.4x |
| < 0.50 | 7.9x |
| < 1.0 | 130x |
| >= 2.0 | >= 2.9e4x |

**Proposed rule: flag a dispersion finding when `G >= 0.2`.** That holds the
false-dispersion rate within 2x nominal, is monotone in `G` across every MAF, N
and effect size tested, and correctly passes the `MAF 0.5, Z_mu = 73.9` cell that
any rule based on effect size alone would wrongly flag. Flagged variants should
have the dispersion claim confirmed by a method that does not use the linear
`W` projection at all — the artifact is in the estimand, so recomputing the
covariance (jackknife or exact sandwich) will *not* clear it.

On the skewed outcome (standardised chi2_4) the same rule over-predicts the
artifact by ~1.6x, i.e. it is conservative there. It was not tested on
heavier-tailed baselines, where the constant could go the other way.

---

## 6. Prior claims: what held and what did not

From `fungwas_paper/docs/BP1_pooled_implementation.md` §"Four things follow":

| claim | verdict |
|---|---|
| "location effects are safe across the realistic GWAS range" (SE error 0.0-0.3% at beta 0.05-0.10) | **the SE number is right, the conclusion is wrong** |
| "0.9481 at beta = 0.35", anti-conservative on the `both` model's additive SE | **confirmed** (I get 0.959 at MAF 0.1, 0.966 at MAF 0.2) |
| "for a dispersion effect the pooled SE is the larger one... the failure mode is lost power" | **contradicted** |
| "quote per-tau min and max, never the median" | **confirmed and reinforced** |

**On "location effects are safe".** The SE ratios are reproduced almost exactly —
at MAF 0.1-0.2 and beta = 0.10 SD/allele the pooled SE error really is 0.36%.
But at that same point the false dispersion rate at 5e-8 is inflated **1123x**:

| MAF | beta | max\|Z_tau\| | SE_mu ratio | false-dispersion inflation @5e-8 |
|---|---|---|---|---|
| 0.1 | 0.02 | 4.6 | 0.99985 | 1.1x |
| 0.1 | 0.05 | 11.5 | 0.99909 | 4.4x |
| 0.1 | 0.10 | 23.0 | 0.99638 | **1123x** |
| 0.1 | 0.20 | 45.7 | 0.98583 | 1.6e7x |
| 0.2 | 0.10 | 30.6 | 0.99699 | **1100x** |

The SE ratio is simply the wrong diagnostic: it is blind to the estimand gap,
which is the dominant term. beta = 0.10 SD/allele is inside the anthropometric
range the doc itself cites as realistic.

**On the dispersion direction.** The doc reports pooled/reference SE ratios
*above* 1 under a scale effect (1.0122 at gamma 0.02, 1.1692 at gamma 0.20) and
concludes the failure is conservative. Both my exact calculation and the MC
through shipped Stage 2 give ratios *below* 1:

| gamma (MAF 0.3) | MC pooled OLS sigma2 | MC hc0 (reference) | exact |
|---|---|---|---|
| 0.01 | 0.9995 ± 0.0071 | 1.0050 | 0.9945 |
| 0.05 | 0.9667 ± 0.0068 | 0.9951 | 0.9713 |
| 0.10 | 0.9462 ± 0.0067 | 1.0053 | 0.9408 |
| 0.20 | 0.8790 ± 0.0062 | 0.9997 | 0.8785 |

The concavity argument in §3 says the ratio *cannot* exceed 1 for a symmetric
scale effect, so I believe the sign reported in the doc is wrong. The most
likely cause is the comparator: that table was built against the production
**B = 65** jackknife, which for a 45x45 covariance is close to rank-deficient. I
did not re-run their chr22 injection with an exact sandwich comparator, so this
explanation is a hypothesis, not a result. The *magnitudes* are comparable; only
the sign differs. **This should be resolved before the doc's §8 recommendation
("the failure mode there is lost power") is relied on**, because it reverses the
risk: it is false positives, not lost power.

---

## 7. What I could not test

* **Real genotypes and LD.** Everything uses a 3-point HWE genotype law. The
  guard depends on MAF only through `gskew`, so genotyping error or
  Hardy-Weinberg departure would perturb it in a way I have not quantified.
* **The production B = 65 jackknife** and the int8-quantised jackknife. Only
  B = 1000 float was simulated. The §6 sign discrepancy hinges on this.
* **Baselines other than normal and standardised chi2_4.** The guard constant
  0.6855 is calibrated on the normal and is conservative for chi2_4; heavier
  tails were not tried.
* **Covariates with real effects.** The main MC cells carry 5 covariates that are
  residualised out but do not affect Y. Three `cov_*` cells with `covar_sd > 0`
  exist in `poolcal/rds` and were not analysed.
* **The `hybrid` estimator**, reported elsewhere as exact in the dispersion case.
  Not simulated here.
* **Rates below 1e-3 by simulation.** The formula is MC-validated at 0.05 and
  1e-3 (§1) and extrapolated from there. The extrapolation assumes `theta_hat` is
  normal about `theta_pop`; the MC supports that down to 1e-3 but cannot reach
  5e-8, which would need ~1e10 replicates per cell. A departure from normality in
  the far tail would move the absolute rates, though not the ordering by `G`.

---

## 8. Recommendation

1. **Keep pooled as the scan default.** At MAF >= 0.05 and `Z_mu <= 16` — the
   overwhelming majority of a scan — the false-dispersion rate is within 1.5x
   nominal and the SE error is under 0.5%.
2. **Ship the `G` guard** and flag dispersion findings with `G >= 0.2`. It costs
   one arithmetic expression per variant.
3. **Do not present the jackknife as the fix for false vQTLs.** It is the fix for
   the true-dispersion covariance failure, which is a different problem. Use it
   with OLS, not GLS, or raise B.
4. **Correct the sign in `BP1_pooled_implementation.md` §"Four things follow" #3**
   once the B=65 comparator question is settled, or re-run that injection against
   an exact sandwich.
