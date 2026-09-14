# Does the calibration measured at N = 30,000 survive to production N = 461,320?

Status: measured 2026-09-14/15 on 36 simulation cells x 20,000 replicates, at
N = 30,000 / 100,000 / 461,320. Companion to `OLS_VS_GLS.md`, which established
the same quantities at N = 30,000 only.

## Bottom line

**Under the null, yes — exactly, and at every tau, MAF and MAC tested. Under a
true scale effect, no.**

The "sag to ~0.94 under a true scale effect" reported at N = 30,000 is not a
constant of the method. It is misspecification bias, and it is measured in
standard errors, which shrink. The absolute misfit `max|beta_tau - W theta|` is
a property of the data-generating process and is **identical at every N** — for
the `ADD` fit under a scale effect it is 0.12162 at all three rungs, to five
decimals. The per-tau standard error falls like `N_eff^-1/2`. So lack of fit
*in SE units* grows like `sqrt(N_eff)`, measured at 10.92 / 19.93 / 42.81 —
ratios 1.83 and 3.92 against `sqrt(N/30,000)` = 1.83 and 3.92. Exact.

Coverage follows it down. The worst case is the dispersion parameter under GLS:

| `LS:beta_sigma2`, true scale effect, MAF 0.30 | N=30,000 | N=100,000 | N=461,320 |
|---|---|---|---|
| OLS | 0.9440 | 0.9434 | 0.9413 |
| **GLS** | **0.9316** | **0.9064** | **0.8437** |

Monte Carlo error 0.0011–0.0018. A 95% interval that was 1.8 points light at
simulation N is **10.6 points light at production N**, on the parameter a vQTL
claim rests on.

**The governing variable is `N_eff = 2Np(1-p)`, not `N`.** Held at constant MAC
(arm C, MAC ~ 922 at every rung) nothing moves at all: lack of fit 2.40 / 2.42 /
2.43, coverage 0.9364 / 0.9360 / 0.9370. Nothing in this pipeline is sensitive
to sample size as such. It is sensitive to how much information the genotype
carries, and more information resolves the model's inadequacy rather than
averaging it away.

Practical reading: **OLS transfers; GLS does not.** OLS loses at most 0.3
coverage points between simulation and production N anywhere tested. GLS loses
up to 10.6. This strengthens, rather than merely repeats, `OLS_VS_GLS.md`'s
recommendation to keep OLS the default — and the lack-of-fit guard it asked for
is now not optional but load-bearing, because the regime it guards against is
the production regime, not the simulation one.

---

## 1. Design

Three rungs x three arms, 20,000 replicates each, so Monte Carlo error on a 0.95
rate is 0.0015 per cell everywhere and the rungs are directly comparable.

| arm | MAF | MAC at 30k / 100k / 461,320 | what it isolates |
|---|---|---|---|
| A | 0.30 | 18,000 / 60,000 / 276,792 | common variant; `N_eff` grows with N |
| B | 0.01 | 600 / 2,000 / 9,226 | constant MAF; `N_eff` grows with N |
| C | 922/(2N) | 922 / 922 / 922 | **constant MAC** — pure N at fixed genotype information |

Arm C is the arm the earlier work did not have. MAC 922 is production's rarest
chr22 variant (MAC 927). Because `SE_tau`, `N_eff` and the number of minor-allele
carriers below any given quantile are all fixed along arm C, anything that moves
along it is an N effect and anything that does not, is not.

DGPs `{null, scale (c=0.05), mixed (a=0.05, c=0.05)}` x shapes `{normal,
gamma(2)}`, each fitted under four models (`ADD`, `LS`, `LSK`, `TRU`) and both
estimators, against exact population truth. Coverage is counted, never inferred.

Pooled rates below carry a **replicate-clustered** Monte Carlo error: the
(fit, parameter) z-values within one replicate are correlated, so the unit of
independence is the replicate, not the z-value.

**Validation.** The N = 30,000 rung reproduces the completed `gls_coverage`
study on an independent seed stream: 160 matched (cell, fit, estimator,
parameter) rows, mean cov95 difference **+0.00056**, largest deviation **2.46**
Monte Carlo standard errors, none beyond 3.

---

## 2. Coverage under the null: transfers exactly

95% intervals, all fits and parameters pooled, MC error +-0.0007 (arms A, B) and
+-0.0010 (arm C).

| arm | est | N=30,000 | N=100,000 | N=461,320 |
|---|---|---|---|---|
| A (MAF 0.30) | OLS | 0.9507 | 0.9505 | 0.9496 |
| A | GLS | 0.9502 | 0.9508 | 0.9497 |
| B (MAF 0.01) | OLS | 0.9489 | 0.9506 | 0.9519 |
| B | GLS | 0.9491 | 0.9500 | 0.9517 |
| C (MAC 922) | OLS | 0.9510 | 0.9492 | 0.9497 |
| C | GLS | 0.9503 | 0.9492 | 0.9491 |

At 99%, every cell lies in 0.9889–0.9904 against a nominal 0.99. Per parameter
(16 fit x parameter x estimator combinations at arm A) every null entry lies in
0.9486–0.9525 at every rung, against MC error 0.0011.

**Type-I error and the genome-wide extrapolation**, from 320,000 null z-values
per cell (160,000 for arm C):

| arm | est | | N=30,000 | N=100,000 | N=461,320 |
|---|---|---|---|---|---|
| A | OLS | sd(z) | 0.9981 | 0.9991 | 1.0004 |
| | | T1 @ 0.05 | 0.0493 | 0.0495 | 0.0504 |
| | | implied p(5e-8) / 5e-8 | 0.94 | 0.97 | 1.01 |
| A | GLS | sd(z) | 0.9995 | 0.9986 | 1.0002 |
| | | implied ratio | 0.99 | 0.96 | 1.01 |
| B | OLS | sd(z) | 1.0038 | 0.9984 | 0.9953 |
| | | implied ratio | 1.12 | 0.95 | 0.87 |
| B | GLS | implied ratio | 1.11 | 0.98 | 0.88 |
| C | OLS | implied ratio | 0.92 | 1.19 | 1.02 |
| C | GLS | implied ratio | 0.91 | 1.17 | 1.02 |

sd(z) is within 0.005 of 1 everywhere, MC error 0.0012–0.0018. The genome-wide
column spans 0.87–1.19 x nominal against the 0.99–1.03 previously reported, but
that is not a degradation with N — it is the extrapolation's own fragility. The
ratio is a 5.45-sigma extrapolation, so `d log(ratio)/d log(sd(z)) ~ 30`: the MC
error on sd(z) alone puts +-0.04 on arm A/B and +-0.06 on arm C, and the
assumption of exact normality in the far tail is doing the rest. **Treat 0.87 to
1.19 as the honest resolution of this method at 20,000 replicates, not as
evidence of an N trend** — there is no monotone pattern in N in any arm.

---

## 3. Coverage under a true scale effect: does not transfer

95% intervals, all fits and parameters pooled.

| arm | est | N=30,000 | N=100,000 | N=461,320 | change |
|---|---|---|---|---|---|
| A (MAF 0.30) | OLS | 0.9445 | 0.9452 | 0.9430 | -0.0015 |
| A | **GLS** | **0.9406** | **0.9353** | **0.9198** | **-0.0208** |
| B (MAF 0.01) | OLS | 0.9382 | 0.9367 | 0.9377 | -0.0005 |
| B | GLS | 0.9336 | 0.9323 | 0.9315 | -0.0021 |
| C (MAC 922) | OLS | 0.9384 | 0.9379 | 0.9385 | +0.0001 |
| C | GLS | 0.9364 | 0.9360 | 0.9370 | +0.0006 |

MC error 0.0007–0.0011. Per parameter at arm A, where the damage concentrates:

| fit:param | est | N=30,000 | N=100,000 | N=461,320 |
|---|---|---|---|---|
| `LS:beta_sigma2` | OLS | 0.9440 | 0.9434 | 0.9413 |
| `LS:beta_sigma2` | **GLS** | 0.9316 | 0.9064 | **0.8437** |
| `ADD:beta_add` | OLS | 0.9458 | 0.9476 | 0.9456 |
| `ADD:beta_add` | **GLS** | 0.9295 | 0.9203 | **0.8988** |
| `TRU:beta_scale` | OLS | 0.9445 | 0.9443 | 0.9414 |
| `TRU:beta_scale` | **GLS** | 0.9413 | 0.9335 | **0.9188** |
| `LSK:beta_kurt` | GLS | 0.9419 | 0.9387 | 0.9192 |
| `LS:beta_mu` | GLS | 0.9465 | 0.9472 | 0.9463 |

MC error 0.0011–0.0018. The location parameter is untouched at every N. Every
degrading entry is a **dispersion or shape** parameter under GLS. Interpolating,
`LS:beta_sigma2` under GLS crosses 0.90 at about N = 120,000 — inside the range
of samples this method is already being applied to.

Under `mixed` the same pattern appears at roughly a third the magnitude
(`LS:beta_sigma2` GLS 0.9412 / 0.9377 / 0.9288).

### Why

The absolute misfit is fixed; the yardstick shrinks.

| arm | fit | `max|resid|` (absolute) | in SE units: 30k / 100k / 461,320 |
|---|---|---|---|
| A | ADD | 0.12162 at every N | 10.92 / 19.93 / **42.81** |
| A | LS | 0.02208 at every N | 2.12 / 3.87 / **8.31** |
| A | TRU | 0.00265 at every N | 0.24 / 0.44 / **0.95** |
| B | ADD | 0.12420 at every N | 2.49 / 4.54 / **9.75** |
| C | ADD | 0.1080 / 0.1081 / 0.1081 | 2.40 / 2.42 / **2.43** |

The arm A ratios are 1.826 and 3.921 against `sqrt(N/30,000)` = 1.826 and 3.921.
Arm C, at constant MAC, does not move. Bias in units of the estimator's own
sampling sd tracks it:

| arm A, scale, `LS:beta_sigma2` | N=30,000 | N=100,000 | N=461,320 |
|---|---|---|---|
| OLS \|bias\|/sd | 0.081 | 0.105 | 0.134 |
| GLS \|bias\|/sd | 0.271 | 0.429 | 0.690 |

and at constant MAC it *falls* (GLS 0.049 / 0.026 / 0.014), which is a small-MAC
finite-sample bias washing out, not an N effect.

This is `OLS_VS_GLS.md`'s finding — GLS concentrates on the directions that
amplify misfit, OLS averages it away — with the N axis attached: the misfit is
invisible at N = 30,000 (under 1 SE for `TRU`, 2 for `LS`) and unmissable at
production. Note that **no fitted model is exactly correct under this scale
DGP**, `TRU` included: its residual is small but non-zero, which is why even it
degrades. That is the honest statement of the regime — "coverage under a true
scale effect" is always coverage under mild misspecification, and mild
misspecification is not scale-free.

---

## 4. The extreme taus improve, or hold — they do not degrade

`BP1_simulation_spec.md` flagged the tail taus as the N-sensitive piece by
construction. Under the null they are not sensitive at all. Stage 1 coverage of
`beta_tau` itself, 95%, MC error 0.0008–0.0016:

| tau | obs below (30k / 100k / 461k) | arm A null | arm C null (MAC 922) |
|---|---|---|---|
| 0.02 | 600 / 2,000 / 9,226 | 0.9498 / 0.9507 / 0.9488 | 0.9505 / 0.9510 / 0.9515 |
| 0.50 | 15,000 / 50,000 / 230,660 | 0.9509 / 0.9495 / 0.9503 | 0.9506 / 0.9498 / 0.9477 |
| 0.98 | 29,400 / 98,000 / 452,094 | 0.9513 / 0.9501 / 0.9507 | 0.9527 / 0.9479 / 0.9525 |

Flat, and at nominal, including at MAC 922 where only ~18 minor-allele carriers
sit below the tau = 0.02 quantile at every rung. Stage 1's own per-tau SE ratio
`mean(SE_tau)/sd(beta_tau_hat)` is 0.988–1.008 at every tau, arm and N.

Under a scale effect the extreme taus degrade with N exactly as section 3
predicts (arm A tau = 0.02: 0.9017 / 0.8550 / **0.7328**), and at constant MAC
they do not (arm C tau = 0.02: 0.9198 / 0.9212 / 0.9225). Same mechanism, not a
separate tail problem.

**The KDE density estimate is not the culprit.** `mean(fhat)/f_true - 1` at
tau = 0.98, arm A: +0.0235 / +0.0138 / +0.0073 (arm C tau=0.02: +0.0309 /
+0.0183 / +0.0100). It shrinks like `N^-0.42`, which is the textbook bias of a
kernel with bandwidth `h ∝ N^-0.2`. It does not reach coverage because the same
`fhat` divides both `beta_tau` and its standard error, and cancels — which is
what the per-tau SE ratio of 1.000 is showing.

---

## 5. SE ratio and efficiency

`mean(reported SE)/empirical sd`, restricted to well-fitting combinations:

| arm | N=30,000 | N=100,000 | N=461,320 |
|---|---|---|---|
| A (MAF 0.30) | 0.9876 | 0.9902 | 0.9996 |
| B (MAF 0.01) | 0.9776 | 0.9784 | 0.9866 |
| C (MAC 922) | 0.9776 | 0.9724 | 0.9763 |

OLS and GLS agree to within 0.002 throughout. The ~2% deficit is a **MAC**
effect, not an N effect: it closes as MAC grows (arm A, arm B) and is flat where
MAC is held fixed (arm C). At production's common variants the reported SE is
essentially exact; at production's rarest it is ~2.4% optimistic, the same as at
N = 30,000.

**Efficiency is N-stable.** Empirical sd(GLS)/sd(OLS), mean [min, max]:

| parameter | N=30,000 | N=100,000 | N=461,320 |
|---|---|---|---|
| `beta_add` | 0.810 [0.53, 1.00] | 0.808 [0.52, 1.00] | 0.807 [0.52, 1.00] |
| `beta_kurt` | 0.848 | 0.874 | 0.880 |
| `beta_scale` | 0.907 | 0.912 | 0.937 |
| `beta_sigma2` | 0.940 | 0.939 | 0.927 |
| `beta_mu` | 0.991 | 0.992 | 0.992 |

GLS keeps its advantage at production N. The advantage is real and it is stable —
it is simply not worth 10 coverage points on the dispersion parameter.

---

## 6. Effective degrees of freedom: an `N_eff` law, not an N law

Harvested from the shipped code via `n_blocks = 1000` and `n_eff = d_resid`:

| arm | parameter | N=30,000 | N=100,000 | N=461,320 |
|---|---|---|---|---|
| A (MAF 0.30) | `beta_mu` | 933.9 | 978.5 | 994.5 |
| A | `beta_scale` | 827.2 | 940.4 | 985.7 |
| B (MAF 0.01) | `beta_mu` | 402.8 | 691.7 | 911.3 |
| B | `beta_scale` | 184.7 | 430.1 | 776.5 |
| **C (MAC 922)** | `beta_mu` | **507.2** | **509.9** | **511.2** |
| **C** | `beta_scale` | **257.0** | **259.1** | **260.1** |

Arm C is flat to three significant figures across a 15x change in N. This
confirms `BP1_effective_df.md`'s law `1/nu = 1/(B-1) + c/N_eff` is entirely an
`N_eff` law: the apparent N dependence is `N_eff = 2Np(1-p)` moving, nothing more.

Production-relevant consequence: at chr22's rarest variants (MAC ~ 927) the
**scale**-parameter effective df is about **260**, not `B - 1 = 999`. The t_260
two-sided critical value at 5e-8 is 5.6165 against the normal 5.4514.

---

## 7. The Q lack-of-fit guard: calibrated, but with two gaps

Q is exactly calibrated under the null at every rung — rejection at p < 0.05 is
0.048–0.053 in all 36 null cells, and 0.0000 at p < 5e-8.

Its power to notice the misfit grows with `N_eff`, as it must, since it measures
the same quantity that does the damage (arm A, scale, OLS, p < 0.05):

| fit | N=30,000 | N=100,000 | N=461,320 |
|---|---|---|---|
| `ADD` | 1.0000 | 1.0000 | 1.0000 |
| `LS` | 0.2982 | 0.5161 | 0.5317 |
| `LSK` | 0.0640 | 0.0965 | 0.2701 |
| `TRU` | 0.0619 | 0.0755 | 0.1413 |

Two gaps, both worth acting on:

1. **Power is only ~53% where the damage is worst.** At arm A, `LS`,
   N = 461,320, GLS coverage of `beta_sigma2` is 0.8437 while Q flags just
   53% of replicates at p < 0.05 (50% at 5e-8). Roughly half of badly affected
   variants pass the guard.
2. **The GLS path emits no Q p-value at all.** Confirmed on a minimal
   reproducer: `param_gwas(..., tau_cov_estimator = "gls")` returns `Q` and `df`
   but no `Q_pval`, `Q_df_eff` or `Q_scale`, where the OLS call returns all of
   them. The guard `OLS_VS_GLS.md` asked for is unavailable on the estimator
   that needs it most. (The GLS `Q` values are also on a visibly different scale
   — 30–49 against OLS's 8–18 on identical near-null inputs at df 43 — so this
   needs its own calibration, not merely plumbing a p-value through.)

---

## 8. What this does not test

* **Only the pooled covariance path.** Every cell uses pooled `Sigma_0` with
  `nu = N`, as the earlier study did. The block-jackknife covariance path was
  not simulated. The section 6 numbers therefore describe what the jackknife df
  law *would* prescribe, not a measured jackknife run; applying t_nu to these
  pooled-path results is conservative (implied Type-I 0.30–0.84 x nominal),
  which is consistent with the package's own statement that the per-SNP
  jackknife df law does not apply to a pooled `Sigma_0`.
* **One covariate model** (intercept only). No LOCO, no principal components,
  no relatedness, no covariate adjustment.
* **One scale effect size** (`c = 0.05`) and one location size (`a = 0.05`).
  Because the damage is driven by lack of fit in SE units, a smaller `c` will
  move the onset to a larger `N_eff`, not remove it; that trade was not mapped.
* **Two error shapes** (normal, gamma(2)); gamma(8) was dropped from the ladder
  to buy replicates. `B = 1000` blocks throughout, so the `1/(B-1)` term was
  held fixed by design.
* **Independent genotypes and no LD**, single variant per replicate. The
  genome-wide figures are extrapolations from sd(z) under normality, not
  measured tail counts; section 2 gives their real resolution.
* **Coverage against each estimator's own population estimand.** GLS is scored
  against the population GLS target, so its degradation is a genuine failure of
  its own interval, not an artefact of comparing it to OLS's target.
* Arms B and C separate MAF from MAC only at the two points tested; the full
  (MAF, N) surface was not filled in.

---

## Reproducing

Workspace `/user/work/fh6520/nscale_run` (raw Stage 1 binaries, per-cell
`out/stage2_cell*.rds`, assembled `out/summary.tsv` and `out/pertau.tsv`).
Scripts in `sim/nscale/src/`: `cells.py` (design), `sim_stage1.py`,
`stage2_driver.R`, `truth.R`, `summarise.R`, `report.R`, `mechanism.R`,
`qguard.R`, `qcheck.R`, `crosscheck.R` (replication against `gls_coverage`),
`run_stage1.sh`, `run_stage2.sh`, `analyse.sh`.

Cost: 720,000 replicates; 1.46 s/replicate at N = 461,320 measured single-process
(3.9 s/replicate at 16-way contention) and 0.56 GB per process; about 470 CPU-hours of Stage 1 across 36 array tasks, Stage 2 a small
fraction of that.
