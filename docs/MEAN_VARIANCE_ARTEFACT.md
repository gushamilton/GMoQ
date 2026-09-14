# How a pure mean effect manufactures a vQTL, and what to do about it

Status: exact population calculation, 2026-09-14. Cross-validated against an
independent implementation. Working code in `sim_meanvar/src/`, package code in
`stage2/R/dispersion-guards.R`.

**Headline.** On a raw, untransformed, right-skewed trait a mean effect of
**|Z_mu| = 9.5 — about 0.02 SD per allele — is already enough to drive the
dispersion parameter past 5e-8 with no true dispersion effect whatsoever.** That is
an ordinary GWAS hit, not an outlier. Log or rank-inverse-normal transformation
removes this mechanism *exactly* and pushes the threshold out to |Z_mu| = 43-93.
A second, unrelated mechanism then takes over, and no transformation touches it.

---

## 1. There are two mechanisms, not one, and the earlier attribution was wrong

`docs/RELATEDNESS.md` reported, for a SNP with a pure mean effect and no variance
effect, at N = 200,000:

| scale | that document | exact, this analysis | true dispersion effect |
|---|---|---|---|
| raw (lognormal) | 28.1 SD | **33.2** | **0.2341 — real, not spurious** |
| log | 2.7 SD | **3.12** | 0 (exactly) |
| RINT | 3.2 SD | **3.09** | -1.7e-4 |

Those three numbers came from a 20-replicate Monte Carlo in
`/user/work/fh6520/scratch_rel/pc.R`, whose "bias in SD units" is `mean/sd` across
replicates, i.e. the per-replicate |Z|. The Monte Carlo error on such a ratio is
about `1/sqrt(2(n-1))` relative, so the published figures carry **+-4.6 on 28.1 and
+-0.49 on 2.7 and 3.2**. All three are therefore consistent with the exact values;
only the raw figure was materially off, and only because 20 replicates cannot pin it
down.

Two of that document's conclusions do not survive.

**"RINT does not remove it" is wrong.** RINT removes the mean-variance mechanism
completely. The clean test is the genotype skewness, `gskew = (1-2p)/sqrt(2p(1-p))`,
which is exactly zero at MAF 0.5:

| analysed scale | MAF 0.30 (gskew 0.617) | MAF 0.50 (gskew 0) |
|---|---|---|
| raw (lognormal) | 33.2 | **34.4** |
| log | 3.12 | **2.8e-13** |
| RINT | 3.09 | **6.4e-14** |

At MAF 0.50 the log and RINT cells are zero to machine precision. Every bit of the
residual that survives transformation is proportional to the *genotype* skewness,
and vanishes when that skewness does. Meanwhile the raw cell is undiminished at
MAF 0.50, because its driver is the *phenotype*, not the genotype.

**"The residual after transformation is W-misspecification" is wrong.** `W` is
exactly as misspecified at MAF 0.50 as at MAF 0.30 — the weight matrix does not know
the allele frequency. A mechanism that disappears at MAF 0.50 cannot be `W`
misspecification. It is the genotype-skewness residue identified independently in
`sim/dispersion` (`docs/FINDINGS_pooled_dispersion.md`), and their guard statistic
predicts it to within 0.4 per cent:

| cell | exact \|Z_sigma2\| | their G | ratio |
|---|---|---|---|
| log, MAF 0.30, Z_mu 57.3 | 3.122 | 3.110 | 1.004 |
| RINT, MAF 0.30, Z_mu 57.3 | 3.088 | 3.110 | 0.993 |

**The raw-scale number is not a bias at all.** On the raw lognormal scale the true
per-allele variance effect *is* 0.2341: a lognormal's conditional SD is proportional
to its conditional mean, so a mean effect really does move the variance. Stage 2
returns 0.1405, which **under**-states the true raw-scale variance effect by 40 per
cent. Calling this a "spurious variance effect" is a category error. It is a real
effect on a scale nobody wants to interpret, estimated with attenuation because the
Gaussian `W` is the wrong family for a lognormal.

### The two mechanisms, side by side

|   | (a) mean-variance coupling | (b) genotype skewness x density curvature |
|---|---|---|
| driver | skewness of the **phenotype** on the analysed scale | skewness of the **genotype**, `(1-2p)/sqrt(2p(1-p))` |
| scales as | **linear** in the mean effect (fitted exponent 1.10) | **quadratic** (fitted exponent 2.00) |
| at MAF 0.5 | unchanged | exactly zero |
| removed by log / RINT | **yes, exactly** | no |
| grows with N | no — the \|Z\| ratio is N-invariant | yes, as sqrt(N) at fixed per-allele effect |
| visible in the Q lack-of-fit | **yes, loudly** (ncp 465 on 43 df) | **no** (ncp 0.00016) |
| common slope across variants | **yes** (5% spread over MAF 0.01-0.5) | **no** (21x spread) |
| right tool | transform the trait, or fit the right `W`; `r_av` as fallback | flag with `G`; nothing removes it |

They are orthogonal. Correcting one does nothing to the other.

---

## 2. The number for the paper

Smallest **pure mean effect** whose spurious dispersion statistic reaches 5e-8, at
N = 461,320. The lognormal here has log-scale SD 0.5, i.e. raw-scale skewness 1.75 —
comparable to untransformed BMI or weight.

| analysed scale | MAF | effect, SD/allele | **\|Z_mu\|** | variance explained |
|---|---|---|---|---|
| **raw (lognormal)** | 0.05 | 0.044 | **9.27** | 0.019% |
| **raw (lognormal)** | 0.30 | 0.021 | **9.48** | 0.019% |
| **raw (lognormal)** | 0.50 | 0.020 | **9.53** | 0.020% |
| log | 0.05 | 0.206 | 43.0 | 0.40% |
| log | 0.30 | 0.213 | 93.3 | 1.89% |
| RINT | 0.05 | 0.206 | 43.1 | 0.40% |
| RINT | 0.30 | 0.214 | 93.9 | 1.91% |

Plainly: **on a raw skewed trait, a mean effect of |Z| = 9.5 manufactures a
genome-wide-significant dispersion hit.** Since 5e-8 corresponds to |Z_mu| = 5.45,
essentially *every* mean-effect locus that clears genome-wide significance by a
factor of two also clears it on the dispersion parameter. A raw-scale vQTL scan on a
skewed trait is, to first order, a re-run of the mean scan.

The raw-scale threshold is almost independent of N (|Z_mu| = 9.0-9.5 at N = 100,000,
200,000 and 461,320) because both statistics grow as sqrt(N) together. Transforming
is not optional.

After transformation, mechanism (b) sets the threshold, and it is MAF-dependent in
the opposite direction — rarer variants are worse:

| MAF | \|Z_sigma2\| at Z_mu 8 | 16 | 30 | 60 |
|---|---|---|---|---|
| 0.005 | 0.64 | 2.54 | 8.86 | 34.04 |
| 0.01 | 0.45 | 1.79 | 6.25 | 24.47 |
| 0.05 | 0.19 | 0.75 | 2.63 | 10.46 |
| 0.10 | 0.12 | 0.48 | 1.70 | 6.77 |
| 0.30 | 0.04 | 0.16 | 0.56 | 2.22 |
| 0.50 | 0.00 | 0.00 | 0.00 | 0.00 |

(N = 461,320, normal phenotype, identity scale; 5.451 is genome-wide.)

At a fixed |Z_mu| the artefact *shrinks* as 1/sqrt(N) — 1.67 at N = 50,000 down to
0.379 at N = 1,000,000 for MAF 0.30, Z_mu 30. At a fixed per-allele effect size it
*grows* as sqrt(N). Both statements are in `G = 0.6855 * gskew * Z_mu^2 / sqrt(N)`;
they differ only in what is held fixed, and the `sim/dispersion` phrasing ("more data
makes it worse") is the fixed-effect-size one.

### Phenotype shape modulates mechanism (b), and the guard does not know it

Pure location effect, MAF 0.30, Z_mu = 30, N = 461,320, **shipped** Gaussian `W`:

| residual law | skewness | excess kurtosis | \|Z_sigma2\| | G | observed / G |
|---|---|---|---|---|---|
| normal | 0.00 | 0.00 | 0.557 | 0.557 | **1.00** |
| chi2_16 | 0.71 | 0.75 | 0.797 | 0.565 | 1.41 |
| chi2_8 | 1.00 | 1.50 | 1.153 | 0.586 | 1.97 |
| chi2_4 | 1.41 | 3.00 | 2.371 | 0.688 | **3.45** |
| chi2_2 | 2.00 | 6.00 | 8.368 | 1.299 | **6.44** |
| t_8 | 0.00 | 1.50 | 0.584 | 0.610 | 0.96 |
| t_6 | 0.00 | 3.00 | 0.623 | 0.652 | 0.96 |

Kurtosis alone is harmless (the symmetric t rows sit at 0.96). **Phenotype skewness
is not**: at chi2_4 the artefact is 3.4x what `G` predicts.

`sim/dispersion` reports `G` *over*-predicting by 1.6x on the same chi2_4 baseline,
i.e. conservative — the opposite sign. **Both are right; the difference is the
estimator.** A full factorial over the three candidate causes settles it:

| chi2_4, MAF 0.30 | W = `z/(2sd)` (shipped) | W = `(q-m)/(2v)` |
|---|---|---|
| OLS | 0.81 | **0.60** |
| GLS | **3.45** | 2.83 |

(observed / `G`, essentially constant over Z_mu = 8, 16, 30; the normal control gives
0.9993-1.0012 across all 24 configurations.) Their headline column is
`both.pooled.**ols**.beta_sigma2` with the distribution-free `W`: ratio 0.60, i.e.
`G` over-predicts by 1/0.60 = **1.67x**, exactly their reported ~1.6x. Whether `W` is
built from the null or the alternative marginal makes no difference at all (agreement
to four decimals).

So the guard is conservative under OLS and **anti-conservative by a factor of 3.4
under GLS on a skewed phenotype**. That matters, because GLS is the recommended
estimator once model fit has been checked. On a normal phenotype every configuration
agrees with `G` to 0.1 per cent, so this is purely a phenotype-skewness effect.
**`G >= 0.2` should be treated as calibrated for OLS and for normal phenotypes only;
under GLS on a skewed trait, use `G >= 0.06`, or recalibrate.**

---

## 3. Model and estimator make no difference to mechanism (b)

MAF 0.30, N = 461,320. `LS` is `make_weights_vqtl`; `TRU` swaps in the
distribution-free scale column; `LSK` adds a `z^2 - 1` kurtosis column.

| scenario | model | est | \|Z_sigma2\| | bias | true | Q ncp (43 df) |
|---|---|---|---|---|---|---|
| pure location, normal | LS | ols | 0.556 | 0.00186 | 0 | — |
| pure location, normal | LS | gls | 0.557 | 0.00186 | 0 | **0.0002** |
| pure location, normal | TRU | gls | 0.557 | 0.00186 | 0 | 0.0002 |
| pure location, normal | LSK | gls | 0.557 | 0.00186 | 0 | 0.0000005 |
| raw lognormal | LS | gls | 50.4 | -0.0936 | 0.2341 | **464.8** |
| raw lognormal | TRU | gls | 54.8 | **+0.0113** | 0.2341 | **2.07** |
| raw lognormal | LSK | gls | 47.5 | -0.0200 | 0.2341 | 31.9 |

Two things follow.

**Nothing in the model tournament sees mechanism (b).** Every fit gives the same
0.557, and the lack-of-fit statistic is 0.0002 on 43 degrees of freedom. The residue
is *perfectly collinear* with the dispersion column, so by construction it leaves no
residual to detect. `Q` is not a guard against it, and must not be presented as one.

**The Q statistic sees mechanism (a), but only well after the damage is done.** At
`Z_mu = 83` the lack-of-fit is 464.8 on 43 df — astronomical, and Stage 2 would refuse
the model. Swapping in the distribution-free scale column (`TRU`) drops it to 2.07 on
43 df and recovers the true raw-scale variance effect to within 5 per cent (+0.0113 on
0.2341) instead of understating it by 40 per cent. So a correctly specified `W` is a
genuine fix, and `Q` genuinely identifies when one is needed.

**But `Q` does not warn in time.** At the effect size where the dispersion artefact
*first* reaches 5e-8, the lack-of-fit is small:

| trait skewness | \|Z_mu\| at 5e-8 | Q ncp (43 df) | **power of Q at 0.05** |
|---|---|---|---|
| 0.30 | 34.7 | 0.20 | **0.053** |
| 0.61 | 19.7 | 0.75 | **0.060** |
| 1.13 | 12.3 | 2.24 | **0.083** |
| 1.75 | 9.48 | 4.69 | **0.131** |

At the very effect size where the dispersion artefact first becomes genome-wide
significant, the lack-of-fit test rejects at 5 to 13 per cent — against a nominal 5
per cent. It is, to a good approximation, blind. **The model tournament is silent at
exactly the effect sizes that matter**, and only becomes loud once the artefact is
already 50-sigma. This is why a correction is needed and not merely tidy: neither `Q`
nor the model choice protects the threshold region, and mechanism (b) is invisible to
`Q` at any effect size (ncp 0.0002).

---

## 4. The fix, and what it is and is not good for

### 4.1 Young et al.'s residualisation is well specified for mechanism (a)

`d_l = alpha_vl - r_av * alpha_l` assumes a single genome-wide slope. On a lognormal
trait with a pure mean effect, the estimated ratio `theta_sigma2 / theta_mu` is:

| MAF | delta 0.01 | 0.03 | 0.10 |
|---|---|---|---|
| 0.01 | 0.7049 | 0.7235 | 0.7847 |
| 0.05 | 0.7047 | 0.7232 | 0.7867 |
| 0.20 | 0.7040 | 0.7220 | 0.7916 |
| 0.30 | 0.7035 | 0.7209 | 0.7922 |
| 0.50 | 0.7025 | 0.7180 | 0.7865 |

Mean 0.738, relative spread **5.1 per cent**, across MAF 0.01-0.50 and Z_mu 1.8-90.
A single slope is a good model. Note the ratio is essentially flat in MAF and drifts
only with effect size — the residual drift is what limits the correction at extreme
effects (§4.3).

`r_av` is a property of the trait, not a constant of nature:

| log-scale SD | raw-scale skewness | r_av (estimated) | r_av (true) |
|---|---|---|---|
| 0.20 | 0.61 | **0.140** | 0.140 |
| 0.35 | 1.13 | 0.377 | 0.464 |
| 0.50 | 1.75 | 0.712 | 1.079 |
| 0.70 | 2.89 | 1.341 | 2.707 |
| 0.90 | 4.75 | 2.321 | 6.270 |

At a mild skewness of 0.61 the implied `r_av` is **0.140**, against Young et al.'s
0.135 for UK Biobank BMI. That is a reassuring external check on the whole picture.
Note also that estimated and true `r_av` diverge as skewness grows — because
`theta_sigma2` is attenuated by the misspecified `W`. This does not matter for the
correction, which estimates `r_av` on the *estimated* scale and applies it on the
*estimated* scale, and is self-consistent there.

### 4.2 It is the wrong tool for mechanism (b), and cannot be made right

Same calculation for a pure location effect on a homoskedastic normal, where the
true `r_av` is zero:

| MAF | Z_mu 16: implied r_av | Z_mu 30: implied r_av |
|---|---|---|
| 0.01 | 0.164 | 0.307 |
| 0.05 | 0.069 | 0.129 |
| 0.30 | 0.0145 | 0.0273 |
| 0.50 | 0.0000 | 0.0000 |

A 21-fold spread, varying with both MAF and effect size. **No single genome-wide
slope can absorb this.** Young's correction is per-trait and linear; this artefact is
per-variant and quadratic. Saying so explicitly matters, because a reader who knows
Young et al. will assume the residualisation covers the whole problem. It does not.

### 4.3 What the residualisation actually achieves

`r_av = 0.699` estimated from a small-effect variant (standing in for an LD-pruned
genome-wide fit), then applied. Raw lognormal, MAF 0.30, N = 461,320:

| scenario | Z_mu | \|Z_sigma2\| before | after | SE ratio |
|---|---|---|---|---|
| pure mean, delta 0.03 | 25.1 | 14.58 | **0.62** | 0.702 |
| pure mean, delta 0.10 | 83.1 | 50.45 | **8.36** | 0.707 |
| pure dispersion, gamma 0.02 | 0.58 | 6.91 | **9.43** | 0.699 |
| pure dispersion, gamma 0.05 | 1.31 | 17.10 | **23.47** | 0.698 |
| both | 81.9 | 55.80 | 17.91 | 0.706 |

- At an ordinary large hit (Z_mu = 25) the artefact is removed completely: 14.6 to
  0.62, comfortably null.
- At an extreme hit (Z_mu = 83) it is reduced 6-fold but a residue of 8.4 survives.
  About half of that is mechanism (b) (`G` = 4.3 there) and about half is the 5 per
  cent drift in `r_av` with effect size. Both are covered by the `G >= 0.2` flag.
- **Genuine dispersion effects are not merely preserved, they gain power** — 6.91 to
  9.43, and 17.10 to 23.47. The baby is not thrown out; it is better fed.

The power gain is real and scale-specific. On the raw lognormal the Gaussian `W`
leaves `theta_mu` and `theta_sigma2` correlated at 0.735, so subtracting
`r_av * theta_mu` removes variance as well as bias and the SE falls 30 per cent.

On a transformed scale the correction is an exact no-op, not a cost:

| analysed scale | r_av implied | corr(mu, sigma2) | corrected SE / raw SE | \|Z\| of a true dispersion effect, before / after |
|---|---|---|---|---|
| raw | 0.699 | 0.735 | **0.699** | 6.91 / **9.43** |
| log | 0.0020 | 7.8e-16 | **1.0000** | 12.07 / 12.07 |
| RINT | 0.0040 | -2.9e-15 | **1.0000** | 12.07 / 12.07 |

On log and RINT the two parameters are orthogonal to machine precision and `r_av` is
zero, so the correction changes nothing whatsoever. It therefore never costs power on
any scale tested. The case for keeping it an extra column rests on the estimand and
on transportability (§5), not on power.

---

## 5. What is implemented

`stage2/R/dispersion-guards.R`, tested by `stage2/tests/test-dispersion-guards.R`.

- **`dispersion_skew_guard(z_mu, N, maf | n_eff)`** — `sim/dispersion`'s
  `G = 0.6855 * gskew * Z_mu^2 / sqrt(N)`, flagged at `G >= 0.2`. Stage 1 does not
  emit an allele frequency, so when `maf` is absent `gskew` is recovered from the
  per-variant `d_resid` that Stage 2 already threads through as `n_eff`, via
  `v = d_resid/N ~= 2p(1-p)` and `gskew = sqrt((1-2v)/v)`. That recovers the
  magnitude but **not the sign**, since `2p(1-p)` cannot separate `p` from `1-p`.
  This is why the guard flags rather than corrects — a threshold on |G| does not
  need the sign. Emitting an allele frequency from Stage 1 would remove the
  restriction and is the obvious follow-up.

- **`estimate_rav(...)`** — the genome-wide coupling, as a method of moments
  corrected for two contaminations an ordinary regression suffers:

  ```
  r_av = sum_l w_l (alpha_l * alpha_vl - c_l) / sum_l w_l (alpha_l^2 - s_al^2)
  ```

  `-c_l` removes the **within-variant sampling covariance** between the two parameter
  estimates, which Stage 2 already reports in `param_cov_upper` and which would
  otherwise masquerade as trait-level coupling. `-s_al^2` removes **regression
  dilution**, which for a scan dominated by null variants is severe. Neither
  correction appears in the original method; both are free here. A delete-one-block
  jackknife gives the SE. Estimate it on an LD-pruned genome-wide set, not on the
  hits — the denominator correction is what makes including null variants safe.

  Measured behaviour, on 50,000 synthetic variants with `r_av` = 0.135 and sampling
  error comparable to the signal (`sim_meanvar/src` trim check):

  | | slope | bias / SE |
  |---|---|---|
  | naive OLS of `alpha_v` on `alpha` | 0.1267 | biased **low** (dilution) |
  | dilution correction only | 0.1997 | biased **high** by 48% (the within-variant covariance) |
  | **both corrections** | **0.1329** | **-0.14** |

  Neither correction alone is enough, and they err in opposite directions. Trimming
  is mild at the default: `trim` = 0, 0.001, 0.01, 0.05 give bias/SE of -0.14, -0.09,
  +0.21, +1.45, so 0.01 is safe and 0.05 is already too aggressive (it removes 40 per
  cent of the denominator mass).

  **The assumption, and how it fails.** Consistency requires genuine dispersion
  effects to be uncorrelated with mean effects — Young et al.'s assumption, not
  testable from these data. With genuine vQTLs present but independent, `r_av` is
  recovered cleanly (0.1328 against 0.135, -1.16 SE). With a real dispersion signal
  genuinely proportional to the mean effect (true total coupling 0.235), the
  estimator returns 0.238 — it **absorbs the real signal into the nuisance
  parameter**. That is the correct behaviour given the assumption, and it is also the
  reason the uncorrected column must stay: if the biology really is "variance scales
  with mean", the residualisation defines that away by construction.

- **`residualise_dispersion(...)`** — the corrected estimate and its SE:

  ```
  Var(d_l) = s_vl^2 + r_av^2 s_al^2 - 2 r_av c_l + Var(r_av) alpha_l^2
  ```

  The fourth term propagates uncertainty in `r_av` itself. It is negligible for null
  variants and matters precisely for the large-|alpha| variants the correction exists
  to adjudicate, so it is on by default; `se_rav = 0` reproduces the original.

- **`dispersion_report(fit, ...)`** — returns uncorrected and corrected statistics
  **side by side**, plus the guard.

### Design decisions, and why

**Extra column, never a silent replacement.** Two reasons, and note that "it costs
power" is *not* one of them — §4.3 shows it costs nothing on any scale tested. First,
the corrected statistic is a *different estimand*: "dispersion beyond what the mean
effect predicts" is a different scientific claim from "dispersion", and silently
swapping it in would misrepresent output that has already been published and
interpreted. Second, `r_av` is a genome-wide nuisance parameter estimated from the
same data, so the corrected statistic does not transport across cohorts, traits or
analysis scales the way the raw one does; a meta-analysis needs the uncorrected
column and each cohort's `r_av` separately.

**Necessary, not ornamental.** The obvious objection is that a transformation or a
correctly specified `W` already solves mechanism (a), so why add anything. §3 answers
it: at the effect size where the artefact first crosses 5e-8, the lack-of-fit
statistic is 0.2 to 9.7 on 43 df — invisible. Nothing in the pipeline tells the
analyst that the fitted family is wrong until the artefact is already enormous. If a
raw-scale analysis is run at all, the residualisation is the only thing standing
between it and a page of false vQTLs. The first-line advice remains to transform;
this is what to do when that advice is not taken, and it is also what makes the
results comparable with Young et al.

**Interaction with the pooled covariance.** `c_l` comes from `param_cov_upper`, which
inherits whichever covariance source is in force. Under the pooled estimator that
covariance carries the working-independence error `sim/dispersion` measured
(`se_ratio < 1` under a true location effect). The residualisation neither fixes nor
worsens it; it passes it through. Nothing here is a substitute for getting the
covariance right.

**Interaction with `tau_cov_estimator = "gls"`.** Both work. But `r_av` **must be
estimated with the same estimator and covariance source used for the scan**, because
the attenuation of `theta_sigma2` is estimator-dependent (0.1405 under LS/GLS versus
0.2341 true, and different again under OLS and under `TRU`). Mixing an `ols`-derived
`r_av` into a `gls` scan is a silent error. This is a real operational constraint.

**GLS's within-variant orthogonalisation is not this.** Stage 2's GLS makes
`theta_mu` and `theta_sigma2` orthogonal in the *sampling* sense within a variant.
`r_av` is their *population* relationship across variants. The tables in §4.1 are
computed under GLS and the coupling is fully present there, so the two are plainly
different things.

**Does the joint `both` fit already absorb it?** No. Every number in §1-4 *is* the
joint location-and-scale fit. Fitting the mean jointly removes the component of
`b(tau)` that is flat in tau; both artefacts are z_tau-shaped, which is exactly what
the scale column is built to pick up. §3 shows LS, LSK and TRU give identical answers
for mechanism (b). Joint fitting is not a fix for either mechanism.

---

## 6. Monte Carlo through the shipped pipeline

Shipped Stage 1 (`core.run_stage1`, `compute_rif`, `pooled_cov.compute_sigma0`) and
shipped Stage 2 (`param_gwas`, `transform = "custom_W"`, `se_mode = "tau_cov"`,
pooled covariance). N = 50,000; 400 replicates; 200 independent variants per
replicate against one phenotype, because the correction is genome-wide and must be
estimated from many variants. 170 null, 20 pure-mean (|Z_mu| 4 to 40, half each
sign), 10 pure-dispersion. `r_av` estimated from **all 200**, with no selection on
truth. That is 68,000 null, 8,000 mean and 4,000 dispersion variant-observations per
cell. Code in `sim_meanvar/src/mv_*.{py,R}`; tables in
`/user/work/fh6520/meanvar/mc/summary/`.

The Stage 2 call was verified against a hand-rolled fit: OLS `theta` max abs error
**exactly 0**, OLS SE 2.2e-16, GLS `theta` 1.7e-12, `cov_vec` round-trip 1.7e-18.

### (a) the artefact is removed, and (b) is not

Rejection of `beta_sigma2` at 0.05 on the **pure-mean** class, where the true
dispersion effect is exactly zero:

| analysed scale | est | uncorrected | residualised |
|---|---|---|---|
| raw | ols | **7545/8000 = 0.943** +-0.003 | **3497/8000 = 0.437** +-0.006 |
| raw | gls | **7723/8000 = 0.965** +-0.002 | **2850/8000 = 0.356** +-0.005 |
| log | gls | 3313/8000 = 0.414 +-0.006 | 3478/8000 = 0.435 +-0.006 |
| rint | gls | 3319/8000 = 0.415 +-0.006 | 3442/8000 = 0.430 +-0.006 |

The residualisation takes `raw` from 0.943 down to 0.437 — **exactly onto the
log/RINT floor of 0.41-0.42**, which is where a scale with no mean-variance coupling
already sits. It removes essentially all of mechanism (a) and, as designed, nothing
at all of mechanism (b). On log and RINT it is a no-op to within Monte Carlo error.

Per effect size (`raw`, OLS), mean |z| on a variant with zero dispersion effect:

| designed Z_mu | 4 | 12 | 20 | 28 | 40 |
|---|---|---|---|---|---|
| uncorrected z | 1.90 | 6.02 | 10.70 | 15.44 | 23.79 |
| residualised z | -0.02 | 0.21 | 1.09 | 2.09 | 4.96 |

(The floor is large here because N = 50,000 with |Z_mu| up to 40 puts `G` as high as
21 for the rarer variants in the mixture — mechanism (b) is deliberately extreme in
this design.)

### Type-I error and the genome-wide extrapolation

On genuinely null variants, corrected and uncorrected alike:

| cell | est | uncorrected k/n | rate | corrected rate | sd(z) null | implied rate @5e-8 |
|---|---|---|---|---|---|---|
| log | ols | 3368/68000 | 0.0495 +-0.0008 | 0.0492 | 0.9963 | 4.46e-08 |
| log | gls | 3319/68000 | 0.0488 +-0.0008 | 0.0486 | 0.9955 | 4.35e-08 |
| raw | ols | 3375/68000 | 0.0496 +-0.0008 | 0.0509 | 0.9994 | 4.91e-08 |
| raw | gls | 3413/68000 | 0.0502 +-0.0008 | 0.0505 | 1.0016 | 5.25e-08 |
| rint | ols | 3413/68000 | 0.0502 +-0.0008 | 0.0502 | 0.9998 | 4.72e-08 |
| rint | gls | 3419/68000 | 0.0503 +-0.0008 | 0.0500 | 0.9980 | 4.71e-08 |

The 5e-8 column uses `2*pnorm(-qnorm(1-5e-8/2)/sd_z)`, the convention in
`gls_coverage/src/summarise.R`; 5e-8 is unreachable by simulation directly. Every
cell lands within 0.87x to 1.05x of nominal. **The null calibration of the dispersion
parameter is sound. The artefact is entirely a mean-effect problem.**

### The baby is not thrown out

Pure-dispersion class. "Sign-aligned bias/SE" is `mean(sign(truth)*(est-truth)/se)`;
negative is shrinkage toward zero. MCSE 0.016.

| cell | est | bias/SE | power uncorrected | power residualised |
|---|---|---|---|---|
| log | gls | -0.13 | 0.9670 +-0.0028 | 0.9673 +-0.0028 |
| raw | ols | -0.07 | 0.9637 +-0.0030 | **0.9802** +-0.0022 |
| raw | gls | **-0.28** | 0.8340 +-0.0059 | **0.9490** +-0.0035 |
| rint | gls | -0.08 | 0.9647 +-0.0029 | 0.9650 +-0.0029 |

The correction **raises** power on the raw scale (0.834 to 0.949) and leaves it
untouched on log and RINT. The estimator is unbiased on log and RINT; on raw under
GLS there is a genuine 0.28 SE shrinkage that grows with effect size, which is `W`
misspecification for a lognormal, not the correction.

### Two real defects the Monte Carlo exposed

**Coverage of the residualised statistic is worse, not better, on the raw scale.**
95% Wald coverage about the true value: null class is nominal everywhere
(0.9491-0.9514), but the `raw` mean class falls to **0.9173** residualised against
0.9324 uncorrected, and the `raw` dispersion class under GLS to **0.9165** against
0.9343. The headline numbers above use `se_rav = 0`, so `r_av`'s own uncertainty is
missing from the interval. Anyone using the residualised column for *estimation*
rather than *testing* must pass `se_rav`.

**The block-jackknife SE for `r_av` is miscalibrated on the raw scale.** Against the
across-replicate SD of `r_av`:

| cell | est | r_av mean | across-replicate sd | mean jackknife SE |
|---|---|---|---|---|
| log | gls | 0.00005 | 0.0378 | 0.0393 |
| rint | gls | 0.0039 | 0.0562 | 0.0696 |
| raw | ols | 1.520 | 0.276 | **0.561** (2x too large) |
| raw | gls | 0.968 | 0.145 | **0.053** (3x too small) |

Well calibrated where the coupling is zero, badly calibrated where it is not. Passing
that SE as `se_rav` over-corrects `raw` OLS into conservatism (mean-class rate 0.437
to 0.0103, but the null rate also drops to 0.0405). The likely cause is that this
design jackknifes 20 blocks of 10 variants inside a single replicate, which is not
what a genome-wide LD-block jackknife over ~1e5 pruned variants does — but that is a
hypothesis, and **the jackknife SE should not be relied on until it is checked at a
realistic number of variants.**

### A design decision confirmed empirically

`r_av` on the raw scale is **1.520 under OLS and 0.968 under GLS** — genuinely
different, because the Gaussian `W` is misspecified for a lognormal and the two
projections of the same population `b(tau)` are different estimands. The GLS value
matches the analytic prediction `2*E[Y]*(exp(s^2)-1)` = 1.13 far better than the OLS
one. This confirms, from simulation, the constraint asserted a priori in §5:
**`r_av` must be estimated with the same estimator and covariance source used for the
scan.** Mixing them is a silent error of 50 per cent here.

---

## 7. What I could not settle

- **A third mechanism exists, and neither study has named it.** On a *skewed*
  phenotype there is a residue at MAF 0.50, where the genotype skewness is exactly
  zero, so it is neither mechanism (a) nor (b):

  | residual law | phenotype skew | \|Z_sigma2\| at MAF 0.50, Z_mu 30 | Q ncp (43 df) |
  |---|---|---|---|
  | normal | 0.00 | 1.3e-13 | 0.0004 |
  | chi2_16 | 0.71 | 0.013 | 0.0013 |
  | chi2_8 | 1.00 | 0.016 | 0.0026 |
  | chi2_4 | 1.41 | 0.167 | 0.013 |
  | chi2_2 | 2.00 | **26.5** | **1725** |

  It is **cubic** in the effect size (fitted log-log slope 3.02 at chi2_4, under both
  weight matrices), against linear for (a) and quadratic for (b). It is negligible
  except at extreme phenotype skewness, and where it is large the lack-of-fit
  statistic is enormous (1725 on 43 df), so `Q` does catch that regime. I have not
  derived it analytically and it is not guarded. It is a third-order term and I do
  not think it is a practical threat, but that is a judgement, not a proof.
- **The sign of `gskew` from `d_resid`.** Recoverable only if Stage 1 emits an allele
  frequency. Currently the guard is magnitude-only.
- **Recalibrating `G` for GLS on skewed phenotypes.** §2 shows the factor is about
  3.4 at chi2_4 and 6.4 at chi2_2 with the shipped `W`, but I have not fitted a
  replacement constant across the MAF x N x skewness grid the way `sim/dispersion`
  fitted 0.6855 for OLS on the normal.
- **The jackknife SE for `r_av` at a realistic number of variants.** §6 shows it is
  2x too large under OLS and 3x too small under GLS on the raw scale, in a design
  with only 200 variants and 20 blocks. Until that is repeated with ~1e5 LD-pruned
  variants, treat `estimate_rav()$se` as indicative.
- **Real genotypes and real traits.** No scan was run, per the brief. Everything is
  simulated with independent variants; LD, genotyping error and Hardy-Weinberg
  departure are untested, and `r_av` for the project's actual `bmi` and `weight`
  phenotypes has not been estimated.
- **Real genotypes, LD, Hardy-Weinberg departure, genotyping error.** Everything here
  uses an exact 3-point HWE genotype law.
- **Heavy-tailed phenotypes beyond t_6.** Moments may not exist; not attempted.
- **Rates below the Monte Carlo floor.** The exact calculations carry no Monte Carlo
  error, but converting a non-centrality to a rate at 5e-8 assumes `theta_hat` is
  normal about `theta_pop` that far into the tail. That assumption is inherited from
  `sim/dispersion` and is supported by their Monte Carlo only down to 1e-3.
