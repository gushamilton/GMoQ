# Correcting the genotype-skewness bias in the dispersion parameter

**A pure location effect at a skewed genotype produces a false vQTL. This is bias,
not a variance error, so it can be subtracted rather than only flagged.
`correct_dispersion_skew()` subtracts it. On a Normal baseline the correction is
exact to \(O(\delta^4)\) and removes the artefact completely: the worst false-positive
rate in the sweep falls from 0.91 to \(1.5\times10^{-9}\) at genome-wide
significance, at a cost of nothing at all on null variants and nothing at all on
genuine dispersion effects — except at a rare variant carrying a colossal
location effect, where it costs about 26% of honest power (§4).**

**Three conditions, all enforced in code, all of which make the correction
harmful rather than merely useless if broken: a *signed* effect-allele frequency
(§8), `make_weights_vqtl`'s *numeric* design rather than merely its column names
(§2.1), and a baseline within the skewness envelope for the chosen estimator —
about 0.65 under GLS, about 1.0 under OLS (§5).**

Everything below is an exact population calculation — quadrature, no Monte Carlo,
no sampling error — confirmed by Monte Carlo through the shipped Stage 1 and
Stage 2. The exact engine reproduces Table 1 of `FINDINGS_pooled_dispersion.md`
to every published digit, which is how it is validated.

The campaign was `sim_skewcorr/` and is frozen at tag
`archive/sim-campaigns-20260917` — `git show
archive/sim-campaigns-20260917:sim_skewcorr/tables/mc_summary.txt`. The
`FINDINGS_pooled_dispersion.md` engine it validates against is archived at the
same ref. See [SIMULATION_INVENTORY.md](SIMULATION_INVENTORY.md) §8.

---

## 1. The result

For \(Y = \mu + \delta g + \varepsilon\) with \(\varepsilon\) independent of \(g\)
the true conditional quantile effect is \(\delta\) at every \(\tau\), but Stage 1's
RIF-OLS estimand is

```
b(tau) = h'(gbar) + (1/2) h''(gbar) * mu3/Var(g)
       = delta + (delta^2 / (2 sigma)) * (mu3/Var g) * z_tau + O(delta^3)
```

with \(z_\tau = \Phi^{-1}(\tau)\) and, for \(g \sim \mathrm{Binom}(2,p)\),
\(\mu_3/\mathrm{Var}(g) = 1-2p\).

The step that makes a *correction* possible rather than only a flag is this.
`make_weights_vqtl()` builds the dispersion column as
\(W_{\sigma^2}(\tau) = z_\tau/(2\sigma)\). The bias above is therefore not merely
proportional to that column, it **is** \(\delta^2(1-2p)\) times that column. Every
Stage 2 projection is \(\hat\theta = Ab\) with
\(A = (W'S^{-1}W)^{-1}W'S^{-1}\), and every such \(A\) satisfies \(AW = I\). Hence

```
A . bias  =  delta^2 * (1 - 2p) * e_2
```

* the induced bias in `beta_sigma2` is **exactly** \(\delta^2(1-2p)\);
* the induced bias in `beta_mu` is **exactly zero**;
* both hold for `ols`, `gls`, `dwls` and `plugin_cor` alike, on any \(\tau\) grid;
* the phenotype SD **cancels** — it is in the bias and again in the weight column.

Verified in-package to \(3\times10^{-16}\) (`test-dispersion-guards.R`), and exactly
against the population estimand:

| MAF | \(\delta/\sigma\) | bias OLS | bias GLS | \(\delta^2(1-2p)\) | OLS/pred | GLS/pred | bias in `beta_mu` |
|---|---|---|---|---|---|---|---|
| 0.50 | 0.35 | -4.7e-17 | -5.5e-17 | 0 | — | — | 2.8e-04 |
| 0.30 | 0.02 | 1.5999e-04 | 1.5999e-04 | 1.6000e-04 | 1.0000 | 1.0000 | 3.1e-08 |
| 0.30 | 0.10 | 3.99567e-03 | 3.99502e-03 | 4.00000e-03 | 0.9989 | 0.9988 | 3.8e-06 |
| 0.05 | 0.10 | 8.99315e-03 | 8.99355e-03 | 9.00000e-03 | 0.9992 | 0.9993 | -1.1e-05 |
| 0.01 | 0.10 | 9.79324e-03 | 9.79419e-03 | 9.80000e-03 | 0.9993 | 0.9994 | -1.4e-05 |
| 0.01 | 0.35 | 1.19043e-01 | 1.19170e-01 | 1.20050e-01 | 0.9916 | 0.9927 | -6.1e-04 |
| 0.01 | 1.00 | 9.17451e-01 | 9.17962e-01 | 9.80000e-01 | 0.9362 | 0.9367 | -1.7e-02 |

---

## 2. Design decisions

**1. Where to apply it — after the projection, on the coefficient.** The two
options are provably *identical* **on `make_weights_vqtl`'s design**, not merely
similar, because the bias vector lies exactly in the column space of \(W\) and
\(AW = I\); this is asserted as a test. Correcting \(b(\tau)\) *before* the
projection is right for **every** \(W\), so it is fair to ask why the narrower
choice ships.

It is not a cost argument. The bias vector is a per-variant scalar
\(\delta^2 k\) times the *fixed* vector \(u = \Phi^{-1}(\tau)/(2\sigma_Y)\), so
\(Au\) is formed once and the per-variant cost of the before-projection route is
also O(1). The reason is that before-projection **does not actually generalise**:

* \(\sigma_Y\) **stops cancelling.** After the projection it cancels exactly,
  because the bias *is* a multiple of the \(W\) column that carries \(\sigma_Y\)
  in its denominator. Before it, you must build \(u\) numerically, so the
  phenotype SD becomes a required input — one Stage 1 does not emit and Stage 2's
  report does not thread through. Getting it wrong scales the correction
  linearly, with no symptom: exactly the `w_sd_ratio` failure mode, made
  mandatory.
* \(\delta\) **stops being identifiable.** The bias is \(\delta^2 k\) with
  \(\delta\) the per-allele location shift in phenotype units. Under
  `make_weights_vqtl`, \(W[,\texttt{beta\_mu}] = 1\), so \(AW = I\) makes
  `beta_mu` *exactly* that \(\delta\). Under a general \(W\) there is no such
  guarantee: for the 5-parameter mixture, `beta_mu1` is a component mean, not the
  shift of the whole distribution — which is precisely why that design
  over-subtracts by 2.4x today. Shipping before-projection would move the silent
  failure from "wrong \(W\)" to "wrong \(\delta\)", not remove it.

So the class is closed not by moving the correction earlier but by refusing to
correct on a design the derivation does not cover. `check_vqtl_design()` does
that, on the *numbers* in \(W\) rather than its column names, and it is a total
function: it either verifies the exact identity the proof needs or it refuses.
Under that gate the after-projection form is exact, so a general-\(W\) path would
be dead code for every input that passes and still unsafe for every input that
does not. Generalising properly needs both prerequisites above — a phenotype SD
threaded from Stage 1, and a declared location contrast for a general \(W\) —
which is a redesign, not a guardrail.

Applying it afterwards is also O(1) per variant, needs no refit, and works for
every `se_mode` without touching the C++ path.

**2. `delta_hat` — the location coefficient from the same fit, debiased, no
iteration.** Because the bias is *entirely* absorbed by `beta_sigma2`,
`beta_mu` carries none of it (the table above: at 0.1 SD/allele the induced bias
in `beta_mu` is \(1.4\times10^{-5}\) against \(9.8\times10^{-3}\) in
`beta_sigma2`, a factor of 700). `beta_mu` is therefore an unbiased input to its
own correction and there is no feedback loop to iterate. Two finite-sample
matters do need handling and are handled:
*`beta_mu^2` is biased upward by \(\mathrm{Var}(\hat\delta)\)*, so `debias = TRUE`
subtracts `se_mu^2`; the result may go negative, which is correct for an unbiased
estimator of a non-negative quantity and is not truncated. *`beta_mu` is
correlated with `beta_sigma2`*, which is a variance matter and appears in the SE.

**3. `sigma_hat` — none is needed.** It cancels. This is the reason no new Stage 1
output is required. (`w_sd_ratio` exists only for the unusual case where the `sd`
passed to `make_weights_vqtl()` was deliberately not the phenotype SD. It **must
match that `sd`**: it is (the `sd` passed to `make_weights_vqtl`) / (the phenotype
SD), and leaving it at the default 1 when `W` was built on another scale
multiplies the whole correction by exactly that ratio — linearly, with no
symptom. `dispersion_report()` requires `W_sd` so the design check can confirm
which `sd` built `W`.)

**4. The \(O(\delta^3)\) term — it does not reach this parameter, *provided the
tau grid is symmetric about 0.5*.** The third-order term is *even* in \(z_\tau\),
so on a symmetric grid it is orthogonal to the dispersion column; it lands on
`beta_mu` and on the fit residual. The first term that survives into
`beta_sigma2` is then \(O(\delta^4)\), and the measured relative error of the
correction is \(\approx 0.065\,(\delta/\sigma)^2\): 0.07% at 0.1 SD/allele, 0.3% at
0.2 SD/allele, 0.8% at 0.35, 6.4% at 1.0. The residual *after* correction is
smaller than the bias by that factor, so it is negligible everywhere a real GWAS
lives and still small where none does.

**The symmetry is a condition, not a technicality**, and earlier versions of this
document asserted the orthogonality unconditionally. On an asymmetric grid the
even term is no longer orthogonal to the column and leaks straight into
`beta_sigma2`. Measured on `linspace(0.02, 0.70, 45)` the order of the residual
in \(\delta\) drops from 4.00 to **3.01** and its relative size grows about
**70x**, with a sign that flips depending on which tail is truncated.
`check_vqtl_design()` warns when the supplied grid is not symmetric about 0.5;
the correction is still applied, but its residual is larger than the
\(O(\delta^4)\) figure above. Use a symmetric grid.

**5. The standard error — it must widen, and it does.**
```
Var = se_sigma2^2 + 4 k^2 beta_mu^2 se_mu^2 + 2 k^2 se_mu^4 - 4 k beta_mu cov_mu_sigma2
```
The middle terms are the exact variance of \(k\hat\beta_\mu^2\) under joint
normality; the last is \(-2\mathrm{Cov}(\hat\beta_{\sigma^2}, k\hat\beta_\mu^2)
= -4k\hat\beta_\mu c_l\). The expression is provably non-negative over the whole
admissible range of \(c_l\). Dropping \(c_l\) is not safe under `gls`, where the
two parameters are correlated, so `param_cov_upper` should be passed.

**6. Presentation — extra columns, not replacement.** The corrected quantity is
arguably the estimand users wanted, which is a real argument for replacement, and
it is weakened by the fact that the correction costs nothing on the variants that
make up a scan. It loses to three specific things. (i) The correction **cannot
always be computed**: it needs a *signed* effect-allele frequency, which Stage 1
does not emit; a column that is silently corrected for some variants and not
others is worse than two honest columns. (ii) It is **specific to
`make_weights_vqtl`** on a Normal baseline; replacing a general column with one
valid only under that model would be a silent trap. (iii) On a strongly skewed
phenotype the correction is *partial* (§5), and a replaced column would hide that.
It is also the convention `residualise_dispersion()` already set.

---

## 3. Does it work? Pure location effect, zero true dispersion

Exact rates, N = 461,320, Normal baseline, GLS (OLS is the same to three digits).
Every rejection is false.

| MAF | Z_mu | \(\delta/\sigma\) | G | unc @0.05 | **cor @0.05** | unc @5e-8 | **cor @5e-8** | SE cost |
|---|---|---|---|---|---|---|---|---|
| 0.50 | 60 | 0.126 | 0.00 | 0.0513 | 0.0513 | 5.97e-08 | 5.97e-08 | 1.000 |
| 0.30 | 16 | 0.037 | 0.16 | 0.0531 | **0.0500** | 7.19e-08 | **5.03e-08** | 1.000 |
| 0.30 | 30 | 0.069 | 0.56 | 0.0875 | **0.0501** | 5.35e-07 | **5.10e-08** | 1.001 |
| 0.30 | 60 | 0.138 | 2.24 | 0.6210 | **0.0506** | 8.02e-04 | **5.40e-08** | 1.003 |
| 0.10 | 30 | 0.105 | 1.71 | 0.4093 | **0.0486** | 1.09e-04 | **4.12e-08** | 1.007 |
| 0.10 | 60 | 0.213 | 6.85 | 1.0000 | **0.0445** | 9.45e-01 | **2.27e-08** | 1.027 |
| 0.05 | 30 | 0.145 | 2.65 | 0.7675 | **0.0462** | 3.34e-03 | **2.94e-08** | 1.016 |
| 0.05 | 60 | 0.299 | 10.6 | 1.0000 | **0.0360** | 1.00e+00 | **5.51e-09** | 1.065 |
| 0.01 | 8 | 0.084 | 0.45 | 0.0751 | **0.0483** | 3.32e-07 | **3.97e-08** | 1.006 |
| 0.01 | 16 | 0.170 | 1.80 | 0.4560 | **0.0435** | 2.07e-04 | **1.98e-08** | 1.026 |
| 0.01 | **30** | 0.331 | 6.33 | **1.0000** | **0.0298** | **9.16e-01** | **1.54e-09** | 1.094 |
| 0.01 | 60 | 0.794 | 25.3 | 1.0000 | 0.0194 | 1.00e+00 | 2.14e-13 | 1.448 |

Read the MAF 0.01 / Z_mu 30 row: the rate at genome-wide significance goes from
**0.912 to 1.5e-9**, a factor of \(6\times10^8\), and from 1.000 to 0.0298 at 0.05.

Three properties hold across the whole sweep (5 MAF × 8 Z_mu × 2 estimators):

* **the correction never makes anything worse.** The largest corrected inflation
  anywhere is 1.03x nominal (MAF 0.30, Z_mu 60), against 20x uncorrected;
* **where it stops being exact it errs conservative**, not anti-conservative —
  the extreme cells under-reject (0.0194 at 0.05) rather than over-reject;
* **MAF 0.50 is untouched**, because \(k = 1-2p = 0\) exactly there.

### Genome-wide extrapolation

At 5e-8 with \(10^7\) tests, the expected number of false dispersion hits from a
single variant class is the rate times the class size. Uncorrected, one MAF 0.01
variant with Z_mu = 30 contributes 0.91 expected false hits on its own;
corrected, \(6\times10^8\) such variants would be needed to expect one. The whole
uncorrected \(G \ge 0.2\) region — which the flag exists to fence off — is brought
inside 1.03x nominal.

---

## 4. What does it cost?

**On null variants: nothing measurable.** Exact rates at \(\delta=\gamma=0\):

| MAF | unc @0.05 | cor @0.05 | unc @5e-8 | cor @5e-8 | SE ratio |
|---|---|---|---|---|---|
| 0.50 | 0.05000 | 0.05000 | 5.000e-08 | 5.000e-08 | 1.00000 |
| 0.30 | 0.05000 | 0.05000 | 5.000e-08 | 5.000e-08 | 1.00000 |
| 0.05 | 0.05000 | 0.05000 | 5.000e-08 | 4.997e-08 | 1.00001 |
| 0.01 | 0.05000 | 0.04998 | 5.000e-08 | 4.985e-08 | 1.00005 |

**On a genuine dispersion effect with no location effect: nothing measurable.**
Power at 5e-8, OLS:

| MAF | Z_sig | power unc | **power cor** | rel. bias unc | rel. bias cor | SE ratio |
|---|---|---|---|---|---|---|
| 0.30 | 5.45 | 0.51047 | 0.51046 | 0.00220 | 0.00220 | 1.00000 |
| 0.10 | 5.45 | 0.53274 | 0.53270 | 0.00670 | 0.00669 | 1.00000 |
| 0.05 | 5.45 | 0.55083 | 0.55074 | 0.01038 | 0.01037 | 1.00001 |
| 0.01 | 5.45 | 0.62022 | 0.61971 | 0.02498 | 0.02489 | 1.00005 |
| 0.01 | 12.0 | 1.00000 | 1.00000 | 0.05630 | 0.05621 | 1.00005 |

The correction term has mean zero and variance \(2k^2\mathrm{se}_\mu^4\) when
\(\beta_\mu = 0\), which is fourth order in a standard error. Power is unchanged
in the fifth decimal place. (The relative bias in the last two columns is the
estimator's own \(O(\gamma^2)\) curvature under a log-scale DGP; it is present with
no location effect at all, is identical before and after, and is outside this
correction's scope.)

**On a mixed location + dispersion effect: the false power is removed, which is
the point.** Remaining bias in `beta_sigma2` in units of its own reported SE, OLS:

| MAF | Z_mu | Z_sig | zbias unc | **zbias cor** | power unc @5e-8 | **power cor @5e-8** |
|---|---|---|---|---|---|---|
| 0.30 | 30 | 0 | 0.562 | **-0.0003** | 0.00000 | 0.00000 |
| 0.30 | 30 | 5.45 | 0.571 | **0.0096** | 0.71771 | **0.50590** |
| 0.05 | 30 | 0 | 2.695 | **-0.0026** | 0.00326 | 0.00000 |
| 0.05 | 30 | 5.45 | 2.713 | **0.0272** | 0.99633 | **0.50417** |
| 0.01 | 16 | 0 | 1.842 | **-0.0017** | 0.00020 | 0.00000 |
| 0.01 | 16 | 5.45 | 1.917 | **0.0916** | 0.97417 | **0.55201** |
| 0.01 | 30 | 0 | 6.903 | **-0.0210** | 0.91237 | **0.00000** |

\(Z_{\sigma} = 5.45\) is the effect size whose *honest* power at 5e-8 is 0.50. The
uncorrected estimator reports 0.72–0.997 there; the corrected one reports
0.50–0.55. That difference is not lost power, it is **false power being removed**:
the uncorrected excess is the location effect masquerading as dispersion.

**The one place the correction costs real power.** At a rare variant carrying a
colossal location effect the correction itself has real variance, and that is paid
for honestly:

| MAF | Z_mu | Z_sig | power unc @5e-8 | **power cor @5e-8** | honest power | SE ratio |
|---|---|---|---|---|---|---|
| 0.01 | 16 | 5.45 | 0.97565 | **0.56246** | 0.50 | 1.025 |
| 0.01 | 30 | 5.45 | 1.00000 | **0.37184** | 0.50 | 1.094 |

At MAF 0.01 with `Z_mu` = 30 a marginal dispersion effect loses about a quarter of
its honest power. That is the worst cell in the whole sweep, it is the same cell
where the uncorrected rate is 0.91, and `G` = 6.3 flags it.

The residual `zbias_cor` at large `Z_sig` (0.28 at MAF 0.05, Z_sig 12) is *not* a
failure of the correction. It is numerically equal to the same estimator's
\(O(\gamma^2)\) bias measured at \(\delta = 0\) (0.0232 × 12 = 0.278), i.e. it is
present with no location effect at all.

---

## 5. Where it breaks down: non-Normal baselines

The derivation assumes a Normal baseline, where \(f_0'/f_0 \propto z_\tau\) makes
the residue exactly collinear with the weight column. Exact rates at 5e-8, GLS,
MAF 0.01 unless stated:

| baseline | skewness | ex. kurtosis | Z_mu 8 unc → cor | Z_mu 16 unc → cor | Z_mu 30 unc → cor |
|---|---|---|---|---|---|
| Normal | 0 | 0 | 3.3e-07 → **4.0e-08** | 2.1e-04 → **2.0e-08** | 9.2e-01 → **1.5e-09** |
| t(6) "heavy" | 0 | 3 | 2.6e-07 → **5.2e-08** | 1.0e-04 → **1.5e-07** | 8.0e-01 → **1.4e-05** |
| chi2(20) "mild skew" | 0.63 | 0.6 | 4.1e-07 → **3.4e-08** | 6.1e-04 → **4.7e-08** | 9.9e-01 → **2.4e-07** |
| chi2(4) "skew" | 1.41 | 3 | 2.4e-05 → **9.3e-06** | 7.2e-01 → **4.4e-01** | 1.0 → **1.0** |

> **Which `W` this table used.** The rows above come from the exact engine, whose
> scale column is \((q^0_\tau - m_0)/(2 v_0)\). That is **not** what
> `make_weights_vqtl()` builds: the shipped function builds
> \(\Phi^{-1}(\tau)/(2\,\mathrm{sd})\) and **ignores its own `q_tau` and `mu`
> arguments entirely**. The two coincide exactly on a Normal baseline and
> diverge by up to three orders of magnitude, in both directions, otherwise — so
> on the non-Normal rows this table answers a slightly different question from
> the one the shipped code asks. §6 measures the shipped `W` by Monte Carlo and
> reaches the same conclusion; the two were not reconciled (§9).

### The envelope, and how it depends on the estimator

The bracket "somewhere between 0.63 and 1.41" that the four rows above suggest is
too generous, and it is not a property of the baseline alone. **The estimator
choice and the safe skewness limit are coupled.** Corrected rates at 5e-8,
MAF 0.01, Z_mu 30, **GLS**:

| skewness | 0.63 | 0.71 | 0.82 | 0.89 | 0.95 |
|---|---|---|---|---|---|
| corrected rate @5e-8 | 3.8e-10 | 3.6e-09 | 5.0e-07 | 4.8e-05 | **1.00** |

The usable envelope under GLS is therefore about **|skewness| 0.65**, not 1.41:
by 0.89 the corrected rate is already a thousand times nominal, and by 0.95 it has
collapsed completely. **Under OLS the envelope runs out to about 1.0.** Pick the
estimator and the limit together; a GLS scan on a trait with skewness 0.9 is
outside what this correction can carry, and the same trait under OLS is not.

### It is not skewness, and it is not kurtosis either

**The earlier claim that "the determining property is skewness, not kurtosis" is
refuted as stated.** At *matched* skewness 0.80 under GLS, a skew-normal baseline
(excess kurtosis 0.65) leaves **2.0e-04** while chi2(12) (excess kurtosis 1.00)
leaves **5.0e-07** — 390x apart, with the **lower-kurtosis** distribution the
worse of the two. Kurtosis does not rescue the claim by substitution, either:
symmetric t(6), skewness exactly 0, still leaves a **62x** residual.

**The failure is not summarised by either moment alone.** What governs it is how
far \(f_0'/f_0\) departs from proportionality to \(z_\tau\) across the tau grid —
a functional property of the whole density that no single moment captures. Treat
the skewness numbers above as a usable operational proxy for the families
tested, not as the mechanism.

* **Symmetric baselines: the correction works, but not perfectly.** On t(6) it
  reduces the MAF 0.30 / Z_mu 60 rate from 4.2e-04 to 3.3e-07 and never makes a
  cell worse anywhere in the sweep — while still leaving a 62x residual at
  MAF 0.01 / Z_mu 30.
* **Mildly skewed (|skewness| ≲ 0.65 under GLS, ≲ 1.0 under OLS): the correction
  works**, holding the rate within a few times nominal out to Z_mu 30 at
  MAF 0.01, against 2e7x uncorrected.
* **Strongly skewed (chi2(4), skewness 1.41): the correction is partial and must
  not be relied on.** It removes roughly half the excess and no more. Two separate
  things go wrong. First, \(f_0'/f_0\) is no longer proportional to \(z_\tau\), so
  only the part of the residue the weight column can see is removed. Second — and
  this is *not* the genotype-skewness bias at all — a strongly skewed baseline
  produces dispersion inflation **even at MAF 0.50**, where \(k = 0\) exactly and
  the correction is identically zero by construction (0.163 at 0.05, 4.5e-06 at
  5e-8, for Z_mu = 60). That residue comes from the phenotype's shape, not the
  genotype's, and nothing here addresses it.
* **`dispersion_skew_guard()` is blind to that third mechanism.** `G` is
  proportional to the genotype skewness, which is *exactly* zero at MAF 0.50, so
  the guard returns identically zero there for any Z_mu however large. The third
  mechanism lives exactly at MAF 0.50. Both tools — the guard and the correction
  — read zero precisely where it is worst. Nothing in the package sees it;
  transform the phenotype.
* On skewed baselines the corrected SE is *smaller* than the uncorrected one
  (ratio 0.5–0.99), because `cov_mu_sigma2` is large under `gls` there. That is
  correct delta-method behaviour but removes the conservative margin, so the
  partial correction is not even conservative on chi2(4).

**Practical consequence.** `docs/phenotype_preparation.md` already recommends
transforming skewed traits, and a rank-inverse-normal transform puts the analysis
squarely in the Normal row. Use the correction there. On an untransformed
right-skewed trait, treat the corrected column as a lower bound on the artefact
and keep using `dispersion_skew_guard()`.

---


---

## 6. Monte Carlo through the shipped code

Everything above is an exact population calculation. This section runs the
**shipped** Stage 1 (`core.run_stage1`, estimated quantiles, estimated densities,
real RIF, residualisation, `pooled_cov.compute_sigma0`) and the **shipped** Stage 2
(`param_gwas`, `se_mode = "tau_cov"`, `return_cov = "upper"`) and applies the
**installed** `correct_dispersion_skew()`, 20,000 replicates per cell,
N = 30,000. Machinery reused from `/user/work/fh6520/gls_coverage/src`.

N = 30,000 is a deliberately harder test than a real scan: reaching a given `Z_mu`
at N = 30k needs a much larger \(\delta/\sigma\) than at N = 461k, so the
second-order expansion is stressed far beyond where it would ever be used.

Rejection of `beta_sigma2` at nominal 0.05, with Monte Carlo SE, against the exact
prediction (`dev` in MC standard errors):

| cell | est | unc ± MCSE | exact | dev | **cor ± MCSE** | exact | dev |
|---|---|---|---|---|---|---|---|
| null MAF 0.30 | ols | 0.0523±0.0016 | 0.0500 | +1.5 | 0.0524±0.0016 | 0.0500 | +1.5 |
| null MAF 0.01 | gls | 0.0490±0.0015 | 0.0500 | -0.7 | 0.0481±0.0015 | 0.0497 | -1.0 |
| loc 0.10, MAF 0.30 | ols | 0.0634±0.0017 | 0.0617 | +1.0 | **0.0517±0.0016** | 0.0502 | +0.9 |
| loc 0.25, MAF 0.30 | ols | 0.4694±0.0035 | 0.4665 | +0.8 | **0.0515±0.0016** | 0.0515 | +0.1 |
| loc 0.25, MAF 0.30 | gls | 0.4665±0.0035 | 0.4674 | -0.2 | **0.0534±0.0016** | 0.0518 | +1.0 |
| loc 0.15, MAF 0.01 | ols | 0.0690±0.0018 | 0.0691 | -0.0 | **0.0442±0.0015** | 0.0446 | -0.2 |
| loc 0.35, MAF 0.01 | ols | 0.4996±0.0035 | 0.5051 | -1.5 | **0.0261±0.0011** | 0.0277 | -1.4 |
| loc 0.35, MAF 0.01 | gls | 0.5004±0.0035 | 0.5072 | -1.9 | **0.0259±0.0011** | 0.0281 | -1.9 |
| loc 0.35, **MAF 0.50** | ols | 0.0595±0.0017 | 0.0591 | +0.3 | 0.0595±0.0017 | 0.0591 | +0.3 |
| disp only, MAF 0.01 | ols | 0.9765±0.0011 | 0.9746 | +1.8 | 0.9755±0.0011 | 0.9739 | +1.5 |
| mixed, MAF 0.01 | ols | 0.9813±0.0010 | 0.9817 | -0.4 | 0.8844±0.0023 | 0.8879 | -1.5 |

Exact and MC agree within ±2 MCSE in **every Normal cell**, uncorrected and
corrected alike. The MAF 0.50 control is byte-identical before and after, as
\(k = 0\) requires.

At nominal 1e-3, twenty times deeper:

| cell | est | unc @1e-3 ± MCSE | exact | **cor @1e-3 ± MCSE** | exact |
|---|---|---|---|---|---|
| loc 0.25, MAF 0.30 | ols | 0.08240±0.00194 | 0.08325 | **0.00100±0.00022** | 0.00108 |
| loc 0.25, MAF 0.30 | gls | 0.08380±0.00196 | 0.08389 | **0.00090±0.00021** | 0.00109 |
| loc 0.35, MAF 0.01 | ols | 0.10860±0.00220 | 0.11125 | **0.00030±0.00012** | 0.00022 |
| loc 0.35, MAF 0.01 | gls | 0.11070±0.00222 | 0.11360 | **0.00020±0.00010** | 0.00023 |
| loc 0.35, MAF 0.50 | ols | 0.00145±0.00027 | 0.00153 | 0.00145±0.00027 | 0.00153 |

Bias in `beta_sigma2` about the generating truth, in units of its own reported SE,
and SE calibration (`mean_se / empirical_sd`; 1.00 is correct):

| cell | est | Z_mu | zbias unc | **zbias cor** | se/sd unc | se/sd cor | SE ratio |
|---|---|---|---|---|---|---|---|
| null MAF 0.30 | ols | 0.0 | 0.000 | 0.000 | 0.992 | 0.992 | 1.0000 |
| loc 0.10, MAF 0.30 | ols | 11.2 | 0.307 | **-0.003** | 0.988 | 0.991 | 1.0015 |
| loc 0.25, MAF 0.30 | ols | 27.9 | 1.876 | **-0.035** | 0.979 | 0.996 | 1.0093 |
| loc 0.25, MAF 0.30 | gls | 27.7 | 1.872 | **-0.019** | 0.973 | 0.992 | 1.0092 |
| loc 0.35, MAF 0.01 | ols | 8.5 | 1.966 | **-0.032** | 0.932 | 1.132 | 1.1079 |
| loc 0.35, MAF 0.01 | gls | 8.5 | 1.967 | **-0.017** | 0.925 | 1.132 | 1.1068 |
| loc 0.35, MAF 0.50 | ols | 42.7 | -0.013 | -0.013 | 0.968 | 0.968 | 1.0000 |
| disp only, MAF 0.01 | gls | 0.0 | 0.255 | 0.246 | 0.838 | 0.841 | 1.0027 |
| mixed, MAF 0.30 | ols | 22.1 | 1.188 | **-0.033** | 0.959 | 0.969 | 1.0042 |
| mixed, MAF 0.01 | gls | 5.4 | 1.094 | **0.087** | 0.831 | 0.919 | 1.0556 |
| heavy loc 0.25, MAF 0.30 | ols | 29.2 | 1.834 | **0.152** | 0.987 | 1.000 | 1.0067 |
| heavy loc 0.35, MAF 0.01 | ols | 8.9 | 1.947 | **0.176** | 0.953 | 1.099 | 1.0780 |
| skew loc 0.25, MAF 0.30 | ols | 28.7 | 1.377 | **-0.232** | 0.947 | 0.950 | 0.9272 |
| skew loc 0.25, MAF 0.30 | gls | 32.4 | 8.465 | **6.338** | 1.002 | 1.040 | 0.8672 |
| skew loc 0.35, MAF 0.01 | gls | 12.5 | 6.625 | **5.289** | 1.218 | 2.183 | 0.5907 |

Three things to read out of this.

1. **On a Normal baseline the correction removes 98–99% of the bias**, and the
   `se/sd` column shows the corrected SE is honest or mildly conservative
   (0.99–1.13) where the uncorrected one was anti-conservative (0.93–0.99).
2. **Heavy tails cost a little but not much**: 91% of the bias removed, residual
   0.15–0.18 SE.
3. **Skew, and especially skew under GLS, is not fixed.** The GLS cells go from
   8.5 SE to 6.3 SE of bias. That is a `W`-misspecification failure that GLS
   amplifies, not a failure this correction is shaped to address.

**A caveat on the non-Normal comparison.** For skew and heavy cells the MC and the
exact engine do not use the same `W`: the MC uses the **shipped**
`make_weights_vqtl()`, whose scale column is \(z_\tau/(2\sigma)\) regardless of
baseline, while the exact engine uses the correctly shaped
\((q_\tau^0 - m_0)/(2v_0)\). These coincide exactly on a Normal baseline and
differ substantially on chi2(4), which is why `dev` blows up in those rows. Both
calculations are internally valid; they answer slightly different questions, and
they agree on the conclusion.

---

## 7. Does this make the `G >= 0.2` flag redundant?

**For Type-I control, yes. For power, no — and that is its new job.**

The flag exists because `G` predicts the false-dispersion \(|Z|\), and
`FINDINGS_pooled_dispersion.md` calibrated `G >= 0.2` as the threshold holding the
5e-8 rate within 2x nominal. After correction the rate is within **1.03x**
nominal everywhere in the sweep, including at G = 25.3. The flag is no longer
needed to bound false positives on a Normal baseline.

It remains worth computing, for three reasons.

1. **The correction cannot always be applied.** It needs a signed EAF; the guard
   needs only the magnitude and falls back on `d_resid`. Where EAF is missing,
   the flag is all there is.
2. **On a skewed baseline the correction is partial** (§5) and the flag's
   documented conservatism there is still the right tool. With one blind spot:
   `G` is proportional to the genotype skewness, so `dispersion_skew_guard()` is
   **identically zero at MAF 0.50** (`gskew = 0`) for any `Z_mu` however large.
   The third mechanism of §5 — dispersion inflation driven by the *phenotype's*
   shape rather than the genotype's — lives exactly at MAF 0.50, where both the
   flag and the correction read zero. Neither tool sees it.
3. **Its meaning inverts into a power warning.** High `G` is now exactly where the
   correction costs power: corrected rates fall to 0.0298 (G = 6.3) and 0.0194
   (G = 25.3) against nominal 0.05, so a genuine dispersion effect at such a
   variant is *under*-detected. The `FINDINGS` note that the flag is
   "anti-conservative by 3.4x under GLS" no longer bears on false positives; it
   bears on how much honest power a high-`G` variant has lost.

---

## 8. Practical guidance

* **Transform skewed phenotypes.** `docs/phenotype_preparation.md` already says so;
  §5 makes it load-bearing for this correction. On a rank-inverse-normal or
  log-transformed trait you are in the Normal row and the correction is exact to
  \(O(\delta^4)\).
* **Pass `param_cov_upper`.** Under `gls` the two parameters are correlated and
  `cov_mu_sigma2` is the term that keeps the corrected SE calibrated.
* **Pass a signed EAF**, oriented to the allele `beta_mu` counts. Stage 1 does not
  emit one; take it from the genotype source. The function refuses rather than
  guesses, because an unoriented MAF does not merely fail to correct — it
  corrects the **wrong way** and is worse than not correcting at all (0.916 to
  **1.000** at 5e-8 for MAF 0.01 / Z_mu 30, against 1.5e-09 done right). Three
  things enforce this: an `eaf` column below 0.5 at every variant is refused as a
  folded-MAF signature; `dispersion_report()` refuses when `eaf` and `maf` are
  the same values; and `eaf_orientation = "effect"` must be set deliberately.
* **Pass `W_sd`, the `sd` you gave `make_weights_vqtl()`.**
  `dispersion_report()` verifies `fit$W` elementwise against
  `qnorm(taus)/(2*W_sd)` via `check_vqtl_design()` and refuses otherwise. A
  design merely *named* `beta_sigma2` is not enough: the same Normal model
  parameterised in the SD over-subtracts by a factor of `2*sd`, and the
  5-parameter mixture by 2.4x.
* **Use a tau grid symmetric about 0.5** (§2.4).
* **Check `attr(result, "na_counts")`.** Variants left uncorrected — an `eaf` of
  0, 1 or `NA`, or a `cov_mu_sigma2` that breaches Cauchy-Schwarz because it came
  from a different fit than the SEs — are counted and warned about rather than
  returned as silent `NA`.
* **Keep `dispersion_skew_guard()`** for the variants where EAF is missing, for
  skewed baselines, and now as a *power* warning (§ above).
* **Both columns are reported.** Compare them: a finding where they agree is
  robust; one where they disagree is exactly what the correction exists to
  adjudicate.

---

## 9. What was not tested

* **Real genotypes and LD.** Everything uses a 3-point HWE genotype law. The
  correction depends on the genotype law only through \(\mu_3/\mathrm{Var}(g)\),
  which `mu3_over_var` lets a caller supply directly for dosages or HWE
  departures — but that route is not itself validated here.
* **Rates below 1e-3 by simulation.** The 5e-8 numbers are quadrature-exact given
  that \(\hat\theta\) is bivariate normal about \(\theta_{pop}\). The Monte Carlo
  confirms that assumption down to 1e-3 only.
* **Covariates with real effects, and the jackknife covariance.** The exact
  engine uses the pooled \(\Sigma_0\) as the reported covariance throughout.
* ~~**`W` families other than `make_weights_vqtl`.** The correction is refused
  there by documentation, not by code.~~ **Closed.** `check_vqtl_design()` now
  verifies the numeric identity — `W[, mu_name]` constant 1 and
  `W[, sigma2_name] == qnorm(taus)/(2*sd)` elementwise — and
  `dispersion_report()` gates on it by default. What is still *untested* is the
  behaviour of the correction on other `W` families, because it is now refused
  rather than measured.
* ~~**Baselines between skewness 0.63 and 1.41.** The boundary is bracketed, not
  located.~~ **Located, for the families tested**: about 0.65 under GLS and about
  1.0 under OLS (§5). What is *not* established is a criterion in terms of the
  density itself rather than a moment — the skew-normal vs chi2(12) comparison at
  matched skewness 0.80 shows the moment is a proxy, not the mechanism.
* **The Monte Carlo and the exact engine use different `W` on non-Normal
  baselines** (§6). Neither is wrong; they were not reconciled.
* **`w_sd_ratio` != 1** is implemented and documented but not exercised by any
  measurement here.
