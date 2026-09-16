# What G-MoQ's standard errors cost, and why

**Status: derived and verified numerically 2026-09-16, then tested against four
real phenotypes. Accuracy floor of the exercise is +/-3%; see caveats.**

## Question

G-MoQ recovers component-specific effects from a mixture phenotype without
knowing membership. A stratified GWAS is told membership. How much must G-MoQ
lose, and does it lose more than that?

For testosterone the comparison is exact: REGENIE female N=199,572 plus male
N=200,138 equals 399,710, which is G-MoQ's N. The same people.

> **Comparator caveat (added 2026-09-16).** Every REGENIE ratio in this document
> is measured against a step 2 run `--ignore-pred`, i.e. with **no LOCO of its
> own** — the sample sizes match but the polygenic adjustment does not. The
> `gamma` comparator below is additionally *not* N-matched (462,869 vs 399,710)
> and is a different sample, not a superset restriction. A fair comparator with
> a proper within-stratum LOCO, and the restated ratios against it, are in
> [COMPARATOR_LEDGER.md](COMPARATOR_LEDGER.md). The bound and the decomposition
> are unaffected; the *observed* ratios they are compared with move.

## The bound

With `Y | g, Z=k ~ N(alpha_k + beta_k g, sigma_k^2)` and
`P(Z=1|g) = expit(eta + gamma g)`, the scores are

```
s_{alpha1} = w1*eps1/sigma1^2      s_eta = w1 - pi
s_{beta1}  = g * s_{alpha1}        s_gamma = g * s_eta
```

Two structural facts follow:

1. **The genotype factors out exactly.** Every g-multiplied score is `g x` a
   non-g score and the residual law does not depend on `g`. Requires `g`
   centred, `g` independent of `Z` (true for autosomes), and
   `|beta1-beta2|*g << |alpha1-alpha2|` -- about 0.01 SD against 2 SD for a
   GWAS variant, so exact to four decimals.

2. **Non-g nuisances are free; the g-block is what costs.** Because `E[g]=0`,
   every cross-term between the g-block and `{alpha1, alpha2, eta, sigma,
   covariates}` vanishes. So estimating `pi`, `sigma` and the intercepts jointly
   costs nothing. What costs is that `beta2` (and `gamma`) must be co-estimated:

```
SE ratio^2 = pi * (J^-1)_11 / sigma1^2        Var(g) cancels identically
```

**The single-parameter version is the wrong quantity.** `1/J_11` -- the
information with everything else known -- saturates at `1/sqrt(pi) = 1.414` as
separation goes to zero. The correct `(J^-1)_11` **diverges**, because at `S=0`
the two component scores become collinear and `J` is singular. The observed
shrink-0.7 ratio of 1.84 exceeds 1.414, so the joint version is the one reality
is testing.

`J` has no closed form -- the integrand is logistic-squared times Gaussian -- but
the limits do, verified to four decimals:

```
S -> 0 :  SE ratio -> 1 / (S * sqrt(2*pi_k));  with gamma free and pi=0.5, 2*sqrt(3)/S^2
S -> inf: SE ratio - 1 = O(S * exp(-S^2/8))
```

### Predicted SE ratio (pi = 0.5, equal variances)

| S | 0.5 | 1.0 | 1.5 | 2.0 | 3.0 | 4.0 | 5.0 | 6.0 |
|---|---|---|---|---|---|---|---|---|
| beta1, beta2 free | 2.239 | 1.443 | 1.285 | 1.261 | 1.202 | 1.096 | 1.034 | 1.010 |
| + gamma free | 15.20 | 4.646 | 2.600 | 1.847 | 1.290 | 1.104 | 1.035 | 1.010 |

At the separation testosterone actually has, the mixture penalty is about **1-2%**.

## The test

Estimated separations, from true sex labels on the covariate-residualised
phenotype (`S = delta_mu / sqrt((sigma_f^2+sigma_m^2)/2)`); the shrink
construction scales the male mean, so `S(x) = (1-x)*S(0)` exactly:

| phenotype | S | (1-x)*S(0) |
|---|---|---|
| main | 5.89 | -- |
| shrink 0.2 | 4.75 | 4.71 |
| shrink 0.5 | 2.99 | 2.94 |
| shrink 0.7 | 1.81 | 1.77 |

Predicted versus observed, using G-MoQ's own fitted skew-normal mixture and its
actual 45x3 weight matrix rather than the idealised Gaussian curve:

| parameter | obs/predicted across the four phenotypes |
|---|---|
| `beta_left` | 1.011, 1.018, 1.021, 1.016 |
| `beta_right` | 0.984, 0.988, 0.981, 0.972 |
| `gamma` | 0.998, 0.998, 0.990, 0.983 |

Twelve cells: mean 0.995, range 0.972-1.021, RMS deviation 1.5%.
**G-MoQ sits on the information bound.**

### Decomposing the observed 12% (female, main)

| source | SE factor |
|---|---|
| mixture information loss at S=5.89 | x1.019 |
| 45 quantiles on (0.02, 0.98) instead of the full likelihood | x1.009 |
| OLS instead of GLS projection | x1.002 |
| **LOCO offset** | **x1.075** |
| unexplained | x1.011 |
| observed | 1.121 |

The information cost is 2%. Most of the rest was the broken LOCO adjustment --
see [LOCO_MIXTURE_HAZARD.md](LOCO_MIXTURE_HAZARD.md), where removing it moved the
observed ratio to 1.0236, confirming this decomposition experimentally.

### gamma

Its comparator is a REGENIE GWAS of sex itself on N=462,869 against G-MoQ's
399,710, carrying a fixed `sqrt(462869/399710) = 1.0761` factor. Divide it out
and the observed main-phenotype ratio is **0.9994**: at full separation G-MoQ
infers membership from the phenotype alone as precisely as a GWAS that is told
it.

## Why the density kernel could not have helped

With `W_{t,j} = -(dF/dtheta_j)(q_t)/f(q_t)` and
`Sigma_{st} = (min(tau_s,tau_t) - tau_s*tau_t)/(f(q_s) f(q_t))`, the density
cancels in the quadratic form:

```
W' Sigma^-1 W  =  D' V^-1 D  =  Fisher information of the multinomial
                                obtained by binning Y at the cut points q_t
```

Verified to 1.6e-14 at T = 45, 199, 999. **The optimally weighted RIF projection
is exactly a grouped-data (minimum-chi-square) estimator, and density estimation
drops out of the efficiency calculus.** This independently explains the
fourth-order kernel result (GitHub issue #1): an 11% change in `fhat` moved the
SE by 0.06% because `fhat` was never in the efficiency term.

Convergence to the MLE, `(W' Sigma^-1 W)^-1` for `beta_left`:
0.4555 (T=45) -> 0.44999 (T=199) -> 0.44762 (T=999) against exact `J^-1` =
0.44742. Full asymptotic efficiency in the limit; 45 quantiles cost 0.9%.

Citable: Zhao & Xiao (2014, *Econometric Theory* 30(6), Thms 8-9) -- the
optimally weighted k-quantile estimator attains the Cramer-Rao bound as
k -> infinity, with an explicit finite-k bound. General-theta version:
Carrasco & Florens (2000, 2014). Finite-T ARE: Ogawa (1951), Mosteller (1946),
Sarhan & Greenberg (1962), Eubank (1981). **Not** Chernoff-Gastwirth-Johns
(1967), which is the T=infinity L-estimator result.

Independent check on the whole W/Sigma machinery: predicted OLS/GLS SE ratio for
`beta_left` main = 1.0081, against the observed arm difference
1.12088/1.11204 = 1.0080.

**Caveat on the OLS arm.** Zhao & Xiao Thm 10 shows equally-weighted *averaging*
over quantiles is asymptotically equivalent to least squares and gains nothing
from more quantiles. G-MoQ's OLS is a projection through a non-constant `W`, not
a flat average, so the theorem does not bite -- and the measured OLS penalty is
0.2% at main separation, rising to 4% at shrink 0.7. But only the GLS arm has an
efficiency guarantee.

**No published efficiency theory exists for RIF/unconditional quantile
regression specifically.** Firpo-Fortin-Lemieux (2009) prove consistency and
asymptotic normality only; no semiparametric bound, no MLE comparison. Firpo
(2007) and Froelich & Melly (2013) cover quantile *treatment* effects, not this.
That gap is why the answer here was computed rather than cited.

## Caveats

- The `beta_right` and `gamma` cells land 1-3% **below** a quantity that is a
  lower bound, which is impossible. That is misspecification in the quadrature
  (the fitted skew-normal mixture is not the true density, and the
  Brownian-bridge Sigma is slightly larger than the covariate-residualised
  pooled Sigma actually used). Treat +/-3% as the accuracy floor.
- Separation was estimated using true sex labels. Unlabelled EM recovers
  S = 6.09, 4.94, 3.01, 1.56 -- good through shrink 0.5, degrading at shrink 0.7
  where the fit destabilises (pi drifts to 0.63).
- `Var(g)` was assumed common across the combined and stratified samples. It
  cancels algebraically; residual mis-cancellation is second order.

## Simulation check against formal mixture models (2026-09-16)

`/user/work/fh6520/gmoq_mixture_sim` -- two-component Gaussian mixture, pi=0.5,
MAF 0.3, N=100,000, 1,000 paired replicates per cell, S in {1,2,3,4} plus
equal-effect and unequal-variance cells. The simulation's Stage 1/Stage 2
reimplementation reproduces the installed package bit-for-bit (max |diff| = 0 on
q_tau, fhat and b(tau); 3.5e-17 on SE), so it exercises the real estimator.

Efficiency against an oracle that knows membership:

| S | G-MoQ GLS | EM MLE | flexmix (informed init) | CR bound |
|---|---|---|---|---|
| 1 | 1.43 / 1.52 | 1.44 / 1.43 | 1.43 / 1.47 | 1.95 / 2.05 |
| 2 | 1.26 / 1.29 | 1.28 / 1.29 | 1.33 / 1.37 | 1.26 / 1.24 |
| 3 | 1.22 / 1.23 | 1.22 / 1.21 | 1.21 / 1.14 | 1.25 / 1.20 |
| 4 | 1.10 / 1.13 | 1.09 / 1.11 | 1.06 / 1.09 | 1.09 / 1.12 |

**A formal mixture model recovers nothing G-MoQ leaves behind.** The curves
coincide at every separation. G-MoQ's disadvantage against the oracle grows as
components merge, but that is the information content of the data -- the MLE
loses exactly as much.

At S=1 nothing attains the asymptotic bound of 1.95-2.05: the information matrix
has condition number 4.2e4 and N=1e5 is not in the asymptotic regime, so every
mixture-based estimator is effectively regularised and lands near 1.43-1.53.
The bound is correct; the regime is not asymptotic.

Confirmed on 100 real chr22 SNPs at N=399,710: a hand-rolled EM is
indistinguishable from G-MoQ on every axis (r 0.977-0.978 vs 0.979-0.981; median
SE ratio 1.127/1.074 vs 1.121/1.068), with 0/100 convergence failures.

### Practical cost of the alternatives

| | per SNP | chr22 (17,550 variants) | cores |
|---|---|---|---|
| G-MoQ Stage 1 | ~43 SNPs/s | ~7 min | **8 threads** |
| G-MoQ Stage 1, single-core equivalent | ~5.4 SNPs/s | ~54 min | 1 |
| EM mixture regression | 1.08 s | 5.3 h | 1 |
| `flexmix` (defaults) | 6.9 s | 34 h | 1 |
| `mixtools::regmixEM` | 112-2,429 s | prohibitive | 1 |

**Compare like with like.** The measured G-MoQ run used `--threads 8`; every
mixture fit was single-core. Against the EM the honest figure is therefore about
**6x**, not the ~45x a naive reading of the wall-clock gives. The durable
advantage is not raw speed: Stage 1 emits all 45 quantiles once, and that single
output serves any distributional model, where each mixture fit answers one
question for one SNP.

The G-MoQ timing also *includes* the per-SNP block jackknife (`--blocks 1000`)
and writing a 66 MB covariance payload, neither of which the production
`pooled` + `ols` path reads. An unprofiled but plausible share of that 7 minutes
is therefore work the default analysis discards. See BACKLOG H8.

`mixtools::regmixEM` timed out on 100% of replicates at S<=2 under a 1,200 s cap.
The hand EM fails to converge on 44.2% of replicates at S=1, with SEs 37-46% too
large there. **`flexmix` at default initialisation converged to a degenerate
one-component solution in 200/200 replicates in all four cells**, with biases of
1-2.5 empirical SDs; informed initialisation fixes it, but that requires the same
marginal mixture fit G-MoQ already performs. On real data at N=400k the same
failure appears as stopping at iteration 2 with `converged = TRUE`, both
component SDs equal to the whole-phenotype SD, and NaN standard errors on about a
third of variants -- flexmix's default `tolerance = 1e-6` is *relative* to a
log-likelihood of ~-5e5, making the effective threshold ~0.5.

### Resolved: the OLS projection penalty is the membership column

The simulation put the cost of the unweighted OLS projection at 4-5% (S=2),
24% (S=3) and **57% (S=4)**, growing with separation, while real testosterone at
S=5.89 showed **0.8%** -- a factor of ~70.

**The cause is neither the distribution family nor the covariance path. It is
that the simulation's W had two columns and testosterone's has three.** The
simulation used `make_weights_normal_mixture(include_membership = FALSE)`;
production's `build_skew_mixture()` always emits `gamma`. At matched separation
and matched K, the skew and normal families agree:

| W | K | S=2 | S=3 | S=4 | S=6 |
|---|---|---|---|---|---|
| normal mixture, means only | 2 | 1.031 | 1.232 | 1.593 | 3.295 |
| normal mixture, `+gamma` | 3 | 1.006 | 1.010 | 1.008 | 1.000 |
| skew mixture, means only | 2 | 1.019/1.086 | 1.213/1.328 | 1.568/1.732 | 3.091/4.404 |
| skew mixture, `+gamma` | 3 | 1.035/1.023 | 1.033/1.022 | 1.022/1.014 | 1.003/1.005 |

Deleting `gamma` from the *actual* shipped testosterone W moves it from 1.013 to
**2.132/2.509**.

The mechanism follows from the same `Sigma` algebra as the density-cancellation
result above. With `Sigma = diag(1/f) B diag(1/f)`, a bimodal `f` has a deep
trough, so `Sigma` collapses onto a single eigendirection -- 98.5% of the
spectrum at S=6. The membership weight
`W_gamma = -(p1(1-p1)/f)(F1-F2)` equals `-p1(1-p1)/f` at the trough, so **the
gamma column is that eigendirection**, sitting 0.011 deg from the leading
eigenvector at S=6 and 0.124 deg for real testosterone. Including it makes
`col(W)` a `Sigma`-invariant subspace, which is exactly the Gauss-Markov
condition, so OLS is already efficient. Omitting it leaves OLS nowhere to put
the trough noise, and it leaks into the component means.

The general diagnostic is `rho_max`, the largest canonical correlation between
the in-subspace and complement noise; `V_gls` is `V_ols` less a Schur
complement, and `1/sqrt(1-rho_max^2)` bounds the penalty. Across 38
configurations Spearman is 0.993 and the bound is never violated. Every shipped
W sits at `rho_max` 0.03-0.29 (penalty < 4%); the means-only W reaches 0.67-0.99.

**Default unchanged: OLS.** Full evidence, the mixture-family prediction and
the `Sigma`-misspecification fragility result are in
`sim/projeff/src`, with aggregate tables in `sim/projeff/tables`, register
`projection-efficiency-20260916-01`.
