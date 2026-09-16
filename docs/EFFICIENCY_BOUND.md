# What G-MoQ's standard errors cost, and why

**Status: derived and verified numerically 2026-09-16, then tested against four
real phenotypes. Accuracy floor of the exercise is +/-3%; see caveats.**

## Question

G-MoQ recovers component-specific effects from a mixture phenotype without
knowing membership. A stratified GWAS is told membership. How much must G-MoQ
lose, and does it lose more than that?

For testosterone the comparison is exact: REGENIE female N=199,572 plus male
N=200,138 equals 399,710, which is G-MoQ's N. The same people.

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
