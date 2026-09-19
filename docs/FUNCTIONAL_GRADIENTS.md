# Functionals of the fitted distribution: why not to average the tau grid

Status: package capability, 2026-09-19. Code in `stage2/R/functional-grad.R`;
unit checks in `tests/tier0/t0_functional_grad.R`.

## The gap

Stage 2 supplies a weight matrix `W` whose columns are the quantile
derivatives of a parametric family,

```
W_j(tau) = d q_tau / d theta_j,
```

and fits `theta` with a covariance. It reports per-parameter effects. It has
no route to a *functional* of the distribution -- a mean, a variance, a tail
probability `P(Y > c)`. The common workaround is to average the tau grid, for
example `colMeans(W) %*% theta`, or an unweighted projection with `W = 1`.

That workaround is only valid when the quantile map is smooth and well
conditioned. When it is not, it is materially wrong, and it is wrong in a way
that no per-tau diagnostic shows.

## The worked example: a near-step mixture

On real testosterone data a two-component skew-normal mixture separates into
male and female components that barely overlap. The fitted `W` then has a
`gamma` column that spikes to `-16.57` at `tau = 0.5`, against about `-1.2`
and `-1.8` at the neighbouring taus: the transition between the components is
almost a step, and the quantile map is near-singular in that direction.

For this fixture (components below):

| route | `d E[Y] / d gamma` |
|---|---|
| tau-grid average of the `gamma` column | **-0.7546** |
| exact analytic value | **-0.6104** |

The grid route inflates the marginal mean effect by about 24 per cent. The
`beta_left` and `beta_right` columns are unaffected: their tau averages
(`0.4938`, `0.5013`) agree with the analytic values (the component weights,
`0.5001`, `0.4999`). Only the near-singular direction is wrong.

The fitted fixture, used verbatim in the tier-0 test:

| component | weight | xi | omega | alpha |
|---|---|---|---|---|
| lower | 0.500095747692 | 0.296217247249 | 0.551173832619 | -0.927107313477 |
| upper | 0.499904252308 | 2.653440657413 | 0.374512611297 | -1.032212367482 |

A milder version of the same problem affects the shifted-lognormal BMI
family: the grid projection gave mean-effect slopes of 0.98 against a REGENIE
reference where the exact route gave 0.85-0.96.

## The correct recipe

Do not average the grid. Differentiate the functional on the parametric
family. The fitted baseline parameters are already stored in the weight
object, so the exact answer is available:

```r
# The functional: the distribution mean, written in the baseline parameters.
f_mean <- function(p) {
  d <- p$alpha / sqrt(1 + p$alpha^2)
  E <- p$xi + p$omega * d * sqrt(2 / pi)          # skew-normal component mean
  p$weight[1] * E[1] + p$weight[2] * E[2]
}

# df/dtheta, in the same units and order as colnames(object$W):
g <- functional_grad(object, f_mean)
g[["gamma"]]        # -0.6104, not the grid value -0.7546

# A Stage 2 theta and its covariance give an estimate and a delta-method SE:
functional_estimate(object, f_mean, theta = theta_hat, cov_theta = cov_theta)
```

`f` is any smooth function of the baseline parameters. The package supplies
the parameterisation (how each Stage 2 column moves the distribution) and
takes a central finite difference; the user supplies only the functional.
The `W` values are not used by the gradient -- only the column names and
their order -- because the gradient is a property of the family, not of the
quantile map that happened to be fitted.

## Parameter semantics by weight family

The Stage 2 parameters are the columns of `W`. The table gives what each one
does to the distribution, its sign convention, and the baseline-parameter
list that `f` receives.

### Skew-normal mixture (`matched_skew_normal_mixture`; stores `components`)

Baseline: `Y ~ w_lower SN(xi_lower, omega_lower, alpha_lower) +
w_upper SN(xi_upper, omega_upper, alpha_upper)`, components ordered by
location. `f` receives `list(weight, xi, omega, alpha)`, each a length-two
vector in `(lower, upper)` order.

| column | distribution parameter moved | sign convention |
|---|---|---|
| `gamma` | log-odds of membership in the **lower** component: `w_lower = plogis(gamma)` | **Increasing `gamma` moves mass into the lower component.** Its `W` column mean is negative whenever the upper component has the larger mean, and its mean-functional gradient is `-w_upper (1 - w_upper) (E_upper - E_lower)`. This sign is counter-intuitive and has been mis-guessed in review. |
| `beta_left` | location `xi_lower` of the lower component | increasing shifts the lower component right; `d E[Y] / d beta_left = w_lower` |
| `beta_right` | location `xi_upper` of the upper component | increasing shifts the upper component right; `d E[Y] / d beta_right = w_upper` |

### Shifted lognormal (`shifted_lognormal`, `matched_shifted_lognormal`)

Baseline: `Y = shift + L`, `L ~ LogNormal(meanlog, sdlog)`, so
`log(Y - shift) ~ N(meanlog, sdlog^2)` for `Y > shift`. `f` receives
`list(shift, meanlog, sdlog)`.

| column | distribution parameter moved | sign convention |
|---|---|---|
| `additive` | the shift | increasing shifts the whole distribution right; `W` column is exactly 1 at every tau, and `d E[Y] / d additive = 1` |
| `multiplicative` | the log-scale location `meanlog` | increasing multiplies the positive part by `exp(.)`; `d E[Y] / d multiplicative = exp(meanlog + sdlog^2 / 2)` |

### Log-scale vQTL (`bmi_lognormal_meanlog_sigma2`, `matched_logscale_vqtl`, `vqtl_lognormal_mu_sig2`)

Baseline: the same shifted lognormal. `f` receives
`list(shift, meanlog, sdlog)`. `sigma2` is the variance of the log,
`sdlog^2`.

| column | distribution parameter moved | sign convention |
|---|---|---|
| `meanlog` / `b1_mu` | the log-scale mean `meanlog` | as `multiplicative` above |
| `sigma2` / `g1_sigma2` | the log-scale variance `sdlog^2` | increasing widens the distribution; `sdlog = sqrt(sdlog0^2 + theta)` |

### Normal mean/variance vQTL (`make_weights_vqtl`; caller attaches `params0`)

Baseline: `Y ~ N(mu, sigma2)`. `make_weights_vqtl()` returns a bare matrix,
so the caller attaches the baseline parameters as
`params0 = list(mu = , sigma2 = )` (or `sd`). `f` receives
`list(mu, sd, sigma2)`.

| column | distribution parameter moved | sign convention |
|---|---|---|
| `beta_mu` | the mean `mu` | increasing shifts the distribution right; `d E[Y] / d beta_mu = 1` |
| `beta_sigma2` | the variance `sigma2` | increasing widens the distribution; `sd = sqrt(sigma2)` |

An object may carry `param_map`, a function `function(theta) -> params`, to
support any other family. This overrides the built-in dispatch.
