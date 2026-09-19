#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 0 -- functional gradients on the parametric family.
#
# A weight matrix gives d q_tau / d theta.  It does NOT give a functional of
# the distribution, and averaging the tau grid is invalid when the quantile
# map is near-singular.  These checks pin the two routes apart on a fitted
# two-component skew-normal mixture whose gamma column spikes to -16.57 at
# tau = 0.5, and confirm the helper reproduces the exact answer and the
# analytic component derivatives.
#
# The fixture is the fitted clinical-testosterone mixture.  It is a handful of
# aggregate distribution parameters, hard-coded here so the test needs no
# private file and no optional package: W is used only for its column names,
# never for its values.
#
# Run:  Rscript --vanilla tests/tier0/t0_functional_grad.R
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier0/t0_functional_grad.R")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
gmoq_begin("tier0", "functional_grad")

f_mean_mix <- function(p) {
  delta <- p$alpha / sqrt(1 + p$alpha^2)
  E <- p$xi + p$omega * delta * sqrt(2 / pi)
  p$weight[1] * E[1] + p$weight[2] * E[2]
}

## =====================================================================
## 1. The testosterone-like skew-normal mixture: the observed bug
## =====================================================================
comp <- data.frame(
  component = c("lower", "upper"),
  weight = c(0.500095747692, 0.499904252308),
  xi     = c(0.296217247249, 2.653440657413),
  omega  = c(0.551173832619, 0.374512611297),
  alpha  = c(-0.927107313477, -1.032212367482)
)
obj_mix <- list(
  W = matrix(0, 3, 3,
             dimnames = list(NULL, c("gamma", "beta_left", "beta_right"))),
  components = comp,
  type = "matched_skew_normal_mixture"
)

g_mix <- functional_grad(obj_mix, f_mean_mix)

delta <- comp$alpha / sqrt(1 + comp$alpha^2)
E <- comp$xi + comp$omega * delta * sqrt(2 / pi)
w_lower <- comp$weight[1]; w_upper <- comp$weight[2]
gamma_analytic <- -w_upper * (1 - w_upper) * (E[2] - E[1])

gmoq_check("T0.FGRAD.mixture_gamma_is_analytic",
           paste("The mixture gamma gradient is the exact analytic",
                 "-w_upper (1 - w_upper) (E_upper - E_lower)."),
           isTRUE(abs(g_mix[["gamma"]] - gamma_analytic) < 1e-6),
           stat = abs(g_mix[["gamma"]] - gamma_analytic),
           criterion = "|grad_gamma - analytic| < 1e-6")

gmoq_close("T0.FGRAD.mixture_gamma_value",
           "The gamma gradient is -0.6104 for the testosterone-like fixture, not the grid value -0.7546.",
           g_mix[["gamma"]], -0.610398843791, 1e-4)

# The measured tau-average of the gamma column for this fixture.  This is the
# value the invalid grid route produces; the exact gradient must not equal it.
grid_gamma_measured <- -0.754629667951
gmoq_check("T0.FGRAD.mixture_gamma_is_not_grid_average",
           "The exact gamma gradient differs from the measured tau-grid average by more than 0.1.",
           isTRUE(abs(g_mix[["gamma"]] - grid_gamma_measured) > 0.1),
           stat = g_mix[["gamma"]] - grid_gamma_measured,
           criterion = "|grad_gamma - mean(W[, gamma])| > 0.1")

gmoq_close("T0.FGRAD.mixture_beta_left_is_weight",
           "The lower-component location gradient is the lower component weight.",
           g_mix[["beta_left"]], w_lower, 1e-8)
gmoq_close("T0.FGRAD.mixture_beta_right_is_weight",
           "The upper-component location gradient is the upper component weight.",
           g_mix[["beta_right"]], w_upper, 1e-8)
gmoq_check("T0.FGRAD.gradient_order_matches_W",
           "The gradient is named and ordered exactly as the W columns.",
           identical(names(g_mix), colnames(obj_mix$W)),
           criterion = "names(grad) == colnames(W)")

## =====================================================================
## 2. Smooth, low-separation case: the grid route is adequate there
## =====================================================================
# Two overlapping Normal components: the quantile map is well conditioned and
# the tau-grid average of the gamma column agrees with the exact gradient.
taus <- seq(0.02, 0.98, length.out = 45)
p1 <- 0.5; mu1 <- 0; sd1 <- 1; mu2 <- 1; sd2 <- 1
mix_cdf <- function(x) p1 * pnorm(x, mu1, sd1) + (1 - p1) * pnorm(x, mu2, sd2)
q_mix <- vapply(taus, function(tt)
  stats::uniroot(function(x) mix_cdf(x) - tt, c(-40, 60),
                 tol = .Machine$double.eps^0.75)$root, numeric(1))
W_smooth <- make_weights_normal_mixture_vqtl(taus, q_mix, p1, mu1, sd1, mu2, sd2)
obj_smooth <- list(
  W = W_smooth,
  components = data.frame(component = c("lower", "upper"),
                          weight = c(p1, 1 - p1),
                          xi = c(mu1, mu2), omega = c(sd1, sd2),
                          alpha = c(0, 0)),
  type = "matched_skew_normal_mixture"
)
g_smooth <- functional_grad(obj_smooth, f_mean_mix)
grid_smooth <- mean(W_smooth[, "gamma"])

gmoq_close("T0.FGRAD.smooth_matches_grid",
           "For a low-separation mixture the exact gradient agrees with the tau-grid average to 2 per cent.",
           g_smooth[["gamma"]], grid_smooth, 0.005,
           criterion = "|grad_gamma - mean(W[, gamma])| <= 0.005")

## =====================================================================
## 3. Shifted lognormal: additive is 1, multiplicative is analytic
## =====================================================================
p0 <- list(shift = 1.0, meanlog = 0.3, sdlog = 0.8)
ln_cdf <- function(y, p) ifelse(y > p$shift,
                                plnorm(y - p$shift, p$meanlog, p$sdlog), 0)
ln_pdf <- function(y, p) ifelse(y > p$shift,
                                dlnorm(y - p$shift, p$meanlog, p$sdlog), 0)
ln_grad_shift <- function(y, p) ifelse(y > p$shift,
                                       -dlnorm(y - p$shift, p$meanlog, p$sdlog), 0)
ln_grad_meanlog <- function(y, p) {
  out <- numeric(length(y)); ok <- y > p$shift
  z <- (log(y[ok] - p$shift) - p$meanlog) / p$sdlog
  out[ok] <- -dnorm(z) / p$sdlog
  out
}
W_ln <- make_weights_generic(
  taus, q_tau = qlnorm(taus, p0$meanlog, p0$sdlog) + p0$shift,
  dist_cdf = ln_cdf, dist_pdf = ln_pdf, params = p0,
  grad_funcs = list(additive = ln_grad_shift, multiplicative = ln_grad_meanlog)
)
obj_ln <- list(W = W_ln, params0 = p0, type = "shifted_lognormal")
f_mean_ln <- function(p) p$shift + exp(p$meanlog + p$sdlog^2 / 2)
g_ln <- functional_grad(obj_ln, f_mean_ln)

gmoq_close("T0.FGRAD.lognormal_additive_is_one",
           "The additive (shift) marginal-mean derivative is exactly 1.",
           g_ln[["additive"]], 1, 1e-8)
gmoq_close("T0.FGRAD.lognormal_multiplicative_analytic",
           "The multiplicative (meanlog) marginal-mean derivative is exp(meanlog + sdlog^2 / 2).",
           g_ln[["multiplicative"]], exp(p0$meanlog + p0$sdlog^2 / 2), 1e-6)
gmoq_close("T0.FGRAD.lognormal_additive_W_is_one",
           "The additive W column is exactly 1 at every tau, so its grid mean is also 1.",
           mean(W_ln[, "additive"]), 1, 1e-10)
gmoq_check("T0.FGRAD.lognormal_grid_gap_reported",
           "For the lognormal the 45-point grid average of the multiplicative column is reported, not asserted.",
           TRUE, stat = mean(W_ln[, "multiplicative"]) -
             g_ln[["multiplicative"]], severity = "report",
           criterion = "mean(W[, multiplicative]) - grad_multiplicative")

## =====================================================================
## 4. Normal mean/variance vQTL: mean gradient is 1 on beta_mu, 0 on variance
## =====================================================================
mu0 <- 2.0; sd0 <- 1.5
W_v <- make_weights_vqtl(taus, q_tau = qnorm(taus, mu0, sd0), mu = mu0, sd = sd0)
obj_v <- list(W = W_v, params0 = list(mu = mu0, sigma2 = sd0^2),
              type = "normal_vqtl")
f_mean_norm <- function(p) p$mu
g_v <- functional_grad(obj_v, f_mean_norm)
gmoq_close("T0.FGRAD.normal_vqtl_mu", "The mean's beta_mu gradient is 1.",
           g_v[["beta_mu"]], 1, 1e-8)
gmoq_close("T0.FGRAD.normal_vqtl_sigma2", "The mean does not move with beta_sigma2.",
           g_v[["beta_sigma2"]], 0, 1e-8)

## =====================================================================
## 5. Delta-method standard error
## =====================================================================
theta_mix <- c(gamma = 0.1, beta_left = 0.02, beta_right = -0.01)
cov_mix <- diag(c(0.01, 0.02, 0.03))
est <- functional_estimate(obj_mix, f_mean_mix, theta = theta_mix,
                           cov_theta = cov_mix)
g_at <- functional_grad(obj_mix, f_mean_mix, theta = theta_mix)
gmoq_check("T0.FGRAD.estimate_at_theta",
           "functional_estimate evaluates f at the supplied theta, not at the baseline.",
           isTRUE(all.equal(est$estimate, local({
             w_lo <- plogis(qlogis(w_lower) + theta_mix[["gamma"]])
             f_mean_mix(list(weight = c(w_lo, 1 - w_lo),
                             xi = c(comp$xi[1] + theta_mix[["beta_left"]],
                                    comp$xi[2] + theta_mix[["beta_right"]]),
                             omega = comp$omega, alpha = comp$alpha))
           }))),
           criterion = "estimate equals f(map(theta))")
gmoq_close("T0.FGRAD.estimate_grad_at_theta",
           "The returned gradient is taken at theta, not at the baseline.",
           est$grad, g_at, 1e-12)
gmoq_check("T0.FGRAD.delta_se_finite_positive",
           "The delta-method SE is finite and strictly positive.",
           isTRUE(is.finite(est$se) && est$se > 0),
           stat = est$se, criterion = "0 < se < Inf")
gmoq_close("T0.FGRAD.delta_se_formula",
           "The delta-method SE is sqrt(g' Cov g).",
           est$se, sqrt(as.numeric(t(g_at) %*% cov_mix %*% g_at)), 1e-12)
est_nocov <- functional_estimate(obj_mix, f_mean_mix, theta = theta_mix)
gmoq_check("T0.FGRAD.delta_se_na_without_cov",
           "Without a covariance the SE is NA rather than a wrong number.",
           isTRUE(is.na(est_nocov$se)),
           criterion = "se is NA when cov_theta is NULL")

## =====================================================================
## 6. param_map override supports any other family
## =====================================================================
obj_pm <- list(W = matrix(0, 1, 1, dimnames = list(NULL, "theta")),
               param_map = function(theta) list(x = 2 + theta[["theta"]]))
g_pm <- functional_grad(obj_pm, function(p) p$x^2)
gmoq_close("T0.FGRAD.param_map_override",
           "An object carrying param_map uses it instead of the built-in families; d/dtheta (2+theta)^2 at 0 is 4.",
           g_pm[["theta"]], 4, 1e-6)

gmoq_finish()