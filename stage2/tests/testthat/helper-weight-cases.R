# Deterministic weight-builder calls shared by test-weight-density.R and by
# fixtures/make_weights_default_fixture.R, which records their output under the
# PREVIOUS package commit. The density = "model" output of the current package
# must match that record. ("model" was the default before 0.6.0; a W with a
# membership column now defaults to "hybrid", so those cases ask for "model".)

weight_cases_taus <- function() seq(0.02, 0.98, length.out = 45)

weight_cases_q <- function(taus) {
  # A separated two-Normal mixture's quantiles, found by root finding, so the
  # grid is fixed without any random numbers.
  cdf <- function(x) 0.45 * stats::pnorm(x, -1.8, 0.6) + 0.55 * stats::pnorm(x, 1.5, 0.8)
  vapply(taus, function(t) stats::uniroot(function(x) cdf(x) - t, c(-10, 10),
                                          tol = 1e-13)$root, numeric(1))
}

weight_cases_default <- function() {
  taus <- weight_cases_taus()
  q <- weight_cases_q(taus)
  mix_cdf <- function(x, p) p$p1 * stats::pnorm(x, p$mu1, p$sd1) +
    (1 - p$p1) * stats::pnorm(x, p$mu2, p$sd2)
  mix_pdf <- function(x, p) p$p1 * stats::dnorm(x, p$mu1, p$sd1) +
    (1 - p$p1) * stats::dnorm(x, p$mu2, p$sd2)
  pm <- list(p1 = 0.45, mu1 = -1.8, sd1 = 0.6, mu2 = 1.5, sd2 = 0.8)
  grad_p <- make_fd_grad("p1", mix_cdf)
  ln_cdf <- function(x, p) { o <- numeric(length(x)); k <- x > p$shift
    o[k] <- stats::plnorm(x[k] - p$shift, p$meanlog, p$sdlog); o }
  ln_pdf <- function(x, p) { o <- numeric(length(x)); k <- x > p$shift
    o[k] <- stats::dlnorm(x[k] - p$shift, p$meanlog, p$sdlog); o }
  pl <- list(shift = -4, meanlog = 1.4, sdlog = 0.35)
  list(
    mix = make_weights_normal_mixture(taus, q, 0.45, -1.8, 0.6, 1.5, 0.8),
    mix_membership = make_weights_normal_mixture(taus, q, 0.45, -1.8, 0.6, 1.5, 0.8,
                                                 include_membership = TRUE,
                                                 density = "model"),
    mix_vqtl = make_weights_normal_mixture_vqtl(taus, q, 0.45, -1.8, 0.6, 1.5, 0.8,
                                                density = "model"),
    generic_mixture = make_weights_generic(
      taus, q, mix_cdf, mix_pdf, pm,
      grad_funcs = list(
        gamma = function(x, p) grad_p(x, p) * p$p1 * (1 - p$p1),
        beta_left = make_fd_grad("mu1", mix_cdf),
        beta_right = make_fd_grad("mu2", mix_cdf)),
      density = "model"),
    generic_lognormal = make_weights_generic(
      taus, q, ln_cdf, ln_pdf, pl,
      grad_funcs = list(
        additive = function(x, p) -ln_pdf(x, p),
        multiplicative = make_fd_grad("meanlog", ln_cdf)))
  )
}
