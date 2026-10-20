# Weight builders: each column of W is the derivative of the quantile function
# with respect to one model parameter, divided by the chosen density.

taus <- seq(0.02, 0.98, length.out = 45)
P <- list(p1 = 0.45, mu1 = -1.8, sd1 = 0.6, mu2 = 1.5, sd2 = 0.8)
mix_q <- function(tt, p1 = P$p1, mu1 = P$mu1, mu2 = P$mu2, sd1 = P$sd1, sd2 = P$sd2) {
  cdf <- function(x) p1 * pnorm(x, mu1, sd1) + (1 - p1) * pnorm(x, mu2, sd2)
  vapply(tt, function(t) uniroot(function(x) cdf(x) - t, c(-15, 15), tol = 1e-13)$root, 0)
}
q <- mix_q(taus)
f_model <- P$p1 * dnorm(q, P$mu1, P$sd1) + (1 - P$p1) * dnorm(q, P$mu2, P$sd2)
# A Stage 1 density that differs from the model, most in the gap.
fhat <- f_model * (1 + 0.8 * exp(-q^2))

mix <- function(...) make_weights_normal_mixture(taus, q, P$p1, P$mu1, P$sd1,
                                                 P$mu2, P$sd2, include_membership = TRUE, ...)

test_that("Normal mean/variance weights are the quantile derivatives", {
  mu <- 0.3; s2 <- 1.7; z <- qnorm(taus); eps <- 1e-6
  qf <- function(mu, s2) mu + sqrt(s2) * z
  W <- make_weights_vqtl(taus, qf(mu, s2), mu = mu, sd = sqrt(s2))
  expect_identical(colnames(W), c("beta_mu", "beta_sigma2"))
  expect_equal(unname(W[, "beta_mu"]), (qf(mu + eps, s2) - qf(mu - eps, s2)) / (2 * eps),
               tolerance = 1e-8)
  expect_equal(unname(W[, "beta_sigma2"]), (qf(mu, s2 + eps) - qf(mu, s2 - eps)) / (2 * eps),
               tolerance = 1e-7)
  expect_error(make_weights_vqtl(taus, qf(mu, s2), mu = mu, sd = 0), "positive")
  expect_silent(check_vqtl_design(W, sd = sqrt(s2)))
  # a design parameterised in the SD, not the variance, is refused
  expect_error(check_vqtl_design(cbind(beta_mu = 1, beta_sigma2 = z / sqrt(s2)),
                                 sd = sqrt(s2), taus = taus))
})

test_that("mixture weights (model density) are the quantile derivatives", {
  W <- mix(density = "model")
  eps <- 1e-5
  d <- function(...) {
    a <- list(...); up <- lo <- a
    up[[1]] <- a[[1]] + eps; lo[[1]] <- a[[1]] - eps
    nm <- names(a)[1]
    (do.call(mix_q, c(list(taus), setNames(up, nm))) -
       do.call(mix_q, c(list(taus), setNames(lo, nm)))) / (2 * eps)
  }
  expect_equal(unname(W[, "beta_1"]), d(mu1 = P$mu1), tolerance = 1e-6)
  expect_equal(unname(W[, "beta_2"]), d(mu2 = P$mu2), tolerance = 1e-6)
  # gamma acts on the log-odds of membership of component 1
  qg <- function(g) mix_q(taus, p1 = plogis(qlogis(P$p1) + g))
  expect_equal(unname(W[, "gamma"]), (qg(eps) - qg(-eps)) / (2 * eps), tolerance = 1e-6)
  # a shift of both components moves every quantile by the same amount
  expect_lt(max(abs(W[, "beta_1"] + W[, "beta_2"] - 1)), 1e-12)
})

test_that("the default is hybrid with a membership column and model without (0.6.0)", {
  expect_identical(mix(stage1_density = fhat), mix(density = "hybrid", stage1_density = fhat))
  expect_error(mix(), "default whenever W has a membership column")
  expect_identical(make_weights_normal_mixture(taus, q, P$p1, P$mu1, P$sd1, P$mu2, P$sd2),
                   make_weights_normal_mixture(taus, q, P$p1, P$mu1, P$sd1, P$mu2, P$sd2,
                                               density = "model"))
  expect_equal(weight_density_mode(
    make_weights_normal_mixture_vqtl(taus, q, P$p1, P$mu1, P$sd1, P$mu2, P$sd2,
                                     stage1_density = fhat)), "hybrid")
  expect_error(make_weights_normal_mixture_vqtl(taus, q, P$p1, P$mu1, P$sd1, P$mu2, P$sd2),
               "needs `stage1_density`")
  expect_error(mix(density = "bogus"), "must be one of")
  # Generic: hybrid only when a gradient is a membership column.
  cdf <- function(x, p) pnorm(x, p$mu, p$sd); pdf <- function(x, p) dnorm(x, p$mu, p$sd)
  Wn <- make_weights_generic(taus, q, cdf, pdf, list(mu = 0, sd = 1),
                             grad_funcs = list(beta_mu = make_fd_grad("mu", cdf)))
  expect_equal(weight_density_mode(Wn), "model")
})

test_that("hybrid gamma column is -dF/dgamma divided by the Stage 1 density", {
  W <- mix(density = "hybrid", stage1_density = fhat)
  dF_dgamma <- P$p1 * (1 - P$p1) * (pnorm(q, P$mu1, P$sd1) - pnorm(q, P$mu2, P$sd2))
  expect_equal(unname(W[, "gamma"]), -dF_dgamma / fhat, tolerance = 1e-14)
  W0 <- mix(density = "model")
  expect_identical(unname(W[, c("beta_1", "beta_2")]), unname(W0[, c("beta_1", "beta_2")]))
  expect_equal(weight_density_mode(W), "hybrid")
  expect_equal(weight_density_mode(W, label = TRUE), "hybrid_stage1_membership_density")
  expect_equal(attr(W, "stage1_density_columns"), "gamma")
  expect_equal(attr(W, "stage1_density"), fhat)
  expect_equal(attr(W, "model_density"), pmax(f_model, 1e-12))

  Wv <- make_weights_normal_mixture_vqtl(taus, q, P$p1, P$mu1, P$sd1, P$mu2, P$sd2,
                                         density = "hybrid", stage1_density = fhat)
  expect_equal(unname(Wv[, "gamma"]), -dF_dgamma / fhat, tolerance = 1e-14)
  Wv0 <- make_weights_normal_mixture_vqtl(taus, q, P$p1, P$mu1, P$sd1, P$mu2, P$sd2,
                                          density = "model")
  expect_identical(unname(Wv[, -1]), unname(Wv0[, -1]))

  # Generic builder with a finite-difference membership gradient.
  cdf <- function(x, p) p$p1 * pnorm(x, p$mu1, p$sd1) + (1 - p$p1) * pnorm(x, p$mu2, p$sd2)
  pdf <- function(x, p) p$p1 * dnorm(x, p$mu1, p$sd1) + (1 - p$p1) * dnorm(x, p$mu2, p$sd2)
  gp <- make_fd_grad("p1", cdf)
  gg <- function(x, p) gp(x, p) * p$p1 * (1 - p$p1)
  Wg <- make_weights_generic(taus, q, cdf, pdf, P,
                             grad_funcs = list(gamma = gg,
                                               beta_left = make_fd_grad("mu1", cdf),
                                               beta_right = make_fd_grad("mu2", cdf)),
                             density = "hybrid", stage1_density = fhat,
                             unit_sum = c("beta_left", "beta_right"))
  dF <- vapply(q, function(x) gg(x, P), numeric(1))
  expect_identical(unname(Wg[, "gamma"]), -dF / fhat)
  expect_equal(unname(Wg[, -1]), unname(W[, -1]), tolerance = 1e-6)
  expect_equal(unname(Wg[, "gamma"]), -dF_dgamma / fhat, tolerance = 1e-6)
})

test_that("identities hold under model and hybrid; stage1 breaks them and warns", {
  for (W in list(mix(density = "model"), mix(density = "hybrid", stage1_density = fhat))) {
    expect_lt(max(abs(W[, "beta_1"] + W[, "beta_2"] - 1)), 1e-12)
  }
  expect_warning(Ws <- mix(density = "stage1", stage1_density = fhat), "sensitivity")
  expect_gt(max(abs(Ws[, "beta_1"] + Ws[, "beta_2"] - 1)), 0.1)
  expect_equal(weight_density_mode(Ws), "stage1")
  expect_equal(attr(Ws, "stage1_density_columns"), c("gamma", "beta_1", "beta_2"))

  # Shifted lognormal: the additive column is exactly 1 under model and hybrid.
  ln_cdf <- function(x, p) plnorm(x - p$shift, p$meanlog, p$sdlog)
  ln_pdf <- function(x, p) dlnorm(x - p$shift, p$meanlog, p$sdlog)
  pl <- list(shift = -4, meanlog = 1.4, sdlog = 0.35)
  gf <- list(additive = function(x, p) -ln_pdf(x, p),
             multiplicative = make_fd_grad("meanlog", ln_cdf))
  Wm <- make_weights_generic(taus, q, ln_cdf, ln_pdf, pl, gf, unit_sum = "additive")
  expect_true(all(Wm[, "additive"] == 1))
  Wh <- make_weights_generic(taus, q, ln_cdf, ln_pdf, pl, gf, density = "hybrid",
                             stage1_density = fhat, stage1_columns = character(0),
                             unit_sum = "additive")
  expect_identical(unname(Wh[, ]), unname(Wm[, ]))
  expect_equal(weight_density_mode(Wh), "hybrid")
  expect_warning(make_weights_generic(taus, q, ln_cdf, ln_pdf, pl, gf, density = "stage1",
                                      stage1_density = fhat, unit_sum = "additive"),
                 "additive")
  # A hybrid call that names no membership column is a mistake, not a no-op.
  expect_error(make_weights_generic(taus, q, ln_cdf, ln_pdf, pl, gf, density = "hybrid",
                                    stage1_density = fhat), "membership column")
  # A violated identity stops under model / hybrid.
  expect_error(make_weights_generic(taus, q, ln_cdf, ln_pdf, pl, gf,
                                    unit_sum = "multiplicative"), "identity violated")
})

test_that("the Stage 1 density is required, checked and readable from JSON", {
  expect_error(mix(density = "hybrid"), "needs `stage1_density`")
  expect_error(mix(density = "hybrid", stage1_density = fhat[-1]), "length")
  expect_error(mix(density = "hybrid", stage1_density = replace(fhat, 3, -1)), "positive")
  expect_error(mix(density = "hybrid",
                   stage1_density = list(fhat = fhat, taus = taus + 0.01)), "tau grid")
  expect_error(mix(density = "hybrid",
                   stage1_density = list(fhat = fhat, q_tau = q + 0.1)), "q_tau")
  js <- tempfile(fileext = ".json")
  fmt <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
  writeLines(sprintf('{"n": 1000, "mean": 0.1, "sd": 1.7, "taus": [%s], "q_tau": [%s], "fhat": [%s]}',
                     fmt(taus), fmt(q), fmt(fhat)), js)
  s1 <- read_stage1_density(js)
  expect_equal(s1$fhat, fhat, tolerance = 1e-15)
  W1 <- mix(density = "hybrid", stage1_density = js)
  W2 <- mix(density = "hybrid", stage1_density = fhat)
  expect_equal(unname(W1[, ]), unname(W2[, ]), tolerance = 1e-15)
  expect_equal(attr(W1, "stage1_density_source"), normalizePath(js))
})

test_that("stage1_kde_density reproduces the Stage 1 Silverman KDE", {
  set.seed(11)
  y <- c(rnorm(700, -2, 0.7), rnorm(300, 2, 1))
  qq <- as.numeric(quantile(y, c(0.1, 0.5, 0.9), type = 8))
  h <- 1.06 * sqrt(mean((y - mean(y))^2)) * length(y)^(-0.2)
  ref <- vapply(qq, function(x) mean(exp(-0.5 * ((y - x) / h)^2)) / (h * sqrt(2 * pi)), 0)
  expect_equal(stage1_kde_density(y, qq), ref, tolerance = 1e-14)
  ref4 <- vapply(qq, function(x) { u <- (y - x) / h
    mean(exp(-0.5 * u^2) * (1.5 - 0.5 * u^2)) / (h * sqrt(2 * pi)) }, 0)
  expect_equal(stage1_kde_density(y, qq, "fourth_order_gaussian"), pmax(ref4, 1e-8),
               tolerance = 1e-14)
})

test_that("provenance: weight_density_mode and the Stage 2 gate", {
  W0 <- mix(density = "model"); Wh <- mix(density = "hybrid", stage1_density = fhat)
  expect_equal(weight_density_mode(W0), "model")
  expect_equal(weight_density_mode(W0, label = TRUE), "parametric_model_density")
  expect_equal(weight_density_mode(list(W = W0,
                 provenance = list(density_weight_mode = "stage1_rif_density"))), "stage1")
  expect_equal(weight_density_mode(list(W = Wh)), "hybrid")
  chk <- gmoq:::.check_require_density_mode
  expect_silent(chk(W0, NULL))
  expect_silent(chk(W0, "model"))
  expect_silent(chk(Wh, "hybrid_stage1_membership_density"))
  expect_error(chk(W0, "hybrid"), "refusing mismatched weights")
  expect_error(chk(Wh, "bogus"), "must be one of")
  expect_true("require_density_mode" %in% names(formals(param_gwas_from_file)))
  expect_true("require_density_mode" %in% names(formals(param_gwas_multi_from_file)))
})

test_that("hybrid removes the gamma attenuation when the fitted density misses the gap", {
  # Truth: two separated Laplace components, the SNP moves membership only.
  # The analyst fits two Normals with the right weights, means and sds, so the
  # model CDF is close but its density in the gap is ~0.56x the truth. Stage 1
  # is reproduced exactly: RIF on the Silverman KDE, OLS slope on dosage.
  set.seed(20260923)
  N <- 200000; delta <- 0.25; b <- 0.5; loc <- 2
  g <- rbinom(N, 2, 0.3)
  p_low <- plogis(delta * (g - mean(g)))
  lower <- runif(N) < p_low
  lap <- function(n) rexp(n, 1 / b) * sample(c(-1, 1), n, replace = TRUE)
  y <- ifelse(lower, -loc, loc) + lap(N)
  tt <- seq(0.05, 0.95, by = 0.05)
  qt <- as.numeric(quantile(y, tt, type = 8))
  fh <- stage1_kde_density(y, qt)
  gc <- g - mean(g)
  bt <- vapply(seq_along(tt), function(k) {
    rif <- qt[k] + (tt[k] - (y <= qt[k])) / fh[k]
    sum(gc * rif) / sum(gc^2)
  }, numeric(1))

  p1 <- mean(lower); s <- b * sqrt(2)
  Wp <- make_weights_normal_mixture(tt, qt, p1, -loc, s, loc, s, include_membership = TRUE,
                                    density = "model")
  Wh <- make_weights_normal_mixture(tt, qt, p1, -loc, s, loc, s, include_membership = TRUE,
                                    density = "hybrid", stage1_density = fh)
  mid <- which.min(abs(tt - 0.5))
  expect_lt(attr(Wh, "model_density")[mid] / fh[mid], 0.75)   # the gap misfit is real

  ols <- function(W) drop(solve(crossprod(W), crossprod(W, bt)))
  gp <- ols(Wp)[["gamma"]] / delta
  gh <- ols(Wh)[["gamma"]] / delta
  expect_gt(gh, 0.85); expect_lt(gh, 1.15)
  expect_lt(gp, 0.8)
  expect_gt(gh - gp, 0.2)
})
