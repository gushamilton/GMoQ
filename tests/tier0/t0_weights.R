#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 0 -- the weight builders against hand-computed W.
#
# A weight column is only correct if it is the derivative of the quantile
# function with respect to its parameter:
#
#     W_j(tau) = dq_tau/dtheta_j = -(dF/dtheta_j)(q_tau) / f(q_tau)
#
# So every column here is checked against a CENTRAL DIFFERENCE OF THE ACTUAL
# QUANTILE FUNCTION, not against a re-typed copy of the algebra. That is the
# difference between testing the code and testing the derivation: a sign error
# in the published formula would pass the second and fail the first.
#
# Run:  Rscript --vanilla tests/tier0/t0_weights.R
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier0/t0_weights.R")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
source(file.path(GMOQ_SUITE_ROOT, "lib", "models.R"))
gmoq_begin("tier0", "weights")

taus <- c(0.02, 0.05, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.95, 0.98)

# A central difference of a quantile function supplied as qfun(tau, params).
fd_dq <- function(qfun, params, name, taus, eps = 1e-5) {
  up <- params; up[[name]] <- up[[name]] + eps
  dn <- params; dn[[name]] <- dn[[name]] - eps
  (qfun(taus, up) - qfun(taus, dn)) / (2 * eps)
}

# =====================================================================
# 1. make_weights_vqtl: Normal mean and VARIANCE
# =====================================================================
mu0 <- 2.0; sd0 <- 1.5
q_norm <- function(taus, p) stats::qnorm(taus, mean = p$mu, sd = sqrt(p$v))
W <- make_weights_vqtl(taus, q_tau = stats::qnorm(taus, mu0, sd0), mu = mu0, sd = sd0)

gmoq_check("T0.W.vqtl_colnames",
           "make_weights_vqtl returns the documented column names in the documented order.",
           identical(colnames(W), c("beta_mu", "beta_sigma2")),
           criterion = "colnames == c(beta_mu, beta_sigma2)")

gmoq_close("T0.W.vqtl_mu",
           "The mean column is dq/dmu = 1, verified by differencing qnorm().",
           W[, "beta_mu"], fd_dq(q_norm, list(mu = mu0, v = sd0^2), "mu", taus), 1e-8)

gmoq_close("T0.W.vqtl_sigma2",
           "The variance column is dq/dsigma^2 = z/(2 sigma), verified by differencing qnorm().",
           W[, "beta_sigma2"], fd_dq(q_norm, list(mu = mu0, v = sd0^2), "v", taus), 1e-7)

gmoq_refuses("T0.W.vqtl_rejects_bad_sd",
             "A non-positive baseline sd is refused rather than producing an infinite column.",
             make_weights_vqtl(taus, stats::qnorm(taus), mu = 0, sd = 0), "sd must be positive")
gmoq_refuses("T0.W.vqtl_length_mismatch",
             "A q_tau of the wrong length is refused rather than recycled.",
             make_weights_vqtl(taus, stats::qnorm(taus)[1:3], mu = 0, sd = 1), ".")

# =====================================================================
# 2. make_weights_normal_mixture_vqtl: five columns of a 2-Normal mixture
#
# The baseline is Y ~ p1 N(mu1, sd1^2) + (1-p1) N(mu2, sd2^2), and the
# membership parameter is the LOG-ODDS of component 1, so p1 = plogis(gamma).
# The quantile function has no closed form, so it is obtained by root-finding
# on the mixture CDF and then differenced.
# =====================================================================
MP <- list(gamma = stats::qlogis(0.4), mu1 = 1.2, sd1 = 0.6, mu2 = 3.0, sd2 = 0.9)

mix_cdf <- function(y, p) {
  p1 <- stats::plogis(p$gamma)
  p1 * stats::pnorm(y, p$mu1, p$sd1) + (1 - p1) * stats::pnorm(y, p$mu2, p$sd2)
}
mix_q <- function(taus, p) {
  vapply(taus, function(tt)
    stats::uniroot(function(y) mix_cdf(y, p) - tt, interval = c(-40, 60),
                   tol = .Machine$double.eps^0.75)$root,
    numeric(1))
}

q_mix <- mix_q(taus, MP)
Wm <- make_weights_normal_mixture_vqtl(
  taus, q_tau = q_mix, p1 = stats::plogis(MP$gamma),
  mu1 = MP$mu1, sd1 = MP$sd1, mu2 = MP$mu2, sd2 = MP$sd2,
  density = "model")   # the closed form below is the model-density W

gmoq_check("T0.W.mix_colnames",
           "The mixture builder returns gamma, beta_mu1, beta_mu2, beta_sigma1, beta_sigma2.",
           identical(colnames(Wm),
                     c("gamma", "beta_mu1", "beta_mu2", "beta_sigma1", "beta_sigma2")),
           criterion = "documented column names and order")

for (pair in list(c("gamma", "gamma"), c("beta_mu1", "mu1"), c("beta_mu2", "mu2"),
                  c("beta_sigma1", "sd1"), c("beta_sigma2", "sd2"))) {
  gmoq_close(paste0("T0.W.mix_", pair[1]),
             sprintf("Mixture column %s is dq/d%s, verified against the root-found quantile function.",
                     pair[1], pair[2]),
             Wm[, pair[1]], fd_dq(mix_q, MP, pair[2], taus, eps = 1e-5), 2e-5)
}

# The membership-only builder must agree with the five-column one on the
# columns they share -- they are documented as the same weights.
Wm2 <- make_weights_normal_mixture(
  taus, q_tau = q_mix, p1 = stats::plogis(MP$gamma),
  mu1 = MP$mu1, sd1 = MP$sd1, mu2 = MP$mu2, sd2 = MP$sd2,
  include_membership = TRUE, density = "model")
gmoq_close("T0.W.mix_two_builders_agree",
           "make_weights_normal_mixture and ..._vqtl agree on the shared columns.",
           as.numeric(Wm2[, c("gamma", "beta_1", "beta_2")]),
           as.numeric(Wm[, c("gamma", "beta_mu1", "beta_mu2")]), 1e-12)

# A structural identity worth pinning: the two mean columns sum to 1 at every
# tau, because p1 f1 + (1-p1) f2 = f. docs/RELATEDNESS.md leans on exactly this
# to show that additive relatedness loads entirely onto the mean parameter.
gmoq_close("T0.W.mix_mean_columns_sum_to_one",
           "beta_mu1 + beta_mu2 = 1 at every tau, the identity the relatedness argument uses.",
           Wm[, "beta_mu1"] + Wm[, "beta_mu2"], rep(1, length(taus)), 1e-12)

# =====================================================================
# 3. make_weights_generic and make_fd_grad
# =====================================================================
dist_cdf <- function(y, p) stats::plnorm(y, meanlog = p$meanlog, sdlog = p$sdlog)
dist_pdf <- function(y, p) stats::dlnorm(y, meanlog = p$meanlog, sdlog = p$sdlog)
LP <- list(meanlog = 0.4, sdlog = 0.7)
q_ln <- stats::qlnorm(taus, LP$meanlog, LP$sdlog)

# analytic dF/dtheta for the lognormal
g_meanlog <- function(y, p) -stats::dnorm((log(y) - p$meanlog) / p$sdlog) / p$sdlog
g_sdlog   <- function(y, p) { z <- (log(y) - p$meanlog) / p$sdlog
                              -z * stats::dnorm(z) / p$sdlog }

Wg <- make_weights_generic(taus, q_ln, dist_cdf, dist_pdf, params = LP,
                           grad_funcs = list(beta_meanlog = g_meanlog, beta_sdlog = g_sdlog))

q_lnf <- function(taus, p) stats::qlnorm(taus, p$meanlog, p$sdlog)
gmoq_close("T0.W.generic_meanlog",
           "The generic builder's meanlog column is dq/dmeanlog of the lognormal quantile.",
           Wg[, "beta_meanlog"], fd_dq(q_lnf, LP, "meanlog", taus), 1e-6)
gmoq_close("T0.W.generic_sdlog",
           "The generic builder's sdlog column is dq/dsdlog of the lognormal quantile.",
           Wg[, "beta_sdlog"], fd_dq(q_lnf, LP, "sdlog", taus), 1e-6)

# The published identity W = -dF/dtheta / f, written out directly.
gmoq_close("T0.W.generic_identity",
           "make_weights_generic implements W_j = -(dF/dtheta_j)/f exactly.",
           Wg[, "beta_meanlog"], -g_meanlog(q_ln, LP) / dist_pdf(q_ln, LP), 1e-12)

# make_fd_grad must reproduce the analytic gradient.
fdg <- make_fd_grad("sdlog", dist_cdf)
gmoq_close("T0.W.fd_grad",
           "make_fd_grad's central difference matches the analytic dF/dtheta.",
           vapply(q_ln, function(q) fdg(q, LP), numeric(1)), g_sdlog(q_ln, LP), 1e-6)

# =====================================================================
# 4. The suite's own sinh-arcsinh columns (tests/lib/models.R)
#
# These generate the truth for tiers 1-3, so they are tested to the same
# standard as the shipped builders: against a difference of the sinh-arcsinh
# quantile function itself. Checked at a SKEWED, HEAVY-TAILED baseline as well
# as at the symmetric one, because tier 2 runs cells at both.
# =====================================================================
sig <- 1.3
sas_q <- function(taus, p) {
  z <- stats::qnorm(taus)
  p$mu + p$sigma * sinh((asinh(z) + p$e) / exp(p$ld))
}
for (base in list(list(mu = 0, sigma = sig, e = 0, ld = 0),
                  list(mu = 0.4, sigma = sig, e = 0.8, ld = 0.25))) {
  tagb <- if (base$e == 0 && base$ld == 0) "symmetric" else "skewed"
  cols <- gmoq_cf_columns(taus, base$sigma, e0 = base$e, d0 = base$ld)
  for (pair in list(c("beta_mu", "mu"), c("beta_sigma", "sigma"),
                    c("beta_skew", "e"), c("beta_kurt", "ld"))) {
    gmoq_close(paste0("T0.W.sas_", pair[1], ".", tagb),
               sprintf("Suite truth column %s is d q/d %s of the sinh-arcsinh quantile function (%s baseline).",
                       pair[1], pair[2], tagb),
               cols[[pair[1]]], fd_dq(sas_q, base, pair[2], taus), 1e-6)
  }
}
cols <- gmoq_cf_columns(taus, sig)

# The shape columns must not be in the span of location and scale, or the
# tournament could never tell them apart even in principle.
Wls <- cbind(cols$beta_mu, cols$beta_sigma)
Mls <- diag(length(taus)) - Wls %*% solve(crossprod(Wls), t(Wls))
res_frac <- vapply(c("beta_skew", "beta_kurt"), function(nm)
  sum((Mls %*% cols[[nm]])^2) / sum(cols[[nm]]^2), numeric(1))
gmoq_check("T0.W.shape_columns_are_identifiable",
           "The skew and tail columns are not IN the span of location and scale, so the models are distinguishable in principle.",
           all(res_frac > 1e-3), stat = min(res_frac),
           criterion = "each shape column keeps more than 0.1 per cent of its norm after projection onto [1, sigma]")
gmoq_check("T0.W.shape_columns_are_weakly_identified",
           "How little of each shape column survives projection onto location and scale: the quantitative reason the tournament cannot separate shape models.",
           TRUE, severity = "report", stat = min(res_frac),
           criterion = "residual norm fraction after projecting out [1, sigma]",
           detail = paste(sprintf("%s=%.4f", names(res_frac), res_frac), collapse = " "))

# The model set must be full rank -- a rank-deficient candidate makes
# solve(crossprod(W)) fail inside param_gwas and the tournament meaningless.
ms <- gmoq_model_set(taus, sig,
                     mix = list(q_tau = q_mix, p1 = stats::plogis(MP$gamma),
                                mu1 = MP$mu1, sd1 = MP$sd1, mu2 = MP$mu2, sd2 = MP$sd2))
ranks <- vapply(ms, function(M) qr(M)$rank == ncol(M), logical(1))
gmoq_check("T0.W.model_set_full_rank",
           "Every candidate model in the tournament has full column rank on the tau grid.",
           all(ranks), criterion = "rank(W) == ncol(W) for all candidates",
           detail = if (!all(ranks)) paste("deficient:", paste(names(ranks)[!ranks], collapse = ", ")) else NULL)

gmoq_check("T0.W.model_set_covers_families",
           "Every simulated truth family names a candidate model that exists in the set.",
           all(vapply(GMOQ_TRUTH_MAP, function(e)
             e$true %in% names(ms) && all(e$wrong %in% names(ms)), logical(1))),
           criterion = "truth map is consistent with the model set")

gmoq_finish()
