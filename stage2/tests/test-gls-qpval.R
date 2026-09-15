# ============================================================
# Test: calibrated lack-of-fit p-value on the GLS path
#
# The GLS Q is NOT the OLS Q and must not borrow its reference.
#
#   OLS:  Q = r' D^-1 r,  r = M b        -> sum_i lambda_i chi2_1 (per-SNP spectrum)
#   GLS:  Q = r' S^-1 r,  r = b - W theta_gls
#           -> chi2_{T-K} when Sigma is known
#           -> (nu m/(nu-m+1)) F_{m, nu-m+1} when Sigma is estimated on nu df
#
# Pins:
#   1. Q_pval / Q_df_eff / Q_scale are present on the GLS path.
#   2. The reference is exact under a CORRECTLY SPECIFIED model: Q_pval is
#      uniform, both with Sigma known and with Sigma estimated at finite nu.
#   3. The chi-square reference that a naive port of the OLS path would give is
#      measurably anti-conservative at finite nu -- i.e. the Hotelling
#      correction is load-bearing, not cosmetic.
#   4. Q_pval falls under misspecification.
#   5. The OLS path is UNCHANGED (no new fields, p-values still reproduce
#      .q_pvalue_calibrated exactly).
#   6. Every inapplicable combination refuses loudly rather than returning a
#      number from the wrong family.
# ============================================================

suppressPackageStartupMessages(library(Rcpp))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg) else "tests/test-gls-qpval.R"
pkg_root <- normalizePath(file.path(dirname(this_file), ".."))

Rcpp::sourceCpp(file.path(pkg_root, "src", "stage2_se_helper.cpp"))
source(file.path(pkg_root, "R", "utils.R"))
source(file.path(pkg_root, "R", "covariance.R"))
source(file.path(pkg_root, "R", "weight-functions.R"))
source(file.path(pkg_root, "R", "tail-reference.R"))
source(file.path(pkg_root, "R", "gwas.R"))

ok <- function(cond, msg) {
  if (isTRUE(cond)) cat("PASS:", msg, "\n") else stop("FAIL: ", msg)
}
quiet <- function(expr) suppressMessages(suppressWarnings(expr))

# ---- fixtures ---------------------------------------------------------
# Brownian-bridge cross-quantile covariance: the structure Stage 1 produces.
Tt <- 45L
taus <- seq(0.02, 0.98, length.out = Tt)
fdens <- dnorm(qnorm(taus))
Sigma <- outer(taus, taus, function(a, b) pmin(a, b) - a * b) / outer(fdens, fdens)
Sigma <- Sigma + diag(1e-8 * mean(diag(Sigma)), Tt)
cholS <- t(chol(Sigma))
q_tau <- qnorm(taus)
W <- cbind(intercept = 1, scale = q_tau - mean(q_tau))
K <- ncol(W)
m <- Tt - K
theta <- c(0.30, 0.12)

upper_row_major <- function(S) {
  out <- numeric(nrow(S) * (nrow(S) + 1) / 2); idx <- 1L
  for (r in seq_len(nrow(S))) for (cc in r:ncol(S)) { out[idx] <- S[r, cc]; idx <- idx + 1L }
  out
}

# nu * V ~ Wishart_T(nu, Sigma) via Bartlett: O(T^2), not O(nu T).
rwish <- function(nu) {
  A <- matrix(0, Tt, Tt)
  A[lower.tri(A)] <- rnorm(Tt * (Tt - 1) / 2)
  diag(A) <- sqrt(rchisq(Tt, df = nu - seq_len(Tt) + 1))
  L <- cholS %*% A
  tcrossprod(L) / nu
}

# P null variants. nu = Inf supplies Sigma exactly; `resid` adds a fixed
# systematic residual (misspecification).
make_data <- function(P, nu = Inf, resid = NULL) {
  mu <- as.numeric(W %*% theta)
  if (!is.null(resid)) mu <- mu + resid
  ncov <- Tt * (Tt + 1) / 2
  b <- matrix(NA_real_, Tt, P); se <- matrix(NA_real_, Tt, P)
  cv <- numeric(P * ncov)
  for (j in seq_len(P)) {
    b[, j] <- mu + cholS %*% rnorm(Tt)
    S <- if (is.finite(nu)) rwish(nu) else Sigma
    cv[((j - 1) * ncov + 1):(j * ncov)] <- upper_row_major(S)
    se[, j] <- sqrt(diag(S))
  }
  colnames(b) <- paste0("snp", seq_len(P))
  list(stage1 = list(taus = taus, q_tau = q_tau, Q_slope = b, SE_tau = se),
       cov_data = list(cov_vec = cv, offsets = seq(0, by = ncov * 4, length.out = P)))
}

fit_gls <- function(d, ...) quiet(param_gwas(
  d$stage1, transform = "custom_W", transform_args = list(W = W),
  se_mode = "tau_cov", cov_data = d$cov_data, tau_cov_estimator = "gls", ...))

# ---- 1. the three fields exist on the GLS path ------------------------
set.seed(20260915)
d0 <- make_data(2000L, nu = 400)
g0 <- fit_gls(d0, nu = 400)
ok(!is.null(g0$Q_pval),   "GLS path returns Q_pval")
ok(!is.null(g0$Q_df_eff), "GLS path returns Q_df_eff")
ok(!is.null(g0$Q_scale),  "GLS path returns Q_scale")
ok(identical(g0$Q_ref_dist, "hotelling_f"), "GLS reference is recorded as hotelling_f")
ok(all(g0$Q_df_eff == m),
   "GLS Q_df_eff is exactly T - K: the GLS quadratic form has a flat spectrum")
ok(abs(g0$Q_scale[1] - 400 / (400 - m - 1)) < 1e-12,
   "GLS Q_scale is the Hotelling null-mean inflation nu/(nu - m - 1)")
ok(abs(mean(g0$Q) / (g0$Q_df_eff[1] * g0$Q_scale[1]) - 1) < 0.02,
   "Q_df_eff * Q_scale is the null mean of Q (calibrated ratio has mean 1)")

# ---- 2. uniform under correct specification ---------------------------
# (a) Sigma known: Q ~ chi2_{T-K} exactly.
set.seed(11)
dk <- make_data(4000L, nu = Inf)
gk <- fit_gls(dk, nu = 1e9)
ok(abs(mean(gk$Q) - m) < 0.6, "Sigma known: mean(Q) is T - K")
ok(abs(var(gk$Q) - 2 * m) < 8, "Sigma known: var(Q) is 2(T - K)")
ok(suppressWarnings(ks.test(gk$Q_pval, "punif")$p.value) > 0.01,
   "Sigma known: Q_pval is uniform (KS)")
ok(abs(mean(gk$Q_pval < 0.05) - 0.05) < 0.012,
   "Sigma known: rejection at 0.05 is nominal")

# (b) Sigma estimated at a finite nu: the Hotelling reference is still uniform.
set.seed(12)
nu_f <- 120
de <- make_data(4000L, nu = nu_f)
ge <- fit_gls(de, nu = nu_f)
ok(suppressWarnings(ks.test(ge$Q_pval, "punif")$p.value) > 0.01,
   "Sigma estimated: Q_pval is uniform (KS)")
ok(abs(mean(ge$Q_pval < 0.05) - 0.05) < 0.012,
   "Sigma estimated: rejection at 0.05 is nominal")
ok(abs(mean(ge$Q_pval < 0.20) - 0.20) < 0.022,
   "Sigma estimated: rejection at 0.20 is nominal")

# ---- 3. the correction is load-bearing --------------------------------
# What a naive port of the OLS path would emit: chi2 on T - K.
p_naive <- pchisq(ge$Q, df = m, lower.tail = FALSE)
ok(mean(p_naive < 0.05) > 0.25,
   paste0("chi2(T - K) on the GLS Q is grossly anti-conservative at nu = ",
          nu_f, " (measured ", round(mean(p_naive < 0.05), 3), " at nominal 0.05)"))
ok(mean(p_naive < 0.05) > 5 * mean(ge$Q_pval < 0.05),
   "the Hotelling correction changes the answer by more than a factor of five")

# ---- 4. power under misspecification ----------------------------------
# A skew-shaped residual, which is what a scale effect on a SKEWED error
# actually produces (docs/N_SCALING.md arm A, gamma2: the symmetric basis
# {1, qnorm(tau)} cannot represent it). Scaled to a target non-centrality
# rather than to a bare constant, so the test does not depend on the arbitrary
# overall scale of Sigma: under Sigma known, Q ~ noncentral chi2_m(lambda) with
#
#   lambda = r' (Sinv - Sinv W (W'Sinv W)^-1 W' Sinv) r
#
# so E[Q] = m + lambda and the test can assert against a number it chose.
set.seed(13)
shape <- (q_tau^2 - mean(q_tau^2)) / sd(q_tau^2)
Sinv <- solve(Sigma)
SinvW <- Sinv %*% W
Mid <- Sinv - SinvW %*% solve(crossprod(W, SinvW), t(SinvW))
lambda_target <- 60
lam1 <- as.numeric(crossprod(shape, Mid %*% shape))
stopifnot(lam1 > 0)
skew_resid <- shape * sqrt(lambda_target / lam1)
ok(abs(as.numeric(crossprod(skew_resid, Mid %*% skew_resid)) - lambda_target) < 1e-6,
   "misspecification residual is scaled to the intended non-centrality")
dm <- make_data(2000L, nu = 400, resid = skew_resid)
gm <- fit_gls(dm, nu = 400)
ok(abs(mean(gm$Q) - (m + lambda_target) * gm$Q_scale[1]) < 0.12 * (m + lambda_target),
   "misspecified: mean(Q) matches the non-central prediction (m + lambda)")
ok(mean(gm$Q_pval < 0.05) > 0.90,
   paste0("misspecified: Q_pval rejects at 0.05 (measured ",
          round(mean(gm$Q_pval < 0.05), 3), ")"))
ok(median(gm$Q_pval) < median(g0$Q_pval),
   "misspecified Q_pval is stochastically smaller than the null one")

# ---- 5. the OLS path is unchanged -------------------------------------
o <- quiet(param_gwas(d0$stage1, transform = "custom_W", transform_args = list(W = W),
                      se_mode = "tau_cov", cov_data = d0$cov_data,
                      tau_cov_estimator = "ols"))
ok(is.null(o$Q_ref_dist) && is.null(o$Q_ref_df) && is.null(o$Q_ref_source),
   "OLS output gains no new fields")
ok(!is.null(o$Q_pval) && !is.null(o$Q_df_eff) && !is.null(o$Q_scale),
   "OLS still returns its own calibrated triple")
# The OLS p-value still comes from the eigenvalue-spectrum machinery, untouched.
ref <- quiet(.q_pvalue_calibrated(Q = o$Q, cov_data = d0$cov_data,
                                  SE_tau = d0$stage1$SE_tau, W = W,
                                  method = "saddlepoint", n_threads = 1L))
ok(identical(o$Q_pval, ref$pval) && identical(o$Q_df_eff, ref$df_eff) &&
     identical(o$Q_scale, ref$scale),
   "OLS Q_pval is bit-identical to .q_pvalue_calibrated (reference untouched)")
ok(median(o$Q_df_eff) < m,
   "OLS effective df is still well below T - K (its spectrum is not flat)")
ok(mean(o$Q) < 0.6 * mean(g0$Q),
   "the two statistics are on different scales, as documented")

# ---- 6. refusals ------------------------------------------------------
errs <- function(expr) inherits(try(quiet(expr), silent = TRUE), "try-error")
ok(errs(param_gwas(d0$stage1, transform = "custom_W", transform_args = list(W = W),
                   se_mode = "tau_cov", cov_data = d0$cov_data,
                   tau_cov_estimator = "gls", q_method = "liu", nu = 400)),
   "q_method = 'liu' refuses on the GLS path")
ok(errs(param_gwas(d0$stage1, transform = "custom_W", transform_args = list(W = W),
                   se_mode = "tau_cov", cov_data = d0$cov_data,
                   tau_cov_estimator = "gls", q_method = "satterthwaite", nu = 400)),
   "q_method = 'satterthwaite' refuses on the GLS path")
ok(errs(param_gwas(d0$stage1, transform = "custom_W", transform_args = list(W = W),
                   se_mode = "tau_cov", cov_data = d0$cov_data,
                   tau_cov_estimator = "ols", q_method = "hotelling")),
   "q_method = 'hotelling' refuses on the OLS path")

# nu below m + 1: the residual block of S is singular, Q is degenerate.
gsing <- fit_gls(d0, nu = m - 5)
ok(all(is.na(gsing$Q_pval)),
   "nu below m + 1 yields no Q p-value (Q is degenerate, not just mis-referenced)")
got_warning <- FALSE
invisible(withCallingHandlers(
  suppressMessages(param_gwas(d0$stage1, transform = "custom_W",
                              transform_args = list(W = W), se_mode = "tau_cov",
                              cov_data = d0$cov_data, tau_cov_estimator = "gls",
                              nu = m - 5)),
  warning = function(w) { got_warning <<- TRUE; invokeRestart("muffleWarning") }))
ok(got_warning, "nu below m + 1 warns rather than failing silently")

# nu unknown: omit by default, chi2 only on explicit request.
gu <- fit_gls(d0)
ok(!is.null(gu$Q) && is.null(gu$Q_pval),
   "nu unknown: Q is still returned, Q_pval is omitted (unknown_nu = 'omit')")
gn <- fit_gls(d0, unknown_nu = "normal")
ok(identical(gn$Q_ref_dist, "chisq"),
   "unknown_nu = 'normal' emits the Sigma-known chi2 limit and records it as such")
ok(errs(param_gwas(d0$stage1, transform = "custom_W", transform_args = list(W = W),
                   se_mode = "tau_cov", cov_data = d0$cov_data,
                   tau_cov_estimator = "gls", unknown_nu = "error")),
   "unknown_nu = 'error' stops")

# chisq_naive keeps the same contract as on the OLS path: no calibrated triple.
gc <- fit_gls(d0, q_method = "chisq_naive", nu = 400)
ok(is.null(gc$Q_pval), "chisq_naive returns no calibrated p-value on the GLS path")
ok(isTRUE(all.equal(gc$Q, g0$Q, tolerance = 0)),
   "the Q STATISTIC is bit-identical across q_method on the GLS path")

cat("\nAll GLS Q p-value tests passed.\n")
