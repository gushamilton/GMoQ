# ============================================================
# Test: calibrated null distribution for the lack-of-fit statistic Q
#
#  1. the DEFAULT behaviour is unchanged (chisq_naive reproduces the old
#     p-value exactly) -- the published results must stay reproducible;
#  2. the C++ helper agrees with a plain-R reference implementation;
#  3. the covariance upper triangle is read in ROW-major order;
#  4. the calibrated p-values are ordered sensibly and lie in [0, 1].
# ============================================================

suppressPackageStartupMessages(library(Rcpp))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg) else "tests/test-q-calibration.R"
pkg_root <- normalizePath(file.path(dirname(this_file), ".."))

Rcpp::sourceCpp(file.path(pkg_root, "src", "stage2_se_helper.cpp"))
source(file.path(pkg_root, "R", "utils.R"))
source(file.path(pkg_root, "R", "covariance.R"))
source(file.path(pkg_root, "R", "weight-functions.R"))
source(file.path(pkg_root, "R", "gwas.R"))

ok <- function(cond, msg) {
  if (isTRUE(cond)) cat("PASS:", msg, "\n")
  else stop("FAIL: ", msg)
}

set.seed(20260911)
Tt <- 12L
P  <- 40L

upper_row_major <- function(S) {
  out <- numeric(nrow(S) * (nrow(S) + 1) / 2); idx <- 1L
  for (r in seq_len(nrow(S))) for (c in r:ncol(S)) { out[idx] <- S[r, c]; idx <- idx + 1L }
  out
}

# A correlated, positive definite Sigma per SNP.
mk_sigma <- function() {
  A <- matrix(rnorm(Tt * Tt), Tt, Tt)
  S <- crossprod(A) / Tt + diag(0.5, Tt)
  d <- sqrt(diag(S)) * runif(Tt, 0.5, 1.5)
  outer(d, d) * cov2cor(S)
}
Sigmas <- replicate(P, mk_sigma(), simplify = FALSE)

taus <- seq(0.1, 0.9, length.out = Tt)
W <- cbind(1, seq_len(Tt) / Tt)
K <- ncol(W)

SE_tau <- vapply(Sigmas, function(S) sqrt(diag(S)), numeric(Tt))   # Tt x P
Q_slope <- vapply(Sigmas, function(S) {
  as.numeric(W %*% c(0.3, -0.2) + t(chol(S)) %*% rnorm(Tt))
}, numeric(Tt))

cov_vec <- unlist(lapply(Sigmas, upper_row_major))
offsets <- seq(0, by = Tt * (Tt + 1) / 2 * 4, length.out = P)
cov_data <- list(cov_vec = cov_vec, offsets = offsets)

stage1 <- list(taus = taus, q_tau = rep(0, Tt), Q_slope = Q_slope, SE_tau = SE_tau)

# ---- 1. default is the calibrated reference; chisq_naive reproduces the old --
fit_new <- param_gwas(stage1, transform = "custom_W", transform_args = list(W = W),
                      se_mode = "tau_cov", cov_data = cov_data, tau_cov_estimator = "ols")
ok(identical(fit_new$q_method, "saddlepoint"), "default q_method is saddlepoint")
ok(!is.null(fit_new$Q_pval), "default q_method returns a calibrated p-value")

# chisq_naive remains reachable and reproduces the pre-2026-09-11 published behaviour
fit_old <- param_gwas(stage1, transform = "custom_W", transform_args = list(W = W),
                      se_mode = "tau_cov", cov_data = cov_data,
                      tau_cov_estimator = "ols", q_method = "chisq_naive")
ok(identical(fit_old$q_method, "chisq_naive"), "chisq_naive is still selectable")
ok(is.null(fit_old$Q_pval), "chisq_naive returns no calibrated p-value, as published")

ok(isTRUE(all.equal(fit_old$Q, fit_new$Q, tolerance = 0)),
   "the Q STATISTIC is bit-identical across q_method (only the reference changes)")
ok(isTRUE(all.equal(fit_old$params, fit_new$params, tolerance = 0)),
   "parameter estimates are unchanged")

# ---- 2. C++ agrees with a plain-R reference ----------------------------
M <- diag(Tt) - W %*% solve(crossprod(W), t(W))
cal <- compute_q_calibration(cov_vec, offsets, fit_old$Q, SE_tau^2, M, Tt, 1L)
ref <- t(vapply(seq_len(P), function(j) {
  d <- diag(Sigmas[[j]])
  B <- diag(1 / sqrt(d)) %*% M %*% Sigmas[[j]] %*% t(diag(1 / sqrt(d)) %*% M)
  ev <- eigen((B + t(B)) / 2, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > 1e-10 * max(ev)]
  c(sum(ev), sum(ev^2), sum(ev^3), sum(ev^4))
}, numeric(4)))
for (i in 1:4) {
  ok(max(abs(cal[, i] - ref[, i]) / pmax(abs(ref[, i]), 1e-12)) < 1e-8,
     sprintf("C++ trace c%d matches the R reference", i))
}
ok(all(cal[, 6] == Tt - K), "number of non-zero eigenvalues equals T - k")

# ---- 3. row-major storage order ----------------------------------------
# se_tau^2 is diag(Sigma); if the triangle were read column-major the
# implied diagonal would differ for every entry after the first.
se2 <- compute_calibrated_se(cov_vec, offsets, diag(Tt), Tt)
ok(max(abs(t(se2)^2 - SE_tau^2) / pmax(SE_tau^2, 1e-12)) < 1e-10,
   "covariance upper triangle is read in row-major order")

# ---- 4. calibrated p-values behave ------------------------------------
for (meth in c("satterthwaite", "liu", "saddlepoint")) {
  f <- param_gwas(stage1, transform = "custom_W", transform_args = list(W = W),
                  se_mode = "tau_cov", cov_data = cov_data,
                  tau_cov_estimator = "ols", q_method = meth)
  ok(all(is.finite(f$Q_pval)) && all(f$Q_pval >= 0 & f$Q_pval <= 1),
     sprintf("%s returns finite p-values in [0, 1]", meth))
  ok(cor(f$Q, f$Q_pval, method = "spearman") < 0,
     sprintf("%s p-value decreases with Q", meth))
}

f <- param_gwas(stage1, transform = "custom_W", transform_args = list(W = W),
                se_mode = "tau_cov", cov_data = cov_data,
                tau_cov_estimator = "ols", q_method = "saddlepoint")
ok(all(f$Q_df_eff > 0) && all(f$Q_df_eff <= Tt - K + 1e-8),
   "effective df is positive and at most T - k")
ok(median(f$Q_df_eff) < (Tt - K),
   "correlated residuals give an effective df below the naive T - k")

cat("\nAll q-calibration tests passed.\n")
