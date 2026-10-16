# ============================================================
# Regression: pooled Sigma_0 uses a shared per-chromosome payload
# ============================================================
# This fixture is deliberately built from independently supplied SPD matrices,
# betas and d_resid values. It does not use the Stage 1 reader or regenerate
# truth from the functions under test. The compact and compatibility-dense
# representations must give the same OLS, GLS, Q-calibration, covariance and
# model-tournament results, including interleaved chromosomes and unusable d.

suppressPackageStartupMessages(library(Rcpp))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg) else
  "tests/test-pooled-memory-repair.R"
pkg_root <- normalizePath(file.path(dirname(this_file), ".."))

source(file.path(pkg_root, "tests", "load-installed.R"))

assert_close <- function(x, y, tol = 2e-10, label = "value") {
  if (is.null(x) || is.null(y)) stop(label, " is unexpectedly NULL")
  if (length(x) != length(y) || !identical(dim(x), dim(y))) {
    stop(label, " dimensions differ")
  }
  ok <- is.na(x) == is.na(y)
  if (any(!ok)) stop(label, " NA pattern differs")
  if (any(ok & !is.na(x))) {
    delta <- max(abs(x[ok & !is.na(x)] - y[ok & !is.na(y)]))
    if (!is.finite(delta) || delta > tol) stop(label, " differs by ", delta)
  }
  invisible(TRUE)
}

set.seed(20260915)
Tt <- 5L
P <- 8L
taus <- c(0.10, 0.25, 0.50, 0.75, 0.90)
chr <- c(21L, 22L, 21L, 22L, 21L, 22L, 21L, 22L)
d <- c(2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0)

# External known-truth covariance matrices: Sigma_v = Sigma_0 / d.
make_spd <- function(seed) {
  set.seed(seed)
  X <- matrix(rnorm(Tt * Tt), Tt, Tt)
  S <- crossprod(X) + diag(seq(0.7, 1.1, length.out = Tt))
  0.5 * (S + t(S))
}
S21 <- make_spd(21L)
S22 <- make_spd(22L)
mats <- list("21" = S21, "22" = S22)

dense <- .pooled_cov_data(mats, chr, d, Tt)
compact <- .pooled_compact_cov_data(mats, chr, d, Tt)
ncov <- Tt * (Tt + 1L) / 2L

# The compact payload is O(chromosomes*T^2), with one scalar per SNP. The dense
# compatibility object is intentionally O(SNP*T^2).
stopifnot(length(dense$cov_vec) == P * ncov)
stopifnot(length(compact$cov_vec) == length(mats) * ncov)
stopifnot(length(compact$cov_vec) < length(dense$cov_vec))
stopifnot(length(compact$offsets) == P, length(compact$cov_scales) == P)
stopifnot(all(compact$offsets[c(TRUE, FALSE)] == 0))
stopifnot(all(compact$offsets[c(FALSE, TRUE)] == ncov * 4))
assert_close(compact$cov_scales, 1 / d, label = "cov_scales")

# Exercise the production sidecar path as well as the pure builder. The
# metadata is a small hand-written fixture with matching provenance, so this
# still does not depend on a Stage 1 reader or generated truth.
tmp <- tempfile("pooled_memory_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
pfx <- file.path(tmp, "fixture")
sigma_file <- paste0(pfx, ".sigma0.tsv.gz")
con <- gzfile(sigma_file, "wt")
writeLines(paste(c("chr", "tau", paste0("s_", sprintf("%.3f", taus))),
                 collapse = "\t"), con)
for (k in names(mats)) for (i in seq_len(Tt)) {
  writeLines(paste(c(k, sprintf("%.3f", taus[i]),
                     sprintf("%.17g", mats[[k]][i, ])), collapse = "\t"), con)
}
close(con)
writeLines(c(
  "{", '  "sigma0_meta_version": 1,',
  '  "package_version": "fixture",',
  '  "estimator": "fixture",',
  '  "phenotype": "fixture",', sprintf('  "n_taus": %d,', Tt),
  sprintf('  "taus": [%s],', paste(sprintf("%.3f", taus), collapse = ", ")),
  '  "n_samples": 100,', '  "n_covariates": 1,',
  '  "covariates": ["age"],', '  "loco_predictor": null,',
  '  "loco_scale_mode": null,', '  "loco_scale_factor": null,',
  '  "chromosomes": ["21", "22"],',
  '  "fingerprint": "fixture-fingerprint"', "}"),
  paste0(pfx, ".sigma0.meta.json"))
writeLines(c(
  "{", '  "cov_meta_version": 4,',
  '  "jackknife": {"n_blocks": 10, "n_blocks_nonempty": 10},',
  '  "analysis": {', '    "phenotype": "fixture",',
  '    "n_samples": 100,', '    "n_covariates": 1,',
  '    "covariates": ["age"],', '    "loco_predictor": null,',
  '    "loco_scale_mode": null,', '    "loco_scale_factor": null,',
  '    "fingerprint": "fixture-fingerprint"', '  }', "}"),
  paste0(pfx, ".cov.meta.json"))
from_file <- .pooled_tau_cov(sigma_file, chr = chr, d_resid = d, taus = taus)
stopifnot(identical(from_file$cov_data$representation, "chromosome_shared"))
assert_close(from_file$cov_data$cov_vec, compact$cov_vec,
             label = "production compact covariance")
assert_close(from_file$cov_data$cov_scales, compact$cov_scales,
             label = "production compact scales")

# Expand compact data independently and compare it with the old dense layout.
expanded <- unlist(lapply(seq_len(P), function(i) {
  j <- compact$offsets[i] / 4 + seq_len(ncov)
  compact$cov_vec[j] * compact$cov_scales[i]
}), use.names = FALSE)
assert_close(expanded, dense$cov_vec, label = "expanded compact covariance")

W <- cbind(location = 1, spread = qnorm(taus) - mean(qnorm(taus)))
A <- solve(crossprod(W), t(W))
beta <- matrix(c(
  0.12, 0.03, -0.01, 0.04, 0.02, 0.06, 0.01, -0.02,
  0.08, 0.00,  0.02, 0.01, 0.04, 0.03, 0.00,  0.01,
  0.05, 0.02,  0.01, 0.00, 0.03, 0.02, 0.01, -0.01,
  0.02, 0.01,  0.00, 0.02, 0.01, 0.00, 0.01,  0.02,
  0.01, 0.00,  0.01, 0.01, 0.00, 0.01, 0.00,  0.01
), nrow = Tt, byrow = TRUE)
SE_tau <- vapply(seq_len(P), function(i) {
  S <- if (chr[i] == 21L) S21 else S22
  sqrt(diag(S) / d[i])
}, numeric(Tt))
stage1 <- list(taus = taus, q_tau = rep(0, Tt),
               Q_slope = beta, SE_tau = SE_tau)

# Every downstream consumer receives the same optional scale vector.
se_dense <- compute_calibrated_se(dense$cov_vec, dense$offsets, A, Tt)
se_compact <- compute_calibrated_se(compact$cov_vec, compact$offsets, A, Tt,
                                    cov_scales = compact$cov_scales)
assert_close(se_dense, se_compact, label = "calibrated SE")

pc_dense <- compute_param_cov(dense$cov_vec, dense$offsets, A, Tt)
pc_compact <- compute_param_cov(compact$cov_vec, compact$offsets, A, Tt,
                                cov_scales = compact$cov_scales)
assert_close(pc_dense, pc_compact, label = "parameter covariance")

gls_dense <- compute_gls_param_gwas(beta, dense$cov_vec, dense$offsets, W,
                                    return_cov = TRUE, n_threads = 1L)
gls_compact <- compute_gls_param_gwas(
  beta, compact$cov_vec, compact$offsets, W,
  return_cov = TRUE, n_threads = 1L, cov_scales = compact$cov_scales
)
assert_close(gls_dense$params, gls_compact$params, label = "GLS parameters")
assert_close(gls_dense$se, gls_compact$se, label = "GLS SE")
assert_close(gls_dense$Q, gls_compact$Q, label = "GLS Q")
assert_close(gls_dense$param_cov_upper, gls_compact$param_cov_upper,
             label = "GLS parameter covariance")

fits <- lapply(list(dense = dense, compact = compact), function(cd) {
  param_gwas(stage1, transform = "custom_W", transform_args = list(W = W),
             se_mode = "tau_cov", cov_data = cd, tau_cov_estimator = "ols",
             return_cov = "upper", q_method = "saddlepoint", nu = 100)
})
assert_close(fits[[1]]$params, fits[[2]]$params, label = "OLS parameters")
assert_close(fits[[1]]$SE_params, fits[[2]]$SE_params, label = "OLS SE")
assert_close(fits[[1]]$Q, fits[[2]]$Q, label = "OLS Q")
assert_close(fits[[1]]$Q_pval, fits[[2]]$Q_pval, label = "OLS Q p-value")
assert_close(fits[[1]]$P_params, fits[[2]]$P_params, label = "OLS parameter p-values")
assert_close(fits[[1]]$param_cov_upper, fits[[2]]$param_cov_upper,
             label = "OLS parameter covariance")

gls_fits <- lapply(list(dense = dense, compact = compact), function(cd) {
  param_gwas(stage1, transform = "custom_W", transform_args = list(W = W),
             se_mode = "tau_cov", cov_data = cd, tau_cov_estimator = "gls",
             return_cov = "upper", q_method = "saddlepoint", nu = 100)
})
assert_close(gls_fits[[1]]$params, gls_fits[[2]]$params, label = "mapped GLS parameters")
assert_close(gls_fits[[1]]$SE_params, gls_fits[[2]]$SE_params, label = "mapped GLS SE")
assert_close(gls_fits[[1]]$Q, gls_fits[[2]]$Q, label = "mapped GLS Q")
assert_close(gls_fits[[1]]$Q_pval, gls_fits[[2]]$Q_pval, label = "mapped GLS Q p-value")
assert_close(gls_fits[[1]]$param_cov_upper, gls_fits[[2]]$param_cov_upper,
             label = "mapped GLS parameter covariance")

W_list <- list(location = W[, 1, drop = FALSE], both = W)
multi_dense <- param_gwas_multi(stage1, W_list, N = 100, cov_data = dense)
multi_compact <- param_gwas_multi(stage1, W_list, N = 100, cov_data = compact)
assert_close(multi_dense, multi_compact, label = "model tournament")

# A bad d_resid is represented by an NA scale and is skipped by all kernels;
# it must not manufacture a dense per-SNP zero matrix.
d_bad <- d
d_bad[3] <- 0
compact_bad <- .pooled_compact_cov_data(mats, chr, d_bad, Tt)
stopifnot(compact_bad$n_bad == 1L, is.na(compact_bad$cov_scales[3]))
se_bad <- compute_calibrated_se(compact_bad$cov_vec, compact_bad$offsets, A, Tt,
                                cov_scales = compact_bad$cov_scales)
stopifnot(all(is.na(se_bad[3, ])), all(is.finite(se_bad[-3, ])))
gls_bad <- compute_gls_param_gwas(
  beta, compact_bad$cov_vec, compact_bad$offsets, W,
  return_cov = TRUE, n_threads = 1L, cov_scales = compact_bad$cov_scales
)
stopifnot(all(is.na(gls_bad$params[, 3])), all(is.finite(gls_bad$params[, -3])))

cat("OK: compact pooled covariance matches dense consumers; payload is chromosome-shared\n")
