# ============================================================
# Test: parameter p-values use a t_nu reference, with nu carried
#       from Stage 1 metadata
#
# Sigma is ESTIMATED by a delete-one-block jackknife, so for a fixed
# contrast a, a' Sigma_jk a / a' Sigma a ~ chi2_nu / nu and therefore
# theta/SE ~ t_nu exactly. Reading it against a normal inflates Type-I
# error 17x at p = 5e-8 when B = 65.
#
# Covered here:
#   1. nu round-trips Stage 1 .cov.meta.json -> Stage 2 output
#   2. <param>_p equals 2 * pt(-|theta/se|, nu) to numerical precision
#   3. a large nu converges to the normal answer
#   4. a covariance file with no block count emits NO p-value (default),
#      or a normal one only when that is asked for explicitly
#   5. the reference distribution actually used is recorded in the output
# ============================================================

suppressPackageStartupMessages(library(fungwasStage2))
suppressPackageStartupMessages(library(data.table))

cat("fungwasStage2 loaded from:", find.package("fungwasStage2"), "\n")

set.seed(20260914)

tmp <- file.path(tempdir(), paste0("tref_", as.integer(runif(1, 1, 1e6))))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

Tt <- 9L
P <- 40L
taus <- seq(0.1, 0.9, length.out = Tt)

# --- a positive-definite tau covariance with a realistic bridge-like shape ---
R_bb <- outer(taus, taus, function(s, t) {
  (pmin(s, t) - s * t) / sqrt(s * (1 - s) * t * (1 - t))
})
sd_tau <- runif(Tt, 0.02, 0.05)
Sigma <- diag(sd_tau) %*% R_bb %*% diag(sd_tau)
Sigma <- 0.5 * (Sigma + t(Sigma))
stopifnot(min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) > 0)

upper_idx <- which(upper.tri(Sigma, diag = TRUE), arr.ind = TRUE)
upper_idx <- upper_idx[order(upper_idx[, "row"], upper_idx[, "col"]), ]
cov_upper <- Sigma[upper_idx]
stopifnot(length(cov_upper) == Tt * (Tt + 1) / 2)

beta_mat <- matrix(rnorm(Tt * P, sd = 0.03), nrow = Tt, ncol = P)
se_mat <- matrix(rep(sqrt(diag(Sigma)), P), nrow = Tt, ncol = P)

write_stage1 <- function(prefix) {
  dt <- data.table(
    snp_id = paste0("rs", seq_len(P)),
    chr = 1L,
    bp = seq_len(P) * 1000L,
    effect_allele = "A",
    other_allele = "G"
  )
  for (i in seq_len(Tt)) {
    dt[[sprintf("beta_tau_%.3f", taus[i])]] <- beta_mat[i, ]
    dt[[sprintf("se_tau_%.3f", taus[i])]] <- se_mat[i, ]
  }
  fwrite(dt, paste0(prefix, ".stage1.tsv.gz"), sep = "\t")
  fwrite(dt[, .(snp_id, chr, bp, effect_allele, other_allele)],
         paste0(prefix, ".cov.ids.tsv.gz"), sep = "\t")

  con <- gzfile(paste0(prefix, ".cov.gz"), "wb")
  for (j in seq_len(P)) writeBin(as.numeric(cov_upper), con, size = 4)
  close(con)
  invisible(prefix)
}

write_meta_v2 <- function(prefix, n_blocks, nu = NULL, nu_source = NULL) {
  nu_txt <- if (is.null(nu)) n_blocks - 1 else nu
  src_txt <- if (is.null(nu_source)) "blocks_minus_one" else nu_source
  writeLines(c(
    "{",
    '  "cov_meta_version": 2,',
    '  "cov_dtype": "float32",',
    sprintf('  "n_taus": %d,', Tt),
    sprintf('  "n_cov_elements": %d,', Tt * (Tt + 1) / 2),
    '  "scale_file": false,',
    '  "jackknife": {',
    '    "method": "delete_one_block",',
    sprintf('    "n_blocks": %d,', n_blocks),
    sprintf('    "n_blocks_nonempty": %d,', n_blocks),
    '    "seed": 42,',
    sprintf('    "nu": %s,', format(nu_txt, scientific = FALSE)),
    sprintf('    "nu_source": "%s"', src_txt),
    "  }",
    "}"
  ), paste0(prefix, ".cov.meta.json"))
}

write_meta_v1 <- function(prefix) {
  # Exactly what Stage 1 used to write: no jackknife block at all.
  writeLines(c(
    "{",
    '  "cov_dtype": "float32",',
    sprintf('  "n_taus": %d,', Tt),
    sprintf('  "n_cov_elements": %d,', Tt * (Tt + 1) / 2),
    sprintf('  "record_bytes": %d,', Tt * (Tt + 1) / 2 * 4),
    '  "scale_file": false',
    "}"
  ), paste0(prefix, ".cov.meta.json"))
}

W <- cbind(intercept = rep(1, Tt), slope = taus - mean(taus))
A <- solve(crossprod(W), t(W))
theta_exp <- A %*% beta_mat
se_exp <- sqrt(diag(A %*% Sigma %*% t(A)))

run_stage2 <- function(prefix, ...) {
  param_gwas_from_file(
    stage1_file = paste0(prefix, ".stage1.tsv.gz"),
    W = W,
    se_mode = "tau_cov",
    cov_file = paste0(prefix, ".cov.gz"),
    tau_cov_estimator = "ols",
    q_method = "chisq_naive",
    ...
  )
}

# ------------------------------------------------------------------
# 1. nu round-trips from Stage 1 metadata, and the p-value is 2*pt()
# ------------------------------------------------------------------
n_blocks <- 65L
pfx <- file.path(tmp, "b65")
write_stage1(pfx)
write_meta_v2(pfx, n_blocks)

res <- suppressMessages(run_stage2(pfx))

stopifnot(all(c("intercept", "intercept_se", "intercept_p",
                "slope", "slope_se", "slope_p",
                "ref_dist", "ref_df") %in% names(res)))
stopifnot(identical(unique(res$ref_dist), "t"))
stopifnot(identical(unique(res$ref_df), as.numeric(n_blocks - 1L)))
cat("OK: nu round-trips from .cov.meta.json (n_blocks =", n_blocks,
    "-> ref_df =", unique(res$ref_df), ")\n")

# the SEs must be the fixed-contrast ones this reference is exact for
stopifnot(max(abs(res$intercept_se - se_exp[1])) < 1e-6)
stopifnot(max(abs(res$slope_se - se_exp[2])) < 1e-6)
stopifnot(max(abs(res$intercept - theta_exp[1, ])) < 1e-8)

nu <- n_blocks - 1L
p_expect_int <- 2 * pt(-abs(res$intercept / res$intercept_se), df = nu)
p_expect_slp <- 2 * pt(-abs(res$slope / res$slope_se), df = nu)
err <- max(abs(res$intercept_p - p_expect_int), abs(res$slope_p - p_expect_slp))
stopifnot(err < 1e-14)
cat("OK: p equals 2*pt(-|theta/se|, nu) to", format(err, digits = 3), "\n")

# and it is NOT the normal answer, in the direction that matters
p_norm <- 2 * pnorm(-abs(res$intercept / res$intercept_se))
stopifnot(all(res$intercept_p >= p_norm - 1e-15))
sharpest <- which.min(p_norm)
stopifnot(res$intercept_p[sharpest] > p_norm[sharpest])
cat("OK: t p-values are uniformly >= normal ones; at the most extreme SNP",
    format(p_norm[sharpest], digits = 3), "->",
    format(res$intercept_p[sharpest], digits = 3), "\n")

# ------------------------------------------------------------------
# 2. Q and Q_pval are untouched by this change
# ------------------------------------------------------------------
stopifnot("Q" %in% names(res), "Q_pval" %in% names(res))
stopifnot(all(is.finite(res$Q)))
cat("OK: Q / Q_pval still present and finite\n")

# ------------------------------------------------------------------
# 3. large nu converges to the normal reference
# ------------------------------------------------------------------
pfx_big <- file.path(tmp, "bbig")
write_stage1(pfx_big)
write_meta_v2(pfx_big, n_blocks = 1000001L)
res_big <- suppressMessages(run_stage2(pfx_big))
stopifnot(identical(unique(res_big$ref_dist), "t"))
conv <- max(abs(res_big$intercept_p - 2 * pnorm(-abs(res_big$intercept / res_big$intercept_se))))
stopifnot(conv < 1e-6)
cat("OK: at nu = 1e6 the t p-value matches the normal one to",
    format(conv, digits = 3), "\n")

# ------------------------------------------------------------------
# 4. old-format covariance: refuse to emit a p-value by default
# ------------------------------------------------------------------
pfx_old <- file.path(tmp, "old")
write_stage1(pfx_old)
write_meta_v1(pfx_old)

warned <- FALSE
res_old <- withCallingHandlers(
  suppressMessages(run_stage2(pfx_old)),
  warning = function(w) {
    if (grepl("Degrees of freedom", conditionMessage(w))) warned <<- TRUE
    invokeRestart("muffleWarning")
  }
)
stopifnot(warned)
stopifnot(!("intercept_p" %in% names(res_old)))
stopifnot(!("slope_p" %in% names(res_old)))
stopifnot(identical(unique(res_old$ref_dist), "none"))
stopifnot(all(is.na(res_old$ref_df)))
# everything else still works
stopifnot(all(is.finite(res_old$intercept)), all(is.finite(res_old$intercept_se)))
cat("OK: pre-v2 covariance yields no p-value, warns, and records ref_dist = none\n")

# a missing sidecar behaves the same way as a pre-v2 one
pfx_none <- file.path(tmp, "nometa")
write_stage1(pfx_none)
res_none <- suppressWarnings(suppressMessages(run_stage2(pfx_none)))
stopifnot(!("intercept_p" %in% names(res_none)))
stopifnot(identical(unique(res_none$ref_dist), "none"))
cat("OK: missing .cov.meta.json behaves identically to a pre-v2 one\n")

# ------------------------------------------------------------------
# 5. the normal fallback exists but must be asked for, and is recorded
# ------------------------------------------------------------------
res_norm <- suppressWarnings(suppressMessages(
  run_stage2(pfx_old, unknown_nu = "normal")
))
stopifnot("intercept_p" %in% names(res_norm))
stopifnot(identical(unique(res_norm$ref_dist), "normal"))
stopifnot(all(!is.finite(res_norm$ref_df)))
err_n <- max(abs(res_norm$intercept_p -
                   2 * pnorm(-abs(res_norm$intercept / res_norm$intercept_se))))
stopifnot(err_n < 1e-14)
cat("OK: unknown_nu = 'normal' emits normal p-values and records ref_dist = normal\n")

# ------------------------------------------------------------------
# 6. unknown_nu = "error" refuses outright
# ------------------------------------------------------------------
failed <- inherits(
  try(suppressMessages(run_stage2(pfx_old, unknown_nu = "error")), silent = TRUE),
  "try-error"
)
stopifnot(failed)
cat("OK: unknown_nu = 'error' stops rather than guessing\n")

# ------------------------------------------------------------------
# 7. an explicit nu overrides the sidecar, for old files
# ------------------------------------------------------------------
res_ovr <- suppressMessages(run_stage2(pfx_old, nu = 199))
stopifnot(identical(unique(res_ovr$ref_dist), "t"))
stopifnot(identical(unique(res_ovr$ref_df), 199))
err_o <- max(abs(res_ovr$intercept_p -
                   2 * pt(-abs(res_ovr$intercept / res_ovr$intercept_se), df = 199)))
stopifnot(err_o < 1e-14)
cat("OK: explicit nu = 199 overrides a missing block count\n")

# an explicit nu also wins over a sidecar that has one
res_ovr2 <- suppressMessages(run_stage2(pfx, nu = 499))
stopifnot(identical(unique(res_ovr2$ref_df), 499))
cat("OK: explicit nu wins over the recorded block count\n")

# ------------------------------------------------------------------
# 8. a measured nu_eff can replace B - 1 with no format change
# ------------------------------------------------------------------
pfx_eff <- file.path(tmp, "eff")
write_stage1(pfx_eff)
write_meta_v2(pfx_eff, n_blocks = 1000L, nu = 812.4, nu_source = "measured_effective")
res_eff <- suppressMessages(run_stage2(pfx_eff))
stopifnot(identical(unique(res_eff$ref_df), 812.4))
err_e <- max(abs(res_eff$intercept_p -
                   2 * pt(-abs(res_eff$intercept / res_eff$intercept_se), df = 812.4)))
stopifnot(err_e < 1e-14)
cat("OK: a measured nu_eff in the same field is honoured (nu = 812.4)\n")

# ------------------------------------------------------------------
# 9. internal resolvers, directly
# ------------------------------------------------------------------
resolve_nu <- fungwasStage2:::.resolve_nu
param_reference <- fungwasStage2:::.param_reference
read_cov_meta <- fungwasStage2:::.read_cov_meta
cov_meta_file <- fungwasStage2:::.cov_meta_file

stopifnot(identical(cov_meta_file("x/y.cov.gz"), "x/y.cov.meta.json"))
stopifnot(is.null(read_cov_meta(file.path(tmp, "does-not-exist.cov.meta.json"))))

m65 <- read_cov_meta(paste0(pfx, ".cov.meta.json"))
stopifnot(isTRUE(all.equal(resolve_nu(m65)$nu, 64)))
stopifnot(identical(resolve_nu(m65)$nu_source, "blocks_minus_one"))

m_old <- read_cov_meta(paste0(pfx_old, ".cov.meta.json"))
stopifnot(is.na(resolve_nu(m_old)$nu))
stopifnot(is.na(resolve_nu(NULL)$nu))

stopifnot(identical(param_reference(nu = 64)$dist, "t"))
stopifnot(identical(param_reference(nu = 64)$df, 64))
stopifnot(identical(
  suppressWarnings(param_reference(nu = NA_real_, unknown_nu = "omit"))$dist, "none"))
stopifnot(identical(
  suppressWarnings(param_reference(nu = NA_real_, unknown_nu = "normal"))$dist, "normal"))
stopifnot(inherits(
  try(param_reference(nu = NA_real_, unknown_nu = "error"), silent = TRUE), "try-error"))
# a non-positive or non-finite nu is treated as unknown, never as a reference
stopifnot(identical(suppressWarnings(param_reference(nu = 0))$dist, "none"))
stopifnot(identical(suppressWarnings(param_reference(nu = -5))$dist, "none"))
cat("OK: .resolve_nu / .param_reference behave as specified\n")

# the base-R JSON fallback must agree with jsonlite where jsonlite exists
if (requireNamespace("jsonlite", quietly = TRUE)) {
  txt <- paste(readLines(paste0(pfx, ".cov.meta.json"), warn = FALSE), collapse = " ")
  stopifnot(identical(fungwasStage2:::.json_number(txt, "n_blocks"), 65))
  stopifnot(identical(fungwasStage2:::.json_number(txt, "nu"), 64))
  stopifnot(identical(fungwasStage2:::.json_string(txt, "nu_source"), "blocks_minus_one"))
  # "nu" must not be matched by the "nu_source" key
  stopifnot(fungwasStage2:::.json_number(txt, "nu") != 0)
  cat("OK: base-R JSON fallback agrees with the jsonlite parse\n")
}

# ------------------------------------------------------------------
# 10. no joint multi-parameter Wald test exists in the package, so there
#     is nothing here for an F reference to apply to. Q is a lack-of-fit
#     statistic on T - k residual df, not a joint test of k parameters,
#     and it is out of scope. If a joint test is added, it must use
#     F_{k, nu-k+1} (Hotelling), not chi-square.
# ------------------------------------------------------------------
stopifnot(!any(grepl("joint", getNamespaceExports("fungwasStage2"), ignore.case = TRUE)))
cat("NOTE: no joint k-parameter test exists; F reference not applicable\n")

cat("OK: all t-reference tests passed\n")
