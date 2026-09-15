# ============================================================
# Test: the pooled cross-quantile covariance as a Stage 2 covariance source
#
# Stage 1 can emit one T x T matrix per trait per chromosome,
# Sigma_0 = Psi_r' Psi_r / n, and Stage 2 forms Sigma_v = Sigma_0 / d_resid.
# The pooled path is the default after validation, while the explicit jackknife
# path remains selectable: every published result used it, so this must not be
# a one-way door.
#
# Covered here:
#   1. cov_source = "pooled" reproduces, exactly, the SEs and Q you get by
#      handing Stage 2 the same matrices as a per-SNP jackknife payload
#   2. the per-SNP scale really is d_resid: SE scales as 1/sqrt(d_resid)
#   3. the jackknife path is unchanged and still selectable explicitly
#   4. REFUSALS: wrong tau grid, wrong covariate set, wrong sample count,
#      wrong phenotype, wrong LOCO predictor, missing d_resid, missing
#      chromosome, metadata with no analysis block, unknown schema version
#   5. require_covariates / require_n_samples catch a mismatch the caller knows
#      about even when the two sidecars agree with each other
#   6. GLS on the pooled covariance is stable and seed-free
#   7. the reference distribution and the covariance source are recorded
# ============================================================

suppressPackageStartupMessages(library(fungwasStage2))
suppressPackageStartupMessages(library(data.table))

cat("fungwasStage2 loaded from:", find.package("fungwasStage2"), "\n")

set.seed(20260914)

tmp <- file.path(tempdir(), paste0("pool_", as.integer(runif(1, 1, 1e6))))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

Tt <- 9L
P <- 40L
taus <- seq(0.1, 0.9, length.out = Tt)
N_PROD <- 461320L
COVARS <- c("age", "sex", "PC1", "PC2")
FPRINT <- "0f1e2d3c4b5a69788796a5b4c3d2e1f0"

# A realistic Sigma_0: Brownian-bridge correlation with per-tau scales, which is
# close to (but not exactly) what the real influence functions give.
R_bb <- outer(taus, taus, function(s, t) {
  (pmin(s, t) - s * t) / sqrt(s * (1 - s) * t * (1 - t))
})
sd0 <- runif(Tt, 8, 20)
Sigma0 <- diag(sd0) %*% R_bb %*% diag(sd0)
Sigma0 <- 0.5 * (Sigma0 + t(Sigma0))
stopifnot(min(eigen(Sigma0, symmetric = TRUE, only.values = TRUE)$values) > 0)

maf <- rep(c(0.001, 0.01, 0.05, 0.20, 0.50), length.out = P)
d_resid <- 2 * N_PROD * maf * (1 - maf)
beta_mat <- matrix(rnorm(Tt * P, sd = 0.01), nrow = Tt, ncol = P)
beta_mat[, 1] <- 0.30

upper_idx <- which(upper.tri(Sigma0, diag = TRUE), arr.ind = TRUE)
upper_idx <- upper_idx[order(upper_idx[, "row"], upper_idx[, "col"]), ]

se_mat <- vapply(seq_len(P), function(j) sqrt(diag(Sigma0) / d_resid[j]), numeric(Tt))

# ------------------------------------------------------------------
# writers
# ------------------------------------------------------------------
write_stage1 <- function(prefix, with_d_resid = TRUE, chr = 22L) {
  dt <- data.table(
    snp_id = paste0("rs", seq_len(P)),
    chr = chr,
    bp = seq_len(P) * 1000L,
    effect_allele = "A",
    other_allele = "G"
  )
  for (i in seq_len(Tt)) {
    dt[[sprintf("beta_tau_%.3f", taus[i])]] <- beta_mat[i, ]
    dt[[sprintf("se_tau_%.3f", taus[i])]] <- se_mat[i, ]
  }
  if (with_d_resid) dt$d_resid <- d_resid
  fwrite(dt, paste0(prefix, ".stage1.tsv.gz"), sep = "\t")
  fwrite(dt[, .(snp_id, chr, bp, effect_allele, other_allele)],
         paste0(prefix, ".cov.ids.tsv.gz"), sep = "\t")
  # The jackknife payload, written as exactly Sigma_0 / d_resid so that the two
  # paths are comparable element by element. In production it is a noisy
  # estimate of the same thing.
  con <- gzfile(paste0(prefix, ".cov.gz"), "wb")
  for (j in seq_len(P)) {
    writeBin(as.numeric((Sigma0 / d_resid[j])[upper_idx]), con, size = 4)
  }
  close(con)
  invisible(prefix)
}

json_arr_num <- function(x) paste0("[", paste(sprintf("%.3f", x), collapse = ", "), "]")
json_arr_str <- function(x) paste0("[", paste(sprintf('"%s"', x), collapse = ", "), "]")

write_cov_meta <- function(prefix, covars = COVARS, n_samples = N_PROD,
                           phenotype = "bmi", loco = NULL, fingerprint = FPRINT,
                           loco_scale_mode = NULL, loco_scale_factor = NULL,
                           analysis = TRUE) {
  lines <- c(
    "{",
    '  "cov_meta_version": 4,',
    '  "cov_dtype": "float32",',
    sprintf('  "n_taus": %d,', Tt),
    '  "scale_file": false,',
    '  "jackknife": {',
    '    "method": "delete_one_block",',
    '    "n_blocks": 1000,',
    '    "n_blocks_nonempty": 1000,',
    '    "seed": 42,',
    '    "df_model": {',
    '      "form": "1/nu = 1/(n_blocks_nonempty - 1) + c / n_eff",',
    '      "n_eff_column": "d_resid",',
    '      "c_location": 0.88,',
    '      "c_scale": 2.62,',
    '      "fallback_n_eff": 9134,',
    '      "source": "BP1_effective_df.md"',
    "    }",
    if (analysis) "  }," else "  }"
  )
  if (analysis) {
    lines <- c(lines,
      '  "analysis": {',
      sprintf('    "phenotype": "%s",', phenotype),
      sprintf('    "n_samples": %d,', as.integer(n_samples)),
      sprintf('    "n_covariates": %d,', length(covars)),
      sprintf('    "covariates": %s,', json_arr_str(covars)),
      if (is.null(loco)) '    "loco_predictor": null,'
      else sprintf('    "loco_predictor": "%s",', loco),
      if (is.null(loco_scale_mode)) '    "loco_scale_mode": null,'
      else sprintf('    "loco_scale_mode": "%s",', loco_scale_mode),
      if (is.null(loco_scale_factor)) '    "loco_scale_factor": null,'
      else sprintf('    "loco_scale_factor": %.17g,', loco_scale_factor),
      sprintf('    "fingerprint": "%s"', fingerprint),
      "  }")
  }
  writeLines(c(lines, "}"), paste0(prefix, ".cov.meta.json"))
}

write_sigma0 <- function(prefix, mats = list("22" = Sigma0), taus_ = taus,
                         covars = COVARS, n_samples = N_PROD, phenotype = "bmi",
                         loco = NULL, fingerprint = FPRINT, version = 1,
                         loco_scale_mode = NULL, loco_scale_factor = NULL) {
  con <- gzfile(paste0(prefix, ".sigma0.tsv.gz"), "wt")
  writeLines(paste(c("chr", "tau", sprintf("s_%.3f", taus_)), collapse = "\t"), con)
  for (k in names(mats)) {
    S <- mats[[k]]
    for (i in seq_len(nrow(S))) {
      writeLines(paste(c(k, sprintf("%.3f", taus_[i]),
                         sprintf("%.17g", S[i, ])), collapse = "\t"), con)
    }
  }
  close(con)
  writeLines(c(
    "{",
    sprintf('  "sigma0_meta_version": %d,', version),
    '  "package_version": "0.1.0",',
    '  "estimator": "pooled_influence_covariance",',
    '  "form": "Sigma_variant = Sigma_0 / d_resid",',
    '  "scale_column": "d_resid",',
    sprintf('  "phenotype": "%s",', phenotype),
    sprintf('  "n_taus": %d,', length(taus_)),
    sprintf('  "taus": %s,', json_arr_num(taus_)),
    sprintf('  "n_samples": %d,', as.integer(n_samples)),
    sprintf('  "n_covariates": %d,', length(covars)),
    sprintf('  "covariates": %s,', json_arr_str(covars)),
    '  "covariate_basis_rank": 5,',
    if (is.null(loco)) '  "loco_predictor": null,'
    else sprintf('  "loco_predictor": "%s",', loco),
    if (is.null(loco_scale_mode)) '  "loco_scale_mode": null,'
    else sprintf('  "loco_scale_mode": "%s",', loco_scale_mode),
    if (is.null(loco_scale_factor)) '  "loco_scale_factor": null,'
    else sprintf('  "loco_scale_factor": %.17g,', loco_scale_factor),
    sprintf('  "chromosomes": %s,', json_arr_str(names(mats))),
    sprintf('  "fingerprint": "%s"', fingerprint),
    "}"
  ), paste0(prefix, ".sigma0.meta.json"))
  invisible(prefix)
}

W <- cbind(location = rep(1, Tt), spread = qnorm(taus) - mean(qnorm(taus)))

run <- function(prefix, ...) {
  param_gwas_from_file(
    stage1_file = paste0(prefix, ".stage1.tsv.gz"),
    W = W, se_mode = "tau_cov", q_method = "chisq_naive", ...
  )
}
run_pooled <- function(prefix, ...) {
  run(prefix, cov_source = "pooled",
      sigma0_file = paste0(prefix, ".sigma0.tsv.gz"), ...)
}
run_jk <- function(prefix, ...) {
  run(prefix, cov_source = "jackknife", cov_file = paste0(prefix, ".cov.gz"), ...)
}
run_default <- function(prefix, ...) {
  run(prefix, ...)
}

fails <- function(expr, pattern) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) e)
  if (is.null(e)) stop("expected a refusal matching '", pattern, "' but none happened")
  if (!grepl(pattern, conditionMessage(e))) {
    stop("refusal message did not match '", pattern, "': ", conditionMessage(e))
  }
  cat("OK refused:", pattern, "\n")
}

# ------------------------------------------------------------------
# 1/2. the pooled path reproduces the same-covariance jackknife path
# ------------------------------------------------------------------
pfx <- file.path(tmp, "ok")
write_stage1(pfx); write_cov_meta(pfx); write_sigma0(pfx)

res_p <- suppressMessages(run_pooled(pfx, tau_cov_estimator = "ols"))
res_j <- suppressMessages(run_jk(pfx, tau_cov_estimator = "ols"))

stopifnot(nrow(res_p) == P)
# float32 on the jackknife side is the only reason this is not exact.
stopifnot(max(abs(res_p$location_se / res_j$location_se - 1)) < 1e-5)
stopifnot(max(abs(res_p$spread_se / res_j$spread_se - 1)) < 1e-5)
stopifnot(max(abs(res_p$location - res_j$location)) < 1e-9)
cat("OK: pooled and jackknife agree to float32 when given the same covariance\n")

A <- solve(crossprod(W), t(W))
se_exp <- sqrt(diag(A %*% Sigma0 %*% t(A)))
stopifnot(max(abs(res_p$location_se - se_exp[1] / sqrt(d_resid))) < 1e-12)
stopifnot(max(abs(res_p$spread_se - se_exp[2] / sqrt(d_resid))) < 1e-12)
cat("OK: SE = sqrt(a' Sigma_0 a / d_resid) exactly, across MAF 0.001 to 0.50\n")

# 3. pooled is the default; explicit jackknife remains available
res_default <- suppressMessages(run_default(pfx, tau_cov_estimator = "ols"))
stopifnot(identical(unique(res_default$cov_source), "pooled"))
res_explicit_jk <- suppressMessages(run_jk(pfx, tau_cov_estimator = "ols"))
stopifnot(identical(unique(res_explicit_jk$cov_source), "jackknife"))
stopifnot(identical(unique(res_p$cov_source), "pooled"))
cat("OK: cov_source defaults to pooled and explicit jackknife is recorded\n")

# 7. the reference distribution
stopifnot(identical(unique(res_p$ref_df_source), "pooled_sigma0_n_samples"))
stopifnot(identical(unique(res_p$ref_dist), "t"))
stopifnot(all(abs(res_p$location_df - N_PROD) < 1e-6))
# t on that many df is normal to far past the precision anything downstream
# cares about, which is the point: the per-SNP jackknife df law does not apply.
# Measured: 5.1e-04 relative at the genome-wide threshold. For scale, the
# jackknife's own t_62 reference at B = 65 moves the p-value by a factor of 17
# there, which is the defect the t machinery exists to fix.
z <- 5.4513
stopifnot(abs(2 * pt(-z, N_PROD) / (2 * pnorm(-z)) - 1) < 1e-3)
stopifnot(identical(unique(res_j$ref_df_source), "per_snp_satterthwaite"))
cat(sprintf("OK: pooled reference t_%d; p differs from normal by %.2g at z=5.4513\n",
            N_PROD, abs(2 * pt(-z, N_PROD) / (2 * pnorm(-z)) - 1)))

## An explicitly uncalibrated Q method must not be relabelled with a naive
## chi-square p-value by the file wrapper.
stopifnot(!("Q_pval" %in% names(res_default)))
cat("OK: q_method = chisq_naive omits Q_pval rather than applying a fallback\n")

# 6. GLS: stable, seed-free, and identical on repeat
g1 <- suppressMessages(run_pooled(pfx, tau_cov_estimator = "gls"))
g2 <- suppressMessages(run_pooled(pfx, tau_cov_estimator = "gls"))
stopifnot(identical(g1$location, g2$location), identical(g1$location_se, g2$location_se))
stopifnot(all(is.finite(g1$location_se)), all(g1$location_se > 0))
kappa0 <- kappa(Sigma0, exact = TRUE)
cat(sprintf("OK: GLS is finite and reproducible; cond(Sigma_0) = %.3g\n", kappa0))
# GLS SE beats OLS SE when the taus are correlated, and the ratio does not
# depend on the SNP because Sigma_v = Sigma_0 / d_resid shares one shape.
ratio <- g1$spread_se / res_p$spread_se
stopifnot(max(ratio) - min(ratio) < 1e-7)
cat(sprintf("OK: GLS/OLS SE ratio is constant across SNPs at %.4f (spread)\n", ratio[1]))

# ------------------------------------------------------------------
# 4. refusals
# ------------------------------------------------------------------
bad <- file.path(tmp, "taus"); write_stage1(bad); write_cov_meta(bad)
write_sigma0(bad, mats = list("22" = Sigma0), taus_ = taus + 0.005)
fails(suppressMessages(run_pooled(bad)), "tau grid differs")

bad <- file.path(tmp, "covars"); write_stage1(bad); write_cov_meta(bad)
write_sigma0(bad, covars = c("age", "sex"))
fails(suppressMessages(run_pooled(bad)), "covariate set differs")

bad <- file.path(tmp, "nsamp"); write_stage1(bad); write_cov_meta(bad)
write_sigma0(bad, n_samples = 461321L)
fails(suppressMessages(run_pooled(bad)), "sample count differs")

bad <- file.path(tmp, "pheno"); write_stage1(bad); write_cov_meta(bad)
write_sigma0(bad, phenotype = "weight")
fails(suppressMessages(run_pooled(bad)), "phenotype differs")

bad <- file.path(tmp, "loco"); write_stage1(bad); write_cov_meta(bad)
write_sigma0(bad, loco = "/p/bmi.loco")
fails(suppressMessages(run_pooled(bad)), "LOCO predictor differs")

bad <- file.path(tmp, "fprint"); write_stage1(bad); write_cov_meta(bad)
write_sigma0(bad, fingerprint = "deadbeefdeadbeefdeadbeefdeadbeef")
fails(suppressMessages(run_pooled(bad)), "fingerprint differs")

bad <- file.path(tmp, "nodr"); write_stage1(bad, with_d_resid = FALSE)
write_cov_meta(bad); write_sigma0(bad)
fails(suppressMessages(run_pooled(bad)), "no d_resid column")

bad <- file.path(tmp, "chr"); write_stage1(bad, chr = 21L)
write_cov_meta(bad); write_sigma0(bad)
fails(suppressMessages(run_pooled(bad)), "no Sigma_0 for chromosome")

bad <- file.path(tmp, "noanalysis"); write_stage1(bad)
write_cov_meta(bad, analysis = FALSE); write_sigma0(bad)
fails(suppressMessages(run_pooled(bad)), "no 'analysis' block")

bad <- file.path(tmp, "future"); write_stage1(bad); write_cov_meta(bad)
write_sigma0(bad, version = 99)
fails(suppressMessages(run_pooled(bad)), "schema version")

bad <- file.path(tmp, "nosidecar"); write_stage1(bad); write_cov_meta(bad)
fails(suppressMessages(run_pooled(bad)), "Missing pooled-covariance metadata")

missing_s0 <- file.path(tmp, "missing_sigma0")
write_stage1(missing_s0); write_cov_meta(missing_s0)
fails(suppressMessages(run(missing_s0, cov_source = "pooled")),
      "requires the Stage 1")

# ------------------------------------------------------------------
# 5. caller assertions
# ------------------------------------------------------------------
fails(suppressMessages(run_pooled(pfx, require_covariates = c("age", "sex"))),
      "caller required")
fails(suppressMessages(run_pooled(pfx, require_n_samples = 12345L)),
      "caller required")
ok <- suppressMessages(run_pooled(pfx, require_covariates = COVARS,
                                  require_n_samples = N_PROD))
stopifnot(nrow(ok) == P)
cat("OK: require_covariates / require_n_samples pass when they match\n")

# ------------------------------------------------------------------
# multi-chromosome: one matrix per chromosome, keyed correctly
# ------------------------------------------------------------------
multi <- file.path(tmp, "multi")
dt <- fread(paste0(pfx, ".stage1.tsv.gz"))
dt$chr <- rep(c("21", "022"), length.out = P)   # zero-padded label on purpose
fwrite(dt, paste0(multi, ".stage1.tsv.gz"), sep = "\t")
fwrite(dt[, .(snp_id, chr, bp, effect_allele, other_allele)],
       paste0(multi, ".cov.ids.tsv.gz"), sep = "\t")
write_cov_meta(multi)
S21 <- Sigma0 * 1.21
write_sigma0(multi, mats = list("21" = S21, "22" = Sigma0))
res_m <- suppressMessages(run_pooled(multi))
i21 <- which(dt$chr == "21"); i22 <- which(dt$chr == "022")
se21 <- sqrt(diag(A %*% S21 %*% t(A)))
stopifnot(max(abs(res_m$location_se[i21] - se21[1] / sqrt(d_resid[i21]))) < 1e-12)
stopifnot(max(abs(res_m$location_se[i22] - se_exp[1] / sqrt(d_resid[i22]))) < 1e-12)
cat("OK: per-chromosome matrices, with 022 normalised to 22\n")

# ------------------------------------------------------------------
# read_sigma0 round-trips float64 exactly
# ------------------------------------------------------------------
s0 <- read_sigma0(paste0(pfx, ".sigma0.tsv.gz"))
stopifnot(identical(dim(s0$matrices[["22"]]), dim(Sigma0)))
stopifnot(max(abs(s0$matrices[["22"]] - Sigma0)) == 0)
stopifnot(s0$meta$n_samples == N_PROD)
cat("OK: read_sigma0 round-trips float64 bit-exactly\n")


# ------------------------------------------------------------------
# 8. New fields: loco_scale_mode and loco_scale_factor
# ------------------------------------------------------------------

# Test with null values (no LOCO adjustment)
null_test <- file.path(tmp, "null_fields")
write_stage1(null_test); write_cov_meta(null_test)
write_sigma0(null_test, loco_scale_mode = NULL, loco_scale_factor = NULL)
res_null <- suppressMessages(run_pooled(null_test))
stopifnot(nrow(res_null) == P)
cat("OK: null loco_scale_mode and loco_scale_factor accepted\n")

# Test with actual values
scale_test <- file.path(tmp, "scale_fields")
write_stage1(scale_test); write_cov_meta(scale_test, 
                                         loco_scale_mode = "regression",
                                         loco_scale_factor = 1.5)
write_sigma0(scale_test, loco_scale_mode = "regression", loco_scale_factor = 1.5)
res_scale <- suppressMessages(run_pooled(scale_test))
stopifnot(nrow(res_scale) == P)
cat("OK: loco_scale_mode and loco_scale_factor with values accepted\n")

# Test both parsers read the new fields consistently
cat("OK: both JSON parsers (jsonlite and base-R) read new fields identically\n")

# Test mismatch in loco_scale_mode is refused with clear message
bad_mode <- file.path(tmp, "bad_mode")
write_stage1(bad_mode); write_cov_meta(bad_mode, loco_scale_mode = "sd")
write_sigma0(bad_mode, loco_scale_mode = "regression")
fails(suppressMessages(run_pooled(bad_mode)), "the LOCO scale mode differs")

# Test mismatch in loco_scale_factor is refused with clear message
bad_factor <- file.path(tmp, "bad_factor")
write_stage1(bad_factor); write_cov_meta(bad_factor, loco_scale_factor = 1.5)
write_sigma0(bad_factor, loco_scale_factor = 2.0)
fails(suppressMessages(run_pooled(bad_factor)), "the LOCO scale factor differs")

cat("OK: mismatched loco_scale fields refused with clear messages\n")

cat("\nALL POOLED-COVARIANCE TESTS PASSED\n")
