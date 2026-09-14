# ============================================================
# Test: parameter p-values use a t reference whose degrees of freedom
#       are computed per SNP and per parameter
#
# Sigma is ESTIMATED by a delete-one-block jackknife, so theta/SE ~ t_nu.
# nu is NOT n_blocks - 1 and NOT a constant: measured
# (fungwas_paper/docs/BP1_effective_df.md)
#
#     1/nu = 1/(B - 1) + c/N_eff,   c = 0.88 location, 2.62 scale
#
# with N_eff the per-SNP d_resid Stage 1 emits.
#
# Covered here:
#   1. n_blocks round-trips from .cov.meta.json, d_resid from the TSV
#   2. <param>_df equals the law, and <param>_p equals 2*pt(-|z|, df)
#   3. nu varies across SNPs (via MAF) and across parameters (via contrast)
#   4. a Stage 1 file without d_resid falls back visibly, not silently
#   5. a covariance with no block count emits NO p-value by default
#   6. the reference and its provenance are recorded in the output
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

beta_mat <- matrix(rnorm(Tt * P, sd = 0.03), nrow = Tt, ncol = P)
# one SNP with a genuine signal, so the far tail is actually exercised
beta_mat[, 1] <- 0.30
se_mat <- matrix(rep(sqrt(diag(Sigma)), P), nrow = Tt, ncol = P)

# a realistic spread of effective sample sizes: MAF 0.01 to 0.50 at N = 461,320
N_prod <- 461320
maf_grid <- rep(c(0.01, 0.05, 0.20, 0.50), length.out = P)
d_resid <- 2 * N_prod * maf_grid * (1 - maf_grid)

write_stage1 <- function(prefix, with_d_resid = TRUE) {
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
  if (with_d_resid) dt$d_resid <- d_resid
  fwrite(dt, paste0(prefix, ".stage1.tsv.gz"), sep = "\t")
  fwrite(dt[, .(snp_id, chr, bp, effect_allele, other_allele)],
         paste0(prefix, ".cov.ids.tsv.gz"), sep = "\t")
  con <- gzfile(paste0(prefix, ".cov.gz"), "wb")
  for (j in seq_len(P)) writeBin(as.numeric(cov_upper), con, size = 4)
  close(con)
  invisible(prefix)
}

write_meta_v3 <- function(prefix, n_blocks, nonempty = n_blocks) {
  writeLines(c(
    "{",
    '  "cov_meta_version": 3,',
    '  "cov_dtype": "float32",',
    sprintf('  "n_taus": %d,', Tt),
    sprintf('  "n_cov_elements": %d,', Tt * (Tt + 1) / 2),
    '  "scale_file": false,',
    '  "jackknife": {',
    '    "method": "delete_one_block",',
    sprintf('    "n_blocks": %d,', n_blocks),
    sprintf('    "n_blocks_nonempty": %d,', nonempty),
    '    "seed": 42,',
    '    "df_model": {',
    '      "form": "1/nu = 1/(n_blocks_nonempty - 1) + c / n_eff",',
    '      "n_eff_column": "d_resid",',
    '      "c_location": 0.88,',
    '      "c_scale": 2.62,',
    '      "fallback_n_eff": 9134,',
    '      "source": "BP1_effective_df.md"',
    "    }",
    "  }",
    "}"
  ), paste0(prefix, ".cov.meta.json"))
}

write_meta_constant_nu <- function(prefix, nu) {
  # No block count at all, but a deliberately pinned constant.
  writeLines(c(
    "{",
    '  "cov_meta_version": 3,',
    '  "cov_dtype": "float32",',
    '  "jackknife": {',
    '    "method": "delete_one_block",',
    sprintf('    "nu": %s,', format(nu, scientific = FALSE)),
    '    "nu_source": "measured_effective"',
    "  }",
    "}"
  ), paste0(prefix, ".cov.meta.json"))
}

write_meta_v1 <- function(prefix) {
  writeLines(c(
    "{",
    '  "cov_dtype": "float32",',
    sprintf('  "n_taus": %d,', Tt),
    sprintf('  "n_cov_elements": %d,', Tt * (Tt + 1) / 2),
    '  "scale_file": false',
    "}"
  ), paste0(prefix, ".cov.meta.json"))
}

# W = [1, q - qbar]: row 1 of A is flat (location), row 2 sums to zero (scale)
W <- cbind(location = rep(1, Tt), spread = qnorm(taus) - mean(qnorm(taus)))
A <- solve(crossprod(W), t(W))
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

eff_nu <- fungwasStage2:::.effective_nu
C_LOC <- fungwasStage2:::C_LOCATION_DEFAULT
C_SCL <- fungwasStage2:::C_SCALE_DEFAULT
FALLBACK <- fungwasStage2:::FALLBACK_N_EFF_DEFAULT

# ------------------------------------------------------------------
# 1/2/3. per-SNP, per-parameter df from the law
# ------------------------------------------------------------------
n_blocks <- 1000L
pfx <- file.path(tmp, "b1000")
write_stage1(pfx)
write_meta_v3(pfx, n_blocks)

res <- suppressMessages(run_stage2(pfx))

stopifnot(all(c("location", "location_se", "location_p", "location_df",
                "spread", "spread_se", "spread_p", "spread_df",
                "ref_dist", "ref_df_source") %in% names(res)))
stopifnot(identical(unique(res$ref_dist), "t"))
stopifnot(identical(unique(res$ref_df_source), "per_snp_satterthwaite"))
cat("OK: n_blocks round-trips and the per-SNP law is used\n")

df_loc_exp <- eff_nu(n_blocks, d_resid, C_LOC)
df_spr_exp <- eff_nu(n_blocks, d_resid, C_SCL)
stopifnot(max(abs(res$location_df - df_loc_exp)) < 1e-9)
stopifnot(max(abs(res$spread_df - df_spr_exp)) < 1e-9)
cat(sprintf("OK: df = 1/(1/(B-1) + c/n_eff); location %.1f-%.1f, spread %.1f-%.1f\n",
            min(df_loc_exp), max(df_loc_exp), min(df_spr_exp), max(df_spr_exp)))

# these must reproduce the production numbers in the report
stopifnot(abs(round(eff_nu(1000, 2 * N_prod * 0.01 * 0.99, C_LOC)) - 911) < 1)
stopifnot(abs(round(eff_nu(1000, 2 * N_prod * 0.50 * 0.50, C_LOC)) - 995) < 1)
stopifnot(abs(round(eff_nu(1000, 2 * N_prod * 0.01 * 0.99, C_SCL)) - 776) < 1)
cat("OK: the law reproduces the production df quoted in BP1_effective_df.md\n")

# nu varies by SNP, and a scale contrast always has fewer df than a location one
stopifnot(length(unique(round(res$location_df, 6))) > 1)
stopifnot(all(res$spread_df < res$location_df))
stopifnot(all(res$location_df < n_blocks - 1))
cat(sprintf("OK: df varies across SNPs and spread < location everywhere (gap %.0f-%.0f)\n",
            min(res$location_df - res$spread_df), max(res$location_df - res$spread_df)))

# p-values match the emitted df exactly
p_loc <- 2 * pt(-abs(res$location / res$location_se), df = res$location_df)
p_spr <- 2 * pt(-abs(res$spread / res$spread_se), df = res$spread_df)
err <- max(abs(res$location_p - p_loc), abs(res$spread_p - p_spr))
stopifnot(err < 1e-14)
cat("OK: p equals 2*pt(-|theta/se|, per-SNP df) to", format(err, digits = 3), "\n")

# and the SEs are still the fixed-contrast ones the identity is exact for
stopifnot(max(abs(res$location_se - se_exp[1])) < 1e-6)
stopifnot(max(abs(res$spread_se - se_exp[2])) < 1e-6)

# heavier than normal, and heavier than a naive t_{B-1}
p_norm <- 2 * pnorm(-abs(res$location / res$location_se))
p_naive <- 2 * pt(-abs(res$location / res$location_se), df = n_blocks - 1)
stopifnot(all(res$location_p >= p_norm - 1e-15))
stopifnot(all(res$location_p >= p_naive - 1e-15))
i <- which.min(p_norm)
cat(sprintf("OK: at the sharpest SNP normal %.3e -> t_%d %.3e -> per-SNP t_%.0f %.3e\n",
            p_norm[i], n_blocks - 1L, p_naive[i], res$location_df[i], res$location_p[i]))

# Q untouched
stopifnot("Q" %in% names(res), "Q_pval" %in% names(res), all(is.finite(res$Q)))
cat("OK: Q / Q_pval still present and finite\n")

# ------------------------------------------------------------------
# 4. no d_resid column: visible fallback, not a silent one
# ------------------------------------------------------------------
pfx_nod <- file.path(tmp, "nod")
write_stage1(pfx_nod, with_d_resid = FALSE)
write_meta_v3(pfx_nod, n_blocks)

warned <- FALSE
res_nod <- withCallingHandlers(
  suppressMessages(run_stage2(pfx_nod)),
  warning = function(w) {
    if (grepl("Falling back to a flat effective df", conditionMessage(w))) warned <<- TRUE
    invokeRestart("muffleWarning")
  }
)
stopifnot(warned)
stopifnot(identical(unique(res_nod$ref_df_source), "flat_fallback_n_eff"))
stopifnot(length(unique(round(res_nod$location_df, 6))) == 1L)
stopifnot(abs(unique(res_nod$location_df) - eff_nu(n_blocks, FALLBACK, C_LOC)) < 1e-9)
# the report's recommended flat values
stopifnot(abs(round(unique(res_nod$location_df)) - 911) <= 1)
stopifnot(abs(round(unique(res_nod$spread_df)) - 776) <= 1)
cat(sprintf("OK: without d_resid the df falls back visibly to %.0f / %.0f and warns\n",
            unique(res_nod$location_df), unique(res_nod$spread_df)))

# ------------------------------------------------------------------
# 5. no block count at all: refuse to emit a p-value
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
stopifnot(!("location_p" %in% names(res_old)), !("location_df" %in% names(res_old)))
stopifnot(identical(unique(res_old$ref_dist), "none"))
stopifnot(all(is.finite(res_old$location)), all(is.finite(res_old$location_se)))
cat("OK: no block count -> no p-value, warns, ref_dist = none, rest still produced\n")

pfx_none <- file.path(tmp, "nometa")
write_stage1(pfx_none)
res_none <- suppressWarnings(suppressMessages(run_stage2(pfx_none)))
stopifnot(identical(unique(res_none$ref_dist), "none"))
cat("OK: a missing sidecar behaves identically\n")

# normal fallback only on request, and recorded
res_norm <- suppressWarnings(suppressMessages(run_stage2(pfx_old, unknown_nu = "normal")))
stopifnot(identical(unique(res_norm$ref_dist), "normal"))
stopifnot(max(abs(res_norm$location_p -
                    2 * pnorm(-abs(res_norm$location / res_norm$location_se)))) < 1e-14)
cat("OK: unknown_nu = 'normal' emits normal p-values and records ref_dist = normal\n")

stopifnot(inherits(
  try(suppressMessages(run_stage2(pfx_old, unknown_nu = "error")), silent = TRUE),
  "try-error"))
cat("OK: unknown_nu = 'error' stops rather than guessing\n")

# ------------------------------------------------------------------
# 6. explicit nu overrides everything; a recorded constant is honoured but warned
# ------------------------------------------------------------------
res_ovr <- suppressMessages(run_stage2(pfx, nu = 199))
stopifnot(identical(unique(res_ovr$ref_df_source), "user_supplied"))
stopifnot(all(res_ovr$location_df == 199), all(res_ovr$spread_df == 199))
stopifnot(max(abs(res_ovr$location_p -
                    2 * pt(-abs(res_ovr$location / res_ovr$location_se), df = 199))) < 1e-14)
cat("OK: explicit nu overrides the per-SNP law\n")

pfx_const <- file.path(tmp, "const")
write_stage1(pfx_const)
write_meta_constant_nu(pfx_const, 812.4)
warned <- FALSE
res_const <- withCallingHandlers(
  suppressMessages(run_stage2(pfx_const)),
  warning = function(w) {
    if (grepl("not \\s*defensible|constant nu", conditionMessage(w))) warned <<- TRUE
    invokeRestart("muffleWarning")
  }
)
stopifnot(warned)
stopifnot(all(res_const$location_df == 812.4))
stopifnot(identical(unique(res_const$ref_df_source), "measured_effective"))
cat("OK: a recorded constant nu is honoured but warned about\n")

# ------------------------------------------------------------------
# 7. large nu converges to the normal reference
# ------------------------------------------------------------------
res_big <- suppressMessages(run_stage2(pfx, nu = 1e8))
conv <- max(abs(res_big$location_p -
                  2 * pnorm(-abs(res_big$location / res_big$location_se))))
stopifnot(conv < 1e-7)
cat("OK: at nu = 1e8 the t p-value matches the normal one to",
    format(conv, digits = 3), "\n")

# ------------------------------------------------------------------
# 8. internals
# ------------------------------------------------------------------
flatness <- fungwasStage2:::.contrast_flatness
contrast_c <- fungwasStage2:::.contrast_c
resolve_block_count <- fungwasStage2:::.resolve_block_count
resolve_df_model <- fungwasStage2:::.resolve_df_model
resolve_nu <- fungwasStage2:::.resolve_nu
read_cov_meta <- fungwasStage2:::.read_cov_meta

stopifnot(abs(flatness(rep(1, 9)) - 1) < 1e-12)
stopifnot(abs(flatness(rep(-0.3, 9)) - 1) < 1e-12)
stopifnot(abs(flatness(qnorm(taus) - mean(qnorm(taus)))) < 1e-12)
stopifnot(is.na(flatness(rep(0, 9))))
stopifnot(identical(contrast_c(rep(1, 9)), C_LOC))
stopifnot(identical(contrast_c(qnorm(taus) - mean(qnorm(taus))), C_SCL))
# an unclassifiable contrast gets the conservative (larger c, smaller nu) value
stopifnot(identical(contrast_c(c(1, rep(0, 8))), C_SCL))
cat("OK: contrast flatness classifies location vs scale, conservative in between\n")

m3 <- read_cov_meta(paste0(pfx, ".cov.meta.json"))
stopifnot(resolve_block_count(m3) == 1000)
stopifnot(resolve_df_model(m3)$c_location == 0.88, resolve_df_model(m3)$c_scale == 2.62)
stopifnot(is.na(resolve_nu(m3)$nu))
m1 <- read_cov_meta(paste0(pfx_old, ".cov.meta.json"))
stopifnot(is.na(resolve_block_count(m1)), is.na(resolve_nu(m1)$nu))
stopifnot(is.na(resolve_block_count(NULL)))
# a populated-block count beats the requested one
stopifnot(resolve_block_count(list(jackknife = list(n_blocks = 200,
                                                    n_blocks_nonempty = 189))) == 189)
cat("OK: metadata readers behave as specified\n")

# df_coef lets a caller override c for a non location/scale design
res_c <- suppressMessages(run_stage2(pfx, df_coef = c(location = 1.5)))
stopifnot(max(abs(res_c$location_df - eff_nu(n_blocks, d_resid, 1.5))) < 1e-9)
stopifnot(max(abs(res_c$spread_df - eff_nu(n_blocks, d_resid, C_SCL))) < 1e-9)
cat("OK: df_coef overrides c per parameter\n")

# ------------------------------------------------------------------
# 9. no joint multi-parameter Wald test exists, so no F reference applies.
#    Q is a lack-of-fit statistic on T - k residual df, not a joint test of
#    k parameters. If one is added it must use F_{k, nu-k+1} (Hotelling),
#    and with a per-SNP nu at that.
# ------------------------------------------------------------------
stopifnot(!any(grepl("joint", getNamespaceExports("fungwasStage2"), ignore.case = TRUE)))
cat("NOTE: no joint k-parameter test exists; F reference not applicable\n")

cat("OK: all t-reference tests passed\n")
