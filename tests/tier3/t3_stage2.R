#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 3, Stage 2 side: time each phase of a Stage 2 pass over a real Stage 1
# artefact set, and the genome-wide scaling that follows from it.
#
#   Rscript --vanilla tests/tier3/t3_stage2.R <stage1_prefix> <out.json>
#
# Phases timed separately, because "Stage 2 took N seconds" is not actionable:
#   read_tsv    parsing the Stage 1 table
#   read_cov    reading and dequantising the covariance payload
#   kernel      the C++ projection / GLS solve and the Q calibration
#   write       emitting the results table
#
# The scaling run repeats the covariance payload to reach larger variant counts.
# That is synthetic and is labelled as such: it measures how Stage 2's cost and
# memory grow with the number of variants, which is the quantity a genome-wide
# plan needs, without pretending to be a genome-wide scan.
# ---------------------------------------------------------------------------
suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
pfx <- args[1L]; outjson <- args[2L]

read_cov <- gmoq_internal(".read_tau_cov_data")
tm <- function(expr) { t0 <- proc.time(); v <- force(expr)
                       list(v = v, sec = as.numeric((proc.time() - t0)[["elapsed"]])) }

r1 <- tm(data.table::fread(paste0(pfx, ".stage1.tsv.gz")))
s1 <- r1$v
beta_cols <- grep("^beta_tau_", names(s1), value = TRUE)
taus <- as.numeric(sub("^beta_tau_", "", beta_cols))
Tt <- length(taus); Mv <- nrow(s1)

r2 <- tm(read_cov(paste0(pfx, ".cov.gz"), Tt, Mv))
cov_data <- r2$v

W <- cbind(beta_mu = rep(1, Tt), beta_sigma = qnorm(taus))
Q_slope <- t(as.matrix(s1[, beta_cols, with = FALSE]))
SE_tau  <- t(as.matrix(s1[, grep("^se_tau_", names(s1)), with = FALSE]))
st <- list(taus = taus, q_tau = rep(0, Tt), Q_slope = Q_slope, SE_tau = SE_tau)

nthr <- as.integer(Sys.getenv("GMOQ_THREADS", "1"))
r3 <- tm(suppressMessages(param_gwas(
  st, transform = "custom_W", transform_args = list(W = W),
  se_mode = "tau_cov", tau_cov_estimator = "ols", cov_data = cov_data,
  q_method = "saddlepoint", n_threads = nthr,
  n_blocks = 1000, n_eff = s1$d_resid)))
fit <- r3$v

r3g <- tm(suppressMessages(param_gwas(
  st, transform = "custom_W", transform_args = list(W = W),
  se_mode = "tau_cov", tau_cov_estimator = "gls", cov_data = cov_data,
  q_method = "saddlepoint", n_threads = nthr,
  n_blocks = 1000, n_eff = s1$d_resid)))

outdt <- data.table::data.table(snp_id = s1$snp_id,
                                beta_mu = fit$params[1, ], beta_mu_se = fit$SE_params[1, ],
                                beta_sigma = fit$params[2, ], beta_sigma_se = fit$SE_params[2, ],
                                Q = fit$Q, Q_pval = fit$Q_pval)
tmpf <- tempfile(fileext = ".tsv.gz")
r4 <- tm(data.table::fwrite(outdt, tmpf, sep = "\t"))
unlink(tmpf)

# --- how the kernel scales with the variant count --------------------------
n_elem <- Tt * (Tt + 1L) / 2L
scale_rows <- list()
for (mult in c(1L, 10L, 40L)) {
  Mx <- Mv * mult
  if (Mx > 400000L) next
  cv <- list(cov_vec = rep(cov_data$cov_vec, mult),
             offsets = seq(0, by = n_elem * 4, length.out = Mx))
  Qx <- Q_slope[, rep(seq_len(Mv), mult), drop = FALSE]
  Sx <- SE_tau[, rep(seq_len(Mv), mult), drop = FALSE]
  stx <- list(taus = taus, q_tau = rep(0, Tt), Q_slope = Qx, SE_tau = Sx)
  tt <- tm(suppressMessages(param_gwas(
    stx, transform = "custom_W", transform_args = list(W = W),
    se_mode = "tau_cov", tau_cov_estimator = "ols", cov_data = cv,
    q_method = "saddlepoint", n_threads = nthr, nu = 900, nu_source = "bench")))
  scale_rows[[length(scale_rows) + 1L]] <- list(
    n_variants = Mx, seconds = tt$sec,
    cov_vec_gb = length(cv$cov_vec) * 8 / 2^30)
  rm(cv, Qx, Sx, stx, tt); invisible(gc(FALSE))
}

res <- list(
  n_variants = Mv, n_taus = Tt, n_threads = nthr,
  stage2_read_tsv_sec = r1$sec,
  stage2_read_cov_sec = r2$sec,
  stage2_kernel_ols_sec = r3$sec,
  stage2_kernel_gls_sec = r3g$sec,
  stage2_write_sec = r4$sec,
  stage2_total_sec = r1$sec + r2$sec + r3$sec + r4$sec,
  gls_solved_fraction = if (is.null(r3g$v$gls_diag)) NA_real_ else
    r3g$v$gls_diag$n_ok / max(r3g$v$gls_diag$n_snps, 1),
  peak_r_mem_gb = as.numeric(sum(gc()[, "used"] * c(8, 8))) / 2^30,
  scaling = scale_rows)

writeLines(jsonlite::toJSON(res, auto_unbox = TRUE, digits = 6, pretty = TRUE), outjson)
cat("Stage 2 timings written to", outjson, "\n")
str(res[setdiff(names(res), "scaling")])
for (r in scale_rows) cat(sprintf("  scaling: %d variants -> %.2f s\n", r$n_variants, r$seconds))
