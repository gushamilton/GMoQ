#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 1 -- Stage 1 to Stage 2, end to end, against a known theta.
#
# The fixtures were produced by tests/tier1/t1_make_fixtures.py, which runs the
# SHIPPED Stage 1 on simulated data and writes the exact artefact set the Stage
# 1 CLI writes. This script drives the SHIPPED Stage 2 entry point
# param_gwas_from_file() over them.
#
# WHAT "KNOWN THETA" MEANS HERE. Two different truths are checked and they are
# deliberately not conflated:
#
#   theta_pop     A b(tau), the projection of Stage 1's own population estimand.
#                 b(tau) is computed exactly from the conditional CDF, not by
#                 Monte Carlo. Recovery of theta_pop is a claim about the
#                 ESTIMATOR and is asserted.
#   theta_latent  the generative (a, b) that produced the data. Recovery of it
#                 also requires the linearisation to hold at this effect size,
#                 so it is measured and reported rather than asserted tightly.
#
# Run:  Rscript --vanilla tests/tier1/t1_end_to_end.R
# (needs $GMOQ_WORK/tier1 populated -- tests/run_tier1.sh does both)
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier1/t1_end_to_end.R")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
suppressPackageStartupMessages(library(data.table))
gmoq_begin("tier1", "end_to_end")

FIX <- Sys.getenv("GMOQ_TIER1_DIR", file.path(Sys.getenv("GMOQ_WORK"), "tier1"))
if (!dir.exists(FIX)) stop("Fixtures not found in ", FIX,
                           ". Run tests/tier1/t1_make_fixtures.py first.")

meta  <- jsonlite::fromJSON(file.path(FIX, "fixtures.json"))
taus  <- as.numeric(meta$taus)
Tt    <- length(taus)
b_true <- as.matrix(data.table::fread(file.path(FIX, "b_true.tsv"), header = FALSE))
causal <- seq.int(meta$n_null + 1L, meta$n_null + meta$n_causal)   # 1-based
lat <- data.table::fread(file.path(FIX, "theta_latent.tsv"))

# The fitted model: location + scale on the Cornish-Fisher columns, which is
# exactly the family the fixtures were simulated under.
W <- cbind(beta_mu = rep(1, Tt), beta_sigma = qnorm(taus))
rownames(W) <- paste0("tau", taus)
A <- solve(crossprod(W), t(W))
theta_pop <- A %*% t(b_true)          # 2 x n_variants

pfx <- function(dir) file.path(FIX, dir, "sim")

# =====================================================================
# 1. The jackknife path, OLS -- the shipped default
# =====================================================================
t0 <- proc.time()[["elapsed"]]
res_ols <- suppressMessages(param_gwas_from_file(
  stage1_file = paste0(pfx("clean"), ".stage1.tsv.gz"),
  W = W, se_mode = "tau_cov", cov_source = "jackknife",
  cov_file = paste0(pfx("clean"), ".cov.gz"),
  tau_cov_estimator = "ols", q_method = "saddlepoint"))
t_ols <- proc.time()[["elapsed"]] - t0

gmoq_check("T1.E2E.output_shape",
           "param_gwas_from_file returns one row per Stage 1 variant with the fitted parameters.",
           nrow(res_ols) == nrow(b_true) &&
             all(c("beta_mu", "beta_mu_se", "beta_mu_p", "beta_mu_df",
                   "beta_sigma", "beta_sigma_se", "Q", "Q_pval",
                   "cov_source", "ref_dist", "ref_df_source") %in% names(res_ols)),
           criterion = "row count and the documented column set")

gmoq_check("T1.E2E.reference_is_t_per_snp",
           "The reference distribution is t with a per-SNP, per-parameter nu, not a normal.",
           all(res_ols$ref_dist == "t") &&
             all(res_ols$ref_df_source == "per_snp_satterthwaite"),
           criterion = "ref_dist == 't', ref_df_source == 'per_snp_satterthwaite'")
gmoq_check("T1.E2E.df_law_applied",
           "The emitted degrees of freedom follow 1/nu = 1/(B-1) + c/d_resid with the shipped c.",
           {
             s1 <- data.table::fread(paste0(pfx("clean"), ".stage1.tsv.gz"))
             B <- meta$blocks
             exp_mu <- 1 / (1 / (B - 1) + 0.88 / s1$d_resid)
             exp_sg <- 1 / (1 / (B - 1) + 2.62 / s1$d_resid)
             max(abs(res_ols$beta_mu_df - exp_mu)) < 1e-6 &&
               max(abs(res_ols$beta_sigma_df - exp_sg)) < 1e-6
           },
           criterion = "per-SNP nu matches the law for both parameters")

# --- recovery of theta_pop -------------------------------------------------
zmat <- rbind((res_ols$beta_mu    - theta_pop[1, ]) / res_ols$beta_mu_se,
              (res_ols$beta_sigma - theta_pop[2, ]) / res_ols$beta_sigma_se)
gmoq_check("T1.E2E.recovers_theta_pop_causal",
           "Every causal variant's estimate is within 5 standard errors of the exact population theta.",
           max(abs(zmat[, causal])) < 5,
           stat = max(abs(zmat[, causal])),
           criterion = "max |theta_hat - theta_pop| / se < 5 over 2 x n_causal")
gmoq_check("T1.E2E.recovers_theta_pop_all",
           "Across all variants the standardised deviation from theta_pop has unit scale and no drift.",
           abs(mean(zmat)) < 4 / sqrt(length(zmat)) && abs(sd(as.numeric(zmat)) - 1) < 0.25,
           stat = sd(as.numeric(zmat)),
           criterion = "|mean z| < 4/sqrt(n) and |sd(z) - 1| < 0.25")

# --- recovery of the latent generative parameter ---------------------------
lat_z <- rbind((res_ols$beta_mu    - lat$a) / res_ols$beta_mu_se,
               (res_ols$beta_sigma - lat$b) / res_ols$beta_sigma_se)
gmoq_check("T1.E2E.recovers_theta_latent",
           "The generative (a, b) is also recovered, so the linearisation holds at this effect size.",
           max(abs(lat_z[, causal])) < 6, stat = max(abs(lat_z[, causal])),
           criterion = "max |theta_hat - theta_latent| / se < 6 over the causal variants")
gmoq_check("T1.E2E.estimand_gap",
           "Size of the gap between the population estimand and the generative parameter -- reported.",
           TRUE, severity = "report",
           stat = max(abs(theta_pop[, causal] - rbind(lat$a[causal], lat$b[causal])) /
                        rbind(res_ols$beta_mu_se[causal], res_ols$beta_sigma_se[causal])),
           criterion = "max |theta_pop - theta_latent| / se, the second-order residue")

# --- the nulls -------------------------------------------------------------
null_idx <- seq_len(meta$n_null)
pnull <- c(res_ols$beta_mu_p[null_idx], res_ols$beta_sigma_p[null_idx])
ks <- suppressWarnings(stats::ks.test(pnull, "punif")$p.value)
gmoq_check("T1.E2E.null_pvalues_uniform",
           "Null-variant parameter p-values are uniform, so the reference distribution is not obviously wrong.",
           ks > 0.001, stat = ks,
           criterion = "Kolmogorov-Smirnov p > 0.001 against Uniform(0,1)",
           detail = "a weak check at this sample size; tier 2 counts the rate properly")

gmoq_check("T1.E2E.Q_pvalue_valid",
           "Every Q p-value is a probability, and the null variants' are not systematically small.",
           all(is.finite(res_ols$Q_pval)) && all(res_ols$Q_pval >= 0 & res_ols$Q_pval <= 1) &&
             mean(res_ols$Q_pval[null_idx] < 0.05) < 0.20,
           stat = mean(res_ols$Q_pval[null_idx] < 0.05),
           criterion = "Q_pval in [0,1]; null rejection at 0.05 below 0.20")

# =====================================================================
# 2. The jackknife path, GLS
# =====================================================================
res_gls <- suppressMessages(param_gwas_from_file(
  stage1_file = paste0(pfx("clean"), ".stage1.tsv.gz"),
  W = W, se_mode = "tau_cov", cov_source = "jackknife",
  cov_file = paste0(pfx("clean"), ".cov.gz"),
  tau_cov_estimator = "gls", q_method = "saddlepoint"))
zg <- rbind((res_gls$beta_mu    - theta_pop[1, ]) / res_gls$beta_mu_se,
            (res_gls$beta_sigma - theta_pop[2, ]) / res_gls$beta_sigma_se)
gmoq_check("T1.E2E.gls_recovers_theta_pop",
           "The GLS estimator also recovers the exact population theta.",
           max(abs(zg[, causal])) < 5, stat = max(abs(zg[, causal])),
           criterion = "max |theta_hat - theta_pop| / se < 5")
gmoq_check("T1.E2E.gls_is_tighter",
           "GLS standard errors are smaller than OLS's, the efficiency claim in docs/OLS_VS_GLS.md.",
           median(res_gls$beta_sigma_se / res_ols$beta_sigma_se) < 1.0,
           stat = median(res_gls$beta_sigma_se / res_ols$beta_sigma_se),
           criterion = "median SE ratio GLS/OLS < 1")

# =====================================================================
# 3. The default really is OLS
# =====================================================================
gmoq_check("T1.E2E.default_estimator_is_ols",
           "tau_cov_estimator defaults to 'ols', which docs/OLS_VS_GLS.md concluded it must.",
           identical(eval(formals(param_gwas_from_file)$tau_cov_estimator)[1], "ols") &&
             identical(eval(formals(param_gwas)$tau_cov_estimator)[1], "ols"),
           criterion = "first choice of tau_cov_estimator is 'ols' in both entry points")
gmoq_check("T1.E2E.default_cov_source_is_jackknife",
           "cov_source defaults to 'jackknife': every published result used it and this is not a one-way door.",
           identical(eval(formals(param_gwas_from_file)$cov_source)[1], "jackknife"),
           criterion = "first choice of cov_source is 'jackknife'")

# =====================================================================
# 4. The pooled covariance path
# =====================================================================
res_pool <- suppressMessages(param_gwas_from_file(
  stage1_file = paste0(pfx("clean"), ".stage1.tsv.gz"),
  W = W, se_mode = "tau_cov", cov_source = "pooled",
  sigma0_file = paste0(pfx("clean"), ".sigma0.tsv.gz"),
  cov_meta_file = paste0(pfx("clean"), ".cov.meta.json"),
  tau_cov_estimator = "ols", q_method = "saddlepoint"))

gmoq_check("T1.E2E.pooled_same_point_estimate",
           "The covariance source changes the standard error and nothing else: OLS theta is identical.",
           max(abs(res_pool$beta_mu - res_ols$beta_mu)) < 1e-12 &&
             max(abs(res_pool$beta_sigma - res_ols$beta_sigma)) < 1e-12,
           stat = max(abs(res_pool$beta_mu - res_ols$beta_mu)),
           criterion = "OLS point estimates bit-identical between cov sources")
gmoq_check("T1.E2E.pooled_reference_recorded",
           "The pooled path records its own df source, so a results file says which covariance made it.",
           all(res_pool$cov_source == "pooled") &&
             all(res_pool$ref_df_source == "pooled_sigma0_n_samples"),
           criterion = "cov_source == 'pooled', ref_df_source == 'pooled_sigma0_n_samples'")
gmoq_check("T1.E2E.pooled_vs_jackknife_se",
           "Pooled and jackknife standard errors agree to within the jackknife's own noise at B = 200.",
           abs(median(res_pool$beta_mu_se / res_ols$beta_mu_se) - 1) < 0.15,
           stat = median(res_pool$beta_mu_se / res_ols$beta_mu_se),
           criterion = "|median SE ratio - 1| < 0.15",
           detail = "the jackknife SE is ~6 per cent conservative by construction; see docs/BP2 section 1")
gmoq_check("T1.E2E.pooled_recovers_theta_pop",
           "The pooled covariance also gives intervals that cover the exact population theta.",
           max(abs((res_pool$beta_sigma - theta_pop[2, causal[1]] * 0 - theta_pop[2, ])[causal] /
                     res_pool$beta_sigma_se[causal])) < 5,
           stat = max(abs((res_pool$beta_sigma - theta_pop[2, ])[causal] /
                            res_pool$beta_sigma_se[causal])),
           criterion = "max |theta_hat - theta_pop| / se < 5 on the scale parameter")

# =====================================================================
# 5. The int8 covariance reader
# =====================================================================
res_i8 <- suppressMessages(param_gwas_from_file(
  stage1_file = paste0(pfx("int8"), ".stage1.tsv.gz"),
  W = W, se_mode = "tau_cov", cov_source = "jackknife",
  cov_file = paste0(pfx("int8"), ".cov.gz"),
  tau_cov_estimator = "ols", q_method = "saddlepoint"))
gmoq_check("T1.E2E.int8_reader_selected",
           "The presence of a .cov.scale.gz sidecar switches the reader to int8 and it dequantises correctly.",
           median(abs(res_i8$beta_mu_se / res_ols$beta_mu_se - 1)) < 0.02,
           stat = median(abs(res_i8$beta_mu_se / res_ols$beta_mu_se - 1)),
           criterion = "median |SE ratio - 1| < 0.02 against the float32 payload")

# =====================================================================
# 6. The parameter covariance output
# =====================================================================
out_cov <- suppressMessages(param_gwas_from_file(
  stage1_file = paste0(pfx("clean"), ".stage1.tsv.gz"),
  W = W, se_mode = "tau_cov", cov_source = "jackknife",
  cov_file = paste0(pfx("clean"), ".cov.gz"),
  tau_cov_estimator = "ols", return_cov = "upper", q_method = "saddlepoint"))
pcu <- out_cov$param_cov_upper
gmoq_check("T1.E2E.param_cov_upper",
           "return_cov='upper' emits the parameter covariance, named and consistent with the SEs.",
           !is.null(pcu) &&
             identical(colnames(pcu),
                       c("cov_beta_mu_beta_mu", "cov_beta_mu_beta_sigma",
                         "cov_beta_sigma_beta_sigma")) &&
             max(abs(sqrt(pcu[, 1]) - res_ols$beta_mu_se)) < 1e-10,
           criterion = "documented names; sqrt of the diagonal equals the reported SE")

gmoq_check("T1.E2E.runtime_seconds",
           "Wall clock for one Stage 2 pass over the fixture -- reported.",
           TRUE, stat = t_ols, severity = "report",
           criterion = sprintf("%d variants x %d taus, jackknife OLS", nrow(b_true), Tt))

gmoq_finish()
