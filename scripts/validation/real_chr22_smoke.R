#!/usr/bin/env Rscript
# Run OLS Stage 2 twice on the bounded real-chr22 Stage 1 smoke output.
# The pooled covariance fit is the primary result; jackknife is a paired check.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(fungwasStage2)
})

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[i], "--")) { i <- i + 1L; next }
    key <- sub("^--", "", x[i])
    if (i < length(x) && !startsWith(x[i + 1L], "--")) {
      out[[key]] <- x[i + 1L]; i <- i + 2L
    } else { out[[key]] <- TRUE; i <- i + 1L }
  }
  out
}
opts <- parse_args(args)
if (is.null(opts[["out"]]) || is.null(opts[["stage1-prefix"]])) {
  stop("usage: real_chr22_smoke.R --stage1-prefix PREFIX --truth truth.json --out SUMMARY.csv")
}

prefix <- opts[["stage1-prefix"]]
truth <- fromJSON(opts[["truth"]])
stage1_file <- paste0(prefix, ".stage1.tsv.gz")
cov_file <- paste0(prefix, ".cov.gz")
sigma0_file <- paste0(prefix, ".sigma0.tsv.gz")
cov_meta_file <- paste0(prefix, ".cov.meta.json")
for (path in c(stage1_file, cov_file, sigma0_file, cov_meta_file)) {
  if (!file.exists(path)) stop("missing Stage 1 output: ", path)
}

stage1 <- fread(stage1_file)
beta_cols <- grep("^beta_tau_", names(stage1), value = TRUE)
taus <- as.numeric(sub("^beta_tau_", "", beta_cols))
W <- make_weights_vqtl(taus, q_tau = qnorm(taus), mu = 0, sd = 1)

fit_one <- function(source) {
  fit <- if (source == "pooled") {
    param_gwas_from_file(
      stage1_file = stage1_file, W = W, se_mode = "tau_cov",
      cov_source = "pooled", sigma0_file = sigma0_file,
      cov_meta_file = cov_meta_file, tau_cov_estimator = "ols",
      q_method = "saddlepoint", n_threads = 2L
    )
  } else {
    param_gwas_from_file(
      stage1_file = stage1_file, W = W, se_mode = "tau_cov",
      cov_source = "jackknife", cov_file = cov_file,
      cov_meta_file = cov_meta_file, tau_cov_estimator = "ols",
      q_method = "saddlepoint", n_threads = 2L
    )
  }
  causal <- fit[fit$snp_id == truth$causal_snp]
  if (nrow(causal) != 1L) stop("causal SNP missing or duplicated in Stage 2 output")
  if (!is.finite(causal$beta_mu) || !is.finite(causal$beta_mu_se)) {
    stop("causal Stage 2 mean estimate is non-finite for ", source)
  }
  if (sign(causal$beta_mu) != sign(truth$expected_beta_mu)) {
    stop("causal Stage 2 mean estimate has the wrong sign for ", source)
  }
  data.table(
    covariance_source = source,
    n_stage1_variants = nrow(stage1),
    causal_snp = truth$causal_snp,
    causal_effect_per_allele = truth$expected_beta_mu,
    causal_beta_mu = causal$beta_mu,
    causal_beta_mu_se = causal$beta_mu_se,
    causal_beta_mu_p = causal$beta_mu_p,
    causal_beta_sigma2 = causal$beta_sigma2,
    causal_beta_sigma2_se = causal$beta_sigma2_se,
    causal_Q = causal$Q,
    causal_Q_p = causal$Q_pval,
    median_beta_mu = median(fit$beta_mu, na.rm = TRUE),
    median_beta_mu_se = median(fit$beta_mu_se, na.rm = TRUE),
    min_Q_p = min(fit$Q_pval, na.rm = TRUE)
  )
}

summary <- rbindlist(list(fit_one("pooled"), fit_one("jackknife")), fill = TRUE)
fwrite(summary, opts[["out"]])
cat("wrote", nrow(summary), "rows to", opts[["out"]], "\n")
