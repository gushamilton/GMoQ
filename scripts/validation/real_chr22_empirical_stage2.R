#!/usr/bin/env Rscript
# Fit the paper additive model to the bounded empirical log_bmi Stage 1 run.
# The pooled covariance is the primary covariance source for this smoke.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(fungwasStage2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument: ", flag)
  args[[i + 1L]]
}

stage1_file <- get_arg("--stage1")
sigma0_file <- get_arg("--sigma0")
cov_meta_file <- get_arg("--cov-meta")
weights_file <- get_arg("--weights")
baseline_file <- if ("--baseline-json" %in% args) get_arg("--baseline-json") else NULL
out_file <- get_arg("--out")

header <- fread(stage1_file, nrows = 0L)
beta_cols <- grep("^beta_tau_", names(header), value = TRUE)
taus <- as.numeric(sub("^beta_tau_", "", beta_cols))
if (!is.null(baseline_file)) {
  baseline <- fromJSON(baseline_file)
  W <- make_weights_vqtl(
    as.numeric(baseline$taus), as.numeric(baseline$q_tau),
    mu = as.numeric(baseline$mean), sd = as.numeric(baseline$sd)
  )
  w_taus <- NULL
} else {
  wobj <- readRDS(weights_file)
  if (is.list(wobj) && !is.null(wobj$W)) {
    W <- wobj$W
    w_taus <- wobj$taus
  } else {
    W <- wobj
    w_taus <- NULL
  }
}
if (!is.null(w_taus)) {
  if (nrow(W) != length(w_taus) && ncol(W) == length(w_taus)) W <- t(W)
  idx <- vapply(taus, function(tau) {
    d <- abs(w_taus - tau)
    j <- which.min(d)
    if (d[j] > 1e-8) NA_integer_ else j
  }, integer(1))
  if (anyNA(idx)) stop("weight taus do not cover Stage 1 taus")
  W <- W[idx, , drop = FALSE]
}
if (is.null(colnames(W))) colnames(W) <- paste0("param", seq_len(ncol(W)))
W <- as.matrix(W)

fit <- param_gwas_from_file(
  stage1_file = stage1_file,
  W = W,
  se_mode = "tau_cov",
  cov_source = "pooled",
  sigma0_file = sigma0_file,
  cov_meta_file = cov_meta_file,
  tau_cov_estimator = "ols",
  q_method = "saddlepoint",
  n_threads = 2L
)
fwrite(fit, out_file, sep = "\t")
cat("wrote", nrow(fit), "rows to", out_file, "with", ncol(W), "parameters\n")
