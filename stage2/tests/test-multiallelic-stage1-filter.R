# ============================================================
# Test: duplicated rsIDs are dropped before Stage 2 fitting
# ============================================================

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
this_file <- sub("^--file=", "", file_arg)
pkg_root <- normalizePath(file.path(dirname(this_file), ".."))
source(file.path(pkg_root, "R", "gwas.R"))

workdir <- tempfile("multiallelic_filter_")
dir.create(workdir)

stage1_file <- file.path(workdir, "test.stage1.tsv")
stage1_dt <- data.table(
  snp_id = c("rs_keep", "rs_multi", "rs_multi"),
  chr = c(1, 1, 1),
  bp = c(100, 200, 200),
  effect_allele = c("A", "C", "C"),
  other_allele = c("G", "A", "T"),
  beta_tau_0.25 = c(0.10, 1.00, -1.00),
  beta_tau_0.75 = c(0.20, 1.20, -1.20),
  se_tau_0.25 = c(0.05, 0.10, 0.10),
  se_tau_0.75 = c(0.05, 0.10, 0.10)
)
fwrite(stage1_dt, stage1_file, sep = "\t")

W <- matrix(c(1, 1), nrow = 2, ncol = 1)
colnames(W) <- "mean"

fit <- param_gwas_from_file(
  stage1_file = stage1_file,
  W = W,
  se_mode = "diagonal"
)

if (!identical(fit$snp_id, "rs_keep")) {
  stop("Expected only the unique rsID to remain, got: ", paste(fit$snp_id, collapse = ", "))
}

fit_unsafe <- param_gwas_from_file(
  stage1_file = stage1_file,
  W = W,
  se_mode = "diagonal",
  drop_multiallelic = FALSE
)

if (nrow(fit_unsafe) != 3L) {
  stop("Expected unsafe mode to keep all three rows")
}

cat("OK: duplicated rsIDs are dropped before Stage 2 fitting\n")
