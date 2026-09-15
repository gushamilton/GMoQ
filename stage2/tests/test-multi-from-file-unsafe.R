# ============================================================
# Test: param_gwas_multi_from_file requires cov IDs unless unsafe
# ============================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rcpp))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
this_file <- sub("^--file=", "", file_arg)
pkg_root <- normalizePath(file.path(dirname(this_file), ".."))
source(file.path(pkg_root, "R", "tail-reference.R"))
source(file.path(pkg_root, "R", "gwas.R"))
source(file.path(pkg_root, "R", "param_gwas.R"))
source(file.path(pkg_root, "R", "RcppExports.R"))
Rcpp::sourceCpp(file.path(pkg_root, "src", "models.cpp"))

workdir <- tempfile("multi_from_file_")
dir.create(workdir)

stage1_file <- file.path(workdir, "test.stage1.tsv")
stage1_dt <- data.table(
  snp_id = c("rs1", "rs2"),
  chr = c(1, 1),
  bp = c(100, 200),
  effect_allele = c("A", "C"),
  other_allele = c("G", "T"),
  beta_tau_0.1 = c(0.0, 0.0),
  se_tau_0.1 = c(1.0, 1.0)
)

fwrite(stage1_dt, stage1_file, sep = "\t")

cov_file <- file.path(workdir, "test.cov.gz")
con <- gzfile(cov_file, "wb")
writeBin(raw(0), con)
close(con)

W_list <- list(add = matrix(1, nrow = 1, ncol = 1))

err <- tryCatch({
  param_gwas_multi_from_file(
    stage1_file = stage1_file,
    W_list = W_list,
    N = 100,
    se_mode = "tau_cov",
    cov_source = "jackknife",
    cov_file = cov_file
  )
  NULL
}, error = function(e) e)

if (is.null(err)) {
  stop("Expected error for missing covariance IDs sidecar, but got none")
}

msg <- conditionMessage(err)
if (!grepl("Missing covariance IDs file", msg, fixed = TRUE)) {
  stop("Unexpected error message: ", msg)
}

# Now proceed with unsafe flag (should warn but not error)
suppressWarnings(
  param_gwas_multi_from_file(
    stage1_file = stage1_file,
    W_list = W_list,
    N = 100,
    se_mode = "tau_cov",
    cov_source = "jackknife",
    cov_file = cov_file,
    unsafe_skip_cov_ids = TRUE
  )
)

cat("OK: unsafe_skip_cov_ids bypasses missing cov IDs check\n")

# A covariance IDs sidecar is indexed against the raw Stage 1 rows. The default
# multiallelic filter drops every row in a duplicated rsID/position group, so
# validation must happen before filtering and the covariance offsets must then
# be subset to the retained raw indices.
stage1_filtered_file <- file.path(workdir, "filtered.stage1.tsv")
stage1_filtered_dt <- data.table(
  snp_id = c("rs_keep", "rs_multi", "rs_multi"),
  chr = c(1L, 1L, 1L),
  bp = c(100L, 200L, 200L),
  effect_allele = c("A", "C", "G"),
  other_allele = c("G", "T", "A"),
  beta_tau_0.1 = c(0.4, 0.2, 0.8),
  beta_tau_0.9 = c(0.1, 0.9, 0.1),
  se_tau_0.1 = c(1, 1, 1),
  se_tau_0.9 = c(1, 1, 1)
)
fwrite(stage1_filtered_dt, stage1_filtered_file, sep = "\t")

filtered_cov_file <- file.path(workdir, "filtered.cov.gz")
con <- gzfile(filtered_cov_file, "wb")
# T = 2, row-major upper triangles (1,1), (1,2), (2,2), one per raw row.
writeBin(as.numeric(c(1, 0, 1, 2, 0, 2, 3, 0, 3)), con, size = 4)
close(con)
fwrite(stage1_filtered_dt[, .(snp_id, chr, bp, effect_allele, other_allele)],
       sub("\\.cov\\.gz$", ".cov.ids.tsv.gz", filtered_cov_file), sep = "\t")

filtered_fit <- param_gwas_multi_from_file(
  stage1_file = stage1_filtered_file,
  W_list = list(add = matrix(c(1, 1), nrow = 2, ncol = 1)),
  N = 100,
  se_mode = "tau_cov",
  cov_source = "jackknife",
  cov_file = filtered_cov_file
)

# Only rs_keep survives. Its beta vector is (0.4, 0.1), giving Q = 0.045
# against the first covariance record and hence the known AIC/BIC below.
if (nrow(filtered_fit) != 1L || !identical(rownames(filtered_fit), NULL)) {
  stop("Expected one retained row in the multi-file jackknife result")
}
if (!identical(colnames(filtered_fit), c("AIC_add", "BIC_add"))) {
  stop("Unexpected model-score columns: ", paste(colnames(filtered_fit), collapse = ", "))
}
expected_scores <- c(AIC_add = 2.045, BIC_add = 0.045 + log(100))
if (max(abs(as.numeric(filtered_fit[1, ]) - expected_scores)) > 1e-6) {
  stop("Retained covariance record was not used: got ",
       paste(format(as.numeric(filtered_fit[1, ]), digits = 12), collapse = ", "))
}
cat("OK: multi-file jackknife validates raw IDs and subsets offsets after filtering\n")

# A mismatched raw-row ID must still refuse before any covariance record is used.
wrong_ids <- copy(stage1_filtered_dt[, .(snp_id, chr, bp, effect_allele, other_allele)])
wrong_ids$snp_id[1L] <- "rs_wrong"
fwrite(wrong_ids,
       sub("\\.cov\\.gz$", ".cov.ids.tsv.gz", filtered_cov_file), sep = "\t")
err_mismatch <- tryCatch({
  param_gwas_multi_from_file(
    stage1_file = stage1_filtered_file,
    W_list = list(add = matrix(c(1, 1), nrow = 2, ncol = 1)),
    N = 100,
    se_mode = "tau_cov",
    cov_source = "jackknife",
    cov_file = filtered_cov_file
  )
  NULL
}, error = function(e) e)
if (is.null(err_mismatch) ||
    !grepl("Covariance IDs mismatch", conditionMessage(err_mismatch), fixed = TRUE)) {
  stop("Expected mismatched raw covariance IDs to refuse")
}
cat("OK: mismatched raw covariance IDs still refuse\n")
