# ============================================================
# Test: duplicated rsID/position rows are dropped before Stage 2
# ============================================================

library(data.table)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
this_file <- sub("^--file=", "", file_arg)
pkg_root <- normalizePath(file.path(dirname(this_file), ".."))
source(file.path(pkg_root, "R", "gwas.R"))

stage1_dt <- data.table(
  snp_id = c("rs_dup", "rs_dup", "rs_pos_a", "rs_pos_b", "rs_keep"),
  chr = c(1L, 1L, 2L, 2L, 3L),
  bp = c(100L, 100L, 200L, 200L, 300L),
  effect_allele = c("A", "A", "C", "G", "T"),
  other_allele = c("C", "G", "A", "A", "C"),
  beta_tau_0.5 = 0,
  se_tau_0.5 = 1
)

filtered <- .filter_multiallelic_stage1(stage1_dt)

if (nrow(filtered$dt) != 1L || filtered$dt$snp_id != "rs_keep") {
  stop("Expected only the unambiguous rs_keep row to remain")
}
if (!identical(filtered$keep_idx, 5L)) {
  stop("Expected keep_idx to preserve the raw Stage 1 row index")
}
if (nrow(filtered$dropped) != 4L) {
  stop("Expected four duplicated/multi-allelic rows to be reported")
}

cat("OK: duplicated rsID and duplicated-position rows are dropped before Stage 2\n")
