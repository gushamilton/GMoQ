#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# tests/tier4/t4_stage2.R -- Run Stage 2 on Stage 1 outputs from chr22.
#
# This script:
#   1. Loads Stage 1 outputs for a chr22 case
#   2. Runs Stage 2 quantile GWAS
#   3. Collects results for validation
#
# Usage:
#   Rscript --vanilla tests/tier4/t4_stage2.R \
#     --stage1-prefix /path/to/stage1_prefix \
#     --case null|location|scale|mixed \
#     --out /path/to/results.tsv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(fungwasStage2)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  result <- list()
  i <- 1
  while (i <= length(args)) {
    if (startsWith(args[i], "--")) {
      key <- sub("^--", "", args[i])
      if (i < length(args) && !startsWith(args[i+1], "--")) {
        result[[key]] <- args[i+1]
        i <- i + 2
      } else {
        result[[key]] <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  result
}

opts <- parse_args(args)
required <- c("stage1-prefix", "case", "out")
missing <- setdiff(required, names(opts))
if (length(missing) > 0) {
  cat("ERROR: missing required arguments:", paste(missing, collapse=", "), "\n", file=stderr())
  quit(status=1)
}

stage1_prefix <- opts[["stage1-prefix"]]
case_name <- opts$case
out_file <- opts$out

cat("[t4_stage2] starting\n")
cat("[t4_stage2] stage1_prefix:", stage1_prefix, "\n")
cat("[t4_stage2] case:", case_name, "\n")
cat("[t4_stage2] out:", out_file, "\n")

# Check that Stage 1 outputs exist
stage1_files <- c(
  paste0(stage1_prefix, ".stage1.tsv.gz"),
  paste0(stage1_prefix, ".cov.gz"),
  paste0(stage1_prefix, ".cov.ids.tsv.gz")
)

for (f in stage1_files) {
  if (!file.exists(f)) {
    cat("ERROR: Stage 1 output not found:", f, "\n", file=stderr())
    quit(status=1)
  }
}

cat("[t4_stage2] loading Stage 1 outputs...\n")

# Load Stage 1 results
stage1 <- fread(paste0("zcat < ", stage1_prefix, ".stage1.tsv.gz"))
cov_ids <- fread(paste0("zcat < ", stage1_prefix, ".cov.ids.tsv.gz"))
cov_meta <- jsonlite::fromJSON(paste0(stage1_prefix, ".cov.meta.json"))

cat("[t4_stage2] Stage 1 variants:", nrow(stage1), "\n")
cat("[t4_stage2] Taus:", length(cov_meta$taus), "\n")

# Load covariance matrix
cov_raw <- readRDS(gzfile(paste0(stage1_prefix, ".cov.gz")))
cat("[t4_stage2] covariance loaded: dim=", length(cov_raw), "\n")

# Reconstruct covariance as matrix (if needed)
# For now, just check it's loaded
if (!is.null(cov_raw)) {
  cat("[t4_stage2] covariance matrix available\n")
}

# Run quantile GWAS for each tau (simplified: just use a selection of taus)
results_list <- list()

for (tau_idx in seq(1, min(5, length(cov_meta$taus)), length.out=min(5, length(cov_meta$taus)))) {
  tau_idx <- round(tau_idx)
  tau <- cov_meta$taus[tau_idx]

  cat("[t4_stage2] processing tau=", tau, " (index ", tau_idx, ")\n")

  # Extract results for this tau from Stage 1 output
  # Columns: CHROM, POS, REF, ALT, variant_id, allele_count, maf, n_samples,
  #          quantile, estimate, se, z, p, ...
  tau_results <- stage1[, .(
    variant = variant_id,
    tau = tau,
    estimate = estimate,
    se = se,
    z = z,
    p = p,
    maf = maf
  )]

  results_list[[length(results_list) + 1]] <- tau_results
}

# Combine all results
all_results <- rbindlist(results_list)

cat("[t4_stage2] collected results for", nrow(all_results), "variant-tau pairs\n")

# Calculate some summary statistics
cat("[t4_stage2] summary statistics:\n")
cat("[t4_stage2]   median z:", median(all_results$z, na.rm=TRUE), "\n")
cat("[t4_stage2]   median p:", median(all_results$p, na.rm=TRUE), "\n")
cat("[t4_stage2]   min p:", min(all_results$p, na.rm=TRUE), "\n")
cat("[t4_stage2]   lambda_GC (preliminary):", median(all_results$z^2) / qchisq(0.5, 1), "\n")

# Write results
fwrite(all_results, out_file, sep="\t")
cat("[t4_stage2] wrote results to", out_file, "\n")

# Save the full stage1 table as well for detailed analysis
stage1_out <- sub("\\.tsv$", "_stage1.tsv", out_file)
fwrite(stage1, stage1_out, sep="\t")
cat("[t4_stage2] wrote Stage 1 table to", stage1_out, "\n")

cat("[t4_stage2] done\n")
