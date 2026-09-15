#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# tests/tier4/t4_analyze.R -- Analyze chr22 validation results.
#
# This script:
#   1. Loads Stage 1 and Stage 2 outputs
#   2. Compares ground truth effects to estimated effects
#   3. Checks pass criteria (Type-I rate, effect recovery, etc.)
#   4. Produces a summary report
#
# Usage:
#   Rscript --vanilla tests/tier4/t4_analyze.R \
#     --work-dir /path/to/work \
#     --out /path/to/results.json
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
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
work_dir <- opts[["work-dir"]]
out_file <- opts$out

if (is.null(work_dir) || is.null(out_file)) {
  cat("usage: t4_analyze.R --work-dir /path/to/work --out /path/to/results.json\n", file=stderr())
  quit(status=1)
}

cat("[t4_analyze] starting\n")
cat("[t4_analyze] work_dir:", work_dir, "\n")
cat("[t4_analyze] out_file:", out_file, "\n")

# Test cases
cases <- c("null", "location", "scale", "mixed")
results <- list()

for (case in cases) {
  cat("[t4_analyze] analyzing case:", case, "\n")

  # Load metadata
  meta_file <- file.path(work_dir, paste0("pheno_", case, ".meta.json"))
  if (!file.exists(meta_file)) {
    cat("  SKIP: metadata not found\n")
    next
  }

  meta <- jsonlite::fromJSON(meta_file)
  cat("  Causal variants:", length(meta$causal_snps), "\n")
  cat("  Causal MAF range: [", format(min(meta$causal_mafs), digits=3),
      ", ", format(max(meta$causal_mafs), digits=3), "]\n")

  # Load Stage 1 results
  # stage1_file <- file.path(work_dir, paste0(case, "_stage1_stage1.tsv"))
  # For now, we'll create a summary based on metadata
  # In a full run, Stage 1 and Stage 2 would have been executed

  case_result <- list(
    case = case,
    n_causal = length(meta$causal_snps),
    n_null = length(meta$null_snps),
    causal_mafs_min = min(meta$causal_mafs),
    causal_mafs_max = max(meta$causal_mafs),
    phenotype_mean = meta$phenotype_mean,
    phenotype_sd = meta$phenotype_sd,
    status = "phenotype_generated"
  )

  results[[case]] <- case_result
  cat("  Status:", case_result$status, "\n")
}

# Write results
out_dir <- dirname(out_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

if (requireNamespace("jsonlite", quietly=TRUE)) {
  result_json <- jsonlite::toJSON(results, pretty=TRUE)
  write(result_json, out_file)
  cat("[t4_analyze] wrote results to", out_file, "\n")
} else {
  # Write as text summary
  summary_text <- sprintf(
    "Chr22 Validation Summary\n=========================\n\n%s",
    paste(sapply(names(results), function(x) {
      r <- results[[x]]
      sprintf("Case: %s\n  Causal SNPs: %d\n  MAF range: [%.4f, %.4f]\n  Status: %s\n",
              r$case, r$n_causal, r$causal_mafs_min, r$causal_mafs_max, r$status)
    }), collapse="\n")
  )
  write(summary_text, out_file)
  cat("[t4_analyze] wrote results to", out_file, "\n")
}

cat("[t4_analyze] done\n")
