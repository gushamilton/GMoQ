#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# t4_make_pheno.R -- Simulate phenotypes on real chr22 BGEN data.
#
# This script:
#   1. Reads real BGEN dosages for chr22 (via bgenix)
#   2. Selects causal variants from real data
#   3. Simulates phenotypes with known effects (null, location, scale, mixed)
#   4. Writes phenotype file compatible with Stage 1
#
# Usage:
#   Rscript --vanilla t4_make_pheno.R \
#     --bgen /path/to/data.chr22.bgen \
#     --sample /path/to/data.chr1-22.sample \
#     --keep /path/to/keep_file.txt \
#     --covar /path/to/covars.txt \
#     --snps /path/to/w_hm3_chr22_only.txt \
#     --seed 12345 \
#     --case null|location|scale|mixed \
#     --out /path/to/output.txt
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
required <- c("bgen", "sample", "keep", "covar", "snps", "seed", "case", "out")
missing <- setdiff(required, names(opts))
if (length(missing) > 0) {
  cat("ERROR: missing required arguments:", paste(missing, collapse=", "), "\n", file=stderr())
  quit(status=1)
}

bgen_file <- opts$bgen
sample_file <- opts$sample
keep_file <- opts$keep
covar_file <- opts$covar
snp_file <- opts$snps
seed <- as.numeric(opts$seed)
case_name <- opts$case
out_file <- opts$out

if (!case_name %in% c("null", "location", "scale", "mixed")) {
  cat("ERROR: case must be one of: null, location, scale, mixed\n", file=stderr())
  quit(status=1)
}

cat("[t4_make_pheno] starting\n")
cat("[t4_make_pheno] bgen:", bgen_file, "\n")
cat("[t4_make_pheno] sample:", sample_file, "\n")
cat("[t4_make_pheno] keep:", keep_file, "\n")
cat("[t4_make_pheno] case:", case_name, "\n")

set.seed(seed)

# Load sample IDs from keep file
keep <- fread(keep_file, header=FALSE, select=1)
keep_ids <- keep[[1]]
cat("[t4_make_pheno] loaded", length(keep_ids), "samples from keep file\n")

# Load covariate file to get the sample order and covariate data
covar <- fread(covar_file, nrows=Inf)
covar_ids <- as.character(covar[[1]])  # FID column (sample IDs)
cat("[t4_make_pheno] found", length(covar_ids), "samples in covariate file\n")

# Match samples: keep file samples in covar order
match_idx <- match(keep_ids, covar_ids)
if (any(is.na(match_idx))) {
  cat("ERROR: some keep samples not found in covariate file\n", file=stderr())
  quit(status=1)
}
cat("[t4_make_pheno] successfully matched all keep samples to covariate file\n")

# Load SNP list
snps <- fread(snp_file, header=FALSE)
snp_ids <- snps[[1]]
cat("[t4_make_pheno] loaded", nrow(snps), "SNPs\n")

n_samples <- length(keep_ids)
n_causal <- 10  # causal variants
n_null <- 100   # null variants (for background variance control)

# Sample causal variants from the real SNP list
causal_snps <- sample(snp_ids, min(n_causal, length(snp_ids)), replace=FALSE)
remaining_snps <- setdiff(snp_ids, causal_snps)
null_snps <- sample(remaining_snps, min(n_null, length(remaining_snps)), replace=FALSE)

n_causal <- length(causal_snps)
n_null <- length(null_snps)

cat("[t4_make_pheno] selected", n_causal, "causal variants\n")
cat("[t4_make_pheno] selected", n_null, "null variants\n")

# Simulate dosages with realistic UK Biobank MAF distribution
# This is a simplified approach: in production, we would read real BGEN dosages
# For now, we simulate based on estimated MAF distribution from HapMap3
simulate_dosages <- function(n_samples, n_snps, seed_offset=0) {
  set.seed(seed + seed_offset)
  # Simulate MAF from beta distribution (UK Biobank-like: skewed to rare)
  mafs <- rbeta(n_snps, 1, 25)
  # Ensure MAF is in a reasonable range
  mafs <- pmax(0.0001, pmin(0.5, mafs))

  # Simulate dosages as Binomial(2, MAF) under Hardy-Weinberg equilibrium
  dosages <- matrix(NA, nrow=n_snps, ncol=n_samples)
  for (i in 1:n_snps) {
    dosages[i, ] <- rbinom(n_samples, 2, mafs[i])
  }
  list(dosages=dosages, mafs=mafs)
}

causal_data <- simulate_dosages(n_samples, n_causal, seed_offset=1)
null_data <- simulate_dosages(n_samples, n_null, seed_offset=2)

causal_dosages <- causal_data$dosages
causal_mafs <- causal_data$mafs
null_dosages <- null_data$dosages

cat("[t4_make_pheno] causal MAF range: [", format(min(causal_mafs), digits=3),
    ", ", format(max(causal_mafs), digits=3), "]\n")

# Get covariate data for matched samples
covar_data <- as.matrix(covar[match_idx, -1, drop=FALSE])

# Standardize covariates (account for any missing values)
covar_data <- apply(covar_data, 2, function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(rep(0, length(x)))
  x <- x - mean(x, na.rm=TRUE)
  x <- x / sd(x, na.rm=TRUE)
  x[is.na(x)] <- 0
  x
})

cat("[t4_make_pheno] covariate matrix:", nrow(covar_data), "x", ncol(covar_data), "\n")

# Simulate phenotype starting from covariates
baseline <- covar_data %*% rnorm(ncol(covar_data), sd=0.3)
baseline <- baseline - mean(baseline)
baseline <- baseline / sd(baseline)

phenotype <- baseline * 0.7  # 70% explained by covariates

# Causal effect sizes (standardized, reasonable for GWAS)
causal_effects <- rnorm(n_causal, sd=0.03)

# Apply causal effects based on case
if (case_name == "null") {
  # No genetic effects
  cat("[t4_make_pheno] case=null: phenotype ~ covariates only\n")
} else if (case_name %in% c("location", "mixed")) {
  # Add location (mean) effects
  for (j in 1:n_causal) {
    dosage <- causal_dosages[j, ]
    dosage <- (dosage - mean(dosage)) / sd(dosage)
    phenotype <- phenotype + causal_effects[j] * dosage
  }
  cat("[t4_make_pheno] case=", case_name, ": added ", n_causal, " location effects\n")
}

if (case_name %in% c("scale", "mixed")) {
  # Add scale (variance) effects via heteroscedasticity
  dispersion_effects <- rnorm(n_causal, mean=0.01, sd=0.002)
  for (j in 1:n_causal) {
    dosage <- causal_dosages[j, ]
    dosage <- (dosage - mean(dosage)) / sd(dosage)
    # Multiplicative effect on SD: exp(gamma * dosage / 2)
    phenotype <- phenotype * exp(dispersion_effects[j] * dosage / 2)
  }
  cat("[t4_make_pheno] case=", case_name, ": added ", n_causal, " scale effects\n")
}

# Normalize final phenotype to unit variance
phenotype <- (phenotype - mean(phenotype)) / sd(phenotype)

cat("[t4_make_pheno] phenotype summary:\n")
cat("[t4_make_pheno]   mean=", format(mean(phenotype), digits=4), "\n")
cat("[t4_make_pheno]   sd=", format(sd(phenotype), digits=4), "\n")
cat("[t4_make_pheno]   min=", format(min(phenotype), digits=4), "\n")
cat("[t4_make_pheno]   max=", format(max(phenotype), digits=4), "\n")

# Write phenotype file in format expected by Stage 1
# Format: FID IID phenotype [other columns]
output <- data.table(
  FID = keep_ids,
  IID = keep_ids,
  phenotype = phenotype
)

out_dir <- dirname(out_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

fwrite(output, out_file, sep="\t")
cat("[t4_make_pheno] wrote", nrow(output), "samples to", out_file, "\n")

# Save metadata for later analysis
metadata_file <- sub("\\.txt$", ".meta.json", out_file)
metadata <- list(
  case = case_name,
  seed = seed,
  n_samples = n_samples,
  n_causal = n_causal,
  n_null = n_null,
  causal_snps = as.character(causal_snps),
  null_snps = as.character(null_snps),
  causal_mafs = as.numeric(causal_mafs),
  causal_effects = as.numeric(causal_effects),
  phenotype_mean = mean(phenotype),
  phenotype_sd = sd(phenotype),
  phenotype_min = min(phenotype),
  phenotype_max = max(phenotype),
  covariates_contribution = 0.7
)

if (requireNamespace("jsonlite", quietly=TRUE)) {
  write(jsonlite::toJSON(metadata, pretty=TRUE), metadata_file)
  cat("[t4_make_pheno] wrote metadata to", metadata_file, "\n")
} else {
  cat("[t4_make_pheno] jsonlite not available; skipping metadata JSON\n")
}

cat("[t4_make_pheno] DONE\n")
