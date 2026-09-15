#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# t4_make_pheno.R -- Build phenotypes on real chr22 BGEN dosages.
#
# This script:
#   1. Extracts REAL dosages from BGEN for requested variants and samples
#   2. Builds phenotypes with known effects on REAL dosages
#   3. Records ground truth (alleles, effects, realized MAFs)
#   4. Validates sample and dosage alignment
#
# Key: some samples in keep file may not be in BGEN - we use only those found
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
cat("[t4_make_pheno] keep:", keep_file, "\n")
cat("[t4_make_pheno] case:", case_name, "\n")

set.seed(seed)

# Load sample IDs from keep file (for filtering)
keep <- fread(keep_file, header=FALSE, select=1)
keep_ids <- keep[[1]]
cat("[t4_make_pheno] loaded", length(keep_ids), "samples from keep file\n")

# Load covariate file to get sample order
covar <- fread(covar_file, nrows=Inf)
covar_ids <- as.character(covar[[1]])
cat("[t4_make_pheno] found", length(covar_ids), "samples in covariate file\n")

# Load SNP list
snps <- fread(snp_file, header=FALSE)
snp_ids <- snps[[1]]
cat("[t4_make_pheno] loaded", nrow(snps), "SNPs\n")

n_samples_keep <- length(keep_ids)
n_causal <- 10
n_null <- 100

# Sample causal variants
causal_snps <- sample(snp_ids, min(n_causal, length(snp_ids)), replace=FALSE)
remaining_snps <- setdiff(snp_ids, causal_snps)
null_snps <- sample(remaining_snps, min(n_null, length(remaining_snps)), replace=FALSE)

n_causal <- length(causal_snps)
n_null <- length(null_snps)

cat("[t4_make_pheno] selected", n_causal, "causal variants\n")
cat("[t4_make_pheno] selected", n_null, "null variants\n")

# ============================================================================
# EXTRACT REAL DOSAGES FROM BGEN
# ============================================================================
cat("[t4_make_pheno] extracting REAL dosages from BGEN...\n")

# Find helper script
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value=TRUE)[1]))
if (length(script_dir) == 0 || nchar(script_dir) == 0) {
  script_dir <- getwd()
}
helper_script <- file.path(dirname(script_dir), "tier4", "extract_dosages.py")
if (!file.exists(helper_script)) {
  helper_script <- "extract_dosages.py"
}

if (!file.exists(helper_script)) {
  cat("ERROR: cannot find extract_dosages.py\n", file=stderr())
  quit(status=1)
}

cat("[t4_make_pheno] using helper script:", helper_script, "\n")

# Create temp files
temp_dir <- tempdir()
causal_snps_file <- file.path(temp_dir, "causal_snps.txt")
null_snps_file <- file.path(temp_dir, "null_snps.txt")
keep_ids_file <- file.path(temp_dir, "keep_ids.txt")
causal_dosage_file <- file.path(temp_dir, "causal_dosages.tsv")
null_dosage_file <- file.path(temp_dir, "null_dosages.tsv")
causal_meta_file <- file.path(temp_dir, "causal_dosages.meta.json")
null_meta_file <- file.path(temp_dir, "null_dosages.meta.json")

# Write IDs to files
writeLines(as.character(causal_snps), causal_snps_file)
writeLines(as.character(null_snps), null_snps_file)
writeLines(as.character(keep_ids), keep_ids_file)

cat("[t4_make_pheno] extracting causal variant dosages...\n")
cmd_causal <- paste(
  "python3", shQuote(helper_script),
  "--bgen", shQuote(bgen_file),
  "--sample", shQuote(sample_file),
  "--variant-ids", shQuote(causal_snps_file),
  "--sample-ids", shQuote(keep_ids_file),
  "--out", shQuote(causal_dosage_file),
  "--meta-out", shQuote(causal_meta_file)
)

if (system(cmd_causal) != 0) {
  cat("ERROR: failed to extract causal dosages\n", file=stderr())
  quit(status=1)
}

cat("[t4_make_pheno] extracting null variant dosages...\n")
cmd_null <- paste(
  "python3", shQuote(helper_script),
  "--bgen", shQuote(bgen_file),
  "--sample", shQuote(sample_file),
  "--variant-ids", shQuote(null_snps_file),
  "--sample-ids", shQuote(keep_ids_file),
  "--out", shQuote(null_dosage_file),
  "--meta-out", shQuote(null_meta_file)
)

if (system(cmd_null) != 0) {
  cat("ERROR: failed to extract null dosages\n", file=stderr())
  quit(status=1)
}

# Load dosages
cat("[t4_make_pheno] loading extracted dosages...\n")
causal_dosage_data <- fread(causal_dosage_file)
causal_var_id <- causal_dosage_data$variant_id
causal_dosages <- as.matrix(causal_dosage_data[, -1])
rownames(causal_dosages) <- causal_var_id

null_dosage_data <- fread(null_dosage_file)
null_var_id <- null_dosage_data$variant_id
null_dosages <- as.matrix(null_dosage_data[, -1])
rownames(null_dosages) <- null_var_id

cat("[t4_make_pheno] causal dosages shape:", nrow(causal_dosages), "x", ncol(causal_dosages), "\n")
cat("[t4_make_pheno] null dosages shape:", nrow(null_dosages), "x", ncol(null_dosages), "\n")

# ============================================================================
# VALIDATE SAMPLE ALIGNMENT
# ============================================================================
cat("[t4_make_pheno] validating sample alignment...\n")

# Get sample IDs from dosages (these are the ones that matched BGEN)
dosage_sample_ids <- as.character(colnames(causal_dosages))
n_matched <- length(dosage_sample_ids)

cat("[t4_make_pheno] successfully matched", n_matched, "out of", length(keep_ids), "requested samples to BGEN\n")
cat("[t4_make_pheno] dosage sample IDs (first 5):", paste(head(dosage_sample_ids, 5), collapse=", "), "\n")

# Match to covariate file: get indices of matched samples in covar file
match_idx <- match(dosage_sample_ids, covar_ids)
if (any(is.na(match_idx))) {
  cat("ERROR: some dosage samples not found in covariate file\n", file=stderr())
  quit(status=1)
}

cat("[t4_make_pheno] sample alignment confirmed\n")

# ============================================================================
# LOAD METADATA
# ============================================================================
cat("[t4_make_pheno] loading BGEN metadata...\n")

causal_meta <- jsonlite::fromJSON(causal_meta_file)
causal_effect_alleles <- causal_meta$effect_alleles
causal_alt_alleles <- sapply(causal_meta$allele_pairs, function(x) x[2])
causal_realized_mafs <- causal_meta$realized_mafs

cat("[t4_make_pheno] causal realized MAF range: [",
    format(min(causal_realized_mafs), digits=3), ", ",
    format(max(causal_realized_mafs), digits=3), "]\n")

# Check dosage distribution
dosage_quantiles <- quantile(as.numeric(causal_dosages), probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
cat("[t4_make_pheno] dosage quantiles:", paste(format(dosage_quantiles, digits=3), collapse=", "), "\n")

n_unique_dosages <- length(unique(as.numeric(causal_dosages)))
cat("[t4_make_pheno] unique dosage values:", n_unique_dosages, "(should be many for real imputed)\n")

# ============================================================================
# BUILD PHENOTYPES FROM REAL DOSAGES
# ============================================================================
cat("[t4_make_pheno] building phenotypes from REAL dosages...\n")

# Get covariate data for matched samples in correct order
covar_data <- as.matrix(covar[match_idx, -1, drop=FALSE])

# Standardize covariates
covar_data <- apply(covar_data, 2, function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(rep(0, length(x)))
  x <- x - mean(x, na.rm=TRUE)
  x <- x / sd(x, na.rm=TRUE)
  x[is.na(x)] <- 0
  x
})

cat("[t4_make_pheno] covariate matrix:", nrow(covar_data), "x", ncol(covar_data), "\n")

# Baseline from covariates
baseline <- covar_data %*% rnorm(ncol(covar_data), sd=0.3)
baseline <- baseline - mean(baseline)
baseline <- baseline / sd(baseline)

phenotype <- baseline * 0.7

# Generate causal effects
causal_effects <- rnorm(n_causal, sd=0.03)

# Apply effects based on case
if (case_name == "null") {
  cat("[t4_make_pheno] case=null: phenotype ~ covariates only (no genetic effects)\n")
} else if (case_name %in% c("location", "mixed")) {
  for (j in 1:n_causal) {
    dosage <- causal_dosages[j, ]
    dosage <- (dosage - mean(dosage)) / sd(dosage)
    phenotype <- phenotype + causal_effects[j] * dosage
  }
  cat("[t4_make_pheno] case=", case_name, ": added ", n_causal, " location effects using REAL dosages\n")
}

if (case_name %in% c("scale", "mixed")) {
  dispersion_effects <- rnorm(n_causal, mean=0.01, sd=0.002)
  for (j in 1:n_causal) {
    dosage <- causal_dosages[j, ]
    dosage <- (dosage - mean(dosage)) / sd(dosage)
    phenotype <- phenotype * exp(dispersion_effects[j] * dosage / 2)
  }
  cat("[t4_make_pheno] case=", case_name, ": added ", n_causal, " scale effects using REAL dosages\n")
}

# Normalize phenotype
phenotype <- (phenotype - mean(phenotype)) / sd(phenotype)

cat("[t4_make_pheno] phenotype summary:\n")
cat("[t4_make_pheno]   mean=", format(mean(phenotype), digits=4), "\n")
cat("[t4_make_pheno]   sd=", format(sd(phenotype), digits=4), "\n")
cat("[t4_make_pheno]   min=", format(min(phenotype), digits=4), "\n")
cat("[t4_make_pheno]   max=", format(max(phenotype), digits=4), "\n")

# ============================================================================
# WRITE OUTPUTS
# ============================================================================
cat("[t4_make_pheno] writing outputs...\n")

output <- data.table(
  FID = dosage_sample_ids,
  IID = dosage_sample_ids,
  phenotype = phenotype
)

out_dir <- dirname(out_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

fwrite(output, out_file, sep="\t")
cat("[t4_make_pheno] wrote", nrow(output), "samples to", out_file, "\n")

# Save metadata
metadata_file <- sub("\\.txt$", ".meta.json", out_file)
metadata <- list(
  case = case_name,
  seed = seed,
  n_samples_requested = n_samples_keep,
  n_samples_matched = n_matched,
  n_causal = n_causal,
  n_null = n_null,
  causal_snps = as.character(causal_var_id),
  causal_effect_alleles = causal_effect_alleles,
  causal_alt_alleles = causal_alt_alleles,
  causal_realized_mafs = causal_realized_mafs,
  causal_effects = as.numeric(causal_effects),
  null_snps = as.character(null_var_id),
  phenotype_mean = mean(phenotype),
  phenotype_sd = sd(phenotype),
  covariates_contribution = 0.7,
  dosage_source = "REAL BGEN (not simulated)",
  n_unique_dosage_values = n_unique_dosages
)

if (case_name %in% c("scale", "mixed")) {
  metadata$dispersion_effects <- as.numeric(dispersion_effects)
}

if (requireNamespace("jsonlite", quietly=TRUE)) {
  write(jsonlite::toJSON(metadata, pretty=TRUE), metadata_file)
  cat("[t4_make_pheno] wrote metadata to", metadata_file, "\n")
}

# Cleanup
unlink(c(causal_snps_file, null_snps_file, keep_ids_file,
         causal_dosage_file, null_dosage_file,
         causal_meta_file, null_meta_file))

cat("[t4_make_pheno] DONE\n")
