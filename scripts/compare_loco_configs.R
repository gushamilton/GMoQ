#!/usr/bin/env Rscript

## Three-way comparison of G-MoQ Stage 2 chromosome results under different
## LOCO configurations (with-LOCO / no-LOCO / two-pass blind LOCO).
##
## Reports, for each configuration:
##   * median SE relative to sex-stratified REGENIE (the precision reference)
##   * correlation of component effects with sex-stratified REGENIE, all
##     variants and on the |z|>=4 subset
##   * genomic inflation lambda for all three Stage 2 parameters
##
## Output is SNP-level aggregate only.
##
## Usage:
##   Rscript scripts/compare_loco_configs.R \
##     --run "with LOCO=/path/chr22.tsv.gz" \
##     --run "no LOCO=/path/chr22.tsv.gz" \
##     --run "two-pass=/path/chr22.tsv.gz" \
##     --regenie-female <file> --regenie-male <file> --chr 22 [--out <tsv>]

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (!is.na(i) && i < length(args)) args[[i + 1L]] else default
}
get_all <- function(flag) {
  i <- which(args == flag); args[i + 1L]
}
r_lib <- get_arg("--r-lib"); if (!is.null(r_lib)) .libPaths(c(r_lib, .libPaths()))
suppressPackageStartupMessages(library(data.table)); setDTthreads(4)

run_specs <- get_all("--run")
if (!length(run_specs)) stop("Need at least one --run 'label=path'.", call. = FALSE)
reg_f <- get_arg("--regenie-female"); reg_m <- get_arg("--regenie-male")
chr <- as.integer(get_arg("--chr", "22"))
out_file <- get_arg("--out")

lam <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  stats::median(stats::qchisq(p, 1, lower.tail = FALSE)) / stats::qchisq(0.5, 1)
}

read_regenie <- function(path, chrom) {
  d <- fread(path, showProgress = FALSE)
  d <- d[CHROM == chrom]
  d[, P := 10^(-LOG10P)]
  d[, .(snp_id = ID, ALLELE0, ALLELE1, BETA, SE, P, N)]
}
rf <- read_regenie(reg_f, chr)
rm_ <- read_regenie(reg_m, chr)

## Harmonise G-MoQ (effect_allele / other_allele) against REGENIE (BETA is per
## ALLELE1).  Flip the REGENIE sign when the alleles are swapped; drop anything
## that is not a clean match.
harmonise <- function(g, r) {
  m <- merge(g, r, by = "snp_id")
  same <- m$effect_allele == m$ALLELE1 & m$other_allele == m$ALLELE0
  swap <- m$effect_allele == m$ALLELE0 & m$other_allele == m$ALLELE1
  m <- m[same | swap]
  sgn <- ifelse(m$effect_allele == m$ALLELE1, 1, -1)
  m[, BETA := BETA * sgn]
  m
}

summarise_run <- function(label, path) {
  s <- fread(cmd = paste("gzip -dc", shQuote(path)), showProgress = FALSE)
  base <- data.table(
    run = label, n_snps = nrow(s),
    lambda_gamma = lam(s$gamma_p),
    lambda_beta_left = lam(s$beta_left_p),
    lambda_beta_right = lam(s$beta_right_p),
    frac_p_lt_0.05_left = mean(s$beta_left_p < 0.05, na.rm = TRUE),
    frac_p_lt_1e3_left = mean(s$beta_left_p < 1e-3, na.rm = TRUE)
  )
  for (cfg in list(list(col = "beta_left", reg = rf, tag = "female"),
                   list(col = "beta_right", reg = rm_, tag = "male"))) {
    g <- s[, .(snp_id, effect_allele, other_allele,
               b = get(cfg$col), se = get(paste0(cfg$col, "_se")))]
    m <- harmonise(g, cfg$reg)
    z <- m$BETA / m$SE
    strong <- is.finite(z) & abs(z) >= 4
    base[[paste0("se_ratio_", cfg$tag)]] <- stats::median(m$se / m$SE, na.rm = TRUE)
    base[[paste0("cor_all_", cfg$tag)]] <- stats::cor(m$b, m$BETA, use = "complete.obs")
    base[[paste0("cor_z4_", cfg$tag)]] <- if (sum(strong) > 2) {
      stats::cor(m$b[strong], m$BETA[strong], use = "complete.obs")
    } else NA_real_
    base[[paste0("n_matched_", cfg$tag)]] <- nrow(m)
    base[[paste0("n_z4_", cfg$tag)]] <- sum(strong)
  }
  base
}

res <- rbindlist(lapply(run_specs, function(spec) {
  label <- sub("=.*$", "", spec); path <- sub("^[^=]*=", "", spec)
  summarise_run(label, path)
}), fill = TRUE)

cols_se <- c("run", "se_ratio_female", "se_ratio_male",
             "cor_all_female", "cor_all_male", "cor_z4_female", "cor_z4_male")
cols_lam <- c("run", "n_snps", "lambda_gamma", "lambda_beta_left",
              "lambda_beta_right", "frac_p_lt_0.05_left", "frac_p_lt_1e3_left")
cat("\n=== Precision and accuracy vs sex-stratified REGENIE ===\n")
print(res[, ..cols_se])
cat("\n=== Genomic inflation ===\n")
print(res[, ..cols_lam])

if (!is.null(out_file)) { fwrite(res, out_file, sep = "\t"); message("wrote ", out_file) }
