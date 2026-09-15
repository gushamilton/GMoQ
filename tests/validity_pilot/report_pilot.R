#!/usr/bin/env Rscript
# Render a compact human-readable pilot report from the reproducible CSVs.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: report_pilot.R <pilot_dir>")
dir <- normalizePath(args[1L])
ds <- read.csv(file.path(dir, "pilot_rows_datasets.csv"), stringsAsFactors = FALSE)
bias <- read.csv(file.path(dir, "pilot_bias.csv"), stringsAsFactors = FALSE)
meta <- jsonlite::fromJSON(file.path(dir, "pilot_meta.json"), simplifyVector = FALSE)

fmt <- function(x) ifelse(is.finite(x), sprintf("%.4f", x), "NA")
fmt_mc <- function(x, s) ifelse(is.finite(x), paste0(sprintf("%.4f", x), " +/- ",
                                                    ifelse(is.finite(s), sprintf("%.4f", s), "NA")), "NA")
lines <- c(
  "G-MoQ compact validity pilot",
  paste0("source_commit: ", meta$source_commit),
  paste0("host: ", meta$host),
  paste0("N=", meta$N, ", reps=", meta$reps, ", taus=", length(meta$taus),
         ", blocks=", meta$blocks, ", threads=", meta$threads),
  paste0("MC unit: ", meta$mc_unit),
  "Interpretation: diagnostic evidence only; this report does not certify scientific validity.",
  "",
  "Dataset-level rates (means over independent phenotype datasets):",
  "cell | path | source | W | n | fit_fail | null_type1 +/- MCSE | null_any +/- MCSE | cover_projection +/- MCSE | cover_scientific +/- MCSE"
)
for (cell in unique(ds$cell)) for (path in unique(ds$path[ds$cell == cell]))
  for (source in unique(ds$source[ds$cell == cell & ds$path == path]))
    for (W in unique(ds$W[ds$cell == cell & ds$path == path & ds$source == source])) {
      z <- ds[ds$cell == cell & ds$path == path & ds$source == source & ds$W == W, , drop = FALSE]
      lines <- c(lines, paste(cell, path, source, W, nrow(z), sum(!z$fit_ok),
                              fmt_mc(mean(z$null_type1, na.rm = TRUE),
                                     sd(z$null_type1, na.rm = TRUE) / sqrt(sum(is.finite(z$null_type1)))),
                              fmt_mc(mean(z$null_any, na.rm = TRUE),
                                     sd(z$null_any, na.rm = TRUE) / sqrt(sum(is.finite(z$null_any)))),
                              fmt_mc(mean(z$causal_cover_proj, na.rm = TRUE),
                                     sd(z$causal_cover_proj, na.rm = TRUE) / sqrt(sum(is.finite(z$causal_cover_proj)))),
                              fmt_mc(mean(z$causal_cover_sci, na.rm = TRUE),
                                     sd(z$causal_cover_sci, na.rm = TRUE) / sqrt(sum(is.finite(z$causal_cover_sci)))), sep = " | "))
    }
lines <- c(lines, "", "Direct-effect causal bias and calibration (one causal SNP per dataset):",
           "cell | source | W | parameter | n | bias_projection +/- MCSE | emp_sd | mean_se | se/emp_sd | bias/emp_sd | cover_projection | cover_scientific")
for (i in seq_len(nrow(bias))) {
  z <- bias[i, ]
  lines <- c(lines, paste(z$cell, z$source, z$W, z$parameter, z$datasets,
                          fmt_mc(z$bias_projection, z$bias_projection_mcse), fmt(z$empirical_sd_projection),
                          fmt(z$mean_reported_se), fmt(z$se_over_empirical_sd),
                          fmt(z$bias_over_empirical_sd), fmt(z$cover_projection),
                          fmt(z$cover_scientific), sep = " | "))
}
pairs_snp1 <- file.path(dir, "pilot_covariate_snp1.csv")
if (file.exists(pairs_snp1)) {
  s1 <- read.csv(pairs_snp1, stringsAsFactors = FALSE)
  lines <- c(lines, "", "Covariate stress: correlated conditional-null SNP1 by parameter:",
             "cell | path | source | W | parameter | n | mean Z(scientific) +/- MCSE | rejection +/- MCSE | effect bias +/- MCSE | mean Q +/- MCSE | Q-p rejection +/- MCSE")
  for (i in seq_len(nrow(s1))) {
    z <- s1[i, ]
    lines <- c(lines, paste(z$cell, z$path, z$source, z$W, z$parameter, z$datasets,
                            fmt_mc(z$mean_z_scientific, z$mean_z_scientific_mcse),
                            fmt_mc(z$rejection_rate, z$rejection_rate_mcse),
                            fmt_mc(z$effect_bias, z$effect_bias_mcse),
                            fmt_mc(z$mean_Q, z$mean_Q_mcse),
                            fmt_mc(z$Q_p_rejection, z$Q_p_rejection_mcse), sep = " | "))
  }
}
pairs_file <- file.path(dir, "pilot_covariate_pairs.csv")
if (file.exists(pairs_file)) {
  pairs <- read.csv(pairs_file, stringsAsFactors = FALSE)
  lines <- c(lines, "", "Paired covariate stress deltas (current - y_first):",
             "cell | source | W | n | current_type1 | y_first_type1 | delta_type1 | current_any | y_first_any | delta_any")
  for (cell in unique(pairs$cell)) for (source in unique(pairs$source[pairs$cell == cell]))
    for (W in unique(pairs$W[pairs$cell == cell & pairs$source == source])) {
      z <- pairs[pairs$cell == cell & pairs$source == source & pairs$W == W, , drop = FALSE]
      lines <- c(lines, paste(cell, source, W, nrow(z), fmt(mean(z$current_null_type1)),
                              fmt(mean(z$y_first_null_type1)), fmt(mean(z$delta_null_type1)),
                              fmt(mean(z$current_null_any)), fmt(mean(z$y_first_null_any)),
                              fmt(mean(z$delta_null_any)), sep = " | "))
    }
}
writeLines(lines, file.path(dir, "pilot_report.txt"))
cat(paste(lines, collapse = "\n"), "\n", sep = "")
