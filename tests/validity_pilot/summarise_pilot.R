#!/usr/bin/env Rscript
# Summarise pilot outputs with the dataset as the Monte Carlo unit.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: summarise_pilot.R <pilot_dir>")
dir <- normalizePath(args[1L])
ds <- read.csv(file.path(dir, "pilot_rows_datasets.csv"), stringsAsFactors = FALSE)
detail <- read.csv(file.path(dir, "pilot_rows.csv"), stringsAsFactors = FALSE)
if (!nrow(ds) || !nrow(detail)) stop("pilot output is empty")

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
sd_or_na <- function(x) if (sum(is.finite(x)) < 2L) NA_real_ else sd(x, na.rm = TRUE)
mcse <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2L) NA_real_ else sd(x) / sqrt(length(x)) }

keys <- unique(ds[c("cell", "path", "source", "W")])
sum_rows <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) {
  z <- ds[ds$cell == keys$cell[i] & ds$path == keys$path[i] &
          ds$source == keys$source[i] & ds$W == keys$W[i], , drop = FALSE]
  ok <- z$fit_ok
  sum_rows[[i]] <- data.frame(
    cell = keys$cell[i], path = keys$path[i], source = keys$source[i], W = keys$W[i],
    datasets = nrow(z), fit_fail = sum(!ok),
    null_type1_mean = mean_or_na(z$null_type1), null_type1_mcse = mcse(z$null_type1),
    null_any_mean = mean_or_na(z$null_any), null_any_mcse = mcse(z$null_any),
    causal_cover_projection_mean = mean_or_na(z$causal_cover_proj),
    causal_cover_projection_mcse = mcse(z$causal_cover_proj),
    causal_cover_scientific_mean = mean_or_na(z$causal_cover_sci),
    causal_cover_scientific_mcse = mcse(z$causal_cover_sci),
    stringsAsFactors = FALSE
  )
}
summary <- do.call(rbind, sum_rows)
write.csv(summary, file.path(dir, "pilot_summary.csv"), row.names = FALSE, na = "")

# Bias and calibration are computed over the causal variant for direct-effect
# cells. Each row below is still a dataset-level statistic, so the standard
# deviation and its MCSE do not count six taus or repeated SNP rows as new
# phenotype replicates.
detail$bias_projection <- detail$estimate - detail$target_projection
detail$bias_scientific <- detail$estimate - detail$target_scientific
direct <- detail[!grepl("^(cov_both_|null_)", detail$cell) & detail$variant == 1, , drop = FALSE]
bkeys <- unique(direct[c("cell", "path", "source", "W", "parameter")])
brows <- vector("list", nrow(bkeys))
for (i in seq_len(nrow(bkeys))) {
  z <- direct[direct$cell == bkeys$cell[i] & direct$path == bkeys$path[i] &
              direct$source == bkeys$source[i] & direct$W == bkeys$W[i] &
              direct$parameter == bkeys$parameter[i], , drop = FALSE]
  cut <- qt(0.975, df = z$param_ref_df)
  z_proj_cover <- abs(z$z_projection) <= cut
  z_sci_cover <- abs(z$z_scientific) <= cut
  brows[[i]] <- data.frame(
    cell = bkeys$cell[i], path = bkeys$path[i], source = bkeys$source[i],
    W = bkeys$W[i], parameter = bkeys$parameter[i], datasets = nrow(z),
    bias_projection = mean_or_na(z$bias_projection),
    bias_projection_mcse = mcse(z$bias_projection),
    empirical_sd_projection = sd_or_na(z$bias_projection),
    mean_reported_se = mean_or_na(z$se),
    se_over_empirical_sd = mean_or_na(z$se) / sd_or_na(z$bias_projection),
    bias_over_empirical_sd = abs(mean_or_na(z$bias_projection)) / sd_or_na(z$bias_projection),
    cover_projection = mean(z_proj_cover, na.rm = TRUE),
    cover_scientific = mean(z_sci_cover, na.rm = TRUE),
    bias_scientific = mean_or_na(z$bias_scientific),
    bias_scientific_mcse = mcse(z$bias_scientific),
    stringsAsFactors = FALSE
  )
}
bias <- do.call(rbind, brows)
write.csv(bias, file.path(dir, "pilot_bias.csv"), row.names = FALSE, na = "")

# Paired covariate comparison. The same dataset IDs and seeds are retained in
# Stage 1 for current and y_first, so deltas are paired diagnostics.
covds <- ds[grepl("^cov_both_", ds$cell), , drop = FALSE]
if (nrow(covds)) {
  pair_keys <- unique(covds[c("cell", "rep", "source", "W")])
  pairs <- vector("list", nrow(pair_keys))
  for (i in seq_len(nrow(pair_keys))) {
    z <- covds[covds$cell == pair_keys$cell[i] & covds$rep == pair_keys$rep[i] &
               covds$source == pair_keys$source[i] & covds$W == pair_keys$W[i], , drop = FALSE]
    if (nrow(z) != 2L) next
    cur <- z[z$path == "current", , drop = FALSE]
    yf <- z[z$path == "y_first", , drop = FALSE]
    if (!nrow(cur) || !nrow(yf)) next
    pairs[[i]] <- data.frame(cell = pair_keys$cell[i], rep = pair_keys$rep[i],
                             source = pair_keys$source[i], W = pair_keys$W[i],
                             current_null_type1 = cur$null_type1,
                             y_first_null_type1 = yf$null_type1,
                             delta_null_type1 = cur$null_type1 - yf$null_type1,
                             current_null_any = cur$null_any,
                             y_first_null_any = yf$null_any,
                             delta_null_any = cur$null_any - yf$null_any,
                             stringsAsFactors = FALSE)
  }
  pairs <- pairs[!vapply(pairs, is.null, logical(1L))]
  if (length(pairs)) write.csv(do.call(rbind, pairs), file.path(dir, "pilot_covariate_pairs.csv"), row.names = FALSE)
}

cat("wrote", nrow(summary), "dataset-arm summaries and", nrow(bias), "causal bias rows\n")
