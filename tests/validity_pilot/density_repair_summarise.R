#!/usr/bin/env Rscript
# Summarise the test-only fourth-order density repair prototype.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: density_repair_summarise.R <probe_dir>")
dir <- normalizePath(args[1L])
detail <- read.csv(file.path(dir, "pilot_rows.csv"), stringsAsFactors = FALSE)
ds <- read.csv(file.path(dir, "pilot_rows_datasets.csv"), stringsAsFactors = FALSE)
aux <- read.csv(file.path(dir, "density_repair_aux.csv"), stringsAsFactors = FALSE)
runtime <- read.csv(file.path(dir, "density_repair_runtime.csv"), stringsAsFactors = FALSE)
if (!nrow(detail) || !nrow(ds) || !nrow(aux) || !nrow(runtime)) stop("repair output is empty")

mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else sd(x) / sqrt(length(x))
}
mean_or_na <- function(x) if (all(!is.finite(x))) NA_real_ else mean(x, na.rm = TRUE)
sd_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else sd(x)
}

direct <- detail[detail$variant == 1L & detail$W %in% c("oracle", "estimated"), , drop = FALSE]
direct$arm <- sub("^repair_v(0|06|24)_", "", direct$cell)
direct$base_cell <- sub("_(A_production_kde|D_fourth_order_kde)$", "", direct$cell)
keys <- unique(direct[c("base_cell", "source", "W", "parameter")])
arm_rows <- list(); ia <- 0L
for (i in seq_len(nrow(keys))) {
  k <- keys[i, ]
  z <- direct[direct$base_cell == k$base_cell & direct$source == k$source &
              direct$W == k$W & direct$parameter == k$parameter, , drop = FALSE]
  for (arm in c("A_production_kde", "D_fourth_order_kde")) {
    za <- z[z$arm == arm, , drop = FALSE]
    if (!nrow(za)) next
    dsz <- ds[ds$cell == za$cell[1L] & ds$path == "current" &
              ds$source == k$source & ds$W == k$W, , drop = FALSE]
    attempted <- nrow(dsz)
    successful <- sum(dsz$fit_ok %in% TRUE)
    ia <- ia + 1L
    proj_bias <- za$estimate - za$target_projection
    sci_bias <- za$estimate - za$target_scientific
    cut <- qt(0.975, df = za$param_ref_df)
    finite_est <- is.finite(za$estimate)
    finite_se <- is.finite(za$se)
    finite_z_proj <- is.finite(za$z_projection) & is.finite(cut)
    finite_z_sci <- is.finite(za$z_scientific) & is.finite(cut)
    proj_cover_n <- sum(finite_z_proj & abs(za$z_projection) <= cut)
    sci_cover_n <- sum(finite_z_sci & abs(za$z_scientific) <= cut)
    ia_rows <- data.frame(
      base_cell = k$base_cell, source = k$source, W = k$W,
      parameter = k$parameter, arm = arm, attempted_datasets = attempted,
      successful_fits = successful, detail_rows = nrow(za),
      nonfinite_estimates = sum(!finite_est), nonfinite_se = sum(!finite_se),
      nonfinite_z_projection = sum(!finite_z_proj), nonfinite_z_scientific = sum(!finite_z_sci),
      bias_projection = mean_or_na(proj_bias), bias_projection_mcse = mcse(proj_bias),
      bias_scientific = mean_or_na(sci_bias), bias_scientific_mcse = mcse(sci_bias),
      empirical_sd_projection = sd_or_na(proj_bias), mean_reported_se = mean_or_na(za$se),
      se_over_empirical_sd = mean_or_na(za$se) / sd_or_na(proj_bias),
      bias_projection_finite_datasets = sum(is.finite(proj_bias)),
      bias_scientific_finite_datasets = sum(is.finite(sci_bias)),
      cover_projection = if (attempted) proj_cover_n / attempted else NA_real_,
      cover_scientific = if (attempted) sci_cover_n / attempted else NA_real_,
      cover_projection_finite = if (sum(finite_z_proj)) proj_cover_n / sum(finite_z_proj) else NA_real_,
      cover_scientific_finite = if (sum(finite_z_sci)) sci_cover_n / sum(finite_z_sci) else NA_real_,
      stringsAsFactors = FALSE)
    arm_rows[[ia]] <- ia_rows
  }
}
arm_summary <- do.call(rbind, arm_rows)
write.csv(arm_summary, file.path(dir, "density_repair_arm_summary.csv"), row.names = FALSE, na = "")

pair_rows <- list(); ip <- 0L
for (i in seq_len(nrow(keys))) {
  k <- keys[i, ]
  z <- direct[direct$base_cell == k$base_cell & direct$source == k$source &
              direct$W == k$W & direct$parameter == k$parameter, , drop = FALSE]
  lhs <- z[z$arm == "A_production_kde", c("rep", "estimate", "se"), drop = FALSE]
  rhs <- z[z$arm == "D_fourth_order_kde", c("rep", "estimate", "se"), drop = FALSE]
  names(lhs)[2:3] <- c("estimate_A", "se_A")
  names(rhs)[2:3] <- c("estimate_D", "se_D")
  zz <- merge(lhs, rhs, by = "rep", all = FALSE)
  if (!nrow(zz)) next
  d <- zz$estimate_D - zz$estimate_A
  dsz <- ds[ds$cell %in% unique(z$cell) & ds$path == "current" &
            ds$source == k$source & ds$W == k$W, , drop = FALSE]
  attempted_pair <- if (nrow(dsz)) max(table(dsz$cell)) else 0L
  ip <- ip + 1L
  pair_rows[[ip]] <- data.frame(
    base_cell = k$base_cell, source = k$source, W = k$W,
    parameter = k$parameter, comparison = "D_fourth_order - A_production",
    attempted_datasets = attempted_pair, paired_datasets = length(d),
    nonfinite_differences = sum(!is.finite(d)),
    mean_difference = mean_or_na(d), difference_mcse = mcse(d),
    empirical_sd_difference = sd_or_na(d), stringsAsFactors = FALSE)
}
pair_summary <- do.call(rbind, pair_rows)
write.csv(pair_summary, file.path(dir, "density_repair_pair_summary.csv"), row.names = FALSE, na = "")

aux$base_cell <- sub("_(A_production_kde|D_fourth_order_kde)$", "", aux$cell)
aux_keys <- unique(aux[c("base_cell", "arm")])
aux_rows <- vector("list", nrow(aux_keys))
for (i in seq_len(nrow(aux_keys))) {
  k <- aux_keys[i, ]
  z <- aux[aux$base_cell == k$base_cell & aux$arm == k$arm, , drop = FALSE]
  aux_rows[[i]] <- data.frame(
    base_cell = k$base_cell, arm = k$arm, datasets = nrow(z),
    mean_raw_nonpositive_tau = mean(z$raw_nonpositive_tau),
    max_raw_nonpositive_tau = max(z$raw_nonpositive_tau),
    nonpositive_dataset_rate = mean(z$raw_nonpositive_tau > 0),
    mean_floor_tau = mean(z$floor_tau), max_floor_tau = max(z$floor_tau),
    floor_dataset_rate = mean(z$floor_tau > 0),
    mean_raw_density_min = mean(z$raw_density_min),
    mean_raw_over_exact = mean(z$raw_over_exact_mean),
    raw_over_exact_mcse = mcse(z$raw_over_exact_mean),
    mean_q_delta_rms = mean(z$q_delta_rms), q_delta_rms_mcse = mcse(z$q_delta_rms),
    mean_bandwidth = mean(z$bandwidth), stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, aux_rows), file.path(dir, "density_repair_aux_summary.csv"),
          row.names = FALSE, na = "")
write.csv(runtime, file.path(dir, "density_repair_runtime_summary.csv"), row.names = FALSE, na = "")

cat("wrote", nrow(arm_summary), "arm rows,", nrow(pair_summary),
    "paired rows and", nrow(aux_keys), "density diagnostics\n")
