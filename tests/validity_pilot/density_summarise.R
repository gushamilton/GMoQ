#!/usr/bin/env Rscript
# Summarise paired density/RIF diagnosis after pilot_stage2.R.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: density_summarise.R <density_probe_dir>")
dir <- normalizePath(args[1L])
detail <- read.csv(file.path(dir, "pilot_rows.csv"), stringsAsFactors = FALSE)
aux <- read.csv(file.path(dir, "density_aux.csv"), stringsAsFactors = FALSE)
if (!nrow(detail) || !nrow(aux)) stop("density probe output is empty")

mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else sd(x) / sqrt(length(x))
}
mean_or_na <- function(x) if (all(!is.finite(x))) NA_real_ else mean(x, na.rm = TRUE)
sd_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else sd(x)
}

detail <- detail[detail$W == "oracle" & detail$variant == 1L, , drop = FALSE]
detail$arm <- sub("^density_v(0|06|24)_", "", detail$cell)
detail$base_cell <- sub("_(A_empirical_kde|B_empirical_exact_f|C_population_exact)$", "", detail$cell)

keys <- unique(detail[c("base_cell", "source", "parameter")])
rows <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) {
  k <- keys[i, ]
  z <- detail[detail$base_cell == k$base_cell & detail$source == k$source &
              detail$parameter == k$parameter, , drop = FALSE]
  arms <- c("A_empirical_kde", "B_empirical_exact_f", "C_population_exact")
  for (arm in arms) {
    za <- z[z$arm == arm, , drop = FALSE]
    if (!nrow(za)) next
    rows[[i]] <- rbind(rows[[i]], data.frame(
      base_cell = k$base_cell, source = k$source, parameter = k$parameter,
      comparison = arm, datasets = nrow(za),
      bias_projection = mean_or_na(za$estimate - za$target_projection),
      bias_projection_mcse = mcse(za$estimate - za$target_projection),
      bias_scientific = mean_or_na(za$estimate - za$target_scientific),
      bias_scientific_mcse = mcse(za$estimate - za$target_scientific),
      empirical_sd_projection = sd_or_na(za$estimate - za$target_projection),
      mean_reported_se = mean_or_na(za$se),
      stringsAsFactors = FALSE))
  }
}
arm_summary <- do.call(rbind, rows)
write.csv(arm_summary, file.path(dir, "density_arm_summary.csv"), row.names = FALSE, na = "")

# Paired differences use the same frozen dataset and replicate number.
pair_defs <- list(c("A_empirical_kde", "B_empirical_exact_f"),
                  c("B_empirical_exact_f", "C_population_exact"))
pair_rows <- list(); ir <- 0L
for (k in seq_len(nrow(keys))) {
  key <- keys[k, ]
  z <- detail[detail$base_cell == key$base_cell & detail$source == key$source &
              detail$parameter == key$parameter, , drop = FALSE]
  for (pd in pair_defs) {
    lhs <- z[z$arm == pd[1L], c("rep", "estimate"), drop = FALSE]
    rhs <- z[z$arm == pd[2L], c("rep", "estimate"), drop = FALSE]
    names(lhs)[2L] <- "lhs"
    names(rhs)[2L] <- "rhs"
    zz <- merge(lhs, rhs, by = "rep", all = FALSE)
    if (!nrow(zz)) next
    d <- zz$lhs - zz$rhs
    ir <- ir + 1L
    pair_rows[[ir]] <- data.frame(
      base_cell = key$base_cell, source = key$source,
      parameter = key$parameter, comparison = paste(pd, collapse = " - "),
      datasets = length(d), mean_difference = mean(d), difference_mcse = mcse(d),
      empirical_sd_difference = sd_or_na(d), stringsAsFactors = FALSE)
  }
}
pair_summary <- do.call(rbind, pair_rows)
write.csv(pair_summary, file.path(dir, "density_pair_summary.csv"), row.names = FALSE, na = "")

aux$arm <- as.character(aux$arm)
aux$base_cell <- sub("_(A_empirical_kde|B_empirical_exact_f|C_population_exact)$", "", aux$cell)
aux_keys <- unique(aux[c("base_cell", "arm")])
aux_rows <- vector("list", nrow(aux_keys))
for (i in seq_len(nrow(aux_keys))) {
  k <- aux_keys[i, ]
  z <- aux[aux$base_cell == k$base_cell & aux$arm == k$arm, , drop = FALSE]
  aux_rows[[i]] <- data.frame(
    base_cell = k$base_cell, arm = k$arm, datasets = nrow(z),
    fhat_over_exact_mean = mean_or_na(z$f_ratio_mean),
    fhat_over_exact_mcse = mcse(z$f_ratio_mean),
    q_delta_rms_mean = mean_or_na(z$q_delta_rms),
    q_delta_rms_mcse = mcse(z$q_delta_rms), stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, aux_rows), file.path(dir, "density_aux_summary.csv"),
          row.names = FALSE, na = "")

cat("wrote", nrow(arm_summary), "arm rows,", nrow(pair_summary),
    "paired rows and", nrow(aux_keys), "auxiliary rows\n")
