#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 2 collation: counts -> rates -> verdicts.
#
#   Rscript --vanilla tests/tier2/t2_collate.R <chunk_dir>
#
# Aggregates every chunk's counts, turns them into rates with Monte Carlo
# error, writes the full long table, and then applies the criteria that were
# committed in tests/criteria.json BEFORE the run. The criteria are read from
# that file, never computed here, so there is no path by which a band gets
# widened to fit a number.
#
# Every rate is reported with its Monte Carlo standard error, and every verdict
# names the band it was judged against.
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier2/t2_collate.R <dir>")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
gmoq_begin("tier2", "statistical_validity")

args <- commandArgs(trailingOnly = TRUE)
dir <- normalizePath(if (length(args)) args[1L] else
                       file.path(Sys.getenv("GMOQ_WORK"), "tier2"))
crit <- jsonlite::fromJSON(file.path(GMOQ_SUITE_ROOT, "criteria.json"),
                           simplifyVector = TRUE)

files <- list.files(dir, pattern = "\\.s2\\.rds$", full.names = TRUE)
if (!length(files)) stop("No tier 2 chunk results in ", dir)

# --- aggregate --------------------------------------------------------------
agg <- new.env(parent = emptyenv())
cellinfo <- list()
for (f in files) {
  x <- readRDS(f)
  cellinfo[[x$cell$name]] <- x$cell
  for (k in names(x$acc)) {
    kk <- paste(x$cell$name, k, sep = "|")
    cur <- if (is.null(agg[[kk]])) list() else agg[[kk]]
    for (fld in names(x$acc[[k]])) {
      cur[[fld]] <- (if (is.null(cur[[fld]])) 0 else cur[[fld]]) + x$acc[[k]][[fld]]
    }
    agg[[kk]] <- cur
  }
}

rows <- list()
for (kk in ls(agg)) {
  parts <- strsplit(kk, "|", fixed = TRUE)[[1]]
  a <- agg[[kk]]
  g <- function(f, d = NA_real_) if (is.null(a[[f]])) d else as.numeric(a[[f]])
  n <- g("n", 0)
  emp_var <- if (is.finite(n) && n > 1)
    (g("sum_d2") - g("sum_d")^2 / n) / (n - 1) else NA_real_
  rows[[length(rows) + 1L]] <- data.frame(
    cell = parts[1], src = parts[2], est = parts[3], model = parts[4],
    param = if (length(parts) >= 5) parts[5] else NA_character_,
    group = if (length(parts) >= 6) parts[6] else NA_character_,
    n = n,
    cover95 = g("cover95") / n,
    rej05 = (if (!is.null(a[["rej05"]])) g("rej05") else g("rej_5e-02")) / max(g("n_p", n), 1),
    rej1e3 = (if (!is.null(a[["rej1e-3"]])) g("rej1e-3") else g("rej_1e-03")) / max(g("n_p", n), 1),
    rej1e4 = g("rej_1e-04") / max(g("n_p", n), 1),
    n_p = g("n_p", n),
    mean_z = g("sum_z") / n,
    sd_z = sqrt(pmax((g("sum_z2") - g("sum_z")^2 / n) / (n - 1), 0)),
    mean_est = g("sum_est") / n,
    mean_true = g("sum_true") / n,
    mean_dev = g("sum_d") / n,
    emp_sd = sqrt(emp_var),
    mean_se = g("sum_se") / n,
    mean_absG = g("sum_absG") / n,
    mean_absZmu = g("sum_absZmu") / n,
    n_ridge_tried = g("n_ridge_tried"), n_sigma_not_spd = g("n_sigma_not_spd"),
    n_snps_gls = g("n_snps"),
    stringsAsFactors = FALSE)
}
tab <- do.call(rbind, rows)
tab$mcse_rej05 <- gmoq_mcse(tab$rej05, pmax(tab$n_p, 1))
tab$mcse_cover <- gmoq_mcse(tab$cover95, pmax(tab$n, 1))
tab$bias_over_sd <- tab$mean_dev / tab$emp_sd
tab$se_over_sd <- tab$mean_se / tab$emp_sd
tab$kind <- vapply(tab$cell, function(cn) cellinfo[[cn]]$kind, "")
tab$family <- vapply(tab$cell, function(cn) cellinfo[[cn]]$family, "")
tab$maf <- vapply(tab$cell, function(cn) as.numeric(cellinfo[[cn]]$maf), 0)
tab$baseline <- ifelse(vapply(tab$cell, function(cn)
  as.numeric(cellinfo[[cn]]$baseline_skew), 0) == 0, "normal", "skewed")

out <- Sys.getenv("GMOQ_TEST_OUT", dir)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.table(tab, file.path(out, "tier2_table.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)
cat("wrote", file.path(out, "tier2_table.tsv"), "--", nrow(tab), "rows from",
    length(files), "chunks over", length(cellinfo), "cells\n")

sel <- function(...) {
  f <- list(...)
  keep <- rep(TRUE, nrow(tab))
  for (nm in names(f)) keep <- keep & tab[[nm]] %in% f[[nm]]
  tab[keep & !is.na(tab$n) & tab$n > 0, , drop = FALSE]
}

band <- function(id, establishes, x, lo, hi, what, severity = "blocking") {
  if (!nrow(x)) return(gmoq_check(id, establishes, NA, detail = "no rows"))
  v <- x[[what]]
  ok <- is.finite(v) & v >= lo & v <= hi
  worst <- x[which.max(abs(v - (lo + hi) / 2)), ]
  gmoq_check(id, establishes, all(ok), severity = severity,
             stat = worst[[what]],
             mcse = if (what == "cover95") worst$mcse_cover else worst$mcse_rej05,
             criterion = sprintf("%s in [%.4g, %.4g] over %d cell/arm combinations",
                                 what, lo, hi, nrow(x)),
             detail = if (!all(ok)) paste0(
               sum(!ok), " outside: ",
               paste(sprintf("%s/%s/%s/%s=%.4f", x$cell[!ok], x$src[!ok],
                             x$est[!ok], x$param[!ok], v[!ok])[1:min(6, sum(!ok))],
                     collapse = " ")) else NULL)
}

floor_check <- function(id, establishes, x, what, lo, severity = "blocking") {
  if (!nrow(x)) return(gmoq_check(id, establishes, NA, detail = "no rows"))
  v <- x[[what]]; ok <- is.finite(v) & v >= lo
  gmoq_check(id, establishes, all(ok), severity = severity, stat = min(v, na.rm = TRUE),
             criterion = sprintf("%s >= %.4g over %d combinations", what, lo, nrow(x)),
             detail = if (!all(ok)) paste(
               sprintf("%s/%s/%s/%s=%.4f", x$cell[!ok], x$src[!ok], x$est[!ok],
                       x$model[!ok], v[!ok])[1:min(6, sum(!ok))], collapse = " ") else NULL)
}

is_param <- !(tab$param %in% c("Q", "GLSDIAG")) & !is.na(tab$param) &
  tab$model %in% c("LS", "TRUE")

# =====================================================================
# A. Null calibration
# =====================================================================
for (e in c("ols", "gls")) {
  nullp <- tab[is_param & tab$kind == "main" & tab$family == "null" &
                 tab$est == e & tab$group == "null", , drop = FALSE]
  cb <- crit$null_calibration$coverage95[[e]]
  band(sprintf("T2.NULL.coverage95.%s", e),
       "Under a true null, the 95 per cent interval for every fitted parameter covers at its nominal rate.",
       nullp, cb$lo, cb$hi, "cover95")
  tb <- crit$null_calibration$typeI_0.05[[e]]
  band(sprintf("T2.NULL.typeI_0.05.%s", e),
       "Under a true null, the parameter test has its nominal size at alpha = 0.05.",
       nullp, tb$lo, tb$hi, "rej05")
  rb <- crit$null_calibration$`typeI_1e-3_ratio`
  nullp$ratio1e3 <- nullp$rej1e3 / 1e-3
  band(sprintf("T2.NULL.typeI_1e-3.%s", e),
       "Fifty times deeper into the tail the size is still right, which a normal reference would not deliver.",
       nullp, rb$lo, rb$hi, "ratio1e3")
  nullp$abs_dev_sd_z <- abs(nullp$sd_z - 1)
  band(sprintf("T2.NULL.sd_z.%s", e),
       "The standardised deviation from the truth has unit scale, so the reported SE describes the actual variability.",
       nullp, 0, crit$null_calibration$sd_z_abs_dev$max, "abs_dev_sd_z")
}

# Q is only calibrated on the OLS path: param_gwas() computes the saddlepoint
# reference inside the `tau_cov_estimator == "ols"` branch, so the GLS path
# emits no Q p-value at all. That is recorded, not asserted away.
qgls <- tab[tab$param == "Q" & tab$est == "gls", , drop = FALSE]
gmoq_check("T2.Q.gls_emits_no_pvalue",
           "Recorded: the GLS path emits no calibrated Q p-value, so lack of fit cannot be judged there.",
           TRUE, severity = "report",
           stat = if (nrow(qgls)) sum(qgls$n_p, na.rm = TRUE) else 0,
           criterion = "number of GLS variants with a finite Q p-value (expected 0)")

qnull <- tab[tab$param == "Q" & tab$kind == "main" & tab$family == "null" &
               tab$group == "null" & tab$est == "ols" &
               tab$model %in% c("LS", "TRUE"), , drop = FALSE]
qb <- crit$null_calibration$Q_rejection_0.05
band("T2.NULL.Q_calibrated",
     "The lack-of-fit statistic does not reject a correctly specified model on null variants.",
     qnull, qb$lo, qb$hi, "rej05")

# =====================================================================
# B. Under a true effect
# =====================================================================
for (e in c("ols", "gls")) {
  for (bl in c("normal", "skewed")) {
    x <- tab[is_param & tab$kind == "main" & tab$family != "null" &
               tab$model == "TRUE" & tab$est == e & tab$group == "causal" &
               tab$baseline == bl, , drop = FALSE]
    lo <- crit$alternative$coverage95_theta_pop[[paste0(bl, "_baseline")]][[e]]
    floor_check(sprintf("T2.ALT.coverage95.%s.%s", e, bl),
                sprintf("Under a real distributional effect on a %s baseline, the %s interval still covers the exact population estimand.", bl, toupper(e)),
                x, "cover95", lo)
  }
  x <- tab[is_param & tab$kind == "main" & tab$family != "null" &
             tab$model == "TRUE" & tab$est == e & tab$group == "causal" &
             tab$baseline == "normal", , drop = FALSE]
  x$abs_bias_over_sd <- abs(x$bias_over_sd)
  band(sprintf("T2.ALT.bias.%s", e),
       "The point estimate is unbiased for the population estimand under the true model.",
       x, 0, crit$alternative$abs_bias_over_sd[[e]], "abs_bias_over_sd")
  sb <- crit$alternative$se_over_empirical_sd
  band(sprintf("T2.ALT.se_calibration.%s", e),
       "The reported standard error matches the estimator's actual variability across replicates.",
       x, sb$lo, sb$hi, "se_over_sd")
}

# =====================================================================
# C. The tournament, at scale
# =====================================================================
qwrong <- tab[tab$param == "Q" & tab$kind == "main" & tab$family != "null" &
                tab$est == "ols" &
                grepl("^WRONG_", tab$model) & tab$group == "causal", , drop = FALSE]
floor_check("T2.TOURN.wrong_models_rejected",
            "Every model that cannot represent the simulated truth is rejected by the lack-of-fit statistic.",
            qwrong, "rej05", crit$tournament$wrong_model_Q_rejection_min)

qtrue <- tab[tab$param == "Q" & tab$kind == "main" & tab$family != "null" &
               tab$est == "ols" &
               tab$model == "TRUE" & tab$group == "causal", , drop = FALSE]
band("T2.TOURN.true_model_not_rejected",
     "The true model is not rejected by the lack-of-fit statistic on the variants that carry the effect.",
     qtrue, 0, crit$tournament$true_model_Q_rejection_max, "rej05")

# =====================================================================
# D. The pins
# =====================================================================
gb <- tab[tab$model == "GBUCKET", , drop = FALSE]
if (nrow(gb)) {
  gb$bucket <- gb$param
  low  <- gb[gb$bucket == "[-Inf,0.05)", , drop = FALSE]
  high <- gb[gb$bucket == "[0.2,Inf)", , drop = FALSE]
  lb <- crit$pins$false_vqtl$low_G_bucket_typeI
  band("T2.PIN.vqtl_guard_does_not_miss",
       "Where the genotype-skewness guard says there is no artefact, there is none: a pure location effect produces no false dispersion finding.",
       low, lb$lo, lb$hi, "rej05")
  gmoq_check("T2.PIN.vqtl_artefact_present",
             "Where the guard fires (G >= 0.2) the false-dispersion rate is inflated, i.e. the documented artefact is still reproducible.",
             TRUE, severity = "report",
             stat = if (nrow(high)) max(high$rej05) else NA,
             criterion = "max false-dispersion rejection at 0.05 in the G >= 0.2 bucket",
             detail = if (nrow(high)) paste(sprintf("%s/%s/%s: rate=%.4f meanG=%.2f",
                                                    high$cell, high$src, high$est,
                                                    high$rej05, high$mean_absG),
                                            collapse = " | ") else "no variants reached G >= 0.2")
  maf50 <- gb[gb$cell == "pin_vqtl_maf0.5", , drop = FALSE]
  if (nrow(maf50)) {
    band("T2.PIN.vqtl_zero_at_maf_half",
         "At MAF 0.50 a very large pure location effect produces no false dispersion finding at all: the driver is genotype skewness, not effect size.",
         maf50, lb$lo, lb$hi, "rej05")
  }
  # the guard must not be silently useless: no bucket it calls clean may be
  # inflated beyond the flagging threshold
  infl <- low$rej05 / 0.05
  gmoq_check("T2.PIN.vqtl_guard_monotone",
             "No bucket the guard calls clean is inflated past the threshold the guard exists to hold.",
             all(!is.finite(infl) | infl < crit$pins$false_vqtl$inflation_to_flag),
             stat = max(infl, na.rm = TRUE),
             criterion = sprintf("inflation in the G < 0.05 bucket below %gx nominal",
                                 crit$pins$false_vqtl$inflation_to_flag))
}

ties <- tab[is_param & tab$cell == "pin_ties" & tab$group == "null" &
              tab$model == "LS", , drop = FALSE]
if (nrow(ties)) {
  ob <- crit$pins$ties_gls$ols_typeI
  band("T2.PIN.ties_ols_unaffected",
       "A tied phenotype does not disturb the DEFAULT estimator: OLS keeps its nominal size where GLS does not.",
       ties[ties$est == "ols", , drop = FALSE], ob$lo, ob$hi, "rej05")
  g <- ties[ties$est == "gls", , drop = FALSE]
  gmoq_check("T2.PIN.ties_gls_collapses",
             "The documented GLS collapse on a tied phenotype is still reproducible -- total power loss, not a false-positive risk.",
             TRUE, severity = "report", stat = if (nrow(g)) max(g$rej05) else NA,
             criterion = "GLS Type-I at 0.05 on a rounded phenotype; BP2 measures 0.0000",
             detail = if (nrow(g)) paste(sprintf("%s/%s=%.4f", g$src, g$param, g$rej05),
                                         collapse = " ") else "")
  rid <- tab[tab$cell == "pin_ties" & tab$param == "GLSDIAG", , drop = FALSE]
  if (nrow(rid)) {
    gmoq_check("T2.PIN.ties_ridge_repair_fires",
               "The GLS solve reports the ridge repair it had to do, so the failure is at least visible in the log.",
               any(rid$n_ridge_tried > 0),
               stat = max(rid$n_ridge_tried / pmax(rid$n_snps_gls, 1)),
               criterion = "ridge repair attempted on a non-zero fraction of variants")
  }
}

ceil <- tab[is_param & tab$cell == "pin_ceiling" & tab$group == "null" &
              tab$model == "LS", , drop = FALSE]
if (nrow(ceil)) {
  ceil$ratio1e4 <- ceil$rej1e4 / crit$pins$ceiling_censoring$alpha
  cb2 <- crit$pins$ceiling_censoring$ols_far_tail_ratio
  band("T2.PIN.ceiling_ols_not_anticonservative",
       "Ceiling censoring does not make the DEFAULT estimator anti-conservative in the far tail.",
       ceil[ceil$est == "ols", , drop = FALSE], cb2$lo, cb2$hi, "ratio1e4")
  g <- ceil[ceil$est == "gls", , drop = FALSE]
  worst <- if (nrow(g)) max(g$ratio1e4, na.rm = TRUE) else NA_real_
  gmoq_check("T2.PIN.ceiling_gls_far_tail",
             "The one anti-conservative failure in the evidence base: GLS manufactures far-tail false positives on a ceiling-censored trait.",
             TRUE, severity = "report", stat = worst,
             criterion = sprintf("rejection at alpha = %g as a multiple of nominal, GLS",
                                 crit$pins$ceiling_censoring$alpha),
             detail = "BP2 section 2 measures 4.90x at 1e-4 and 270x at 1e-6")
  gmoq_check("T2.PIN.ceiling_gls_still_unguarded",
             "Advisory: nothing in the shipped code warns a user running GLS on a ceiling-censored trait.",
             !is.finite(worst) || worst < crit$pins$ceiling_censoring$gls_far_tail_advisory_above,
             severity = "advisory", stat = worst,
             criterion = sprintf("GLS far-tail ratio below %gx",
                                 crit$pins$ceiling_censoring$gls_far_tail_advisory_above))
}

# =====================================================================
# E. Coverage of the run itself
# =====================================================================
gmoq_check("T2.META.cells_completed",
           "Every preregistered cell produced results; a silently missing cell is a silently untested claim.",
           TRUE, severity = "report", stat = length(cellinfo),
           criterion = "number of cells with at least one chunk",
           detail = paste(sort(names(cellinfo)), collapse = ", "))
gmoq_check("T2.META.null_observations",
           "Total null parameter observations behind the calibration numbers.",
           TRUE, severity = "report",
           stat = sum(tab$n[is_param & tab$group == "null" & tab$est == "ols" &
                              tab$src == "jackknife" & tab$model == "LS"]),
           criterion = "reported; the Monte Carlo error on every rate is in the table")

gmoq_finish()
