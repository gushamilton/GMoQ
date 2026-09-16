## Render the Monte Carlo arm in the SAME format as the published
## FINDINGS_covsource_default.md tables, so the two can be read line by line.
WS <- "/user/work/fh6520/gmoq_covmixW"
tab <- read.delim(file.path(WS, "out", "summary_mix.tsv"), stringsAsFactors = FALSE)
tab <- tab[tab$est == "ols", ]          # the shipped estimator
cat("cells:", length(unique(tab$cell)), " reps/cell:", paste(unique(tab$nrep), collapse=","),
    " (Monte Carlo, OLS only)\n")

wide <- function(sub, valcol, msecol = NULL) {
  p <- sub[sub$cov_source == "pooled", ]; j <- sub[sub$cov_source == "jackknife", ]
  m <- merge(p, j, by = c("cell", "fit", "param"), suffixes = c(".p", ".j"))
  m[order(m$maf.p, m$fit, m$param), ]
}

cat("\n\n########## MC1. 95% COVERAGE and SE RATIO, pooled+OLS vs jackknife+OLS ##########\n")
cat("SE ratio = mean reported SE / empirical sd; < 1 is anti-conservative.\n")
cat("MC error in parentheses.  K is the column count of W.\n\n")
for (dg in c("null", "loc", "disp")) {
  s <- tab[tab$dgp == dg, ]; if (!nrow(s)) next
  cat(sprintf("\n--- dgp = %s ---\n", toupper(dg)))
  cat(sprintf("%-28s %-6s %-2s %-12s | %-17s %-17s | %-17s %-17s\n",
              "cell", "W", "K", "param", "cov95 pooled", "cov95 JACKKNIFE",
              "SEratio pooled", "SEratio JACKKNIFE"))
  m <- wide(s)
  for (i in seq_len(nrow(m))) with(m[i, ], cat(sprintf(
    "%-28s %-6s %-2d %-12s | %8.4f (%.4f) %8.4f (%.4f) | %8.4f (%.4f) %8.4f (%.4f)\n",
    cell, fit, K.p, param, cov95_t.p, cov95_t_mcse.p, cov95_t.j, cov95_t_mcse.j,
    se_ratio.p, se_ratio_mcse.p, se_ratio.j, se_ratio_mcse.j)))
}

cat("\n\n########## MC2. TYPE-I at nominal 0.05 and extrapolated to 5e-8 ##########\n")
cat("Meaningful ONLY for parameters whose true value is zero (|Z_true| < 0.75,\n")
cat("marked NULL).  5e-8 is extrapolated from the measured tail scale of the\n")
cat("Wald statistic, exactly as sim/covsource/summarise_cs.R does.\n\n")
for (dg in c("null", "loc", "disp")) {
  s <- tab[tab$dgp == dg, ]; if (!nrow(s)) next
  cat(sprintf("\n--- dgp = %s ---\n", toupper(dg)))
  cat(sprintf("%-28s %-6s %-2s %-12s %5s | %-17s %-17s | %11s %11s\n",
              "cell", "W", "K", "param", "NULL?", "t1@0.05 pooled", "t1@0.05 JACKKNIFE",
              "infl5e8 P", "infl5e8 JK"))
  m <- wide(s)
  for (i in seq_len(nrow(m))) with(m[i, ], cat(sprintf(
    "%-28s %-6s %-2d %-12s %5s | %8.4f (%.4f) %8.4f (%.4f) | %11.4g %11.4g\n",
    cell, fit, K.p, param, ifelse(abs(Z_true.p) < 0.75, "YES", ""),
    t1_05.p, t1_05_mcse.p, t1_05.j, t1_05_mcse.j,
    infl_5e8_from_sd.p, infl_5e8_from_sd.j)))
}

cat("\n\n########## MC3. PAIRED jackknife - pooled, and the estimate identity ##########\n")
pr <- read.delim(file.path(WS, "out", "paired_mix.tsv"), stringsAsFactors = FALSE)
pr <- pr[pr$est == "ols", ]
cat(sprintf("%-28s %-6s %-12s | %9s %11s | %9s %9s | %14s\n", "cell", "W", "param",
            "d_cov95", "d_se_ratio", "d_t1@.05", "+-paired", "sd_jk / sd_pool"))
pr <- pr[order(pr$cell, pr$fit, pr$param), ]
for (i in seq_len(nrow(pr))) with(pr[i, ], cat(sprintf(
  "%-28s %-6s %-12s | %9.4f %11.4f | %9.4f %9.4f | %14.8f\n",
  cell, fit, param, d_cov95_t, d_se_ratio, d_t1_05, d_t1_05_paired_mcse,
  emp_sd_jk_over_pooled)))
cat(sprintf("\nMAX |sd_jk/sd_pooled - 1| over all rows = %.3e  (OLS point estimates are\n",
            max(abs(pr$emp_sd_jk_over_pooled - 1))))
cat("bit-identical between covariance sources, as A is formed before the se branch)\n")

cat("\n\n########## MC4. jackknife noise cost ##########\n")
cat(sprintf("%-28s %-6s %-12s | %12s %12s\n", "cell", "W", "param",
            "cv(SE) pooled", "cv(SE) jk"))
m <- wide(tab)
for (i in seq_len(nrow(m))) with(m[i, ], cat(sprintf(
  "%-28s %-6s %-12s | %12.5f %12.5f\n", cell, fit, param, cv_se.p, cv_se.j)))
