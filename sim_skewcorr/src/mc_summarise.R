## Summarise the Monte Carlo through the SHIPPED Stage 1 + Stage 2, against the
## exact population prediction.  Every rate carries its Monte Carlo SE.
WS <- "/user/work/fh6520/skewcorr"
.libPaths(file.path(WS, "rlib"))
stopifnot(identical(.libPaths()[1], file.path(WS, "rlib")))

ex <- read.delim(file.path(WS, "out", "mc_exact.tsv"), stringsAsFactors = FALSE)
files <- sort(list.files(file.path(WS, "s2"), pattern = "^mc_cell[0-9]+\\.rds$",
                         full.names = TRUE))
cat("cells found:", length(files), "\n\n")

mcse <- function(p, n) sqrt(p * (1 - p) / n)
rows <- list()

for (f in files) {
  o <- readRDS(f)
  ci <- o$cell_idx
  for (est in c("ols", "gls")) {
    r <- o$res[[est]]
    n <- length(r$bs2)
    e <- ex[ex$cell == ci & ex$est == est, ]
    zu <- r$bs2 / r$ss2
    zc <- r$cb  / r$cs
    ru05 <- mean(abs(zu) > qnorm(0.975));  rc05 <- mean(abs(zc) > qnorm(0.975))
    ru3  <- mean(abs(zu) > qnorm(1 - 5e-4)); rc3 <- mean(abs(zc) > qnorm(1 - 5e-4))
    rows[[length(rows) + 1L]] <- data.frame(
      cell = ci, name = o$cell$name, est = est, n = n,
      maf = o$cell$maf, shape = o$cell$shape,
      delta = o$cell$delta, gamma = o$cell$gamma,
      z_mu_mc  = mean(r$bmu) / sd(r$bmu),
      tgt_s2   = e$tgt_s2,
      ## bias of the two estimators about the generating truth
      bias_unc = mean(r$bs2) - e$tgt_s2,
      bias_cor = mean(r$cb)  - e$tgt_s2,
      zbias_unc = (mean(r$bs2) - e$tgt_s2) / mean(r$ss2),
      zbias_cor = (mean(r$cb)  - e$tgt_s2) / mean(r$cs),
      ## the estimate the exact engine predicts the MC should reproduce
      pop_s2_ex = e$pop_s2, e_cor_ex = e$e_cor,
      mean_unc = mean(r$bs2), mean_cor = mean(r$cb),
      ## SE calibration
      sesd_unc = mean(r$ss2) / sd(r$bs2),
      sesd_cor = mean(r$cs)  / sd(r$cb),
      se_ratio_mc = mean(r$cs) / mean(r$ss2),
      se_ratio_ex = e$sd_ratio_cor,
      ## rejection rates with MC error, against the exact prediction
      r05_unc = ru05, r05_unc_se = mcse(ru05, n), r05_unc_ex = e$unc_05,
      r05_cor = rc05, r05_cor_se = mcse(rc05, n), r05_cor_ex = e$cor_05,
      r1e3_unc = ru3, r1e3_unc_se = mcse(ru3, n), r1e3_unc_ex = e$unc_1e3,
      r1e3_cor = rc3, r1e3_cor_se = mcse(rc3, n), r1e3_cor_ex = e$cor_1e3,
      stringsAsFactors = FALSE)
  }
}
tb <- do.call(rbind, rows)
write.table(tb, file.path(WS, "out", "mc_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

fmt <- function(p, s) sprintf("%.4f+-%.4f", p, s)
dev <- function(p, s, e) if (s > 0) sprintf("%+.1f", (p - e) / s) else "  -"

cat("=================================================================================================================\n")
cat("MC through the SHIPPED Stage 1 + Stage 2, vs the EXACT prediction.  N =", 30000, " reps per cell as shown.\n")
cat("dev is (MC - exact) in MC standard errors.\n")
cat("=================================================================================================================\n")
cat(sprintf("%-22s %4s %6s | %20s %9s %4s | %20s %9s %4s\n",
            "cell", "est", "n", "unc@0.05 +- MCSE", "exact", "dev", "cor@0.05 +- MCSE", "exact", "dev"))
for (i in seq_len(nrow(tb))) {
  x <- tb[i, ]
  cat(sprintf("%-22s %4s %6d | %20s %9.4f %4s | %20s %9.4f %4s\n",
              x$name, x$est, x$n,
              fmt(x$r05_unc, x$r05_unc_se), x$r05_unc_ex, dev(x$r05_unc, x$r05_unc_se, x$r05_unc_ex),
              fmt(x$r05_cor, x$r05_cor_se), x$r05_cor_ex, dev(x$r05_cor, x$r05_cor_se, x$r05_cor_ex)))
}

cat("\n\n=================================================================================================================\n")
cat("Deep tail: nominal 1e-3, twenty times deeper than 0.05\n")
cat("=================================================================================================================\n")
cat(sprintf("%-22s %4s | %22s %10s %5s | %22s %10s %5s\n",
            "cell", "est", "unc@1e-3 +- MCSE", "exact", "dev", "cor@1e-3 +- MCSE", "exact", "dev"))
for (i in seq_len(nrow(tb))) {
  x <- tb[i, ]
  cat(sprintf("%-22s %4s | %22s %10.5f %5s | %22s %10.5f %5s\n",
              x$name, x$est,
              sprintf("%.5f+-%.5f", x$r1e3_unc, x$r1e3_unc_se), x$r1e3_unc_ex,
              dev(x$r1e3_unc, x$r1e3_unc_se, x$r1e3_unc_ex),
              sprintf("%.5f+-%.5f", x$r1e3_cor, x$r1e3_cor_se), x$r1e3_cor_ex,
              dev(x$r1e3_cor, x$r1e3_cor_se, x$r1e3_cor_ex)))
}

cat("\n\n=================================================================================================================\n")
cat("Bias in beta_sigma2 about the generating truth, in units of its own reported SE, and SE calibration\n")
cat("(mean_se / empirical_sd; 1.00 is correct)\n")
cat("=================================================================================================================\n")
cat(sprintf("%-22s %4s %7s | %10s %10s | %9s %9s | %9s %9s\n",
            "cell", "est", "Z_mu", "zbias_unc", "zbias_cor",
            "sesd_unc", "sesd_cor", "seratio", "ex"))
for (i in seq_len(nrow(tb))) {
  x <- tb[i, ]
  cat(sprintf("%-22s %4s %7.1f | %10.3f %10.3f | %9.4f %9.4f | %9.4f %9.4f\n",
              x$name, x$est, x$z_mu_mc, x$zbias_unc, x$zbias_cor,
              x$sesd_unc, x$sesd_cor, x$se_ratio_mc, x$se_ratio_ex))
}
cat("\nwrote", file.path(WS, "out", "mc_summary.tsv"), "\n")
