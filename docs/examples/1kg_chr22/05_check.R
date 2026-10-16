# Judge Stage 2 against the simulated truth. Writes recovery.tsv, null_calibration.tsv and qq.png
# and prints PASS/FAIL against these thresholds:
#   recovery     every causal SNP: |est - truth| / SE < 3, and its target parameter p < 1e-4
#                (judged on null/shift/variance and on mixture_hybrid, the default; the
#                 density = "model" fit is reported but known to attenuate gamma here)
#   calibration  per parameter over null SNPs: share with p < 0.05 within 0.05 +/- 3 binomial SDs
#                (about +/-0.015 for 2,000 SNPs, +/-0.027 for 600). lambda_GC is reported
#                but not judged: a median over a few hundred SNPs is too noisy to gate on.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
work <- if (length(args)) args[1] else "work"
truth <- fread(file.path(work, "truth.tsv"))
ld <- fread(file.path(work, "ld_to_causal.tsv"))
params <- list(null = c("beta_mu", "beta_sigma2"), shift = c("beta_mu", "beta_sigma2"),
               variance = c("beta_mu", "beta_sigma2"), mixture = c("gamma", "beta_1", "beta_2"),
               mixture_hybrid = c("gamma", "beta_1", "beta_2"))

rec <- list(); cal <- list(); qq <- list()
for (ph in names(params)) for (src in c("pooled", "jackknife")) {
  r <- fread(file.path(work, sprintf("s2_%s_%s.tsv.gz", ph, src)))
  m <- sub("_hybrid", "", ph)
  causal <- truth[model == m]
  # Null SNPs: r^2 < 0.01 with every causal SNP of the simulated model (none for "null").
  null_i <- if (m == "null") rep(TRUE, nrow(r)) else ld[[m]][match(r$snp_id, ld$snp_id)] < 0.01
  for (pa in params[[ph]]) {
    est <- r[[pa]]; se <- r[[paste0(pa, "_se")]]; p <- r[[paste0(pa, "_p")]]
    # Recovery: every causal SNP of this model, on every parameter (truth 0 off-target).
    if (nrow(causal)) {
      i <- match(causal$snp_id, r$snp_id)
      tr <- ifelse(causal$param == pa, causal$truth, 0)
      rec[[length(rec) + 1]] <- data.table(model = ph, cov_source = src, snp_id = causal$snp_id,
        maf = round(causal$maf, 3), param = pa, truth = round(tr, 3), est = round(est[i], 3),
        se = round(se[i], 3), z_vs_truth = round((est[i] - tr) / se[i], 2), p = signif(p[i], 3),
        target = causal$param == pa)
    }
    chi <- qchisq(p[null_i], 1, lower.tail = FALSE)
    cal[[length(cal) + 1]] <- data.table(model = ph, cov_source = src, param = pa,
      n_null = sum(null_i), lambda_gc = round(median(chi, na.rm = TRUE) / qchisq(0.5, 1), 3),
      p_lt_05 = round(mean(p[null_i] < 0.05, na.rm = TRUE), 4),
      p_lt_001 = round(mean(p[null_i] < 0.001, na.rm = TRUE), 4))
    if (src == "pooled") qq[[paste(ph, pa)]] <- sort(p[null_i])
  }
  if (src == "pooled") cal[[length(cal)]][, Q_reject_05 := round(mean(r$Q_pval[null_i] < 0.05, na.rm = TRUE), 4)]
}
rec <- rbindlist(rec); cal <- rbindlist(cal, fill = TRUE)
fwrite(rec, file.path(work, "recovery.tsv"), sep = "\t")
fwrite(cal, file.path(work, "null_calibration.tsv"), sep = "\t")

png(file.path(work, "qq.png"), width = 1900, height = 700, res = 120)
par(mfrow = c(2, 7), mar = c(4, 4, 2, 1))
for (nm in names(qq)) {
  o <- -log10(qq[[nm]]); e <- -log10(ppoints(length(o)))
  plot(e, o, pch = 20, cex = 0.4, main = nm, xlab = "expected -log10 p", ylab = "observed")
  abline(0, 1, col = "red")
}
invisible(dev.off())

ok_rec <- rec[model != "mixture", all(abs(z_vs_truth) < 3) && all(p[target] < 1e-4)]
ok_cal <- cal[model != "mixture", all(abs(p_lt_05 - 0.05) <= 3 * sqrt(0.05 * 0.95 / n_null))]
options(width = 160)
cat("\n== Causal SNPs, target parameter (pooled)\n")
print(rec[target == TRUE & cov_source == "pooled", .(model, snp_id, maf, param, truth, est, se, z_vs_truth, p)])
cat("\n== Off-target |z| max per model (truth 0; pooled)\n")
print(rec[target == FALSE & cov_source == "pooled", .(max_abs_z = max(abs(z_vs_truth))), by = .(model, param)])
cat("\n== Null-SNP calibration\n"); print(cal)
cat(sprintf("\nRECOVERY    %s\nCALIBRATION %s\n", if (ok_rec) "PASS" else "FAIL", if (ok_cal) "PASS" else "FAIL"))
