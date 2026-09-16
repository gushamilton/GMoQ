## gamma's exact calibration as a function of mixture separation.
## The protection is not a property of "having a membership column"; it is a
## property of W_gamma being Sigma's leading eigendirection, which requires the
## density trough to be deep enough for Sigma to collapse onto it.
t <- read.delim("/user/work/fh6520/gmoq_covmixW/out/exact_mix_a.tsv", stringsAsFactors = FALSE)
s <- t[t$W == "mix3" & t$param == "gamma" & t$dgp == "disp" & t$maf == 0.01, ]
s <- s[order(s$sep, s$Ztarget), ]
cat("EXACT, mix3 gamma, MAF 0.01, true dispersion\n\n")
cat(sprintf("%-8s %5s %4s %9s | %10s %10s %11s %12s\n",
            "shape", "sep", "Ztg", "c", "vEv/vBv", "se_ratio", "cov95", "infl@5e-8"))
for (i in seq_len(nrow(s))) with(s[i, ], cat(sprintf(
  "%-8s %5g %4g %9.5f | %10.5f %10.5f %11.5f %12.4g\n",
  shape, sep, Ztarget, c, vEv_over_vBv, se_ratio_pooled, cov95_pooled, infl_5e8_pooled)))

cat("\nEXACT, mix3 gamma, worst se_ratio over the WHOLE sweep, by separation\n\n")
g <- t[t$W == "mix3" & t$param == "gamma", ]
cat(sprintf("%-8s %5s | %10s %10s %11s %12s %s\n",
            "shape", "sep", "min se_rat", "max se_rat", "cov95 @min", "infl @min", "where"))
for (sh in unique(g$shape)) {
  gg <- g[g$shape == sh, ]; i <- which.min(gg$se_ratio_pooled)
  cat(sprintf("%-8s %5g | %10.5f %10.5f %11.5f %12.4g %s maf=%g Ztg=%g\n",
              sh, gg$sep[1], gg$se_ratio_pooled[i], max(gg$se_ratio_pooled),
              gg$cov95_pooled[i], gg$infl_5e8_pooled[i], gg$dgp[i], gg$maf[i], gg$Ztarget[i]))
}
