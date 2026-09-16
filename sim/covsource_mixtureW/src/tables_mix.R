WS <- "/user/work/fh6520/gmoq_covmixW"
tab <- read.delim(file.path(WS, "out", "exact_mix_a.tsv"), stringsAsFactors = FALSE)
pub <- read.delim(file.path(WS, "out", "exact_table.tsv"), stringsAsFactors = FALSE)
fmt <- function(x, d = 4) formatC(x, format = "f", digits = d)
g <- function(x) formatC(x, format = "g", digits = 3)

cat("\n########## T1. REPRODUCTION: 2-column vqtl W, normal, TRUE DISPERSION ##########\n")
cat("Published FINDINGS_covsource_default.md table is 15,000-rep Monte Carlo;\n")
cat("this column is the EXACT asymptotic value from the same engine.\n\n")
p <- pub[pub$dgp == "disp" & pub$shape == "normal" & pub$param == "beta_sigma2", ]
cat(sprintf("%6s %5s | %-22s %-22s\n", "MAF", "Z", "cov95 pooled+OLS", "SE ratio pooled+OLS"))
cat(sprintf("%6s %5s | %-10s %-11s %-10s %-11s\n", "", "", "exact", "published MC", "exact", "published MC"))
ref <- list("0.3|16"=c(0.9467,0.9900), "0.05|16"=c(0.9285,0.9206), "0.01|16"=c(0.9027,0.8442),
            "0.01|30"=c(0.8475,0.7339), "0.005|16"=c(0.8742,0.7862), "0.005|30"=c(0.7919,0.6519))
for (k in names(ref)) {
  s <- strsplit(k, "|", fixed = TRUE)[[1]]; mf <- as.numeric(s[1]); zz <- as.numeric(s[2])
  r <- p[abs(p$maf - mf) < 1e-9 & abs(p$Z_true_ols - zz) < 0.5, ]
  cat(sprintf("%6g %5g | %-10s %-11s %-10s %-11s\n", mf, zz,
              fmt(r$cov95_pooled_ols[1]), fmt(ref[[k]][1]),
              fmt(r$se_ratio_pooled_ols[1]), fmt(ref[[k]][2])))
}
cat("\n-- Type-I inflation at 5e-8 on beta_mu under a TRUE DISPERSION effect --\n")
pm <- pub[pub$dgp == "disp" & pub$shape == "normal" & pub$param == "beta_mu", ]
cat(sprintf("%6s %7s | %-12s %-14s\n", "MAF", "sd rat", "exact", "published MC"))
for (k in list(c(0.3, 1.06, "1.26x"), c(0.01, 1.67, "295x"), c(0.005, 2.17, "2086x"))) {
  r <- pm[abs(pm$maf - as.numeric(k[1])) < 1e-9 & abs(pm$sd_ratio_hom - as.numeric(k[2])) < 0.02, ]
  cat(sprintf("%6g %7s | %-12s %-14s\n", as.numeric(k[1]), k[2],
              paste0(g(r$infl_5e8_pooled_ols[1]), "x"), k[3]))
}

cat("\n\n########## T2. SAME DGP, DIFFERENT W: does the 3-column W inherit the failure? ##########\n")
cat("shape mix5.89 = testosterone-like separation.  All rows are EXACT (no MC error).\n")
cat("se_ratio = reported pooled SE / true SE.  <1 is anti-conservative.\n\n")
sel <- tab[tab$shape == "mix5.89" & tab$dgp == "disp", ]
sel <- sel[order(sel$maf, sel$Ztarget, sel$W), ]
cat(sprintf("%6s %4s %7s | %-6s %-2s %-12s %8s %9s %9s %12s\n",
            "MAF", "Ztg", "sd rat", "W", "K", "param", "Z_true", "se_ratio", "cov95", "infl@5e-8"))
for (i in seq_len(nrow(sel))) with(sel[i, ] , cat(sprintf(
  "%6g %4g %7.3f | %-6s %-2d %-12s %8.1f %9.4f %9.4f %12s\n",
  maf, Ztarget, sd_ratio_hom, W, K, param, Z_true, se_ratio_pooled, cov95_pooled,
  paste0(g(infl_5e8_pooled), "x"))))

cat("\n\n########## T3. NULL and PURE LOCATION, mix5.89 ##########\n")
sel <- tab[tab$shape == "mix5.89" & tab$dgp %in% c("null", "loc"), ]
sel <- sel[order(sel$dgp, sel$maf, sel$Ztarget, sel$W), ]
cat(sprintf("%-5s %6s %4s | %-6s %-2s %-12s %8s %9s %9s %12s\n",
            "dgp", "MAF", "Ztg", "W", "K", "param", "Z_true", "se_ratio", "cov95", "infl@5e-8"))
for (i in seq_len(nrow(sel))) with(sel[i, ], cat(sprintf(
  "%-5s %6g %4g | %-6s %-2d %-12s %8.1f %9.4f %9.4f %12s\n",
  dgp, maf, Ztarget, W, K, param, Z_true, se_ratio_pooled, cov95_pooled,
  paste0(g(infl_5e8_pooled), "x"))))

cat("\n\n########## T4. SEPARATION LADDER, MAF 0.01, TRUE DISPERSION, Ztarget 30 ##########\n")
sel <- tab[tab$dgp == "disp" & tab$maf == 0.01 & tab$Ztarget == 30, ]
sel <- sel[order(sel$sep, sel$W), ]
cat(sprintf("%-8s %5s | %-6s %-2s %-12s %8s %9s %9s %12s %8s %8s\n",
            "shape", "sep", "W", "K", "param", "Z_true", "se_ratio", "cov95", "infl@5e-8",
            "tau_bar", "partic"))
for (i in seq_len(nrow(sel))) with(sel[i, ], cat(sprintf(
  "%-8s %5s | %-6s %-2d %-12s %8.1f %9.4f %9.4f %12s %8.3f %8.3f\n",
  shape, ifelse(is.na(sep), "-", sep), W, K, param, Z_true, se_ratio_pooled,
  cov95_pooled, paste0(g(infl_5e8_pooled), "x"), tau_bar, participation)))

cat("\n\n########## T5. MECHANISM: where each column's v_j sits, and Delta ##########\n")
cat("v_j = diag(1/f) W (W'W)^-1 e_j is the direction the SE ratio integrates over.\n")
cat("(se_ratio)^2 = 1 / (1 + v'Ev / v'Bv);  E[s,t] = D(min) - t_s D_t - t_t D_s.\n")
cat("Delta(tau) is antisymmetric about tau=0.5 under a scale effect on a symmetric\n")
cat("phenotype, so a v_j concentrated at the median sees v'Ev ~ 0.\n\n")
sel <- tab[tab$dgp == "disp" & tab$maf == 0.01 & tab$Ztarget == 30 &
           tab$shape %in% c("normal", "mix2", "mix4", "mix5.89", "mix6"), ]
sel <- sel[order(sel$sep, sel$W), ]
cat(sprintf("%-8s %5s | %-6s %-12s %9s %9s %10s %10s\n",
            "shape", "sep", "W", "param", "tau_bar", "partic", "vEv/vBv", "se_ratio"))
for (i in seq_len(nrow(sel))) with(sel[i, ], cat(sprintf(
  "%-8s %5s | %-6s %-12s %9.4f %9.4f %10.5f %10.4f\n",
  shape, ifelse(is.na(sep), "-", sep), W, param, tau_bar, participation,
  vEv_over_vBv, se_ratio_pooled)))

cat("\n\n########## T6. WORST se_ratio anywhere in the exact sweep, per (W, dgp) ##########\n")
for (dg in c("null", "loc", "disp")) {
  s <- tab[tab$dgp == dg, ]
  for (wn in c("vqtl2", "mix2", "mix3")) {
    ss <- s[s$W == wn, ]; if (!nrow(ss)) next
    i <- which.min(ss$se_ratio_pooled)
    cat(sprintf("%-5s %-6s K=%d  worst se_ratio %.4f  (%s maf=%g Ztg=%g param=%s)  cov95=%.4f  infl@5e-8=%sx\n",
                dg, wn, ss$K[i], ss$se_ratio_pooled[i], ss$shape[i], ss$maf[i],
                ss$Ztarget[i], ss$param[i], ss$cov95_pooled[i], g(ss$infl_5e8_pooled[i])))
  }
}
cat("\n(exact sweep:", nrow(tab), "rows,", length(unique(paste(tab$shape, tab$maf, tab$dgp, tab$Ztarget))), "cells)\n")
