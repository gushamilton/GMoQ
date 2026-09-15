## Render the head-to-head tables in markdown from summary_cs.tsv / power_cs.tsv.
WS <- "/user/work/fh6520/covsource_run"
outdir <- file.path(WS, "out")
tab <- read.delim(file.path(outdir, "summary_cs.tsv"), stringsAsFactors = FALSE)
pw  <- read.delim(file.path(outdir, "power_cs.tsv"), stringsAsFactors = FALSE)
ex  <- read.delim(file.path(outdir, "exact_table.tsv"), stringsAsFactors = FALSE)
pr  <- read.delim(file.path(outdir, "paired_cs.tsv"), stringsAsFactors = FALSE)

tab <- tab[tab$fit == "LS", ]
tab$cfg <- paste0(tab$cov_source, "+", tab$est)
CFG <- c("pooled+ols", "jackknife+ols", "pooled+gls", "jackknife+gls")

fm <- function(v, d = 4) formatC(v, format = "f", digits = d)
sg <- function(v, d = 3) formatC(v, format = "g", digits = d)

wide <- function(d, val, mcse = NULL, digits = 4, keys = c("cell","maf","Ztarget","param")) {
  u <- unique(d[, keys, drop = FALSE])
  out <- u
  for (cf in CFG) {
    v <- rep(NA_character_, nrow(u))
    for (i in seq_len(nrow(u))) {
      sel <- d$cell == u$cell[i] & d$cfg == cf
      if ("param" %in% keys) sel <- sel & d$param == u$param[i]
      s <- d[sel, ]
      if (nrow(s) == 1) v[i] <- if (is.null(mcse)) fm(s[[val]], digits)
                                else paste0(fm(s[[val]], digits), " (", fm(s[[mcse]], digits), ")")
    }
    out[[cf]] <- v
  }
  out
}

mdtab <- function(df) {
  cat("|", paste(names(df), collapse = " | "), "|\n")
  cat("|", paste(rep("---", ncol(df)), collapse = " | "), "|\n")
  for (i in seq_len(nrow(df)))
    cat("|", paste(as.character(unlist(df[i, ])), collapse = " | "), "|\n")
  cat("\n")
}

sect <- function(t) cat("\n\n### ", t, "\n\n", sep = "")

## ---------------------------------------------------------------- 1. coverage
sect("1. 95% CI coverage, PURE DISPERSION, normal phenotype (t reference, the package's own)")
d <- tab[tab$dgp == "disp" & tab$shape == "normal", ]
d <- d[order(d$maf, d$Ztarget, d$param), ]
mdtab(wide(d, "cov95_t", "cov95_t_mcse"))

sect("2. 95% CI coverage, NULL, normal")
d <- tab[tab$dgp == "null" & tab$shape == "normal", ]
mdtab(wide(d[order(d$maf, d$param), ], "cov95_t", "cov95_t_mcse"))

sect("3. 95% CI coverage, PURE LOCATION and MIXED, normal")
d <- tab[tab$dgp %in% c("loc", "mixed") & tab$shape == "normal", ]
mdtab(wide(d[order(d$dgp, d$maf, d$Ztarget, d$param), ], "cov95_t", "cov95_t_mcse",
           keys = c("cell","dgp","maf","Ztarget","param")))

sect("4. Skewed phenotype (standardised gamma(2)) -- coverage")
d <- tab[tab$shape == "gamma2", ]
mdtab(wide(d[order(d$dgp, d$maf, d$param), ], "cov95_t", "cov95_t_mcse",
           keys = c("cell","dgp","maf","param")))

## -------------------------------------------------------------- 5. SE ratio
sect("5. SE ratio (mean reported SE / empirical sd). <1 = ANTI-CONSERVATIVE")
d <- tab[tab$shape == "normal", ]
mdtab(wide(d[order(d$dgp, d$maf, d$Ztarget, d$param), ], "se_ratio", "se_ratio_mcse",
           keys = c("cell","dgp","maf","Ztarget","param")))

## -------------------------------------------------------------- 6. Type-I
sect("6. Type-I at nominal 0.05 -- NULL cells (both parameters truly zero)")
d <- tab[tab$dgp == "null", ]
mdtab(wide(d[order(d$shape, d$maf, d$param), ], "t1_05", "t1_05_mcse",
           keys = c("cell","shape","maf","param")))

sect("7. Type-I at nominal 0.05 on beta_mu -- under a TRUE DISPERSION effect (beta_mu truly 0)")
d <- tab[tab$dgp == "disp" & tab$param == "beta_mu", ]
mdtab(wide(d[order(d$shape, d$maf, d$Ztarget), ], "t1_05", "t1_05_mcse",
           keys = c("cell","shape","maf","Ztarget","param")))

sect("8. FALSE-vQTL rate at nominal 0.05: beta_sigma2 under a PURE LOCATION effect (zero true dispersion)")
d <- tab[tab$dgp == "loc" & tab$param == "beta_sigma2", ]
mdtab(wide(d[order(d$shape, d$maf, d$Ztarget), ], "t1_05", "t1_05_mcse",
           keys = c("cell","shape","maf","Ztarget","param")))

## ------------------------------------------------- 9. extrapolated to 5e-8
sect("9. Extrapolated Type-I at 5e-8, beta_mu under a TRUE DISPERSION effect (inflation x nominal)")
d <- tab[tab$dgp == "disp" & tab$param == "beta_mu", ]
cat("Two extrapolations per cell: from sd(T) and from the 99.9th percentile of |T|.\n\n")
mdtab(wide(d[order(d$maf, d$Ztarget), ], "infl_5e8_from_sd", NULL, 2,
           keys = c("cell","shape","maf","Ztarget")))
cat("From the 99.9th percentile:\n\n")
mdtab(wide(d[order(d$maf, d$Ztarget), ], "infl_5e8_from_q999", NULL, 2,
           keys = c("cell","shape","maf","Ztarget")))

sect("10. MC validation of the extrapolation: measured Type-I at nominal 1e-3 (NULL cells)")
d <- tab[tab$dgp == "null", ]
mdtab(wide(d[order(d$shape, d$maf, d$param), ], "t1_1e3", "t1_1e3_mcse", 5,
           keys = c("cell","shape","maf","param")))

## -------------------------------------------------------------- 11. power
sect("11. Power to detect a true dispersion effect at MATCHED ACTUAL Type-I = 0.05")
pw$cfg <- paste0(pw$cov_source, "+", pw$est)
p5 <- pw[pw$alpha_actual == 0.05 & pw$fit == "LS", ]
u <- unique(p5[, c("cell","maf","shape","Ztarget")])
u <- u[order(u$shape, u$maf, u$Ztarget), ]
for (cf in CFG) u[[cf]] <- vapply(seq_len(nrow(u)), function(i) {
  s <- p5[p5$cell == u$cell[i] & p5$cfg == cf, ]
  if (nrow(s) == 1) paste0(fm(s$power, 4), " (", fm(s$power_mcse, 4), ")") else NA_character_
}, character(1))
mdtab(u)

sect("12. Power at MATCHED ACTUAL Type-I = 0.01")
p1 <- pw[pw$alpha_actual == 0.01 & pw$fit == "LS", ]
u <- unique(p1[, c("cell","maf","shape","Ztarget")]); u <- u[order(u$shape, u$maf, u$Ztarget), ]
for (cf in CFG) u[[cf]] <- vapply(seq_len(nrow(u)), function(i) {
  s <- p1[p1$cell == u$cell[i] & p1$cfg == cf, ]
  if (nrow(s) == 1) paste0(fm(s$power, 4), " (", fm(s$power_mcse, 4), ")") else NA_character_
}, character(1))
mdtab(u)

## ------------------------------------------------ 13. paired jk - pooled
sect("13. PAIRED jackknife - pooled, same replicates (OLS only), coverage and Type-I")
q <- pr[pr$fit == "LS" & pr$est == "ols", ]
q <- q[order(q$cell, q$param), ]
q$d_cov95_t <- fm(q$d_cov95_t, 4); q$d_se_ratio <- fm(q$d_se_ratio, 4)
q$d_t1_05 <- paste0(fm(q$d_t1_05, 4), " (", fm(q$d_t1_05_paired_mcse, 4), ")")
q$emp_sd_jk_over_pooled <- fm(q$emp_sd_jk_over_pooled, 5)
mdtab(q[, c("cell","param","d_cov95_t","d_se_ratio","d_t1_05","emp_sd_jk_over_pooled")])

## ------------------------------------------- 14. MC vs exact, pooled only
sect("14. MC vs EXACT, pooled+OLS (validates the engine and the shipped path against each other)")
d <- tab[tab$cfg == "pooled+ols" & tab$shape == "normal", ]
m <- merge(d[, c("cell","param","cov95_t","cov95_t_mcse","se_ratio","se_ratio_mcse","maf","Ztarget","dgp")],
           ex[, c("cell","param","cov95_pooled_ols","se_ratio_pooled_ols")],
           by = c("cell","param"))
m <- m[order(m$dgp, m$maf, m$Ztarget, m$param), ]
m$dev_sigma <- fm((m$cov95_t - m$cov95_pooled_ols) / m$cov95_t_mcse, 2)
m$cov95_MC <- paste0(fm(m$cov95_t, 4), " (", fm(m$cov95_t_mcse, 4), ")")
m$cov95_exact <- fm(m$cov95_pooled_ols, 4)
m$seR_MC <- fm(m$se_ratio, 4); m$seR_exact <- fm(m$se_ratio_pooled_ols, 4)
mdtab(m[, c("cell","param","cov95_MC","cov95_exact","dev_sigma","seR_MC","seR_exact")])

sect("15. Realism: what effect each dispersion cell actually requires")
e2 <- ex[ex$dgp == "disp" & ex$param == "beta_sigma2", c("cell","maf","MAC","c","sd_ratio_hom","Z_true_ols")]
e2$c <- fm(e2$c, 4); e2$sd_ratio_hom <- fm(e2$sd_ratio_hom, 2); e2$Z_true_ols <- fm(e2$Z_true_ols, 1)
mdtab(e2[order(e2$maf, e2$Z_true_ols), ])

sect("16. FALSE-vQTL rate extrapolated to 5e-8 (beta_sigma2 under a pure location effect)")
cat("The jackknife is NOT a fix here: the artefact is in the estimand, not the covariance.\n\n")
d <- tab[tab$dgp == "loc" & tab$param == "beta_sigma2", ]
mdtab(wide(d[order(d$shape, d$maf, d$Ztarget), ], "infl_5e8_from_sd", NULL, 2,
           keys = c("cell","shape","maf","Ztarget","param")))

sect("17. Relative noise in the reported SE, sd(SE)/mean(SE) -- the jackknife's only real cost")
d <- tab[tab$param == "beta_sigma2", ]
mdtab(wide(d[order(d$dgp, d$shape, d$maf, d$Ztarget), ], "cv_se", NULL, 4,
           keys = c("cell","dgp","maf","Ztarget","param")))

sect("18. Deep coverage: 99.9% nominal intervals")
d <- tab[tab$dgp %in% c("disp","null") & tab$shape == "normal", ]
mdtab(wide(d[order(d$dgp, d$maf, d$Ztarget, d$param), ], "cov999_t", NULL, 4,
           keys = c("cell","dgp","maf","Ztarget","param")))
