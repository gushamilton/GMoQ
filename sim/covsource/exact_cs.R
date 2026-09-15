## EXACT engine for the cov_source head-to-head.
##
## The pooled estimator's error is a FIXED working-independence bias: it reports
## sqrt(diag(A Sigma_0 A' / d)) when the truth is sqrt(diag(A Sigma_true A' / d)).
## That ratio is available in closed form, with no Monte Carlo error at all, so
## the whole pooled surface -- SE ratio, its direction, coverage, Type-I at 0.05
## and at 5e-8 -- can be computed exactly.
##
## The JACKKNIFE has no such closed form: its error is estimation noise in a
## 45x45 matrix from B = 1000 blocks. Nothing in this file says anything about
## the jackknife; that is what the Monte Carlo is for.
##
## Assumption for the exact Type-I: theta_hat is normal about theta_pop with the
## true SE. That is the same assumption the prior FINDINGS document's exact
## column rests on, and the MC checks it at 0.05 and 1e-3.
WS <- "/user/work/fh6520/covsource_run"
.libPaths(c(file.path(WS, "lib")))
stopifnot(identical(.libPaths()[1], file.path(WS, "lib")))
library(fungwasStage2)
source(file.path(WS, "src", "truth_cs.R"))

TAUS <- seq(0.02, 0.98, length.out = 45)

cells_txt <- paste(readLines(file.path(WS, "cells.json"), warn = FALSE), collapse = "")
blocks <- regmatches(cells_txt, gregexpr("\\{[^{}]*\\}", cells_txt))[[1]]
getf <- function(b, k) {
  m <- regmatches(b, regexpr(paste0("\"", k, "\"[[:space:]]*:[[:space:]]*[^,}]+"), b))
  gsub("\"", "", trimws(sub(paste0("\"", k, "\"[[:space:]]*:[[:space:]]*"), "", m)))
}

zc05 <- qnorm(1 - 0.05 / 2); zc8 <- qnorm(1 - 5e-8 / 2)
rows <- list()
for (i in seq_along(blocks)) {
  b <- blocks[i]
  nm <- getf(b, "name"); dgp <- getf(b, "dgp")
  maf <- as.numeric(getf(b, "maf")); shape <- getf(b, "shape")
  a <- as.numeric(getf(b, "a")); cc <- as.numeric(getf(b, "c"))
  N <- as.integer(getf(b, "N"))
  tr <- pop_truth(TAUS, maf, shape, a, cc)
  W  <- make_weights_vqtl(TAUS, tr$q_tau, mu = tr$EY, sd = tr$sdY)
  d  <- N * tr$Vg
  tt <- cs_theta_truth(W, tr, d)
  ## per-tau pooled/true SE ratio
  ptr <- sqrt(diag(tr$Sigma0) / diag(tr$Sigma_true))
  for (p in c("beta_mu", "beta_sigma2")) {
    rat_ols <- unname(tt$ratio_pooled_ols[p])       # reported / actual, pooled+OLS
    rat_gls <- unname(tt$ratio_pooled_gls[p])       # reported / actual, pooled+GLS
    ## Type-I of the Wald test on this parameter when the truth is zero for it.
    ## |theta|/SE_reported > z  <=>  |Z_true| > z * rat, so rate = 2*pnorm(-z*rat).
    t1_05_ols <- 2 * pnorm(-zc05 * rat_ols); t1_8_ols <- 2 * pnorm(-zc8 * rat_ols)
    t1_05_gls <- 2 * pnorm(-zc05 * rat_gls); t1_8_gls <- 2 * pnorm(-zc8 * rat_gls)
    rows[[length(rows) + 1L]] <- data.frame(
      cell_idx = i - 1L, cell = nm, dgp = dgp, maf = maf, shape = shape,
      a = a, c = cc, sd_ratio_hom = 1 + 2 * cc, MAC = 2 * maf * N, param = p,
      theta_ols = unname(tt$th_ols[p]),
      theta_gls_true = unname(tt$th_gls[p]),
      theta_gls_pooledW = unname(tt$th_gls_pooledW[p]),
      Z_true_ols = unname(tt$th_ols[p] / tt$se_ols_true[p]),
      se_ratio_pooled_ols = rat_ols, se_ratio_pooled_gls = rat_gls,
      cov95_pooled_ols = pnorm(zc05 * rat_ols) - pnorm(-zc05 * rat_ols),
      cov95_pooled_gls = pnorm(zc05 * rat_gls) - pnorm(-zc05 * rat_gls),
      t1_05_pooled_ols = t1_05_ols, t1_5e8_pooled_ols = t1_8_ols,
      infl_5e8_pooled_ols = t1_8_ols / 5e-8,
      t1_05_pooled_gls = t1_05_gls, t1_5e8_pooled_gls = t1_8_gls,
      infl_5e8_pooled_gls = t1_8_gls / 5e-8,
      eff_gls_over_ols_true = unname(tt$se_gls_true[p] / tt$se_ols_true[p]),
      per_tau_ratio_min = min(ptr), per_tau_ratio_max = max(ptr),
      max_resid_ols = max(abs(tt$resid_ols)),
      stringsAsFactors = FALSE)
  }
  cat(".", if (i %% 20 == 0) "\n" else "")
}
cat("\n")
tab <- do.call(rbind, rows)
write.table(tab, file.path(WS, "out", "exact_table.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("wrote", file.path(WS, "out", "exact_table.tsv"), nrow(tab), "rows\n")

## headline: pure-dispersion cells, both parameters
sel <- tab[tab$dgp == "disp" & tab$shape == "normal", ]
sel <- sel[order(sel$maf, sel$Z_true_ols), ]
cat("\n== EXACT: pooled under a TRUE DISPERSION effect (normal) ==\n")
cat(sprintf("%6s %8s %6s %-12s %10s %10s %10s %12s\n",
            "maf", "Zs2", "sdrat", "param", "seR_ols", "cov95_ols", "seR_gls", "infl5e8_ols"))
for (r in seq_len(nrow(sel))) with(sel[r, ], cat(sprintf(
  "%6g %8.1f %6.2f %-12s %10.5f %10.5f %10.5f %12.3g\n",
  maf, Z_true_ols, sd_ratio_hom, param,
  se_ratio_pooled_ols, cov95_pooled_ols, se_ratio_pooled_gls, infl_5e8_pooled_ols)))
