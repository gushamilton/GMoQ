## Assemble the per-cell Stage 2 output of the N-scaling study into two tables:
##   summary.tsv  -- one row per (cell, fit, estimator, parameter)
##   pertau.tsv   -- one row per (cell, tau), Stage 1 calibration
args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]; outdir <- args[2]
.libPaths(c(LIB)); stopifnot(identical(.libPaths()[1], LIB))

files <- sort(list.files(outdir, pattern = "^stage2_cell[0-9]+\\.rds$", full.names = TRUE))
rows <- list(); zrows <- list(); ptrows <- list()
zcrit <- qnorm(1 - 5e-8 / 2)   # 5.4514, the normal genome-wide threshold

for (f in files) {
  x <- readRDS(f)
  medSEtau <- median(sqrt(diag(x$truth$Sigma0) / median(x$d_resid)))
  mac <- 2 * x$N * x$cell$maf
  sds <- list()
  for (key in names(x$res)) {
    fn <- sub("_.*", "", key); est <- sub(".*_", "", key)
    th <- x$res[[key]]$theta; se <- x$res[[key]]$se
    n <- ncol(th)
    for (p in seq_len(nrow(th))) {
      e <- th[p, ]; s <- se[p, ]
      pname <- rownames(th)[p]
      tt_own <- x$theta_true[[fn]][[est]][p]
      tt_ols <- x$theta_true[[fn]][["ols"]][p]
      emp <- sd(e)
      sds[[paste(fn, est, p)]] <- emp
      zz <- (e - tt_own) / s
      zc <- (e - mean(e)) / s
      cov_at <- function(lv, zv = zz) mean(abs(zv) < qnorm(1 - (1 - lv) / 2))
      ## degrees of freedom from the block-jackknife law, 1/nu = 1/(B-1) + c/n_eff
      nu_jk <- if (!is.null(x$ref_df_med[[fn]])) unname(x$ref_df_med[[fn]][p]) else NA_real_
      ## coverage when the reference is t_nu_jk rather than the pooled path's normal
      cov95_t <- if (is.finite(nu_jk)) mean(abs(zz) < qt(0.975, nu_jk)) else NA_real_
      cov99_t <- if (is.finite(nu_jk)) mean(abs(zz) < qt(0.995, nu_jk)) else NA_real_
      ## |z| that the t_nu reference calls genome-wide significant
      zcrit_t <- if (is.finite(nu_jk)) -qt(5e-8 / 2, nu_jk) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        cell = x$cell$name, cell_idx = x$cell_idx, arm = x$cell$arm,
        dgp = x$cell$dgp, maf = x$cell$maf, mac = mac, shape = x$cell$shape,
        N = x$N, n_blocks = x$n_blocks, fit = fn, est = est,
        param = pname, nrep = n,
        theta_true_own = tt_own, theta_true_ols = tt_ols,
        mean_hat = mean(e), emp_sd = emp, mean_se = mean(s), sd_se = sd(s),
        se_ratio = mean(s) / emp,
        bias_sd = (mean(e) - tt_own) / emp,
        cov95 = cov_at(0.95), cov99 = cov_at(0.99), cov999 = cov_at(0.999),
        cov95_centred = cov_at(0.95, zc), cov99_centred = cov_at(0.99, zc),
        nu_jk = nu_jk, cov95_t = cov95_t, cov99_t = cov99_t, zcrit_t = zcrit_t,
        sd_z = sd(zz), sd_z_centred = sd(zc),
        p_gt3 = mean(abs(zz) > 3), p_gt4 = mean(abs(zz) > 4),
        t1_05 = mean(abs(e / s) > qnorm(0.975)),
        t1_1e3 = mean(abs(e / s) > qnorm(1 - 5e-4)),
        t1_1e4 = mean(abs(e / s) > qnorm(1 - 5e-5)),
        max_resid_ols = max(abs(x$theta_true[[fn]]$resid_ols)),
        max_resid_ols_per_SEtau = max(abs(x$theta_true[[fn]]$resid_ols)) / medSEtau,
        an_se_ratio_gls_ols = x$theta_true[[fn]]$se_gls_unit[p] /
          x$theta_true[[fn]]$se_ols_unit[p],
        stringsAsFactors = FALSE)
    }
  }
  for (fn in unique(sub("_.*", "", names(x$res)))) {
    K <- nrow(x$res[[paste0(fn, "_ols")]]$theta)
    for (p in seq_len(K)) {
      zrows[[length(zrows) + 1L]] <- data.frame(
        cell = x$cell$name, fit = fn,
        param = rownames(x$res[[paste0(fn, "_ols")]]$theta)[p],
        eff_sd_gls_over_ols = sds[[paste(fn, "gls", p)]] / sds[[paste(fn, "ols", p)]],
        stringsAsFactors = FALSE)
    }
  }
  pt <- x$pertau
  pt$cell <- x$cell$name; pt$arm <- x$cell$arm; pt$dgp <- x$cell$dgp
  pt$maf <- x$cell$maf; pt$mac <- mac; pt$shape <- x$cell$shape; pt$N <- x$N
  ## observations below the tau-th quantile, the quantity that makes the
  ## extreme taus N-sensitive by construction
  pt$n_below <- pt$tau * x$N
  ptrows[[length(ptrows) + 1L]] <- pt
}

tab <- do.call(rbind, rows)
eff <- do.call(rbind, zrows)
tab <- merge(tab, eff, by = c("cell", "fit", "param"), all.x = TRUE)
tab$mc_se_cov95 <- sqrt(tab$cov95 * (1 - tab$cov95) / tab$nrep)
tab$mc_se_cov99 <- sqrt(tab$cov99 * (1 - tab$cov99) / tab$nrep)
tab$mc_se_se_ratio <- tab$se_ratio / sqrt(2 * (tab$nrep - 1))
tab$mc_se_sd_z <- tab$sd_z / sqrt(2 * (tab$nrep - 1))
## genome-wide extrapolation: normal reference (what the pooled path uses)
tab$t1_5e8_implied <- 2 * pnorm(-zcrit / tab$sd_z)
tab$t1_5e8_ratio <- tab$t1_5e8_implied / 5e-8
## genome-wide extrapolation: t_nu reference from the block-jackknife law
tab$t1_5e8_implied_t <- 2 * pnorm(-tab$zcrit_t / tab$sd_z)
tab$t1_5e8_ratio_t <- tab$t1_5e8_implied_t / 5e-8

ptab <- do.call(rbind, ptrows)
ptab$mc_se_cov95 <- sqrt(ptab$cov95 * (1 - ptab$cov95) / ptab$nrep)

write.table(tab, file.path(outdir, "summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(ptab, file.path(outdir, "pertau.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("summary rows:", nrow(tab), " pertau rows:", nrow(ptab), "->", outdir, "\n")
