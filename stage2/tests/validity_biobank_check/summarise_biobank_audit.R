## Summariser for the bounded OLS covariance audit.
##
## Every rate here is COUNTED from the shipped Stage 2 output. Monte Carlo error
## is attached to all of them. Comparisons between cov_source values are PAIRED
## (same replicates), so the paired MC SE on a difference is reported alongside
## the unpaired SE on each level.
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "."
WS <- if (length(args) >= 2) args[2] else "/user/work/fh6520/covsource_run"
LIB <- if (length(args) >= 3) args[3] else "/user/work/fh6520/gmoq_finish_20260915/inference/lib"
.libPaths(c(LIB, .libPaths()))
stopifnot(identical(.libPaths()[1], LIB))

files <- sort(list.files(outdir, pattern = "^stage2_cell[0-9]+\\.rds$", full.names = TRUE))
cat("cells:", length(files), "\n")

z975  <- qnorm(0.975)
zc8   <- qnorm(1 - 5e-8 / 2)

## Which population estimand does this configuration target?
estimand_of <- function(tt, est, cs, p) {
  if (est == "ols") return(unname(tt$th_ols[p]))
  if (cs == "jackknife") return(unname(tt$th_gls[p]))       # weights -> Sigma_true
  unname(tt$th_gls_pooledW[p])                              # weights = Sigma_0
}

boot_se_ratio <- function(th, se, nboot = 400L) {
  n <- length(th)
  v <- numeric(nboot)
  for (b in seq_len(nboot)) { i <- sample.int(n, n, TRUE); v[b] <- mean(se[i]) / sd(th[i]) }
  sd(v)
}

rows <- list(); stats_store <- list()
for (f in files) {
  x <- readRDS(f)
  cn <- x$cell$name
  for (key in names(x$res)) {
    kk <- strsplit(key, "|", fixed = TRUE)[[1]]
    fn <- kk[1]; est <- kk[2]; cs <- kk[3]
    r  <- x$res[[key]]
    tt <- x$theta_true[[fn]]
    pn <- rownames(r$theta)
    for (pi in seq_along(pn)) {
      p  <- pn[pi]
      th <- r$theta[pi, ]; se <- r$se[pi, ]; pv <- r$pval[pi, ]; df <- r$df[pi, ]
      ok <- is.finite(th) & is.finite(se) & se > 0
      nfail <- sum(!ok)
      th <- th[ok]; se <- se[ok]; pv <- pv[ok]; df <- df[ok]
      n <- length(th)
      if (n < 100) next
      tgt <- estimand_of(tt, est, cs, p)
      zz  <- (th - tgt) / se                     # pivot for coverage
      TT  <- th / se                             # Wald statistic for Type-I / power
      cov_t <- mean(abs(zz) < qt(0.975, df))
      cov_z <- mean(abs(zz) < z975)
      emp <- sd(th); mse_ <- mean(se)
      ## tail scale of the Wald statistic relative to the package's own reference
      nu0 <- median(df)
      lam_sd  <- sd(TT) / sqrt(nu0 / (nu0 - 2))
      lam_q   <- as.numeric(quantile(abs(TT), 0.999)) / qt(1 - 0.0005, nu0)
      tq <- qt(1 - 2.5e-8, nu0)
      ex5e8 <- function(lam) 2 * pt(-tq / lam, nu0)
      t1_05  <- mean(pv < 0.05)
      t1_1e3 <- mean(pv < 1e-3)
      stats_store[[paste(cn, key, p, sep = "@@")]] <- list(TT = TT, ok = ok)
      rows[[length(rows) + 1L]] <- data.frame(
        cell = cn, cell_idx = x$cell_idx, dgp = x$cell$dgp, maf = x$cell$maf,
        shape = x$cell$shape, Ztarget = x$cell$Ztarget, N = x$N,
        fit = fn, est = est, cov_source = cs, param = p, nrep = n, nfail = nfail,
        estimand = tgt, Z_true = tgt / unname(tt$se_ols_true[p]),
        mean_hat = mean(th), emp_sd = emp, mean_se = mse_,
        se_ratio = mse_ / emp,
        se_ratio_mcse = boot_se_ratio(th, se),
        bias_in_sd = (mean(th) - tgt) / emp,
        cv_se = sd(se) / mse_,
        cov95_t = cov_t, cov95_t_mcse = sqrt(cov_t * (1 - cov_t) / n),
        cov99_t = mean(abs(zz) < qt(0.995, df)),
        cov999_t = mean(abs(zz) < qt(0.9995, df)),
        cov95_z = cov_z, cov95_z_mcse = sqrt(cov_z * (1 - cov_z) / n),
        t1_05 = t1_05, t1_05_mcse = sqrt(t1_05 * (1 - t1_05) / n),
        t1_1e3 = t1_1e3, t1_1e3_mcse = sqrt(max(t1_1e3, 1 / n) * (1 - t1_1e3) / n),
        nu_med = nu0, lambda_sd = lam_sd, lambda_q999 = lam_q,
        t1_5e8_from_sd = ex5e8(lam_sd), t1_5e8_from_q999 = ex5e8(lam_q),
        infl_5e8_from_sd = ex5e8(lam_sd) / 5e-8,
        infl_5e8_from_q999 = ex5e8(lam_q) / 5e-8,
        stringsAsFactors = FALSE)
    }
  }
  cat(".", sep = "")
}
cat("\n")
tab <- do.call(rbind, rows)

## ---- paired pooled vs jackknife -------------------------------------------
pair <- list()
kk <- unique(tab[, c("cell", "fit", "est", "param")])
for (i in seq_len(nrow(kk))) {
  a <- stats_store[[paste(kk$cell[i], paste(kk$fit[i], kk$est[i], "pooled", sep = "|"), kk$param[i], sep = "@@")]]
  b <- stats_store[[paste(kk$cell[i], paste(kk$fit[i], kk$est[i], "jackknife", sep = "|"), kk$param[i], sep = "@@")]]
  if (is.null(a) || is.null(b)) next
  ok <- a$ok & b$ok
  ra <- tab[tab$cell == kk$cell[i] & tab$fit == kk$fit[i] & tab$est == kk$est[i] &
            tab$param == kk$param[i] & tab$cov_source == "pooled", ]
  rb <- tab[tab$cell == kk$cell[i] & tab$fit == kk$fit[i] & tab$est == kk$est[i] &
            tab$param == kk$param[i] & tab$cov_source == "jackknife", ]
  if (!nrow(ra) || !nrow(rb)) next
  ## paired MC SE on the Type-I difference, using the common replicates
  Ta <- a$TT[ok[a$ok]]; Tb <- b$TT[ok[b$ok]]
  m <- min(length(Ta), length(Tb)); Ta <- Ta[seq_len(m)]; Tb <- Tb[seq_len(m)]
  ia <- as.numeric(abs(Ta) > qt(0.975, ra$nu_med)); ib <- as.numeric(abs(Tb) > qt(0.975, rb$nu_med))
  pair[[length(pair) + 1L]] <- data.frame(
    cell = kk$cell[i], fit = kk$fit[i], est = kk$est[i], param = kk$param[i],
    d_cov95_t = rb$cov95_t - ra$cov95_t,
    d_se_ratio = rb$se_ratio - ra$se_ratio,
    d_t1_05 = mean(ib) - mean(ia),
    d_t1_05_paired_mcse = sd(ib - ia) / sqrt(m),
    emp_sd_jk_over_pooled = rb$emp_sd / ra$emp_sd,
    stringsAsFactors = FALSE)
}
ptab <- do.call(rbind, pair)

## ---- power at MATCHED ACTUAL Type-I ----------------------------------------
## For each (maf, shape, fit, est, cov_source) take the matching null cell's
## |Wald| distribution for beta_sigma2, read off the empirical alpha-quantile,
## and apply THAT threshold in the dispersion cells. Equalises the actual, not
## the nominal, Type-I error across configurations.
pw <- list()
nullmap <- tab[tab$dgp == "null", ]
for (i in seq_len(nrow(tab))) {
  rr <- tab[i, ]
  if (rr$dgp != "disp" || rr$param != "beta_sigma2") next
  nc <- nullmap[nullmap$maf == rr$maf & nullmap$shape == rr$shape &
                nullmap$fit == rr$fit & nullmap$est == rr$est &
                nullmap$cov_source == rr$cov_source & nullmap$param == "beta_sigma2", ]
  if (!nrow(nc)) next
  Tn <- stats_store[[paste(nc$cell[1], paste(rr$fit, rr$est, rr$cov_source, sep = "|"), "beta_sigma2", sep = "@@")]]$TT
  Td <- stats_store[[paste(rr$cell, paste(rr$fit, rr$est, rr$cov_source, sep = "|"), "beta_sigma2", sep = "@@")]]$TT
  for (al in c(0.05, 0.01)) {
    thr <- as.numeric(quantile(abs(Tn), 1 - al))
    po  <- mean(abs(Td) > thr)
    pw[[length(pw) + 1L]] <- data.frame(
      cell = rr$cell, maf = rr$maf, shape = rr$shape, Ztarget = rr$Ztarget,
      fit = rr$fit, est = rr$est, cov_source = rr$cov_source, alpha_actual = al,
      thresh = thr, nominal_t1_at_thresh = 2 * pt(-thr, rr$nu_med),
      power = po, power_mcse = sqrt(po * (1 - po) / length(Td)),
      stringsAsFactors = FALSE)
  }
}
ptab_pw <- do.call(rbind, pw)

write.table(tab, file.path(outdir, "summary_cs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(ptab)) write.table(ptab, file.path(outdir, "paired_cs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(ptab_pw)) write.table(ptab_pw, file.path(outdir, "power_cs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote summary_cs.tsv (", nrow(tab), " rows), paired_cs.tsv, power_cs.tsv\n", sep = "")
