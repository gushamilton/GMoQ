## Summaries for the mean-variance (dispersion) Monte Carlo.
##
## Rates are counted, not modelled.  Two Monte Carlo SEs are given: the binomial
## one over variant-observations (the convention the brief asked for) and a
## replicate-clustered one (SD of the per-replicate rate / sqrt(n_rep)), because
## the 200 variants inside a replicate share a phenotype, a Sigma_0 and an r_av.
##
## The 5e-8 extrapolation follows /user/work/fh6520/gls_coverage/src/summarise.R:
##   rate_5e8 = 2 * pnorm(-qnorm(1 - 5e-8/2) / sd_z)
## with sd_z the SD of the NULL-class z statistics.
.libPaths(c("/user/work/fh6520/gls_coverage/lib", "/user/work/fh6520/r_packages"))
a <- commandArgs(trailingOnly = TRUE)
s2dir <- a[1]; outdir <- a[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
ZC  <- qnorm(0.975)
Z58 <- qnorm(1 - 5e-8 / 2)
CLSNAME <- c("null", "mean", "disp")

load_cell <- function(dir) {
  fs <- sort(list.files(dir, pattern = "\\.rds$", full.names = TRUE))
  if (!length(fs)) return(NULL)
  parts <- lapply(fs, readRDS)
  Rtot <- 0L
  g <- function(f) unlist(lapply(parts, f))
  off <- cumsum(c(0L, vapply(parts, function(p) p$R, 0L)))
  rep_id <- unlist(lapply(seq_along(parts), function(i) parts[[i]]$rep_id + off[i]))
  cat2 <- function(nm) do.call(rbind, lapply(parts, function(p) p$out[[nm]]))
  list(cls = parts[[1]]$cls, M = parts[[1]]$M, N = parts[[1]]$N,
       nrep = off[length(off)], rep_id = rep_id,
       var_id = g(function(p) p$var_id), maf = g(function(p) p$maf),
       delta = g(function(p) p$delta), gamma = g(function(p) p$gamma),
       mom = do.call(rbind, lapply(parts, function(p) p$mom)),
       fhat = do.call(rbind, lapply(parts, function(p) p$fhat)),
       fpop = do.call(rbind, lapply(parts, function(p) p$fpop)),
       verify = g(function(p) p$verify_ols_max_err),
       out = setNames(lapply(names(parts[[1]]$out), cat2), names(parts[[1]]$out)))
}

rate_row <- function(hit, rep_id) {
  n <- length(hit); k <- sum(hit); p <- k / n
  pr <- tapply(hit, rep_id, mean)
  c(n = n, k = k, rate = p, mcse_bin = sqrt(p * (1 - p) / n),
    mcse_clus = sd(pr) / sqrt(length(pr)))
}

summarise_cell <- function(x, cell) {
  rows <- list(); lev <- list(); rav <- list()
  for (est in c("ols", "gls")) {
    o <- x$out
    b <- o$beta_s2[, est]; s <- o$se_s2[, est]
    d <- o$d_res[, est];  sd_ <- o$se_d[, est]; sdr <- o$se_d_rav[, est]
    tb <- o$true_s2[, est]; td <- o$true_d[, est]
    z <- b / s; zd <- d / sd_; zdr <- d / sdr
    zc <- (b - tb) / s; zdc <- (d - td) / sd_
    cls <- x$cls[x$var_id]
    sdz_null <- sd(z[cls == 0], na.rm = TRUE)
    sdz_null_d <- sd(zd[cls == 0], na.rm = TRUE)
    rv <- tapply(o$r_av[, est], x$rep_id, function(v) v[1])
    rav[[est]] <- data.frame(cell = cell, est = est, nrep = length(rv),
                             rav_mean = mean(rv), rav_sd = sd(rv),
                             rav_q05 = quantile(rv, .05), rav_q95 = quantile(rv, .95),
                             rav_jk_se_mean = mean(tapply(o$se_rav[, est], x$rep_id,
                                                          function(v) v[1]), na.rm = TRUE),
                             row.names = NULL)
    for (k in 0:2) {
      i <- which(cls == k)
      ru <- rate_row(abs(z[i]) > ZC, x$rep_id[i])
      rr <- rate_row(abs(zd[i]) > ZC, x$rep_id[i])
      rrr <- rate_row(abs(zdr[i]) > ZC, x$rep_id[i])
      cu <- rate_row(abs(zc[i]) < ZC, x$rep_id[i])
      cc <- rate_row(abs(zdc[i]) < ZC, x$rep_id[i])
      rows[[length(rows) + 1L]] <- data.frame(
        cell = cell, est = est, class = CLSNAME[k + 1], nobs = ru["n"],
        nrep = x$nrep,
        rej_u_k = ru["k"], rej_u = ru["rate"], rej_u_mcse = ru["mcse_bin"],
        rej_u_mcse_clus = ru["mcse_clus"],
        rej_r_k = rr["k"], rej_r = rr["rate"], rej_r_mcse = rr["mcse_bin"],
        rej_r_mcse_clus = rr["mcse_clus"],
        rej_r_ravse_k = rrr["k"], rej_r_ravse = rrr["rate"],
        mean_est_u = mean(b[i]), mean_true_u = mean(tb[i]),
        bias_u = mean(b[i] - tb[i]), bias_u_per_se = mean((b[i] - tb[i]) / s[i]),
        mean_est_r = mean(d[i]), mean_true_r = mean(td[i]),
        bias_r = mean(d[i] - td[i]), bias_r_per_se = mean((d[i] - td[i]) / sd_[i]),
        cov95_u_k = cu["k"], cov95_u = cu["rate"], cov95_u_mcse = cu["mcse_bin"],
        cov95_u_mcse_clus = cu["mcse_clus"],
        cov95_r_k = cc["k"], cov95_r = cc["rate"], cov95_r_mcse = cc["mcse_bin"],
        cov95_r_mcse_clus = cc["mcse_clus"],
        se_over_sd_u = mean(s[i]) / sd(b[i]),
        se_over_sd_u_ctr = mean(s[i]) / sd(b[i] - tb[i]),
        se_over_sd_r = mean(sd_[i]) / sd(d[i]),
        se_over_sd_r_ctr = mean(sd_[i]) / sd(d[i] - td[i]),
        sd_z_null = sdz_null, sd_z_null_resid = sdz_null_d,
        rate_5e8_u = 2 * pnorm(-Z58 / sdz_null),
        rate_5e8_r = 2 * pnorm(-Z58 / sdz_null_d),
        row.names = NULL)
      if (k > 0) {
        ## variant index encodes the designed effect level (fixed across reps)
        grp <- x$var_id[i]
        lev[[length(lev) + 1L]] <- do.call(rbind, lapply(split(seq_along(i), grp), function(j) {
          jj <- i[j]
          data.frame(cell = cell, est = est, class = CLSNAME[k + 1],
                     var = grp[j][1], n = length(jj),
                     mean_true_z = mean(tb[jj] / s[jj]),
                     mean_z_u = mean(z[jj]), mean_z_r = mean(zd[jj]),
                     pow_u = mean(abs(z[jj]) > ZC), pow_r = mean(abs(zd[jj]) > ZC),
                     mean_est = mean(b[jj]), mean_true = mean(tb[jj]),
                     bias_per_se = mean((b[jj] - tb[jj]) / s[jj]),
                     mean_est_r = mean(d[jj]), mean_true_r = mean(td[jj]),
                     bias_r_per_se = mean((d[jj] - td[jj]) / sd_[jj]),
                     cov_u = mean(abs((b[jj] - tb[jj]) / s[jj]) < ZC),
                     cov_r = mean(abs((d[jj] - td[jj]) / sd_[jj]) < ZC),
                     row.names = NULL)
        }))
      }
    }
  }
  list(main = do.call(rbind, rows), lev = do.call(rbind, lev),
       rav = do.call(rbind, rav))
}

cells <- list.dirs(s2dir, recursive = FALSE)
MAIN <- LEV <- RAV <- DIAG <- list()
for (cd in cells) {
  cell <- basename(cd)
  x <- load_cell(cd)
  if (is.null(x)) next
  cat("cell", cell, ": reps =", x$nrep, " verify_ols_max_err =",
      max(x$verify), "\n")
  s <- summarise_cell(x, cell)
  MAIN[[cell]] <- s$main; LEV[[cell]] <- s$lev; RAV[[cell]] <- s$rav
  DIAG[[cell]] <- data.frame(cell = cell, nrep = x$nrep,
                             mean_y = mean(x$mom[, 1]), sd_y = mean(x$mom[, 2]),
                             mean_y_pop = mean(x$mom[, 3]), sd_y_pop = mean(x$mom[, 4]),
                             sd_latent = mean(x$mom[, 5]),
                             fhat_over_fpop_med = median(x$fhat / x$fpop),
                             fhat_over_fpop_min = min(colMeans(x$fhat / x$fpop)),
                             fhat_over_fpop_max = max(colMeans(x$fhat / x$fpop)),
                             rav_target_lognormal = 2 * mean(x$mom[, 1]) * (exp(0.25) - 1),
                             verify_ols_max_err = max(x$verify))
}
wr <- function(l, nm) {
  d <- do.call(rbind, l)
  write.table(d, file.path(outdir, nm), sep = "\t", quote = FALSE, row.names = FALSE)
  d
}
m <- wr(MAIN, "main.tsv"); wr(LEV, "bylevel.tsv"); rv <- wr(RAV, "rav.tsv")
dg <- wr(DIAG, "diag.tsv")
cat("\n== DIAGNOSTICS ==\n"); print(dg, digits = 4)
cat("\n== r_av ==\n"); print(rv, digits = 4)
cat("\n== MAIN ==\n")
print(m[, c("cell", "est", "class", "nobs", "rej_u_k", "rej_u", "rej_u_mcse",
            "rej_r_k", "rej_r", "rej_r_mcse", "cov95_u", "cov95_r",
            "se_over_sd_u", "se_over_sd_r", "sd_z_null", "rate_5e8_u",
            "rate_5e8_r")], digits = 4)
cat("\nwrote to", outdir, "\n")
