#!/usr/bin/env Rscript
# Aggregate the per-chunk RDS files into calibration tables.
suppressWarnings(suppressMessages({
  .libPaths(c(Sys.getenv("POOLCAL_RLIB"), "/user/work/fh6520/r_packages", .libPaths()))
}))

args <- commandArgs(trailingOnly = TRUE)
rds_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
cells <- unique(sub("\\.chunk[0-9]+\\.rds$", "", basename(files)))
cat(length(files), "chunk files,", length(cells), "cells\n")

mcse_p <- function(p, n) sqrt(pmax(p * (1 - p), 0) / n)

pool_cell <- function(cf) {
  parts <- lapply(cf, readRDS)
  tr <- parts[[1]]$truth
  Tt <- parts[[1]]$meta$Tt
  N <- parts[[1]]$meta$N
  beta <- do.call(rbind, lapply(parts, `[[`, "beta"))
  se_tau <- lapply(c("pooled", "jackknife", "hc0"), function(s)
    do.call(rbind, lapply(parts, function(p) p$se_tau[[s]])))
  names(se_tau) <- c("pooled", "jackknife", "hc0")
  keys <- names(parts[[1]]$fits)
  fits <- lapply(keys, function(k) {
    list(params = do.call(cbind, lapply(parts, function(p) p$fits[[k]]$params)),
         se     = do.call(cbind, lapply(parts, function(p) p$fits[[k]]$se)),
         p      = do.call(cbind, lapply(parts, function(p) p$fits[[k]]$p)),
         df     = do.call(cbind, lapply(parts, function(p) p$fits[[k]]$df)))
  })
  names(fits) <- keys
  list(truth = tr, Tt = Tt, N = N, beta = beta, se_tau = se_tau, fits = fits,
       d = unlist(lapply(parts, `[[`, "d")))
}

par_rows <- list(); tau_rows <- list(); pop_rows <- list()

for (cl in cells) {
  cf <- files[sub("\\.chunk[0-9]+\\.rds$", "", basename(files)) == cl]
  X <- pool_cell(cf)
  tr <- X$truth; Tt <- X$Tt; N <- X$N
  Rn <- nrow(X$beta)

  ## ---- exact population covariances --------------------------------------
  Sp <- tr$S0pop / (N * tr$var_g)             # pooled, per variant
  Sh <- tr$Shcpop / (N * tr$var_g^2)          # exact HC0 sandwich, per variant
  dratio <- diag(Sp) / diag(Sh)               # per-tau variance ratio, exact
  # worst case over ALL linear contrasts: eigenvalues of Sh^-1/2 Sp Sh^-1/2
  eh <- eigen(Sh, symmetric = TRUE)
  ev <- pmax(eh$values, 0)
  keep <- ev > max(ev) * 1e-10
  Hh <- eh$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(ev[keep]), sum(keep))
  lam <- eigen(t(Hh) %*% Sp %*% Hh, symmetric = TRUE, only.values = TRUE)$values
  pop_rows[[length(pop_rows) + 1]] <- data.frame(
    cell = cl, N = N, maf = tr$maf, delta = tr$delta, gamma = tr$gamma,
    diag_ratio_med = median(dratio), diag_ratio_min = min(dratio),
    diag_ratio_max = max(dratio),
    worst_se_ratio_lo = sqrt(min(lam)), worst_se_ratio_hi = sqrt(max(lam)),
    stringsAsFactors = FALSE)

  ## ---- per-tau Monte Carlo ------------------------------------------------
  emp_sd_tau <- apply(X$beta, 2, sd)
  for (s in names(X$se_tau)) {
    ms <- colMeans(X$se_tau[[s]])
    tau_rows[[length(tau_rows) + 1]] <- data.frame(
      cell = cl, source = s, tau = tr$taus, tau_idx = seq_len(Tt),
      mean_se = ms, emp_sd = emp_sd_tau, se_ratio = ms / emp_sd_tau,
      se_ratio_mc = (ms / emp_sd_tau) / sqrt(2 * (Rn - 1)),
      pop_var_ratio = if (s == "pooled") dratio else NA_real_,
      mean_beta = colMeans(X$beta), b_pop = tr$b_pop, b_dgp = tr$b_dgp,
      n = Rn, stringsAsFactors = FALSE)
  }

  ## ---- parameters ---------------------------------------------------------
  for (k in names(X$fits)) {
    sp <- strsplit(k, ".", fixed = TRUE)[[1]]
    wname <- sp[1]; src <- sp[2]; est <- sp[3]
    W <- if (wname == "both") tr$W_both else tr$W_add
    # The estimator's own projection: OLS is fixed; GLS weights by the
    # covariance the chosen source supplies (pooled uses Sigma_0, the
    # jackknife and HC0 both converge on the true sandwich Sh).
    Swt <- if (src == "pooled") Sp else Sh
    Aeff <- if (est == "gls")
      solve(t(W) %*% solve(Swt, W), t(W) %*% solve(Swt, diag(nrow(W))))
    else solve(crossprod(W), t(W))
    th_pop <- as.numeric(Aeff %*% tr$b_pop)   # exact estimand of THIS estimator
    th_dgp <- as.numeric(Aeff %*% tr$b_dgp)   # first-order DGP truth
    # What the source REPORTS as the variance, against the TRUE variance Sh.
    vrep <- if (src == "pooled") diag(Aeff %*% Sp %*% t(Aeff))
            else diag(Aeff %*% Sh %*% t(Aeff))
    vtrue <- diag(Aeff %*% Sh %*% t(Aeff))
    vp <- vrep; vh <- vtrue

    f <- X$fits[[k]]
    if (is.null(f$params)) next
    K <- nrow(f$params)
    for (i in seq_len(K)) {
      th <- f$params[i, ]; se <- f$se[i, ]
      dfv <- if (!is.null(f$df)) f$df[i, ] else rep(Inf, length(th))
      ok <- is.finite(th) & is.finite(se) & se > 0
      th <- th[ok]; se <- se[ok]; dfv <- dfv[ok]; n <- length(th)
      crit95 <- qt(0.975, dfv); crit99 <- qt(0.995, dfv)
      pv <- if (!is.null(f$p)) f$p[i, ][ok] else 2 * pt(-abs(th / se), dfv)
      esd <- sd(th); msd <- mean(se)
      cov95p <- mean(abs(th - th_pop[i]) <= crit95 * se)
      cov99p <- mean(abs(th - th_pop[i]) <= crit99 * se)
      cov95d <- mean(abs(th - th_dgp[i]) <= crit95 * se)
      cov99d <- mean(abs(th - th_dgp[i]) <= crit99 * se)
      # z relative to the DGP truth, for the false-positive extrapolation
      zc <- (th - th_dgp[i]) / se
      par_rows[[length(par_rows) + 1]] <- data.frame(
        cell = cl, model = wname, source = src, estimator = est,
        param = rownames(f$params)[i], n = n, N = N,
        theta_dgp = th_dgp[i], theta_pop = th_pop[i],
        mean_theta = mean(th), bias_dgp = mean(th) - th_dgp[i],
        bias_pop = mean(th) - th_pop[i],
        emp_sd = esd, mean_se = msd, se_ratio = msd / esd,
        se_ratio_mc = (msd / esd) / sqrt(2 * (n - 1)),
        pop_se_ratio_vs_hc0 = sqrt(vp[i] / vh[i]),
        cov95_pop = cov95p, cov95_pop_mc = mcse_p(cov95p, n),
        cov99_pop = cov99p, cov99_pop_mc = mcse_p(cov99p, n),
        cov95_dgp = cov95d, cov95_dgp_mc = mcse_p(cov95d, n),
        cov99_dgp = cov99d, cov99_dgp_mc = mcse_p(cov99d, n),
        rej05 = mean(pv < 0.05), rej01 = mean(pv < 0.01),
        rej001 = mean(pv < 1e-3), rej1e4 = mean(pv < 1e-4),
        z_mean = mean(zc), z_sd = sd(zc),
        mean_z_abs = mean(abs(th / se)),
        stringsAsFactors = FALSE)
    }
  }
  cat("done", cl, "R =", Rn, "\n")
}

pr <- do.call(rbind, par_rows)
tr2 <- do.call(rbind, tau_rows)
write.table(pr, file.path(out_dir, "params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tr2, file.path(out_dir, "tau.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pp <- do.call(rbind, pop_rows)
write.table(pp, file.path(out_dir, "popexact.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", nrow(pr), "param rows,", nrow(tr2), "tau rows,", nrow(pp), "pop rows\n")
