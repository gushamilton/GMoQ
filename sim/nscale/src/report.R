## Report tables for the N-scaling study.
##
## Works from the per-cell .rds files rather than summary.tsv, so that pooled
## rates carry an HONEST Monte Carlo error: the (fit, parameter) z-values within
## one replicate are correlated, so the unit of independence is the REPLICATE.
## Every pooled rate below is therefore the mean over replicates of that
## replicate's mean indicator, and its MC se is sd(per-replicate mean)/sqrt(R).
args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]; outdir <- args[2]
.libPaths(c(LIB)); stopifnot(identical(.libPaths()[1], LIB))

files <- sort(list.files(outdir, pattern = "^stage2_cell[0-9]+\\.rds$", full.names = TRUE))
Z <- list()      # per cell/est: matrix R x K of z-values, plus labels
meta <- list()

for (f in files) {
  x <- readRDS(f)
  for (est in c("ols", "gls")) {
    zl <- list(); lab <- character(0); lof <- numeric(0)
    for (fn in names(x$W)) {
      key <- paste(fn, est, sep = "_")
      th <- x$res[[key]]$theta; se <- x$res[[key]]$se
      medSEtau <- median(sqrt(diag(x$truth$Sigma0) / median(x$d_resid)))
      l <- max(abs(x$theta_true[[fn]]$resid_ols)) / medSEtau
      for (p in seq_len(nrow(th))) {
        tt <- x$theta_true[[fn]][[est]][p]
        zl[[length(zl) + 1L]] <- (th[p, ] - tt) / se[p, ]
        lab <- c(lab, paste(fn, rownames(th)[p], sep = ":"))
        lof <- c(lof, l)
      }
    }
    Zm <- do.call(cbind, zl); colnames(Zm) <- lab
    id <- paste(x$cell$name, est, sep = "|")
    Z[[id]] <- Zm
    meta[[id]] <- list(cell = x$cell$name, arm = x$cell$arm, dgp = x$cell$dgp,
                       shape = x$cell$shape, maf = x$cell$maf, N = x$N,
                       mac = 2 * x$N * x$cell$maf, est = est, lof = lof,
                       reps = nrow(Zm))
  }
}

## pooled rate over a set of ids and a column filter, with replicate-clustered MC se
pooled <- function(ids, fun, colkeep = NULL) {
  per_rep <- unlist(lapply(ids, function(id) {
    Zm <- Z[[id]]
    if (!is.null(colkeep)) {
      k <- colkeep(meta[[id]])
      if (!any(k)) return(numeric(0))
      Zm <- Zm[, k, drop = FALSE]
    }
    rowMeans(fun(Zm))
  }))
  if (!length(per_rep)) return(c(est = NA, se = NA, n = 0))
  c(est = mean(per_rep), se = sd(per_rep) / sqrt(length(per_rep)),
    n = length(per_rep))
}
ind95 <- function(Zm) abs(Zm) < qnorm(0.975)
ind99 <- function(Zm) abs(Zm) < qnorm(0.995)

Ns <- sort(unique(vapply(meta, function(m) m$N, 0)))
sel <- function(...) {
  cond <- list(...)
  names(Filter(function(m) all(vapply(names(cond), function(k)
    m[[k]] %in% cond[[k]], TRUE)), meta))
}
fmt <- function(v) if (is.na(v["est"])) "     --      " else
  sprintf("%.4f+-%.4f", v["est"], v["se"])

hdr <- function(s) cat("\n\n============ ", s, " ============\n")

## ---- 1. headline coverage by N x dgp x estimator ---------------------------
for (lv in c(95, 99)) {
  hdr(sprintf("COVERAGE OF %d%% INTERVALS, all fits and parameters pooled", lv))
  fn <- if (lv == 95) ind95 else ind99
  cat(sprintf("%-8s %-6s %-6s %16s %16s\n", "dgp", "arm", "est", "", ""))
  cat(sprintf("%-8s %-6s %-6s", "dgp", "arm", "est"))
  for (N in Ns) cat(sprintf(" %16s", paste0("N=", N)))
  cat("\n")
  for (arm in c("A", "B", "C")) {
    for (dg in c("null", "scale", "mixed")) {
      for (est in c("ols", "gls")) {
        ids0 <- sel(arm = arm, dgp = dg, est = est)
        if (!length(ids0)) next
        cat(sprintf("%-8s %-6s %-6s", dg, arm, est))
        for (N in Ns) {
          ids <- intersect(ids0, sel(N = N))
          cat(sprintf(" %16s", fmt(pooled(ids, fn))))
        }
        cat("\n")
      }
    }
  }
}

## ---- 1b. coverage per model PARAMETER --------------------------------------
allpar <- sort(unique(unlist(lapply(Z, colnames))))
for (dg in c("null", "scale", "mixed")) {
  hdr(sprintf("COVERAGE 95%% PER PARAMETER, arm A (MAF 0.30), dgp = %s", dg))
  cat(sprintf("%-22s %-5s", "fit:param", "est"))
  for (N in Ns) cat(sprintf(" %16s", paste0("N=", N)))
  cat("\n")
  for (pm in allpar) {
    for (est in c("ols", "gls")) {
      ids0 <- sel(arm = "A", dgp = dg, est = est)
      if (!length(ids0)) next
      line <- sprintf("%-22s %-5s", pm, est); any_ok <- FALSE
      for (N in Ns) {
        ids <- intersect(ids0, sel(N = N))
        v <- pooled(ids, ind95, function(m) colnames(Z[[paste(m$cell, m$est,
                                                              sep = "|")]]) == pm)
        if (!is.na(v["est"])) any_ok <- TRUE
        line <- paste0(line, sprintf(" %16s", fmt(v)))
      }
      if (any_ok) cat(line, "\n")
    }
  }
}

## ---- 2. coverage restricted to WELL-FITTING model/cell combinations --------
hdr("COVERAGE 95%, restricted to (fit,param) with lack of fit < 0.5")
keep_good <- function(m) m$lof < 0.5
cat(sprintf("%-8s %-6s %-6s", "dgp", "arm", "est"))
for (N in Ns) cat(sprintf(" %16s", paste0("N=", N)))
cat("\n")
for (arm in c("A", "B", "C")) for (dg in c("null", "scale", "mixed"))
  for (est in c("ols", "gls")) {
    ids0 <- sel(arm = arm, dgp = dg, est = est)
    if (!length(ids0)) next
    cat(sprintf("%-8s %-6s %-6s", dg, arm, est))
    for (N in Ns) cat(sprintf(" %16s",
      fmt(pooled(intersect(ids0, sel(N = N)), ind95, keep_good))))
    cat("\n")
  }

## ---- 3. Type-I and the genome-wide extrapolation under the null ------------
hdr("NULL CELLS: sd(z) and the extrapolation to the 5e-8 threshold")
zc <- qnorm(1 - 5e-8 / 2)
cat(sprintf("%-6s %-6s %-9s %10s %10s %10s %10s %10s\n",
            "arm", "est", "N", "n_z", "sd(z)", "mc_se", "t1@0.05", "p5e8/5e-8"))
for (arm in c("A", "B", "C")) for (est in c("ols", "gls")) for (N in Ns) {
  ids <- sel(arm = arm, dgp = "null", est = est, N = N)
  if (!length(ids)) next
  zz <- unlist(lapply(ids, function(id) as.vector(Z[[id]])))
  s <- sd(zz); n <- length(zz)
  cat(sprintf("%-6s %-6s %-9d %10d %10.4f %10.4f %10.4f %10.3f\n",
              arm, est, N, n, s, s / sqrt(2 * (n - 1)),
              mean(abs(zz) > qnorm(0.975)), 2 * pnorm(-zc / s) / 5e-8))
}

## ---- 4. SE ratio (mean reported SE / empirical sd) -------------------------
hdr("SE RATIO  mean(reported SE)/sd(theta_hat), by N")
tab <- read.delim(file.path(outdir, "summary.tsv"), stringsAsFactors = FALSE)
for (arm in c("A", "B", "C")) for (est in c("ols", "gls")) {
  cat(sprintf("%-6s %-6s", arm, est))
  for (N in Ns) {
    s <- tab[tab$arm == arm & tab$est == est & tab$N == N &
               tab$max_resid_ols_per_SEtau < 0.5, "se_ratio"]
    cat(sprintf(" %16s", if (!length(s)) "  --  " else
      sprintf("%.4f (n=%d)", mean(s), length(s))))
  }
  cat("\n")
}

## ---- 5. efficiency sd(GLS)/sd(OLS) by N ------------------------------------
hdr("EFFICIENCY  empirical sd(GLS)/sd(OLS), by N and parameter")
for (pm in sort(unique(tab$param))) {
  cat(sprintf("%-14s", pm))
  for (N in Ns) {
    s <- tab[tab$param == pm & tab$N == N & tab$est == "ols" &
               tab$max_resid_ols_per_SEtau < 0.5, "eff_sd_gls_over_ols"]
    cat(sprintf(" %16s", if (!length(s)) "  --  " else
      sprintf("%.3f [%.3f,%.3f]", mean(s), min(s), max(s))))
  }
  cat("\n")
}

## ---- 6. per-tau Stage 1 calibration ----------------------------------------
hdr("PER-TAU Stage 1 coverage of beta_tau (95%), extreme taus vs the middle")
pt <- read.delim(file.path(outdir, "pertau.tsv"), stringsAsFactors = FALSE)
show <- c(0.02, 0.06, 0.14, 0.50, 0.86, 0.94, 0.98)
for (arm in c("A", "B", "C")) for (dg in c("null", "scale")) {
  ids <- pt$arm == arm & pt$dgp == dg
  if (!any(ids)) next
  cat(sprintf("\n-- arm %s, dgp %s --\n%-8s", arm, dg, "tau"))
  for (N in Ns) cat(sprintf(" %22s", paste0("N=", N)))
  cat("\n")
  for (tv in show) {
    cat(sprintf("%-8.2f", tv))
    for (N in Ns) {
      r <- pt[ids & pt$N == N & abs(pt$tau - tv) < 0.011, ]
      cat(sprintf(" %22s", if (!nrow(r)) "  --  " else
        sprintf("%.4f+-%.4f n<=%.0f", mean(r$cov95),
                mean(r$mc_se_cov95) / sqrt(nrow(r)), mean(r$n_below))))
    }
    cat("\n")
  }
}

hdr("PER-TAU KDE density bias  mean(fhat)/f_true - 1")
for (arm in c("A", "B", "C")) {
  ids <- pt$arm == arm & pt$dgp == "null"
  if (!any(ids)) next
  cat(sprintf("\n-- arm %s --\n%-8s", arm, "tau"))
  for (N in Ns) cat(sprintf(" %14s", paste0("N=", N)))
  cat("\n")
  for (tv in show) {
    cat(sprintf("%-8.2f", tv))
    for (N in Ns) {
      r <- pt[ids & pt$N == N & abs(pt$tau - tv) < 0.011, ]
      cat(sprintf(" %14s", if (!nrow(r)) " -- " else
        sprintf("%+.4f", mean(r$fhat_bias))))
    }
    cat("\n")
  }
}

hdr("PER-TAU SE ratio at Stage 1  mean(SE_tau)/sd(beta_tau_hat)")
for (arm in c("A", "B", "C")) {
  ids <- pt$arm == arm & pt$dgp == "null"
  if (!any(ids)) next
  cat(sprintf("\n-- arm %s --\n%-8s", arm, "tau"))
  for (N in Ns) cat(sprintf(" %14s", paste0("N=", N)))
  cat("\n")
  for (tv in show) {
    cat(sprintf("%-8.2f", tv))
    for (N in Ns) {
      r <- pt[ids & pt$N == N & abs(pt$tau - tv) < 0.011, ]
      cat(sprintf(" %14s", if (!nrow(r)) " -- " else
        sprintf("%.4f", mean(r$se_ratio_tau))))
    }
    cat("\n")
  }
}

## ---- 7. the block-jackknife df law across N --------------------------------
hdr("EFFECTIVE DF from 1/nu = 1/(B-1) + c/n_eff, and what it does at 5e-8")
cat(sprintf("%-6s %-12s %-9s %10s %10s %12s\n",
            "arm", "param", "N", "nu_jk", "zcrit_t", "t1ratio_t"))
for (arm in c("A", "B", "C")) for (pm in c("beta_mu", "beta_scale")) for (N in Ns) {
  s <- tab[tab$arm == arm & tab$param == pm & tab$N == N & tab$dgp == "null" &
             tab$est == "ols", ]
  if (!nrow(s)) next
  cat(sprintf("%-6s %-12s %-9d %10.1f %10.4f %12.3f\n", arm, pm, N,
              mean(s$nu_jk), mean(s$zcrit_t), mean(s$t1_5e8_ratio_t)))
}
cat("\n(for comparison the pooled path uses nu = N, i.e. zcrit = 5.4514)\n")
