## Stage 2 driver for the N-scaling study.
##
## Runs the SHIPPED param_gwas() (OLS and GLS) on simulated Stage 1 output
## carrying the POOLED covariance, and stores per-replicate theta_hat / SE so
## that coverage can be COUNTED against the exact population truth.
##
## Extends gls_coverage/src/stage2_driver.R with:
##   * workspace / library paths as arguments;
##   * harvest of the block-jackknife effective-df law (n_blocks + n_eff), so
##     the t_nu reference can be compared against the pooled path's nu = N
##     ("numerically normal") reference at every rung of the N ladder;
##   * per-tau Stage 1 calibration accumulators (coverage of beta_tau itself
##     against the exact population beta_tau), which is where an N effect at
##     the extreme taus would first appear;
##   * KDE density accumulators (fhat against the exact population f_tau).
args     <- commandArgs(trailingOnly = TRUE)
WS       <- args[1]
LIB      <- args[2]
SRC      <- args[3]
cell_idx <- as.integer(args[4])
indir    <- args[5]
outdir   <- args[6]
subn     <- if (length(args) >= 7) as.integer(args[7]) else 2500L

.libPaths(c(LIB))
stopifnot(identical(.libPaths()[1], LIB))
library(fungwasStage2)
source(file.path(SRC, "truth.R"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cells_txt <- paste(readLines(file.path(WS, "cells.json"), warn = FALSE), collapse = "")
blocks <- regmatches(cells_txt, gregexpr("\\{[^{}]*\\}", cells_txt))[[1]]
getf <- function(b, k) {
  m <- regmatches(b, regexpr(paste0("\"", k, "\"[[:space:]]*:[[:space:]]*[^,}]+"), b))
  v <- sub(paste0("\"", k, "\"[[:space:]]*:[[:space:]]*"), "", m)
  gsub("\"", "", trimws(v))
}
b <- blocks[cell_idx + 1L]
cell <- list(name = getf(b, "name"), arm = getf(b, "arm"), dgp = getf(b, "dgp"),
             maf = as.numeric(getf(b, "maf")), shape = getf(b, "shape"),
             a = as.numeric(getf(b, "a")), c = as.numeric(getf(b, "c")),
             N = as.integer(getf(b, "N")), reps = as.integer(getf(b, "reps")),
             n_blocks = as.integer(getf(b, "n_blocks")))
N <- cell$N; NB <- cell$n_blocks
cat("CELL", cell_idx, cell$name, "N =", N, "B =", NB, "\n")

TAUS <- seq(0.02, 0.98, length.out = 45); Tt <- length(TAUS)

## ---- exact population truth -------------------------------------------------
tr <- pop_truth(TAUS, cell$maf, cell$shape, cell$a, cell$c)
z  <- qnorm(TAUS)
qe <- eps_quantile(TAUS, cell$shape)

W_ADD <- matrix(1, Tt, 1, dimnames = list(paste0("tau", TAUS), "beta_add"))
W_LS  <- make_weights_vqtl(TAUS, tr$q_tau, mu = tr$EY, sd = tr$sdY)
W_LSK <- cbind(beta_mu = rep(1, Tt), beta_sigma2 = z / (2 * tr$sdY), beta_kurt = z^2 - 1)
W_TRU <- cbind(beta_mu = rep(1, Tt), beta_scale = qe)
rownames(W_LSK) <- rownames(W_TRU) <- paste0("tau", TAUS)
fits <- list(ADD = W_ADD, LS = W_LS, LSK = W_LSK, TRU = W_TRU)

S0inv <- solve(tr$Sigma0)
truth_tbl <- list()
for (fn in names(fits)) {
  W <- fits[[fn]]
  A <- solve(crossprod(W), t(W))
  th_ols <- drop(A %*% tr$beta)
  th_gls <- drop(solve(t(W) %*% S0inv %*% W, t(W) %*% S0inv %*% tr$beta))
  names(th_ols) <- names(th_gls) <- colnames(W)
  truth_tbl[[fn]] <- list(
    ols = th_ols, gls = th_gls,
    resid_ols = tr$beta - drop(W %*% th_ols),
    resid_gls = tr$beta - drop(W %*% th_gls),
    se_ols_unit = sqrt(diag(A %*% tr$Sigma0 %*% t(A))),
    se_gls_unit = sqrt(diag(solve(t(W) %*% S0inv %*% W))))
}

## ---- per-tau accumulators (Stage 1 calibration) -----------------------------
pt <- list(n = 0L, sz = numeric(Tt), sz2 = numeric(Tt),
           c95 = numeric(Tt), c99 = numeric(Tt),
           sb = numeric(Tt), sb2 = numeric(Tt), sse = numeric(Tt),
           sf = numeric(Tt), sf2 = numeric(Tt), sq = numeric(Tt))
ref_c <- list(); ref_df_med <- list()

acc <- list(); done <- 0L
buf_beta <- list(); buf_d <- list(); buf_s0 <- list()
buf_f <- list(); buf_q <- list(); nbuf <- 0L

flush_buf <- function() {
  if (nbuf == 0L) return(invisible(NULL))
  Bm <- do.call(rbind, buf_beta); dv <- unlist(buf_d); S0l <- do.call(c, buf_s0)
  Fm <- do.call(rbind, buf_f); Qm <- do.call(rbind, buf_q)
  R <- nrow(Bm); names(S0l) <- as.character(seq_len(R))
  cd <- fungwasStage2:::.pooled_cov_data(S0l, chr = seq_len(R), d_resid = dv, n_taus = Tt)
  cov_data <- list(cov_vec = cd$cov_vec, offsets = cd$offsets)
  SE_tau <- vapply(seq_len(R), function(i) sqrt(diag(S0l[[i]]) / dv[i]), numeric(Tt))

  ## ---- per-tau Stage 1 calibration -----------------------------------------
  Zt <- (t(Bm) - tr$beta) / SE_tau                       # Tt x R
  pt$n   <<- pt$n + R
  pt$sz  <<- pt$sz  + rowSums(Zt)
  pt$sz2 <<- pt$sz2 + rowSums(Zt^2)
  pt$c95 <<- pt$c95 + rowSums(abs(Zt) < qnorm(0.975))
  pt$c99 <<- pt$c99 + rowSums(abs(Zt) < qnorm(0.995))
  pt$sb  <<- pt$sb  + colSums(Bm)
  pt$sb2 <<- pt$sb2 + colSums(Bm^2)
  pt$sse <<- pt$sse + rowSums(SE_tau)
  pt$sf  <<- pt$sf  + colSums(Fm)
  pt$sf2 <<- pt$sf2 + colSums(Fm^2)
  pt$sq  <<- pt$sq  + colSums(Qm)

  Qsl <- t(Bm); colnames(Qsl) <- as.character(seq_len(R)); rownames(Qsl) <- paste0("tau", TAUS)
  stage1 <- list(taus = TAUS, q_tau = tr$q_tau, Q_slope = Qsl, SE_tau = SE_tau)
  for (fn in names(fits)) {
    for (est in c("ols", "gls")) {
      r <- param_gwas(stage1, transform = "custom_W",
                      transform_args = list(W = fits[[fn]]),
                      se_mode = "tau_cov", cov_data = cov_data,
                      tau_cov_estimator = est, nu = N,
                      nu_source = "pooled_sigma0_n_samples", q_method = "saddlepoint")
      key <- paste(fn, est, sep = "_"); prev <- acc[[key]]
      acc[[key]] <<- list(theta = c(prev$theta, list(r$params)),
                          se    = c(prev$se,    list(r$SE_params)),
                          Q     = c(prev$Q,     list(r$Q)),
                          Qp    = c(prev$Qp,    list(r$Q_pval)))
    }
    ## block-jackknife effective-df law, harvested once per fit from the
    ## shipped code (theta/SE are identical to the nu = N call; only the
    ## reference distribution differs).
    if (is.null(ref_c[[fn]])) {
      rj <- param_gwas(stage1, transform = "custom_W",
                       transform_args = list(W = fits[[fn]]),
                       se_mode = "tau_cov", cov_data = cov_data,
                       tau_cov_estimator = "ols", nu = NULL,
                       n_blocks = NB, n_eff = dv, q_method = "saddlepoint")
      ref_c[[fn]] <<- rj$param_ref_c
      ref_df_med[[fn]] <<- apply(rj$param_ref_df, 1, median)
    }
  }
  done <<- done + R
  buf_beta <<- list(); buf_d <<- list(); buf_s0 <<- list()
  buf_f <<- list(); buf_q <<- list(); nbuf <<- 0L
  cat("  ...", done, "reps\n"); flush.console()
}

rdb <- function(p, n) readBin(p, "double", n = n, size = 8, endian = "little")
tags <- sub("\\.beta\\.bin$", "",
            sort(list.files(indir, pattern = sprintf("^cell%02d_chunk[0-9]+\\.beta\\.bin$", cell_idx))))
cat("chunks:", length(tags), "\n")
d_all <- numeric(0); mac_all <- numeric(0)
for (tg in tags) {
  txt <- paste(readLines(file.path(indir, paste0(tg, ".json")), warn = FALSE), collapse = "")
  nrep <- as.integer(sub(".*\"reps\": *([0-9]+).*", "\\1", txt))
  bb <- matrix(rdb(file.path(indir, paste0(tg, ".beta.bin")), nrep * Tt),
               nrow = nrep, ncol = Tt, byrow = TRUE)
  dd <- rdb(file.path(indir, paste0(tg, ".dresid.bin")), nrep)
  mm <- rdb(file.path(indir, paste0(tg, ".mac.bin")), nrep)
  ff <- matrix(rdb(file.path(indir, paste0(tg, ".fhat.bin")), nrep * Tt),
               nrow = nrep, ncol = Tt, byrow = TRUE)
  qq <- matrix(rdb(file.path(indir, paste0(tg, ".qtau.bin")), nrep * Tt),
               nrow = nrep, ncol = Tt, byrow = TRUE)
  ss <- rdb(file.path(indir, paste0(tg, ".sigma0.bin")), nrep * Tt * Tt)
  sl <- lapply(seq_len(nrep), function(i)
    matrix(ss[((i - 1L) * Tt * Tt + 1L):(i * Tt * Tt)], nrow = Tt, ncol = Tt, byrow = TRUE))
  rm(ss)
  buf_beta[[length(buf_beta) + 1L]] <- bb; buf_d[[length(buf_d) + 1L]] <- dd
  buf_s0[[length(buf_s0) + 1L]] <- sl
  buf_f[[length(buf_f) + 1L]] <- ff; buf_q[[length(buf_q) + 1L]] <- qq
  d_all <- c(d_all, dd); mac_all <- c(mac_all, mm); nbuf <- nbuf + nrep
  if (nbuf >= subn) flush_buf()
}
flush_buf()

res <- list()
for (key in names(acc))
  res[[key]] <- list(theta = do.call(cbind, acc[[key]]$theta),
                     se    = do.call(cbind, acc[[key]]$se),
                     Q     = unlist(acc[[key]]$Q),
                     Qp    = unlist(acc[[key]]$Qp))

pertau <- data.frame(
  tau = TAUS, nrep = pt$n,
  beta_true = tr$beta, f_true = tr$f, q_true = tr$q_tau,
  mean_z = pt$sz / pt$n,
  sd_z = sqrt(pmax(pt$sz2 / pt$n - (pt$sz / pt$n)^2, 0) * pt$n / (pt$n - 1)),
  cov95 = pt$c95 / pt$n, cov99 = pt$c99 / pt$n,
  mean_beta = pt$sb / pt$n,
  emp_sd_beta = sqrt(pmax(pt$sb2 / pt$n - (pt$sb / pt$n)^2, 0) * pt$n / (pt$n - 1)),
  mean_se_tau = pt$sse / pt$n,
  mean_fhat = pt$sf / pt$n,
  sd_fhat = sqrt(pmax(pt$sf2 / pt$n - (pt$sf / pt$n)^2, 0) * pt$n / (pt$n - 1)),
  mean_qhat = pt$sq / pt$n,
  stringsAsFactors = FALSE)
pertau$se_ratio_tau <- pertau$mean_se_tau / pertau$emp_sd_beta
pertau$fhat_bias <- pertau$mean_fhat / pertau$f_true - 1

saveRDS(list(cell = cell, cell_idx = cell_idx, N = N, n_blocks = NB, reps = done,
             taus = TAUS, truth = tr, W = fits, theta_true = truth_tbl,
             d_resid = d_all, mac = mac_all, res = res,
             pertau = pertau, ref_c = ref_c, ref_df_med = ref_df_med),
        file.path(outdir, sprintf("stage2_cell%02d.rds", cell_idx)))
cat("WROTE stage2_cell", cell_idx, " reps =", done, "\n")
