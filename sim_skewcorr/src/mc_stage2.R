## Stage 2 driver for the skew-correction validation.
##
## Adapted from /user/work/fh6520/gls_coverage/src/stage2_driver.R: same chunked
## read of the Stage 1 binaries, same .pooled_cov_data() call, same SHIPPED
## param_gwas() for ols and gls.  Added: return_cov = "upper" so the within-
## variant covariance is available, and the INSTALLED correct_dispersion_skew()
## applied per replicate.
WS <- "/user/work/fh6520/skewcorr"
.libPaths(file.path(WS, "rlib"))
stopifnot(identical(.libPaths()[1], file.path(WS, "rlib")))
library(fungwasStage2)

args     <- commandArgs(trailingOnly = TRUE)
cell_idx <- as.integer(args[1])
indir    <- args[2]
outdir   <- args[3]
subn     <- if (length(args) >= 4) as.integer(args[4]) else 2000L
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- cell definition (minimal JSON parse; no jsonlite under --vanilla) ------
cells_txt <- paste(readLines(file.path(WS, "cells.json"), warn = FALSE), collapse = "")
blocks <- regmatches(cells_txt, gregexpr("\\{[^{}]*\\}", cells_txt))[[1]]
getf <- function(b, k) {
  m <- regmatches(b, regexpr(paste0("\"", k, "\"[[:space:]]*:[[:space:]]*[^,}]+"), b))
  gsub("\"", "", trimws(sub(paste0("\"", k, "\"[[:space:]]*:[[:space:]]*"), "", m)))
}
b <- blocks[cell_idx + 1L]
cell <- list(name = getf(b, "name"), maf = as.numeric(getf(b, "maf")),
             shape = getf(b, "shape"), delta = as.numeric(getf(b, "delta")),
             gamma = as.numeric(getf(b, "gamma")), N = as.integer(getf(b, "N")))
N <- cell$N
cat("CELL", cell_idx, cell$name, "N =", N, "\n")

TAUS <- seq(0.02, 0.98, length.out = 45); Tt <- length(TAUS)

## ---- the baseline W, built from the NULL population law --------------------
## For a normal baseline this is exactly make_weights_vqtl(); for skew/heavy the
## null quantiles come from the exact engine's dump (see mc_truth.py).
ln <- readLines(file.path(WS, "truth", sprintf("truth_cell%02d.txt", cell_idx)))
mv <- as.numeric(strsplit(trimws(ln[1]), "[[:space:]]+")[[1]])
tr <- list(m0 = mv[1], v0 = mv[2],
           q0 = as.numeric(strsplit(trimws(ln[2]), "[[:space:]]+")[[1]]))
stopifnot(length(tr$q0) == Tt)
W <- make_weights_vqtl(TAUS, q_tau = tr$q0, mu = tr$m0, sd = sqrt(tr$v0))
## make_weights_vqtl's sigma2 column is z_tau/(2 sd), which equals
## (q0 - m0)/(2 v0) exactly for a normal baseline. For skew/heavy baselines the
## SHIPPED weight is still z_tau/(2 sd) -- that is what a user would build,
## knowing only mu and sd -- so leave it as the package produced it.

acc <- list(); done <- 0L
buf_beta <- list(); buf_d <- list(); buf_s0 <- list(); nbuf <- 0L

flush_buf <- function() {
  if (nbuf == 0L) return(invisible(NULL))
  Bm <- do.call(rbind, buf_beta); dv <- unlist(buf_d); S0l <- do.call(c, buf_s0)
  R <- nrow(Bm); names(S0l) <- as.character(seq_len(R))
  cd <- fungwasStage2:::.pooled_cov_data(S0l, chr = seq_len(R), d_resid = dv, n_taus = Tt)
  cov_data <- list(cov_vec = cd$cov_vec, offsets = cd$offsets)
  SE_tau <- vapply(seq_len(R), function(i) sqrt(diag(S0l[[i]]) / dv[i]), numeric(Tt))
  Qsl <- t(Bm); colnames(Qsl) <- as.character(seq_len(R)); rownames(Qsl) <- paste0("tau", TAUS)
  stage1 <- list(taus = TAUS, q_tau = tr$q0, Q_slope = Qsl, SE_tau = SE_tau)
  for (est in c("ols", "gls")) {
    r <- param_gwas(stage1, transform = "custom_W", transform_args = list(W = W),
                    se_mode = "tau_cov", cov_data = cov_data,
                    tau_cov_estimator = est, return_cov = "upper",
                    nu = N, nu_source = "pooled_sigma0_n_samples",
                    q_method = "saddlepoint")
    cn <- colnames(r$param_cov_upper)
    ci <- which(cn %in% c("cov_beta_mu_beta_sigma2", "cov_beta_sigma2_beta_mu"))
    cl <- as.numeric(r$param_cov_upper[, ci[1]])
    cc <- correct_dispersion_skew(
      beta_mu = r$params["beta_mu", ], beta_sigma2 = r$params["beta_sigma2", ],
      se_mu = r$SE_params["beta_mu", ], se_sigma2 = r$SE_params["beta_sigma2", ],
      cov_mu_sigma2 = cl, eaf = cell$maf)
    prev <- acc[[est]]
    acc[[est]] <<- list(
      bmu  = c(prev$bmu,  r$params["beta_mu", ]),
      bs2  = c(prev$bs2,  r$params["beta_sigma2", ]),
      smu  = c(prev$smu,  r$SE_params["beta_mu", ]),
      ss2  = c(prev$ss2,  r$SE_params["beta_sigma2", ]),
      cov  = c(prev$cov,  cl),
      cb   = c(prev$cb,   cc$beta_sigma2_skewcorr),
      cs   = c(prev$cs,   cc$se_sigma2_skewcorr))
  }
  done <<- done + R
  buf_beta <<- list(); buf_d <<- list(); buf_s0 <<- list(); nbuf <<- 0L
  cat("  ...", done, "reps\n"); flush.console()
}

tags <- sub("\\.beta\\.bin$", "",
            sort(list.files(indir, pattern = sprintf("^cell%02d_chunk[0-9]+\\.beta\\.bin$", cell_idx))))
cat("chunks:", length(tags), "\n")
for (tg in tags) {
  txt <- paste(readLines(file.path(indir, paste0(tg, ".json")), warn = FALSE), collapse = "")
  nrep <- as.integer(sub(".*\"reps\": *([0-9]+).*", "\\1", txt))
  bb <- matrix(readBin(file.path(indir, paste0(tg, ".beta.bin")), "double", n = nrep*Tt,
                       size = 8, endian = "little"), nrow = nrep, ncol = Tt, byrow = TRUE)
  dd <- readBin(file.path(indir, paste0(tg, ".dresid.bin")), "double", n = nrep, size = 8, endian = "little")
  ss <- readBin(file.path(indir, paste0(tg, ".sigma0.bin")), "double", n = nrep*Tt*Tt, size = 8, endian = "little")
  sl <- lapply(seq_len(nrep), function(i)
    matrix(ss[((i-1L)*Tt*Tt + 1L):(i*Tt*Tt)], nrow = Tt, ncol = Tt, byrow = TRUE))
  rm(ss)
  buf_beta[[length(buf_beta)+1L]] <- bb; buf_d[[length(buf_d)+1L]] <- dd
  buf_s0[[length(buf_s0)+1L]] <- sl; nbuf <- nbuf + nrep
  if (nbuf >= subn) flush_buf()
}
flush_buf()

saveRDS(list(cell = cell, cell_idx = cell_idx, reps = done, taus = TAUS,
             truth = tr, W = W, res = acc),
        file.path(outdir, sprintf("mc_cell%02d.rds", cell_idx)))
cat("WROTE mc_cell", cell_idx, " reps =", done, "\n")
