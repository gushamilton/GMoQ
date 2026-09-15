## Stage 2 driver for the bounded biobank covariance audit.
##
## Runs the shipped param_gwas() over OLS with cov_source in {pooled, jackknife}
## on the SAME Stage 1 replicates, so every comparison is paired and the
## Monte Carlo error on a difference is far smaller than on either rate.
##
## Nothing here reimplements an estimator.  The pooled covariance is assembled
## by the package's own .pooled_cov_data(); the jackknife covariance is the
## float32 upper triangle Stage 1 produced, handed over in exactly the layout
## and precision .read_tau_cov_data() would have produced from a .cov.gz.
##
## Degrees of freedom follow the shipped conventions for each source:
##   pooled    -> nu = n_samples, nu_source = "pooled_sigma0_n_samples"
##                (what param_gwas_from_file(cov_source="pooled") does)
##   jackknife -> nu = NULL, n_blocks = 1000, n_eff = d_resid per SNP
##                (the per-SNP Satterthwaite law, rung 2 of .param_reference)
args     <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3L)
cell_idx <- as.integer(args[1])
indir    <- args[2]
outdir   <- args[3]
subn     <- if (length(args) >= 4) as.integer(args[4]) else 2500L
WS       <- if (length(args) >= 5) args[5] else "/user/work/fh6520/covsource_run"
LIB      <- if (length(args) >= 6) args[6] else "/user/work/fh6520/gmoq_finish_20260915/inference/lib"
.libPaths(c(LIB, .libPaths()))
stopifnot(identical(.libPaths()[1], LIB))
library(fungwasStage2)
source(file.path(WS, "src", "truth_cs.R"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- cell (minimal JSON parse; no jsonlite under --vanilla) -----------------
cells_txt <- paste(readLines(file.path(WS, "cells.json"), warn = FALSE), collapse = "")
blocks <- regmatches(cells_txt, gregexpr("\\{[^{}]*\\}", cells_txt))[[1]]
getf <- function(b, k) {
  m <- regmatches(b, regexpr(paste0("\"", k, "\"[[:space:]]*:[[:space:]]*[^,}]+"), b))
  gsub("\"", "", trimws(sub(paste0("\"", k, "\"[[:space:]]*:[[:space:]]*"), "", m)))
}
b <- blocks[cell_idx + 1L]
cell <- list(name = getf(b, "name"), dgp = getf(b, "dgp"),
             maf = as.numeric(getf(b, "maf")), shape = getf(b, "shape"),
             a = as.numeric(getf(b, "a")), c = as.numeric(getf(b, "c")),
             N = as.integer(getf(b, "N")), reps = as.integer(getf(b, "reps")),
             Ztarget = suppressWarnings(as.numeric(getf(b, "Ztarget"))))
N <- cell$N
cat("CELL", cell_idx, cell$name, "N =", N, "\n")

TAUS <- seq(0.02, 0.98, length.out = 45); Tt <- length(TAUS)
NUP  <- Tt * (Tt + 1L) / 2L
N_BLOCKS <- 1000L

## ---- exact population truth -------------------------------------------------
tr   <- pop_truth(TAUS, cell$maf, cell$shape, cell$a, cell$c)
fits <- cs_weights(TAUS, tr, make_weights_vqtl)
d_pop <- N * tr$Vg
truth_tbl <- lapply(fits, function(W) cs_theta_truth(W, tr, d_pop))

## row-major upper-triangle index, the order Stage 1 wrote and C++ reads
iu <- which(upper.tri(matrix(0, Tt, Tt), diag = TRUE))
ord <- order(row(matrix(0, Tt, Tt))[iu], col(matrix(0, Tt, Tt))[iu])
iu <- iu[ord]
upper_to_full <- function(v) { M <- matrix(0, Tt, Tt); M[iu] <- v; M[lower.tri(M)] <- t(M)[lower.tri(M)]; M }

## ---- accumulate -------------------------------------------------------------
acc <- new.env(parent = emptyenv()); done <- 0L
buf <- new.env(parent = emptyenv())
reset_buf <- function() { buf$beta <- list(); buf$se <- list(); buf$d <- list()
                          buf$s0 <- list(); buf$jk <- list(); buf$n <- 0L }
reset_buf()

push <- function(key, r) {
  prev <- acc[[key]]
  acc[[key]] <- list(theta = c(prev$theta, list(r$params)),
                     se    = c(prev$se,    list(r$SE_params)),
                     pval  = c(prev$pval,  list(r$P_params)),
                     df    = c(prev$df,    list(r$param_ref_df)),
                     Q     = c(prev$Q,     list(r$Q)),
                     Qp    = c(prev$Qp,    list(if (is.null(r$Q_pval)) rep(NA_real_, length(r$Q)) else r$Q_pval)))
}

flush_buf <- function() {
  if (buf$n == 0L) return(invisible(NULL))
  Bm  <- do.call(rbind, buf$beta)
  SEm <- do.call(rbind, buf$se)
  dv  <- unlist(buf$d)
  S0m <- do.call(rbind, buf$s0)          # R x NUP, float64 upper triangle
  JKm <- do.call(rbind, buf$jk)          # R x NUP, float32-valued upper triangle
  R   <- nrow(Bm)

  ## --- pooled: assembled by the package itself -------------------------------
  S0l <- lapply(seq_len(R), function(i) upper_to_full(S0m[i, ]))
  names(S0l) <- as.character(seq_len(R))
  cd <- fungwasStage2:::.pooled_cov_data(S0l, chr = seq_len(R), d_resid = dv, n_taus = Tt)
  cov_pooled <- list(cov_vec = cd$cov_vec, offsets = cd$offsets)
  rm(S0l); gc(FALSE)

  ## --- jackknife: the Stage 1 payload, same layout the .cov.gz reader yields --
  cov_jk <- list(cov_vec = as.numeric(t(JKm)),
                 offsets = seq(0, by = NUP * 4, length.out = R))

  Qsl <- t(Bm); colnames(Qsl) <- as.character(seq_len(R)); rownames(Qsl) <- paste0("tau", TAUS)
  SEt <- t(SEm); colnames(SEt) <- colnames(Qsl); rownames(SEt) <- rownames(Qsl)
  stage1 <- list(taus = TAUS, q_tau = tr$q_tau, Q_slope = Qsl, SE_tau = SEt)

  ## Bounded audit configuration: the shipped OLS transform with pooled and
  ## jackknife covariance, paired on the unchanged Stage 1 payload.
  for (fn in "LS") for (est in "ols") for (cs in c("pooled", "jackknife")) {
    cdta <- if (cs == "pooled") cov_pooled else cov_jk
    r <- if (cs == "pooled") {
      param_gwas(stage1, transform = "custom_W", transform_args = list(W = fits[[fn]]),
                 se_mode = "tau_cov", cov_data = cdta, tau_cov_estimator = est,
                 nu = N, nu_source = "pooled_sigma0_n_samples", q_method = "saddlepoint")
    } else {
      param_gwas(stage1, transform = "custom_W", transform_args = list(W = fits[[fn]]),
                 se_mode = "tau_cov", cov_data = cdta, tau_cov_estimator = est,
                 nu = NULL, n_blocks = N_BLOCKS, n_eff = dv, q_method = "saddlepoint")
    }
    push(paste(fn, est, cs, sep = "|"), r)
  }
  done <<- done + R
  reset_buf()
  cat("  ...", done, "reps\n"); flush.console()
}

tags <- sub("\\.beta\\.bin$", "",
            sort(list.files(indir, pattern = sprintf("^cell%03d_chunk[0-9]+\\.beta\\.bin$", cell_idx))))
cat("chunks:", length(tags), "\n")
stopifnot(length(tags) > 0)
d_all <- numeric(0); mac_all <- numeric(0)
rdm <- function(p, n, ncol, size = 8) matrix(readBin(file.path(indir, p), "double", n = n * ncol,
                                                     size = size, endian = "little"),
                                             nrow = n, ncol = ncol, byrow = TRUE)
for (tg in tags) {
  txt <- paste(readLines(file.path(indir, paste0(tg, ".json")), warn = FALSE), collapse = "")
  nrep <- as.integer(sub(".*\"reps\": *([0-9]+).*", "\\1", txt))
  buf$beta[[length(buf$beta)+1L]] <- rdm(paste0(tg, ".beta.bin"),   nrep, Tt)
  buf$se[[length(buf$se)+1L]]     <- rdm(paste0(tg, ".setau.bin"),  nrep, Tt)
  buf$s0[[length(buf$s0)+1L]]     <- rdm(paste0(tg, ".sigma0u.bin"), nrep, NUP)
  buf$jk[[length(buf$jk)+1L]]     <- rdm(paste0(tg, ".jkcovu.bin"), nrep, NUP, size = 4)
  dd <- readBin(file.path(indir, paste0(tg, ".dresid.bin")), "double", n = nrep, size = 8, endian = "little")
  mm <- readBin(file.path(indir, paste0(tg, ".mac.bin")), "double", n = nrep, size = 8, endian = "little")
  buf$d[[length(buf$d)+1L]] <- dd
  d_all <- c(d_all, dd); mac_all <- c(mac_all, mm); buf$n <- buf$n + nrep
  if (buf$n >= subn) flush_buf()
}
flush_buf()

res <- list()
for (key in ls(acc)) {
  A <- acc[[key]]
  res[[key]] <- list(theta = do.call(cbind, A$theta),
                     se    = do.call(cbind, A$se),
                     pval  = do.call(cbind, A$pval),
                     df    = do.call(cbind, A$df),
                     Q     = unlist(A$Q), Qp = unlist(A$Qp))
}
saveRDS(list(cell = cell, cell_idx = cell_idx, N = N, reps = done, taus = TAUS,
             n_blocks = N_BLOCKS, truth = tr, W = fits, theta_true = truth_tbl,
             d_pop = d_pop, d_resid = d_all, mac = mac_all, res = res),
        file.path(outdir, sprintf("stage2_cell%03d.rds", cell_idx)))
cat("WROTE stage2_cell", cell_idx, " reps =", done, "\n")
