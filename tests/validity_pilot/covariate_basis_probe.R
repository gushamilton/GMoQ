#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: covariate_basis_probe.R <probe_dir> <out.csv>")
dir <- normalizePath(args[1]); out <- args[2]
test_lib <- Sys.getenv("GMOQ_TEST_LIB"); dep_lib <- Sys.getenv("GMOQ_DEP_LIB")
.libPaths(unique(c(test_lib, dep_lib[nzchar(dep_lib)], .libPaths())))
suppressPackageStartupMessages(library(fungwasStage2))
pcd <- get(".pooled_cov_data", envir = asNamespace("fungwasStage2"))
read_arr <- function(f, dims) array(readBin(f, "double", n = prod(dims), size = 8, endian = "little"), dim = dims)
tr_files <- list.files(file.path(dir, "truth"), full.names = TRUE, pattern = "\\.json$")
rows <- list(); ir <- 0L
add <- function(x) { ir <<- ir + 1L; rows[[ir]] <<- x }
for (tf in tr_files) {
  tr <- jsonlite::fromJSON(tf, simplifyVector = FALSE); cell <- tr$cell$name; kind <- tr$kind
  for (path in c("current", "y_first")) {
    pfx <- file.path(dir, "stage1", paste0(cell, "__", kind, "__", path))
    jm <- jsonlite::fromJSON(paste0(pfx, ".json"), simplifyVector = FALSE)
    R <- as.integer(jm$reps); T <- as.integer(jm$T); M <- as.integer(jm$M); B <- as.integer(jm$blocks); N <- as.integer(jm$N); NE <- as.integer(jm$ne)
    beta <- read_arr(paste0(pfx, ".beta.bin"), c(T, M, R)); se <- read_arr(paste0(pfx, ".se.bin"), c(T, M, R))
    cov <- read_arr(paste0(pfx, ".cov.bin"), c(NE, M, R)); dres <- read_arr(paste0(pfx, ".dresid.bin"), c(M, R))
    s0 <- read_arr(paste0(pfx, ".sigma0.bin"), c(T, T, R)); qt <- read_arr(paste0(pfx, ".q_tau.bin"), c(T, R))
    ym <- read_arr(paste0(pfx, ".y_mean.bin"), R); ys <- read_arr(paste0(pfx, ".y_sd.bin"), R)
    taus <- as.numeric(unlist(tr$taus)); W0 <- make_weights_vqtl(taus, qnorm(taus), 0, 1)
    for (r in seq_len(R)) {
      Qs <- beta[, , r, drop = FALSE]; dim(Qs) <- c(T, M); Ss <- se[, , r, drop = FALSE]; dim(Ss) <- c(T, M)
      colnames(Qs) <- colnames(Ss) <- paste0("v", seq_len(M)); st <- list(taus=taus, q_tau=qt[,r], Q_slope=Qs, SE_tau=Ss)
      cvj <- list(cov_vec=as.numeric(cov[,,r]), offsets=seq(0, by=NE*4, length.out=M)); S <- s0[,,r,drop=FALSE]; dim(S) <- c(T,T)
      cvp <- pcd(list(ALL=S), rep("1", M), dres[,r], T); Ws <- list(oracle=W0, estimated=make_weights_vqtl(taus, qt[,r], ym[r], ys[r]))
      for (src in c("jackknife", "pooled")) for (wn in names(Ws)) {
        cvd <- if (src == "jackknife") cvj else list(cov_vec=cvp$cov_vec, offsets=cvp$offsets)
        fit <- tryCatch(suppressWarnings(suppressMessages(param_gwas(st, transform="custom_W", transform_args=list(W=Ws[[wn]]), se_mode="tau_cov", tau_cov_estimator="ols", cov_data=cvd, q_method="saddlepoint", n_blocks=if(src=="jackknife") B else NULL, n_eff=if(src=="jackknife") dres[,r] else NULL, nu=if(src=="pooled") N else NULL, nu_source=if(src=="pooled") "pooled_sigma0_n_samples" else NULL))), error=function(e) list(error=conditionMessage(e)))
        if (!is.null(fit$error)) next
        for (k in seq_len(nrow(fit$params))) { nm <- rownames(fit$params)[k]; est <- fit$params[k,1]; se0 <- fit$SE_params[k,1]; pval <- fit$P_params[k,1]
          add(data.frame(cell=cell, basis=kind, path=path, source=src, W=wn, parameter=nm, rep=r, estimate=est, se=se0, p=pval, z=est/se0, Q=fit$Q[1], Q_p=if(is.null(fit$Q_pval)) NA else fit$Q_pval[1]))
        }
      }
    }
  }
}
d <- do.call(rbind, rows); write.csv(d, out, row.names=FALSE, na="")
keys <- unique(d[c("cell","basis","path","source","W","parameter")]); sm <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) { z <- d[d$cell==keys$cell[i] & d$basis==keys$basis[i] & d$path==keys$path[i] & d$source==keys$source[i] & d$W==keys$W[i] & d$parameter==keys$parameter[i],]
 sm[[i]] <- data.frame(keys[i,], datasets=nrow(z), mean_z=mean(z$z), mcse_z=sd(z$z)/sqrt(nrow(z)), rejection=mean(z$p<.05), mcse_rejection=sd(z$p<.05)/sqrt(nrow(z)), bias=mean(z$estimate), mcse_bias=sd(z$estimate)/sqrt(nrow(z)), mean_Q=mean(z$Q), mcse_Q=sd(z$Q)/sqrt(nrow(z)), Q_rejection=mean(z$Q_p<.05,na.rm=TRUE)) }
write.csv(do.call(rbind, sm), sub("\\.csv$", "_summary.csv", out), row.names=FALSE, na="")
cat("wrote", nrow(d), "rows\n")
