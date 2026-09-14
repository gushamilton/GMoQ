#!/usr/bin/env Rscript
# Drive the SHIPPED Stage 2 param_gwas over one collated cell,
# for jackknife/pooled x ols/gls, optionally on a subset of taus.
.libPaths(c("/user/work/fh6520/kde_lumpy/rlib", "/user/work/fh6520/r_packages"))
suppressMessages(library(fungwasStage2))
suppressMessages(library(jsonlite))

args <- commandArgs(TRUE)
dir <- args[1]
tag <- if (length(args) >= 2) args[2] else "full"      # output suffix
keepf <- if (length(args) >= 3) args[3] else ""        # json file with 1-based keep idx

meta <- fromJSON(file.path(dir, "meta.json"))
taus0 <- as.numeric(meta$taus); T0 <- length(taus0); nel0 <- T0 * (T0 + 1) / 2
P <- meta$n_snp
W0 <- matrix(as.numeric(unlist(meta$W)), nrow = T0)
colnames(W0) <- c("shift", "scale")

keep <- seq_len(T0)
if (nzchar(keepf)) keep <- as.integer(fromJSON(keepf)$keep)
Tt <- length(keep); nel <- Tt * (Tt + 1) / 2
stopifnot(Tt >= 3)

rd <- function(f, n) { con <- file(file.path(dir, f), "rb"); v <- readBin(con, "numeric", n = n, size = 8); close(con); v }
beta <- matrix(rd("beta.bin", P * T0), nrow = T0)[keep, , drop = FALSE]
se   <- matrix(rd("se.bin",   P * T0), nrow = T0)[keep, , drop = FALSE]
snpm <- matrix(rd("meta_snp.bin", P * 4), nrow = 4)
W <- W0[keep, , drop = FALSE]
taus <- taus0[keep]

# indices, within a row-major packed upper triangle of size T0, of the
# sub-triangle for `keep`
pos <- matrix(0L, T0, T0); k <- 0L
for (i in seq_len(T0)) for (j in i:T0) { k <- k + 1L; pos[i, j] <- k }
sub <- integer(nel); k <- 0L
for (ii in seq_len(Tt)) for (jj in ii:Tt) { k <- k + 1L; sub[k] <- pos[keep[ii], keep[jj]] }

subset_cov <- function(v) {
  m <- matrix(v, nrow = nel0)          # nel0 x P
  as.numeric(m[sub, , drop = FALSE])
}
covv  <- subset_cov(rd("cov.bin",    P * nel0))
poolv <- subset_cov(rd("pooled.bin", P * nel0))
offsets <- seq(0, by = nel * 4, length.out = P)
stage1 <- list(taus = taus, q_tau = rep(0, Tt), Q_slope = beta, SE_tau = se)

cat("tag=", tag, " taus kept:", Tt, "of", T0, "\n")
for (src in c("jackknife", "pooled")) {
  cd <- list(cov_vec = if (src == "jackknife") covv else poolv, offsets = offsets)
  for (est in c("ols", "gls")) {
    a <- list(stage1 = stage1, transform = "custom_W",
              transform_args = list(W = W), se_mode = "tau_cov",
              cov_data = cd, tau_cov_estimator = est, n_threads = 4L,
              unknown_nu = "omit", q_method = "saddlepoint")
    if (src == "pooled") {
      a$nu <- as.numeric(meta$N); a$nu_source <- "pooled_sigma0_n_samples"
    } else {
      a$n_blocks <- as.numeric(meta$blocks); a$n_eff <- as.numeric(snpm[4, ])
    }
    r <- withCallingHandlers(do.call(param_gwas, a),
          message = function(m) { cat("[msg]", conditionMessage(m)); invokeRestart("muffleMessage") })
    key <- paste0(src, "_", est)
    nz <- function(x) if (is.null(x)) rep(NA_real_, P) else as.numeric(x)
    M8 <- cbind(as.numeric(r$params[1, ]), as.numeric(r$params[2, ]),
                as.numeric(r$SE_params[1, ]), as.numeric(r$SE_params[2, ]),
                nz(if (is.null(r$P_params)) NULL else r$P_params[1, ]),
                nz(if (is.null(r$P_params)) NULL else r$P_params[2, ]),
                nz(r$param_ref_df[1, ]), nz(r$param_ref_df[2, ]),
                nz(r$Q), nz(r$Q_pval))
    con <- file(file.path(dir, paste0("s2_", tag, "_", key, ".bin")), "wb")
    writeBin(as.numeric(t(M8)), con, size = 8); close(con)
    cat(key, "ref=", r$param_ref_dist, "\n")
  }
}

# covariance conditioning diagnostics, first 200 SNPs
Sig <- matrix(0, Tt, Tt); iu <- which(upper.tri(Sig, diag = TRUE))
iu <- iu[order(row(Sig)[iu], col(Sig)[iu])]
cond <- dup <- mineig <- numeric(0)
for (j in seq_len(min(200L, P))) {
  v <- covv[((j - 1) * nel + 1):(j * nel)]
  S <- matrix(0, Tt, Tt); S[iu] <- v; S[lower.tri(S)] <- t(S)[lower.tri(S)]
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  cond <- c(cond, max(ev) / max(min(ev), .Machine$double.eps))
  mineig <- c(mineig, min(ev) / max(ev))
  R <- suppressWarnings(cov2cor(S)); diag(R) <- 0
  dup <- c(dup, sum(abs(R) > 1 - 1e-6) / 2)
}
write_json(list(tag = tag, taus_kept = keep, n_taus = Tt,
                cond_median = median(cond), mineig_ratio_median = median(mineig),
                dup_pairs_median = median(dup)),
           file.path(dir, paste0("cond_", tag, ".json")), auto_unbox = TRUE)
cat("cond median", signif(median(cond), 3), " dup pairs", median(dup), "\n")
