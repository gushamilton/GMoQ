args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]; outdir <- args[2]
.libPaths(c(LIB)); stopifnot(identical(.libPaths()[1], LIB))
x <- readRDS(file.path(outdir, "stage2_cell25.rds"))
for (key in names(x$res)) {
  r <- x$res[[key]]
  cat(sprintf("%-10s Q: %-10s n=%-7s   Qp: %-10s n=%s\n", key,
              class(r$Q)[1], length(r$Q), class(r$Qp)[1], length(r$Qp)))
}
cat("\n-- direct param_gwas probe, same inputs, ols vs gls --\n")
library(fungwasStage2)
set.seed(1)
Tt <- 45; TAUS <- seq(0.02, 0.98, length.out = Tt)
S0 <- x$truth$Sigma0; d <- median(x$d_resid)
R <- 5
S0l <- setNames(replicate(R, S0, simplify = FALSE), as.character(seq_len(R)))
cd <- fungwasStage2:::.pooled_cov_data(S0l, chr = seq_len(R),
                                       d_resid = rep(d, R), n_taus = Tt)
B <- replicate(R, x$truth$beta + drop(chol(S0 / d) %*% rnorm(Tt)))
rownames(B) <- paste0("tau", TAUS); colnames(B) <- as.character(seq_len(R))
st <- list(taus = TAUS, q_tau = x$truth$q_tau, Q_slope = B,
           SE_tau = matrix(sqrt(diag(S0) / d), Tt, R))
for (est in c("ols", "gls")) {
  r <- param_gwas(st, transform = "custom_W",
                  transform_args = list(W = x$W$LS), se_mode = "tau_cov",
                  cov_data = list(cov_vec = cd$cov_vec, offsets = cd$offsets),
                  tau_cov_estimator = est, nu = x$N,
                  nu_source = "pooled_sigma0_n_samples", q_method = "saddlepoint")
  cat(sprintf("est=%-4s  names: %s\n", est, paste(names(r), collapse = ", ")))
  cat(sprintf("          Q = %s\n", paste(signif(r$Q, 4), collapse = " ")))
  cat(sprintf("          Q_pval = %s\n",
              if (is.null(r$Q_pval)) "NULL" else paste(signif(r$Q_pval, 3), collapse = " ")))
  cat(sprintf("          df = %s\n",
              if (is.null(r$df)) "NULL" else paste(r$df, collapse = " ")))
}
