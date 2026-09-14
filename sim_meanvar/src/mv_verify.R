## Independent, hand-rolled verification of the Stage 2 call used by the
## mean-variance Monte Carlo.  Reproduces param_gwas()'s OLS and GLS output from
## first principles on one replicate and reports the max errors.
.libPaths(c("/user/work/fh6520/gls_coverage/lib", "/user/work/fh6520/r_packages"))
stopifnot("fungwasStage2" %in% rownames(installed.packages()))
suppressPackageStartupMessages(library(fungwasStage2))

a <- commandArgs(trailingOnly = TRUE)
indir <- a[1]; tag <- a[2]
TAUS <- seq(0.02, 0.98, length.out = 45); Tt <- 45L

meta <- paste(readLines(file.path(indir, paste0(tag, ".json")), warn = FALSE), collapse = "")
gi <- function(k) as.integer(sub(paste0(".*\"", k, "\": *([0-9]+).*"), "\\1", meta))
R <- gi("reps"); M <- gi("M"); N <- gi("N")
rd <- function(k, dims) {
  v <- readBin(file.path(indir, paste0(tag, ".", k, ".bin")), "double",
               n = prod(dims), size = 8, endian = "little")
  array(aperm(array(v, rev(dims))), dims)
}
r <- 1L
B  <- rd("beta", c(R, M, Tt)); D <- rd("dresid", c(R, M))
S0 <- rd("sigma0", c(R, Tt, Tt)); QT <- rd("qtau", c(R, Tt)); MOM <- rd("mom", c(R, 5L))
S0r <- matrix(S0[r, , ], Tt, Tt); d <- D[r, ]; Bm <- matrix(B[r, , ], M, Tt)

cat("Sigma0 symmetric:", max(abs(S0r - t(S0r))) < 1e-14, "\n")
Sl <- list(S0r); names(Sl) <- "1"
cd <- fungwasStage2:::.pooled_cov_data(Sl, chr = rep(1L, M), d_resid = d, n_taus = Tt)
m0 <- matrix(0, Tt, Tt); iu <- which(upper.tri(m0, diag = TRUE))
iu <- iu[order(row(m0)[iu], col(m0)[iu])]; ne <- Tt * (Tt + 1) / 2
rt <- sapply(c(1L, 57L, M), function(i) {
  v <- cd$cov_vec[((i - 1) * ne + 1):(i * ne)]
  S <- matrix(0, Tt, Tt); S[iu] <- v; S <- S + t(S) - diag(diag(S))
  max(abs(S - S0r / d[i]))
})
cat("cov_vec round-trip max abs err:", format(rt, digits = 3), "\n")

W <- make_weights_vqtl(TAUS, QT[r, ], mu = MOM[r, 1], sd = MOM[r, 2])
cat("colnames(W):", colnames(W), "\n")
A <- solve(crossprod(W), t(W))
SE_tau <- outer(sqrt(diag(S0r)), 1 / sqrt(d))
Qsl <- t(Bm); colnames(Qsl) <- paste0("v", seq_len(M)); rownames(Qsl) <- paste0("tau", TAUS)
dimnames(SE_tau) <- dimnames(Qsl)
s1 <- list(taus = TAUS, q_tau = QT[r, ], Q_slope = Qsl, SE_tau = SE_tau)
cdl <- list(cov_vec = cd$cov_vec, offsets = cd$offsets)

ro <- param_gwas(s1, "custom_W", list(W = W), se_mode = "tau_cov", cov_data = cdl,
                 tau_cov_estimator = "ols", nu = N,
                 nu_source = "pooled_sigma0_n_samples", return_cov = "upper")
rg <- param_gwas(s1, "custom_W", list(W = W), se_mode = "tau_cov", cov_data = cdl,
                 tau_cov_estimator = "gls", nu = N,
                 nu_source = "pooled_sigma0_n_samples", return_cov = "upper")

man_o_th <- A %*% Qsl
man_o_se <- vapply(seq_len(M), function(i) sqrt(diag(A %*% (S0r / d[i]) %*% t(A))), numeric(2))
man_o_cov <- vapply(seq_len(M), function(i) (A %*% (S0r / d[i]) %*% t(A))[1, 2], 0.0)
man_g <- vapply(seq_len(M), function(i) {
  Si <- solve(S0r / d[i]); Mx <- solve(t(W) %*% Si %*% W)
  V <- Mx
  c(drop(Mx %*% t(W) %*% Si %*% Bm[i, ]), sqrt(diag(V)), V[1, 2])
}, numeric(5))

cat("params rownames:", rownames(ro$params), "\n")
cat("param_cov_upper colnames:", colnames(ro$param_cov_upper), "\n")
cat("OLS theta  max abs err :", max(abs(ro$params - man_o_th)), "\n")
cat("OLS SE     max rel err :", max(abs(ro$SE_params / man_o_se - 1)), "\n")
cat("OLS cov_mu_sigma2 max rel err:",
    max(abs(ro$param_cov_upper[, "cov_beta_mu_beta_sigma2"] / man_o_cov - 1)), "\n")
cat("OLS cov diag vs SE^2 max rel err:",
    max(abs(ro$param_cov_upper[, "cov_beta_mu_beta_mu"] / man_o_se[1, ]^2 - 1)),
    max(abs(ro$param_cov_upper[, "cov_beta_sigma2_beta_sigma2"] / man_o_se[2, ]^2 - 1)), "\n")
cat("GLS theta  max rel err :", max(abs(rg$params / man_g[1:2, ] - 1)), "\n")
cat("GLS SE     max rel err :", max(abs(rg$SE_params / man_g[3:4, ] - 1)), "\n")
cat("GLS cov_mu_sigma2 max rel err:",
    max(abs(rg$param_cov_upper[, "cov_beta_mu_beta_sigma2"] / man_g[5, ] - 1)), "\n")
cat("param_ref_dist:", ro$param_ref_dist, " df[1,1]:", ro$param_ref_df[1, 1], "\n")
