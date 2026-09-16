#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE)
S1 <- a[[1]]; SG0 <- a[[2]]; WTS <- a[[3]]; S2 <- a[[4]]; OUT <- a[[5]]

wl <- readRDS(WTS); W <- wl$W
s1 <- fread(S1); bc <- grep("^beta_tau_", names(s1), value = TRUE)
B <- as.matrix(s1[, bc, with = FALSE]); D <- s1$d_resid
sg <- fread(SG0); Sig0 <- as.matrix(sg[, grep("^s_", names(sg), value = TRUE), with = FALSE])
stopifnot(nrow(Sig0) == 45, ncol(Sig0) == 45)
s2 <- fread(S2)
stopifnot(identical(s1$snp_id, s2$snp_id) || all(s1$snp_id == s2$snp_id))

A  <- solve(crossprod(W), t(W))            # 3 x 45
TH <- B %*% t(A)                           # P x 3
V1 <- diag(A %*% Sig0 %*% t(A))            # per-parameter unit variance
SE <- sqrt(outer(1/D, V1))
colnames(TH) <- colnames(SE) <- colnames(W)
NS <- s2$gamma_df[1]
cat("N (df) =", NS, "  n snps =", nrow(s1), "\n\n")
for (p in colnames(W)) {
  th <- TH[, p]; se <- SE[, p]
  tt <- th/se; pv <- 2 * pt(-abs(tt), df = NS)
  cat(sprintf("%-11s theta max|rel err| %.3e   se max|rel err| %.3e   p max|rel err| %.3e\n",
      p, max(abs(th - s2[[p]])/pmax(abs(s2[[p]]),1e-300)),
      max(abs(se - s2[[paste0(p,"_se")]])/s2[[paste0(p,"_se")]]),
      max(abs(pv - s2[[paste0(p,"_p")]])/s2[[paste0(p,"_p")]])))
}
# unit sd of the contrast in RIF units: sqrt(a' Sigma0 a) == sd_u
cat("\nsqrt(diag(A Sig0 A')):", sprintf("%.6g", sqrt(V1)), "\n")
saveRDS(list(A = A, W = W, Sig0 = Sig0, sd_u = sqrt(V1), N = NS), file.path(OUT, "contrasts.rds"))
write.table(A, file.path(OUT, "A_contrasts.tsv"), sep = "\t", quote = FALSE, col.names = NA)
