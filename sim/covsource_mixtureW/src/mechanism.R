## MECHANISM + the one thing that can break it.
##
## Claim to be tested:
##   (1) Sigma_true - Sigma_pooled = D E D, E[s,t] = Delta(min) - t_s D_t - t_t D_s.
##   (2) On its own diagonal, E[t,t] = Delta(t) * (1 - 2t).  This vanishes at
##       tau = 0.5 FOR ANY Delta.
##   (3) For a symmetric bimodal phenotype the density trough is at tau = 0.5,
##       Sigma = D B D collapses onto the trough spike (projeff: 98.5% of the
##       spectrum at S=6), and W_gamma IS that spike.
##   (4) Therefore the leading eigenvalue of Sigma is IDENTICAL under pooled and
##       true covariance, and the column that loads on it -- gamma -- is exactly
##       calibrated.  beta_left/beta_right are the residual directions, where E
##       does not vanish, so they inherit the full pooled misspecification.
##   (5) PREDICTION: move the trough off the median (p1 != 0.5) and gamma's
##       protection must degrade in proportion to |1 - 2*tau_trough|.
WS <- "/user/work/fh6520/gmoq_covmixW"
.libPaths(c(file.path(WS, "lib")))
stopifnot(identical(.libPaths()[1], file.path(WS, "lib")))
library(fungwasStage2)
source(file.path(WS, "src", "truth_mix.R"))
TAUS <- seq(0.02, 0.98, length.out = 45); N <- 461320L

cat("\n########## M1. E on its own diagonal is Delta(tau)*(1-2tau) ##########\n")
for (sh in c("normal", "mix5.89")) {
  tr <- pop_truth(TAUS, 0.01, sh, 0, 0.2)
  Ed <- diag(tr$Bw - tr$Bbridge)
  pred <- tr$Delta * (1 - 2 * TAUS)
  cat(sprintf("  %-9s max|diag(E) - Delta*(1-2tau)| = %.3e   (max|diag E| = %.3e)\n",
              sh, max(abs(Ed - pred)), max(abs(Ed))))
  i <- which.min(abs(TAUS - 0.5))
  cat(sprintf("             at tau=%.3f : Delta=%.3e  diag(E)=%.3e\n",
              TAUS[i], tr$Delta[i], Ed[i]))
}

cat("\n########## M2. spectrum of Sigma_pooled vs Sigma_true ##########\n")
cat("If the leading eigen-DIRECTION is the tau=0.5 trough spike, its eigenVALUE\n")
cat("must agree between the two covariance sources while the rest diverge.\n\n")
cat(sprintf("%-10s %6s | %10s %10s %9s | %10s %10s %9s | %8s\n",
            "shape", "share1", "lam1 pool", "lam1 true", "ratio", "lam2 pool",
            "lam2 true", "ratio", "ang(v1)"))
for (sh in c("normal", "mix2", "mix4", "mix5.89", "mix6")) {
  tr <- pop_truth(TAUS, 0.01, sh, 0, 0.2)
  e0 <- eigen(tr$Sigma0, symmetric = TRUE); eT <- eigen(tr$Sigma_true, symmetric = TRUE)
  share <- e0$values[1] / sum(e0$values)
  ang <- NA
  if (grepl("^mix", sh)) {
    W3 <- build_Ws(TAUS, tr)$mix3
    u <- W3[, "gamma"]
    ang <- acos(min(1, abs(sum(u * e0$vectors[, 1])) / sqrt(sum(u^2)))) * 180 / pi
  }
  cat(sprintf("%-10s %6.4f | %10.3e %10.3e %9.5f | %10.3e %10.3e %9.5f | %8s\n",
              sh, share, e0$values[1], eT$values[1], e0$values[1] / eT$values[1],
              e0$values[2], eT$values[2], e0$values[2] / eT$values[2],
              ifelse(is.na(ang), "-", sprintf("%.3f", ang))))
}

cat("\n########## M3. where each column's u sits relative to the leading eigvec ##########\n")
tr <- pop_truth(TAUS, 0.01, "mix5.89", 0, 0.2)
Ws <- build_Ws(TAUS, tr); e0 <- eigen(tr$Sigma0, symmetric = TRUE)
cat(sprintf("%-6s %-12s %10s %10s %10s\n", "W", "param",
            "ang(u,v1)", "share on v1", "se_ratio"))
for (wn in names(Ws)) {
  W <- Ws[[wn]]; A <- solve(crossprod(W), t(W)); ex <- exact_W(W, tr, N * tr$Vg)
  for (j in seq_len(ncol(W))) {
    u <- A[j, ]
    cs <- abs(sum(u * e0$vectors[, 1])) / sqrt(sum(u^2))
    ## fraction of u'Sigma u carried by the leading eigendirection
    pr <- e0$values[1] * sum(u * e0$vectors[, 1])^2 / as.numeric(t(u) %*% tr$Sigma0 %*% u)
    cat(sprintf("%-6s %-12s %10.3f %10.5f %10.4f\n", wn, colnames(W)[j],
                acos(min(1, cs)) * 180 / pi, pr, ex$ratio[j]))
  }
}

cat("\n########## M4. BREAKING IT: move the trough off the median (p1 != 0.5) ##########\n")
cat("Ztarget 30 under a TRUE DISPERSION effect at MAF 0.01, S = 5.89 throughout.\n")
cat("tau_trough = F(argmin f).  Protection of gamma should track |1 - 2*tau_trough|.\n\n")
maxZ <- function(maf, sh, a, cc, wn) {
  t2 <- pop_truth(TAUS, maf, sh, a, cc)
  max(abs(exact_W(build_Ws(TAUS, t2)[[wn]], t2, N * t2$Vg)$Z_true))
}
solve_c <- function(maf, sh, Zt) {
  hi <- 1; while (maxZ(maf, sh, 0, hi, "mix3") < Zt) { hi <- hi * 2; if (hi > 512) return(NA) }
  uniroot(function(v) maxZ(maf, sh, 0, v, "mix3") - Zt, c(1e-9, hi), tol = 1e-10)$root
}
cat(sprintf("%-14s %5s %10s %10s | %-12s %9s %9s %12s\n",
            "shape", "p1", "tau_trough", "|1-2t_tr|", "param", "se_ratio", "cov95", "infl@5e-8"))
for (p1 in c(0.50, 0.40, 0.30, 0.20, 0.10)) {
  sh <- if (p1 == 0.5) "mix5.89" else sprintf("mix5.89p%g", p1)
  cc <- solve_c(0.01, sh, 30); if (!is.finite(cc)) { cat(sh, " Z=30 unreachable\n"); next }
  tr <- pop_truth(TAUS, 0.01, sh, 0, cc)
  ## trough of the density, on a fine grid, expressed as a tau
  gr <- seq(min(tr$q_tau), max(tr$q_tau), length.out = 4001)
  fg <- vapply(gr, function(y) sum(dbinom(0:2, 2, 0.01) *
        eps_pdf((y - 0) / (1 + cc * (0:2)), tr$shape) / (1 + cc * (0:2))), 0)
  inner <- gr > quantile(tr$q_tau, .15) & gr < quantile(tr$q_tau, .85)
  ytr <- gr[inner][which.min(fg[inner])]
  ttr <- sum(dbinom(0:2, 2, 0.01) * eps_cdf((ytr) / (1 + cc * (0:2)), tr$shape))
  ex <- exact_W(build_Ws(TAUS, tr)$mix3, tr, N * tr$Vg)
  for (j in seq_along(ex$ratio)) {
    cat(sprintf("%-14s %5.2f %10.4f %10.4f | %-12s %9.4f %9.4f %12.4g\n",
                sh, p1, ttr, abs(1 - 2 * ttr), names(ex$ratio)[j],
                ex$ratio[j], ex$cov95[j], ex$infl_5e8[j]))
  }
}
cat("\nDONE\n")
