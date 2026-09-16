## Does the membership column's protection survive a W built from a MIS-FITTED
## mixture?  Production estimates (p1, mu1, sd1, mu2, sd2) by EM on the observed
## phenotype; the exact engine above hands the W the population values.  Perturb
## each parameter and re-measure.  Also: how far is the W's own f from Stage 1's.
WS <- "/user/work/fh6520/gmoq_covmixW"
.libPaths(c(file.path(WS, "lib")))
stopifnot(identical(.libPaths()[1], file.path(WS, "lib")))
library(fungwasStage2)
source(file.path(WS, "src", "truth_mix.R"))
TAUS <- seq(0.02, 0.98, length.out = 45); N <- 461320L

## reference cell: testosterone-like separation, MAF 0.01, true dispersion Z30
SH <- "mix5.89"
maxZ <- function(cc) {
  t2 <- pop_truth(TAUS, 0.01, SH, 0, cc)
  max(abs(exact_W(build_Ws(TAUS, t2)$mix3, t2, N * t2$Vg)$Z_true))
}
hi <- 1; while (maxZ(hi) < 30) hi <- hi * 2
CC <- uniroot(function(v) maxZ(v) - 30, c(1e-9, hi), tol = 1e-10)$root
tr <- pop_truth(TAUS, 0.01, SH, 0, CC)
P  <- eps_mix_params(SH)
base <- list(p1 = P$p1,
             mu1 = tr$EY + tr$sdY * P$mu1, sd1 = tr$sdY * P$sd1,
             mu2 = tr$EY + tr$sdY * P$mu2, sd2 = tr$sdY * P$sd2)

meas <- function(pp) {
  W <- make_weights_normal_mixture(TAUS, tr$q_tau, p1 = pp$p1, mu1 = pp$mu1,
        sd1 = pp$sd1, mu2 = pp$mu2, sd2 = pp$sd2, include_membership = TRUE)
  exact_W(W, tr, N * tr$Vg)
}
b <- meas(base)
cat("\n########## R1. W built from MIS-FITTED mixture parameters ##########\n")
cat(sprintf("cell: shape=%s MAF=0.01 c=%.5f (Ztarget 30), N=%d\n", SH, CC, N))
cat("Each row perturbs ONE fitted parameter of the W by the stated relative\n")
cat("amount, holding the DGP fixed.  se_ratio = pooled reported SE / true SE.\n\n")
cat(sprintf("%-22s | %9s %9s %9s | %9s %9s\n", "W perturbation",
            "gamma", "beta_1", "beta_2", "gam cov95", "gam infl"))
cat(sprintf("%-22s | %9.4f %9.4f %9.4f | %9.4f %9.4f\n", "none (population W)",
            b$ratio[1], b$ratio[2], b$ratio[3], b$cov95[1], b$infl_5e8[1]))
pert <- list()
for (d in c(-0.20, -0.05, 0.05, 0.20)) {
  pert[[sprintf("p1 x (1%+.0f%%)", d * 100)]] <- within(base, p1 <- min(0.99, max(0.01, p1 * (1 + d))))
  pert[[sprintf("sd1,sd2 x (1%+.0f%%)", d * 100)]] <- within(base, { sd1 <- sd1 * (1 + d); sd2 <- sd2 * (1 + d) })
  pert[[sprintf("separation x (1%+.0f%%)", d * 100)]] <- within(base, {
    m <- (mu1 + mu2) / 2; mu1 <- m + (mu1 - m) * (1 + d); mu2 <- m + (mu2 - m) * (1 + d) })
  pert[[sprintf("mu1 shift %+.2f sd", d)]] <- within(base, mu1 <- mu1 + d * tr$sdY)
}
for (nm in names(pert)) {
  e <- meas(pert[[nm]])
  cat(sprintf("%-22s | %9.4f %9.4f %9.4f | %9.4f %9.4f\n", nm,
              e$ratio[1], e$ratio[2], e$ratio[3], e$cov95[1], e$infl_5e8[1]))
}

cat("\n########## R2. same, at the REAL-DATA operating point (|Z| ~ 11, MAF 0.30) ##########\n")
maxZ2 <- function(cc) {
  t2 <- pop_truth(TAUS, 0.30, SH, 0, cc)
  max(abs(exact_W(build_Ws(TAUS, t2)$mix3, t2, N * t2$Vg)$Z_true))
}
hi <- 1e-3; while (maxZ2(hi) < 11) hi <- hi * 2
CC2 <- uniroot(function(v) maxZ2(v) - 11, c(1e-9, hi), tol = 1e-10)$root
tr <- pop_truth(TAUS, 0.30, SH, 0, CC2)
base2 <- list(p1 = P$p1, mu1 = tr$EY + tr$sdY * P$mu1, sd1 = tr$sdY * P$sd1,
              mu2 = tr$EY + tr$sdY * P$mu2, sd2 = tr$sdY * P$sd2)
meas2 <- function(pp) {
  W <- make_weights_normal_mixture(TAUS, tr$q_tau, p1 = pp$p1, mu1 = pp$mu1,
        sd1 = pp$sd1, mu2 = pp$mu2, sd2 = pp$sd2, include_membership = TRUE)
  exact_W(W, tr, N * tr$Vg)
}
cat(sprintf("%-22s | %9s %9s %9s\n", "W perturbation", "gamma", "beta_1", "beta_2"))
e <- meas2(base2)
cat(sprintf("%-22s | %9.5f %9.5f %9.5f\n", "none (population W)", e$ratio[1], e$ratio[2], e$ratio[3]))
for (d in c(-0.20, 0.20)) {
  for (nm in c("p1", "sd", "sep")) {
    pp <- base2
    if (nm == "p1")  pp$p1 <- min(0.99, max(0.01, pp$p1 * (1 + d)))
    if (nm == "sd")  { pp$sd1 <- pp$sd1 * (1 + d); pp$sd2 <- pp$sd2 * (1 + d) }
    if (nm == "sep") { m <- (pp$mu1 + pp$mu2) / 2
                       pp$mu1 <- m + (pp$mu1 - m) * (1 + d); pp$mu2 <- m + (pp$mu2 - m) * (1 + d) }
    e <- meas2(pp)
    cat(sprintf("%-22s | %9.5f %9.5f %9.5f\n", sprintf("%s x (1%+.0f%%)", nm, d * 100),
                e$ratio[1], e$ratio[2], e$ratio[3]))
  }
}
cat("\nDONE\n")
