## Controls that must pass before any number is drawn.
WS <- "/user/work/fh6520/gmoq_covmixW"
.libPaths(c(file.path(WS, "lib")))
stopifnot(identical(.libPaths()[1], file.path(WS, "lib")))
library(fungwasStage2)
source(file.path(WS, "src", "truth_cs.R"))   # PUBLISHED engine -> pop_truth etc
pop_truth_pub <- pop_truth
cs_theta_truth_pub <- cs_theta_truth
source(file.path(WS, "src", "truth_mix.R"))  # this study's engine (overrides)

TAUS <- seq(0.02, 0.98, length.out = 45)
ok <- function(lbl, cond, extra = "") cat(sprintf("%-62s %s %s\n", lbl,
  if (isTRUE(cond)) "PASS" else "**FAIL**", extra))

cat("== C1. new pop_truth reproduces the PUBLISHED one for shape=normal ==\n")
for (cell in list(c(0.01, 0, 0.3), c(0.005, 0, 1.6), c(0.30, 0.05, 0.1))) {
  A <- pop_truth_pub(TAUS, cell[1], "normal", cell[2], cell[3])
  B <- pop_truth(TAUS,     cell[1], "normal", cell[2], cell[3])
  ok(sprintf("  maf=%g a=%g c=%g  Sigma0 / Sigma_true identical", cell[1], cell[2], cell[3]),
     max(abs(A$Sigma0 - B$Sigma0)) < 1e-12 &&
     max(abs(A$Sigma_true - B$Sigma_true)) / max(abs(A$Sigma_true)) < 1e-12 &&
     max(abs(A$beta - B$beta)) < 1e-12,
     sprintf("relerr=%.3g", max(abs(A$Sigma_true - B$Sigma_true)) / max(abs(A$Sigma_true))))
}

cat("\n== C2. mixture eps law is self-consistent (cdf/pdf/quantile) ==\n")
for (S in c(2, 4, 5.89, 6)) {
  sh <- paste0("mix", S)
  m  <- integrate(function(e) e * eps_pdf(e, sh), -20, 20, rel.tol = 1e-12)$value
  v  <- integrate(function(e) e^2 * eps_pdf(e, sh), -20, 20, rel.tol = 1e-12)$value - m^2
  qq <- eps_quantile(c(0.1, 0.5, 0.9), sh)
  fd <- (eps_cdf(qq + 1e-6, sh) - eps_cdf(qq - 1e-6, sh)) / 2e-6
  ok(sprintf("  S=%-5g mean 0, var 1, pdf = d/de cdf", S),
     abs(m) < 1e-10 && abs(v - 1) < 1e-9 && max(abs(fd - eps_pdf(qq, sh))) < 1e-6,
     sprintf("mean=%.2e var-1=%.2e", m, v - 1))
}

cat("\n== C3. mix3 W == build_skew_mixture()'s W at alpha = 0 ==\n")
## build_skew_mixture() forms gamma / beta_left / beta_right by finite-difference
## gradients of the skew-normal mixture CDF, divided by -f.  At alpha = 0 the
## skew-normal IS the normal, so this must equal make_weights_normal_mixture().
if (!requireNamespace("sn", quietly = TRUE)) {
  cat("  sn not installed; using pnorm as the alpha=0 skew-normal (identical by definition)\n")
  psn_ <- function(x, xi, omega, alpha) pnorm(x, xi, omega)
  dsn_ <- function(x, xi, omega, alpha) dnorm(x, xi, omega)
} else {
  psn_ <- sn::psn; dsn_ <- sn::dsn
}
p <- list(p1 = 0.5, xi1 = -0.8, omega1 = 0.55, alpha1 = 0,
          xi2 = 0.8, omega2 = 0.55, alpha2 = 0)
cdf <- function(x, z) z$p1 * psn_(x, z$xi1, z$omega1, z$alpha1) +
                      (1 - z$p1) * psn_(x, z$xi2, z$omega2, z$alpha2)
pdf <- function(x, z) z$p1 * dsn_(x, z$xi1, z$omega1, z$alpha1) +
                      (1 - z$p1) * dsn_(x, z$xi2, z$omega2, z$alpha2)
qt <- vapply(TAUS, function(t) uniroot(function(x) cdf(x, p) - t, c(-8, 8), tol = 1e-13)$root, 0)
grad_p <- make_fd_grad("p1", cdf)
W_prod <- make_weights_generic(
  taus = TAUS, q_tau = qt, dist_cdf = cdf, dist_pdf = pdf, params = p,
  grad_funcs = list(gamma = function(x, z) grad_p(x, z) * z$p1 * (1 - z$p1),
                    beta_left  = make_fd_grad("xi1", cdf),
                    beta_right = make_fd_grad("xi2", cdf)))
W_pkg <- make_weights_normal_mixture(TAUS, qt, p1 = p$p1, mu1 = p$xi1, sd1 = p$omega1,
                                     mu2 = p$xi2, sd2 = p$omega2, include_membership = TRUE)
rel <- max(abs(W_prod - W_pkg)) / max(abs(W_pkg))
ok("  build_skew_mixture(alpha=0) == make_weights_normal_mixture(TRUE)", rel < 1e-6,
   sprintf("max rel diff = %.3g ; dim %s / %s", rel,
           paste(dim(W_prod), collapse = "x"), paste(dim(W_pkg), collapse = "x")))
ok("  column names match production", identical(colnames(W_prod), c("gamma", "beta_left", "beta_right")),
   paste(colnames(W_prod), collapse = ","))

cat("\n== C4. column counts of the Ws actually used ==\n")
tr <- pop_truth(TAUS, 0.01, "mix5.89", 0, 0)
Ws <- build_Ws(TAUS, tr)
for (nm in names(Ws)) cat(sprintf("  %-6s ncol=%d  cols=%s\n", nm, ncol(Ws[[nm]]),
                                  paste(colnames(Ws[[nm]]), collapse = ",")))
ok("  vqtl2 has 2 columns", ncol(Ws$vqtl2) == 2L)
ok("  mix2  has 2 columns", ncol(Ws$mix2) == 2L)
ok("  mix3  has 3 columns (PRODUCTION)", ncol(Ws$mix3) == 3L)

cat("\n== C5. Sigma_0 == Sigma_true under the NULL (the reason nulls miss this) ==\n")
for (sh in c("normal", "mix5.89")) {
  tr0 <- pop_truth(TAUS, 0.01, sh, 0, 0)
  ok(sprintf("  shape=%-8s Delta(tau) == 0 and Sigma0 == Sigma_true", sh),
     max(abs(tr0$Delta)) < 1e-10 &&
     max(abs(tr0$Sigma0 - tr0$Sigma_true)) / max(abs(tr0$Sigma0)) < 1e-10,
     sprintf("max|Delta|=%.2e", max(abs(tr0$Delta))))
}

cat("\n== C6. first-order law  B_w - B = Delta(min) - t_s D_t - t_t D_s ==\n")
for (cc in c(0.02, 0.1, 0.3)) {
  trc <- pop_truth(TAUS, 0.01, "mix5.89", 0, cc)
  Dl <- trc$Delta
  Emat <- outer(seq_along(Dl), seq_along(Dl), Vectorize(function(s, t) {
    mn <- if (trc$q_tau[s] <= trc$q_tau[t]) Dl[s] else Dl[t]
    mn - trc$taus[s] * Dl[t] - trc$taus[t] * Dl[s]
  }))
  num <- max(abs(trc$Bw - trc$Bbridge - Emat))
  cat(sprintf("  c=%-5g max|B_w - B - E| = %.3e   max|B_w - B| = %.3e   (ratio %.4f)\n",
              cc, num, max(abs(trc$Bw - trc$Bbridge)),
              num / max(abs(trc$Bw - trc$Bbridge))))
}
cat("  (ratio -> 0 as c -> 0 confirms E is the exact first-order term)\n")

cat("\n== C7. Sigma-invariance of col(W): angle of W_gamma to leading eigvec ==\n")
for (S in c(2, 4, 5.89, 6)) {
  trS <- pop_truth(TAUS, 0.01, paste0("mix", S), 0, 0)
  WS3 <- build_Ws(TAUS, trS)$mix3
  ev <- eigen(trS$Sigma_true, symmetric = TRUE)
  lead <- ev$vectors[, 1]
  spec <- ev$values[1] / sum(ev$values)
  ang <- function(u, v) acos(min(1, abs(sum(u * v)) / sqrt(sum(u^2) * sum(v^2)))) * 180 / pi
  ## how much of the leading eigenvector lies in col(W)
  proj <- function(W, x) { P <- W %*% solve(crossprod(W), t(W)); as.numeric(P %*% x) }
  r3 <- sqrt(sum(proj(WS3, lead)^2) / sum(lead^2))
  r2 <- sqrt(sum(proj(WS3[, -1, drop = FALSE], lead)^2) / sum(lead^2))
  cat(sprintf("  S=%-5g lead eigval share %.4f  angle(W_gamma, v1)=%7.3f deg  |P_W3 v1|=%.5f  |P_W2 v1|=%.5f\n",
              S, spec, ang(WS3[, "gamma"], lead), r3, r2))
}
cat("\nDONE\n")
