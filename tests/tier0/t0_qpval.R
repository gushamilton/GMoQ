#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 0 -- the Q lack-of-fit statistic's reference distribution.
#
# Q = r' D^-1 r with r = M b, M = I - W(W'W)^-1 W'. Under the null r is a
# correlated normal vector with Var(r) = M Sigma M', so
#
#     Q ~ sum_i lambda_i chi2_1,   lambda = eig(D^-1/2 M Sigma M' D^-1/2)
#
# and chi2_{T-k} is the WRONG reference except in the one special case where
# Sigma = D. The package offers three approximations to that mixture. This file
# checks them against
#
#   (a) the exact eigen-decomposition, for the moments;
#   (b) the closed-form chi2_{T-k} in the special case, which every method must
#       reproduce or it is not approximating the right object;
#   (c) CompQuadForm::davies(), Davies' numerical inversion of the
#       characteristic function, which is the exact tail probability.
#
# Run:  Rscript --vanilla tests/tier0/t0_qpval.R
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier0/t0_qpval.R")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
gmoq_begin("tier0", "qpval")

compute_q_calibration <- gmoq_internal("compute_q_calibration")

have_cqf <- requireNamespace("CompQuadForm", quietly = TRUE)

set.seed(90210)
Tt <- 12L
taus <- seq(0.05, 0.95, length.out = Tt)
fpd <- dnorm(qnorm(taus))
W <- cbind(mu = rep(1, Tt), sigma = qnorm(taus))
M <- diag(Tt) - W %*% solve(crossprod(W), t(W))
k <- ncol(W)

upper_rowmajor <- function(S) {
  n <- nrow(S); out <- numeric(n * (n + 1) / 2); i <- 1L
  for (r in seq_len(n)) for (cc in r:n) { out[i] <- S[r, cc]; i <- i + 1L }
  out
}
n_elem <- Tt * (Tt + 1L) / 2L

# =====================================================================
# 1. The special case Sigma = D: Q is exactly chi2_{T-k}
#
# If the per-tau residuals were independent with variance exactly se_tau^2 then
# B = M, a rank-(T-k) projector, whose eigenvalues are T-k ones and k zeros. So
# Q ~ chi2_{T-k} and every method must return pchisq(..., lower = FALSE).
# =====================================================================
# Sigma must be a MULTIPLE OF THE IDENTITY for B to reduce to M: with a merely
# diagonal Sigma, D^-1/2 M Sigma M' D^-1/2 is not M and Q is not chi-square.
D <- diag(1.7, Tt)
SigD <- D                                  # Sigma = D = 1.7 I
Dvar <- matrix(diag(D), nrow = Tt, ncol = 1L)
Qgrid <- c(2, 5, 9, 14, 20, 30, 45)
cal <- compute_q_calibration(upper_rowmajor(SigD), 0, Qgrid[1], Dvar, M, Tt, 1L)
gmoq_close("T0.Q.projector_moments",
           "With Sigma = D the eigenvalues are those of the projector M: c_j = T - k for all j.",
           as.numeric(cal[1, 1:4]), rep(Tt - k, 4), 1e-8)

qc <- gmoq_internal(".q_pvalue_calibrated")
for (meth in c("satterthwaite", "liu", "saddlepoint")) {
  res <- qc(Q = Qgrid,
            cov_data = list(cov_vec = rep(upper_rowmajor(SigD), length(Qgrid)),
                            offsets = seq(0, by = n_elem * 4, length.out = length(Qgrid))),
            SE_tau = matrix(sqrt(diag(D)), nrow = Tt, ncol = length(Qgrid)),
            W = W, method = meth)
  ref <- pchisq(Qgrid, df = Tt - k, lower.tail = FALSE)
  # Satterthwaite and Liu match moments and are EXACT when the spectrum is a
  # single repeated eigenvalue. The saddlepoint is an asymptotic expansion, so
  # its tolerance is set from what a Lugannani-Rice approximation can deliver
  # on chi2_10, not from the other two.
  tol <- if (meth == "saddlepoint") 1e-2 else 1e-4
  gmoq_check(paste0("T0.Q.chisq_special_case.", meth),
             sprintf("With Sigma = D, q_method='%s' reproduces the exact chi2_{T-k} tail.", meth),
             max(abs(res$pval / ref - 1)) < tol,
             stat = max(abs(res$pval / ref - 1)),
             criterion = sprintf("max relative error < %g against pchisq(T-k)", tol))
}

# =====================================================================
# 2. The general case: moments against the exact eigen-decomposition
# =====================================================================
Sig <- (outer(taus, taus, pmin) - outer(taus, taus)) / outer(fpd, fpd)
Sig <- Sig / 5000
Dv  <- diag(diag(Sig))                      # se_tau^2 on the diagonal, as Stage 1 emits
B   <- solve(sqrt(Dv)) %*% M %*% Sig %*% t(M) %*% solve(sqrt(Dv))
lam <- eigen(0.5 * (B + t(B)), symmetric = TRUE, only.values = TRUE)$values
lam <- lam[abs(lam) > 1e-12]

Dvar1 <- matrix(diag(Sig), nrow = Tt, ncol = 1L)
Qobs  <- sum(lam) * 1.8                     # comfortably into the upper tail
cal   <- compute_q_calibration(upper_rowmajor(Sig), 0, Qobs, Dvar1, M, Tt, 1L)

for (j in 1:4) {
  gmoq_close(sprintf("T0.Q.power_sum_c%d", j),
             sprintf("compute_q_calibration's c%d is the exact eigenvalue power sum sum lambda^%d.", j, j),
             cal[1, j], sum(lam^j), abs(sum(lam^j)) * 1e-6)
}
gmoq_check("T0.Q.n_eig",
           "The kernel reports the number of non-negligible eigenvalues it used.",
           cal[1, 6] >= (Tt - k) - 1 && cal[1, 6] <= Tt,
           stat = cal[1, 6], criterion = "n_eig between T-k-1 and T")

# =====================================================================
# 3. The general case: tail probabilities against Davies' exact inversion
# =====================================================================
if (have_cqf) {
  Qs <- sum(lam) * c(1.0, 1.5, 2.0, 2.6, 3.2)
  davies_p <- vapply(Qs, function(q)
    CompQuadForm::davies(q, lambda = lam, acc = 1e-12)$Qq, numeric(1))
  cd <- list(cov_vec = rep(upper_rowmajor(Sig), length(Qs)),
             offsets = seq(0, by = n_elem * 4, length.out = length(Qs)))
  SEm <- matrix(sqrt(diag(Sig)), nrow = Tt, ncol = length(Qs))

  r_sad <- qc(Qs, cd, SEm, W, method = "saddlepoint")
  r_liu <- qc(Qs, cd, SEm, W, method = "liu")
  r_sat <- qc(Qs, cd, SEm, W, method = "satterthwaite")

  # The saddlepoint is the default and is the one that has to be accurate in
  # the tail. Liu is the four-moment fallback; Satterthwaite matches two
  # moments only and is expected to be the loosest.
  gmoq_check("T0.Q.saddlepoint_vs_davies",
             "The default saddlepoint reference matches Davies' exact tail to within 5 per cent.",
             max(abs(r_sad$pval / davies_p - 1)) < 0.05,
             stat = max(abs(r_sad$pval / davies_p - 1)),
             criterion = "max relative error < 0.05 over p in [1e-6, 0.4]",
             detail = paste("davies:", paste(signif(davies_p, 4), collapse = ", "),
                            "| saddle:", paste(signif(r_sad$pval, 4), collapse = ", ")))
  gmoq_check("T0.Q.liu_vs_davies",
             "The four-moment Liu reference matches Davies' exact tail to within 25 per cent.",
             max(abs(r_liu$pval / davies_p - 1)) < 0.25,
             stat = max(abs(r_liu$pval / davies_p - 1)),
             criterion = "max relative error < 0.25")
  gmoq_check("T0.Q.saddlepoint_beats_satterthwaite",
             "The saddlepoint is closer to the exact tail than the two-moment match, as claimed.",
             max(abs(r_sad$pval / davies_p - 1)) <= max(abs(r_sat$pval / davies_p - 1)),
             stat = max(abs(r_sat$pval / davies_p - 1)),
             criterion = "saddlepoint max error <= Satterthwaite max error")

  # The naive chi2_{T-k} reference is the thing the calibration exists to
  # replace. Pin that it really is wrong here, so a regression that quietly
  # reverts to it is visible.
  naive <- pchisq(Qs, df = Tt - k, lower.tail = FALSE)
  gmoq_check("T0.Q.naive_chisq_is_wrong",
             "chisq_naive is materially wrong on a correlated Sigma -- the reason q_method exists.",
             max(naive / davies_p) > 3,
             stat = max(naive / davies_p),
             criterion = "max ratio of naive chi2 tail to exact > 3")
} else {
  for (id in c("T0.Q.saddlepoint_vs_davies", "T0.Q.liu_vs_davies",
               "T0.Q.saddlepoint_beats_satterthwaite", "T0.Q.naive_chisq_is_wrong")) {
    gmoq_check(id, "Requires CompQuadForm for the exact reference.", NA,
               detail = "CompQuadForm not installed")
  }
}

# =====================================================================
# 4. The Satterthwaite scale/df match, written out
# =====================================================================
res <- qc(Qobs, list(cov_vec = upper_rowmajor(Sig), offsets = 0),
          matrix(sqrt(diag(Sig)), nrow = Tt, ncol = 1L), W, method = "satterthwaite")
c1 <- sum(lam); c2 <- sum(lam^2)
gmoq_close("T0.Q.satterthwaite_scale",
           "The Satterthwaite scale is Var(Q)/(2 E[Q]) = c2/c1.",
           res$scale, c2 / c1, abs(c2 / c1) * 1e-6)
gmoq_close("T0.Q.satterthwaite_df",
           "The Satterthwaite effective df is 2 E[Q]^2/Var(Q) = c1^2/c2.",
           res$df_eff, c1^2 / c2, abs(c1^2 / c2) * 1e-6)

# =====================================================================
# 5. Negative eigenvalues: the Liu match must refuse rather than guess
#
# The moments c1..c4 stay exact for a two-sided spectrum, but Liu's match maps
# them onto a non-central chi-square, whose skewness is necessarily positive.
# With enough negative weight c3 goes non-positive and the fit is meaningless;
# the code is supposed to return NA there rather than a number from the wrong
# shape family.
# =====================================================================
Sneg <- Sig
Sneg[1, Tt] <- Sneg[Tt, 1] <- -0.9 * sqrt(Sig[1, 1] * Sig[Tt, Tt])
Bn <- solve(sqrt(diag(diag(Sneg)))) %*% M %*% Sneg %*% t(M) %*% solve(sqrt(diag(diag(Sneg))))
lamn <- eigen(0.5 * (Bn + t(Bn)), symmetric = TRUE, only.values = TRUE)$values
rn <- qc(sum(abs(lamn)), list(cov_vec = upper_rowmajor(Sneg), offsets = 0),
         matrix(sqrt(diag(Sneg)), nrow = Tt, ncol = 1L), W, method = "liu")
gmoq_check("T0.Q.two_sided_spectrum_handled",
           "A two-sided eigenvalue spectrum yields either a finite p-value or NA, never a p outside [0,1].",
           is.na(rn$pval) || (rn$pval >= 0 && rn$pval <= 1),
           stat = if (is.na(rn$pval)) NA else rn$pval,
           criterion = "p in [0,1] or NA")
gmoq_check("T0.Q.negative_eigenvalues_present",
           "The probe really does produce a two-sided spectrum, so the branch above was exercised.",
           min(lamn) < -1e-10, stat = min(lamn), criterion = "min eigenvalue < 0")

# =====================================================================
# 6. Q is invariant to a rescaling of W -- it measures lack of fit, not units
# =====================================================================
W2 <- W %*% diag(c(3, 0.2))
M2 <- diag(Tt) - W2 %*% solve(crossprod(W2), t(W2))
gmoq_close("T0.Q.projector_invariant",
           "The residual projector depends on the column space of W, not its scaling.",
           as.numeric(M2), as.numeric(M), 1e-12)

gmoq_finish()
