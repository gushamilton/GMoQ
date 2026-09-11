# ============================================================
# Test: the Q calibration on a spectrum with KNOWN negative eigenvalues.
#
# The stored int8 covariance is not positive definite for essentially every
# SNP, so B = D^-1/2 M Sigma M' D^-1/2 legitimately has negative eigenvalues
# that carry real mass in Q = sum_i lambda_i chi2_1.  Two things have to hold:
#
#  1. the eigenvalue retention rule keeps them -- it must select by MAGNITUDE,
#     so that only the ncol(W) structural zeros are dropped;
#  2. the saddlepoint root-find stays inside the CGF's domain of convergence,
#     which a two-sided spectrum bounds BELOW as well as above, at
#     1/(2 lambda_min).  Leaving it does not crash: 1 - 2 lambda_i t merely
#     goes negative for some i and K'/K'' return values off the wrong branch.
#
# Every reference here is computed independently in plain R (exact power sums,
# a bracketed Lugannani-Rice, and Imhof numerical inversion).
# ============================================================

suppressPackageStartupMessages(library(Rcpp))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg) else "tests/test-q-negative-eigenvalues.R"
pkg_root <- normalizePath(file.path(dirname(this_file), ".."))

Rcpp::sourceCpp(file.path(pkg_root, "src", "stage2_se_helper.cpp"))
source(file.path(pkg_root, "R", "utils.R"))
source(file.path(pkg_root, "R", "covariance.R"))
source(file.path(pkg_root, "R", "weight-functions.R"))
source(file.path(pkg_root, "R", "gwas.R"))

ok <- function(cond, msg) {
  if (isTRUE(cond)) cat("PASS:", msg, "\n") else stop("FAIL: ", msg)
}

upper_row_major <- function(S) {
  out <- numeric(nrow(S) * (nrow(S) + 1) / 2); idx <- 1L
  for (r in seq_len(nrow(S))) for (c in r:ncol(S)) { out[idx] <- S[r, c]; idx <- idx + 1L }
  out
}

# Build a symmetric matrix with an exactly prescribed spectrum.
with_spectrum <- function(lam, seed) {
  set.seed(seed)
  n <- length(lam)
  U <- qr.Q(qr(matrix(rnorm(n * n), n, n)))
  S <- U %*% diag(lam) %*% t(U)
  (S + t(S)) / 2
}

# ---- independent references -------------------------------------------
kprime <- function(t, lam) sum(lam / (1 - 2 * lam * t))

# Lugannani-Rice with a bracket derived from the domain, solved by uniroot.
lr_ref <- function(q, lam) {
  hi_lim <- 0.5 / max(lam)
  stopifnot(q > sum(lam))                      # upper tail
  hi <- hi_lim * 0.5
  while (kprime(hi, lam) <= q) hi <- hi + 0.5 * (hi_lim - hi)
  z <- stats::uniroot(function(t) kprime(t, lam) - q, c(0, hi),
                      tol = .Machine$double.eps^0.9)$root
  den <- 1 - 2 * lam * z
  K  <- -0.5 * sum(log(den))
  K2 <- sum(2 * lam^2 / den^2)
  w <- sign(z) * sqrt(2 * (z * q - K)); v <- z * sqrt(K2)
  list(zeta = z, p = stats::pnorm(w + log(v / w) / w, lower.tail = FALSE))
}

# Exponentially tilted importance sampling: an EXACT reference that does not
# reuse any of the machinery under test beyond the tilt point.  Tilting a
# chi2_1 by theta rescales it by 1/(1 - 2 theta), so under the tilt at t
# X_i ~ chi2_1/(1 - 2 lam_i t) and P(Q>q) = E_tilt[exp(-t Q + K(t)) 1{Q>q}]
# for any t strictly inside the CGF domain.  Fixed seed, so deterministic.
# (Imhof/Davies numerical inversion was tried and is not usable here: the
# tail probability is a ~1e-7 residue after cancelling a 0.5-sized integral.)
is_ref <- function(q, lam, t, n = 2e6, seed = 1234L) {
  set.seed(seed)
  den <- 1 - 2 * lam * t
  stopifnot(all(den > 0))
  Kt <- -0.5 * sum(log(den))
  k <- length(lam); nb <- 20L; per <- n %/% nb
  s <- 0; s2 <- 0
  for (b in seq_len(nb)) {
    z <- matrix(stats::rnorm(per * k), per, k)
    Qs <- (z * z) %*% (lam / den)
    w <- ifelse(Qs > q, exp(-t * Qs + Kt), 0)
    s <- s + sum(w); s2 <- s2 + sum(w * w)
  }
  m <- s / (nb * per)
  list(p = m, rse = sqrt(max(s2 / (nb * per) - m * m, 0) / (nb * per)) / m)
}

# ============================================================
# 1. Negative eigenvalues are RETAINED, and the power sums are exact.
# ============================================================
Tt <- 20L
# Most negative eigenvalue -0.9: the CGF's lower domain edge is then
# 1/(2 * -0.9) = -0.556, so the pre-fix hard-coded starting bracket of -1.0
# lies OUTSIDE the domain of convergence.  This is not a contrived regime:
# the real int8 bmi spectra reach -0.84 against a largest of 8-14.
lam_true <- c(12, 7, 4, 2.5, 1.6, 1.1, 0.7, 0.45, 0.3, 0.2, 0.12, 0.07, 0.04,
              -0.9, -0.55, -0.3, -0.16, -0.08, -0.03, -0.01)
ok(sum(lam_true < 0) == 7L, "test spectrum really has negative eigenvalues")
ok(min(lam_true) < -0.5,
   "test spectrum puts the pre-fix starting bracket lo = -1 outside the CGF domain")

S <- with_spectrum(lam_true, seed = 11)
cov_vec <- upper_row_major(S)
offsets <- 0
Dvar <- matrix(1, Tt, 1L)            # D = I
M_id  <- diag(Tt)                    # so B = Sigma exactly
q_obs <- 6 * sum(lam_true)

cal <- compute_q_calibration(cov_vec, offsets, q_obs, Dvar, M_id, Tt, 1L)

ok(cal[1, 6] == Tt,
   "no eigenvalue is discarded when M is full rank (nothing is structurally zero)")
for (k in 1:4) {
  ok(abs(cal[1, k] - sum(lam_true^k)) / abs(sum(lam_true^k)) < 1e-10,
     sprintf("power sum c%d equals the exact value including the negative eigenvalues", k))
}
# The point of the fix, stated as a number: a one-sided rule would return the
# positive part of the trace instead, which is wrong by this much.
c1_pos_only <- sum(lam_true[lam_true > 0])
ok(abs(cal[1, 1] - c1_pos_only) / c1_pos_only > 0.05,
   "c1 differs materially from the positive-only sum a one-sided rule would give")

# ============================================================
# 2. The structural zeros -- and only those -- are still dropped.
# ============================================================
W <- cbind(1, seq_len(Tt) / Tt, (seq_len(Tt) / Tt)^2)
kW <- ncol(W)
Mp <- diag(Tt) - W %*% solve(crossprod(W), t(W))
# A genuinely indefinite Sigma with a wide dynamic range.
S2 <- with_spectrum(c(seq(9, 0.02, length.out = 15), -c(0.8, 0.4, 0.2, 0.05, 0.01)),
                    seed = 23)
cv2 <- upper_row_major(S2)
d2  <- pmax(diag(S2), 1e-6)
B2  <- diag(1 / sqrt(d2)) %*% Mp %*% S2 %*% t(diag(1 / sqrt(d2)) %*% Mp)
ev2 <- eigen((B2 + t(B2)) / 2, symmetric = TRUE, only.values = TRUE)$values
ev2_keep <- ev2[abs(ev2) > 1e-8 * max(abs(ev2))]
ok(length(ev2_keep) == Tt - kW,
   "the R reference retains exactly T - ncol(W) eigenvalues")
ok(sum(ev2_keep < 0) > 0, "and some of those retained are negative")

cal2 <- compute_q_calibration(cv2, 0, 3 * sum(ev2_keep), matrix(d2, ncol = 1),
                             Mp, Tt, 1L)
ok(cal2[1, 6] == Tt - kW, "C++ drops exactly the ncol(W) structural zeros")
for (k in 1:4) {
  ok(abs(cal2[1, k] - sum(ev2_keep^k)) / abs(sum(ev2_keep^k)) < 1e-8,
     sprintf("c%d matches the R reference on the projected indefinite case", k))
}

# ============================================================
# 3. The saddlepoint root-find stays in the domain of convergence.
# ============================================================
for (mult in c(1.5, 3, 6, 12, 30)) {
  qq <- mult * sum(lam_true)
  cc <- compute_q_calibration(cov_vec, 0, qq, Dvar, M_id, Tt, 1L)
  p_cpp <- cc[1, 5]
  ref <- lr_ref(qq, lam_true)
  ok(p_cpp >= 0, sprintf("saddlepoint converges at q/E[Q] = %g (two-sided spectrum)", mult))
  ok(abs(log(p_cpp) - log(ref$p)) < 1e-8,
     sprintf("saddlepoint matches the bracketed R reference at q/E[Q] = %g", mult))
  # The root itself must be strictly inside (1/(2 lam_min), 1/(2 lam_max)).
  ok(ref$zeta > 0.5 / min(lam_true) && ref$zeta < 0.5 / max(lam_true),
     sprintf("the saddlepoint root is inside the CGF domain at q/E[Q] = %g", mult))
}

# Against an entirely different method: exact importance-sampled Monte Carlo.
for (mult in c(3, 6, 20)) {
  qq <- mult * sum(lam_true)
  p_cpp <- compute_q_calibration(cov_vec, 0, qq, Dvar, M_id, Tt, 1L)[1, 5]
  ref <- is_ref(qq, lam_true, lr_ref(qq, lam_true)$zeta)
  ok(ref$rse < 0.02, sprintf("IS reference is precise at q/E[Q] = %g (rse %.3f%%)",
                             mult, 100 * ref$rse))
  ok(abs(p_cpp / ref$p - 1) < 0.12,
     sprintf("saddlepoint within 12%% of exact Monte Carlo at q/E[Q] = %g (%.3e vs %.3e, ratio %.3f)",
             mult, p_cpp, ref$p, p_cpp / ref$p))
}

# Plain Monte Carlo too, with no tilting at all, where it can reach.
set.seed(99)
qq <- 3 * sum(lam_true)
p_cpp <- compute_q_calibration(cov_vec, 0, qq, Dvar, M_id, Tt, 1L)[1, 5]
nmc <- 4e6L; hits <- 0L
for (b in 1:20) {
  z <- matrix(stats::rnorm((nmc %/% 20L) * Tt), nmc %/% 20L, Tt)
  hits <- hits + sum(((z * z) %*% lam_true) > qq)
}
p_mc <- hits / nmc
ok(abs(p_cpp / p_mc - 1) < 0.12,
   sprintf("saddlepoint within 12%% of plain Monte Carlo (%.3e vs %.3e, %d hits)",
           p_cpp, p_mc, hits))

# A far more hostile spectrum: the negative side dominates the domain, so the
# admissible interval for t is narrow and asymmetric.
lam_hard <- c(3, 1.2, 0.6, 0.2, 0.05, -2.5, -1.1, -0.4, -0.09)
ok(sum(lam_hard) > 0, "hard spectrum still has a positive mean")
Sh <- with_spectrum(lam_hard, seed = 31); Th <- length(lam_hard)
for (mult in c(2, 5, 15, 50)) {
  qq <- mult * sum(lam_hard)
  p_cpp <- compute_q_calibration(upper_row_major(Sh), 0, qq,
                                 matrix(1, Th, 1L), diag(Th), Th, 1L)[1, 5]
  ok(p_cpp >= 0 && p_cpp <= 1,
     sprintf("saddlepoint converges on the hostile spectrum at q/E[Q] = %g", mult))
  ok(abs(log(p_cpp) - log(lr_ref(qq, lam_hard)$p)) < 1e-8,
     sprintf("and matches the bracketed reference at q/E[Q] = %g", mult))
}

# ============================================================
# 4. No regression on a positive definite spectrum.
# ============================================================
lam_pd <- c(12, 7, 4, 2.5, 1.6, 1.1, 0.7, 0.45, 0.3, 0.2, 0.12, 0.07, 0.04,
            0.02, 0.01, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4)
Spd <- with_spectrum(lam_pd, seed = 17)
cpd <- compute_q_calibration(upper_row_major(Spd), 0, 5 * sum(lam_pd),
                             Dvar, M_id, Tt, 1L)
ok(cpd[1, 6] == Tt, "PD case retains the whole spectrum")
for (k in 1:4) {
  ok(abs(cpd[1, k] - sum(lam_pd^k)) / sum(lam_pd^k) < 1e-10,
     sprintf("PD power sum c%d is exact", k))
}
ok(abs(log(cpd[1, 5]) - log(lr_ref(5 * sum(lam_pd), lam_pd)$p)) < 1e-9,
   "PD saddlepoint matches the reference (lo = 0 bracket loses nothing)")
# The smallest PD eigenvalue here is 2e-4 = 1.7e-5 of the largest, well above
# the 1e-8 retention tolerance: real structure is never at risk of being cut.
ok(min(lam_pd) / max(lam_pd) > 1e-6,
   "smallest genuine eigenvalue sits far above the retention tolerance")

# ============================================================
# 5. Liu's four-moment fit refuses a non-positive skewness.
# ============================================================
# Enough negative weight to drive c3 = sum lambda^3 below zero: the moments
# are still exact, but there is no (non-central) chi-square with negative
# skewness to match them to, so the approximation must return NA rather than
# a number from the wrong shape family.
lam_skew <- c(1.0, 0.9, 0.8, 0.7, 0.6, -1.6, -0.5, -0.2)
ok(sum(lam_skew) > 0 && sum(lam_skew^3) < 0,
   "constructed spectrum has positive mean but negative third power sum")
Ts <- length(lam_skew) + 1L
# Put that spectrum on the range of the residual projector exactly: with
# Sigma = V diag(lam_skew) V' and V an orthonormal basis of the complement of
# W, B = D^-1/2 M Sigma M' D^-1/2 has spectrum lam_skew plus one structural
# zero.  This exercises the real projector path, not a shortcut.
set.seed(53)
Wk <- matrix(rnorm(Ts), Ts, 1L)
Mk <- diag(Ts) - Wk %*% solve(crossprod(Wk), t(Wk))
Vk <- qr.Q(qr(Mk))[, seq_len(Ts - 1L), drop = FALSE]
Ss <- Vk %*% diag(lam_skew) %*% t(Vk); Ss <- (Ss + t(Ss)) / 2
q_skew <- 8 * sum(lam_skew)
cs <- compute_q_calibration(upper_row_major(Ss), 0, q_skew,
                            matrix(1, Ts, 1L), Mk, Ts, 1L)
ok(cs[1, 6] == Ts - 1L, "one structural zero dropped, the rest kept")
ok(abs(cs[1, 1] - sum(lam_skew)) / abs(sum(lam_skew)) < 1e-8, "c1 exact here")
ok(cs[1, 3] < 0 && abs(cs[1, 3] - sum(lam_skew^3)) / abs(sum(lam_skew^3)) < 1e-8,
   "c3 is exact and negative")
cal_l <- .q_pvalue_calibrated(
  Q = q_skew,
  cov_data = list(cov_vec = upper_row_major(Ss), offsets = 0),
  SE_tau = matrix(1, Ts, 1L), W = Wk,
  method = "liu", n_threads = 1L)
ok(is.na(cal_l$p_liu[1]),
   "liu returns NA when the skewness a chi-square could match is not available")
ok(is.finite(cal_l$p_saddlepoint[1]) && cal_l$p_saddlepoint[1] >= 0,
   "the saddlepoint, which makes no shape assumption, still converges there")
ok(abs(log(cal_l$p_saddlepoint[1]) - log(lr_ref(q_skew, lam_skew)$p)) < 1e-8,
   "and agrees with the bracketed R reference on that spectrum")

cat("\nAll negative-eigenvalue tests passed.\n")
