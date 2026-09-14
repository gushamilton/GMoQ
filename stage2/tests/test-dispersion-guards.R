#!/usr/bin/env Rscript
## test-dispersion-guards.R
##
## Covers stage2/R/dispersion-guards.R: the genotype-skewness guard, the
## genome-wide mean-variance coupling estimator, and the residualised dispersion
## statistic.  Self-contained, base R only, no network, a few seconds.

suppressMessages(library(fungwasStage2))
set.seed(20260914)
ok <- function(msg) cat(sprintf("  ok  %s\n", msg))
near <- function(a, b, tol) all(abs(a - b) <= tol)

cat("== dispersion_skew_guard ==\n")

## gskew from MAF, against the closed form
p <- c(0.005, 0.01, 0.05, 0.1, 0.3, 0.5)
gsk <- (1 - 2 * p) / sqrt(2 * p * (1 - p))
G <- dispersion_skew_guard(z_mu = rep(10, length(p)), N = 461320, maf = p)
stopifnot(near(G, 0.6855 * gsk * 100 / sqrt(461320), 1e-12))
ok("G matches 0.6855 * gskew * Z^2 / sqrt(N) from MAF")

## MAF 0.5 has zero genotype skewness, so no artefact however large the mean effect
stopifnot(abs(dispersion_skew_guard(z_mu = 200, N = 461320, maf = 0.5)) < 1e-12)
ok("G is exactly zero at MAF 0.5, for any Z_mu")

## recovery from d_resid: v = d_resid/N ~= 2p(1-p), gskew magnitude = sqrt((1-2v)/v)
N <- 461320
d_resid <- 2 * p * (1 - p) * N
Gn <- dispersion_skew_guard(z_mu = rep(10, length(p)), N = N, n_eff = d_resid)
stopifnot(near(abs(Gn), abs(G), 1e-9))
ok("G recovered from d_resid agrees with G from MAF (magnitude)")

## it is a magnitude: p and 1-p are indistinguishable from d_resid alone
stopifnot(near(dispersion_skew_guard(10, N = N, n_eff = 2 * 0.3 * 0.7 * N),
               dispersion_skew_guard(10, N = N, n_eff = 2 * 0.7 * 0.3 * N), 1e-12))
ok("d_resid route cannot distinguish p from 1-p, as documented")

stopifnot(inherits(try(dispersion_skew_guard(1), silent = TRUE), "try-error"))
ok("refuses when neither maf nor (n_eff, N) is supplied")


cat("\n== estimate_rav ==\n")

## Simulate summary statistics with a KNOWN coupling and known sampling error.
M <- 20000
r_true <- 0.135                                   # Young et al.'s UK Biobank BMI value
alpha  <- rnorm(M, 0, 0.02)                       # true mean effects
d_true <- rnorm(M, 0, 0.004)                      # genuine dispersion effects, independent of alpha
alpha_v <- r_true * alpha + d_true
se_a <- runif(M, 0.010, 0.020)                    # sampling error comparable to the signal
se_v <- runif(M, 0.004, 0.008)
rho  <- 0.3                                       # within-variant sampling correlation
cov_l <- rho * se_a * se_v
e1 <- rnorm(M); e2 <- rho * e1 + sqrt(1 - rho^2) * rnorm(M)
a_hat <- alpha   + se_a * e1
v_hat <- alpha_v + se_v * e2
blk <- rep(1:22, length.out = M)

fit <- estimate_rav(a_hat, v_hat, se_a, se_v, cov_l, block = blk, trim = 0)
cat(sprintf("   r_av = %.4f (se %.4f, %d blocks), truth %.4f\n",
            fit$r_av, fit$se, fit$n_blocks, r_true))
stopifnot(abs(fit$r_av - r_true) < 4 * fit$se)
ok("recovers a known r_av within 4 jackknife SE")
stopifnot(is.finite(fit$se), fit$se > 0, fit$n_blocks == 22L)
ok("block jackknife SE is finite and positive")

## The two corrections each matter.  Naive OLS of v_hat on a_hat is contaminated by
## BOTH the within-variant sampling covariance and regression dilution.
naive <- sum(a_hat * v_hat) / sum(a_hat^2)
no_cov <- sum(a_hat * v_hat) / sum(a_hat^2 - se_a^2)
cat(sprintf("   naive OLS slope          %.4f  (truth %.4f)\n", naive, r_true))
cat(sprintf("   dilution-corrected only  %.4f\n", no_cov))
cat(sprintf("   both corrections         %.4f\n", fit$r_av))
stopifnot(abs(fit$r_av - r_true) < abs(naive - r_true))
ok("the corrected estimator beats naive OLS")

## With zero within-variant covariance the cov correction is a no-op
f0 <- estimate_rav(a_hat, v_hat, se_a, se_v, rep(0, M), trim = 0)
stopifnot(near(f0$r_av, sum(a_hat * v_hat) / sum(a_hat^2 - se_a^2), 1e-12))
ok("with cov = 0 it reduces to the dilution-corrected slope")

stopifnot(inherits(try(estimate_rav(1:5, 1:5, rep(1, 5), rep(1, 5), rep(0, 5)),
                       silent = TRUE), "try-error"))
ok("refuses fewer than 20 usable variants")

## missing covariance warns rather than silently assuming zero
w <- tryCatch({ estimate_rav(a_hat, v_hat, se_a, se_v, NULL, trim = 0); "nowarn" },
              warning = function(e) "warned")
stopifnot(identical(w, "warned"))
ok("warns when the within-variant covariance is not supplied")


cat("\n== residualise_dispersion ==\n")

r_hat <- fit$r_av; se_r <- fit$se
res <- residualise_dispersion(a_hat, v_hat, se_a, se_v, cov_l, r_av = r_hat, se_rav = 0)
stopifnot(near(res$beta_sigma2_resid, v_hat - r_hat * a_hat, 1e-12))
ok("point estimate is alpha_v - r_av * alpha")
stopifnot(near(res$se_sigma2_resid,
               sqrt(se_v^2 + r_hat^2 * se_a^2 - 2 * r_hat * cov_l), 1e-12))
ok("SE is the delta method including the within-variant covariance")

## The residualised statistic is calibrated when the ONLY dispersion signal present
## is mean-variance coupling.  Use a strong coupling so the contrast is unambiguous:
## this mimics a raw lognormal, where the measured r_av is 0.7, not 0.135.
r_strong <- 0.70
alpha2  <- rnorm(M, 0, 0.02)
v2_hat  <- r_strong * alpha2 + se_v * rnorm(M)    # NO genuine dispersion effect at all
a2_hat  <- alpha2 + se_a * rnorm(M)
r2 <- estimate_rav(a2_hat, v2_hat, se_a, se_v, rep(0, M), trim = 0)$r_av
before <- mean(abs(v2_hat / se_v) > 1.96)
rr <- residualise_dispersion(a2_hat, v2_hat, se_a, se_v, rep(0, M), r_av = r2)
after <- mean(abs(rr$z_sigma2_resid) > 1.96)
cat(sprintf("   r_av recovered %.4f (truth %.2f)\n", r2, r_strong))
cat(sprintf("   rejection at 0.05 before correction %.4f, after %.4f (nominal 0.05)\n",
            before, after))
stopifnot(abs(r2 - r_strong) < 0.05)
stopifnot(before > 0.20, abs(after - 0.05) < 0.01)
ok("a pure mean-variance coupling is removed: rejection returns to nominal")

## including se_rav inflates the SE, and only for large |alpha|
r3 <- residualise_dispersion(a_hat, v_hat, se_a, se_v, cov_l, r_av = r_hat, se_rav = se_r)
stopifnot(all(r3$se_sigma2_resid >= res$se_sigma2_resid - 1e-15))
stopifnot(cor(r3$se_sigma2_resid - res$se_sigma2_resid, abs(a_hat)) > 0.5)
ok("se_rav inflates the SE, and does so in proportion to |alpha|")

stopifnot(inherits(try(residualise_dispersion(a_hat, v_hat, se_a, se_v, cov_l,
                                              r_av = c(1, 2)), silent = TRUE), "try-error"))
ok("refuses a non-scalar r_av")


cat("\n== dispersion_report ==\n")

P <- 50
fake <- list(
  params    = rbind(beta_mu = a_hat[1:P], beta_sigma2 = v_hat[1:P]),
  SE_params = rbind(beta_mu = se_a[1:P], beta_sigma2 = se_v[1:P]),
  param_cov_upper = cbind(cov_beta_mu_beta_mu = se_a[1:P]^2,
                          cov_beta_mu_beta_sigma2 = cov_l[1:P],
                          cov_beta_sigma2_beta_sigma2 = se_v[1:P]^2))
colnames(fake$params) <- colnames(fake$SE_params) <- paste0("rs", 1:P)
rep1 <- dispersion_report(fake, r_av = r_hat, se_rav = se_r, N = 461320, maf = 0.3)
stopifnot(nrow(rep1) == P)
stopifnot(all(c("beta_sigma2", "z_sigma2", "beta_sigma2_resid", "z_sigma2_resid",
                "G_skew_guard", "flag_skew") %in% names(rep1)))
ok("returns uncorrected and corrected columns side by side, plus the guard")
stopifnot(near(rep1$cov_mu_sigma2, cov_l[1:P], 1e-12))
ok("picks the mu/sigma2 covariance out of param_cov_upper by name")
stopifnot(near(rep1$G_skew_guard,
               dispersion_skew_guard(rep1$z_mu, N = 461320, maf = 0.3), 1e-12))
ok("guard column matches dispersion_skew_guard")

## the uncorrected column must be untouched by the correction
stopifnot(near(rep1$beta_sigma2, v_hat[1:P], 1e-12))
ok("the uncorrected estimate is preserved verbatim, not overwritten")

## omitting r_av omits the corrected columns rather than defaulting to some value
rep2 <- dispersion_report(fake, r_av = NULL, N = 461320, maf = 0.3)
stopifnot(!("beta_sigma2_resid" %in% names(rep2)))
ok("no r_av supplied => no corrected columns invented")

w2 <- tryCatch({ dispersion_report(list(params = fake$params, SE_params = fake$SE_params),
                                   r_av = r_hat); "nowarn" },
               warning = function(e) "warned")
stopifnot(identical(w2, "warned"))
ok("warns when param_cov_upper is absent")



cat("\n== dispersion_skew_bias / correct_dispersion_skew ==\n")

## --- the identity the whole correction rests on: A W = I, so the bias vector
## --- (which IS delta^2 k times the dispersion column) projects onto
## --- beta_sigma2 in full and onto beta_mu not at all -- for ols AND gls.
TAUS <- seq(0.02, 0.98, length.out = 45)
zt <- qnorm(TAUS); sdY <- 1.3
Wv <- make_weights_vqtl(TAUS, q_tau = sdY * zt, mu = 0, sd = sdY)
Sig <- outer(TAUS, TAUS, pmin) - outer(TAUS, TAUS) + diag(1e-8, length(TAUS))
delta <- 0.2; pfreq <- 0.01; kk <- 1 - 2 * pfreq
bias_b <- (delta^2 / (2 * sdY)) * kk * zt          # the derived bias in b(tau)
for (nm in c("ols", "gls")) {
  A <- if (nm == "ols") solve(crossprod(Wv), t(Wv)) else
       solve(t(Wv) %*% solve(Sig) %*% Wv, t(Wv) %*% solve(Sig))
  stopifnot(near(A %*% Wv, diag(2), 1e-10))
  proj <- drop(A %*% bias_b)
  stopifnot(abs(proj[1]) < 1e-12)                 # nothing lands on the location term
  stopifnot(near(proj[2], delta^2 * kk, 1e-12))   # all of it lands on the scale term
}
ok("A W = I, so the bias projects onto beta_sigma2 in full: identical for ols and gls")

## correcting b(tau) BEFORE the projection == correcting the coefficient AFTER
b_obs <- rep(delta, length(TAUS)) + bias_b
for (nm in c("ols", "gls")) {
  A <- if (nm == "ols") solve(crossprod(Wv), t(Wv)) else
       solve(t(Wv) %*% solve(Sig) %*% Wv, t(Wv) %*% solve(Sig))
  before <- drop(A %*% (b_obs - bias_b))
  after  <- drop(A %*% b_obs) - c(0, delta^2 * kk)
  stopifnot(near(before, after, 1e-12))
}
ok("correcting b(tau) first and correcting the coefficient after are identical")

## --- the point estimate ---------------------------------------------------
stopifnot(near(dispersion_skew_bias(0.2, 0.0, eaf = 0.01), 0.04 * 0.98, 1e-12))
ok("bias = (1 - 2p) * beta_mu^2")
stopifnot(near(dispersion_skew_bias(0.2, 0.01, eaf = 0.01),
               0.98 * (0.04 - 1e-4), 1e-12))
ok("debias subtracts Var(beta_mu_hat) from beta_mu_hat^2")
stopifnot(near(dispersion_skew_bias(0.2, 0.01, eaf = 0.01, debias = FALSE),
               0.98 * 0.04, 1e-12))
ok("debias = FALSE leaves the raw square")

stopifnot(abs(dispersion_skew_bias(5, 0.01, eaf = 0.5)) < 1e-15)
ok("exactly zero at EAF 0.5, for any beta_mu")

## a dispersion effect is invariant under an allele flip, and so is the correction:
## flipping sends p -> 1-p and beta -> -beta, and both b(tau) and the bias change sign.
stopifnot(near(dispersion_skew_bias(0.2, 0.01, eaf = 0.01),
               -dispersion_skew_bias(-0.2, 0.01, eaf = 0.99), 1e-14))
ok("covariant under an allele flip")

## the sign is not recoverable from an unoriented frequency, so it is refused
stopifnot(inherits(try(dispersion_skew_bias(0.2, 0.01), silent = TRUE), "try-error"))
ok("refuses to guess: an effect-allele frequency or mu3_over_var is mandatory")

## dosages / non-HWE: mu3_over_var overrides the 1 - 2p identity
stopifnot(near(dispersion_skew_bias(0.2, 0.0, mu3_over_var = 0.37), 0.04 * 0.37, 1e-12))
ok("mu3_over_var overrides the Hardy-Weinberg identity")

## --- the standard error ---------------------------------------------------
bm <- 0.2; sm <- 0.01; ss <- 0.012; cl <- 1.3e-5; kk <- 1 - 2 * 0.01
cc <- correct_dispersion_skew(bm, 0.05, sm, ss, cl, eaf = 0.01)
v_expect <- ss^2 + 4 * kk^2 * bm^2 * sm^2 + 2 * kk^2 * sm^4 - 4 * kk * bm * cl
stopifnot(near(cc$se_sigma2_skewcorr, sqrt(v_expect), 1e-14))
ok("SE is the stated delta-method variance including the within-variant covariance")
stopifnot(cc$se_sigma2_skewcorr > ss)
ok("the SE widens: the correction is estimated, and that is paid for")

## non-negative for any admissible covariance (|c| <= se_mu * se_sigma2)
cg <- seq(-1, 1, by = 0.05) * sm * ss
vv <- correct_dispersion_skew(rep(bm, length(cg)), rep(0.05, length(cg)),
                              sm, ss, cg, eaf = 0.01)
stopifnot(all(is.finite(vv$se_sigma2_skewcorr)), all(vv$se_sigma2_skewcorr > 0))
ok("variance stays positive across the whole admissible covariance range")

## length-1 arguments recycle against beta_mu
r1 <- correct_dispersion_skew(c(0.1, 0.2), 0.05, sm, ss, cl, eaf = 0.01)
r2 <- correct_dispersion_skew(c(0.1, 0.2), c(0.05, 0.05), c(sm, sm),
                              c(ss, ss), c(cl, cl), eaf = c(0.01, 0.01))
stopifnot(near(unlist(r1), unlist(r2), 1e-14), nrow(r1) == 2L)
ok("length-1 arguments recycle to length(beta_mu)")

## on a null variant the correction is a no-op to within its own noise
nn <- correct_dispersion_skew(0, 0.001, sm, ss, 0, eaf = 0.01)
## the whole correction is then -k*se_mu^2, which is a hundredth of a standard error
stopifnot(abs(nn$bias_skew) / ss < 0.01, nn$se_sigma2_skewcorr / ss < 1.0001)
ok("costs essentially nothing when there is no location effect")

## --- end to end: a pure location effect at a skewed genotype --------------
## RIF-OLS b(tau) generated at the POPULATION estimand, plus sampling noise, so
## every rejection of beta_sigma2 is false by construction.
set.seed(11)
Nrep <- 20000
kk <- 1 - 2 * 0.01
A <- solve(crossprod(Wv), t(Wv))
V <- A %*% Sig %*% t(A) / 4000          # a plausible parameter covariance
th <- c(0.2, 0.2^2 * kk)                 # truth: location 0.2, dispersion ZERO + bias
L <- chol(V)
Z <- matrix(rnorm(2 * Nrep), nrow = 2)
hat <- th + t(L) %*% Z
se1 <- sqrt(V[1, 1]); se2 <- sqrt(V[2, 2]); c12 <- V[1, 2]
r_unc <- mean(abs(hat[2, ] / se2) > 1.959964)
cc <- correct_dispersion_skew(hat[1, ], hat[2, ], se1, se2, c12, eaf = 0.01)
r_cor <- mean(abs(cc$z_sigma2_skewcorr) > 1.959964)
cat(sprintf("   false dispersion rate at 0.05: %.4f uncorrected, %.4f corrected (nominal 0.05)\n",
            r_unc, r_cor))
stopifnot(r_unc > 0.20)                  # the artefact is there
stopifnot(abs(r_cor - 0.05) < 0.01)      # and the correction removes it
ok("a pure location effect's false dispersion rate returns to nominal")

## --- dispersion_report wiring ---------------------------------------------
rep3 <- dispersion_report(fake, N = 461320, maf = 0.3, eaf = 0.3)
stopifnot(all(c("bias_skew", "beta_sigma2_skewcorr", "se_sigma2_skewcorr",
                "z_sigma2_skewcorr", "p_sigma2_skewcorr") %in% names(rep3)))
ok("dispersion_report adds the corrected columns when eaf is supplied")
stopifnot(near(rep3$beta_sigma2, v_hat[1:P], 1e-12))
ok("and still preserves the uncorrected estimate verbatim")
rep4 <- dispersion_report(fake, N = 461320, maf = 0.3)
stopifnot(!("beta_sigma2_skewcorr" %in% names(rep4)))
ok("no eaf supplied => no corrected columns invented")

cat("\nAll dispersion-guard tests passed.\n")
