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

## the residualised statistic is calibrated on truly-null dispersion effects
idx <- abs(d_true) < 1e-9   # none exactly, so build an explicitly null set instead
alpha2  <- rnorm(M, 0, 0.02)
v2_hat  <- r_true * alpha2 + se_v * rnorm(M)      # NO genuine dispersion effect at all
a2_hat  <- alpha2 + se_a * rnorm(M)
r2 <- estimate_rav(a2_hat, v2_hat, se_a, se_v, rep(0, M), trim = 0)$r_av
before <- mean(abs(v2_hat / se_v) > 1.96)
rr <- residualise_dispersion(a2_hat, v2_hat, se_a, se_v, rep(0, M), r_av = r2)
after <- mean(abs(rr$z_sigma2_resid) > 1.96)
cat(sprintf("   rejection at 0.05 before correction %.4f, after %.4f (nominal 0.05)\n",
            before, after))
stopifnot(before > 0.10, abs(after - 0.05) < 0.01)
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

cat("\nAll dispersion-guard tests passed.\n")
