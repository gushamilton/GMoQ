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

## Trimming is defined on each variant's ratio contribution.  Give a small
## group an ordinary ratio but a much larger corrected denominator: numerator-
## only trimming discards that group, while ratio trimming retains it.
Mtrim <- 1000L
ratio_target <- c(seq(0.5, 1.5, length.out = 900L),
                  seq(1.0, 1.1, length.out = 100L))
a_trim <- c(rep(1, 900L), rep(10, 100L))
sa_trim <- rep(0.1, Mtrim)
den_trim <- a_trim^2 - sa_trim^2
v_trim <- ratio_target * den_trim / a_trim
sv_trim <- rep(0.1, Mtrim)
trimmed <- estimate_rav(a_trim, v_trim, sa_trim, sv_trim,
                        cov_mu_sigma2 = rep(0, Mtrim), trim = 0.1)
num_trim <- a_trim * v_trim
ratio_trim <- num_trim / den_trim
lim_trim <- quantile(ratio_trim, c(0.1, 0.9))
keep_trim <- is.finite(ratio_trim) & ratio_trim >= lim_trim[1] &
  ratio_trim <= lim_trim[2]
expected_trim <- sum(num_trim[keep_trim]) / sum(den_trim[keep_trim])
old_numerator_trim <- {
  lim_old <- quantile(num_trim, c(0.1, 0.9))
  keep_old <- num_trim >= lim_old[1] & num_trim <= lim_old[2]
  sum(num_trim[keep_old]) / sum(den_trim[keep_old])
}
stopifnot(abs(trimmed$r_av - expected_trim) < 1e-12,
          abs(trimmed$r_av - old_numerator_trim) > 1e-3)
ok("trim removes extreme ratio contributions, not merely extreme numerators")

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

## The design every later section uses.  make_weights_vqtl on a Normal baseline
## with a tau grid symmetric about 0.5 -- the one case the correction is exact
## for, and the one the guardrails let through.
TAUS <- seq(0.02, 0.98, length.out = 45)
zt <- qnorm(TAUS); sdY <- 1.3
Wv <- make_weights_vqtl(TAUS, q_tau = sdY * zt, mu = 0, sd = sdY)

P <- 50
fake <- list(
  W         = Wv,
  params    = rbind(beta_mu = a_hat[1:P], beta_sigma2 = v_hat[1:P]),
  SE_params = rbind(beta_mu = se_a[1:P], beta_sigma2 = se_v[1:P]),
  param_cov_upper = cbind(cov_beta_mu_beta_mu = se_a[1:P]^2,
                          cov_beta_mu_beta_sigma2 = cov_l[1:P],
                          cov_beta_sigma2_beta_sigma2 = se_v[1:P]^2))
df_report <- rep(37, P)
fake$param_ref_dist <- "t"
fake$param_ref_df <- rbind(beta_mu = df_report, beta_sigma2 = df_report)
fake$P_params <- rbind(beta_mu = 2 * pt(-abs(a_hat[1:P] / se_a[1:P]), df_report),
                       beta_sigma2 = 2 * pt(-abs(v_hat[1:P] / se_v[1:P]), df_report))
colnames(fake$params) <- colnames(fake$SE_params) <- paste0("rs", 1:P)
colnames(fake$param_ref_df) <- colnames(fake$P_params) <- paste0("rs", 1:P)
rep1 <- dispersion_report(fake, r_av = r_hat, se_rav = se_r, N = 461320, maf = 0.3)
stopifnot(nrow(rep1) == P)
stopifnot(all(c("beta_sigma2", "z_sigma2", "beta_sigma2_resid", "z_sigma2_resid",
                "G_skew_guard", "flag_skew") %in% names(rep1)))
ok("returns uncorrected and corrected columns side by side, plus the guard")
stopifnot(near(rep1$cov_mu_sigma2, cov_l[1:P], 1e-12))
ok("picks the mu/sigma2 covariance out of param_cov_upper by name")
stopifnot(near(rep1$p_sigma2,
               2 * pt(-abs(v_hat[1:P] / se_v[1:P]), df_report), 1e-15),
          near(rep1$p_sigma2_resid,
               2 * pt(-abs(rep1$z_sigma2_resid), df_report), 1e-15))
ok("dispersion_report reuses the fit's calibrated t reference for raw and residual p-values")
fake_no_ref <- fake
fake_no_ref$P_params <- NULL
fake_no_ref$param_ref_df <- NULL
fake_no_ref$param_ref_dist <- "none"
rep_no_ref <- dispersion_report(fake_no_ref, r_av = r_hat, N = 461320, maf = 0.3)
stopifnot(all(is.na(rep_no_ref$p_sigma2)),
          all(is.na(rep_no_ref$p_sigma2_resid)))
ok("dispersion_report omits p-values when the fit has no supported reference")
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
## TAUS / zt / sdY / Wv are defined above, with the dispersion_report fixture.
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
stopifnot(near(dispersion_skew_bias(0.2, 0.0, eaf = 0.01, eaf_orientation = "effect"), 0.04 * 0.98, 1e-12))
ok("bias = (1 - 2p) * beta_mu^2")
stopifnot(near(dispersion_skew_bias(0.2, 0.01, eaf = 0.01, eaf_orientation = "effect"),
               0.98 * (0.04 - 1e-4), 1e-12))
ok("debias subtracts Var(beta_mu_hat) from beta_mu_hat^2")
stopifnot(near(dispersion_skew_bias(0.2, 0.01, eaf = 0.01, eaf_orientation = "effect", debias = FALSE),
               0.98 * 0.04, 1e-12))
ok("debias = FALSE leaves the raw square")

stopifnot(abs(dispersion_skew_bias(5, 0.01, eaf = 0.5)) < 1e-15)
ok("exactly zero at EAF 0.5, for any beta_mu")

## a dispersion effect is invariant under an allele flip, and so is the correction:
## flipping sends p -> 1-p and beta -> -beta, and both b(tau) and the bias change sign.
stopifnot(near(dispersion_skew_bias(0.2, 0.01, eaf = 0.01, eaf_orientation = "effect"),
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
cc <- correct_dispersion_skew(bm, 0.05, sm, ss, cl, eaf = 0.01, eaf_orientation = "effect")
v_expect <- ss^2 + 4 * kk^2 * bm^2 * sm^2 + 2 * kk^2 * sm^4 - 4 * kk * bm * cl
stopifnot(near(cc$se_sigma2_skewcorr, sqrt(v_expect), 1e-14))
ok("SE is the stated delta-method variance including the within-variant covariance")
stopifnot(cc$se_sigma2_skewcorr > ss)
ok("the SE widens: the correction is estimated, and that is paid for")

## non-negative for any admissible covariance (|c| <= se_mu * se_sigma2)
cg <- seq(-1, 1, by = 0.05) * sm * ss
vv <- correct_dispersion_skew(rep(bm, length(cg)), rep(0.05, length(cg)),
                              sm, ss, cg, eaf = 0.01, eaf_orientation = "effect")
stopifnot(all(is.finite(vv$se_sigma2_skewcorr)), all(vv$se_sigma2_skewcorr > 0))
ok("variance stays positive across the whole admissible covariance range")

## length-1 arguments recycle against beta_mu
r1 <- correct_dispersion_skew(c(0.1, 0.2), 0.05, sm, ss, cl, eaf = 0.01, eaf_orientation = "effect")
r2 <- correct_dispersion_skew(c(0.1, 0.2), c(0.05, 0.05), c(sm, sm),
                              c(ss, ss), c(cl, cl), eaf = c(0.01, 0.01), eaf_orientation = "effect")
stopifnot(near(unlist(r1), unlist(r2), 1e-14), nrow(r1) == 2L)
ok("length-1 arguments recycle to length(beta_mu)")

## on a null variant the correction is a no-op to within its own noise
nn <- correct_dispersion_skew(0, 0.001, sm, ss, 0, eaf = 0.01, eaf_orientation = "effect")
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
cc <- correct_dispersion_skew(hat[1, ], hat[2, ], se1, se2, c12, eaf = 0.01, eaf_orientation = "effect")
r_cor <- mean(abs(cc$z_sigma2_skewcorr) > 1.959964)
cat(sprintf("   false dispersion rate at 0.05: %.4f uncorrected, %.4f corrected (nominal 0.05)\n",
            r_unc, r_cor))
stopifnot(r_unc > 0.20)                  # the artefact is there
stopifnot(abs(r_cor - 0.05) < 0.01)      # and the correction removes it
ok("a pure location effect's false dispersion rate returns to nominal")

## --- dispersion_report wiring ---------------------------------------------
## maf is FOLDED (0.3); eaf is the SIGNED effect-allele frequency (0.7 -- the
## effect allele is the major one here).  They must differ, and the orientation
## must be asserted, and W_sd must match the sd that built fake$W.
rep3 <- dispersion_report(fake, N = 461320, maf = 0.3, eaf = 0.7,
                          eaf_orientation = "effect", W_sd = sdY)
stopifnot(all(c("bias_skew", "beta_sigma2_skewcorr", "se_sigma2_skewcorr",
                "z_sigma2_skewcorr", "p_sigma2_skewcorr") %in% names(rep3)))
ok("dispersion_report adds the corrected columns when eaf is supplied")
stopifnot(near(rep3$beta_sigma2, v_hat[1:P], 1e-12))
ok("and still preserves the uncorrected estimate verbatim")
rep4 <- dispersion_report(fake, N = 461320, maf = 0.3)
stopifnot(!("beta_sigma2_skewcorr" %in% names(rep4)))
ok("no eaf supplied => no corrected columns invented")

## ===========================================================================
## GUARDRAILS.  Each block below must fail if its guardrail is removed.
## ===========================================================================

cat("\n== guardrail 1: EAF orientation ==\n")

errs <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
warns <- function(expr) {
  got <- FALSE
  withCallingHandlers(
    tryCatch(force(expr), error = function(e) NULL),
    warning = function(w) { got <<- TRUE; invokeRestart("muffleWarning") })
  got
}

## A realistic FOLDED MAF column: below 0.5 at every variant, many distinct
## values.  This is what a user gets from almost any frequency source, and it is
## the single most likely way to point the correction the wrong way.
maf_col <- seq(0.01, 0.49, length.out = 40)
eaf_col <- ifelse(seq_along(maf_col) %% 2 == 0, maf_col, 1 - maf_col)  # genuine, straddles
bm_v <- rep(0.2, length(maf_col)); sm_v <- rep(0.01, length(bm_v))
sv_v <- rep(0.012, length(bm_v)); cl_v <- rep(0, length(bm_v))

stopifnot(errs(correct_dispersion_skew(bm_v, 0.05, sm_v, sv_v, cl_v, eaf = maf_col)))
ok("refuses an eaf column that lies below 0.5 at EVERY variant (folded signature)")

stopifnot(errs(dispersion_skew_bias(bm_v, sm_v, eaf = maf_col)))
ok("dispersion_skew_bias refuses the same folded signature")

## the refusal is on the DATA, so asserting the orientation downgrades it to a
## loud warning rather than silently accepting it
stopifnot(warns(correct_dispersion_skew(bm_v, 0.05, sm_v, sv_v, cl_v,
                                        eaf = maf_col,
                                        eaf_orientation = "effect")))
ok("eaf_orientation = 'effect' downgrades the folded signature to a warning, not silence")

## a genuine EAF column straddles 0.5 and passes without a murmur
stopifnot(!warns(correct_dispersion_skew(bm_v, 0.05, sm_v, sv_v, cl_v,
                                         eaf = eaf_col)))
ok("a genuine EAF column straddling 0.5 passes silently")

## saying outright that the column is folded is always refused
stopifnot(errs(correct_dispersion_skew(bm_v, 0.05, sm_v, sv_v, cl_v,
                                       eaf = eaf_col, eaf_orientation = "folded")))
ok("eaf_orientation = 'folded' is refused outright")
stopifnot(errs(correct_dispersion_skew(bm_v, 0.05, sm_v, sv_v, cl_v,
                                       eaf = eaf_col, eaf_orientation = "yes")))
ok("an unrecognised eaf_orientation is refused")

## too few distinct sub-0.5 values to be conclusive: warn rather than accept
stopifnot(warns(correct_dispersion_skew(0.2, 0.05, 0.01, 0.012, 0, eaf = 0.01)))
ok("an unverifiable (scalar, sub-0.5) eaf warns when the orientation is not asserted")
stopifnot(!warns(correct_dispersion_skew(0.2, 0.05, 0.01, 0.012, 0, eaf = 0.01,
                                         eaf_orientation = "effect")))
ok("...and is silent once the orientation is asserted")

## WHY it is refused: a folded k does not merely fail to correct, it corrects
## the WRONG WAY and is worse than not correcting at all.
k_true <- 1 - 2 * eaf_col
k_fold <- 1 - 2 * pmin(eaf_col, 1 - eaf_col)      # what a folded column gives
delta_v <- 0.2
obs_v <- delta_v^2 * k_true                        # pure artefact, truth zero
z_unc  <- abs(obs_v / sv_v)
z_rgt  <- abs(correct_dispersion_skew(rep(delta_v, length(k_true)), obs_v, sm_v,
                                      sv_v, cl_v, debias = FALSE,
                                      mu3_over_var = k_true)$z_sigma2_skewcorr)
z_wrg  <- abs(correct_dispersion_skew(rep(delta_v, length(k_true)), obs_v, sm_v,
                                      sv_v, cl_v, debias = FALSE,
                                      mu3_over_var = k_fold)$z_sigma2_skewcorr)
flip <- k_true != k_fold                           # effect allele is the major one
cat(sprintf("   |Z| at the flipped variants: %.2f uncorrected, %.2f right way, %.2f folded\n",
            mean(z_unc[flip]), mean(z_rgt[flip]), mean(z_wrg[flip])))
stopifnot(all(z_rgt[flip] < 1e-6))                 # the right k removes it exactly
stopifnot(mean(z_wrg[flip]) > 1.8 * mean(z_unc[flip]))
ok("a folded k roughly doubles the false |Z| instead of removing it -- worse than no correction")

## dispersion_report takes BOTH frequencies, which is what invites the confusion
maf50 <- seq(0.05, 0.45, length.out = P)
stopifnot(errs(dispersion_report(fake, N = 461320, maf = maf50, eaf = maf50,
                                 eaf_orientation = "effect", W_sd = sdY)))
ok("dispersion_report refuses when eaf and maf are the SAME OBJECT")

maf_named <- stats::setNames(maf50, paste0("v", seq_len(P)))  # not identical(); same values
stopifnot(!identical(maf_named, maf50))
stopifnot(errs(dispersion_report(fake, N = 461320, maf = maf50, eaf = maf_named,
                                 eaf_orientation = "effect", W_sd = sdY)))
ok("dispersion_report refuses when eaf and maf are ELEMENT-WISE identical")

## a scalar maf recycled against a constant eaf is the same mistake
stopifnot(errs(dispersion_report(fake, N = 461320, maf = 0.3,
                                 eaf = rep(0.3, P),
                                 eaf_orientation = "effect", W_sd = sdY)))
ok("...including when one of them is a length-1 value recycled against the other")

## the orientation must be stated deliberately, with no working default
stopifnot(errs(dispersion_report(fake, N = 461320, maf = 0.3, eaf = 0.7,
                                 W_sd = sdY)))
ok("dispersion_report refuses an eaf with no eaf_orientation: there is no default")


cat("\n== guardrail 2: the W identity, not the column names ==\n")

## the shipped design passes
chk <- check_vqtl_design(Wv, sd = sdY)
stopifnot(isTRUE(chk$ok), chk$dev_mu < 1e-12, chk$dev_sigma2 < 1e-12,
          isTRUE(chk$taus_symmetric))
stopifnot(near(chk$sd_implied, sdY, 1e-9))
ok("check_vqtl_design accepts make_weights_vqtl's design and recovers its sd")

## the SAME Normal model built through make_weights_generic, parameterised in the
## SD rather than the variance.  Column names identical; numbers are z_tau, not
## z_tau/(2 sd).  A name check accepts it; so does any scale-free shape check,
## because it is still proportional to z_tau.
dist_cdf <- function(y, pr) pnorm(y, mean = pr$mu, sd = pr$sd)
dist_pdf <- function(y, pr) dnorm(y, mean = pr$mu, sd = pr$sd)
grad_mu  <- function(y, pr) -dnorm((y - pr$mu) / pr$sd) / pr$sd
grad_sd  <- function(y, pr) { z <- (y - pr$mu) / pr$sd; -(z * dnorm(z)) / pr$sd }
Wsd <- make_weights_generic(TAUS, q_tau = sdY * zt, dist_cdf, dist_pdf,
                            params = list(mu = 0, sd = sdY),
                            grad_funcs = list(beta_mu = grad_mu,
                                              beta_sigma2 = grad_sd))
stopifnot(identical(colnames(Wsd), colnames(Wv)))         # names agree...
stopifnot(near(Wsd[, "beta_sigma2"], zt, 1e-12))          # ...numbers do not
stopifnot(near(Wsd[, "beta_sigma2"] - 2 * sdY * Wv[, "beta_sigma2"], 0, 1e-12))
ok("an SD-parameterised Normal design carries the SAME column names, ratio 2*sd")
stopifnot(errs(check_vqtl_design(Wsd, sd = sdY)))
ok("check_vqtl_design refuses it even though the names match")

## and it refuses for the right reason: the message names the implied sd
m_sd <- tryCatch(check_vqtl_design(Wsd, sd = sdY), error = conditionMessage)
stopifnot(grepl("proportional to qnorm", m_sd), grepl("standard deviation", m_sd))
ok("...with a message that diagnoses the SD-vs-variance parameterisation")

## how much damage it would do: exactly the 1.6x sign flip the review measured
Ag   <- solve(crossprod(Wsd), t(Wsd))
kk01 <- 1 - 2 * 0.01
bias_b_sd <- (delta_v^2 / (2 * sdY)) * kk01 * zt          # bias in b(tau), as derived
proj_sd   <- drop(Ag %*% bias_b_sd)                       # what beta_sigma2 actually carries
stopifnot(abs(proj_sd[1]) < 1e-12)
stopifnot(near(proj_sd[2], delta_v^2 * kk01 / (2 * sdY), 1e-12))
wrong <- correct_dispersion_skew(delta_v, proj_sd[2], 1e-6, 0.012, 0,
                                 eaf = 0.01, eaf_orientation = "effect")
cat(sprintf("   truth 0: uncorrected %+.3e, 'corrected' on the wrong W %+.3e (%.2fx, sign flipped)\n",
            proj_sd[2], wrong$beta_sigma2_skewcorr,
            abs(wrong$beta_sigma2_skewcorr) / proj_sd[2]))
stopifnot(proj_sd[2] > 0, wrong$beta_sigma2_skewcorr < 0)
stopifnot(abs(abs(wrong$beta_sigma2_skewcorr) / proj_sd[2] - 1.6) < 0.01)
ok("on that W the correction would over-subtract by 2*sd: +1.51e-02 -> -2.41e-02")

## the 5-parameter mixture design, which passes a name check via mu_name
Wmix <- make_weights_normal_mixture_vqtl(TAUS, q_tau = sdY * zt, p1 = 0.5,
                                         mu1 = -1, sd1 = 1, mu2 = 1, sd2 = 1,
                                         density = "model")
stopifnot(all(c("beta_mu1", "beta_sigma2") %in% colnames(Wmix)))
stopifnot(errs(check_vqtl_design(Wmix, sd = sdY, mu_name = "beta_mu1")))
ok("check_vqtl_design refuses the 5-parameter mixture design under mu_name = 'beta_mu1'")

## sd is mandatory: without it the check is not identifiable
stopifnot(errs(check_vqtl_design(Wv, sd = NULL)))
stopifnot(errs(check_vqtl_design(Wv)))
ok("check_vqtl_design refuses without the sd that built W")
stopifnot(errs(check_vqtl_design(Wv, sd = 2 * sdY)))
ok("...and refuses the right W paired with the wrong sd")

## an asymmetric tau grid warns: the O(delta^3) orthogonality needs symmetry
TA <- seq(0.02, 0.70, length.out = 45)
Wa <- make_weights_vqtl(TA, q_tau = sdY * qnorm(TA), mu = 0, sd = sdY)
stopifnot(warns(check_vqtl_design(Wa, sd = sdY)))
ok("warns on a tau grid that is not symmetric about 0.5")
stopifnot(!warns(check_vqtl_design(Wv, sd = sdY)))
ok("...and does not warn on a symmetric one")

## dispersion_report gates on it
fake_bad <- fake; fake_bad$W <- Wsd
stopifnot(errs(dispersion_report(fake_bad, N = 461320, maf = 0.3, eaf = 0.7,
                                 eaf_orientation = "effect", W_sd = sdY)))
ok("dispersion_report refuses a fit whose W is merely NAMED beta_sigma2")

fake_noW <- fake; fake_noW$W <- NULL
stopifnot(errs(dispersion_report(fake_noW, N = 461320, maf = 0.3, eaf = 0.7,
                                 eaf_orientation = "effect", W_sd = sdY)))
ok("dispersion_report refuses a fit that carries no W at all")

stopifnot(errs(dispersion_report(fake, N = 461320, maf = 0.3, eaf = 0.7,
                                 eaf_orientation = "effect")))
ok("dispersion_report refuses without W_sd, the sd that built W")

## the opt-out is deliberate and loud, not silent
stopifnot(warns(dispersion_report(fake_bad, N = 461320, maf = 0.3, eaf = 0.7,
                                  eaf_orientation = "effect",
                                  check_design = FALSE)))
ok("check_design = FALSE proceeds but warns prominently")

## mu3_over_var is gated too -- it is the same correction by another route
stopifnot(errs(dispersion_report(fake_bad, N = 461320, maf = 0.3,
                                 mu3_over_var = 0.98, W_sd = sdY)))
ok("the mu3_over_var route is gated on the design as well")


cat("\n== NA corrections are counted, not returned silently ==\n")

## eaf of 0, 1 or NA gives no k
eaf_bad <- c(0.7, 0, 1, NA, 0.6)
nb <- length(eaf_bad)
cnt <- NULL
res_na <- withCallingHandlers(
  correct_dispersion_skew(rep(0.2, nb), 0.05, 0.01, 0.012, 0, eaf = eaf_bad),
  warning = function(w) { cnt <<- conditionMessage(w); invokeRestart("muffleWarning") })
nac <- attr(res_na, "na_counts")
stopifnot(!is.null(nac))
stopifnot(nac$n == nb, nac$n_eaf_na == 1L, nac$n_eaf_out_of_range == 2L,
          nac$n_k_unusable == 3L, nac$n_na == 3L)
ok("an eaf of 0, 1 or NA is counted on attr(, 'na_counts')")
stopifnot(!is.null(cnt), grepl("have NO corrected estimate", cnt))
stopifnot(sum(is.na(res_na$beta_sigma2_skewcorr)) == 3L)
ok("...and warns rather than returning the NAs silently")

## a cov_mu_sigma2 from a different source than the SEs breaches Cauchy-Schwarz
## and yields a negative delta-method variance
bmC <- 0.2; smC <- 0.01; kC <- 1 - 2 * 0.01
ssC <- 2 * kC * bmC * smC                # AM-GM equality: the tightest admissible point
clC <- 1.01 * smC * ssC                  # 1% past the Cauchy-Schwarz bound
cnt2 <- NULL
res_cs <- withCallingHandlers(
  correct_dispersion_skew(bmC, 0.05, smC, ssC, clC, eaf = 0.01,
                          eaf_orientation = "effect"),
  warning = function(w) { cnt2 <<- conditionMessage(w); invokeRestart("muffleWarning") })
nac2 <- attr(res_cs, "na_counts")
stopifnot(nac2$n_cov_cauchy_schwarz == 1L, nac2$n_var_negative == 1L,
          nac2$n_na == 1L, is.na(res_cs$se_sigma2_skewcorr))
ok("a covariance past |c| = se_mu * se_sigma2 gives a negative variance, and is counted")
stopifnot(!is.null(cnt2), grepl("Cauchy-Schwarz", cnt2))
ok("...and the warning names Cauchy-Schwarz and the mismatched covariance file")

## an admissible covariance does none of this
stopifnot(!warns(correct_dispersion_skew(bmC, 0.05, smC, ssC, 0.99 * smC * ssC,
                                         eaf = 0.01, eaf_orientation = "effect")))
ok("an admissible covariance produces no NA and no warning")

## and the counts survive the dispersion_report wrapper
eaf_rep <- c(0.7, 0, NA, rep(0.6, P - 3))
rep_na <- suppressWarnings(
  dispersion_report(fake, N = 461320, maf = 0.3, eaf = eaf_rep,
                    eaf_orientation = "effect", W_sd = sdY))
nac3 <- attr(rep_na, "na_counts")
stopifnot(!is.null(nac3), nac3$n == P, nac3$n_k_unusable == 2L, nac3$n_na == 2L)
ok("dispersion_report carries the na_counts attribute through")


cat("\nAll dispersion-guard tests passed.\n")
