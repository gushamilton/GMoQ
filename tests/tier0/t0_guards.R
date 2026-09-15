#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 0 -- the two dispersion guards, against their published closed forms.
#
# Two independent mechanisms turn a large MEAN effect into a false dispersion
# finding, and they need different tools (stage2/R/dispersion-guards.R):
#
#   (a) mean-variance coupling on a skewed analysed scale -- LINEAR in the mean
#       effect, common slope across variants, removed exactly by log/RINT.
#       Tool: estimate_rav() + residualise_dispersion().
#   (b) genotype skewness x density curvature -- QUADRATIC in the mean effect,
#       per-variant, zero at MAF 0.5, survives every phenotype transform.
#       Tool: dispersion_skew_guard(), which flags and does not correct.
#
# This file checks the arithmetic of both. Whether the ARTEFACTS themselves are
# still present at the documented magnitudes is tier 2 (T2.PIN.*): that needs
# simulation, not algebra.
#
# Run:  Rscript --vanilla tests/tier0/t0_guards.R
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier0/t0_guards.R")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
gmoq_begin("tier0", "guards")

set.seed(20260914)
N <- 461320

# =====================================================================
# 1. dispersion_skew_guard:  G = 0.6855 * gskew * Z_mu^2 / sqrt(N)
# =====================================================================
p <- c(0.005, 0.01, 0.05, 0.1, 0.3, 0.5)
gskew <- (1 - 2 * p) / sqrt(2 * p * (1 - p))
G <- dispersion_skew_guard(z_mu = rep(10, length(p)), N = N, maf = p)
gmoq_close("T0.G.closed_form",
           "G is 0.6855 * gskew * Z_mu^2 / sqrt(N) with gskew = (1-2p)/sqrt(2p(1-p)).",
           G, 0.6855 * gskew * 100 / sqrt(N), 1e-14)

gmoq_check("T0.G.zero_at_maf_half",
           "G is exactly zero at MAF 0.5 for any Z_mu: the driver is skewness, not effect size.",
           abs(dispersion_skew_guard(z_mu = 200, N = N, maf = 0.5)) < 1e-12,
           criterion = "|G| < 1e-12 at MAF 0.5, Z_mu = 200")

gmoq_check("T0.G.quadratic_in_z",
           "G scales as Z_mu^2, the quadratic signature of mechanism (b).",
           abs(dispersion_skew_guard(20, N, maf = 0.1) /
                 dispersion_skew_guard(10, N, maf = 0.1) - 4) < 1e-12,
           criterion = "G(2Z)/G(Z) == 4")
gmoq_check("T0.G.shrinks_as_sqrt_N",
           "At fixed Z_mu, G falls as 1/sqrt(N); at fixed per-allele effect it grows as sqrt(N).",
           abs(dispersion_skew_guard(10, 4 * N, maf = 0.1) /
                 dispersion_skew_guard(10, N, maf = 0.1) - 0.5) < 1e-12,
           criterion = "G(4N)/G(N) == 1/2 at fixed Z")

# recovery of |gskew| from d_resid when no allele frequency is available
d_resid <- 2 * p * (1 - p) * N
Gn <- dispersion_skew_guard(z_mu = rep(10, length(p)), N = N, n_eff = d_resid)
gmoq_close("T0.G.recovered_from_d_resid",
           "Without a MAF, |gskew| is recovered from d_resid/N = 2p(1-p); only the sign is lost.",
           abs(Gn), abs(G), 1e-9)
gmoq_refuses("T0.G.requires_inputs",
             "The guard refuses rather than guessing when neither maf nor (n_eff, N) is given.",
             dispersion_skew_guard(z_mu = 10), "Supply either maf")

# The published values, reproduced. The Monte Carlo table in
# docs/FINDINGS_pooled_dispersion.md section 1 is at N = 100,000 (its separate
# "bigN" rows are the N = 461,320 ones), which is what the formula reproduces:
#   MAF 0.30, Z_mu 40.1 -> G = 2.155;  MAF 0.01, Z_mu 8.7 -> G = 1.146;
#   MAF 0.50, Z_mu 73.9 -> G = 0.000 exactly.
N_MC <- 100000
gmoq_check("T0.G.reproduces_documented_cells",
           "G reproduces the published values for the three decisive cells.",
           abs(dispersion_skew_guard(40.1, N_MC, maf = 0.30) - 2.155) < 0.02 &&
             abs(dispersion_skew_guard(8.7, N_MC, maf = 0.01) - 1.146) < 0.02 &&
             abs(dispersion_skew_guard(73.9, N_MC, maf = 0.50)) < 1e-10,
           stat = dispersion_skew_guard(40.1, N_MC, maf = 0.30),
           criterion = "G(40.1, MAF 0.30) = 2.155; G(8.7, MAF 0.01) = 1.146; G(73.9, MAF 0.50) = 0, all +- 0.02")

# =====================================================================
# 2. estimate_rav: both moment corrections, against a known truth
#
# Reproduces docs/MEAN_VARIANCE_ARTEFACT.md section 5: with true r_av = 0.135,
# a naive OLS of alpha_v on alpha is biased LOW by regression dilution, the
# dilution correction alone is biased HIGH by the within-variant sampling
# covariance, and only both together recover the truth.
# =====================================================================
n_var <- 50000
r_true <- 0.135

#' Build a synthetic scan with a known coupling and a known sampling covariance.
make_scan <- function(sd_alpha, s_a, s_v, rho_within, seed) {
  set.seed(seed)
  alpha <- rnorm(n_var, 0, sd_alpha)
  d_true <- rnorm(n_var, 0, sd_alpha / 5)
  c_l <- rho_within * s_a * s_v
  L <- chol(matrix(c(s_a^2, c_l, c_l, s_v^2), 2, 2))
  e <- matrix(rnorm(2 * n_var), n_var, 2) %*% L
  list(bm = alpha + e[, 1], bv = r_true * alpha + d_true + e[, 2],
       s_a = rep(s_a, n_var), s_v = rep(s_v, n_var), c_l = rep(c_l, n_var))
}

# (A) REGRESSION DILUTION, isolated: no within-variant covariance at all, and a
# sampling error comparable with the signal. A naive regression must come out
# too small, and the denominator correction must fix it.
A_ <- make_scan(sd_alpha = 0.02, s_a = 0.02, s_v = 0.02, rho_within = 0, seed = 1)
naive_A <- unname(stats::coef(stats::lm(A_$bv ~ A_$bm))[2])
fixed_A <- estimate_rav(A_$bm, A_$bv, A_$s_a, A_$s_v,
                        cov_mu_sigma2 = A_$c_l, trim = 0)$r_av
gmoq_check("T0.RAV.naive_is_diluted",
           "A naive regression of the dispersion effect on the mean effect is biased LOW by dilution.",
           naive_A < r_true * 0.8, stat = naive_A,
           criterion = sprintf("naive < 0.8 x %.3f (attenuation by Var/(Var+s_a^2) = 0.5)", r_true))
gmoq_check("T0.RAV.dilution_correction_works",
           "The denominator correction -s_a^2 removes the dilution and recovers the truth.",
           abs(fixed_A - r_true) < 0.01, stat = fixed_A,
           criterion = sprintf("|r_av - %.3f| < 0.01", r_true))

# (B) WITHIN-VARIANT SAMPLING COVARIANCE, isolated: negligible dilution, but the
# two parameter estimates are correlated within a variant. Correcting only the
# denominator must now be biased HIGH, and only subtracting c_l fixes it.
B_ <- make_scan(sd_alpha = 0.05, s_a = 0.005, s_v = 0.02, rho_within = 0.6, seed = 2)
dilution_only <- sum(B_$bm * B_$bv) / sum(B_$bm^2 - B_$s_a^2)
both <- estimate_rav(B_$bm, B_$bv, B_$s_a, B_$s_v,
                     cov_mu_sigma2 = B_$c_l, trim = 0)$r_av
gmoq_check("T0.RAV.dilution_only_is_biased_high",
           "Correcting dilution alone is biased HIGH: the within-variant covariance must go too.",
           dilution_only > r_true * 1.05, stat = dilution_only,
           criterion = "dilution-only > 1.05 x true r_av")
gmoq_check("T0.RAV.both_corrections_recover_truth",
           "With both corrections the moment estimator recovers the true coupling.",
           abs(both - r_true) < 0.005, stat = both,
           criterion = sprintf("|r_av - %.3f| < 0.005", r_true))

# The trimming default costs a little accuracy; the documentation says 0.01 is
# safe and 0.05 already too aggressive. Measured here rather than assumed.
tr_bias <- vapply(c(0, 0.001, 0.01, 0.05), function(tt)
  estimate_rav(A_$bm, A_$bv, A_$s_a, A_$s_v, cov_mu_sigma2 = A_$c_l, trim = tt)$r_av - r_true,
  numeric(1))
gmoq_check("T0.RAV.trim_cost",
           "The default trim = 0.01 biases r_av by less than 10 per cent; trim = 0.05 is worse.",
           abs(tr_bias[3]) < 0.1 * r_true && abs(tr_bias[4]) >= abs(tr_bias[3]),
           stat = tr_bias[3],
           criterion = "|bias at trim=0.01| < 0.0135 and monotone in trim")

bm <- A_$bm; bv <- A_$bv; s_a <- A_$s_a; s_v <- A_$s_v; c_l <- A_$c_l

gmoq_check("T0.RAV.warns_without_covariance",
           "Omitting the within-variant covariance warns: it is an assumption, not a default.",
           tryCatch({ estimate_rav(bm, bv, s_a, s_v); FALSE },
                    warning = function(w) grepl("cov_mu_sigma2 not supplied", conditionMessage(w))),
           criterion = "a warning is raised naming cov_mu_sigma2")
gmoq_refuses("T0.RAV.refuses_too_few_variants",
             "Fewer than 20 usable variants is refused rather than producing a number.",
             estimate_rav(bm[1:5], bv[1:5], s_a[1:5], s_v[1:5], cov_mu_sigma2 = c_l[1:5]),
             "Fewer than 20 usable variants")

jk <- estimate_rav(bm, bv, s_a, s_v, cov_mu_sigma2 = c_l,
                   block = rep(1:22, length.out = n_var))
gmoq_check("T0.RAV.jackknife_se_present",
           "A block vector with >= 10 blocks produces a delete-one-block SE.",
           is.finite(jk$se) && jk$n_blocks == 22L,
           stat = jk$se, criterion = "finite SE, 22 blocks")
gmoq_check("T0.RAV.jackknife_se_needs_ten_blocks",
           "Fewer than 10 blocks yields NA rather than an SE from too few replicates.",
           is.na(estimate_rav(bm, bv, s_a, s_v, cov_mu_sigma2 = c_l,
                              block = rep(1:5, length.out = n_var))$se),
           criterion = "SE is NA below 10 blocks")

# =====================================================================
# 3. residualise_dispersion: the delta-method variance, written out
# =====================================================================
rd <- residualise_dispersion(bm[1:5], bv[1:5], s_a[1:5], s_v[1:5], c_l[1:5],
                             r_av = 0.699, se_rav = 0.03)
gmoq_close("T0.RAV.residual_estimate",
           "The corrected estimate is alpha_v - r_av alpha.",
           rd$beta_sigma2_resid, bv[1:5] - 0.699 * bm[1:5], 1e-14)
v_hand <- s_v[1:5]^2 + 0.699^2 * s_a[1:5]^2 - 2 * 0.699 * c_l[1:5] + 0.03^2 * bm[1:5]^2
gmoq_close("T0.RAV.residual_variance",
           "Its variance is s_v^2 + r^2 s_a^2 - 2 r c_l + Var(r) alpha^2, including the r_av term.",
           rd$se_sigma2_resid, sqrt(v_hand), 1e-14)
rd0 <- residualise_dispersion(bm[1:5], bv[1:5], s_a[1:5], s_v[1:5], c_l[1:5], r_av = 0.699)
gmoq_check("T0.RAV.se_rav_default_zero",
           "se_rav defaults to 0, reproducing the uncorrected-r_av version of the published method.",
           all(rd0$se_sigma2_resid <= rd$se_sigma2_resid),
           criterion = "se_rav = 0 gives the smaller SE")
gmoq_check("T0.RAV.noop_when_rav_zero",
           "With r_av = 0 the correction is exactly a no-op on the point estimate.",
           max(abs(residualise_dispersion(bm[1:5], bv[1:5], s_a[1:5], s_v[1:5], c_l[1:5],
                                          r_av = 0)$beta_sigma2_resid - bv[1:5])) < 1e-15,
           criterion = "r_av = 0 leaves beta_sigma2 unchanged")

# =====================================================================
# 4. dispersion_report: the wiring, and the threshold it hard-codes
# =====================================================================
Tt <- 9L; taus <- seq(0.1, 0.9, length.out = Tt)
Wv <- make_weights_vqtl(taus, qnorm(taus), mu = 0, sd = 1)
P <- 6L
Sig <- diag(1e-4, Tt)
upper_rowmajor <- function(S) { n <- nrow(S); o <- numeric(n*(n+1)/2); i <- 1L
  for (r in seq_len(n)) for (cc in r:n) { o[i] <- S[r, cc]; i <- i + 1L }; o }
n_elem <- Tt * (Tt + 1L) / 2L
bmat <- matrix(rnorm(Tt * P, 0, 0.001), Tt, P,
               dimnames = list(paste0("tau", taus), paste0("v", seq_len(P))))
st1 <- list(taus = taus, q_tau = qnorm(taus), Q_slope = bmat,
            SE_tau = matrix(sqrt(diag(Sig)), Tt, P))
cdat <- list(cov_vec = rep(upper_rowmajor(Sig), P),
             offsets = seq(0, by = n_elem * 4, length.out = P))
fit <- param_gwas(st1, transform = "custom_W", transform_args = list(W = Wv),
                  se_mode = "tau_cov", tau_cov_estimator = "ols", cov_data = cdat,
                  return_cov = "upper", q_method = "chisq_naive",
                  nu = 900, nu_source = "test")
rep_df <- dispersion_report(fit, r_av = 0.135, se_rav = 0.01, N = N, maf = 0.05)
gmoq_check("T0.G.report_columns",
           "dispersion_report returns the uncorrected, residualised and guard statistics side by side.",
           all(c("z_sigma2", "beta_sigma2_resid", "G_skew_guard", "flag_skew") %in% names(rep_df)),
           criterion = "all three families of column present")
gmoq_close("T0.G.report_guard_matches",
           "The report's G column is dispersion_skew_guard() applied to its own z_mu.",
           rep_df$G_skew_guard,
           dispersion_skew_guard(rep_df$z_mu, N = N, maf = 0.05), 1e-14)
gmoq_check("T0.G.report_flag_threshold",
           "The shipped flag fires at |G| >= 0.2, the threshold calibrated for OLS on a normal trait.",
           identical(rep_df$flag_skew, !is.na(rep_df$G_skew_guard) &
                       abs(rep_df$G_skew_guard) >= 0.2),
           criterion = "flag_skew == (|G| >= 0.2)")
gmoq_check("T0.G.threshold_is_not_gls_calibrated",
           "Recorded, not a pass/fail: 0.2 is the OLS/normal threshold. docs/MEAN_VARIANCE_ARTEFACT.md says GLS on a skewed trait needs 0.06, and the package hard-codes only 0.2.",
           TRUE, stat = 0.2,
           criterion = "known gap, reported in docs/TEST_SUITE.md",
           detail = "dispersion_report() offers no way to set the threshold")

gmoq_finish()
