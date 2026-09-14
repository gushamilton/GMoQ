# d3_ravstructure.R -- does Young et al.'s genome-wide residualisation have the right
# structural form for the mean-variance mechanism, and does it destroy genuine dispersion?
# Exact population calculation throughout.

.libPaths("/user/work/fh6520/r_packages")
source(Sys.getenv("ENGINE", "engine.R"))
options(width = 210, scipen = 5)
E <- function(...) gmoq_exact(..., ngrid = 200001)

cat("\n############################################################################\n")
cat("## 1.  IS THE COUPLING A SINGLE GENOME-WIDE CONSTANT?\n")
cat("## Young's d_l = alpha_vl - r_av * alpha_l assumes one slope r_av for all variants.\n")
cat("## Here: raw lognormal, pure mean effect.  If theta_sigma2/theta_mu is constant\n")
cat("## across MAF and effect size, a single r_av is the right model.\n")
cat("############################################################################\n")
R1 <- NULL
for (mf in c(0.01, 0.05, 0.2, 0.3, 0.5)) for (dd in c(0.01, 0.03, 0.10)) {
  e <- E(maf = mf, delta = dd, s = 0.5, mu0 = 0.5, base = make_base("normal"),
         transform = "exp", N = 461320, cov_mode = "pooled", fit = "gls", W_mode = "vqtl")
  R1 <- rbind(R1, data.frame(maf = mf, delta = dd, Z_mu = e$Z_mu,
                             theta_mu = e$theta[1], theta_sig2 = e$theta[2],
                             ratio_est = e$theta[2] / e$theta[1],
                             ratio_true = e$true_sigma2 / e$true_mu))
}
print(R1, row.names = FALSE, digits = 5)
cat(sprintf("\n  estimated ratio: mean %.5f, sd %.6f, relative spread %.3f%%\n",
            mean(R1$ratio_est), sd(R1$ratio_est), 100 * sd(R1$ratio_est) / mean(R1$ratio_est)))
cat(sprintf("  true ratio     : mean %.5f, sd %.6f\n", mean(R1$ratio_true), sd(R1$ratio_true)))
cat("  A near-constant estimated ratio => a single genome-wide r_av is well specified.\n")

cat("\n## 1b.  r_av is a property of the TRAIT, not a universal constant.\n")
cat("##      Varying the lognormal's log-scale sd s (i.e. the trait's own skewness):\n")
R1b <- NULL
for (ss in c(0.2, 0.35, 0.5, 0.7, 0.9)) {
  e <- E(maf = 0.3, delta = 0.02, s = ss, mu0 = 0.5, base = make_base("normal"),
         transform = "exp", N = 461320, cov_mode = "pooled", fit = "gls", W_mode = "vqtl")
  # skewness of a lognormal on the raw scale
  sk <- (exp(ss^2) + 2) * sqrt(exp(ss^2) - 1)
  R1b <- rbind(R1b, data.frame(s_logscale = ss, trait_skewness = sk,
                               r_av_est = e$theta[2] / e$theta[1],
                               r_av_true = e$true_sigma2 / e$true_mu))
}
print(R1b, row.names = FALSE, digits = 5)

cat("\n\n############################################################################\n")
cat("## 2.  DOES THE RESIDUALISATION ACTUALLY REMOVE IT, AND DOES IT SPARE A\n")
cat("##     GENUINE DISPERSION EFFECT?\n")
cat("## r_av is estimated ONCE (below, from a small-effect variant, standing in for a\n")
cat("## genome-wide LD-pruned fit) and then applied to every row.\n")
cat("############################################################################\n")
ref <- E(maf = 0.3, delta = 0.005, s = 0.5, mu0 = 0.5, base = make_base("normal"),
         transform = "exp", N = 461320, cov_mode = "pooled", fit = "gls", W_mode = "vqtl")
r_av <- ref$theta[2] / ref$theta[1]
cat(sprintf("\n  r_av estimated at a Z_mu = %.1f variant:  r_av = %.5f\n", ref$Z_mu, r_av))

R2 <- NULL
for (sc in list(list(d = 0.10, g = 0.00, lab = "PURE MEAN effect (truth: no dispersion)"),
                list(d = 0.03, g = 0.00, lab = "PURE MEAN effect, smaller"),
                list(d = 0.00, g = 0.02, lab = "PURE DISPERSION effect"),
                list(d = 0.00, g = 0.05, lab = "PURE DISPERSION effect, larger"),
                list(d = 0.10, g = 0.02, lab = "BOTH"))) {
  e <- E(maf = 0.3, delta = sc$d, gamma = sc$g, s = 0.5, mu0 = 0.5,
         base = make_base("normal"), transform = "exp", N = 461320,
         cov_mode = "pooled", fit = "gls", W_mode = "vqtl")
  # Young residual, and its SE including the within-variant covariance and r_av's own error.
  # Within-variant Cov(theta_mu, theta_sigma2) is what stage2's param_cov_upper carries.
  Si <- solve((outer(e$taus, e$taus, pmin) - outer(e$taus, e$taus)) /
                outer(e$fy, e$fy) / (e$N * e$Vg))
  Vth <- solve(t(e$W) %*% Si %*% e$W)
  d_l  <- e$theta[2] - r_av * e$theta[1]
  se_d <- sqrt(Vth[2,2] + r_av^2 * Vth[1,1] - 2 * r_av * Vth[1,2])
  R2 <- rbind(R2, data.frame(scenario = sc$lab, delta = sc$d, gamma = sc$g,
                             Z_mu = e$Z_mu,
                             absZ_raw = abs(e$Z_sigma2), absZ_resid = abs(d_l / se_d),
                             se_inflation = se_d / e$se_rep[2],
                             cov_mu_sig2 = Vth[1,2] / sqrt(Vth[1,1] * Vth[2,2])))
}
print(R2, row.names = FALSE, digits = 4)
cat("\n  absZ_raw vs absZ_resid: the correction should collapse the PURE MEAN rows to ~0\n")
cat("  and leave the PURE DISPERSION rows essentially untouched.\n")

cat("\n\n############################################################################\n")
cat("## 3.  DOES THE RESIDUALISATION TOUCH THE SKEWNESS MECHANISM?  (It should not.)\n")
cat("## Pure location effect, homoskedastic normal, identity scale: there is no\n")
cat("## mean-variance coupling at all, so the true r_av here is ZERO and the\n")
cat("## residualisation is a no-op.  The artefact survives it.\n")
cat("############################################################################\n")
R3 <- NULL
for (mf in c(0.01, 0.05, 0.3, 0.5)) for (zt in c(16, 30)) {
  dd <- zt / sqrt(461320 * 2 * mf * (1 - mf))
  e <- E(maf = mf, delta = dd, s = 1, mu0 = 0, base = make_base("normal"),
         transform = "identity", N = 461320, cov_mode = "pooled", fit = "gls", W_mode = "vqtl")
  R3 <- rbind(R3, data.frame(maf = mf, Z_mu_target = zt,
                             theta_mu = e$theta[1], theta_sig2 = e$theta[2],
                             implied_rav = e$theta[2] / e$theta[1],
                             absZ_sig2 = abs(e$Z_sigma2), G_guard = e$G))
}
print(R3, row.names = FALSE, digits = 4)
cat("\n  'implied_rav' is NOT constant here -- it scales with delta and with MAF through\n")
cat("  the genotype skewness.  A single genome-wide slope cannot absorb it, which is\n")
cat("  exactly why Young's correction is the wrong tool for this second mechanism.\n")
