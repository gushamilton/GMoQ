# d4_followup.R -- three open questions left by d2/d3.
#  1. On a SKEWED phenotype the G guard under-predicted by up to 6x (d2 part C3),
#     where sim/dispersion reports it over-predicts by 1.6x. Is the difference the
#     weight matrix (shipped Gaussian z/(2sd) vs their distribution-free (q-m)/(2v))?
#  2. Is there a residual artefact at MAF 0.5 on a skewed phenotype? If so it is a
#     THIRD mechanism, since genotype skewness is exactly zero there.
#  3. What does the r_av correction cost in power on a scale with NO coupling?

.libPaths("/user/work/fh6520/r_packages")
source(Sys.getenv("ENGINE", "engine.R"))
options(width = 210, scipen = 5)
E <- function(...) gmoq_exact(..., ngrid = 200001)
NN <- 461320

cat("\n############################################################################\n")
cat("## 1 + 2.  Skewed phenotype, PURE LOCATION effect, Z_mu = 30.\n")
cat("##   Wvqtl = shipped make_weights_vqtl, z_tau/(2 sd)  [what the pipeline uses]\n")
cat("##   Wdf   = distribution-free (q_tau - m)/(2 v)      [what sim/dispersion uses]\n")
cat("##   MAF 0.50 has genotype skewness exactly 0, so ANY signal there is a\n")
cat("##   mechanism neither study has named.\n")
cat("############################################################################\n")
S <- NULL
for (bs in list(make_base("normal"), make_base("chi2", 16), make_base("chi2", 8),
                make_base("chi2", 4), make_base("chi2", 2))) {
  for (mf in c(0.30, 0.50)) {
    dd <- 30 / sqrt(NN * 2 * mf * (1 - mf))
    a <- E(maf = mf, delta = dd, s = 1, base = bs, transform = "identity", N = NN,
           cov_mode = "pooled", fit = "gls", W_mode = "vqtl")
    b <- E(maf = mf, delta = dd, s = 1, base = bs, transform = "identity", N = NN,
           cov_mode = "pooled", fit = "gls", W_mode = "distfree")
    S <- rbind(S, data.frame(residual = bs$label, ph_skew = bs$skew, maf = mf,
                             gskew = a$gskew,
                             absZ_Wvqtl = abs(a$Z_sigma2), absZ_Wdf = abs(b$Z_sigma2),
                             G = a$G,
                             ratio_vqtl = abs(a$Z_sigma2) / ifelse(a$G == 0, NA, a$G),
                             ratio_df   = abs(b$Z_sigma2) / ifelse(a$G == 0, NA, a$G),
                             Q_ncp_vqtl = a$Q_ncp, Q_ncp_df = b$Q_ncp))
  }
}
print(S, row.names = FALSE, digits = 4)
cat("\n  If absZ_Wdf/G ~ 1 while absZ_Wvqtl/G >> 1, the excess is the shipped Gaussian W,\n")
cat("  not the mechanism, and the guard is anti-conservative for the SHIPPED pipeline.\n")
cat("  Any non-zero entry in the maf = 0.50 rows is a mechanism driven by PHENOTYPE\n")
cat("  skewness rather than genotype skewness.\n")

cat("\n## 1b.  Scaling of the MAF-0.50 skewed-phenotype term in the effect size.\n")
for (bs in list(make_base("chi2", 4), make_base("chi2", 2))) {
  ds <- 2^seq(log2(0.005), log2(0.08), length.out = 8)
  v <- vapply(ds, function(d) abs(E(maf = 0.5, delta = d, s = 1, base = bs,
                transform = "identity", N = NN, cov_mode = "pooled", fit = "gls",
                W_mode = "vqtl")$theta[2]), 0.0)
  cat(sprintf("   %-8s  W=vqtl  log-log slope in delta = %.3f\n",
              bs$label, coef(lm(log(v) ~ log(ds)))[2]))
  v2 <- vapply(ds, function(d) abs(E(maf = 0.5, delta = d, s = 1, base = bs,
                transform = "identity", N = NN, cov_mode = "pooled", fit = "gls",
                W_mode = "distfree")$theta[2]), 0.0)
  cat(sprintf("   %-8s  W=df    log-log slope in delta = %.3f\n",
              bs$label, coef(lm(log(v2) ~ log(ds)))[2]))
}

cat("\n\n############################################################################\n")
cat("## 3.  How skewed must a trait be before the mean-variance mechanism bites?\n")
cat("## Lognormal with log-scale sd s.  Reporting the mean effect |Z_mu| at which the\n")
cat("## spurious dispersion statistic reaches 5e-8, and r_av.\n")
cat("############################################################################\n")
S3 <- NULL
for (ss in c(0.10, 0.20, 0.35, 0.50, 0.70)) {
  dd <- try(solve_delta(qnorm(1 - 5e-8 / 2), maf = 0.3, s = ss, mu0 = 0.5,
                        base = make_base("normal"), transform = "exp", N = NN,
                        cov_mode = "pooled", fit = "gls", W_mode = "vqtl",
                        ngrid = 200001, hi = 1), silent = TRUE)
  if (inherits(dd, "try-error") || is.na(dd)) next
  e <- E(maf = 0.3, delta = dd, s = ss, mu0 = 0.5, base = make_base("normal"),
         transform = "exp", N = NN, cov_mode = "pooled", fit = "gls", W_mode = "vqtl")
  S3 <- rbind(S3, data.frame(s_logscale = ss,
                             trait_skew = (exp(ss^2) + 2) * sqrt(exp(ss^2) - 1),
                             r_av = e$theta[2] / e$theta[1],
                             Z_mu_at_5e8 = e$Z_mu,
                             beta_SD_per_allele = e$true_mu / e$sdY,
                             Q_ncp = e$Q_ncp, Q_df = e$df_Q))
}
print(S3, row.names = FALSE, digits = 4)
cat("\n  Q_ncp on Q_df df is the population lack-of-fit of the Gaussian vQTL model.\n")
cat("  Where Q_ncp is large the model tournament already refuses the fitted family.\n")

cat("\n\n############################################################################\n")
cat("## 4.  What does the r_av correction COST when there is no coupling to remove?\n")
cat("## Ratio of the residualised SE to the uncorrected SE, on each analysed scale.\n")
cat("## r_av is set to the value the trait actually implies on that scale.\n")
cat("############################################################################\n")
S4 <- NULL
for (tr in c("exp", "identity", "rint")) {
  ref <- E(maf = 0.3, delta = 0.005, s = 0.5, mu0 = 0.5, base = make_base("normal"),
           transform = tr, N = NN, cov_mode = "pooled", fit = "gls", W_mode = "vqtl")
  rv <- ref$theta[2] / ref$theta[1]
  e  <- E(maf = 0.3, delta = 0, gamma = 0.02, s = 0.5, mu0 = 0.5,
          base = make_base("normal"), transform = tr, N = NN,
          cov_mode = "pooled", fit = "gls", W_mode = "vqtl")
  Si <- solve((outer(e$taus, e$taus, pmin) - outer(e$taus, e$taus)) /
                outer(e$fy, e$fy) / (e$N * e$Vg))
  V <- solve(t(e$W) %*% Si %*% e$W)
  se_d <- sqrt(V[2,2] + rv^2 * V[1,1] - 2 * rv * V[1,2])
  S4 <- rbind(S4, data.frame(scale = c(exp="raw",identity="log",rint="RINT")[tr],
                             r_av_implied = rv,
                             corr_mu_sigma2 = V[1,2] / sqrt(V[1,1] * V[2,2]),
                             se_ratio_corrected_over_raw = se_d / e$se_rep[2],
                             absZ_true_disp_raw = abs(e$Z_sigma2),
                             absZ_true_disp_corrected =
                               abs((e$theta[2] - rv * e$theta[1]) / se_d)))
}
print(S4, row.names = FALSE, digits = 4)
cat("\n  se_ratio < 1 means the correction GAINS power; > 1 means it costs power.\n")
