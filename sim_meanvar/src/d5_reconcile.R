# d5_reconcile.R -- reconcile the chi2_4 sign disagreement with sim/dispersion.
# They report G is CONSERVATIVE on a standardised chi2_4 phenotype (over-predicts 1.6x).
# d4 found it ANTI-conservative (under-predicts 2.3-3.4x).  Three candidate causes:
#   (i)  the weight column: shipped z_tau/(2 sd) vs their (q_tau - m)/(2 v)   [tested in d4]
#   (ii) the law W is built from: the ALTERNATIVE marginal (what Stage 1 sees) vs the
#        NULL marginal (what fdr_exact.core uses)
#   (iii) OLS vs GLS -- their headline column is both.pooled.OLS.beta_sigma2.

.libPaths("/user/work/fh6520/r_packages")
source(Sys.getenv("ENGINE", "engine.R"))
options(width = 210, scipen = 5)
E <- function(...) gmoq_exact(..., ngrid = 200001)
NN <- 461320

cat("\n## Full factorial over the three candidate causes.\n")
cat("## Standardised chi2_4 phenotype, pure location effect, MAF 0.30.\n")
cat("## 'ratio' = exact |Z_sigma2| / G.  >1 anti-conservative, <1 conservative.\n\n")
R <- NULL
for (zt in c(8, 16, 30)) {
  dd <- zt / sqrt(NN * 2 * 0.3 * 0.7)
  for (wm in c("vqtl", "distfree")) for (wf in c("alt", "null")) for (fi in c("ols", "gls")) {
    e <- E(maf = 0.3, delta = dd, s = 1, base = make_base("chi2", 4),
           transform = "identity", N = NN, cov_mode = "pooled", fit = fi,
           W_mode = wm, W_from = wf)
    R <- rbind(R, data.frame(Z_mu_target = zt, W_col = wm, W_from = wf, est = fi,
                             absZ = abs(e$Z_sigma2), G = e$G, ratio = abs(e$Z_sigma2) / e$G))
  }
}
print(R, row.names = FALSE, digits = 4)

cat("\n## Same, for the NORMAL phenotype, as a control: every configuration should\n")
cat("## give ratio ~ 1 there, since that is where G was calibrated.\n\n")
R2 <- NULL
for (zt in c(8, 16, 30)) {
  dd <- zt / sqrt(NN * 2 * 0.3 * 0.7)
  for (wm in c("vqtl", "distfree")) for (wf in c("alt", "null")) for (fi in c("ols", "gls")) {
    e <- E(maf = 0.3, delta = dd, s = 1, base = make_base("normal"),
           transform = "identity", N = NN, cov_mode = "pooled", fit = fi,
           W_mode = wm, W_from = wf)
    R2 <- rbind(R2, data.frame(Z_mu_target = zt, W_col = wm, W_from = wf, est = fi,
                               absZ = abs(e$Z_sigma2), G = e$G, ratio = abs(e$Z_sigma2) / e$G))
  }
}
print(R2, row.names = FALSE, digits = 4)
cat(sprintf("\n  normal control: ratio range %.4f to %.4f across all 24 configurations\n",
            min(R2$ratio), max(R2$ratio)))
cat(sprintf("  chi2_4        : ratio range %.4f to %.4f\n", min(R$ratio), max(R$ratio)))
cat("\n  If NO configuration reproduces a ratio below 1 on chi2_4, the disagreement is\n")
cat("  not explained by any of (i)-(iii) and one of the two calculations is wrong.\n")

cat("\n\n## What does the r_av correction COST where there is no coupling to remove?\n")
cat("## (this is d4 section 4, rerun because the log was truncated)\n\n")
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
                             se_corrected_over_raw = se_d / e$se_rep[2],
                             absZ_true_disp_before = abs(e$Z_sigma2),
                             absZ_true_disp_after =
                               abs((e$theta[2] - rv * e$theta[1]) / se_d)))
}
print(S4, row.names = FALSE, digits = 4)
cat("\n  se ratio < 1: the correction gains power.  > 1: it costs power.\n")

cat("\n\n## Lack-of-fit Q at the effect size where the artefact FIRST reaches 5e-8.\n")
cat("## This is the practical question: does the model tournament warn you in time?\n\n")
S5 <- NULL
for (ss in c(0.10, 0.20, 0.35, 0.50)) {
  dd <- try(solve_delta(qnorm(1 - 5e-8/2), maf = 0.3, s = ss, mu0 = 0.5,
                        base = make_base("normal"), transform = "exp", N = NN,
                        cov_mode = "pooled", fit = "gls", W_mode = "vqtl",
                        ngrid = 200001, hi = 1), silent = TRUE)
  if (inherits(dd, "try-error") || is.na(dd)) next
  e <- E(maf = 0.3, delta = dd, s = ss, mu0 = 0.5, base = make_base("normal"),
         transform = "exp", N = NN, cov_mode = "pooled", fit = "gls", W_mode = "vqtl")
  S5 <- rbind(S5, data.frame(trait_skew = (exp(ss^2) + 2) * sqrt(exp(ss^2) - 1),
                             Z_mu = e$Z_mu, absZ_sig2 = abs(e$Z_sigma2),
                             Q_ncp = e$Q_ncp, Q_df = e$df_Q,
                             Q_power_at_0.05 = 1 - pchisq(qchisq(0.95, e$df_Q),
                                                          e$df_Q, ncp = e$Q_ncp)))
}
print(S5, row.names = FALSE, digits = 4)
cat("\n  Q_power_at_0.05 is the chance the lack-of-fit test flags the model at the very\n")
cat("  effect size where the dispersion artefact first becomes genome-wide significant.\n")
