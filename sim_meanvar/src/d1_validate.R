# d1_validate.R -- two independent validations of engine.R.
#
# V1  Against sim/dispersion's exact calculator (a separate implementation, in Python)
#     on the cells it covers: normal baseline, pure location effect, no transform.
# V2  Against the Monte Carlo in /user/work/fh6520/scratch_rel/pc.R section E, which is
#     the source of the three numbers in docs/RELATEDNESS.md.

.libPaths("/user/work/fh6520/r_packages")
source(Sys.getenv("ENGINE", "engine.R"))
options(width = 200, scipen = 6)

EX <- Sys.getenv("EXACTDIR", "exact")

cat("==========================================================================\n")
cat("V1  engine.R (R)  vs  sim/dispersion fdr_exact.py (Python), independent code\n")
cat("    normal baseline, pure location effect, identity transform.\n")
cat("    Their design_W is (q-m)/(2v) built at the NULL cell -> W_mode='distfree',\n")
cat("    W_from='null'.  Comparing theta_pop for beta_sigma2 and the pooled se_ratio.\n")
cat("==========================================================================\n")

cells <- expand.grid(maf = c(0.005, 0.01, 0.05, 0.2, 0.3, 0.5),
                     N = c(50000, 461320), stringsAsFactors = FALSE)
res <- list()
for (i in seq_len(nrow(cells))) {
  f <- file.path(EX, sprintf("z_normal_maf%s_N%d.tsv", cells$maf[i], cells$N[i]))
  if (!file.exists(f)) { cat("MISSING", f, "\n"); next }
  d <- read.delim(f, check.names = FALSE)
  d <- d[d$kind == "loc", ]
  for (r in seq_len(nrow(d))) {
    e <- gmoq_exact(maf = cells$maf[i], delta = d$delta[r], gamma = 0, s = 1,
                    base = make_base("normal"), transform = "identity",
                    N = cells$N[i], cov_mode = "pooled", fit = "gls",
                    W_mode = "distfree", W_from = "null", ngrid = 60001)
    res[[length(res) + 1]] <- data.frame(
      maf = cells$maf[i], N = cells$N[i], target_z = d$target_z[r], delta = d$delta[r],
      py_theta2  = d[["both.pooled.gls.beta_sigma2.theta_pop"]][r],
      R_theta2   = e$theta[2],
      py_zrep    = d[["both.pooled.gls.beta_sigma2.z_rep"]][r],
      R_zrep     = e$Z_sigma2,
      py_seratio = d[["both.pooled.gls.beta_sigma2.se_ratio"]][r],
      R_seratio  = e$se_ratio[2],
      py_gskew   = d$geno_skew[r], R_gskew = e$gskew)
  }
}
res <- do.call(rbind, res)
res$rel_theta2 <- (res$R_theta2 - res$py_theta2) / abs(res$py_theta2)
res$rel_zrep   <- (res$R_zrep   - res$py_zrep)   / abs(res$py_zrep)
res$rel_seratio<- (res$R_seratio- res$py_seratio)/ abs(res$py_seratio)

cat(sprintf("\n n cells compared: %d\n", nrow(res)))
cat(sprintf(" max |rel err| theta_sigma2 : %.3e\n", max(abs(res$rel_theta2))))
cat(sprintf(" max |rel err| z_rep        : %.3e\n", max(abs(res$rel_zrep))))
cat(sprintf(" max |rel err| se_ratio     : %.3e\n", max(abs(res$rel_seratio))))
cat(sprintf(" max |rel err| geno_skew    : %.3e\n",
            max(abs((res$R_gskew - res$py_gskew) / res$py_gskew))))
cat("\n worst 6 rows by |rel err| in z_rep:\n")
print(head(res[order(-abs(res$rel_zrep)),
               c("maf","N","target_z","py_zrep","R_zrep","rel_zrep")], 6), row.names = FALSE)
cat("\n a representative slice (maf 0.30, N 461320):\n")
print(res[res$maf == 0.3 & res$N == 461320,
          c("target_z","delta","py_theta2","R_theta2","py_zrep","R_zrep")], row.names = FALSE)

cat("\n\n==========================================================================\n")
cat("V2  The three numbers in docs/RELATEDNESS.md, reproduced exactly.\n")
cat("    pc.R section E DGP:  y0 = exp(0.5 + 0.10*g + N(0, 0.5^2)),  g~Bin(2,0.3), N=2e5.\n")
cat("    'bias in SD units' there = mean(theta_sigma2)/sd(theta_sigma2) over 20 reps,\n")
cat("    i.e. the per-replicate |Z| of the spurious dispersion estimate.\n")
cat("    NOTE pc.R residualises on an intercept only and uses the ALTERNATIVE marginal\n")
cat("    for W (it computes sd(y) from the simulated y), so W_from='alt'.\n")
cat("==========================================================================\n")

tab <- NULL
for (tr in c("exp", "identity", "rint")) {
  e <- gmoq_exact(maf = 0.3, delta = 0.10, gamma = 0, s = 0.5, mu0 = 0.5,
                  base = make_base("normal"), transform = tr, N = 200000,
                  cov_mode = "pooled", fit = "gls", W_mode = "vqtl", W_from = "alt")
  tab <- rbind(tab, data.frame(
    scale = c(exp = "raw (lognormal)", identity = "log", rint = "RINT")[tr],
    theta_sigma2 = e$theta[2], true_sigma2 = e$true_sigma2,
    bias = e$bias_sigma2, se_rep = e$se_rep[2], sd_true = e$sd_true[2],
    Z_sigma2 = e$Z_sigma2, ncp_bias = e$ncp_sigma2, Z_mu = e$Z_mu,
    G = e$G, se_ratio = e$se_ratio[2]))
}
print(tab, row.names = FALSE, digits = 5)
cat("\n  doc claims:  raw 28.1 SD,  log 2.7 SD,  RINT 3.2 SD\n")
cat("  pc.R MC used 20 replicates, so the MC error on each of those ratios is\n")
cat(sprintf("  about 1/sqrt(2*(20-1)) = %.3f relative, i.e. +-%.2f on 28.1 and +-%.2f on 3.0.\n",
            1 / sqrt(2 * 19), 28.1 / sqrt(2 * 19), 3.0 / sqrt(2 * 19)))
cat("  (the ratio mean/sd is itself estimated; its relative SE is ~1/sqrt(2(n-1)) for the sd\n")
cat("   plus 1/(sqrt(n)*ratio) for the mean, the former dominating here.)\n")
