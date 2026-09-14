# d2_characterise.R -- Task 1.  Exact characterisation of the spurious dispersion artefact.
# Everything here is a population calculation; there is no Monte Carlo error anywhere.

.libPaths("/user/work/fh6520/r_packages")
source(Sys.getenv("ENGINE", "engine.R"))
options(width = 210, scipen = 5)
ZGW <- qnorm(1 - 5e-8 / 2)   # 5.451

E <- function(...) gmoq_exact(..., ngrid = 200001)

cat("\n############################################################################\n")
cat("## A.  THE CLEAN 2x2:  the two mechanisms separated.\n")
cat("##\n")
cat("## DGP:  X | g ~ N(0.5 + delta*g, 0.5^2);  Y = T(X).  N = 200,000.\n")
cat("## 'raw'  = exp(X): a lognormal trait whose conditional SD is proportional to its\n")
cat("##          conditional mean -- a GENUINE mean-variance relationship.\n")
cat("## 'log'  = X itself: exactly homoskedastic normal, TRUE dispersion effect = 0.\n")
cat("## MAF 0.50 has genotype skewness EXACTLY zero, killing the skewness mechanism.\n")
cat("##\n")
cat("##   (log, MAF .50) = neither mechanism      (log, MAF .30) = skewness only\n")
cat("##   (raw, MAF .50) = mean-variance only     (raw, MAF .30) = both\n")
cat("############################################################################\n")

A <- NULL
for (tr in c("exp", "identity", "rint")) for (mf in c(0.30, 0.50)) {
  e <- E(maf = mf, delta = 0.10, s = 0.5, mu0 = 0.5, base = make_base("normal"),
         transform = tr, N = 200000, cov_mode = "pooled", fit = "gls", W_mode = "vqtl")
  A <- rbind(A, data.frame(
    scale = c(exp = "raw", identity = "log", rint = "RINT")[tr], maf = mf,
    gskew = e$gskew, Z_mu = e$Z_mu,
    true_sig2 = e$true_sigma2, est_sig2 = e$theta[2], bias = e$bias_sigma2,
    absZ_sig2 = abs(e$Z_sigma2), G_guard = e$G,
    p_sig2 = 2 * pnorm(-abs(e$Z_sigma2))))
}
print(A, row.names = FALSE, digits = 4)

cat("\n## A2.  How each mechanism scales in the effect size delta.\n")
cat("## log-log slope of |theta_sigma2| against delta (delta 0.01 -> 0.16):\n")
ds <- 2^seq(log2(0.01), log2(0.16), length.out = 9)
sl <- NULL
for (tr in c("exp", "identity")) for (mf in c(0.30, 0.50)) {
  v <- vapply(ds, function(d) abs(E(maf = mf, delta = d, s = 0.5, mu0 = 0.5,
                base = make_base("normal"), transform = tr, N = 200000,
                cov_mode = "pooled", fit = "gls", W_mode = "vqtl")$theta[2]), 0.0)
  co <- coef(lm(log(v) ~ log(ds)))
  sl <- rbind(sl, data.frame(scale = c(exp = "raw", identity = "log")[tr], maf = mf,
                             loglog_slope = co[2],
                             interpretation = ifelse(co[2] < 1.4, "LINEAR in delta",
                                              ifelse(co[2] < 2.4, "QUADRATIC in delta", "higher order"))))
}
print(sl, row.names = FALSE, digits = 4)
cat("## linear-in-delta  => genuine mean-variance relationship (Young et al. mechanism)\n")
cat("## quadratic        => genotype-skewness x density-curvature (sim/dispersion mechanism)\n")

cat("\n\n############################################################################\n")
cat("## B.  THE NUMBER FOR THE PAPER.\n")
cat("## Smallest PURE MEAN effect whose spurious dispersion statistic reaches\n")
cat("## genome-wide significance, |Z_sigma2| = 5.451 (p = 5e-8).\n")
cat("## beta_SD/allele is on the ANALYSED scale.  Z_mu is the mean effect's own Z.\n")
cat("############################################################################\n")

B <- NULL
for (NN in c(1e5, 2e5, 461320)) for (tr in c("exp", "identity", "rint")) for (mf in c(0.05, 0.30, 0.50)) {
  d <- try(solve_delta(ZGW, maf = mf, s = 0.5, mu0 = 0.5, base = make_base("normal"),
                       transform = tr, N = NN, cov_mode = "pooled", fit = "gls",
                       W_mode = "vqtl", ngrid = 200001, hi = 3), silent = TRUE)
  if (inherits(d, "try-error") || is.na(d)) { next }
  e <- E(maf = mf, delta = d, s = 0.5, mu0 = 0.5, base = make_base("normal"),
         transform = tr, N = NN, cov_mode = "pooled", fit = "gls", W_mode = "vqtl")
  B <- rbind(B, data.frame(N = NN, scale = c(exp="raw",identity="log",rint="RINT")[tr], maf = mf,
                           delta_logscale = d, beta_SD_per_allele = e$true_mu / e$sdY,
                           Z_mu_at_that_point = e$Z_mu, varexpl_pct = 100 * e$Z_mu^2 / NN,
                           G_guard = e$G))
}
print(B, row.names = FALSE, digits = 4)

cat("\n\n############################################################################\n")
cat("## C.  MAPS.  Pure LOCATION effect, homoskedastic normal (no mean-variance\n")
cat("## coupling at all), identity scale, N = 461,320.  So every number below is\n")
cat("## the SKEWNESS mechanism in isolation.\n")
cat("############################################################################\n")

cat("\n## C1.  MAF x Z_mu.  Entries are the spurious |Z_sigma2|.\n")
mafs <- c(0.005, 0.01, 0.05, 0.10, 0.30, 0.50); zs <- c(8, 16, 30, 60)
C1 <- matrix(NA, length(mafs), length(zs), dimnames = list(paste0("maf ", mafs), paste0("Zmu ", zs)))
for (i in seq_along(mafs)) for (j in seq_along(zs)) {
  dd <- zs[j] / sqrt(461320 * 2 * mafs[i] * (1 - mafs[i]))
  C1[i, j] <- abs(E(maf = mafs[i], delta = dd, s = 1, mu0 = 0, base = make_base("normal"),
                    transform = "identity", N = 461320, cov_mode = "pooled",
                    fit = "gls", W_mode = "vqtl")$Z_sigma2)
}
print(round(C1, 2))
cat("   (5.451 = genome-wide.  MAF 0.50 row should be ~0: genotype skewness is zero there.)\n")

cat("\n## C2.  N, at fixed Z_mu = 30, MAF 0.30.  Does more data help or hurt?\n")
C2 <- NULL
for (NN in c(5e4, 1e5, 2e5, 461320, 1e6)) {
  dd <- 30 / sqrt(NN * 2 * 0.3 * 0.7)
  e <- E(maf = 0.3, delta = dd, s = 1, base = make_base("normal"), transform = "identity",
         N = NN, cov_mode = "pooled", fit = "gls", W_mode = "vqtl")
  C2 <- rbind(C2, data.frame(N = NN, beta_SD_per_allele = dd, absZ_sig2 = abs(e$Z_sigma2),
                             G = e$G, ratio_to_sqrtN = abs(e$Z_sigma2) / sqrt(NN)))
}
print(C2, row.names = FALSE, digits = 4)

cat("\n## C3.  Phenotype shape (skewness / excess kurtosis of the residual),\n")
cat("##      MAF 0.30, Z_mu = 30, N = 461,320, identity scale.\n")
C3 <- NULL
for (bs in list(make_base("normal"), make_base("chi2", 16), make_base("chi2", 8),
                make_base("chi2", 4), make_base("chi2", 2), make_base("t", 8), make_base("t", 6))) {
  dd <- 30 / sqrt(461320 * 2 * 0.3 * 0.7)
  e <- try(E(maf = 0.3, delta = dd, s = 1, base = bs, transform = "identity",
             N = 461320, cov_mode = "pooled", fit = "gls", W_mode = "vqtl"), silent = TRUE)
  if (inherits(e, "try-error")) next
  C3 <- rbind(C3, data.frame(residual = bs$label, skew = bs$skew, exkurt = bs$kurt,
                             absZ_sig2 = abs(e$Z_sigma2), G = e$G,
                             ratio_obs_over_G = abs(e$Z_sigma2) / e$G))
}
print(C3, row.names = FALSE, digits = 4)
cat("   ratio_obs_over_G = 1 means the G guard (calibrated on the normal) is exactly right.\n")

cat("\n## C4.  Fitted model and estimator.  MAF 0.30, N = 461,320.\n")
cat("##      LS  = shipped make_weights_vqtl, [1, z_tau/(2 sd)]\n")
cat("##      TRU = distribution-free scale column [1, (q_tau - m)/(2 v)]\n")
cat("##      LSK = LS plus a kurtosis column (z^2 - 1)\n")
C4 <- NULL
for (sc in list(list(tr = "identity", lab = "pure location, normal (skewness mechanism)"),
                list(tr = "exp",      lab = "raw lognormal (mean-variance mechanism)"))) {
  for (mo in list(list(m = "vqtl", x = "none", n = "LS"),
                  list(m = "distfree", x = "none", n = "TRU"),
                  list(m = "vqtl", x = "kurt", n = "LSK"))) {
    for (fi in c("ols", "gls")) {
      dd <- if (sc$tr == "identity") 30 / sqrt(461320 * 2 * 0.3 * 0.7) else 0.10
      e <- E(maf = 0.3, delta = dd, s = if (sc$tr == "identity") 1 else 0.5,
             mu0 = if (sc$tr == "identity") 0 else 0.5,
             base = make_base("normal"), transform = sc$tr, N = 461320,
             cov_mode = "pooled", fit = fi, W_mode = mo$m, W_extra = mo$x)
      C4 <- rbind(C4, data.frame(scenario = sc$lab, model = mo$n, est = fi,
                                 Z_mu = e$Z_mu, absZ_sig2 = abs(e$Z_sigma2),
                                 bias = e$bias_sigma2, true_sig2 = e$true_sigma2,
                                 se_ratio = e$se_ratio[2],
                                 Q_ncp = e$Q_ncp, Q_df = e$df_Q))
    }
  }
}
print(C4, row.names = FALSE, digits = 4)
cat("\n   Q_ncp is the population lack-of-fit statistic against a chi-square on Q_df df.\n")
cat("   If Q_ncp is huge the model tournament already knows the fitted family is wrong.\n")
