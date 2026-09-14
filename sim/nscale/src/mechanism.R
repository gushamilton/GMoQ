## Why GLS coverage falls with N under a true scale effect, and OLS does not.
##
## Hypothesis: the misspecification residual beta_tau - W theta is a property of
## the DGP, FIXED in absolute units and independent of N.  The per-tau SE shrinks
## like N^-1/2.  So lack of fit measured in SE units grows like N^1/2, and any
## bias it induces, measured in units of the estimator's own sampling sd, grows
## the same way.  Coverage then falls monotonically in N.
args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]; outdir <- args[2]
.libPaths(c(LIB)); stopifnot(identical(.libPaths()[1], LIB))
tab <- read.delim(file.path(outdir, "summary.tsv"), stringsAsFactors = FALSE)

cat("=== lack of fit is N-relative: same absolute residual, shrinking SE_tau ===\n")
cat(sprintf("%-6s %-7s %-5s %12s %12s %12s\n", "arm", "dgp", "fit",
            "N=30000", "N=1e+05", "N=461320"))
for (arm in c("A", "B", "C")) for (dg in c("null", "scale")) for (fi in c("ADD", "LS", "TRU")) {
  s <- tab[tab$arm == arm & tab$dgp == dg & tab$fit == fi & tab$est == "ols", ]
  if (!nrow(s)) next
  cat(sprintf("%-6s %-7s %-5s", arm, dg, fi))
  for (N in c(30000, 1e5, 461320))
    cat(sprintf(" %12.3f", mean(s[s$N == N, "max_resid_ols_per_SEtau"])))
  cat("   (abs resid:")
  for (N in c(30000, 1e5, 461320))
    cat(sprintf(" %.5f", mean(s[s$N == N, "max_resid_ols"])))
  cat(")\n")
}

cat("\n=== bias in units of the estimator OWN sampling sd, dgp = scale ===\n")
cat(sprintf("%-6s %-5s %-12s %-5s %10s %10s %10s\n", "arm", "fit", "param", "est",
            "N=30000", "N=1e+05", "N=461320"))
sc <- tab[tab$dgp == "scale", ]
for (arm in c("A", "B", "C"))
  for (k in unique(paste(sc$fit, sc$param))) {
    fi <- sub(" .*", "", k); pm <- sub(".* ", "", k)
    for (est in c("ols", "gls")) {
      s <- sc[sc$arm == arm & sc$fit == fi & sc$param == pm & sc$est == est, ]
      if (!nrow(s)) next
      v <- vapply(c(30000, 1e5, 461320),
                  function(N) mean(abs(s[s$N == N, "bias_sd"])), 0)
      if (max(v, na.rm = TRUE) < 0.05) next
      cat(sprintf("%-6s %-5s %-12s %-5s %10.3f %10.3f %10.3f\n",
                  arm, fi, pm, est, v[1], v[2], v[3]))
    }
  }

cat("\n=== does |bias_sd| grow like sqrt(N)?  ratio to the N=30,000 value ===\n")
cat("   (sqrt(N/30000) = 1.00, 1.83, 3.92 at the three rungs)\n")
for (arm in c("A", "B", "C")) {
  s <- sc[sc$arm == arm & sc$est == "gls" & sc$fit == "LS" & sc$param == "beta_sigma2", ]
  if (!nrow(s)) next
  v <- vapply(c(30000, 1e5, 461320), function(N) mean(abs(s[s$N == N, "bias_sd"])), 0)
  cat(sprintf("arm %s  GLS LS:beta_sigma2  %.3f %.3f %.3f   -> x%.2f x%.2f\n",
              arm, v[1], v[2], v[3], v[2] / v[1], v[3] / v[1]))
  s2 <- sc[sc$arm == arm & sc$est == "ols" & sc$fit == "LS" & sc$param == "beta_sigma2", ]
  v2 <- vapply(c(30000, 1e5, 461320), function(N) mean(abs(s2[s2$N == N, "bias_sd"])), 0)
  cat(sprintf("arm %s  OLS LS:beta_sigma2  %.3f %.3f %.3f   -> x%.2f x%.2f\n",
              arm, v2[1], v2[2], v2[3], v2[2] / v2[1], v2[3] / v2[1]))
}

cat("\n=== N_eff = 2 N p(1-p), which is what actually moves ===\n")
for (arm in c("A", "B", "C")) {
  s <- tab[tab$arm == arm, ]
  cat(sprintf("arm %s  MAC:", arm))
  for (N in c(30000, 1e5, 461320)) cat(sprintf(" %9.0f", mean(s[s$N == N, "mac"])))
  cat("\n")
}
