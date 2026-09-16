## WHERE DOES DIVERGENCE BEGIN?
## For the production 3-column W at testosterone separation, find the |Z| at
## which the pooled arm's reported SE first falls 1% / 5% / 10% below the truth,
## per parameter and MAF.  Everything exact.
WS <- "/user/work/fh6520/gmoq_covmixW"
.libPaths(c(file.path(WS, "lib")))
stopifnot(identical(.libPaths()[1], file.path(WS, "lib")))
library(fungwasStage2)
source(file.path(WS, "src", "truth_mix.R"))
TAUS <- seq(0.02, 0.98, length.out = 45); N <- 461320L
SH <- "mix5.89"

ex_at <- function(maf, cc, wn) {
  tr <- pop_truth(TAUS, maf, SH, 0, cc)
  ex <- exact_W(build_Ws(TAUS, tr)[[wn]], tr, N * tr$Vg)
  list(Z = max(abs(ex$Z_true)), r = ex$ratio, cov = ex$cov95, infl = ex$infl_5e8)
}
## invert: smallest c whose worst-parameter se_ratio drops below (1 - thr)
find_c <- function(maf, wn, thr, parm = NULL) {
  f <- function(cc) {
    e <- ex_at(maf, cc, wn)
    rr <- if (is.null(parm)) min(e$r) else e$r[parm]
    (1 - thr) - rr
  }
  lo <- 1e-6; hi <- 1e-3
  while (f(hi) < 0) { hi <- hi * 2; if (hi > 64) return(NA_real_) }
  uniroot(f, c(lo, hi), tol = 1e-9)$root
}

cat("\n########## O1. onset of pooled mis-calibration, production W (mix3, K=3) ##########\n")
cat("shape = mix5.89 (testosterone-like), TRUE DISPERSION effect.\n")
cat("Entry is the |Z| the SHIPPED 3-column W would report on that variant at the\n")
cat("effect size where the stated SE shortfall is first reached.  '-' = never.\n\n")
cat(sprintf("%6s %-12s | %14s %14s %14s\n", "MAF", "param",
            "|Z| @ 1% short", "|Z| @ 5% short", "|Z| @ 10% short"))
for (maf in c(0.005, 0.01, 0.05, 0.30)) {
  for (parm in c("gamma", "beta_1")) {
    out <- character(0)
    for (thr in c(0.01, 0.05, 0.10)) {
      cc <- find_c(maf, "mix3", thr, parm)
      out <- c(out, if (is.finite(cc)) sprintf("%14.1f", ex_at(maf, cc, "mix3")$Z) else
                    sprintf("%14s", "-"))
    }
    cat(sprintf("%6g %-12s | %s\n", maf, parm, paste(out, collapse = "")))
  }
  ## the published 2-column vqtl W on the SAME phenotype, for contrast
  out <- character(0)
  for (thr in c(0.01, 0.05, 0.10)) {
    cc <- find_c(maf, "vqtl2", thr, "beta_mu")
    out <- c(out, if (is.finite(cc)) sprintf("%14.1f", ex_at(cc, maf, "vqtl2")$Z) else
                  sprintf("%14s", "-"))
  }
}

cat("\n########## O2. the chr22-testosterone regime ##########\n")
cat("Real chr22 testosterone: N=399,710, max |z| ~ 11, common variants.\n")
cat("Ratio of the two arms' reported SEs there, exactly:\n\n")
cat(sprintf("%6s %5s | %-12s %10s %10s %10s\n", "MAF", "|Z|", "param",
            "se_ratio", "cov95", "infl@5e-8"))
for (maf in c(0.05, 0.30)) for (Zt in c(5, 8, 11)) {
  f <- function(cc) ex_at(maf, cc, "mix3")$Z - Zt
  hi <- 1e-3; while (f(hi) < 0) hi <- hi * 2
  cc <- uniroot(f, c(1e-9, hi), tol = 1e-10)$root
  e <- ex_at(maf, cc, "mix3")
  for (j in seq_along(e$r))
    cat(sprintf("%6g %5g | %-12s %10.5f %10.5f %10.4f\n",
                maf, Zt, names(e$r)[j], e$r[j], e$cov[j], e$infl[j]))
}
cat("\nDONE\n")
