## Follow-up diagnostics: (i) how close is (d) to (b)?  (ii) is the residual
## lambda signal or inflation?  SNP-level aggregates only.
suppressPackageStartupMessages(library(data.table)); setDTthreads(4)
options(width = 210)

runs <- list(
  a  = "/user/work/fh6520/Gmoq_paper/.local/runs/paper-stage2-testosterone-grid-all-20260916-01/testosterone/chr22/pooled_ols/chr22.tsv.gz",
  b  = "/user/work/fh6520/gmoq_noloco/stage2/pooled_ols/chr22.tsv.gz",
  c  = "/user/work/fh6520/gmoq_twopass/stage2/pooled_ols_resid_mean/chr22.tsv.gz",
  d  = "/user/work/fh6520/gmoq_locovar/stage2/covariate/chr22.tsv.gz"
)
dat <- lapply(runs, function(p) fread(cmd = paste("gzip -dc", shQuote(p)), showProgress = FALSE))
pars <- c("beta_left", "beta_right", "gamma")

cat("=== between-run agreement of point estimates (Pearson r on shared SNPs) ===\n")
cmp <- rbindlist(lapply(pars, function(par) {
  rbindlist(lapply(list(c("a","b"), c("a","d"), c("b","d"), c("b","c"), c("c","d")), function(pp) {
    m <- merge(dat[[pp[1]]][, .(snp_id, x = get(par), xs = get(paste0(par,"_se")))],
               dat[[pp[2]]][, .(snp_id, y = get(par), ys = get(paste0(par,"_se")))], by = "snp_id")
    data.table(param = par, pair = paste(pp, collapse = " vs "),
               r_beta = round(cor(m$x, m$y), 5),
               r_z = round(cor(m$x/m$xs, m$y/m$ys), 5),
               med_se_ratio = round(median(m$xs/m$ys), 4))
  }))
}))
print(cmp)

## Signal vs inflation.  The chr22 testosterone signal is dominated by SHBG.
## Use the sex-stratified REGENIE scans (themselves run with --ignore-pred) to
## define signal windows, drop everything within 1 Mb of a |z| >= 5 variant, and
## recompute lambda on what is left.
reg_dir <- "/user/work/fh6520/fungwas_paper/results/regenie/combined"
rr <- rbindlist(lapply(c("testosterone_female_testosterone_all.regenie",
                         "testosterone_male_testosterone_all.regenie",
                         "sex_all.regenie"), function(f) {
  d <- fread(file.path(reg_dir, f), showProgress = FALSE)[CHROM == 22]
  d[, .(ID, GENPOS, z = BETA/SE)]
}))
hits <- rr[abs(z) >= 5, unique(GENPOS)]
cat("\nREGENIE |z|>=5 positions on chr22:", length(hits), "\n")

lam <- function(p) { p <- p[is.finite(p) & p > 0 & p <= 1]
  median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }

near <- function(bp, hits, win = 1e6) {
  if (!length(hits)) return(rep(FALSE, length(bp)))
  vapply(bp, function(x) any(abs(hits - x) <= win), TRUE)
}
out <- rbindlist(lapply(names(runs), function(k) {
  s <- dat[[k]]
  drop <- near(s$bp, hits)
  rbindlist(lapply(pars, function(par) {
    data.table(run = k, param = par, n_all = nrow(s), n_kept = sum(!drop),
               lambda_all = round(lam(s[[paste0(par,"_p")]]), 4),
               lambda_nosignal = round(lam(s[[paste0(par,"_p")]][!drop]), 4))
  }))
}))
cat("\n=== lambda before and after masking 1 Mb around REGENIE |z|>=5 ===\n")
print(dcast(out, param + n_all + n_kept ~ run, value.var = c("lambda_all","lambda_nosignal")))

perm <- "/user/work/fh6520/gmoq_perm/stage2/pooled_ols/chr22.tsv.gz"
if (file.exists(perm)) {
  s <- fread(cmd = paste("gzip -dc", shQuote(perm)), showProgress = FALSE)
  drop <- near(s$bp, hits)
  cat("\npermuted null, all SNPs   :",
      paste(sprintf("%s=%.4f", pars, vapply(pars, function(p) lam(s[[paste0(p,"_p")]]), 0)), collapse="  "))
  cat("\npermuted null, masked     :",
      paste(sprintf("%s=%.4f", pars, vapply(pars, function(p) lam(s[[paste0(p,"_p")]][!drop]), 0)), collapse="  "), "\n")
}

## Where did (d) move gamma relative to (b)?  If this is deconfounding the shift
## should be a shrinkage of the estimates, not a rescaling of the SEs.
m <- merge(dat$b[, .(snp_id, bb = gamma, bs = gamma_se)],
           dat$d[, .(snp_id, db = gamma, ds = gamma_se)], by = "snp_id")
cat("\n=== gamma: (d) against (b) ===\n")
cat(sprintf("median |beta| ratio d/b : %.4f\n", median(abs(m$db)/abs(m$bb))))
cat(sprintf("regression slope d ~ b  : %.4f\n", coef(lm(db ~ bb, data = m))[[2]]))
cat(sprintf("median SE ratio  d/b    : %.4f\n", median(m$ds/m$bs)))
cat(sprintf("median chi2 d / median chi2 b : %.4f\n",
            median((m$db/m$ds)^2)/median((m$bb/m$bs)^2)))
fwrite(out, "/user/work/fh6520/gmoq_locovar/analysis/lambda_masked.tsv", sep = "\t")
fwrite(cmp, "/user/work/fh6520/gmoq_locovar/analysis/between_run.tsv", sep = "\t")
