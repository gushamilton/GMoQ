## Four-way (plus a same-code control) comparison of LOCO handling.
## testosterone, chr22, skew family, pooled covariance + OLS.
## Writes SNP-level aggregate summaries only.
suppressPackageStartupMessages(library(data.table)); setDTthreads(4)

reg_dir <- "/user/work/fh6520/fungwas_paper/results/regenie/combined"
read_reg <- function(f) {
  d <- fread(file.path(reg_dir, f), showProgress = FALSE)[CHROM == 22]
  d[, .(ID, A0 = ALLELE0, A1 = ALLELE1, BETA, SE)]
}
reg <- list(
  beta_left  = read_reg("testosterone_female_testosterone_all.regenie"),
  beta_right = read_reg("testosterone_male_testosterone_all.regenie"),
  gamma      = read_reg("sex_all.regenie")
)
## REGENIE codes sex Male = 1; the mixture's gamma is signed the other way.
sign_flip <- c(beta_left = 1, beta_right = 1, gamma = -1)
reg_label <- c(beta_left = "female", beta_right = "male", gamma = "sex")

lam <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
}
lam_z <- function(z) {
  z <- z[is.finite(z)]
  median(z^2) / qchisq(0.5, 1)
}

harmonise <- function(s, par) {
  r <- reg[[par]]
  m <- merge(s[, .(snp_id, ea = effect_allele, oa = other_allele,
                   b = get(par), se = get(paste0(par, "_se")),
                   p = get(paste0(par, "_p")))],
             r, by.x = "snp_id", by.y = "ID")
  same <- m$ea == m$A1 & m$oa == m$A0
  swap <- m$ea == m$A0 & m$oa == m$A1
  m <- m[same | swap]
  sgn <- ifelse(m$ea == m$A1 & m$oa == m$A0, 1, -1) * sign_flip[[par]]
  m[, bh := b * sgn]
  m[is.finite(bh) & is.finite(se) & se > 0 & is.finite(BETA) & is.finite(SE) & SE > 0]
}

runs <- list(
  "(a)  subtract-from-Y, production"   =
    "/user/work/fh6520/Gmoq_paper/.local/runs/paper-stage2-testosterone-grid-all-20260916-01/testosterone/chr22/pooled_ols/chr22.tsv.gz",
  "(a') subtract-from-Y, same code"    =
    "/user/work/fh6520/gmoq_locovar/stage2/subtract/chr22.tsv.gz",
  "(b)  no LOCO at all"                =
    "/user/work/fh6520/gmoq_noloco/stage2/pooled_ols/chr22.tsv.gz",
  "(c)  two-pass blind LOCO"           =
    "/user/work/fh6520/gmoq_twopass/stage2/pooled_ols_resid_mean/chr22.tsv.gz",
  "(d)  LOCO as RIF-column covariate"  =
    "/user/work/fh6520/gmoq_locovar/stage2/covariate/chr22.tsv.gz"
)
have <- vapply(runs, file.exists, TRUE)
if (any(!have)) message("MISSING: ", paste(names(runs)[!have], collapse = "; "))
runs <- runs[have]

dat <- lapply(runs, function(p) fread(cmd = paste("gzip -dc", shQuote(p)), showProgress = FALSE))

res <- rbindlist(lapply(names(runs), function(k) {
  s <- dat[[k]]
  rbindlist(lapply(names(reg), function(par) {
    m <- harmonise(s, par)
    sub <- m[abs(BETA / SE) >= 4]
    data.table(
      run = k, param = par, regenie = reg_label[[par]], n = nrow(m),
      se_ratio  = median(m$se / m$SE),
      absb_ratio = median(abs(m$bh) / abs(m$BETA)),
      cor_all   = cor(m$bh, m$BETA),
      n_hiz     = nrow(sub),
      cor_hiz   = if (nrow(sub) > 2) cor(sub$bh, sub$BETA) else NA_real_,
      slope_hiz = if (nrow(sub) > 2) coef(lm(bh ~ BETA, data = sub))[[2]] else NA_real_,
      lambda    = lam(s[[paste0(par, "_p")]])
    )
  }))
}))

## Is the subtract-LOCO lambda advantage structure control, or just bigger SEs?
## Take each run's own point estimates and give them the SE of run (a'), then
## recompute lambda.  If a run's lambda lands on (a')'s under that swap, the
## whole difference was SE inflation.
ref_key <- if ("(a') subtract-from-Y, same code" %in% names(runs)) {
  "(a') subtract-from-Y, same code"
} else names(runs)[[1]]
swap <- rbindlist(lapply(names(runs), function(k) {
  rbindlist(lapply(names(reg), function(par) {
    a <- dat[[ref_key]][, .(snp_id, se_ref = get(paste0(par, "_se")))]
    m <- merge(dat[[k]][, .(snp_id, b = get(par), se = get(paste0(par, "_se")))],
               a, by = "snp_id")
    m <- m[is.finite(b) & is.finite(se) & se > 0 & is.finite(se_ref) & se_ref > 0]
    data.table(run = k, param = par, n = nrow(m),
               lambda_own      = lam_z(m$b / m$se),
               lambda_ref_se   = lam_z(m$b / m$se_ref),
               lambda_own_beta_ref_se_ratio = lam_z(m$b / m$se_ref) / lam_z(m$b / m$se),
               med_se_vs_ref   = median(m$se / m$se_ref))
  }))
}))

options(width = 210)
for (par in names(reg)) {
  cat("\n=========  ", par, "   (REGENIE reference: ", reg_label[[par]], ")\n", sep = "")
  print(res[param == par, .(run, n,
                            se_ratio = round(se_ratio, 4),
                            absb_ratio = round(absb_ratio, 3),
                            cor_all = round(cor_all, 4),
                            n_hiz, cor_hiz = round(cor_hiz, 4),
                            slope_hiz = round(slope_hiz, 3),
                            lambda = round(lambda, 4))])
}
cat("\n\n=========  SE-swap diagnostic (reference SEs = ", ref_key, ")\n", sep = "")
print(swap[, .(run, param, med_se_vs_ref = round(med_se_vs_ref, 4),
               lambda_own = round(lambda_own, 4),
               lambda_with_ref_se = round(lambda_ref_se, 4))])

perm <- "/user/work/fh6520/gmoq_perm/stage2/pooled_ols/chr22.tsv.gz"
if (file.exists(perm)) {
  s <- fread(cmd = paste("gzip -dc", shQuote(perm)), showProgress = FALSE)
  cat("\npermuted true null lambda:",
      paste(sprintf("%s=%.4f", names(reg),
                    vapply(names(reg), function(p) lam(s[[paste0(p, "_p")]]), 0)),
            collapse = "  "), "\n")
}

fwrite(res,  "/user/work/fh6520/gmoq_locovar/analysis/fourway.tsv", sep = "\t")
fwrite(swap, "/user/work/fh6520/gmoq_locovar/analysis/se_swap.tsv", sep = "\t")
cat("\nwrote /user/work/fh6520/gmoq_locovar/analysis/{fourway,se_swap}.tsv\n")
