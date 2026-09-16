#!/usr/bin/env Rscript
# Task 1: did the LOCO change buy discoveries?
suppressPackageStartupMessages({library(data.table)})
args <- commandArgs(TRUE)
CHR   <- as.integer(args[[1]])
SUB   <- args[[2]]   # subtract-mode stage2 tsv.gz
COV   <- args[[3]]   # covariate-mode stage2 tsv.gz
RGF   <- args[[4]]   # regenie female
RGM   <- args[[5]]   # regenie male
RGS   <- args[[6]]   # regenie sex
OUT   <- args[[7]]

rd_g <- function(f) {
  d <- fread(f)
  d[, .(snp_id, chr, bp, effect_allele, other_allele,
        gamma, gamma_se, gamma_p, gamma_df,
        beta_left, beta_left_se, beta_left_p,
        beta_right, beta_right_se, beta_right_p)]
}
rd_r <- function(f, tag) {
  d <- fread(f)
  d <- d[TEST == "ADD"]
  d[, p := 10^(-LOG10P)]
  keep <- c("ID","ALLELE1","ALLELE0","BETA","SE","N","A1FREQ","INFO","p")
  d <- d[, keep, with = FALSE]
  setnames(d, keep, paste0(c("snp_id","a1","a0","beta","se","n","af","info","p"), "_", tag))
  d
}

S <- rd_g(SUB); C <- rd_g(COV)
F_ <- rd_r(RGF, "f"); M_ <- rd_r(RGM, "m"); X_ <- rd_r(RGS, "s")
setnames(F_, "snp_id_f", "snp_id"); setnames(M_, "snp_id_m", "snp_id"); setnames(X_, "snp_id_s", "snp_id")

D <- merge(merge(merge(S, F_, by = "snp_id"), M_, by = "snp_id"), X_, by = "snp_id")
setnames(C, setdiff(names(C), "snp_id"), paste0(setdiff(names(C), "snp_id"), ".cov"))
D <- merge(D, C, by = "snp_id")

## --- allele harmonisation against REGENIE ALLELE1 ------------------------
D[, align := fifelse(effect_allele == a1_f & other_allele == a0_f,  1,
             fifelse(effect_allele == a0_f & other_allele == a1_f, -1, NA_real_))]
cat("harmonisation: aligned", sum(D$align == 1, na.rm=TRUE),
    " flipped", sum(D$align == -1, na.rm=TRUE),
    " unmatched", sum(is.na(D$align)), "\n")
stopifnot(all(D$a1_f == D$a1_m), all(D$a1_f == D$a1_s))
D <- D[!is.na(align)]

for (mode in c("", ".cov")) {
  D[, (paste0("bL", mode)) := get(paste0("beta_left",  mode)) * align]
  D[, (paste0("bR", mode)) := get(paste0("beta_right", mode)) * align]
  D[, (paste0("bG", mode)) := -get(paste0("gamma",     mode)) * align]  # gamma negated: REGENIE Male=1
  D[, (paste0("sL", mode)) := get(paste0("beta_left_se",  mode))]
  D[, (paste0("sR", mode)) := get(paste0("beta_right_se", mode))]
  D[, (paste0("sG", mode)) := get(paste0("gamma_se",      mode))]
  D[, (paste0("pL", mode)) := get(paste0("beta_left_p",  mode))]
  D[, (paste0("pR", mode)) := get(paste0("beta_right_p", mode))]
  D[, (paste0("pG", mode)) := get(paste0("gamma_p",      mode))]
}
D[, `:=`(rb_f = beta_f, rs_f = se_f, rp_f = p_f,
         rb_m = beta_m, rs_m = se_m, rp_m = p_m,
         rb_s = beta_s, rs_s = se_s, rp_s = p_s)]

PARAMS <- list(beta_left = list(g = "L", r = "f"),
               beta_right = list(g = "R", r = "m"),
               gamma = list(g = "G", r = "s"))
THR <- c(`5e-8` = 5e-8, `1e-5` = 1e-5)

se_tab <- rbindlist(lapply(names(PARAMS), function(pn) {
  g <- PARAMS[[pn]]$g; r <- PARAMS[[pn]]$r
  data.table(chr = CHR, param = pn,
             se_ratio_subtract  = median(D[[paste0("s", g)]]        / D[[paste0("rs_", r)]]),
             se_ratio_covariate = median(D[[paste0("s", g, ".cov")]] / D[[paste0("rs_", r)]]),
             beta_shift_med_in_SE = median(abs(D[[paste0("b", g, ".cov")]] - D[[paste0("b", g)]]) /
                                           D[[paste0("s", g)]]),
             beta_shift_max_in_SE = max(abs(D[[paste0("b", g, ".cov")]] - D[[paste0("b", g)]]) /
                                        D[[paste0("s", g)]]))
}))

## --- crude independent-locus count: greedy 500kb clumping ---------------
n_loci <- function(bp, p, thr, win = 5e5) {
  idx <- which(p < thr); if (!length(idx)) return(0L)
  o <- idx[order(p[idx])]; keep <- integer(0)
  for (i in o) if (!any(abs(bp[i] - bp[keep]) < win)) keep <- c(keep, i)
  length(keep)
}

rows <- list(); chg <- list()
for (pn in names(PARAMS)) {
  g <- PARAMS[[pn]]$g; r <- PARAMS[[pn]]$r
  rp <- D[[paste0("rp_", r)]]
  for (tn in names(THR)) {
    thr <- THR[[tn]]
    rhit <- rp < thr
    for (mode in c("subtract", "covariate")) {
      sfx <- if (mode == "subtract") "" else ".cov"
      gp <- D[[paste0("p", g, sfx)]]
      ghit <- gp < thr
      ov <- sum(ghit & rhit)
      rows[[length(rows)+1]] <- data.table(
        chr = CHR, param = pn, thr = tn, loco_mode = mode,
        gmoq_hits = sum(ghit), regenie_hits = sum(rhit), overlap = ov,
        precision = ifelse(sum(ghit) > 0, ov/sum(ghit), NA_real_),
        recall    = ifelse(sum(rhit) > 0, ov/sum(rhit), NA_real_),
        gmoq_loci = n_loci(D$bp, gp, thr), regenie_loci = n_loci(D$bp, rp, thr),
        sign_conc_at_regenie_hit = if (sum(rhit) > 0)
          mean(sign(D[[paste0("b", g, sfx)]][rhit]) == sign(D[[paste0("rb_", r)]][rhit])) else NA_real_)
    }
    zcf <- D[[paste0("b", g)]] / D[[paste0("s", g, ".cov")]]
    pcf <- 2 * pt(-abs(zcf), df = D$gamma_df)
    zbb <- D[[paste0("b", g, ".cov")]] / D[[paste0("s", g)]]
    pbb <- 2 * pt(-abs(zbb), df = D$gamma_df)
    a <- D[[paste0("p", g)]] < thr; b <- D[[paste0("p", g, ".cov")]] < thr
    chg[[length(chg)+1]] <- data.table(
      chr = CHR, param = pn, thr = tn,
      n_sub = sum(a), n_cov = sum(b),
      gained = sum(!a & b), lost = sum(a & !b), changed = sum(a != b),
      mcnemar_p = if (sum(a != b) > 0)
        stats::binom.test(sum(!a & b), sum(a != b), 0.5)$p.value else NA_real_,
      n_se_only   = sum(pcf < thr),   # subtract beta + covariate SE
      n_beta_only = sum(pbb < thr),   # covariate beta + subtract SE
      gained_true_pos = sum(!a & b & rhit), lost_true_pos = sum(a & !b & rhit))
  }
}
res <- rbindlist(rows); chgt <- rbindlist(chg)
fwrite(res,  file.path(OUT, sprintf("discovery_chr%d.tsv", CHR)), sep = "\t")
fwrite(chgt, file.path(OUT, sprintf("statuschange_chr%d.tsv", CHR)), sep = "\t")
fwrite(se_tab, file.path(OUT, sprintf("seratio_chr%d.tsv", CHR)), sep = "\t")
cat("\n=== n variants analysed:", nrow(D), "===\n")
print(se_tab); cat("\n"); print(res); cat("\n"); print(chgt)
saveRDS(D, file.path(OUT, sprintf("merged_chr%d.rds", CHR)))
