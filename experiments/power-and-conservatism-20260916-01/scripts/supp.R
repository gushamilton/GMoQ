#!/usr/bin/env Rscript
# Task 1 supplements: inflation on the real scans, what the G-MoQ-only hits look like
# in REGENIE, and how much of the recall is predicted by the SE ratio alone.
suppressPackageStartupMessages(library(data.table))
OUT <- "/user/work/fh6520/gmoq_power_conserv/out"
set.seed(2)
PAR <- list(beta_left = list(g="L", r="f"), beta_right = list(g="R", r="m"), gamma = list(g="G", r="s"))

allD <- list(); for (ch in c(22, 21)) allD[[as.character(ch)]] <- readRDS(file.path(OUT, sprintf("merged_chr%d.rds", ch)))

## 1. genomic inflation on the real scans (dominated by true signal, so read as a ceiling)
inf <- rbindlist(lapply(names(allD), function(ch) {
  D <- allD[[ch]]
  rbindlist(lapply(names(PAR), function(pn) {
    g <- PAR[[pn]]$g; r <- PAR[[pn]]$r
    lam <- function(p) median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
    data.table(chr = as.integer(ch), param = pn,
               lambda_gmoq_subtract  = lam(D[[paste0("p", g)]]),
               lambda_gmoq_covariate = lam(D[[paste0("p", g, ".cov")]]),
               lambda_regenie        = lam(D[[paste0("rp_", r)]]))
  }))
}))
cat("=== genomic control on the REAL scans (signal included; not a null) ===\n"); print(inf)

## 2. where do the G-MoQ-only hits sit in REGENIE?
own <- rbindlist(lapply(names(allD), function(ch) {
  D <- allD[[ch]]
  rbindlist(lapply(c("beta_left", "beta_right"), function(pn) {
    g <- PAR[[pn]]$g; r <- PAR[[pn]]$r
    rbindlist(lapply(c(5e-8, 1e-5), function(thr) {
      only <- D[[paste0("p", g, ".cov")]] < thr & !(D[[paste0("rp_", r)]] < thr)
      rp <- D[[paste0("rp_", r)]][only]
      data.table(chr = as.integer(ch), param = pn, thr = thr, n_gmoq_only = sum(only),
                 regenie_p_median = if (sum(only)) median(rp) else NA_real_,
                 frac_regenie_p_lt_1e3 = if (sum(only)) mean(rp < 1e-3) else NA_real_,
                 frac_regenie_p_lt_0.05 = if (sum(only)) mean(rp < 0.05) else NA_real_,
                 expected_frac_lt_0.05_if_null = 0.05)
    }))
  }))
}))
cat("\n=== G-MoQ-only hits: their REGENIE p-values ===\n"); print(own)

## 3. recall predicted by the SE ratio alone, against observed
pred <- rbindlist(lapply(names(allD), function(ch) {
  D <- allD[[ch]]
  rbindlist(lapply(names(PAR), function(pn) {
    g <- PAR[[pn]]$g; r <- PAR[[pn]]$r
    rbindlist(lapply(c(5e-8, 1e-5), function(thr) {
      zc <- qnorm(thr/2, lower.tail = FALSE)
      rhit <- D[[paste0("rp_", r)]] < thr
      if (!sum(rhit)) return(NULL)
      zR <- qnorm(D[[paste0("rp_", r)]][rhit]/2, lower.tail = FALSE)
      out <- data.table(chr = as.integer(ch), param = pn, thr = thr, n_regenie = sum(rhit))
      for (m in c("", ".cov")) {
        k <- median(D[[paste0("s", g, m)]] / D[[paste0("rs_", r)]])
        nm <- if (m == "") "subtract" else "covariate"
        out[[paste0("se_ratio_", nm)]] <- k
        out[[paste0("recall_pred_", nm)]] <- mean(zR / k > zc)
        out[[paste0("recall_obs_", nm)]]  <- mean(D[[paste0("p", g, m)]][rhit] < thr)
      }
      out
    }))
  }))
}))
cat("\n=== recall: predicted from the SE ratio alone vs observed ===\n"); print(pred)

fwrite(inf, file.path(OUT, "inflation_real.tsv"), sep = "\t")
fwrite(own, file.path(OUT, "gmoq_only_hits.tsv"), sep = "\t")
fwrite(pred, file.path(OUT, "recall_pred.tsv"), sep = "\t")

## 4. pooled chr21+chr22 discovery table with a paired bootstrap CI on the gain
d22 <- fread(file.path(OUT, "discovery_chr22.tsv")); d21 <- fread(file.path(OUT, "discovery_chr21.tsv"))
D <- rbind(d22, d21)
pooled <- D[, .(chr = "21+22", gmoq_hits = sum(gmoq_hits), regenie_hits = sum(regenie_hits),
                overlap = sum(overlap)), by = .(param, thr, loco_mode)]
pooled[, `:=`(precision = overlap/gmoq_hits, recall = overlap/regenie_hits)]
cat("\n=== pooled chr21+chr22 ===\n"); print(pooled)
fwrite(pooled, file.path(OUT, "discovery_pooled.tsv"), sep = "\t")

## 5. block-bootstrap CI on the discovery gain (blocks of 250 contiguous SNPs)
bs <- rbindlist(lapply(names(allD), function(ch) {
  Dt <- copy(allD[[ch]]); setorder(Dt, bp); P <- nrow(Dt)
  Dt[, blk := (seq_len(P) - 1L) %/% 250L]
  nb <- max(Dt$blk) + 1L
  rbindlist(lapply(names(PAR), function(pn) {
    g <- PAR[[pn]]$g
    rbindlist(lapply(c(5e-8, 1e-5), function(thr) {
      a <- Dt[[paste0("p", g)]] < thr; b <- Dt[[paste0("p", g, ".cov")]] < thr
      bb <- replicate(2000, { s <- sample(0:(nb-1), nb, TRUE)
        i <- unlist(lapply(s, function(k) which(Dt$blk == k)))
        sum(b[i]) - sum(a[i]) })
      data.table(chr = as.integer(ch), param = pn, thr = thr,
                 gain = sum(b) - sum(a), gain_lo = quantile(bb, .025), gain_hi = quantile(bb, .975))
    }))
  }))
}))
cat("\n=== block-bootstrap 95% CI on the discovery gain (covariate - subtract) ===\n"); print(bs)
fwrite(bs, file.path(OUT, "gain_bootstrap.tsv"), sep = "\t")
