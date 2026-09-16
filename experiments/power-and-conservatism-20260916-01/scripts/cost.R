#!/usr/bin/env Rscript
# What a conservatism of factor k would actually cost, on the real chr21+chr22 scans.
suppressPackageStartupMessages(library(data.table))
OUT <- "/user/work/fh6520/gmoq_power_conserv/out"
PAR <- list(beta_left = list(g="L", r="f"), beta_right = list(g="R", r="m"), gamma = list(g="G", r="s"))
D <- rbindlist(lapply(c(22, 21), function(ch) {
  d <- readRDS(file.path(OUT, sprintf("merged_chr%d.rds", ch)))
  d[, chrn := ch][]
}), fill = TRUE)

res <- rbindlist(lapply(names(PAR), function(pn) {
  g <- PAR[[pn]]$g; r <- PAR[[pn]]$r
  z <- abs(qnorm(D[[paste0("p", g, ".cov")]] / 2))
  rbindlist(lapply(c(1, 1.5, 2.5, 5), function(k) {
    zc <- qnorm(5e-8 / (2 * k), lower.tail = FALSE)
    data.table(param = pn, k = k, crit = zc,
               hits_5e8 = sum(z > qnorm(2.5e-8, lower.tail = FALSE)),
               hits_under_k = sum(z > zc),
               lost = sum(z > qnorm(2.5e-8, lower.tail = FALSE)) - sum(z > zc),
               regenie_hits = sum(D[[paste0("rp_", r)]] < 5e-8),
               recall_nominal = mean(D[[paste0("p", g, ".cov")]][D[[paste0("rp_", r)]] < 5e-8] < 5e-8),
               recall_under_k = mean(z[D[[paste0("rp_", r)]] < 5e-8] > zc))
  }))
}))
print(res)
fwrite(res, file.path(OUT, "cost_of_conservatism.tsv"), sep = "\t")
