#!/usr/bin/env Rscript
# What the permuted-phenotype null can and cannot establish, with LD-aware uncertainty.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); PERM <- a[[1]]; OUT <- a[[2]]
set.seed(1)
d <- fread(PERM)
setorder(d, bp)
P <- nrow(d)
cat("permuted-null variants:", P, "\n")
pars <- c("gamma", "beta_left", "beta_right")
NB <- 100L                      # contiguous LD blocks for the block bootstrap
blk <- cut(seq_len(P), breaks = NB, labels = FALSE)
NBOOT <- 2000L
alphas <- c(0.05, 0.01, 1e-3, 1e-4, 1e-5)

boot_idx <- lapply(seq_len(NBOOT), function(b) {
  bs <- sample(seq_len(NB), NB, replace = TRUE)
  unlist(lapply(bs, function(k) which(blk == k)), use.names = FALSE)
})

res <- list()
for (pn in pars) {
  pv <- d[[paste0(pn, "_p")]]
  z2 <- qchisq(pv, 1, lower.tail = FALSE)
  lam <- median(z2) / qchisq(0.5, 1)
  # scale factor c such that z_obs = c * z_null, from the median
  cmed <- sqrt(lam)
  lam_b <- vapply(boot_idx, function(i) median(z2[i]) / qchisq(0.5, 1), 0)
  for (al in alphas) {
    obs <- sum(pv < al)
    ob <- vapply(boot_idx, function(i) sum(pv[i] < al), 0)
    res[[length(res) + 1]] <- data.table(
      param = pn, alpha = al, expected = al * P, observed = obs,
      ratio = obs / (al * P),
      se_indep = sqrt(al * (1 - al) * P) / (al * P),
      se_block_boot = sd(ob) / (al * P),
      lambda = lam, lambda_se_block = sd(lam_b),
      c_scale_from_median = cmed)
  }
}
R <- rbindlist(res)
# how many independent tests the block bootstrap implies
R[, m_eff_implied := (se_indep / se_block_boot)^2 * P]
fwrite(R, file.path(OUT, "perm_null_rates.tsv"), sep = "\t")
print(R)

# Is the observed tail consistent with a pure scale deflation fitted at the median?
cat("\n--- pure-scale test: predicted counts from lambda alone ---\n")
for (pn in pars) {
  pv <- d[[paste0(pn, "_p")]]; z2 <- qchisq(pv, 1, lower.tail = FALSE)
  cs <- sqrt(median(z2)/qchisq(0.5, 1))
  pred <- vapply(alphas, function(al) {
    zc <- qnorm(al/2, lower.tail = FALSE)
    P * 2 * pnorm(zc / cs, lower.tail = FALSE)
  }, 0)
  cat(sprintf("%-11s  ", pn))
  for (j in seq_along(alphas))
    cat(sprintf("a=%-6g obs=%-5d pred=%.1f | ", alphas[j], sum(pv < alphas[j]), pred[j]))
  cat("\n")
}
cat("\nc_scale extrapolated rejection rate at 5e-8 (pure-scale assumption):\n")
for (pn in pars) {
  pv <- d[[paste0(pn, "_p")]]; cs <- sqrt(median(qchisq(pv,1,lower.tail=FALSE))/qchisq(0.5,1))
  r <- 2*pnorm(qnorm(2.5e-8, lower.tail = FALSE)/cs, lower.tail = FALSE)
  cat(sprintf("  %-11s c=%.4f  rate=%.3e  (nominal 5e-8, factor %.2f)\n", pn, cs, r, r/5e-8))
}
