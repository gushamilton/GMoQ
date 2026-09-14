d <- "/user/work/fh6520/nscale_run/smoke/out"
t <- read.delim(file.path(d, "summary.tsv"))
cat("--- df law columns (unique) ---\n")
print(unique(t[, c("cell", "param", "nu_jk", "zcrit_t", "t1_5e8_ratio", "t1_5e8_ratio_t")]))
cat("\n--- any NA in key cols? ---\n")
print(colSums(is.na(t[, c("nu_jk", "zcrit_t", "eff_sd_gls_over_ols", "cov95_t")])))
p <- read.delim(file.path(d, "pertau.tsv"))
cat("\n--- per-tau, extreme taus ---\n")
print(p[p$tau < 0.07 | p$tau > 0.93,
        c("cell", "tau", "n_below", "cov95", "se_ratio_tau", "fhat_bias", "sd_z")])
