WS <- "/user/work/fh6520/gmoq_covmixW"
tab <- read.delim(file.path(WS, "out", "summary_mix.tsv"), stringsAsFactors = FALSE)
tab <- tab[tab$est == "ols", ]
cat("cells present:", length(unique(tab$cell)), "  reps/cell:", unique(tab$nrep), "\n\n")

show <- function(sub, title) {
  if (!nrow(sub)) return(invisible())
  cat("\n##", title, "\n")
  cat(sprintf("%-28s %-6s %-2s %-12s %-10s | %8s %7s | %8s %7s | %10s %10s | %9s\n",
              "cell", "W", "K", "param", "cov_source", "cov95", "+-mcse",
              "se_rat", "+-mcse", "t1@0.05", "+-mcse", "infl@5e-8"))
  sub <- sub[order(sub$cell, sub$fit, sub$param, sub$cov_source), ]
  for (i in seq_len(nrow(sub))) with(sub[i, ], cat(sprintf(
    "%-28s %-6s %-2d %-12s %-10s | %8.4f %7.4f | %8.4f %7.4f | %10.4f %10.4f | %9.3g\n",
    cell, fit, K, param, cov_source, cov95_t, cov95_t_mcse,
    se_ratio, se_ratio_mcse, t1_05, t1_05_mcse, infl_5e8_from_sd)))
}
for (cl in unique(tab$cell)) show(tab[tab$cell == cl, ], cl)

pr <- file.path(WS, "out", "paired_mix.tsv")
if (file.exists(pr)) {
  p <- read.delim(pr, stringsAsFactors = FALSE)
  p <- p[p$est == "ols", ]
  cat("\n## PAIRED jackknife - pooled (same replicates)\n")
  cat(sprintf("%-28s %-6s %-12s | %10s %10s | %10s %10s | %12s\n",
              "cell", "W", "param", "d_cov95", "d_se_ratio", "d_t1@.05",
              "+-paired", "sd_jk/sd_pool"))
  for (i in seq_len(nrow(p))) with(p[i, ], cat(sprintf(
    "%-28s %-6s %-12s | %10.4f %10.4f | %10.4f %10.4f | %12.6f\n",
    cell, fit, param, d_cov95_t, d_se_ratio, d_t1_05, d_t1_05_paired_mcse,
    emp_sd_jk_over_pooled)))
}
