## Predicted-vs-measured OLS/GLS projection penalty.
options(width = 200)
outdir <- "/user/work/fh6520/gmoq_projeff/out"
CELLS <- c("N2","N3","N4","N5","N6","K4","K6","MIX4","MIX6")
LAB <- c(N2="normal S=2", N3="normal S=3", N4="normal S=4", N5="normal S=5",
         N6="normal S=6", K4="skew a=-1 S=4", K6="skew a=-1 S=6",
         MIX4="skew a=+2 S=4", MIX6="skew a=+2 S=6")
rows <- NULL
for (C in CELLS) {
  f <- file.path(outdir, paste0("sim_", C, ".rds")); if (!file.exists(f)) next
  z <- readRDS(f); r <- z$res; n <- nrow(r)
  osd1 <- sd(r$or_b1); osd2 <- sd(r$or_b2)
  fin <- function(x) x[is.finite(x)]
  mk <- function(tag, b1, b2, s1, s2) data.frame(
    cell = C, model = LAB[[C]], arm = tag, n = n,
    sd_b1 = sd(fin(b1)), sd_b2 = sd(fin(b2)),
    eff1 = sd(fin(b1))/osd1, eff2 = sd(fin(b2))/osd2,
    meanSE1 = mean(fin(s1)), meanSE2 = mean(fin(s2)),
    cal1 = mean(fin(s1))/sd(fin(b1)), cal2 = mean(fin(s2))/sd(fin(b2)))
  a <- rbind(mk("2col OLS", r$o2_b1, r$o2_b2, r$o2_s1, r$o2_s2),
             mk("2col GLS", r$g2_b1, r$g2_b2, r$g2_s1, r$g2_s2),
             mk("3col OLS", r$o3_b1, r$o3_b2, r$o3_s1, r$o3_s2),
             mk("3col GLS", r$g3_b1, r$g3_b2, r$g3_s1, r$g3_s2))
  # predicted asymptotic penalty: median over reps of the per-rep analytic ratio
  a$pred2_1 <- median(fin(r$pr2_1)); a$pred2_2 <- median(fin(r$pr2_2))
  a$pred3_1 <- median(fin(r$pr3_1)); a$pred3_2 <- median(fin(r$pr3_2))
  rows <- rbind(rows, a)
}
saveRDS(rows, file.path(outdir, "sim_summary.rds"))

cat("\n================ MEASURED empirical SD ratio OLS/GLS vs PREDICTED ================\n")
cat(sprintf("%-16s %-9s %7s %8s %8s | %8s %8s %8s\n",
            "model","W","n","measOLS/GLS_b1","pred_b1","meas_b2","pred_b2",""))
for (C in unique(rows$cell)) {
  a <- rows[rows$cell == C, ]
  o2 <- a[a$arm=="2col OLS",]; g2 <- a[a$arm=="2col GLS",]
  o3 <- a[a$arm=="3col OLS",]; g3 <- a[a$arm=="3col GLS",]
  cat(sprintf("%-16s %-9s %7d %14.4f %8.4f | %8.4f %8.4f\n",
      a$model[1], "2col", a$n[1], o2$sd_b1/g2$sd_b1, a$pred2_1[1],
      o2$sd_b2/g2$sd_b2, a$pred2_2[1]))
  cat(sprintf("%-16s %-9s %7d %14.4f %8.4f | %8.4f %8.4f\n",
      "", "3col(+gamma)", a$n[1], o3$sd_b1/g3$sd_b1, a$pred3_1[1],
      o3$sd_b2/g3$sd_b2, a$pred3_2[1]))
}

cat("\n================ EFFICIENCY vs ORACLE (empirical SD / oracle SD) ================\n")
cat(sprintf("%-16s %-10s %8s %8s %8s %8s\n","model","arm","eff_b1","eff_b2","cal_b1","cal_b2"))
for (i in seq_len(nrow(rows))) with(rows[i,],
  cat(sprintf("%-16s %-10s %8.3f %8.3f %8.3f %8.3f\n", model, arm, eff1, eff2, cal1, cal2)))

cat("\n============ ABSOLUTE empirical SD: does +gamma cost anything? ============\n")
cat(sprintf("%-16s %10s %10s %10s %10s\n","model","2colOLS","2colGLS","3colOLS","3colGLS"))
for (C in unique(rows$cell)) {
  a <- rows[rows$cell == C, ]
  cat(sprintf("%-16s %10.6f %10.6f %10.6f %10.6f\n", a$model[1],
      a$sd_b1[a$arm=="2col OLS"], a$sd_b1[a$arm=="2col GLS"],
      a$sd_b1[a$arm=="3col OLS"], a$sd_b1[a$arm=="3col GLS"]))
}
