suppressPackageStartupMessages(library(data.table))
OUT <- "/user/work/fh6520/gmoq_power_conserv/out"
for (ch in c(22, 21)) {
  D <- readRDS(file.path(OUT, sprintf("merged_chr%d.rds", ch)))
  cat("chr", ch, " n=", nrow(D), " aligned=", sum(D$align == 1), " flipped=", sum(D$align == -1), "\n")
}
D <- readRDS(file.path(OUT, "merged_chr21.rds"))
miss <- D[rp_f < 5e-8 & pL.cov >= 5e-8]
cat("\nchr21 beta_left: REGENIE 5e-8 hits G-MoQ(covariate) misses:", nrow(miss), "\n")
print(miss[, .(snp_id, bp, regenie_p = rp_f, regenie_z = abs(qnorm(rp_f/2)),
               gmoq_p = pL.cov, gmoq_z = abs(qnorm(pL.cov/2)),
               se_ratio = sL.cov / rs_f)])
D22 <- readRDS(file.path(OUT, "merged_chr22.rds"))
m2 <- D22[rp_m < 5e-8 & pR.cov >= 5e-8]
cat("\nchr22 beta_right misses:", nrow(m2), "\n")
print(m2[, .(snp_id, bp, regenie_z = abs(qnorm(rp_m/2)), gmoq_z = abs(qnorm(pR.cov/2)))])
