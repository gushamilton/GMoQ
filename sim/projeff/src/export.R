o <- "/user/work/fh6520/gmoq_projeff/out"
d <- readRDS(file.path(o, "diag.rds"))
write.table(d, file.path(o, "table1_diagnostic_by_model.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
p <- as.data.frame(readRDS(file.path(o, "pred_mixture.rds")))
write.table(p, file.path(o, "table2_mixture_family_prediction.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
s <- readRDS(file.path(o, "sim_summary.rds"))
write.table(s, file.path(o, "table3_simulation_measured.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("rows:", nrow(d), nrow(p), nrow(s), "\n")
cat("spearman(rho_max, worst):",
    { a <- d[d$algebra_check < 1e-6, ]; a$w <- ave(a$se_ratio, a$model, FUN = max)
      u <- unique(a[, c("rho_max","w")]); round(cor(u$rho_max, u$w, method="spearman"), 4) }, "\n")
