## Does the shipped Q lack-of-fit statistic flag the cells whose coverage
## degrades?  If the damage grows with N_eff, the diagnostic that measures it
## should grow too -- i.e. the guard should get MORE reliable at production N,
## not less.  Q_pval is per replicate, so this is a power curve.
args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]; outdir <- args[2]
.libPaths(c(LIB)); stopifnot(identical(.libPaths()[1], LIB))

files <- sort(list.files(outdir, pattern = "^stage2_cell[0-9]+\\.rds$", full.names = TRUE))
rows <- list()
for (f in files) {
  x <- readRDS(f)
  for (fn in names(x$W)) for (est in c("ols", "gls")) {
    qp <- x$res[[paste(fn, est, sep = "_")]]$Qp
    if (is.null(qp)) next
    rows[[length(rows) + 1L]] <- data.frame(
      arm = x$cell$arm, dgp = x$cell$dgp, shape = x$cell$shape, N = x$N,
      mac = 2 * x$N * x$cell$maf, fit = fn, est = est,
      rej05 = mean(qp < 0.05, na.rm = TRUE),
      rej1e3 = mean(qp < 1e-3, na.rm = TRUE),
      rej5e8 = mean(qp < 5e-8, na.rm = TRUE),
      nrep = sum(!is.na(qp)), stringsAsFactors = FALSE)
  }
}
q <- do.call(rbind, rows)
Ns <- c(30000, 1e5, 461320)

cat("=== Q rejection rate at p < 0.05.  dgp = null is the FALSE-POSITIVE rate;\n")
cat("    dgp = scale is the POWER to notice the model does not fit. ===\n\n")
for (est in c("ols", "gls")) {
  cat("---- estimator:", est, "----\n")
  cat(sprintf("%-6s %-7s %-5s %11s %11s %11s\n", "arm", "dgp", "fit",
              "N=30000", "N=1e+05", "N=461320"))
  for (arm in c("A", "B", "C")) for (dg in c("null", "scale")) for (fi in c("ADD", "LS", "LSK", "TRU")) {
    s <- q[q$arm == arm & q$dgp == dg & q$fit == fi & q$est == est, ]
    if (!nrow(s)) next
    cat(sprintf("%-6s %-7s %-5s", arm, dg, fi))
    for (N in Ns) cat(sprintf(" %11.4f", mean(s[s$N == N, "rej05"])))
    cat("\n")
  }
  cat("\n")
}
cat("=== same, at the genome-wide threshold p < 5e-8 (est = ols) ===\n")
cat(sprintf("%-6s %-7s %-5s %11s %11s %11s\n", "arm", "dgp", "fit",
            "N=30000", "N=1e+05", "N=461320"))
for (arm in c("A", "B", "C")) for (dg in c("null", "scale")) for (fi in c("ADD", "LS", "TRU")) {
  s <- q[q$arm == arm & q$dgp == dg & q$fit == fi & q$est == "ols", ]
  if (!nrow(s)) next
  cat(sprintf("%-6s %-7s %-5s", arm, dg, fi))
  for (N in Ns) cat(sprintf(" %11.4f", mean(s[s$N == N, "rej5e8"])))
  cat("\n")
}
