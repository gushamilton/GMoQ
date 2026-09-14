## External replication check: our N = 30,000 rung against the completed
## gls_coverage study (different seed stream, same machinery and same truth).
## gls_coverage is read-only here.
args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]; outdir <- args[2]
.libPaths(c(LIB)); stopifnot(identical(.libPaths()[1], LIB))

new <- read.delim(file.path(outdir, "summary.tsv"), stringsAsFactors = FALSE)
old <- read.delim("/user/work/fh6520/gls_coverage/out/summary.tsv", stringsAsFactors = FALSE)

## map our cell names onto theirs: dgp_maf<maf>_<shape>_N30k
new30 <- new[new$N == 30000 & new$arm %in% c("A", "B"), ]
new30$oldcell <- paste0(new30$dgp, "_maf", ifelse(new30$maf == 0.3, "0.3", "0.01"),
                        "_", new30$shape, "_N30k")
m <- merge(new30, old[old$N == 30000, c("cell", "fit", "est", "param",
                                        "cov95", "se_ratio", "sd_z", "nrep")],
           by.x = c("oldcell", "fit", "est", "param"),
           by.y = c("cell", "fit", "est", "param"),
           suffixes = c(".new", ".old"))
cat("matched rows:", nrow(m), "\n\n")
d <- m$cov95.new - m$cov95.old
mc <- sqrt(m$cov95.new * (1 - m$cov95.new) / m$nrep.new +
           m$cov95.old * (1 - m$cov95.old) / m$nrep.old)
cat(sprintf("cov95  new-old: mean %+.5f  sd %.5f  max|diff| %.4f\n",
            mean(d), sd(d), max(abs(d))))
cat(sprintf("       standardised by MC error: mean %+.3f  sd %.3f  max|z| %.2f\n",
            mean(d / mc), sd(d / mc), max(abs(d / mc))))
cat(sprintf("       fraction |z| > 3: %.4f (expect ~0.003)\n", mean(abs(d / mc) > 3)))
ds <- m$se_ratio.new - m$se_ratio.old
cat(sprintf("se_ratio new-old: mean %+.5f  sd %.5f  max|diff| %.4f\n",
            mean(ds), sd(ds), max(abs(ds))))
dz <- m$sd_z.new - m$sd_z.old
cat(sprintf("sd(z)    new-old: mean %+.5f  sd %.5f  max|diff| %.4f\n",
            mean(dz), sd(dz), max(abs(dz))))

cat("\n-- worst 8 rows by |cov95 difference| --\n")
o <- order(-abs(d))[1:min(8, nrow(m))]
print(m[o, c("oldcell", "fit", "est", "param", "cov95.new", "cov95.old")],
      row.names = FALSE)
