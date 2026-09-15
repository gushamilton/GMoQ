source("R/pooled-covariance.R")
run_id <- "pooled-memory-repair-20260915-01"
Tt <- 45L
set.seed(20260915)
X <- matrix(rnorm(Tt * Tt), Tt, Tt)
S <- crossprod(X) + diag(Tt)
mats <- list("22" = S)
rows <- lapply(c(1000L, 10000L), function(P) {
  chr <- rep(22L, P)
  d <- seq(10, 10 + P - 1)
  gc()
  td <- system.time(dense <- .pooled_cov_data(mats, chr, d, Tt))["elapsed"]
  dense_bytes <- as.numeric(object.size(dense$cov_vec)) +
    as.numeric(object.size(dense$offsets))
  rm(dense)
  gc()
  tc <- system.time(compact <- .pooled_compact_cov_data(mats, chr, d, Tt))["elapsed"]
  compact_bytes <- as.numeric(object.size(compact$cov_vec)) +
    as.numeric(object.size(compact$offsets)) +
    as.numeric(object.size(compact$cov_scales))
  out <- data.frame(run_id = run_id, tau_count = Tt, snps = P,
                    dense_bytes = dense_bytes,
                    compact_bytes = compact_bytes,
                    dense_bytes_per_snp = dense_bytes / P,
                    compact_bytes_per_snp = compact_bytes / P,
                    dense_elapsed_sec = unname(td),
                    compact_elapsed_sec = unname(tc),
                    reduction = dense_bytes / compact_bytes)
  rm(compact)
  gc()
  out
})
profile <- do.call(rbind, rows)
write.table(profile, "/user/work/fh6520/gmoq_finish_20260915/inputs/pooled-memory-repair-20260915-01/memory-profile.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(profile, stdout(), sep = "\t", row.names = FALSE, quote = FALSE)
