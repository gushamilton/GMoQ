.libPaths(c('/user/work/fh6520/gmoq_finish_20260915/root/lib',
            '/user/work/fh6520/r_packages', .libPaths()))
library(fungwasStage2)
f <- get('.pooled_cov_data', asNamespace('fungwasStage2'))
T <- 45L
S <- diag(T)
for (m in c(1000L, 10000L)) {
  tm <- system.time(cd <- f(list(ALL=S), rep('1', m), rep(1000, m), T))
  cat(sprintf('variants=%d elapsed=%.3f cov_payload_bytes=%d total_object_bytes=%d\n',
              m, tm[['elapsed']], as.numeric(object.size(cd$cov_vec)),
              as.numeric(object.size(cd))))
  rm(cd); gc()
}
