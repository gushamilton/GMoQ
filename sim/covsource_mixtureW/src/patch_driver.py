## Derive stage2_driver_mix.R from the PUBLISHED sim/covsource driver by five
## localised substitutions.  Nothing else changes: the shipped param_gwas(),
## the package's own .pooled_cov_data(), and the pooled/jackknife df
## conventions are all untouched.
import io
src = "/user/work/fh6520/covsource_run/src/stage2_driver_cs.R"
dst = "/user/work/fh6520/gmoq_covmixW/src/stage2_driver_mix.R"
s = io.open(src).read()
subs = [
    ('WS <- "/user/work/fh6520/covsource_run"',
     'WS <- "/user/work/fh6520/gmoq_covmixW"'),
    ('source(file.path(WS, "src", "truth_cs.R"))',
     'source(file.path(WS, "src", "truth_mix.R"))'),
    ('readLines(file.path(WS, "cells.json"), warn = FALSE)',
     'readLines(file.path(WS, "cells_mix.json"), warn = FALSE)'),
    ('fits <- cs_weights(TAUS, tr, make_weights_vqtl)',
     'fits <- build_Ws(TAUS, tr)   # vqtl2 (2 col) + mix2 (2 col) + mix3 (3 col, PRODUCTION)\n'
     'cat("W set:", paste(sprintf("%s(K=%d)", names(fits), vapply(fits, ncol, 0L)), collapse=" "), "\\n")'),
    ('truth_tbl <- lapply(fits, function(W) cs_theta_truth(W, tr, d_pop))',
     'truth_tbl <- lapply(fits, function(W) exact_W(W, tr, d_pop))'),
    ('N_CHUNKS_EXPECTED <- 20L', 'N_CHUNKS_EXPECTED <- 10L'),
]
for a, b in subs:
    assert s.count(a) == 1, (a, s.count(a))
    s = s.replace(a, b)
io.open(dst, "w").write(s)
print("wrote", dst)
