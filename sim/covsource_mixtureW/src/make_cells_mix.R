## Monte Carlo cell grid.  Effect sizes CALIBRATED by root-finding on the exact
## true Z (same convention as sim/covsource/make_cells.R), where "true Z" is
## max_j |theta_j / SE_true_j| over the columns of the W that SHIPS for that
## phenotype shape: vqtl2 for the unimodal normal cells (the published
## campaign's convention), mix3 for the bimodal mixture cells (production).
WS <- "/user/work/fh6520/gmoq_covmixW"
.libPaths(c(file.path(WS, "lib")))
stopifnot(identical(.libPaths()[1], file.path(WS, "lib")))
library(fungwasStage2)
source(file.path(WS, "src", "truth_mix.R"))

TAUS <- seq(0.02, 0.98, length.out = 45)
N    <- 461320L
REPS <- 5000L          # see README: MC error budget

maxZ <- function(maf, shape, a, cc, wname) {
  tr <- pop_truth(TAUS, maf, shape, a, cc)
  ex <- exact_W(build_Ws(TAUS, tr)[[wname]], tr, N * tr$Vg)
  max(abs(ex$Z_true))
}
solve_par <- function(f, Zt, hi0 = 1) {
  hi <- hi0
  while (f(hi) < Zt) { hi <- hi * 2; if (hi > 512) stop("unreachable Z=", Zt) }
  uniroot(function(v) f(v) - Zt, lower = 1e-9, upper = hi, tol = 1e-10)$root
}

cells <- list()
add <- function(name, dgp, maf, shape, a, cc, Zt = NA, ref = "mix3")
  cells[[length(cells) + 1L]] <<- list(name = name, dgp = dgp, maf = maf, shape = shape,
                                       a = a, c = cc, N = N, reps = REPS,
                                       Ztarget = Zt, ref = ref)

## ---- R. reproduction link to the PUBLISHED campaign (2-column vqtl W) ------
for (maf in c(0.01, 0.005)) {
  cc <- solve_par(function(v) maxZ(maf, "normal", 0, v, "vqtl2"), 30)
  add(sprintf("R_disp_normal_maf%g_Z30", maf), "disp", maf, "normal", 0, cc, 30, "vqtl2")
}

## ---- M. production regime: testosterone-like separation S = 5.89 ----------
SEP <- "mix5.89"
for (maf in c(0.005, 0.01, 0.05, 0.30)) {
  add(sprintf("M_null_maf%g", maf), "null", maf, SEP, 0, 0, NA)
  a16 <- solve_par(function(v) maxZ(maf, SEP, v, 0, "mix3"), 16)
  add(sprintf("M_loc_maf%g_Z16", maf), "loc", maf, SEP, a16, 0, 16)
  c30 <- solve_par(function(v) maxZ(maf, SEP, 0, v, "mix3"), 30)
  add(sprintf("M_disp_maf%g_Z30", maf), "disp", maf, SEP, 0, c30, 30)
}
## where does divergence begin: low-Z dispersion at the two extreme MAFs
for (maf in c(0.01, 0.30)) {
  c5 <- solve_par(function(v) maxZ(maf, SEP, 0, v, "mix3"), 5)
  add(sprintf("M_disp_maf%g_Z5", maf), "disp", maf, SEP, 0, c5, 5)
}

## ---- S. separation ladder at MAF 0.01, true dispersion ---------------------
## Matched on the PHYSICAL effect (the same c as the S=5.89 / MAF 0.01 / Z30
## cell), not on Z: at S=2 the 3-column W cannot reach Z=30 under a pure scale
## effect at any c <= 512, because it has no variance column.  Holding c fixed
## makes separation the only thing that varies.
c_ref <- cells[[which(vapply(cells, function(k) k$name, "") ==
                      "M_disp_maf0.01_Z30")]]$c
for (sh in c("mix2", "mix4", "mix6")) {
  add(sprintf("S_disp_%s_maf0.01_cmatched", sh), "disp", 0.01, sh, 0, c_ref, NA)
}

esc <- function(s) paste0("\"", s, "\"")
js <- vapply(cells, function(k) paste0(
  "{\"name\": ", esc(k$name), ", \"dgp\": ", esc(k$dgp),
  ", \"maf\": ", format(k$maf, digits = 10), ", \"shape\": ", esc(k$shape),
  ", \"a\": ", format(k$a, digits = 15), ", \"c\": ", format(k$c, digits = 15),
  ", \"N\": ", k$N, ", \"reps\": ", k$reps, ", \"ref\": ", esc(k$ref),
  ", \"Ztarget\": ", ifelse(is.na(k$Ztarget), "null", k$Ztarget), "}"), character(1))
writeLines(paste0("[\n", paste(js, collapse = ",\n"), "\n]"),
           file.path(WS, "cells_mix.json"))

cat(sprintf("%3s %-28s %-5s %-8s %6s %10s %10s %7s %8s\n",
            "idx", "name", "dgp", "shape", "maf", "a", "c", "sdrat", "MAC"))
for (i in seq_along(cells)) with(cells[[i]], {
  tr <- pop_truth(TAUS, maf, shape, a, c); Ws <- build_Ws(TAUS, tr); d <- N * tr$Vg
  zs <- unlist(lapply(names(Ws), function(nm) {
    e <- exact_W(Ws[[nm]], tr, d); setNames(e$Z_true, paste0(nm, ":", colnames(Ws[[nm]]))) }))
  cat(sprintf("%3d %-28s %-5s %-8s %6g %10.5f %10.5f %7.2f %8.0f  |Z| %s\n",
              i - 1L, name, dgp, shape, maf, a, c, 1 + 2 * c, 2 * maf * N,
              paste(sprintf("%s=%.1f", names(zs), zs), collapse = " ")))
})
cat("n_cells =", length(cells), " reps/cell =", REPS,
    " total reps =", length(cells) * REPS, "\n")
