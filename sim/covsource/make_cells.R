## Build the cell grid for the cov_source head-to-head.
##
## Effect sizes are CALIBRATED, not guessed: for each (MAF, shape, target Z) the
## script solves for the a or c that makes the relevant parameter's TRUE Z equal
## the target, where "true Z" = theta_pop / (exact asymptotic SE of the shipped
## OLS estimator under the TRUE covariance), at d_resid = N * Var(g).
## Everything in this file is exact (root-finder tolerance); no Monte Carlo.
WS <- "/user/work/fh6520/covsource_run"
.libPaths(c(file.path(WS, "lib")))
stopifnot(identical(.libPaths()[1], file.path(WS, "lib")))
library(fungwasStage2)
source(file.path(WS, "src", "truth_cs.R"))

TAUS <- seq(0.02, 0.98, length.out = 45)
N    <- 461320L

Zof <- function(maf, shape, a, cc, which) {
  tr <- pop_truth(TAUS, maf, shape, a, cc)
  W  <- make_weights_vqtl(TAUS, tr$q_tau, mu = tr$EY, sd = tr$sdY)
  d  <- N * tr$Vg
  tt <- cs_theta_truth(W, tr, d)
  unname(tt$th_ols[which] / tt$se_ols_true[which])
}

## Expanding-bracket root find; Z is increasing but saturates, so the upper
## bound has to grow until it brackets (or we declare the target unreachable).
solve_par <- function(f, Ztarget, hi0 = 1) {
  hi <- hi0
  while (f(hi) < Ztarget) {
    hi <- hi * 2
    if (hi > 512) stop("Z target ", Ztarget, " unreachable below par = 512")
  }
  uniroot(function(v) f(v) - Ztarget, lower = 1e-8, upper = hi, tol = 1e-12)$root
}
solve_a <- function(maf, shape, Ztarget)
  solve_par(function(a) Zof(maf, shape, a, 0, "beta_mu"), Ztarget)
solve_c <- function(maf, shape, Ztarget)
  solve_par(function(cc) Zof(maf, shape, 0, cc, "beta_sigma2"), Ztarget)

cells <- list()
add <- function(name, dgp, maf, shape, a, cc, reps, Ztarget = NA) {
  cells[[length(cells) + 1L]] <<- list(name = name, dgp = dgp, maf = maf,
                                       shape = shape, a = a, c = cc,
                                       N = N, reps = reps, Ztarget = Ztarget)
}

MAFS <- c(0.30, 0.05, 0.01, 0.005)
REPS <- 15000L

## ---- A. main grid, normal, Z = 16 -----------------------------------------
for (maf in MAFS) {
  a16 <- solve_a(maf, "normal", 16)
  c16 <- solve_c(maf, "normal", 16)
  add(sprintf("A_null_maf%g_normal", maf),  "null",  maf, "normal", 0,   0,   REPS)
  add(sprintf("A_loc_maf%g_normal_Z16", maf),  "loc",   maf, "normal", a16, 0,   REPS, 16)
  add(sprintf("A_disp_maf%g_normal_Z16", maf), "disp",  maf, "normal", 0,   c16, REPS, 16)
  add(sprintf("A_mixed_maf%g_normal_Z16", maf),"mixed", maf, "normal", a16, c16, REPS, 16)
}

## ---- B. pure-dispersion effect ladder, normal ------------------------------
for (maf in MAFS) for (Z in c(5, 10, 30, 60)) {
  add(sprintf("B_disp_maf%g_normal_Z%d", maf, Z), "disp", maf, "normal",
      0, solve_c(maf, "normal", Z), REPS, Z)
}

## ---- C. pure-location ladder (the false-vQTL scenario), normal -------------
for (maf in c(0.30, 0.005)) for (Z in c(5, 30, 60)) {
  add(sprintf("C_loc_maf%g_normal_Z%d", maf, Z), "loc", maf, "normal",
      solve_a(maf, "normal", Z), 0, REPS, Z)
}

## ---- D. skewed phenotype (standardised gamma(2), skew 1.41) ----------------
for (maf in c(0.30, 0.01)) {
  a16 <- solve_a(maf, "gamma2", 16); c16 <- solve_c(maf, "gamma2", 16)
  add(sprintf("D_null_maf%g_gamma2", maf), "null", maf, "gamma2", 0, 0, REPS)
  add(sprintf("D_loc_maf%g_gamma2_Z16", maf), "loc", maf, "gamma2", a16, 0, REPS, 16)
  add(sprintf("D_disp_maf%g_gamma2_Z16", maf), "disp", maf, "gamma2", 0, c16, REPS, 16)
}

## ---- E. low-Z dispersion cells: the only regime where POWER discriminates.
## At Z >= 5 every configuration has power ~1 at nominal 0.05, so a power
## comparison there measures nothing.
for (maf in MAFS) for (Z in c(2, 3)) {
  add(sprintf("E_disp_maf%g_normal_Z%d", maf, Z), "disp", maf, "normal",
      0, solve_c(maf, "normal", Z), REPS, Z)
}

## ---- write ------------------------------------------------------------------
esc <- function(s) paste0("\"", s, "\"")
js <- vapply(cells, function(k) paste0(
  "{\"name\": ", esc(k$name), ", \"dgp\": ", esc(k$dgp),
  ", \"maf\": ", format(k$maf, digits = 10),
  ", \"shape\": ", esc(k$shape),
  ", \"a\": ", format(k$a, digits = 15),
  ", \"c\": ", format(k$c, digits = 15),
  ", \"N\": ", k$N, ", \"reps\": ", k$reps,
  ", \"Ztarget\": ", ifelse(is.na(k$Ztarget), "null", k$Ztarget), "}"), character(1))
## atomic write: jobs already running read cells.json by index, and the first
## 44 entries are reproduced bit-for-bit by this script, so appending is safe
## provided the replacement is atomic.
tmp <- file.path(WS, "cells.json.tmp")
writeLines(paste0("[\n", paste(js, collapse = ",\n"), "\n]"), tmp)
stopifnot(file.rename(tmp, file.path(WS, "cells.json")))

cat(sprintf("%-32s %-6s %6s %-7s %9s %9s %8s %8s %6s %8s\n",
            "name", "dgp", "maf", "shape", "a", "c", "Zmu", "Zs2", "sdrat", "MAC"))
for (i in seq_along(cells)) {
  k <- cells[[i]]
  zm <- Zof(k$maf, k$shape, k$a, k$c, "beta_mu")
  zs <- Zof(k$maf, k$shape, k$a, k$c, "beta_sigma2")
  cat(sprintf("%3d %-32s %-6s %6g %-7s %9.5f %9.5f %8.2f %8.2f %6.2f %8.0f\n",
              i - 1L, k$name, k$dgp, k$maf, k$shape, k$a, k$c, zm, zs,
              1 + 2 * k$c, 2 * k$maf * N))
}
cat("n_cells =", length(cells), " total reps =", sum(vapply(cells, function(k) k$reps, 1L)), "\n")
