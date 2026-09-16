## EXACT pooled-vs-true surface for the 2-column vqtl W and the 3-column
## production mixture W, on a common DGP.  Zero Monte Carlo error.
##
## Identity this rests on (proved in verify_W.R C6, exact to 4e-16):
##     Sigma_true - Sigma_pooled = diag(1/f) E diag(1/f)
##     E[s,t] = Delta(min(t_s,t_t)) - t_s Delta(t_t) - t_t Delta(t_s)
##     Delta(tau) = sum_g w_g (F_g(q_tau) - tau),  w_g = p_g (g-Eg)^2 / Var(g)
## so for column j, with v_j = diag(1/f) W (W'W)^-1 e_j,
##     (SE_pooled / SE_true)^2 = 1 / (1 + v_j' E v_j / v_j' B v_j).
WS <- "/user/work/fh6520/gmoq_covmixW"
.libPaths(c(file.path(WS, "lib")))
stopifnot(identical(.libPaths()[1], file.path(WS, "lib")))
library(fungwasStage2)
source(file.path(WS, "src", "truth_mix.R"))

TAUS <- seq(0.02, 0.98, length.out = 45)
N    <- 461320L
args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args)) args[1] else file.path(WS, "out", "exact_mix.tsv")

cache <- new.env(parent = emptyenv())
PT <- function(maf, shape, a, cc) {
  k <- paste(maf, shape, a, cc)
  if (is.null(cache[[k]])) cache[[k]] <- pop_truth(TAUS, maf, shape, a, cc)
  cache[[k]]
}

## max |Z| over the columns of the reference W: the Z a real scan would report
maxZ <- function(maf, shape, a, cc, wname) {
  tr <- PT(maf, shape, a, cc)
  W  <- build_Ws(TAUS, tr)[[wname]]
  ex <- exact_W(W, tr, N * tr$Vg)
  max(abs(ex$Z_true))
}
solve_par <- function(f, Zt, hi0 = 1) {
  hi <- hi0
  while (f(hi) < Zt) { hi <- hi * 2; if (hi > 512) return(NA_real_) }
  uniroot(function(v) f(v) - Zt, lower = 1e-9, upper = hi, tol = 1e-10)$root
}

SHAPES <- c("normal", "mix2", "mix4", "mix5.89", "mix6")
MAFS   <- c(0.005, 0.01, 0.05, 0.30)
Z_DISP <- c(5, 10, 16, 30)
Z_LOC  <- c(5, 16, 30)

cells <- list()
addc <- function(...) cells[[length(cells) + 1L]] <<- list(...)
for (sh in SHAPES) {
  ref <- if (sh == "normal") "vqtl2" else "mix3"   # calibrate on the W that ships for that shape
  for (maf in MAFS) {
    addc(shape = sh, maf = maf, dgp = "null", a = 0, c = 0, Zt = 0, ref = ref)
    for (Z in Z_LOC) {
      a <- solve_par(function(v) maxZ(maf, sh, v, 0, ref), Z)
      if (is.finite(a)) addc(shape = sh, maf = maf, dgp = "loc", a = a, c = 0, Zt = Z, ref = ref)
    }
    for (Z in Z_DISP) {
      cc <- solve_par(function(v) maxZ(maf, sh, 0, v, ref), Z)
      if (is.finite(cc)) addc(shape = sh, maf = maf, dgp = "disp", a = 0, c = cc, Zt = Z, ref = ref)
    }
    cat("calibrated", sh, maf, "\n"); flush.console()
  }
}
cat("n cells =", length(cells), "\n")

rows <- list()
for (k in cells) {
  tr <- PT(k$maf, k$shape, k$a, k$c)
  Ws <- build_Ws(TAUS, tr)
  d  <- N * tr$Vg
  ## per-tau pooled/true ratio: the raw material, W-free
  ptr <- sqrt(diag(tr$Sigma0) / diag(tr$Sigma_true))
  for (wn in names(Ws)) {
    W  <- Ws[[wn]]
    ex <- exact_W(W, tr, d)
    ## exactness control: the identity above must reproduce `ratio` exactly
    V  <- ex$V
    vBv <- diag(t(V) %*% tr$Bbridge %*% V); vEv <- diag(t(V) %*% ex$Emat %*% V)
    rat_id <- sqrt(1 / (1 + vEv / vBv))
    for (j in seq_len(ncol(W))) {
      rows[[length(rows) + 1L]] <- data.frame(
        shape = k$shape, sep = .mix_sep(k$shape), maf = k$maf, dgp = k$dgp,
        Ztarget = k$Zt, a = k$a, c = k$c, sd_ratio_hom = 1 + 2 * k$c,
        MAC = 2 * k$maf * N, W = wn, K = ncol(W), param = colnames(W)[j],
        theta = ex$theta[j], Z_true = ex$Z_true[j],
        se_ratio_pooled = ex$ratio[j],
        se_ratio_identity = rat_id[j],
        cov95_pooled = ex$cov95[j],
        t1_05_pooled = ex$t1_05[j],
        t1_5e8_pooled = ex$t1_5e8[j],
        infl_5e8_pooled = ex$infl_5e8[j],
        tau_bar = ex$tau_bar[j], participation = ex$part[j],
        vEv_over_vBv = (vEv / vBv)[j],
        per_tau_ratio_min = min(ptr), per_tau_ratio_max = max(ptr),
        max_abs_Delta = max(abs(tr$Delta)),
        Delta_at_median = tr$Delta[which.min(abs(TAUS - 0.5))],
        f_min = min(tr$f), f_max = max(tr$f),
        stringsAsFactors = FALSE)
    }
  }
}
tab <- do.call(rbind, rows)
stopifnot(max(abs(tab$se_ratio_pooled - tab$se_ratio_identity)) < 1e-9)
cat("IDENTITY CONTROL PASS: max|ratio - identity| =",
    max(abs(tab$se_ratio_pooled - tab$se_ratio_identity)), "\n")
write.table(tab, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", OUT, nrow(tab), "rows\n")
