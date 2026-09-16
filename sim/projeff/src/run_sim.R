## Simulation confirmation of the analytic OLS/GLS projection penalty.
## Reuses the validated Stage1/Stage2 reimplementation in core.R (copied from
## /user/work/fh6520/gmoq_mixture_sim, which reproduces the installed package
## bit-for-bit). Adds: (i) a GLS arm for the 3-column (+membership) W,
## (ii) a skew-normal generative family, (iii) higher separations.
source("/user/work/fh6520/gmoq_projeff/core_base.R")
suppressMessages(library(parallel))

args  <- commandArgs(trailingOnly = TRUE)
CELL  <- args[1]; NREP <- as.integer(args[2]); NCORE <- as.integer(args[3])

## ---- skew-normal sampling (stochastic representation; no `sn` needed) ----
rsn_ <- function(n, xi, omega, alpha) {
  d <- alpha/sqrt(1 + alpha^2)
  u0 <- rnorm(n); v <- rnorm(n)
  xi + omega * (d * abs(u0) + sqrt(1 - d^2) * v)
}

CELLS <- list(
  N2   = list(fam="norm", N=1e5, S=2, b1=0.07, b2=0.035),
  N3   = list(fam="norm", N=1e5, S=3, b1=0.07, b2=0.035),
  N4   = list(fam="norm", N=1e5, S=4, b1=0.07, b2=0.035),
  N5   = list(fam="norm", N=1e5, S=5, b1=0.07, b2=0.035),
  N6   = list(fam="norm", N=1e5, S=6, b1=0.07, b2=0.035),
  # skew family at matched separation: tests family-independence of the penalty
  K4   = list(fam="skew", N=1e5, S=4, b1=0.07, b2=0.035, al=-1.0),
  K6   = list(fam="skew", N=1e5, S=6, b1=0.07, b2=0.035, al=-1.0),
  # "mixture family" cell: bimodal antibody-like, skewed components, the shape
  # the shipped build_skew_mixture() W is built for
  MIX4 = list(fam="skew", N=1e5, S=4, b1=0.07, b2=0.035, al=2.0),
  MIX6 = list(fam="skew", N=1e5, S=6, b1=0.07, b2=0.035, al=2.0)
)
cf <- CELLS[[CELL]]; stopifnot(!is.null(cf))
PI1 <- 0.5; MAF <- 0.3
AL  <- if (is.null(cf$al)) 0 else cf$al

sim_one <- function(N) {
  g <- rbinom(N, 2L, MAF)
  z <- as.integer(runif(N) < PI1)
  if (cf$fam == "norm") {
    mu <- ifelse(z == 1L, 0 + cf$b1*g, cf$S + cf$b2*g)
    y  <- rnorm(N, mu, 1)
  } else {
    d <- AL/sqrt(1 + AL^2); sdsn <- sqrt(1 - 2*d^2/pi)   # omega = 1
    sh <- ifelse(z == 1L, cf$b1*g, cf$S*sdsn + cf$b2*g)
    y  <- sh + rsn_(N, 0, 1, AL)
  }
  list(g = g, z = z, y = y)
}

one_rep <- function(r) {
  set.seed(4200000L + 1000L*which(names(CELLS) == CELL) + r)
  d <- sim_one(cf$N); y <- d$y; g <- d$g
  o  <- fit_oracle(y, g, d$z)
  s1 <- stage1_rif(y, g, TAUS)
  mm <- marginal_mix(y)                      # marginal 2-normal EM, as production
  W2 <- make_W(TAUS, s1$q_tau, mm$p1, mm$mu1, mm$sd1, mm$mu2, mm$sd2, FALSE)
  W3 <- make_W(TAUS, s1$q_tau, mm$p1, mm$mu1, mm$sd1, mm$mu2, mm$sd2, TRUE)
  p2 <- gmoq_project(s1, W2)
  p3 <- gmoq_project(s1, W3)
  data.frame(rep = r,
    or_b1 = o[["b1"]], or_b2 = o[["b2"]], or_se1 = o[["se1"]], or_se2 = o[["se2"]],
    o2_b1 = p2$ols[1], o2_b2 = p2$ols[2], o2_s1 = p2$ols_se[1], o2_s2 = p2$ols_se[2],
    g2_b1 = p2$gls[1], g2_b2 = p2$gls[2], g2_s1 = p2$gls_se[1], g2_s2 = p2$gls_se[2],
    o3_gm = p3$ols[1], o3_b1 = p3$ols[2], o3_b2 = p3$ols[3],
    o3_sg = p3$ols_se[1], o3_s1 = p3$ols_se[2], o3_s2 = p3$ols_se[3],
    g3_gm = p3$gls[1], g3_b1 = p3$gls[2], g3_b2 = p3$gls[3],
    g3_sg = p3$gls_se[1], g3_s1 = p3$gls_se[2], g3_s2 = p3$gls_se[3],
    # per-replicate ANALYTIC ratio from this rep's own W and Sigma: the exact
    # asymptotic prediction, computed on exactly the matrices the estimator used
    pr2_1 = p2$ols_se[1]/p2$gls_se[1], pr2_2 = p2$ols_se[2]/p2$gls_se[2],
    pr3_1 = p3$ols_se[2]/p3$gls_se[2], pr3_2 = p3$ols_se[3]/p3$gls_se[3])
}

res <- mclapply(seq_len(NREP), function(r) try(one_rep(r), silent = TRUE), mc.cores = NCORE)
bad <- vapply(res, function(z) inherits(z, "try-error"), logical(1))
cat("cell", CELL, "failures:", sum(bad), "of", NREP, "\n")
res <- do.call(rbind, res[!bad])
saveRDS(list(res = res, cf = cf, AL = AL),
        file.path("/user/work/fh6520/gmoq_projeff/out", paste0("sim_", CELL, ".rds")))
cat("done", CELL, nrow(res), "reps\n")
