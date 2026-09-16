## The DGPs the 3-column W actually EXISTS for.
##
## Y | g  ~  p1(g) N(mu1 + a1 g, (s (1+c g))^2) + (1-p1(g)) N(mu2 + a2 g, (s (1+c g))^2)
##          logit p1(g) = logit p1_0 + gam * g
##
## Five regimes, and which of the W's columns are NULL in each:
##   null   a1=a2=gam=c=0          all three
##   loc    a1=a2=a                gamma
##   disp   c>0                    gamma
##   comp   a1=0, a2=a             gamma AND beta_left   <- component-specific
##   memb   gam!=0                 beta_left AND beta_right <- membership effect
##
## The last two matter because the exact sweep shows gamma is protected and the
## component means are not.  If the null parameter is a component mean, the
## pooled arm's mis-calibration lands on a false-positive rate.
WS <- "/user/work/fh6520/gmoq_covmixW"
.libPaths(c(file.path(WS, "lib")))
stopifnot(identical(.libPaths()[1], file.path(WS, "lib")))
library(fungwasStage2)
source(file.path(WS, "src", "truth_mix.R"))
TAUS <- seq(0.02, 0.98, length.out = 45); N <- 461320L

## exact population truth for the 2-component DGP above
pop_truth_mix2 <- function(taus, maf, S, p1_0 = 0.5, a1 = 0, a2 = 0, gam = 0, cc = 0) {
  gv <- 0:2; pg <- dbinom(gv, 2, maf)
  Eg <- sum(pg * gv); Vg <- sum(pg * gv^2) - Eg^2; gc <- gv - Eg
  ## baseline (g = 0) components, standardised to mean 0 / var 1
  p2_0 <- 1 - p1_0
  k <- sqrt(1 + p1_0 * p2_0 * S^2)
  m1 <- -p2_0 * S / k; m2 <- p1_0 * S / k; s <- 1 / k
  ## per-genotype parameters
  ## c scales the WHOLE phenotype (both component means and both sds), which is
  ## exactly Y = a g + (1 + c g) eps in the published covsource DGP; a1/a2 shift
  ## the component means; gam tilts the membership log-odds.
  pig  <- plogis(qlogis(p1_0) + gam * gv)
  sc   <- 1 + cc * gv
  mu1g <- m1 * sc + a1 * gv; mu2g <- m2 * sc + a2 * gv
  sdg  <- s * sc
  Fcond <- function(y, i) pig[i] * pnorm(y, mu1g[i], sdg[i]) +
                          (1 - pig[i]) * pnorm(y, mu2g[i], sdg[i])
  fcond <- function(y, i) pig[i] * dnorm(y, mu1g[i], sdg[i]) +
                          (1 - pig[i]) * dnorm(y, mu2g[i], sdg[i])
  Fy <- function(y) sum(pg * vapply(seq_along(gv), function(i) Fcond(y, i), 0))
  fy <- function(y) sum(pg * vapply(seq_along(gv), function(i) fcond(y, i), 0))
  lo <- min(mu1g, mu2g) - 14 * max(sdg); hi <- max(mu1g, mu2g) + 14 * max(sdg)
  q <- vapply(taus, function(tt) uniroot(function(y) Fy(y) - tt, c(lo, hi), tol = 1e-13)$root, 0)
  f <- vapply(q, fy, 0)
  Fg <- outer(seq_along(gv), seq_along(taus), Vectorize(function(i, j) Fcond(q[j], i)))
  beta <- vapply(seq_along(taus), function(j) {
    -(sum(pg * gv * Fg[, j]) - Eg * taus[j]) / (f[j] * Vg) }, 0)
  ## marginal moments
  EYg  <- pig * mu1g + (1 - pig) * mu2g
  EY2g <- pig * (mu1g^2 + sdg^2) + (1 - pig) * (mu2g^2 + sdg^2)
  EY <- sum(pg * EYg); VY <- sum(pg * EY2g) - EY^2
  Tt <- length(taus)
  Bbr <- outer(taus, taus, pmin) - outer(taus, taus)
  S0 <- Bbr / outer(f, f)
  w <- pg * gc^2 / Vg
  Bw <- matrix(0, Tt, Tt)
  for (a in seq_len(Tt)) for (b in a:Tt) {
    mn <- if (q[a] <= q[b]) Fg[, a] else Fg[, b]
    e <- sum(w * (mn - taus[a] * Fg[, b] - taus[b] * Fg[, a] + taus[a] * taus[b]))
    Bw[a, b] <- e; Bw[b, a] <- e
  }
  Delta <- as.numeric(crossprod(w, Fg)) - taus
  list(taus = taus, q_tau = q, f = f, beta = beta,
       Sigma0 = S0, Sigma_true = Bw / outer(f, f),
       Bbridge = Bbr, Bw = Bw, Delta = Delta,
       EY = EY, sdY = sqrt(VY), Vg = Vg,
       base = list(p1 = p1_0, mu1 = m1, sd1 = s, mu2 = m2, sd2 = s))
}
## W from the baseline components, rescaled onto the marginal (what production
## gets by fitting the mixture to the observed phenotype)
Wmix3 <- function(tr, K3 = TRUE) {
  b <- tr$base
  make_weights_normal_mixture(tr$taus, tr$q_tau, p1 = b$p1,
    mu1 = tr$EY + tr$sdY * b$mu1, sd1 = tr$sdY * b$sd1,
    mu2 = tr$EY + tr$sdY * b$mu2, sd2 = tr$sdY * b$sd2, include_membership = K3)
}

## control: with a1=a2=0, gam=0 this must reproduce the mix<S> engine exactly
t1 <- pop_truth_mix2(TAUS, 0.01, 5.89, 0.5, 0, 0, 0, 0.2)
t2 <- pop_truth(TAUS, 0.01, "mix5.89", 0, 0.2)
cat(sprintf("CONTROL: pop_truth_mix2 == pop_truth(mix5.89)  rel err Sigma_true = %.3e, beta = %.3e\n",
            max(abs(t1$Sigma_true - t2$Sigma_true)) / max(abs(t2$Sigma_true)),
            max(abs(t1$beta - t2$beta)) / max(abs(t2$beta))))
stopifnot(max(abs(t1$Sigma_true - t2$Sigma_true)) / max(abs(t2$Sigma_true)) < 1e-10)

ZC8 <- qnorm(1 - 5e-8 / 2); ZC05 <- qnorm(1 - 0.05 / 2)
meas <- function(tr, K3 = TRUE) {
  W <- Wmix3(tr, K3); A <- solve(crossprod(W), t(W))
  d <- N * tr$Vg
  th <- drop(A %*% tr$beta)
  st <- sqrt(diag(A %*% (tr$Sigma_true / d) %*% t(A)))
  sp <- sqrt(diag(A %*% (tr$Sigma0 / d) %*% t(A)))
  r <- sp / st
  data.frame(param = colnames(W), K = ncol(W), theta = th, Z_true = th / st,
             se_ratio = r, cov95 = pnorm(ZC05 * r) - pnorm(-ZC05 * r),
             t1_05 = 2 * pnorm(-ZC05 * r), infl_5e8 = 2 * pnorm(-ZC8 * r) / 5e-8,
             stringsAsFactors = FALSE)
}
solve_for_Z <- function(maf, S, which, Zt) {
  f <- function(v) {
    ar <- list(a1 = 0, a2 = 0, gam = 0, cc = 0); ar[[which]] <- v
    m <- meas(pop_truth_mix2(TAUS, maf, S, 0.5, ar$a1, ar$a2, ar$gam, ar$cc))
    max(abs(m$Z_true))
  }
  hi <- 1e-3; while (f(hi) < Zt) { hi <- hi * 2; if (hi > 512) return(NA_real_) }
  uniroot(function(v) f(v) - Zt, c(1e-9, hi), tol = 1e-10)$root
}

cat("\n########## D1. component-specific and membership effects, S=5.89 ##########\n")
cat("EXACT.  'NULL?' marks parameters whose true value is zero: those are the\n")
cat("ones on which se_ratio < 1 becomes a FALSE-POSITIVE rate.\n\n")
cat(sprintf("%-6s %6s %4s | %-12s %5s %9s %9s %9s %12s\n",
            "dgp", "MAF", "Ztg", "param", "NULL?", "Z_true", "se_ratio", "cov95", "infl@5e-8"))
for (dg in c("loc", "disp", "comp", "memb")) {
  wh <- switch(dg, loc = "a1", disp = "cc", comp = "a2", memb = "gam")
  for (maf in c(0.01, 0.30)) for (Zt in c(16, 30)) {
    v <- solve_for_Z(maf, 5.89, wh, Zt)
    if (!is.finite(v)) { cat(sprintf("%-6s %6g %4g | unreachable\n", dg, maf, Zt)); next }
    ar <- list(a1 = 0, a2 = 0, gam = 0, cc = 0); ar[[wh]] <- v
    if (dg == "loc") ar$a2 <- v          # location = BOTH components shift
    tr <- pop_truth_mix2(TAUS, maf, 5.89, 0.5, ar$a1, ar$a2, ar$gam, ar$cc)
    m <- meas(tr)
    for (i in seq_len(nrow(m))) with(m[i, ], cat(sprintf(
      "%-6s %6g %4g | %-12s %5s %9.1f %9.4f %9.4f %12.4g\n",
      dg, maf, Zt, param, ifelse(abs(Z_true) < 0.75, "YES", ""), Z_true,
      se_ratio, cov95, infl_5e8)))
  }
}

cat("\n########## D2. the worst false-positive rate the production W can carry ##########\n")
cat("For each DGP, the largest 5e-8 inflation among the NULL parameters only.\n\n")
cat(sprintf("%-6s %6s %4s | %-12s %9s %12s | %-12s %9s %12s\n",
            "dgp", "MAF", "Ztg", "worst null par", "se_ratio", "infl@5e-8",
            "(K=2 W)", "se_ratio", "infl@5e-8"))
for (dg in c("loc", "disp", "comp", "memb")) {
  wh <- switch(dg, loc = "a1", disp = "cc", comp = "a2", memb = "gam")
  for (maf in c(0.005, 0.01, 0.05, 0.30)) {
    v <- solve_for_Z(maf, 5.89, wh, 30); if (!is.finite(v)) next
    ar <- list(a1 = 0, a2 = 0, gam = 0, cc = 0); ar[[wh]] <- v
    if (dg == "loc") ar$a2 <- v
    tr <- pop_truth_mix2(TAUS, maf, 5.89, 0.5, ar$a1, ar$a2, ar$gam, ar$cc)
    m3 <- meas(tr, TRUE); m2 <- meas(tr, FALSE)
    n3 <- m3[abs(m3$Z_true) < 0.75, ]; n2 <- m2[abs(m2$Z_true) < 0.75, ]
    r3 <- if (nrow(n3)) n3[which.max(n3$infl_5e8), ] else NULL
    r2 <- if (nrow(n2)) n2[which.max(n2$infl_5e8), ] else NULL
    cat(sprintf("%-6s %6g %4g | %-12s %9.4f %12.4g | %-12s %9s %12s\n", dg, maf, 30,
        if (is.null(r3)) "-" else r3$param, if (is.null(r3)) NA else r3$se_ratio,
        if (is.null(r3)) NA else r3$infl_5e8,
        if (is.null(r2)) "-" else r2$param,
        if (is.null(r2)) "-" else sprintf("%.4f", r2$se_ratio),
        if (is.null(r2)) "-" else sprintf("%.4g", r2$infl_5e8)))
  }
}
cat("\nDONE\n")
