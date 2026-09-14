# engine.R -- EXACT population analysis of the spurious-dispersion artefact in G-MoQ Stage 2.
#
# No Monte Carlo anywhere in this file.  Everything is a population calculation:
# the only numerical approximations are (i) root-finding for the marginal quantiles
# and (ii) a fine trapezoid grid for the marginal moments.  Both are pushed to ~1e-10.
#
# Data-generating process, on a BASE scale X:
#     X | g = k   ~   mu0 + delta*k  +  s*exp(gamma*k) * E
# with E a standardised residual (normal / standardised chi2_df / t_df),
# g ~ Binomial(2, p) under HWE.
#   delta = per-allele LOCATION effect on the base scale
#   gamma = per-allele LOG-SD effect on the base scale (gamma = 0 => NO true dispersion effect)
#
# The analysed phenotype is Y = T(X) for a monotone increasing T:
#   "identity"  -- analyse X itself
#   "exp"       -- analyse exp(X); with normal E this is the raw lognormal
#   "rint"      -- rank-inverse-normal transform of X (identical to RINT of exp(X),
#                  since RINT is invariant to any monotone transform)
#
# KEY IDENTITY that makes this cheap and exact.  The RIF-regression population
# estimand at quantile tau is
#     b(tau) = -Cov(1{Y <= q_tau}, g) / ( f_Y(q_tau) * Var(g) )
# Under a monotone T, q_tau^Y = T(q_tau^X) and F_k(q_tau^Y) = F_k(q_tau^X) exactly,
# while f_Y(q^Y) = f_X(q^X) / T'(q^X).  Hence
#     b_Y(tau) = b_X(tau) * T'(q_tau^X).
# So a single base-scale solve serves every transform.

suppressWarnings(suppressMessages({}))
options(warn = 1)

## ----------------------------------------------------------------- residual laws
make_base <- function(kind = c("normal", "chi2", "t"), df = 4) {
  kind <- match.arg(kind)
  if (kind == "normal") {
    return(list(F = pnorm, f = dnorm, lo = -40, hi = 40,
                skew = 0, kurt = 0, label = "normal"))
  }
  if (kind == "chi2") {
    # standardised chi-squared with df degrees of freedom: mean 0, sd 1
    m <- df; v <- 2 * df; sdv <- sqrt(v)
    list(F = function(e) pchisq(e * sdv + m, df = df),
         f = function(e) dchisq(e * sdv + m, df = df) * sdv,
         lo = -m / sdv, hi = 60,
         skew = sqrt(8 / df), kurt = 12 / df,
         label = paste0("chi2_", df))
  } else {
    # standardised Student t with df degrees of freedom (df > 2)
    sdv <- sqrt(df / (df - 2))
    list(F = function(e) pt(e * sdv, df = df),
         f = function(e) dt(e * sdv, df = df) * sdv,
         lo = -80, hi = 80,
         skew = 0, kurt = if (df > 4) 6 / (df - 4) else Inf,
         label = paste0("t_", df))
  }
}

## --------------------------------- baseline (delta = gamma = 0) marginal of Y = T(X)
gmoq_null_marginal <- function(maf, s, mu0, base, transform, taus, ngrid) {
  lo <- mu0 + s * base$lo; hi <- mu0 + s * base$hi
  Fm <- function(x) base$F((x - mu0) / s)
  fm <- function(x) base$f((x - mu0) / s) / s
  qx <- vapply(taus, function(tt)
    stats::uniroot(function(x) Fm(x) - tt, lower = lo, upper = hi,
                   tol = .Machine$double.eps^0.75)$root, 0.0)
  xs <- seq(lo, hi, length.out = ngrid); dx <- xs[2] - xs[1]
  Fx <- Fm(xs)
  Ty <- switch(transform,
               identity = xs, exp = exp(xs),
               rint = qnorm(pmin(pmax(Fx, 1e-300), 1 - 1e-16)))
  d <- fm(xs)
  trap <- function(v) { v[1] <- v[1] / 2; v[ngrid] <- v[ngrid] / 2; sum(v) * dx }
  EY <- trap(Ty * d); EY2 <- trap(Ty * Ty * d)
  qy <- switch(transform, identity = qx, exp = exp(qx), rint = qnorm(taus))
  list(qy = qy, EY = EY, VY = EY2 - EY^2)
}

## ----------------------------------------------------------------- the engine
gmoq_exact <- function(maf, delta, gamma = 0, s = 1, mu0 = 0,
                       base = make_base("normal"),
                       transform = c("identity", "exp", "rint"),
                       N = 200000,
                       taus = seq(0.02, 0.98, length.out = 45),
                       cov_mode = c("pooled", "true"),
                       fit = c("gls", "ols"),
                       W_mode = c("vqtl", "distfree"),
                       W_from = c("alt", "null"),
                       W_extra = c("none", "kurt"),
                       ngrid = 400001) {
  transform <- match.arg(transform)
  fit       <- match.arg(fit)
  W_mode    <- match.arg(W_mode)
  W_from    <- match.arg(W_from)
  W_extra   <- match.arg(W_extra)
  cov_mode  <- match.arg(cov_mode)

  p  <- maf
  w  <- c((1 - p)^2, 2 * p * (1 - p), p^2)     # HWE weights for g = 0,1,2
  gv <- 0:2
  gb <- sum(w * gv); gc <- gv - gb
  Vg <- sum(w * gc^2)                          # = 2p(1-p)
  gskew <- sum(w * gc^3) / Vg^1.5              # = (1-2p)/sqrt(2p(1-p))

  mk <- mu0 + delta * gv                       # per-genotype base-scale location
  sk <- s * exp(gamma * gv)                    # per-genotype base-scale scale

  Fk <- function(x) sapply(1:3, function(k) base$F((x - mk[k]) / sk[k]))
  fk <- function(x) sapply(1:3, function(k) base$f((x - mk[k]) / sk[k]) / sk[k])
  Fmix <- function(x) as.numeric(Fk(x) %*% w)
  fmix <- function(x) as.numeric(fk(x) %*% w)

  # --- base-scale marginal quantiles by bisection (vectorised, ~1e-12)
  lo <- min(mk + sk * base$lo); hi <- max(mk + sk * base$hi)
  qx <- vapply(taus, function(tt)
    stats::uniroot(function(x) Fmix(x) - tt, lower = lo, upper = hi,
                   tol = .Machine$double.eps^0.75)$root, 0.0)

  Fk_q <- t(vapply(qx, function(x) as.numeric(Fk(x)), numeric(3)))  # T x 3
  fx_q <- vapply(qx, fmix, 0.0)

  # --- transform and its derivative at the quantiles
  Tfun <- switch(transform,
                 identity = function(x) x,
                 exp      = function(x) exp(x),
                 rint     = function(x) qnorm(Fmix(x)))
  dT_q <- switch(transform,
                 identity = rep(1, length(qx)),
                 exp      = exp(qx),
                 rint     = fx_q / dnorm(qnorm(vapply(qx, Fmix, 0.0))))

  qy   <- vapply(qx, Tfun, 0.0)
  fy_q <- fx_q / dT_q

  # --- POPULATION RIF-regression estimand  b(tau)   (exact)
  covterm <- as.numeric(Fk_q %*% (w * gc))      # Cov(1{Y<=q}, g)
  b <- -covterm / (fy_q * Vg)

  # --- marginal + conditional moments of Y, by fine trapezoid on the base scale
  xs <- seq(lo, hi, length.out = ngrid); dx <- xs[2] - xs[1]
  # Fmix is vectorised over x (Fk returns length(x) x 3), so call it once on the whole
  # grid.  Calling it pointwise here costs ~1e6 R-level closure invocations and was the
  # single dominant cost in the RINT cells.
  Ty <- switch(transform,
               identity = xs,
               exp      = exp(xs),
               rint     = qnorm(pmin(pmax(Fmix(xs), 1e-300), 1 - 1e-16)))
  dens_k <- vapply(1:3, function(k) base$f((xs - mk[k]) / sk[k]) / sk[k], numeric(ngrid))
  trap <- function(v) { v[1] <- v[1] / 2; v[ngrid] <- v[ngrid] / 2; sum(v) * dx }
  EY_k  <- vapply(1:3, function(k) trap(Ty      * dens_k[, k]), 0.0)
  EY2_k <- vapply(1:3, function(k) trap(Ty * Ty * dens_k[, k]), 0.0)
  VY_k  <- EY2_k - EY_k^2
  EY    <- sum(w * EY_k)
  VY    <- sum(w * EY2_k) - EY^2
  sdY   <- sqrt(VY)

  # --- TRUTH on the Y scale: the linear-in-g contrasts of the conditional moments
  true_mu     <- sum(w * gc * EY_k) / Vg
  true_sigma2 <- sum(w * gc * VY_k) / Vg

  # --- the law the design matrix is built from.
  # "alt"  = the realised phenotype marginal (what Stage 1 actually sees: it computes
  #          quantiles and densities from the observed phenotype, SNP effect included).
  # "null" = the delta=gamma=0 baseline (what sim/dispersion's design_W uses).
  if (W_from == "alt") {
    qW <- qy; mW <- EY; vW <- VY; sdW <- sdY
  } else {
    n0 <- gmoq_null_marginal(maf, s, mu0, base, transform, taus, ngrid)
    qW <- n0$qy; mW <- n0$EY; vW <- n0$VY; sdW <- sqrt(n0$VY)
  }

  # --- weight matrix.
  # "vqtl"     = exactly stage2/R/weight-functions.R::make_weights_vqtl, i.e. z_tau/(2 sd).
  # "distfree" = the distribution-free generalisation (q_tau - m)/(2 v); identical for a
  #              normal baseline, different once the analysed scale is skewed.
  if (W_mode == "vqtl") {
    W <- cbind(beta_mu = rep(1, length(taus)), beta_sigma2 = qnorm(taus) / (2 * sdW))
  } else {
    W <- cbind(beta_mu = rep(1, length(taus)), beta_sigma2 = (qW - mW) / (2 * vW))
  }
  # LSK adds the kurtosis column z^2 - 1, exactly as gls_coverage/src/stage2_driver.R does
  if (W_extra == "kurt") W <- cbind(W, beta_kurt = qnorm(taus)^2 - 1)

  # --- covariance of the Stage 1 tau-slopes
  qmin <- outer(taus, taus, pmin)
  Sp <- (qmin - outer(taus, taus)) / outer(fy_q, fy_q) / (N * Vg)
  idx <- outer(seq_along(taus), seq_along(taus), function(i, j) pmin(i, j))
  St <- matrix(0, length(taus), length(taus))
  for (k in 1:3) {
    Fq <- Fk_q[, k]
    St <- St + w[k] * gc[k]^2 * (matrix(Fq[idx], length(taus)) - outer(Fq, Fq))
  }
  St <- St / outer(fy_q, fy_q) / (N * Vg^2)

  Sfit <- if (cov_mode == "pooled") Sp else St

  if (fit == "gls") {
    Si  <- solve(Sfit)
    A   <- solve(t(W) %*% Si %*% W)
    P   <- A %*% t(W) %*% Si
    Vrep <- A
  } else {
    P    <- solve(t(W) %*% W) %*% t(W)
    Vrep <- P %*% Sfit %*% t(P)
  }
  Vtrue <- P %*% St %*% t(P)

  theta   <- as.numeric(P %*% b)
  se_rep  <- sqrt(diag(Vrep))
  sd_true <- sqrt(diag(Vtrue))

  bias2 <- theta[2] - true_sigma2
  bias1 <- theta[1] - true_mu

  # lack-of-fit statistic Q (population non-centrality), GLS only
  r  <- as.numeric(b - W %*% theta)
  Qn <- if (fit == "gls") as.numeric(t(r) %*% solve(Sfit) %*% r) else NA_real_

  Zmu <- theta[1] / se_rep[1]
  G   <- 0.6855 * gskew * Zmu^2 / sqrt(N)

  list(taus = taus, b = b, W = W, qy = qy, fy = fy_q,
       theta = theta, se_rep = se_rep, sd_true = sd_true,
       true_mu = true_mu, true_sigma2 = true_sigma2,
       bias_mu = bias1, bias_sigma2 = bias2,
       ncp_sigma2 = bias2 / se_rep[2], ncp_mu = bias1 / se_rep[1],
       Z_sigma2 = theta[2] / se_rep[2], Z_mu = Zmu,
       se_ratio = se_rep / sd_true, EY = EY, sdY = sdY,
       VY_k = VY_k, EY_k = EY_k, Vg = Vg, gskew = gskew, G = G,
       Q_ncp = Qn, N = N, df_Q = length(taus) - ncol(W))
}

## rejection rate for beta_sigma2 at level alpha, exact given (ncp, ratio)
rej_rate <- function(theta2_pop, se_rep, sd_true, alpha) {
  z <- qnorm(1 - alpha / 2)
  ncp <- theta2_pop / se_rep; r <- sd_true / se_rep
  pnorm((-z - ncp) / r) + pnorm(-(z - ncp) / r)
}

## solve for the delta at which the spurious |Z_sigma2| hits a target
solve_delta <- function(target_absZ = 5.4513, ..., lo = 1e-5, hi = 2) {
  f <- function(d) abs(gmoq_exact(delta = d, ...)$Z_sigma2) - target_absZ
  if (f(hi) < 0) return(NA_real_)
  if (f(lo) > 0) return(lo)
  stats::uniroot(f, c(lo, hi), tol = 1e-9)$root
}
