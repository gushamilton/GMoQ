## Extension of sim/covsource/truth_cs.R to a BIMODAL 2-component-normal eps,
## and to the 3-column production mixture/skew W.
##
## Everything here is exact (root-finder / closed form).  No Monte Carlo.
##
## DGP is UNCHANGED from the published covsource campaign:
##     Y = a*g + (1 + c*g) * eps ,   g ~ Binom(2, maf),  eps standardised.
## Only the LAW of eps changes: "normal", "gamma2" (as published), or
## "mixS<sep>" = standardised 50/50 two-normal mixture with component
## separation <sep> in units of the component sd.
##
## The W is what the whole measurement is about:
##   vqtl2 : make_weights_vqtl                       (2 cols, PUBLISHED campaign)
##   mix2  : make_weights_normal_mixture(FALSE)      (2 cols, means only)
##   mix3  : make_weights_normal_mixture(TRUE)       (3 cols, PRODUCTION)
## mix3 is analytically identical to scripts/stage2/prepare_matched_weights.R's
## build_skew_mixture() at alpha = 0; verified numerically in verify_W.R.

GAMMA_SHAPE <- c(gamma2 = 2, gamma8 = 8)

## Shape grammar for the bimodal law:
##   "mix<S>"          50/50 two-normal mixture, component separation S
##   "mix<S>p<p1>"     the same with mixing proportion p1 for the LOWER component
## Standardised: X ~ p1 N(-p2 S, 1) + p2 N(p1 S, 1)  (mean 0),
## Var(X) = 1 + p1 p2 S^2, eps = X / kappa.  The p1 != 0.5 cells exist to move
## the density trough OFF the median, which is the one thing that can break the
## membership column's protection (see mechanism.R).
.mix_parse <- function(shape) {
  if (!grepl("^mix", shape)) return(NULL)
  body <- sub("^mix", "", shape)
  if (grepl("p", body, fixed = TRUE)) {
    pp <- strsplit(body, "p", fixed = TRUE)[[1]]
    list(S = as.numeric(pp[1]), p1 = as.numeric(pp[2]))
  } else list(S = as.numeric(body), p1 = 0.5)
}
.mix_sep <- function(shape) { m <- .mix_parse(shape); if (is.null(m)) NA_real_ else m$S }
.mix_p1  <- function(shape) { m <- .mix_parse(shape); if (is.null(m)) NA_real_ else m$p1 }
.mix_comp <- function(shape) {
  m <- .mix_parse(shape); p1 <- m$p1; p2 <- 1 - p1; S <- m$S
  k <- sqrt(1 + p1 * p2 * S^2)
  list(p1 = p1, m1 = -p2 * S / k, m2 = p1 * S / k, s = 1 / k, S = S, kappa = k)
}

eps_cdf <- function(e, shape) {
  if (shape == "normal") return(pnorm(e))
  if (grepl("^mix", shape)) { z <- .mix_comp(shape)
    return(z$p1 * pnorm(e, z$m1, z$s) + (1 - z$p1) * pnorm(e, z$m2, z$s)) }
  kk <- GAMMA_SHAPE[[shape]]
  pgamma(kk + sqrt(kk) * e, shape = kk, scale = 1)
}
eps_pdf <- function(e, shape) {
  if (shape == "normal") return(dnorm(e))
  if (grepl("^mix", shape)) { z <- .mix_comp(shape)
    return(z$p1 * dnorm(e, z$m1, z$s) + (1 - z$p1) * dnorm(e, z$m2, z$s)) }
  kk <- GAMMA_SHAPE[[shape]]
  sqrt(kk) * dgamma(kk + sqrt(kk) * e, shape = kk, scale = 1)
}
eps_quantile <- function(p, shape) {
  if (shape == "normal") return(qnorm(p))
  if (grepl("^mix", shape)) { z <- .mix_comp(shape)
    lo <- min(z$m1, z$m2) - 14 * z$s; hi <- max(z$m1, z$m2) + 14 * z$s
    return(vapply(p, function(pp) uniroot(function(e) eps_cdf(e, shape) - pp,
                                          lower = lo, upper = hi, tol = 1e-13)$root,
                  numeric(1))) }
  kk <- GAMMA_SHAPE[[shape]]
  (qgamma(p, shape = kk, scale = 1) - kk) / sqrt(kk)
}

## standardised-eps component parameters (used to build the mixture W)
eps_mix_params <- function(shape) {
  z <- .mix_comp(shape)
  list(p1 = z$p1, mu1 = z$m1, sd1 = z$s, mu2 = z$m2, sd2 = z$s, S = z$S)
}

## ---------------------------------------------------------------------------
## pop_truth: identical algebra to truth_cs.R, with the generalised eps.
## Returns Sigma0 (what POOLED targets) and Sigma_true (what the JACKKNIFE
## targets), both on the "N*Vg" unit scale, plus the bridge decomposition
## B / B_w and the Delta(tau) profile that drives their difference.
## ---------------------------------------------------------------------------
pop_truth <- function(taus, maf, shape, a, cc) {
  gv <- 0:2
  pg <- dbinom(gv, 2, maf)
  sg <- 1 + cc * gv
  Eg <- sum(pg * gv); Vg <- sum(pg * gv^2) - Eg^2
  gc <- gv - Eg
  Fy <- function(y) sum(pg * eps_cdf((y - a * gv) / sg, shape))
  fy <- function(y) sum(pg * eps_pdf((y - a * gv) / sg, shape) / sg)
  lo <- min(a * gv + sg * eps_quantile(1e-10, shape)) - 1
  hi <- max(a * gv + sg * eps_quantile(1 - 1e-10, shape)) + 1
  q <- vapply(taus, function(tt)
    uniroot(function(y) Fy(y) - tt, lower = lo, upper = hi, tol = 1e-13)$root, numeric(1))
  f <- vapply(q, fy, numeric(1))
  Fg <- outer(seq_along(gv), seq_along(taus),
              Vectorize(function(i, j) eps_cdf((q[j] - a * gv[i]) / sg[i], shape)))
  beta <- vapply(seq_along(taus), function(i) {
    Eg1 <- sum(pg * gv * Fg[, i])
    -(Eg1 - Eg * taus[i]) / (f[i] * Vg)
  }, numeric(1))
  EY  <- a * Eg
  EY2 <- sum(pg * (a^2 * gv^2 + sg^2))
  VY  <- EY2 - EY^2

  Tt <- length(taus)
  ## bridge (pooled target, unit scale before the 1/f sandwich)
  Bbr <- outer(taus, taus, pmin) - outer(taus, taus)
  S0  <- Bbr / outer(f, f)

  ## gc^2-weighted bridge (true target)
  w <- pg * gc^2 / Vg                      # sums to 1
  Bw <- matrix(0, Tt, Tt)
  for (s in seq_len(Tt)) for (t in s:Tt) {
    mn <- if (q[s] <= q[t]) Fg[, s] else Fg[, t]
    e  <- sum(w * (mn - taus[s] * Fg[, t] - taus[t] * Fg[, s] + taus[s] * taus[t]))
    Bw[s, t] <- e; Bw[t, s] <- e
  }
  Strue <- Bw / outer(f, f)

  ## Delta(tau) = sum_g w_g (F_g(q_tau) - tau): the first-order driver of
  ## B_w - B.  Exactly:  B_w - B = Delta(min) - tau_s Delta_t - tau_t Delta_s
  ## + O(Delta^2); Delta == 0 under the null, so the two sources agree there.
  Delta <- as.numeric(crossprod(w, Fg)) - taus

  dimnames(S0) <- dimnames(Strue) <- list(paste0("tau", taus), paste0("tau", taus))
  list(taus = taus, q_tau = q, f = f, beta = beta, Sigma0 = S0, Sigma_true = Strue,
       Bbridge = Bbr, Bw = Bw, Delta = Delta,
       EY = EY, sdY = sqrt(VY), Eg = Eg, Vg = Vg, Fg = Fg,
       maf = maf, shape = shape, a = a, c = cc)
}

## ---------------------------------------------------------------------------
## the three weight matrices, all built at the MARGINAL q_tau (as production
## does: the mixture is fitted to the observed phenotype, not to a g-stratum)
## ---------------------------------------------------------------------------
build_Ws <- function(taus, tr) {
  out <- list()
  out$vqtl2 <- fungwasStage2::make_weights_vqtl(taus, tr$q_tau, mu = tr$EY, sd = tr$sdY)
  if (grepl("^mix", tr$shape)) {
    P <- eps_mix_params(tr$shape)
    ## rescale the standardised-eps components onto the marginal Y scale
    m1 <- tr$EY + tr$sdY * P$mu1; s1 <- tr$sdY * P$sd1
    m2 <- tr$EY + tr$sdY * P$mu2; s2 <- tr$sdY * P$sd2
    out$mix2 <- fungwasStage2::make_weights_normal_mixture(
      taus, tr$q_tau, p1 = P$p1, mu1 = m1, sd1 = s1, mu2 = m2, sd2 = s2,
      include_membership = FALSE)
    out$mix3 <- fungwasStage2::make_weights_normal_mixture(
      taus, tr$q_tau, p1 = P$p1, mu1 = m1, sd1 = s1, mu2 = m2, sd2 = s2,
      include_membership = TRUE)
  }
  out
}

## ---------------------------------------------------------------------------
## exact per-column surface for one W.  Under OLS the point estimate is
## bit-identical between covariance sources (A is formed before the se branch),
## so theta and Z_true are shared; only the reported SE differs.
## ---------------------------------------------------------------------------
ZC05 <- qnorm(1 - 0.05 / 2)
ZC8  <- qnorm(1 - 5e-8 / 2)

exact_W <- function(W, tr, d_resid) {
  A  <- solve(crossprod(W), t(W))
  S0 <- tr$Sigma0 / d_resid
  ST <- tr$Sigma_true / d_resid
  th <- drop(A %*% tr$beta)
  se_true   <- sqrt(diag(A %*% ST %*% t(A)))
  se_pooled <- sqrt(diag(A %*% S0 %*% t(A)))
  rat <- se_pooled / se_true                      # reported / actual
  ## v_j = diag(1/f) %*% A'_j : the direction in bridge coordinates.
  V <- diag(1 / tr$f) %*% t(A)
  ## concentration diagnostics for the mechanism
  wts <- V^2; wts <- sweep(wts, 2, colSums(wts), "/")
  tau_bar <- as.numeric(crossprod(wts, tr$taus))            # |v|^2-weighted mean tau
  part    <- 1 / colSums(wts^2) / length(tr$taus)           # participation ratio (0,1]
  ## alignment of v_j with the Delta-driven perturbation, exactly
  Dl <- tr$Delta
  Emat <- outer(seq_along(Dl), seq_along(Dl), Vectorize(function(s, t) {
    mn <- if (tr$q_tau[s] <= tr$q_tau[t]) Dl[s] else Dl[t]
    mn - tr$taus[s] * Dl[t] - tr$taus[t] * Dl[s]
  }))
  vBv  <- diag(t(V) %*% tr$Bbridge %*% V)
  vEv  <- diag(t(V) %*% Emat %*% V)
  ## GLS estimands, so the PUBLISHED summarise_cs.R runs against this object
  ## unchanged.  Sigma is near-singular for a bimodal f, so guard the solves.
  S0i <- tryCatch(solve(S0), error = function(e) NULL)
  STi <- tryCatch(solve(ST), error = function(e) NULL)
  nmz <- setNames(rep(NA_real_, ncol(W)), colnames(W))
  th_gls <- th_glsP <- se_gls_true <- nmz
  if (!is.null(STi)) {
    th_gls <- tryCatch(drop(solve(t(W) %*% STi %*% W, t(W) %*% STi %*% tr$beta)),
                       error = function(e) nmz)
    se_gls_true <- tryCatch(sqrt(diag(solve(t(W) %*% STi %*% W))), error = function(e) nmz)
  }
  if (!is.null(S0i)) {
    th_glsP <- tryCatch(drop(solve(t(W) %*% S0i %*% W, t(W) %*% S0i %*% tr$beta)),
                        error = function(e) nmz)
  }
  names(th_gls) <- names(th_glsP) <- names(se_gls_true) <- colnames(W)
  list(theta = th, se_true = se_true, se_pooled = se_pooled, ratio = rat,
       ## aliases consumed by sim/covsource/summarise_cs.R
       th_ols = th, se_ols_true = se_true, se_ols_pooled = se_pooled,
       th_gls = th_gls, th_gls_pooledW = th_glsP, se_gls_true = se_gls_true,
       ratio_pooled_ols = rat,
       Z_true = th / se_true,
       cov95 = pnorm(ZC05 * rat) - pnorm(-ZC05 * rat),
       t1_05 = 2 * pnorm(-ZC05 * rat),
       t1_5e8 = 2 * pnorm(-ZC8 * rat),
       infl_5e8 = 2 * pnorm(-ZC8 * rat) / 5e-8,
       tau_bar = tau_bar, part = part,
       first_order = -0.5 * vEv / vBv,   # predicted (ratio - 1) to first order
       V = V, Emat = Emat)
}
