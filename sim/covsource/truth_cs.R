## Exact population truth for the cov_source head-to-head.
##
## Extends gls_coverage/src/truth.R with the TRUE asymptotic covariance of the
## Stage 1 beta_tau vector, which is what "pooled" approximates and what the
## block jackknife estimates.  Everything here is exact to the root-finder /
## closed-form tolerance -- no Monte Carlo.
##
## Model:  Y = a*g + (1 + c*g)*eps ,  g ~ Binom(2, maf) ,  eps standardised.
##
## Stage 1 estimand (RIF / unconditional quantile partial effect):
##   beta_tau = Cov(g, RIF_tau(Y)) / Var(g)
##
## Influence function of beta_tau (ignoring the O_p(N^-1) contribution from
## estimating q_tau and f_tau, which is the same object Stage 1's own jackknife
## picks up):
##   IF_i(tau) = gc_i * psi_tau(Y_i) / Vg ,
##   psi_tau(y) = -(1{y <= q_tau} - tau) / f_tau
## so
##   N * Cov(beta_s, beta_t) = E[ gc^2 psi_s psi_t ] / Vg^2 .
##
## The POOLED estimator instead uses Sigma_0 / d_resid with
##   Sigma_0[s,t] = E[psi_s psi_t] = (min(t_s,t_t) - t_s t_t)/(f_s f_t)
##   d_resid      = sum_i gc_i^2  ->  N * Vg
## i.e. it replaces the gc^2-weighted expectation by the unweighted one.  Under
## the null those agree exactly; under a dispersion effect they do not, and the
## ratio of the two is the whole of the pooled estimator's working-independence
## error.

GAMMA_SHAPE <- c(gamma2 = 2, gamma8 = 8)

eps_cdf <- function(e, shape) {
  if (shape == "normal") return(pnorm(e))
  k <- GAMMA_SHAPE[[shape]]
  pgamma(k + sqrt(k) * e, shape = k, scale = 1)
}
eps_pdf <- function(e, shape) {
  if (shape == "normal") return(dnorm(e))
  k <- GAMMA_SHAPE[[shape]]
  sqrt(k) * dgamma(k + sqrt(k) * e, shape = k, scale = 1)
}
eps_quantile <- function(p, shape) {
  if (shape == "normal") return(qnorm(p))
  k <- GAMMA_SHAPE[[shape]]
  (qgamma(p, shape = k, scale = 1) - k) / sqrt(k)
}

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
  ## per-genotype CDF at each population quantile: Fg[g, tau]
  Fg <- outer(seq_along(gv), seq_along(taus),
              Vectorize(function(i, j) eps_cdf((q[j] - a * gv[i]) / sg[i], shape)))
  beta <- vapply(seq_along(taus), function(i) {
    Eg1 <- sum(pg * gv * Fg[, i])
    -(Eg1 - Eg * taus[i]) / (f[i] * Vg)
  }, numeric(1))
  EY  <- a * Eg
  EY2 <- sum(pg * (a^2 * gv^2 + sg^2))
  VY  <- EY2 - EY^2

  ## pooled population target: Sigma_0 (unit scale; divided by d_resid downstream)
  S0 <- (outer(taus, taus, pmin) - outer(taus, taus)) / outer(f, f)

  ## TRUE covariance, unit scale: N * Vg * Cov(beta) = E[gc^2 psi psi'] / Vg
  ## E[(1{Y<=q_s}-t_s)(1{Y<=q_t}-t_t) | g] = Fg(min) - t_s Fg(q_t) - t_t Fg(q_s) + t_s t_t
  Tt <- length(taus)
  Strue <- matrix(0, Tt, Tt)
  for (s in seq_len(Tt)) for (t in s:Tt) {
    mn <- if (q[s] <= q[t]) Fg[, s] else Fg[, t]
    e  <- sum(pg * gc^2 * (mn - taus[s] * Fg[, t] - taus[t] * Fg[, s] + taus[s] * taus[t]))
    v  <- e / (f[s] * f[t] * Vg)
    Strue[s, t] <- v; Strue[t, s] <- v
  }
  dimnames(S0) <- dimnames(Strue) <- list(paste0("tau", taus), paste0("tau", taus))

  list(taus = taus, q_tau = q, f = f, beta = beta, Sigma0 = S0, Sigma_true = Strue,
       EY = EY, sdY = sqrt(VY), Eg = Eg, Vg = Vg, Fg = Fg,
       maf = maf, shape = shape, a = a, c = cc)
}

## Weight matrices used throughout.  W_LS is the SHIPPED make_weights_vqtl
## (normal baseline location+scale); the package's default vQTL model.
cs_weights <- function(taus, tr, make_weights_vqtl) {
  Tt <- length(taus)
  W_ADD <- matrix(1, Tt, 1, dimnames = list(paste0("tau", taus), "beta_add"))
  W_LS  <- make_weights_vqtl(taus, tr$q_tau, mu = tr$EY, sd = tr$sdY)
  list(ADD = W_ADD, LS = W_LS)
}

## Population theta and SEs for one W, both estimators, both covariance sources.
## d_resid enters only as a scalar: Sigma_pooled = Sigma0 / d_resid,
## Sigma_true_v = Sigma_true / d_resid  (both unit matrices are "N*Vg" scaled).
cs_theta_truth <- function(W, tr, d_resid) {
  A  <- solve(crossprod(W), t(W))
  S0 <- tr$Sigma0 / d_resid
  ST <- tr$Sigma_true / d_resid
  S0i <- solve(S0); STi <- solve(ST)
  th_ols <- drop(A %*% tr$beta)
  th_gls <- drop(solve(t(W) %*% STi %*% W, t(W) %*% STi %*% tr$beta))
  th_gls_pooledW <- drop(solve(t(W) %*% S0i %*% W, t(W) %*% S0i %*% tr$beta))
  names(th_ols) <- names(th_gls) <- names(th_gls_pooledW) <- colnames(W)
  ## OLS: theta = A beta, so Var = A Sigma A'
  se_ols_true   <- sqrt(diag(A %*% ST %*% t(A)))
  se_ols_pooled <- sqrt(diag(A %*% S0 %*% t(A)))          # what pooled+OLS reports
  ## GLS with the correct Sigma
  se_gls_true   <- sqrt(diag(solve(t(W) %*% STi %*% W)))
  ## GLS run with Sigma_0 as the weight: theta = B beta, B = (W'S0^-1 W)^-1 W'S0^-1
  B <- solve(t(W) %*% S0i %*% W, t(W) %*% S0i)
  se_gls_pooled_rep <- sqrt(diag(solve(t(W) %*% S0i %*% W)))   # reported
  se_gls_pooled_act <- sqrt(diag(B %*% ST %*% t(B)))           # actual
  list(th_ols = th_ols, th_gls = th_gls, th_gls_pooledW = th_gls_pooledW,
       se_ols_true = se_ols_true, se_ols_pooled = se_ols_pooled,
       se_gls_true = se_gls_true,
       se_gls_pooled_rep = se_gls_pooled_rep,
       se_gls_pooled_act = se_gls_pooled_act,
       ## exact SE-ratio (reported / actual) for each configuration
       ratio_pooled_ols = se_ols_pooled / se_ols_true,
       ratio_pooled_gls = se_gls_pooled_rep / se_gls_pooled_act,
       resid_ols = tr$beta - drop(W %*% th_ols))
}
