## Exact population truth for the GLS/OLS coverage measurement.
##
## Y = a*g + (1 + c*g)*eps ,  g ~ Binom(2, maf) ,  eps standardised (mean 0, var 1)
##
## Stage 1s estimand (RIF / unconditional quantile partial effect):
##   beta_tau = Cov(g, RIF_tau(Y)) / Var(g)
##            = -[ E(g 1{Y<=q_tau}) - E(g) tau ] / ( f(q_tau) Var(g) )
## q_tau, f are the MARGINAL population quantile and density.  Exact to the
## root-finder tolerance -- no Monte Carlo anywhere in this file.
##
## Population pooled covariance (Psi residualised on the intercept only):
##   Psi_t = -(1{Y<=q_t} - tau_t)/f_t  =>  Sigma_0[s,t] = (min(t_s,t_t) - t_s t_t)/(f_s f_t)

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
  Fy <- function(y) sum(pg * eps_cdf((y - a * gv) / sg, shape))
  fy <- function(y) sum(pg * eps_pdf((y - a * gv) / sg, shape) / sg)
  lo <- min(a * gv + sg * eps_quantile(1e-9, shape)) - 1
  hi <- max(a * gv + sg * eps_quantile(1 - 1e-9, shape)) + 1
  q <- vapply(taus, function(tt)
    uniroot(function(y) Fy(y) - tt, lower = lo, upper = hi, tol = 1e-12)$root, numeric(1))
  f <- vapply(q, fy, numeric(1))
  beta <- vapply(seq_along(taus), function(i) {
    Eg1 <- sum(pg * gv * eps_cdf((q[i] - a * gv) / sg, shape))
    -(Eg1 - Eg * taus[i]) / (f[i] * Vg)
  }, numeric(1))
  EY  <- a * Eg
  EY2 <- sum(pg * (a^2 * gv^2 + sg^2))
  VY  <- EY2 - EY^2
  S0 <- (outer(taus, taus, pmin) - outer(taus, taus)) / outer(f, f)
  list(taus = taus, q_tau = q, f = f, beta = beta, Sigma0 = S0,
       EY = EY, sdY = sqrt(VY), Eg = Eg, Vg = Vg,
       maf = maf, shape = shape, a = a, c = cc)
}
