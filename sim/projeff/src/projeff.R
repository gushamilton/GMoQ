# ---------------------------------------------------------------------------
# projeff.R -- exact OLS-vs-GLS projection efficiency for G-MoQ Stage 2,
# across weight-matrix families, plus a Gauss-Markov departure diagnostic.
#
# Algebra (exact; no simulation needed for the asymptotic ratio):
#   P = orthonormal basis of col(W) (T x K), Q = orthonormal complement
#   W = P R,  R = P'W  (K x K invertible)
#   A = P'Sig P,  B = P'Sig Q,  D = Q'Sig Q
#   V_ols = (W'W)^-1 W' Sig W (W'W)^-1 = R^-1 A R^-T
#   V_gls = (W' Sig^-1 W)^-1           = R^-1 (A - B D^-1 B') R^-T
# The GLS gain is EXACTLY a Schur complement: the part of the in-subspace noise
# linearly predictable from the orthogonal-complement noise.
#   R2_j       = 1 - V_gls[j,j]/V_ols[j,j]
#   SEratio_j  = 1/sqrt(1 - R2_j)
#   rho_i      = canonical correlations between P-block and Q-block noise
#   sup over directions of SEratio = 1/sqrt(1 - rho_max^2)
# Gauss-Markov holds iff B = 0 iff all rho_i = 0 iff col(W) is Sigma-invariant.
# ---------------------------------------------------------------------------

TAUS <- seq(0.02, 0.98, length.out = 45)

## ---- skew normal, self-contained (no `sn` dependency) ---------------------
dsn_ <- function(x, xi, omega, alpha)
  2/omega * dnorm((x - xi)/omega) * pnorm(alpha * (x - xi)/omega)

psn_ <- function(x, xi, omega, alpha) {
  z <- (x - xi)/omega
  if (abs(alpha) < 1e-12) return(pnorm(z))
  ng <- 4001L
  u  <- seq(0, 1, length.out = ng) * alpha
  du <- u[2] - u[1]
  Tz <- vapply(z, function(zz) {
    g <- exp(-0.5 * zz^2 * (1 + u^2))/(1 + u^2)
    sum((g[-1] + g[-ng])/2) * du/(2 * pi)
  }, numeric(1))
  pnorm(z) - 2 * Tz
}

qfun_ <- function(pf, taus, lo, hi)
  vapply(taus, function(tt)
    uniroot(function(x) pf(x) - tt, lower = lo, upper = hi,
            tol = .Machine$double.eps^0.75)$root, numeric(1))

## ---- weight-matrix builders (mirror stage2/R/weight-functions.R) ----------
## All are W_{t,j} = -(dF/dtheta_j)(q_t)/f_den(q_t).
## f_den defaults to the model pdf; real testosterone uses the Stage-1 KDE.

W_norm_mix <- function(q, p1, mu1, sd1, mu2, sd2, membership, f_den = NULL) {
  F1 <- pnorm(q, mu1, sd1); f1 <- dnorm(q, mu1, sd1)
  F2 <- pnorm(q, mu2, sd2); f2 <- dnorm(q, mu2, sd2)
  f  <- if (is.null(f_den)) pmax(p1*f1 + (1-p1)*f2, 1e-300) else f_den
  W1 <- p1*f1/f; W2 <- (1-p1)*f2/f
  if (membership) cbind(gamma = -(p1*(1-p1)/f)*(F1-F2), beta_1 = W1, beta_2 = W2)
  else            cbind(beta_1 = W1, beta_2 = W2)
}

W_skew_mix <- function(q, p1, xi1, om1, al1, xi2, om2, al2, membership, f_den = NULL) {
  F1 <- psn_(q, xi1, om1, al1); f1 <- dsn_(q, xi1, om1, al1)
  F2 <- psn_(q, xi2, om2, al2); f2 <- dsn_(q, xi2, om2, al2)
  f  <- if (is.null(f_den)) pmax(p1*f1 + (1-p1)*f2, 1e-300) else f_den
  W1 <- p1*f1/f; W2 <- (1-p1)*f2/f
  if (membership) cbind(gamma = -(p1*(1-p1)/f)*(F1-F2), beta_left = W1, beta_right = W2)
  else            cbind(beta_left = W1, beta_right = W2)
}

W_vqtl <- function(taus, sdv)
  cbind(beta_mu = rep(1, length(taus)), beta_sigma2 = qnorm(taus)/(2*sdv))

W_mix_vqtl <- function(q, p1, mu1, sd1, mu2, sd2) {
  F1 <- pnorm(q, mu1, sd1); f1 <- dnorm(q, mu1, sd1)
  F2 <- pnorm(q, mu2, sd2); f2 <- dnorm(q, mu2, sd2)
  f  <- pmax(p1*f1 + (1-p1)*f2, 1e-300)
  cbind(gamma       = -(p1*(1-p1)/f)*(F1-F2),
        beta_mu1    = p1*f1/f,
        beta_mu2    = (1-p1)*f2/f,
        beta_sigma1 = (p1*f1*(q-mu1)/sd1)/f,
        beta_sigma2 = ((1-p1)*f2*(q-mu2)/sd2)/f)
}

## ---- Brownian-bridge Sigma (asymptotic RIF tau-covariance) ---------------
## Sigma_st = (min - product)/(f_s f_t), up to 1/(N Var(g)) which cancels.
bridge_Sigma <- function(taus, f) (outer(taus, taus, pmin) - outer(taus, taus))/outer(f, f)

## ---- the diagnostic ------------------------------------------------------
proj_diag <- function(W, Sig, label = "") {
  Tt <- nrow(W); K <- ncol(W)
  qrW <- qr(W); Qf <- qr.Q(qrW, complete = TRUE)
  P  <- Qf[, seq_len(K), drop = FALSE]
  Qc <- Qf[, (K+1):Tt, drop = FALSE]
  R  <- crossprod(P, W); Rinv <- solve(R)
  A  <- crossprod(P, Sig %*% P)
  Bm <- crossprod(P, Sig %*% Qc)
  D  <- crossprod(Qc, Sig %*% Qc)
  Schur <- A - Bm %*% solve(D, t(Bm))

  Ainv  <- solve(crossprod(W), t(W))
  V_ols <- Ainv %*% Sig %*% t(Ainv)
  V_gls <- solve(crossprod(W, solve(Sig, W)))
  V_ols_s <- Rinv %*% A %*% t(Rinv)
  V_gls_s <- Rinv %*% Schur %*% t(Rinv)
  chk <- max(abs(V_ols_s - V_ols)/pmax(abs(V_ols), 1e-300),
             abs(V_gls_s - V_gls)/pmax(abs(V_gls), 1e-300))

  # canonical correlations between in-subspace and complement noise
  eA <- eigen(A, symmetric = TRUE)
  Amh <- eA$vectors %*% diag(1/sqrt(pmax(eA$values, 0)), K) %*% t(eA$vectors)
  G <- Amh %*% (A - Schur) %*% Amh
  rho2 <- pmin(pmax(sort(eigen(G, symmetric = TRUE)$values, decreasing = TRUE), 0), 1-1e-15)

  lam <- eigen(Sig, symmetric = TRUE)$values
  kant <- (max(lam) + min(lam))^2 / (4 * max(lam) * min(lam))   # variance-scale bound

  se_ratio <- sqrt(diag(V_ols)/diag(V_gls))
  R2 <- 1 - diag(V_gls)/diag(V_ols)

  # spectral profile: which bridge modes carry the OLS variance
  eS <- eigen(Sig, symmetric = TRUE)
  keff <- vapply(seq_len(K), function(j) {
    a <- Ainv[j, ]
    contrib <- eS$values * as.numeric(crossprod(eS$vectors, a))^2
    sum(seq_along(contrib) * contrib)/sum(contrib)
  }, numeric(1))

  data.frame(model = label, param = colnames(W), K = K,
             se_ratio = se_ratio, R2 = R2,
             rho_max = sqrt(rho2[1]),
             bound_rho = 1/sqrt(1 - rho2[1]),
             bound_kant = sqrt(kant),
             kappa_Sig = max(lam)/min(lam),
             k_eff = keff,
             algebra_check = chk,
             stringsAsFactors = FALSE)
}
