source("/user/work/fh6520/gmoq_projeff/projeff.R")
options(width = 200)

## If Stage 2 weights by Sw but the truth is St, then for a projection
## A (= (W'W)^-1W' for OLS, or (W'Sw^-1W)^-1W'Sw^-1 for GLS):
##   true Var  = A St A'          reported Var = A Sw A'  (GLS: = (W'Sw^-1W)^-1)
## calibration = sqrt(reported/true). OLS's POINT estimate does not depend on
## Sigma at all; GLS's does, so GLS can be both inefficient AND miscalibrated.
arm <- function(W, Sw, St) {
  Ao <- solve(crossprod(W), t(W))
  Ag <- solve(crossprod(W, solve(Sw, W)), crossprod(W, solve(Sw)))
  list(ols_true = sqrt(diag(Ao %*% St %*% t(Ao))),
       ols_rep  = sqrt(diag(Ao %*% Sw %*% t(Ao))),
       gls_true = sqrt(diag(Ag %*% St %*% t(Ag))),
       gls_rep  = sqrt(diag(Ag %*% Sw %*% t(Ag))))
}
show <- function(lab, W, Sw, St) {
  a <- arm(W, Sw, St)
  for (j in seq_len(ncol(W)))
    cat(sprintf("  %-34s %-11s OLS cal %6.3f eff %6.3f | GLS cal %6.3f eff %6.3f\n",
        if (j == 1) lab else "", colnames(W)[j],
        a$ols_rep[j]/a$ols_true[j], a$ols_true[j]/min(a$gls_true[j], a$ols_true[j]),
        a$gls_rep[j]/a$gls_true[j], a$gls_true[j]/min(a$gls_true[j], a$ols_true[j])))
}

cat("cal = reported SE / true SE (1 = calibrated, <1 = anti-conservative)\n")
cat("eff = true SE / best true SE in that row pair\n")

## ---- (1) real testosterone: KDE vs model density is a REAL 16% Sigma error
z <- readRDS("/user/work/fh6520/gmoq_noloco/grid/testosterone/weights/weights_skew.rds")
Sk <- bridge_Sigma(z$taus, z$density_stage1)   # Sigma implied by the Stage-1 KDE
Sm <- bridge_Sigma(z$taus, z$density_model)    # Sigma implied by the fitted model
cat("\n## REAL testosterone W (3col shipped). Truth = KDE-Sigma, Stage 2 weights by model-Sigma\n")
show("shipped 3col", z$W, Sm, Sk)
cat("\n## same, correctly specified (Sw = St = KDE)\n")
show("shipped 3col", z$W, Sk, Sk)
cat("\n## means-only 2col W, correctly specified -- where GLS IS worth a lot\n")
show("2col nogamma", z$W[, c("beta_left","beta_right")], Sk, Sk)
cat("\n## means-only 2col W, Sigma misspecified by the same 16% density error\n")
show("2col nogamma", z$W[, c("beta_left","beta_right")], Sm, Sk)

## ---- (2) synthetic: graded Sigma perturbation at normal-mixture S=4 and S=6
for (S in c(4, 6)) {
  p1 <- .5; pf <- function(x) p1*pnorm(x,0,1) + (1-p1)*pnorm(x,S,1)
  q <- qfun_(pf, TAUS, -40, S+40); f <- p1*dnorm(q,0,1) + (1-p1)*dnorm(q,S,1)
  St <- bridge_Sigma(TAUS, f)
  W2 <- W_norm_mix(q,p1,0,1,S,1,FALSE); W3 <- W_norm_mix(q,p1,0,1,S,1,TRUE)
  for (eps in c(0, 0.02, 0.05, 0.10)) {
    set.seed(11); fp <- f * exp(eps * sin(3*pi*TAUS))   # smooth density error
    Sw <- bridge_Sigma(TAUS, fp)
    cat(sprintf("\n## normal-mix S=%g, density error %.0f%%\n", S, 100*eps))
    show("2col means", W2, Sw, St)
    show("3col +gamma", W3, Sw, St)
  }
}
