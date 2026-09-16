source("/user/work/fh6520/gmoq_projeff/projeff.R")
options(width = 200)
ALL <- NULL
add <- function(x) ALL <<- rbind(ALL, x)

## =========================================================================
## ENDPOINT B: the REAL testosterone W (45x3 skew family), exactly as shipped
## =========================================================================
z <- readRDS("/user/work/fh6520/gmoq_noloco/grid/testosterone/weights/weights_skew.rds")
Wreal <- z$W; taus_r <- z$taus; q_r <- z$q_tau
fk <- z$density_stage1        # Stage-1 KDE: the density that scaled the RIF
fm <- z$density_model         # fitted skew-normal mixture pdf
cat("== real testosterone W ==\n")
print(z$components)
add(proj_diag(Wreal, bridge_Sigma(taus_r, fk), "REAL testo skew W (Sig from KDE)"))
add(proj_diag(Wreal, bridge_Sigma(taus_r, fm), "REAL testo skew W (Sig from model pdf)"))
# real W restricted to the two mean columns (drop membership), to isolate K
add(proj_diag(Wreal[, c("beta_left","beta_right")], bridge_Sigma(taus_r, fk),
              "REAL testo skew W, membership DROPPED"))

## =========================================================================
## ENDPOINT A: the simulation's normal-mixture W, 2 columns, no membership
## =========================================================================
for (S in c(1, 2, 3, 4, 5, 6)) {
  p1 <- 0.5; mu1 <- 0; mu2 <- S; s1 <- s2 <- 1
  pf <- function(x) p1*pnorm(x, mu1, s1) + (1-p1)*pnorm(x, mu2, s2)
  q  <- qfun_(pf, TAUS, mu1 - 40*s1, mu2 + 40*s2)
  f  <- p1*dnorm(q, mu1, s1) + (1-p1)*dnorm(q, mu2, s2)
  Sg <- bridge_Sigma(TAUS, f)
  add(proj_diag(W_norm_mix(q, p1, mu1, s1, mu2, s2, FALSE), Sg,
                sprintf("normal-mix 2col  S=%g", S)))
  add(proj_diag(W_norm_mix(q, p1, mu1, s1, mu2, s2, TRUE), Sg,
                sprintf("normal-mix 3col  S=%g", S)))
  add(proj_diag(W_mix_vqtl(q, p1, mu1, s1, mu2, s2), Sg,
                sprintf("normal-mix 5col  S=%g", S)))
}

## =========================================================================
## SKEW family at MATCHED separations, using testosterone's fitted shapes
## =========================================================================
cp <- z$components
al1 <- cp$alpha[1]; al2 <- cp$alpha[2]
om1 <- cp$omega[1]; om2 <- cp$omega[2]
# skew-normal moments, to set separation on a comparable scale
sn_m  <- function(xi, om, al) { d <- al/sqrt(1+al^2); xi + om*d*sqrt(2/pi) }
sn_sd <- function(om, al)     { d <- al/sqrt(1+al^2); om*sqrt(1 - 2*d^2/pi) }
sdpool <- sqrt((sn_sd(om1,al1)^2 + sn_sd(om2,al2)^2)/2)
for (S in c(1, 2, 3, 4, 5, 6, 6.21)) {
  xi1 <- 0; xi2 <- S*sdpool - (sn_m(0,om2,al2) - sn_m(0,om1,al1))
  pf <- function(x) 0.5*psn_(x, xi1, om1, al1) + 0.5*psn_(x, xi2, om2, al2)
  q  <- qfun_(pf, TAUS, xi1 - 40*om1, xi2 + 40*om2)
  f  <- 0.5*dsn_(q, xi1, om1, al1) + 0.5*dsn_(q, xi2, om2, al2)
  Sg <- bridge_Sigma(TAUS, f)
  add(proj_diag(W_skew_mix(q, .5, xi1,om1,al1, xi2,om2,al2, FALSE), Sg,
                sprintf("skew-mix 2col    S=%.2f", S)))
  add(proj_diag(W_skew_mix(q, .5, xi1,om1,al1, xi2,om2,al2, TRUE), Sg,
                sprintf("skew-mix 3col    S=%.2f", S)))
}

## =========================================================================
## Reference families: plain Normal mean/variance (vqtl)
## =========================================================================
q0 <- qnorm(TAUS); f0 <- dnorm(q0)
add(proj_diag(W_vqtl(TAUS, 1), bridge_Sigma(TAUS, f0), "normal vqtl 2col"))

## =========================================================================
## Does the KDE-vs-model density mismatch explain testosterone? Rebuild the
## skew W at testosterone's own parameters with each density in the denominator
## =========================================================================
xi1 <- cp$xi[1]; xi2 <- cp$xi[2]
add(proj_diag(W_skew_mix(q_r, .5, xi1,om1,al1, xi2,om2,al2, TRUE, f_den = fm),
              bridge_Sigma(taus_r, fm), "testo params, model pdf everywhere"))
add(proj_diag(W_skew_mix(q_r, .5, xi1,om1,al1, xi2,om2,al2, TRUE, f_den = fk),
              bridge_Sigma(taus_r, fk), "testo params, KDE everywhere"))
add(proj_diag(W_skew_mix(q_r, .5, xi1,om1,al1, xi2,om2,al2, FALSE, f_den = fk),
              bridge_Sigma(taus_r, fk), "testo params, KDE, 2col"))

saveRDS(ALL, "/user/work/fh6520/gmoq_projeff/out/diag.rds")
fm2 <- ALL
fm2[, c("se_ratio","R2","rho_max","bound_rho","bound_kant","k_eff")] <-
  round(fm2[, c("se_ratio","R2","rho_max","bound_rho","bound_kant","k_eff")], 4)
fm2$kappa_Sig <- signif(fm2$kappa_Sig, 3)
fm2$algebra_check <- signif(fm2$algebra_check, 2)
print(fm2, row.names = FALSE)
