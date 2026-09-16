source("/user/work/fh6520/gmoq_projeff/projeff.R")
options(width = 210, digits = 6)

## Mechanism: is the gamma column the missing Sigma-eigendirection?
## For each model we report
##   - principal angles between col(W) and the span of Sigma's top-m eigenvectors
##   - absolute V (so 2col and 3col models can be compared on the same footing)
##   - pooled-scale invariance of GLS

ang <- function(W, V) {   # principal angles (deg) between col(W) and col(V)
  qa <- qr.Q(qr(W)); qb <- qr.Q(qr(V))
  s <- svd(crossprod(qa, qb))$d
  round(acos(pmin(pmax(s, -1), 1)) * 180/pi, 3)
}

cell <- function(label, q, f, taus, Wlist) {
  Sg <- bridge_Sigma(taus, f)
  eS <- eigen(Sg, symmetric = TRUE)
  cat("\n##", label, "\n")
  cat("  Sigma top-5 eigenvalue share:",
      paste(round(eS$values[1:5]/sum(eS$values), 4), collapse = " "), "\n")
  for (nm in names(Wlist)) {
    W <- Wlist[[nm]]
    Ainv <- solve(crossprod(W), t(W))
    Vo <- Ainv %*% Sg %*% t(Ainv)
    Vg <- solve(crossprod(W, solve(Sg, W)))
    cat(sprintf("  %-14s K=%d | angle(colW, v1)=%7.3f deg | angles to top-K eigspace: %s\n",
                nm, ncol(W), ang(W, eS$vectors[, 1, drop = FALSE])[1],
                paste(ang(W, eS$vectors[, seq_len(ncol(W)), drop = FALSE]), collapse = ",")))
    cat(sprintf("     SEols %s\n     SEgls %s\n     ratio %s\n",
                paste(signif(sqrt(diag(Vo)), 5), collapse = " "),
                paste(signif(sqrt(diag(Vg)), 5), collapse = " "),
                paste(round(sqrt(diag(Vo)/diag(Vg)), 4), collapse = " ")))
  }
}

## ---- normal mixture, matched separations --------------------------------
for (S in c(2, 4, 6)) {
  p1 <- .5; pf <- function(x) p1*pnorm(x,0,1) + (1-p1)*pnorm(x,S,1)
  q <- qfun_(pf, TAUS, -40, S+40)
  f <- p1*dnorm(q,0,1) + (1-p1)*dnorm(q,S,1)
  cell(sprintf("normal mixture S=%g (model pdf denominator)", S), q, f, TAUS,
       list(`2col means` = W_norm_mix(q,p1,0,1,S,1,FALSE),
            `3col +gamma` = W_norm_mix(q,p1,0,1,S,1,TRUE),
            `1col gammaonly` = W_norm_mix(q,p1,0,1,S,1,TRUE)[,1,drop=FALSE]))
}

## ---- real testosterone W ------------------------------------------------
z <- readRDS("/user/work/fh6520/gmoq_noloco/grid/testosterone/weights/weights_skew.rds")
cell("REAL testosterone skew W, Sigma from Stage-1 KDE", z$q_tau, z$density_stage1, z$taus,
     list(`3col shipped` = z$W, `2col nogamma` = z$W[, c("beta_left","beta_right")],
          `1col gammaonly` = z$W[, "gamma", drop = FALSE]))

## ---- GLS invariance to the pooled scale factor --------------------------
cat("\n## pooled-covariance scale invariance\n")
p1 <- .5; S <- 4; pf <- function(x) p1*pnorm(x,0,1) + (1-p1)*pnorm(x,S,1)
q <- qfun_(pf, TAUS, -40, S+40); f <- p1*dnorm(q,0,1) + (1-p1)*dnorm(q,S,1)
Sg <- bridge_Sigma(TAUS, f); W <- W_norm_mix(q,p1,0,1,S,1,FALSE)
g1 <- solve(crossprod(W, solve(Sg, W)), crossprod(W, solve(Sg)))
g2 <- solve(crossprod(W, solve(7.3*Sg, W)), crossprod(W, solve(7.3*Sg)))
cat("  max |A_gls(Sigma) - A_gls(c*Sigma)| =", max(abs(g1 - g2)), "\n")
cat("  -> the GLS projection matrix is scale-free, so a pooled Sigma_0 shared\n")
cat("     across SNPs delivers exactly the same reweighting for every SNP.\n")

## ---- mixture-family prediction: shipped 3col skew W, model pdf ----------
## Bimodal antibody phenotypes (jcv/mcv) use the SAME build_skew_mixture()
## 3-column W, with the MODEL pdf in the denominator. Sweep plausible shapes.
cat("\n## PREDICTION for the mixture family (3col skew W, model pdf denominator)\n")
cat(sprintf("%6s %6s %6s %7s %9s %9s %9s %9s\n",
            "S","alpha","p1","rho_max","gamma","beta_L","beta_R","worstbound"))
pred <- NULL
for (S in c(2,3,4,5,6)) for (al in c(0,-1,2)) for (p1 in c(0.3,0.5,0.7)) {
  om1 <- 1; om2 <- 1
  d <- al/sqrt(1+al^2); msn <- function(om) om*d*sqrt(2/pi); sdsn <- om1*sqrt(1-2*d^2/pi)
  xi1 <- 0; xi2 <- S*sdsn
  pf <- function(x) p1*psn_(x,xi1,om1,al) + (1-p1)*psn_(x,xi2,om2,al)
  q <- qfun_(pf, TAUS, xi1-50, xi2+50)
  f <- p1*dsn_(q,xi1,om1,al) + (1-p1)*dsn_(q,xi2,om2,al)
  if (any(!is.finite(f)) || min(f) <= 0) next
  Sg <- bridge_Sigma(TAUS, f)
  W3 <- W_skew_mix(q,p1,xi1,om1,al,xi2,om2,al,TRUE)
  W2 <- W_skew_mix(q,p1,xi1,om1,al,xi2,om2,al,FALSE)
  d3 <- proj_diag(W3, Sg, "x"); d2 <- proj_diag(W2, Sg, "x")
  cat(sprintf("%6.1f %6.1f %6.1f %7.4f %9.4f %9.4f %9.4f %9.4f   [2col: %.3f %.3f]\n",
              S, al, p1, d3$rho_max[1], d3$se_ratio[1], d3$se_ratio[2], d3$se_ratio[3],
              d3$bound_rho[1], d2$se_ratio[1], d2$se_ratio[2]))
  pred <- rbind(pred, cbind(S=S, alpha=al, p1=p1,
                            r3_gamma=d3$se_ratio[1], r3_bL=d3$se_ratio[2], r3_bR=d3$se_ratio[3],
                            rho3=d3$rho_max[1], r2_bL=d2$se_ratio[1], r2_bR=d2$se_ratio[2],
                            rho2=d2$rho_max[1]))
}
saveRDS(pred, "/user/work/fh6520/gmoq_projeff/out/pred_mixture.rds")
