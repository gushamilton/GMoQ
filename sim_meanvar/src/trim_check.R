# check the default trim=0.01 does not materially attenuate r_av
.libPaths(c("/user/work/fh6520/meanvar/rlib","/user/work/fh6520/r_packages"))
suppressMessages(library(fungwasStage2)); set.seed(7)
M <- 50000; r_true <- 0.135
alpha <- rnorm(M,0,0.02); av <- r_true*alpha + rnorm(M,0,0.004)
sa <- runif(M,0.010,0.020); sv <- runif(M,0.004,0.008); rho <- 0.3; cl <- rho*sa*sv
e1 <- rnorm(M); e2 <- rho*e1+sqrt(1-rho^2)*rnorm(M)
a <- alpha+sa*e1; v <- av+sv*e2; blk <- rep(1:22,length.out=M)
for (tr in c(0, 0.001, 0.01, 0.05)) {
  f <- estimate_rav(a, v, sa, sv, cl, block=blk, trim=tr)
  cat(sprintf("trim %.3f  r_av %.5f  se %.5f  n_used %6d  denom %.4f  bias/se %+.2f\n",
      tr, f$r_av, f$se, f$n_used, f$denominator, (f$r_av-r_true)/f$se))
}
cat("\n# with genuine vQTLs present and INDEPENDENT of alpha (Young's assumption)\n")
v2 <- r_true*alpha + rnorm(M,0,0.004) + ifelse(runif(M)<0.002, rnorm(M,0,0.05), 0)
v2h <- v2 + sv*e2
f <- estimate_rav(a, v2h, sa, sv, cl, block=blk, trim=0.01)
cat(sprintf("  r_av %.5f (se %.5f), truth %.3f -> bias/se %+.2f\n", f$r_av, f$se, r_true, (f$r_av-r_true)/f$se))
cat("\n# with genuine vQTLs CORRELATED with alpha (assumption violated)\n")
v3 <- r_true*alpha + rnorm(M,0,0.004) + 0.10*alpha
v3h <- v3 + sv*e2
f <- estimate_rav(a, v3h, sa, sv, cl, block=blk, trim=0.01)
cat(sprintf("  r_av %.5f vs true coupling %.3f -- absorbs the real signal, as expected\n", f$r_av, r_true+0.10))
