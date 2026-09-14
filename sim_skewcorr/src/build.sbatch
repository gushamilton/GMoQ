#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=skewbuild
#SBATCH --time=0:50:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o /user/work/fh6520/skewcorr/logs/build.%j.out
WS=/user/work/fh6520/skewcorr
TREE=/user/work/fh6520/wt_gls-default
LIB=$WS/rlib
DOCLIB=/user/work/fh6520/r_packages     # roxygen2 lives here; added EXPLICITLY, never silently
mkdir -p $LIB
module load languages/R/4.5.1

echo "=== 1. roxygenise (regenerates NAMESPACE from fungwasStage2-package.R) ==="
Rscript --vanilla -e "
.libPaths(c('$LIB', '$DOCLIB', .libPaths()))
stopifnot('$LIB' %in% .libPaths())
cat('roxygen2:', as.character(packageVersion('roxygen2')), '\n')
roxygen2::roxygenise('$TREE/stage2')
" 2>&1 | tail -30
echo "=== roxygenise rc ${PIPESTATUS[0]} ==="

echo "=== 2. NAMESPACE guard ==="
grep -n useDynLib $TREE/stage2/NAMESPACE || { echo "!!! NO useDynLib !!!"; exit 9; }
if [ -e $TREE/stage2/NAMESPACE.backup ]; then echo "!!! NAMESPACE.backup PRESENT !!!"; exit 9; fi
grep -n "dispersion_skew_bias\|correct_dispersion_skew" $TREE/stage2/NAMESPACE

echo "=== 3. install into the private lib ==="
Rscript --vanilla -e "
.libPaths('$LIB')
stopifnot(identical(.libPaths()[1], '$LIB'))
install.packages('$TREE/stage2', repos=NULL, type='source', lib='$LIB')
" 2>&1 | tail -25
echo "=== install rc ${PIPESTATUS[0]} ==="

echo "=== 4. load + the identity the correction rests on ==="
Rscript --vanilla -e "
.libPaths('$LIB')
stopifnot(identical(.libPaths()[1], '$LIB'))
library(fungwasStage2)
cat('LOADED OK  version', as.character(packageVersion('fungwasStage2')), '\n')
taus <- seq(0.02, 0.98, length.out=45)
z <- qnorm(taus); sdY <- 1.3
W <- make_weights_vqtl(taus, q_tau = sdY*z, mu = 0, sd = sdY)
S <- outer(taus,taus,pmin) - outer(taus,taus) + diag(1e-8, 45)
delta <- 0.2; k <- 1 - 2*0.01
bias_b <- (delta^2/(2*sdY)) * k * z          # the derived bias in b(tau)
for (nm in c('ols','gls')) {
  A <- if (nm=='ols') solve(crossprod(W), t(W)) else solve(t(W)%*%solve(S)%*%W, t(W)%*%solve(S))
  cat(sprintf('%s: max|AW - I| = %.3e   projected bias = (%.10f, %.10f)  predicted = (0, %.10f)\n',
      nm, max(abs(A%*%W - diag(2))), (A%*%bias_b)[1], (A%*%bias_b)[2], delta^2*k))
}
cat('\n--- correct_dispersion_skew ---\n')
print(correct_dispersion_skew(beta_mu=0.2, beta_sigma2=0.0392, se_mu=0.01,
                              se_sigma2=0.012, cov_mu_sigma2=0, eaf=0.01))
cat('bias helper:', dispersion_skew_bias(0.2, 0.01, eaf=0.01), '\n')
cat('allele-flip covariance (must be equal):',
    dispersion_skew_bias(0.2, 0.01, eaf=0.01), 'vs',
    -dispersion_skew_bias(-0.2, 0.01, eaf=0.99), '\n')
cat('MAF 0.5 is exactly zero:', dispersion_skew_bias(0.2, 0.01, eaf=0.5), '\n')
cat('refuses an unsigned frequency: ')
print(tryCatch(dispersion_skew_bias(0.2, 0.01), error=function(e) 'ERROR (correct)'))
"
echo "=== smoke rc $? ==="

echo "=== 5. package tests ==="
for t in test-namespace-dynlib.R test-dispersion-guards.R test-vQTL.R test-mixture-vQTL.R; do
  echo "--- $t ---"
  Rscript --vanilla -e ".libPaths('$LIB'); source('$TREE/stage2/tests/$t')" 2>&1 | tail -20
  echo "rc=$?"
done
