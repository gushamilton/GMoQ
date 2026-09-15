#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=cssetup
#SBATCH --time=0:40:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH -o /user/work/fh6520/covsource_run/logs/setup.%j.out
set -x
WS=/user/work/fh6520/covsource_run
module load languages/R/4.5.1
R_LIBS=$WS/lib R CMD INSTALL --no-docs --no-html --library=$WS/lib $WS/pkgsrc/stage2 2>&1 | tail -25
echo "=== install done ==="
Rscript --vanilla -e ".libPaths(c(\"$WS/lib\")); print(.libPaths()); library(fungwasStage2); cat(\"LOADED OK\\n\"); print(packageVersion(\"fungwasStage2\")); ns<-asNamespace(\"fungwasStage2\"); cat(\"pooled_cov_data:\", exists(\".pooled_cov_data\", envir=ns), \" param_gwas:\", exists(\"param_gwas\", envir=ns), \"\\n\")"
echo "=== PYTHON TIMING PILOT N=461320 ==="
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd $WS/pkgsrc/stage1
python -c "
import numpy as np, sys, time, resource
sys.path.insert(0,\".\")
from fungwas_stage1 import core, pooled_cov
N=461320; T=45; taus=np.linspace(0.02,0.98,45)
rng=np.random.default_rng(7)
for maf in (0.30,0.005):
    g=rng.binomial(2,maf,N).astype(float); y=rng.normal(size=N)
    t0=time.time(); r=core.run_stage1(g.reshape(-1,1), y, taus, n_blocks=1000); t1=time.time()
    rif,q,f=core.compute_rif(y,taus)
    X=np.ones((N,1)); Qb,_=np.linalg.qr(X); rr=rif-Qb@(Qb.T@rif)
    S0=pooled_cov.compute_sigma0(rr); t2=time.time()
    cv=r[\"covariances\"][0]
    print(\"maf\",maf,\"run_stage1 %.2fs\"%(t1-t0),\"sigma0 %.2fs\"%(t2-t1),\"d_resid\",r[\"d_resid\"][0],\"cov ok\",cv.shape, np.isfinite(cv).all())
    print(\"   jk_se[0:3]\",np.sqrt(np.diag(cv))[:3],\" pooled_se[0:3]\",np.sqrt(np.diag(S0)/r[\"d_resid\"][0])[:3])
print(\"maxRSS GB\", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6)
"
