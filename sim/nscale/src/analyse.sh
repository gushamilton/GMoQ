#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=test
#SBATCH --job-name=nsan
#SBATCH --time=0:55:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH -o /user/work/fh6520/nscale_run/logs/analyse.%j.out
WS=/user/work/fh6520/nscale_run
SRC=/user/work/fh6520/wt_nscale/sim/nscale/src
module load languages/R/4.5.1
Rscript --vanilla $SRC/summarise.R $WS/lib $WS/out
echo "############ CROSSCHECK vs gls_coverage (N = 30,000) ############"
Rscript --vanilla $SRC/crosscheck.R $WS/lib $WS/out
echo "############ REPORT ############"
Rscript --vanilla $SRC/report.R $WS/lib $WS/out
