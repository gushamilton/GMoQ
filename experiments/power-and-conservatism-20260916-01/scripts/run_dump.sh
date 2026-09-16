#!/usr/bin/env bash
#SBATCH --job-name=pc-psi
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/user/work/fh6520/gmoq_power_conserv/logs/psi_%j.log
set -euo pipefail
S=/user/work/fh6520/gmoq_power_conserv
DATA=/user/work/fh6520/fungwas_paper

source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
T=${SLURM_CPUS_PER_TASK:-8}
export PYTHONPATH=/user/work/fh6520/gmoq_locovar/pkg:/user/work/fh6520/gmoq_finish_20260915/inputs/pydeps
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T
TAUS=$(python -c 'import numpy as np; print(",".join("%.17g" % t for t in np.linspace(0.02,0.98,45)))')

BGEN=/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr22.bgen
SAMPLE=/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample
PHENO=$DATA/data/derived/paper/phenotypes_paper.txt
COVAR=$DATA/data/derived/paper/covars_no_sex_paper.txt
KEEP=$DATA/data/derived/paper/keep_files/keep_testosterone_series.txt
PRED=$DATA/results/paper/regenie_step1/full/testosterone/testosterone_step1_pred.list

mkdir -p $S/psi $S/tmp
head -200 $DATA/inputs/data/w_hm3_chr22_only.txt > $S/tmp/snps_small.txt

export PSI_OUT=$S/psi
python $S/scripts/dump_psi.py \
  --bgen "$BGEN" --sample "$SAMPLE" \
  --pheno "$PHENO" --pheno-col testosterone \
  --covar "$COVAR" --keep "$KEEP" --snps "$S/tmp/snps_small.txt" \
  --loco-preds "$PRED" --loco-scale regression --loco-mode covariate \
  --taus "$TAUS" --blocks 1000 --batch-size 100 --cov-dtype float32 \
  --seed 42 --threads $T \
  --out "$S/tmp/probe_chr22"

ls -la $S/psi
cat $S/psi/meta_chr22.txt
echo PSIDONE
