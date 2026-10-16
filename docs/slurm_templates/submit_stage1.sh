#!/bin/bash
#SBATCH --job-name=fungwas_stage1
#SBATCH --partition=short
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --array=1-22
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/stage1_chr%a_%j.log

# FungWas Stage 1 SLURM Template
# Runs one chromosome per array task

# Load modules (adjust for your cluster)
module load lang/python/3.10

# Configuration
BGEN_DIR="/path/to/bgen/files"
PHENO_FILE="inputs/regenie_pheno.txt"
PHENO_COL="your_phenotype"   # column in PHENO_FILE (first column must match BGEN sample IDs)
COVAR_FILE="inputs/regenie_covar.txt"
HM3_DIR="inputs"  # Contains chr${CHR}_hm3.txt files
OUT_DIR="results/stage1"

# Get chromosome from array task
CHR=${SLURM_ARRAY_TASK_ID}
CHR_PAD=$(printf "%02d" $CHR)

# Create output directory
mkdir -p ${OUT_DIR}
mkdir -p logs

# Run Stage 1
echo "Running Stage 1 for chromosome ${CHR}"
fungwas-stage1 \
    --bgen "${BGEN_DIR}/data.chr${CHR_PAD}.bgen" \
    --pheno "${PHENO_FILE}" \
    --pheno-col "${PHENO_COL}" \
    --covar "${COVAR_FILE}" \
    --snps "${HM3_DIR}/chr${CHR}_hm3.txt" \
    --batch-size 1000 \
    --threads 4 \
    --output-rtau \
    --out "${OUT_DIR}/chr${CHR}"

echo "Chromosome ${CHR} complete"
