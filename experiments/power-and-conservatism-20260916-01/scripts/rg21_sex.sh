#!/bin/bash
#SBATCH --job-name=rg21_sex
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --output=/user/work/fh6520/gmoq_power_conserv/logs/rg21_sex.out
#SBATCH --error=/user/work/fh6520/gmoq_power_conserv/logs/rg21_sex.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
set -euo pipefail
source /etc/profile.d/modules.sh
module purge
module load apps/regenie/4.1
cd /user/work/fh6520/fungwas_paper
regenie --step 2 --bgen /bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr21.bgen --sample /bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample   --keep /user/work/fh6520/gmoq_fair_comparator/keep_step2_both.txt   --phenoFile data/processed/phenotypes_final.txt --phenoCol sex   --covarFile data/processed/covars_no_sex.txt   --extract inputs/data/w_hm3_rsids_only.txt   --bsize 1000 --bt --threads 16   --pred /user/work/fh6520/gmoq_fair_comparator/step1/sex/sex_step1_pred.list   --out /user/work/fh6520/gmoq_power_conserv/regenie21/sex_chr21
echo DONE_loco_sex
