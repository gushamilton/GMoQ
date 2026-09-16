#!/bin/bash
#SBATCH --job-name=rg21_female
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --output=/user/work/fh6520/gmoq_power_conserv/logs/rg21_female.out
#SBATCH --error=/user/work/fh6520/gmoq_power_conserv/logs/rg21_female.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
set -euo pipefail
source /etc/profile.d/modules.sh
module purge
module load apps/regenie/4.1
cd /user/work/fh6520/fungwas_paper
regenie --step 2 --bgen /bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr21.bgen --sample /bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample   --keep /user/work/fh6520/gmoq_fair_comparator/keep_step2_female.txt   --phenoFile data/processed/phenotypes_final.txt   --phenoCol testosterone   --covarFile data/processed/covars_no_sex.txt   --extract inputs/data/w_hm3_rsids_only.txt   --bsize 1000 --threads 16   --pred /user/work/fh6520/gmoq_fair_comparator/step1/female/testosterone_female_step1_pred.list   --out /user/work/fh6520/gmoq_power_conserv/regenie21/testosterone_female_chr21
echo DONE_loco_female
