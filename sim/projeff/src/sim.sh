#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=projeff
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G
#SBATCH --time=0-08:00:00
#SBATCH --output=/user/work/fh6520/gmoq_projeff/logs/%x_%a.out
#SBATCH --array=1-9

source /etc/profile.d/modules.sh
module load languages/R/4.5.1
export R_PROFILE_USER=/dev/null
export R_LIBS=/user/work/fh6520/gmoq_verify_work/deplib
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

CELLS=(N2 N3 N4 N5 N6 K4 K6 MIX4 MIX6)
C=${CELLS[$((SLURM_ARRAY_TASK_ID-1))]}
cd /user/work/fh6520/gmoq_projeff
echo "cell=$C host=$(hostname) start=$(date)"
Rscript run_sim.R $C 1000 20
echo "end=$(date)"
