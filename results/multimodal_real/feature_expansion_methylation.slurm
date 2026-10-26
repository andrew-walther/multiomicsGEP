#!/bin/bash
#SBATCH --job-name=methylation-count
#SBATCH --array=1-30
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=methylation-count-%A_%a.log
# Prepare the frozen training universe and design before submission.
# Run from the project root with all local ignored RDS inputs available.
# Each array task owns one terminal cell; do not duplicate a live worker.
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
Rscript results/multimodal_real/run_feature_expansion_multimodal.R "${SLURM_ARRAY_TASK_ID}"
