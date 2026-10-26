#!/bin/bash
#SBATCH --job-name=feature-expansion
#SBATCH --array=1-24
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=feature-expansion-%A_%a.log
# Prepare inputs once from the repository root before submission:
# Rscript results/benchmark_sim/run_feature_expansion.R prepare
# Set PDAC_DATA_ROOT for preparation. Fits read the frozen local training RDS.
# Each task is independent; never launch two workers for the same cell.
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
Rscript results/benchmark_sim/run_feature_expansion.R "${SLURM_ARRAY_TASK_ID}"
