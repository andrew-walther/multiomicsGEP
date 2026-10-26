#!/usr/bin/env bash
#SBATCH --job-name=round2-collapse
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=round2-collapse-%A_%a.log
# Run from the repository root after installing the project's R dependencies.
# Each array cell is one fixed training-only diagnosis arm, with no nested
# workers. 8GB covers dense 144 x 13000 data, moments and temporary projections.
# TCGA training cache must already exist at data/round2_tcga_training.rds;
# transfer it securely, never commit controlled data. No jobs are submitted here.
# Do not edit the R runner or fitter while an array is running.
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
Rscript results/multimodal_real/run_round2_collapse_diagnosis.R "${SLURM_ARRAY_TASK_ID}"
