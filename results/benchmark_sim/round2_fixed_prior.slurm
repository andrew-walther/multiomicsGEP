#!/usr/bin/env bash
#SBATCH --job-name=round2-fixed-prior
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=round2-fixed-prior-%A_%a.log
# From repo root: Rscript results/benchmark_sim/run_round2_fixed_prior.R prepare
# Then Andrew may submit this array; after every cell finishes:
# Rscript results/benchmark_sim/run_round2_fixed_prior.R finalize
# No jobs are submitted automatically. Transfer the approved TCGA-only cache
# securely to data/round2_tcga_training.rds; never commit controlled data.
# Each independent fit is serial.8GB covers dense TCGA folds and moments;
# simulation cells are much smaller. Never edit running scripts or fitters.
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
Rscript results/benchmark_sim/run_round2_fixed_prior.R cell "${SLURM_ARRAY_TASK_ID}"
