#!/usr/bin/env bash
#SBATCH --job-name=round2-profiled-cox
#SBATCH --array=1-82,85-86,89-90,93-94,97-98,101-102
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=round2-profiled-cox-%A_%a.log
# Run from repo root after the fixed-prior R1 aggregate is complete:
# Rscript results/benchmark_sim/run_round2_profiled_cox.R prepare
# This freezes10 unchanged framework CV references and creates102-cell design.
# Andrew may submit these92 NEW fits; reference cells are deliberately omitted.
# After all finish: Rscript results/benchmark_sim/run_round2_profiled_cox.R finalize
# No job submission here. Transfer only approved TCGA training cache securely,
# never commit it.8GB per serial fit covers dense TCGA moments and temporaries.
# No nested parallelism. Never edit a running runner or fitter.
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
Rscript results/benchmark_sim/run_round2_profiled_cox.R cell "${SLURM_ARRAY_TASK_ID}"
