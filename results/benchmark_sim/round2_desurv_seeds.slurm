#!/usr/bin/env bash
#SBATCH --job-name=round2-desurv
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=round2-desurv-%A_%a.log
# From the repo root, first run run_desurv_h2h.R prepare with the controlled
# PDAC data and pinned DeSurv package installed for this platform.
# Then submit this array yourself. Each seed uses production4000/4000 limits;
# it writes an independent fit and status, never overwriting another seed.
# After every seed status exists, run run_desurv_h2h.R finalize to perform the
# unchanged consensus fit and external report. No seed is selected by external C.
# Memory: <4GB expected for ~273 patients and ~2000 genes; no nested workers.
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
Rscript -e '
if(dir.exists("data/round2_R_library")) .libPaths(c(normalizePath("data/round2_R_library"),.libPaths()))
stopifnot(packageDescription("DeSurv")$RemoteSha=="cc9d4c2486b3125704ab71eea61f35e040a40080")
source("results/benchmark_sim/run_desurv_h2h.R")
desurv_seed_fit(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")),
 readRDS("results/benchmark_sim/outputs/desurv_h2h/round2/training.rds"),
 "results/benchmark_sim/outputs/desurv_h2h/round2/seeds")'
