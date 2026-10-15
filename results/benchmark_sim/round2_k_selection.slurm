#!/usr/bin/env bash
#SBATCH --job-name=round2-rank
#SBATCH --array=1-800
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=round2-rank-%A_%a.log
# Run from repository root with the project's R dependencies. No jobs submitted.
# Each cell uses the unchanged committed 800-cell design and one serial fit.
# The six-worker local panel took over one hour; observed workers used about
# 700MB RSS. Reserve 4GB per array cell for R packages/compilation and refits.
# This fresh output path intentionally differs from every local attempt.
# Do not edit the runner or its fitters while the array is running.
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
Rscript - <<'RS'
source("results/benchmark_sim/run_k_selection_panel.R")
i <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
stopifnot(length(i)==1L,!is.na(i),i>=1L,i<=800L)
grid <- k_panel_grid()
cfg <- yaml::read_yaml("config/globals.yml")
out_dir <- "results/benchmark_sim/outputs/k_selection/round2/attempt3_longleaf"
path <- file.path(out_dir,"cells",sprintf("cell_%04d.csv",i))
if(file.exists(path))stop("Existing checkpoint; do not overwrite: ",path)
dir.create(file.path(out_dir,"cells"),recursive=TRUE,showWarnings=FALSE)
load_panel_code()
print(panel_cell(i))
RS
# Once all 800 checkpoints exist, use k_panel_summary() on their row-bound CSVs
# to assess the unchanged rule; keep errors in the denominator. Never select on
# a partial array. Local attempt3 remains the primary report for Round2 Task1.
