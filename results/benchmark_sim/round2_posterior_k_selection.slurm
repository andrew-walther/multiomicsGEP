#!/usr/bin/env bash
#SBATCH --job-name=round2-posterior-rank
#SBATCH --array=1-160
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=round2-posterior-rank-%A_%a.log
# Run from repository root with project R dependencies. One serial protected
# single-modality cell per task; same fixed data/criteria as local attempt5.
# 4GB allows R packages, posterior moments and recursive refits. No submission.
# After completion compare all160 cells with the unchanged Task1-3 reference
# rows using k_panel_summary(); never assess adoption from a partial array.
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
Rscript - <<'RS'
source("results/benchmark_sim/run_k_selection_panel.R")
grid <- k_panel_grid()
ids <- which(grid$model=="single" & grid$arm=="protected")
j <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
stopifnot(length(j)==1L,!is.na(j),j>=1L,j<=length(ids))
i <- ids[j]
cfg <- yaml::read_yaml("config/globals.yml")
out_dir <- "results/benchmark_sim/outputs/k_selection/round2/attempt5_posterior_longleaf"
p <- file.path(out_dir,"cells",sprintf("cell_%04d.csv",i))
if(file.exists(p))stop("Existing checkpoint; do not overwrite: ",p)
dir.create(file.path(out_dir,"cells"),recursive=TRUE,showWarnings=FALSE)
load_panel_code()
print(panel_cell(i,posterior_start=TRUE))
RS
