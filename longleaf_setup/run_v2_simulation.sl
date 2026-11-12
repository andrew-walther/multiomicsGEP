#!/bin/bash
#SBATCH --job-name=mgep_v2sim
#SBATCH -n 1
#SBATCH --mem=4g
#SBATCH -t 01:00:00
#SBATCH --output=longleaf_setup/logs/%x_%j.out
#SBATCH --error=longleaf_setup/logs/%x_%j.err

# multiomicsGEP — Run V2 simulation (code/Supervised_Bayesian_MF_V2.R)
# Produces: results/history/early-benchmarks/tables/full_sim/ (7 CSVs) + results/history/early-benchmarks/figures/full_sim/ (8 PDF+PNG)
# Submit from repo root:  sbatch longleaf_setup/run_v2_simulation.sl

module purge
module load r/4.4.0

export REPO_ROOT=$(pwd)

mkdir -p results/history/early-benchmarks/tables/full_sim results/history/early-benchmarks/figures/full_sim

Rscript results/legacy/full_sim/run_simulation.R
