#!/bin/bash
#SBATCH --job-name=mgep_modsim
#SBATCH -n 1
#SBATCH --mem=4g
#SBATCH -t 01:00:00
#SBATCH --output=longleaf_setup/logs/%x_%j.out
#SBATCH --error=longleaf_setup/logs/%x_%j.err

# multiomicsGEP — Run modular simulation (update_L/F/beta/tau modules)
# Produces: results/history/early-benchmarks/tables/modular_sim/ (7 CSVs) + results/history/early-benchmarks/figures/modular_sim/ (8 PDF+PNG)
# Submit from repo root:  sbatch longleaf_setup/run_modular_simulation.sl

module purge
module load r/4.4.0

export REPO_ROOT=$(pwd)

mkdir -p results/history/early-benchmarks/tables/modular_sim results/history/early-benchmarks/figures/modular_sim

Rscript results/legacy/modular_sim_block/run_modular_simulation.R
