#!/bin/bash
#SBATCH --job-name=mgep_demos
#SBATCH -n 1
#SBATCH --mem=2g
#SBATCH -t 00:30:00
#SBATCH --output=longleaf_setup/logs/%x_%j.out
#SBATCH --error=longleaf_setup/logs/%x_%j.err

# multiomicsGEP — Run all 4 interactive demo scripts
# Submit from repo root:  sbatch longleaf_setup/run_demos.sl

module purge
module load r/4.4.0

export REPO_ROOT=$(pwd)

for demo in demos/demo_update_beta.R demos/demo_update_L.R demos/demo_update_F.R demos/demo_update_tau.R; do
  echo "=== Running $demo ==="
  Rscript "$demo"
  echo ""
done
