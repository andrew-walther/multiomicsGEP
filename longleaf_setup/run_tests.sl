#!/bin/bash
#SBATCH --job-name=mgep_tests
#SBATCH -n 1
#SBATCH --mem=2g
#SBATCH -t 00:10:00
#SBATCH --output=longleaf_setup/logs/%x_%j.out
#SBATCH --error=longleaf_setup/logs/%x_%j.err

# multiomicsGEP — Run test suite (105 tests)
# Submit from repo root:  sbatch longleaf_setup/run_tests.sl

module purge
module load r/4.4.0

export REPO_ROOT=$(pwd)

Rscript tests/run_tests.R
