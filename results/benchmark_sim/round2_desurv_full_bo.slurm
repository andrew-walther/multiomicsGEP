#!/usr/bin/env bash
#SBATCH --job-name=round2-desurv-bo
#SBATCH --cpus-per-task=30
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=round2-desurv-bo-%j.log
# Run from repository root. Set PDAC_DATA_ROOT to the controlled data folder
# containing original/ and the unchanged author load_data_internal.R.
# Install pinned DeSurv plus DiceKriging/lhs for the platform before submission.
# Copy unchanged paper R/bo_helpers.R to ignored data/round2_desurv_paper_bo_helpers.R
# (MD5 is checked by the runner). Outputs include the paper03 one-SE rank assessment.
# This is one adaptive BO search, not independent hyperparameter array cells:
# 50 initial+100 BO points,30 starts,5 folds, no external selection.
# Author-reported4-8h on30 cores; request48h to expose, not truncate, long fits.
# Memory: budget48GB for concurrent grid workers and raw gene matrices.
# The independent100-seed final-fit array is round2_desurv_seeds.slurm.
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
test -f data/round2_desurv_paper_bo_helpers.R
desurv_helpers_md5=$(md5sum data/round2_desurv_paper_bo_helpers.R | cut -d ' ' -f 1)
if [[ "$desurv_helpers_md5" != "6f02c072e5b64c5d7582ff76773bee8b" ]]; then
  echo "Pinned author BO helper is missing or changed" >&2
  exit 1
fi
Rscript results/benchmark_sim/run_desurv_h2h.R bo "${SLURM_CPUS_PER_TASK}"
