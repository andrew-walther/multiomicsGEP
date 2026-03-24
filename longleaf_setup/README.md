# Longleaf Execution Guide — multiomicsGEP

Run the Supervised Bayesian Matrix Factorization project on UNC's Longleaf HPC cluster.

**Requirements:** R 4.4.0, two packages (`ebnm`, `survival`).

---

## What's in This Directory

```
longleaf_setup/
├── README.md                  # This file — the guide you're reading
├── install_packages.R         # One-time R package installer (ebnm + survival)
├── run_tests.sl               # SLURM job: run the 105-test suite
├── run_v2_simulation.sl       # SLURM job: run the V2 monolithic simulation
├── run_modular_simulation.sl  # SLURM job: run the modular simulation
├── run_demos.sl               # SLURM job: run all 4 demo scripts
└── logs/                      # SLURM writes job output here (stdout + stderr)
```

### What are `.sl` files?

`.sl` files are **SLURM job scripts** — shell scripts with special `#SBATCH` comment
lines at the top that tell Longleaf how much memory, time, and CPUs your job needs.
When you run `sbatch some_script.sl`, Longleaf reads those directives, puts your job
in a queue, and runs it on a compute node when resources are available. You don't need
to wait at your terminal — the output goes to a log file you can check later.

### What are the two simulations?

The project has two ways to run the same algorithm:

- **V2 simulation** (`run_v2_simulation.sl`): Runs `code/Supervised_Bayesian_MF_V2.R`,
  the single monolithic script that contains everything — the CAVI loop, all update
  equations, helpers, and visualization. This is the original implementation.

- **Modular simulation** (`run_modular_simulation.sl`): Runs
  `results/run_modular_simulation.R`, which builds the CAVI loop from the four
  standalone update modules (`code/update_L.R`, `update_F.R`, `update_beta.R`,
  `update_tau.R`). Same algorithm, same math, but using the modular code that has
  its own test suite.

Both produce the same outputs (7 CSV tables + 8 figure pairs in PDF and PNG).
You can run either or both — they use different random seeds so results will
differ slightly, but the algorithm behavior should be consistent.

### What do the other jobs do?

- **Tests** (`run_tests.sl`): Runs the 105-test suite that validates the four
  modular update scripts. Good to run first to confirm everything works on Longleaf.

- **Demos** (`run_demos.sl`): Runs the four educational demo scripts that illustrate
  the mathematical properties of each update (signal recovery, shrinkage behavior,
  error-in-variables effects, etc.). Output is narrative text — no figures are saved.

---

## How Longleaf Works (Quick Primer)

Longleaf is a shared compute cluster. You **SSH in** to a login node, but you
**don't run heavy computation there** — instead you either:

1. **`srun`**: Request an interactive session on a compute node (like getting a
   temporary terminal on a powerful machine). Good for installing packages and
   quick test runs.

2. **`sbatch`**: Submit a job script to the queue. Longleaf runs it when resources
   are free and saves the output to a log file. Good for anything that takes more
   than a few minutes.

Your files live on a shared filesystem, so both the login node and compute nodes
can see the same `/work/users/...` directory.

---

## Phase 1 — First-Time Setup

Do this once when you first set up the project on Longleaf.

```bash
# 1. SSH into Longleaf from your Mac terminal
ssh <onyen>@longleaf.unc.edu

# 2. Go to your work directory
#    Replace the letters and onyen with yours.
#    Example: if your onyen is "awalther", the path would be
#    /work/users/a/w/awalther
cd /work/users/<first_letter>/<second_letter>/<onyen>

# 3. Clone the repo (one-time)
git clone https://github.com/andrew-walther/multiomicsGEP.git
cd multiomicsGEP

# 4. Create the directories that jobs will write into
#    (git doesn't track empty directories, so we create them manually)
mkdir -p longleaf_setup/logs \
         results/tables/full_sim results/tables/modular_sim \
         results/figures/full_sim results/figures/modular_sim

# 5. Install R packages in an interactive session
#    srun gives you a shell on a compute node (the login node is not
#    meant for running R). The flags request 4 GB RAM for 30 minutes.
srun -p interact -n 1 --mem=4G -t 00:30:00 --pty bash

# Now you're on a compute node. Load R and install packages:
module add r/4.4.0
Rscript longleaf_setup/install_packages.R

# When it finishes, leave the compute node:
exit
```

After this, `ebnm` and `survival` are installed in your personal R library on
Longleaf. You won't need to do this again unless you delete your library or
the R module version changes.

---

## Phase 2 — Test Run (Verify Before Full Submission)

Before submitting batch jobs, confirm that R and the code work interactively.

```bash
# 6. Start an interactive session and run the test suite
srun -p interact --mem=2G --time=00:10:00 --pty bash
module add r/4.4.0
export REPO_ROOT=$(pwd)
Rscript tests/run_tests.R

# You should see:
#   --- tests/test_update_beta.R: 24/24 tests passed ---
#   --- tests/test_update_L.R:    28/28 tests passed ---
#   --- tests/test_update_F.R:    26/26 tests passed ---
#   --- tests/test_update_tau.R:  27/27 tests passed ---
#   FINAL: ... All tests PASSED.

exit
```

If all 105 tests pass, you're good to go. If something fails, check that you
ran steps 4 and 5 from Phase 1.

---

## Phase 3 — Submit Jobs

Now you can submit batch jobs. These run in the background — you can log out
and come back later.

```bash
# 7. Submit whichever jobs you want (from the repo root directory)
sbatch longleaf_setup/run_tests.sl               # ~30 seconds
sbatch longleaf_setup/run_v2_simulation.sl        # ~5 minutes
sbatch longleaf_setup/run_modular_simulation.sl   # ~5 minutes
sbatch longleaf_setup/run_demos.sl                # ~2 minutes

# Each sbatch command prints a job ID like:
#   Submitted batch job 12345678

# 8. Check on your jobs
squeue -u <onyen>              # shows your running/queued jobs
squeue -u <onyen> | wc -l      # quick count (subtract 1 for the header line)
```

**Tip:** You don't need to submit all four. If you just want the modular
simulation, run only `sbatch longleaf_setup/run_modular_simulation.sl`.

---

## Phase 4 — Check Results

When `squeue` shows no more jobs (or you get an email if configured), check
the output.

```bash
# 9. View job logs (stdout from your R scripts)
cat longleaf_setup/logs/mgep_tests_*.out         # test results
cat longleaf_setup/logs/mgep_v2sim_*.out          # V2 simulation console output
cat longleaf_setup/logs/mgep_modsim_*.out         # modular simulation console output
cat longleaf_setup/logs/mgep_demos_*.out          # demo narrative output

# Check for errors (empty = good)
cat longleaf_setup/logs/mgep_v2sim_*.err

# Verify output files were created
ls results/tables/full_sim/       # should see 7 CSVs
ls results/figures/full_sim/      # should see 8 PDF + 8 PNG files
ls results/tables/modular_sim/    # should see 7 CSVs
ls results/figures/modular_sim/   # should see 8 PDF + 8 PNG files
```

### Expected outputs per simulation

Each simulation produces:

| Directory | Contents |
|-----------|----------|
| `results/tables/*/factor_summary_table.csv` | Per-factor beta, log-rank p, sparsity, PVE |
| `results/tables/*/beta_comparison_table.csv` | Estimated vs. true betas with posterior SD |
| `results/tables/*/cindex_comparison.csv` | C-index: supervised loadings vs. PCA |
| `results/tables/*/convergence_history.csv` | RMSE and ELBO per iteration |
| `results/tables/*/top_features_GEP[1-5].csv` | Top 10 features per factor |
| `results/tables/*/loading_correlation_matrix.csv` | True vs. estimated loading correlations |
| `results/tables/*/ph_test_results.csv` | Proportional hazards test |
| `results/figures/*/fig[1-8]_*.pdf` + `.png` | RMSE trace, ELBO, beta comparison, heatmap, KM curves, signal recovery, correlation heatmap, tau distribution |

---

## Phase 5 — Get Results to Your Mac

Run this **on your Mac** (not on Longleaf) to copy the results back:

```bash
# 10. From your Mac terminal:
scp -r <onyen>@longleaf.unc.edu:/work/users/<first>/<second>/<onyen>/multiomicsGEP/results/ \
       ~/GithubProjects/multiomicsGEP/results/
```

This copies the entire `results/` directory (tables, figures, logs) to your
local repo, where you can open the PDFs and CSVs directly.

---

## Quick Reference

### SLURM scripts at a glance

| Script | Job Name | Mem | Time | What it runs |
|--------|----------|-----|------|--------------|
| `run_tests.sl` | mgep_tests | 2g | 10 min | Test suite (105 tests) |
| `run_v2_simulation.sl` | mgep_v2sim | 4g | 1 hr | V2 monolithic simulation |
| `run_modular_simulation.sl` | mgep_modsim | 4g | 1 hr | Modular simulation |
| `run_demos.sl` | mgep_demos | 2g | 30 min | All 4 demo scripts |

All scripts log to `longleaf_setup/logs/` and must be submitted from the repo root.

### Common commands

| What | Command |
|------|---------|
| Submit a job | `sbatch longleaf_setup/run_tests.sl` |
| Check your queue | `squeue -u <onyen>` |
| Cancel a job | `scancel <job_id>` |
| Cancel all your jobs | `scancel -u <onyen>` |
| View job output | `cat longleaf_setup/logs/mgep_tests_*.out` |
| View job errors | `cat longleaf_setup/logs/mgep_tests_*.err` |
| Interactive session | `srun -p interact --mem=2G -t 00:10:00 --pty bash` |
| Check available R | `module avail r` |
| Load R | `module add r/4.4.0` |

### Troubleshooting

- **"Cannot find repo root"** error in R: You submitted the job from the wrong
  directory. `cd` to `multiomicsGEP/` before running `sbatch`.
- **Package not found**: Re-run Phase 1, step 5 (install packages).
- **Job stuck in queue**: Longleaf is busy. Check `squeue -u <onyen>` — status
  `PD` means pending (waiting for resources), `R` means running.
- **Out of memory**: Increase `--mem` in the `.sl` file (e.g., `--mem=8g`).
- **Out of time**: Increase `-t` (e.g., `-t 02:00:00` for 2 hours).
