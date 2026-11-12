# Codex prompt: rerun task B and task C's expression comparator at K = 5

Date: 2026-11-12. Author: Andrew Walther.

Paste into one Codex session. The two parts run **in order in that one session**: B first,
then C once task C's current K = 7 runs have finished. If the Codex session running task C
is still active, this can start now in a separate session, because Part 1 doesn't touch C's
files.

---

Work in `~/GithubProjects/multiomicsGEP` on `main`. Do not create a worktree or change
branches. Read AGENTS.md, then the 2026-11-12 "Task F" entry in DECISIONS.md. Andrew
adopted K = 5 for the working YFB model. New work reads `benchmark$k_working_yfb` (5) from
`config/globals.yml`. **Do not change `benchmark$k_merged_yfb_desurv`**, which stays 7 for
historical runners and for task C's running K = 7 jobs. Then read
`docs/findings/2026-11-12-beta-prior-comparison-findings.md` (task B at K = 7) and
`docs/findings/2026-11-12-training-cv-redo-findings.md` (task F).

Everything here is the working model except where stated: `fit_ssbmf_yfb(canonicalize = TRUE,
w_S = 1)`, point-Laplace L and F, per-cohort z-score input on the 2,218 D4 genes, pooled
TCGA and CPTAC without strata, no pruning, config `max_iter`/`tol`/`alpha`. The primary
metric is within-cohort concordance.

### Part 1: task B at K = 5 (real data only)

1. In `results/benchmark_sim/run_training_cv_redo.R`, add a mode `brerun` next to the
   existing `arerun` mode. Follow `arerun`'s pattern, and don't change any existing mode
   or default.
   - Arms: shared point-Laplace, shared point-normal and shared Normal β
     (`prior_beta` with `beta_prior_mode = "shared"`).
   - Input: `D4_Z`. K from the command line. Folds 0 to 5.
   - Extend `f_fit()` only with optional arguments whose defaults keep current behavior.
   - Record point-normal inclusion probabilities with `lfp_pip()` from
     `run_lf_prior_comparison.R` (already sourced).
   - The Normal (per-program) reference at K = 5 is the existing rank-sweep checkpoint
     `rank_D4_Z_K05_f*.rds`. Do not refit it.
2. Run `Rscript results/benchmark_sim/run_training_cv_redo.R brerun 3 5`, with at most 3
   workers while task C is running and `OMP_NUM_THREADS=1`. That is 18 fits, about 1 to 2
   minutes each.
3. Write `b_rerun_paired.csv` to `results/benchmark_sim/outputs/training_cv_redo_2026-11-12/`.
   Per arm against Normal on the same folds, include:
   - mean within-cohort and pooled C, with the paired difference and its SE
   - folds converged
   - mean survival-active programs
   - folds with every β zero
   - the full-fit inclusion probabilities.
4. Apply part (c) of B's registered rule (CV not lower than Normal by more than one SE) at
   K = 5. Parts (a) and (b) are simulation criteria and are not rerun. State that plainly.

### Part 2: task C's expression comparator at K = 5

1. First check that task C's K = 7 runs have finished: no `run_methylation_single_modality.R`
   process (`ps aux`), and its findings note no longer says "in progress". If they're still
   running, stop here, say so, and resume Part 2 when they finish. Do not run C fits
   alongside them.
2. In `results/multimodal_real/run_methylation_single_modality.R`, add a K option that
   defaults to the current `k_merged_yfb_desurv` (7), so every existing job and checkpoint
   is unchanged. Write K = 5 checkpoints under a separate key or folder so they can't
   overwrite or resume K = 7 checkpoints.
3. At K = 5, fit the expression comparator (D4 genes, the same matched TCGA patients and
   saved folds) under both setups, PL and PE: the full fit and five folds each, 12 fits.
   Then, for each methylation count whose CV completed at K = 7, project the K = 5 CV time
   from the K = 7 timings. Run the K = 5 methylation full fit and CV for those counts only
   if the projected total is under 3 hours. Otherwise record them as deferred with the
   projection.
4. Recompute C's comparisons at K = 5:
   - CV concordance against expression
   - programs with nonzero loadings
   - survival-active programs
   - methylation-only programs (score correlation below 0.5 with every expression program,
     and survival-active in the full fit and in at least 3 of 5 folds).
   Report K = 5 next to the existing K = 7 results without replacing them.

### Writing and rules

- Add a section "Rerun at K = 5" to each findings note: B's note
  (`docs/findings/2026-11-12-beta-prior-comparison-findings.md`) and C's note
  (`docs/findings/2026-11-12-methylation-single-modality-findings.md`). Keep the K = 7
  results and label both ranks clearly.
- Do not edit DECISIONS.md, PROJECT_STATUS.qmd, the week plan or the attempts log. Put
  proposed text for them in a final "Text for the shared status files (K = 5)" section of
  each note.
- Patient-level fits and scores go under the ignored `data/` directory. Commit only
  aggregate CSVs.
- Script conventions per AGENTS.md: roxygen on every new or changed function, author
  Andrew Walther, no tool names.
- Run `Rscript tests/run_tests.R` once at the end if you changed anything under `code/`,
  or touched a runner that tests cover, and report the count.
- Stage only your own paths (`git add <paths>`, never `git add -A`). Commit each part
  separately with a detailed message and no tool attribution, then push.

**Done when:** Part 1's 18 fits and paired CSV are committed with B's K = 5 section, and
Part 2 is either complete (12 expression fits, any qualifying methylation fits, C's K = 5
section) or explicitly paused because C's K = 7 runs were still going. Stop and report
both outcomes.
