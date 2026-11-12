# Codex prompt: independent review of task A (L/F priors)

Date: 2026-11-12. Author: Andrew Walther.

Paste into a new Codex session. It can run alongside the other Codex sessions; it only
reads, reruns a few fits and writes one review note.

---

Work in `~/GithubProjects/multiomicsGEP` on `main`. Do not create a worktree or change
branches. You are an independent reviewer of task A. Form your own judgment from the
code, data and plan. Do not read the task A findings note or the 2026-11-12 task A
DECISIONS.md entry, A result attempts-log entry or PROJECT_STATUS section until you have
written your own verdict. Then compare and note any disagreement.

**What to read.**

1. The plan: "Reporting rules for all tasks" and "Task A" in
   `docs/plans/2026-11-09-week-plan.md`, including the dated 2026-11-12 revision note
   (pooled training, no strata, within-cohort concordance as the primary CV metric).
2. The registration: the 2026-11-12 "A: point-exponential vs point-Laplace priors on L and F
   (registered before execution)" entry in `docs/plans/2026-10-16-round2-attempts-log.md`,
   and the "Task A fits pooled training data" entry in DECISIONS.md.
3. The code: `results/benchmark_sim/run_lf_prior_comparison.R` and what it calls
   (`run_yfb_survival_weight_p2000.R`: `s9_fit`, `s9_score`, `s9_zscore`;
   `run_yfb_near_threshold.R`: `s14_data`; `run_yfb_supervision_decision.R`:
   `yfb_decision_data`; `code/fit_ssbmf_yfb.R` for `center_projection` and
   `prior_LF = "point_exponential"`).
4. The outputs: `results/benchmark_sim/outputs/lf_prior_2026-11-09/` (all CSVs).
   Patient-level checkpoints are in `data/lf_prior_2026-11-12/` and the prepared input in
   `data/lf_prior_real_2026-11-12.rds`; read them locally if needed, never commit them.

**Check.**

- Does the code implement the registered arms, cells, folds, metrics and criteria? In
  particular: PE input is X/sd within cohort without centering and nonnegative; PL is the
  working-model call; `center_projection = TRUE` only for PE and RANK; the S-11 folds are
  used unchanged; within-cohort concordance pairs patients only within the same cohort
  and the score orientation is not chosen on held-out data.
- Is the recovery, false-active and calibration scoring correct and identical across arms?
  Is a failed or missing fit counted as a failure, not dropped?
- Recompute the four criteria from `sim_fits.csv` and `real_fits.csv` yourself and compare
  with `criteria.csv` and `real_paired.csv`.
- Rerun two fits from scratch in a scratch directory (one PE simulation fit and one PE
  real-data fold) with the runner's functions and confirm the saved metrics.
- Check the cohort-label and mean-expression correlations in `real_programs.csv` and
  whether they support or contradict any explanation for PE's real-data behavior.
- Check the PL reuse comparisons (`pl_reuse_check.csv`, `pl_s11_reuse_check.csv`) and
  whether any disagreement is a fit difference or a scoring-definition difference.
- Note any deviation from the plan, any undisclosed limitation, and anything that would
  change the decision.

**Output.** Write `docs/reviews/2026-11-12-task-a-independent-review.md` with: your
verdict on the decision (keep point-Laplace or not), each check and its result, the two
reruns, any errors with file and line, and then the comparison with the findings note.
Add a row to `docs/reviews/README.md`. Do not edit DECISIONS.md, PROJECT_STATUS.qmd, the
week plan or the attempts log; put any proposed changes in a final section of your review.
Stage only your two files (`git add <paths>`, never `git add -A`), commit with a detailed
message and no tool attribution, and push. At most 2 CPU workers. Stop and report.
