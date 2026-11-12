# Codex prompt: independent review of task A (L/F priors)

Date: 2026-11-12. Author: Andrew Walther.

Paste into a new Codex session. It can run alongside the other Codex sessions; it only
reads, reruns a few fits and writes one review note covering task A and follow-up A2.

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

**Follow-up A2 (review it too, after task A).** The same blind-reading rule applies: do not
read the A2 findings note (`docs/findings/2026-11-12-lf-prior-followup-findings.md`), the
"A2 result" attempts-log entry or the A2 paragraph in DECISIONS.md until you have your own
verdict. Read the registration, the 2026-11-12 "A2: point-exponential with cohort
differences removed (registered before execution)" entry in the attempts log, the runner
`results/benchmark_sim/run_lf_prior_followup.R` (it sources the task A runner and reuses its
helpers), and the outputs in `results/benchmark_sim/outputs/lf_prior_2026-11-09/a2/`
(checkpoints in `data/lf_prior_2026-11-12/a2/`). Two setups, both point-exponential L and F:
a shifted per-cohort z-score (SZ), pooled, and X/sd with `strata_id` = cohort (PES).

- Does the code implement the registered setups? SZ adds one per-gene constant
  c_j = −min_i z_ij to the per-cohort z-score. In the synthetic data the training z-score is
  shifted by its own minimum and the validation z-score by the training constants. PES passes
  the training cohort labels as `strata_id` with `center_projection = TRUE`.
- The SZ constants are computed once from all 273 patients, held-out folds included. Decide
  whether this can affect any held-out metric.
- PES was not refitted in simulation. Its criteria 1 to 3 are taken from task A's PE
  simulation fits, on the argument that one stratum with `center_projection = TRUE` is the
  same fit. Check that argument in `code/fit_ssbmf_yfb.R`.
- Recompute the four criteria for each setup from `a2_sim_fits.csv`, `a2_real_fits.csv` and
  task A's `real_fits.csv` (arm PL), and compare with `a2_criteria.csv` and
  `a2_real_paired.csv`.
- Rerun one SZ simulation fit and one PES real-data fold from scratch and confirm the saved
  metrics.
- Check the false-active coefficient ratios (|E β_j| / |E β_k| to the matched program) in
  `a2_sim_fits.csv`. They are a description made after seeing the results, not a criterion.
- Check whether the PES programs match task A's PE programs, using the saved E[F] in both
  checkpoint folders. Compare both the directions and the column norms.

**Output.** Write `docs/reviews/2026-11-12-task-a-independent-review.md` with: your
verdict on the decision (keep point-Laplace or not), each check and its result, the four
reruns (two for task A, two for A2), any errors with file and line, and then the comparison
with each findings note. Give a separate verdict for A2: does either setup qualify?
Add a row to `docs/reviews/README.md`. Do not edit DECISIONS.md, PROJECT_STATUS.qmd, the
week plan or the attempts log; put any proposed changes in a final section of your review.
Stage only your two files (`git add <paths>`, never `git add -A`), commit with a detailed
message and no tool attribution, and push. At most 2 CPU workers. Stop and report.
