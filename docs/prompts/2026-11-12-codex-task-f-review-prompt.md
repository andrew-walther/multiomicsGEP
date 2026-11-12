# Codex prompt: independent review of task F (training-CV redo)

Date: 2026-11-12. Author: Andrew Walther.

Paste into a new Codex session. It reads, reruns two fits and writes one review note.

---

Work in `~/GithubProjects/multiomicsGEP` on `main`. Do not create a worktree or change
branches. You are an independent reviewer of task F. Form your own verdict before reading
`docs/findings/2026-11-12-training-cv-redo-findings.md`, the 2026-11-12 "Task F" entry in
DECISIONS.md, the "F result" attempts-log entry or the PROJECT_STATUS section. Then compare.

**Read.** "Task F" and "Reporting rules for all tasks" in
`docs/plans/2026-11-09-week-plan.md`; the "F: training-CV redo ... (registered before
execution)" entry in `docs/plans/2026-10-16-round2-attempts-log.md`; the 2026-05-27 D1–D5
definitions in `docs/findings/2026-05-27-desurv-alignment-report.qmd` and
`results/benchmark_sim/run_desurv_comparison.R`; the runner
`results/benchmark_sim/run_training_cv_redo.R` (it sources
`run_lf_prior_comparison.R`); outputs in
`results/benchmark_sim/outputs/training_cv_redo_2026-11-12/`. Patient-level files are in
`data/training_cv_redo_2026-11-12/` and `data/training_cv_redo_inputs_2026-11-12.rds`;
never commit them.

**Check.**

- Is reducing D1–D5 to the D2 and D4 screens under YFB justified by the plan and the
  2026-05-27 definitions? Does `f_prepare()` reproduce the D2 preprocessing exactly, and
  are all screens outcome-free?
- Are the three scaling inputs built as registered, and is `center_projection` used only
  for the uncentered X/sd input?
- Does reusing task A's point-Laplace checkpoints as the D4/z-score/K = 7 reference pair
  correctly with the new fits (same folds, same call)?
- Recompute the screen and scaling decisions and the one-SE K selection from `fits.csv`.
  Is the one-SE rule implemented as registered? Comment on unpaired vs paired SE and on
  how decisive K = 5 is.
- Is the task A rerun at K = 5 the registered one?
- Rerun one rank-sweep fit (K = 5, one fold) and the D2 fit for one fold from scratch and
  confirm the saved metrics.
- Note any deviation, undisclosed limitation, or anything that would change a choice.

**Output.** `docs/reviews/2026-11-12-task-f-independent-review.md` with your verdict on
each of the three choices, each check, the reruns, errors with file and line, then the
comparison with the findings. Add a row to `docs/reviews/README.md`. Do not edit
DECISIONS.md, PROJECT_STATUS.qmd, the week plan or the attempts log; put proposed changes in
a final section. Stage only your two files, commit with a detailed message and no tool
attribution, push. At most 2 CPU workers. Stop and report.
