# Codex prompt: independent review of task B (β prior)

Date: 2026-11-12. Author: Andrew Walther.

Paste into a new Codex session. It can run alongside other sessions; it reads, reruns two
fits and writes one review note.

---

Work in `~/GithubProjects/multiomicsGEP` on `main`. Do not create a worktree or change
branches. You are an independent reviewer of task B. Form your own verdict from the plan,
code and outputs before reading the task B findings note
(`docs/findings/2026-11-12-beta-prior-comparison-findings.md`), the 2026-11-12 "Task B:
Normal β stays" DECISIONS.md entry, the "B result" attempts-log entry or the PROJECT_STATUS
section. Then compare and record any disagreement.

**Read.** "Task B. β prior" and "Reporting rules for all tasks" in
`docs/plans/2026-11-09-week-plan.md`; the "B: β prior (registered before execution)" entry
in `docs/plans/2026-10-16-round2-attempts-log.md`; the 2026-11-12 correction entry for
Normal β in DECISIONS.md; `results/benchmark_sim/run_lf_prior_comparison.R` (`lfp_fit`,
`lfp_pip`, `lfp_sim_task`, `lfp_real_task`, `lfp_summary_beta`) and the β update in
`code/fit_ssbmf_yfb.R` and `code/update_beta.R` for `beta_prior_mode = "shared"`; outputs
`beta_*.csv` and, for the Normal arm, A's `sim_fits.csv`/`real_fits.csv` in
`results/benchmark_sim/outputs/lf_prior_2026-11-09/`. Patient-level checkpoints are in
`data/lf_prior_2026-11-12/`; never commit them.

**Check.**

- The arms differ only in the β prior and mode. Is the choice of shared mode for the
  sparse priors justified, and is it disclosed?
- The inclusion probability formula in `lfp_pip()` matches the point-normal EBNM
  posterior for the prior used in the final sweep. Recompute it for one fit from its
  saved prior trace.
- Recompute the rule's three parts per sparse prior from the CSVs, including the version
  without the S-14 cell and whether Normal meets the false-activity limit there.
- Rerun one shared point-normal simulation fit and one real fold from scratch and confirm
  the saved metrics.
- Note any deviation from the plan, undisclosed limitation, or anything that would change
  the decision.

**Output.** `docs/reviews/2026-11-12-task-b-independent-review.md` with your verdict, each
check, the reruns, errors with file and line, and then the comparison with the findings.
Add a row to `docs/reviews/README.md`. Do not edit DECISIONS.md, PROJECT_STATUS.qmd, the
week plan or the attempts log; put proposed changes in a final section. Stage only your
two files, commit with a detailed message and no tool attribution, push. At most 2 CPU
workers. Stop and report.
