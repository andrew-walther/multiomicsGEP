# Codex prompt: survival weight set by a target share (S-16)

Date: 2026-11-13. Author: Andrew Walther.

Andrew approved this plan on 2026-11-13:
[`docs/plans/2026-11-13-share-targeted-survival-weight-plan.md`](../plans/2026-11-13-share-targeted-survival-weight-plan.md).
Read it in full before starting. It is the registration. Do not change its fixed choices,
targets, selection rule or pass criteria. If one turns out to be impossible, stop, write
the review handoff explaining why, and push. Andrew will not be available to answer
questions overnight.

## Background in brief

w_S is a power on the Cox partial likelihood. Along a program's own loading direction,
survival's share of the curvature scales like w_S · b² / (p · PVE). On the real PDAC data it
is about 10⁻⁴ at w_S = 1, and the S-15 grid (w_S up to 50) never moved the loadings. The
plan sets w_S from a target directional share s\* computed from a w_S = 1 pilot fit on the
same training data (with recalibration), makes large weights converge, then tests whether
supervision that strong changes the programs and the training CV C. External scoring
waits for review.

Read first: the plan; `docs/findings/2026-11-04-survival-weight-p2000-findings.md` (S-9);
`docs/findings/2026-11-05-supervision-strength-findings.md` (S-13);
`docs/findings/2026-11-05-near-threshold-findings.md` (S-14);
`docs/findings/2026-11-06-supervised-real-refit-findings.md` (S-15);
`docs/findings/2026-11-06-desurv-style-analyses-findings.md`;
`docs/findings/2026-11-12-k5-exhibits-findings.md` (the K = 5 working fit);
the survival-weight box and Eq. AB-F of
`derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.tex`;
the F update at `code/fit_ssbmf_yfb.R` (canonical sequential sweep) and the
`survival_block_weight` option.

Existing runners to reuse rather than duplicate: `run_yfb_survival_weight_p2000.R` (S-9),
`run_yfb_supervision_strength.R` (S-13), `run_yfb_near_threshold.R` (S-14),
`run_yfb_supervised_real_refit.R` (S-15 folds and real-data CV),
`score_yfb_supervised_external.R` (external scoring), and the S-11 DeSurv-style scripts
listed in their findings note. Extend them through arguments where possible. Add a new
runner only where extending would be awkward.

## Order and priorities

Andrew wants results that can go into prelim Chapter 4 when he reviews it on
Saturday 2026-11-14 or Sunday 2026-11-15. Work in this order and commit and push after each
step, so partial results are usable:

1. Step 1 (rule check and `share_targeted_ws()`). No new fits.
2. Step 2 (diagnose and fix large-weight non-convergence).
3. Step 3 (simulation screen and null check).
4. Step 4 (real-data training CV) and Step 5 (program structure).
5. Step 7 only if time remains.
6. Write the findings note and review handoff, push, and stop. **Do not run Step 6
   (external scoring).** It runs only after the reviewing session and Andrew approve,
   because the one allowed external scoring should not be spent on an unreviewed fitter
   change.

Stop conditions are in the plan. In short: if Step 2 finds no fix, continue only at the
targets whose weights converge without one. A target that fails the Step 3 null check does
not go to real data. Do not tune s\*, the cap, the selection rule or the program rule after
seeing results. Put any idea for doing so in the findings note as a proposed follow-up.

## Model routing

- **Planner and orchestrator: GPT-6.1 Sol** (Opus 5.5 if run in Claude Code). It reads the
  plan, keeps the step order and stop conditions, makes every judgment call, does the
  derivation for any Step 2 fix, checks each subagent's result before relying on it, and
  writes the findings note and review handoff.
- **Lower-level execution: GPT-6.1 Terra subagents** (Sonnet 5.5 in Claude Code). Use them
  for bounded, well-specified work: extending runners, writing tests, launching fit grids,
  checking that outputs are complete, building aggregate tables and figures, and running
  the test suite.
- Give each subagent a complete brief: the goal, the files, the fixed settings from the
  plan, what is already ruled out and what done looks like. Subagents do not choose
  targets, change thresholds, decide stop conditions or interpret results for the chapter.
- The independent reviews of the code and of the findings use a fresh Sol (Opus) agent that
  sees the plan and the work but not your conclusions.

## Rules

- Script headers and roxygen per AGENTS.md, author Andrew Walther, no tool names.
- Work in this checkout on a new branch `s16-share-targeted-weight` created from `main`.
  No worktree. Push the branch. Do not merge it into `main`: AGENTS.md requires review of
  supervision work before merging.
- Every new fitter control is opt-in. Fitter defaults do not change. At w_S = 1 with new
  controls off, results must match the current fitter exactly (add a test).
- Run `Rscript tests/run_tests.R` after any change to an update or the fitter, and
  `verify_derivation.R` in the YFB derivation folder if an update's algebra changes. All
  tests must pass.
- Use `timefix = FALSE` in every `coxph` call on simulated times. Evaluate held-out Cox
  scores with `frozen_reverse_cindex()` and the documented direction. Keep all records from
  each PACA-AU donor together in the pooled external interval.
- Weights and pilots come from training data only (the fold's training data inside CV).
- Save E[L²] and E[F²] in every new checkpoint so precision shares are exact.
- Patient-level inputs, outputs and fit checkpoints go under the ignored `data/` folder.
  Commit only aggregate outputs.
- At most 2 CPU workers. Record wall-clock times so the remaining steps can be scheduled.
- Set `max_iter` to 500 for every S-16 fit by overriding `cfg$cavi$max_iter` in the
  runners, as the plan says.
- Do not edit DECISIONS.md, PROJECT_STATUS.qmd, ROADMAP.md, the week plan, the progress
  book or the prelim chapter. Put proposed text for them in the findings note.
- Index every new document in its folder's README (date, link, purpose, status).
- Stage only your own paths (`git add <paths>`, never `git add -A`). Commit each step
  separately with a detailed message and no tool attribution, then push.
- Get an independent review of the code changes and of the findings before the final
  commit, and report its result in the findings note.
- Write prose plainly: no semicolons or dashes as punctuation where a period or comma
  works, no metaphors, headers that say what follows.

## What to hand back

When you finish, or when you stop at a stop condition, write two documents and push them:

1. **Findings:** `docs/findings/2026-11-13-share-targeted-weight-findings.md`. Summary
   first (what was run, what passed, what was selected), then one section per step with
   its tables, the exact commands to reproduce each output, deviations from the plan with
   reasons, limits, and a final "Text for the shared status files" section. Include a
   short draft paragraph that could go into the progress-book chapter and Chapter 4,
   labelled as a draft.
2. **Review handoff:** `docs/reviews/2026-11-13-share-targeted-weight-handoff.md`, addressed
   to the reviewing Claude Code session. It must contain: the commits (hashes) per step;
   what to check first; each pass or fail verdict against the plan's registered criteria,
   with the number it rests on; every place you used judgment the plan did not fix; any
   numbers you are unsure of; open questions for Andrew; and which results are ready for
   Chapter 4 and which are not.

Then stop. Do not start follow-up studies. The reviewing session will check the findings
before anything is adopted or written into the chapter.
