# Codex prompt: stronger supervision with a fixed prior on β (S-17)

Date: 2026-11-14. Author: Andrew Walther.

Andrew approved this study on 2026-11-14. The registration is
[`docs/plans/2026-11-14-bounded-risk-supervision-plan.md`](../plans/2026-11-14-bounded-risk-supervision-plan.md).
Read it in full first. Do not change its fixed choices, grids, gates or selection rule.
If one is impossible, write the review handoff explaining why, push, and stop. Andrew will
not be available to answer questions while this runs.

## Background in brief

S-16 (merged into `main` at `df32eb1d`) found that survival weights large enough to move
the real-data loadings never converge: the ELBO keeps rising while |η| and |β| grow.
Canonicalization makes β a log hazard ratio per SD of the program score, and its prior
variance is re-estimated every sweep, so nothing bounds it. S-17 fixes the β prior (ridge
by default, Laplace as a comparison), the Bayesian counterpart of DeSurv's fixed
elastic-net penalty, and repeats the share-targeted weights under that bound. Read the
reviewing session notes at the end of
`docs/reviews/2026-11-13-share-targeted-weight-handoff.md` before starting. They explain
why S-16's null gate is replaced.

Reuse the S-16 code: `code/share_targeted_ws.R`, `results/benchmark_sim/run_share_targeted_weight.R`,
`run_share_targeted_diagnostics.R`, `share_targeted_structure.R` and their tests. Extend
them through arguments (β prior mode, fixed prior, σ_b). Add new files only where extending
would be awkward. The fixed prior already exists in the fitter: `beta_prior_mode = "shared"`,
`beta_prior_fix = TRUE`, `beta_prior_init`. These pass through `fit_control`, warm starts
keep the fixed prior, and checkpoints already save the needed moments. A review of the
plan found these runner gaps, which you need to close with opt-in arguments and tests:

- `s16_fit()` and `s16_diagnostic_fit()` hardcode `prior_beta = "normal"`. The Laplace arm
  needs `prior_beta = "point_laplace"`.
- One `fit_control` applies to every arm of a dataset or fold. The working-model reference
  (estimated per-program prior) needs its own control.
- `s16_calibrate()` does not trace, and the diagnostics runner has no recalibration, a fixed
  2026-11-13 output folder, and a stop when a checkpoint already exists. Step 1 needs
  traced, recalibrated fits written to the S-17 output folder.
- The trace needs the per-program score SD before canonicalization and β before and after
  it, plus the stable-convergence check defined in Step 1.

## Order

1. Step 1, the convergence screen. Commit and push.
2. Step 2, real-data training CV and program structure, for the Normal σ_b values that
   passed. Commit and push.
3. Step 3, the simulation screen and null check. Commit and push.
4. Write the findings note and review handoff, push, and stop. **Do not run Step 4
   (external scoring).** Do not start follow-up studies.

Real data comes before simulation so that Andrew has the real-data answer first. The
adoption rule in the plan still requires the Step 3 null check, so do not call any arm
adopted or recommended before Step 3 is done.

## Model routing

- **Planner and orchestrator: GPT-6.1 Sol** (Opus 5.5 in Claude Code). It keeps the step
  order and gates, makes every judgment call, checks each subagent's result before relying
  on it, and writes the findings and handoff.
- **Execution: GPT-6.1 Terra subagents** (Sonnet 5.5 in Claude Code), with complete
  briefs, for runner extensions, tests, fit grids, completeness checks and tables.
  Subagents do not change grids, gates or rules or interpret results.
- Independent reviews of the code and of the findings use a fresh Sol (Opus) agent that
  sees the plan and the work but not your conclusions.

## Rules

- Script headers and roxygen per AGENTS.md, author Andrew Walther, no tool names.
- Create the branch `s17-bounded-risk-supervision` from `main` in this checkout. No
  worktree. I authorize pushing commits on this branch to origin
  (`https://github.com/andrew-walther/multiomicsGEP.git`, branch
  `s17-bounded-risk-supervision`) after each step. Do not merge into `main`.
- Fitter defaults do not change. At the working configuration with new controls off,
  results must match the current fitter exactly (keep the S-16 baseline test passing).
- Run `Rscript tests/run_tests.R` after any change to the fitter or an update, and
  `derivations/SSBMF-YFB-model-derivations/verify_derivation.R` if an update's algebra
  changes.
- `max_iter = 500` for every S-17 fit. `timefix = FALSE` in every `coxph` call on simulated
  times. Weights and pilots from training data only.
- Save τ, E[L²], E[F²], E[β²] and the final Cox working weights in every new checkpoint.
- Patient-level inputs, outputs and checkpoints under the ignored `data/` folder. Commit
  only aggregate outputs.
- At most 2 CPU workers, one BLAS thread each. Record wall-clock times.
- Do not edit DECISIONS.md, PROJECT_STATUS.qmd, ROADMAP.md, the week plan, the progress book
  or the prelim chapter. Put proposed text in the findings note.
- Index every new document in its folder's README.
- Stage only your own paths. Commit each step with a detailed message and no tool
  attribution.
- Write prose plainly: no semicolons or dashes as punctuation where a period or comma
  works, no metaphors, headers that say what follows.

## What to hand back

1. `docs/findings/2026-11-14-bounded-risk-findings.md`: summary first (what ran, which
   gates passed, what was selected, whether the programs changed), one section per step
   with tables and reproduction commands, deviations with reasons, limits, "Text for the
   shared status files", and a draft paragraph for the progress book and Chapter 4 labelled
   as a draft. The paragraph must say what the arms mean in terms of supervision strength,
   not only their labels.
2. `docs/reviews/2026-11-14-bounded-risk-handoff.md`, addressed to the reviewing Claude
   Code session: commits per step, what to check first, each gate's verdict with the
   number it rests on, every judgment call the plan did not fix, uncertain numbers, open
   questions for Andrew, and which results are ready for Chapter 4.

Then stop.
