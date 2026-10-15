# Round 2 attempt history

Authority: `docs/plans/Codex_Handoff_Round2_10_16_26.md` and Andrew's request.
Log each attempt here and commit its design before running it. Append outcomes;
never remove or rewrite a failed attempt. Defaults, thresholds and test beds stay
fixed. Selection uses simulation and training data only. New real-data tables
must include the immutable benchmark snapshot rows.

## Benchmark and integration record

Before Round 2, all ten SHA256 entries passed in
`results/benchmark_snapshot_2026-10-14/`. The annotated benchmark tag resolves to
`a352f17ca9cff523caa80756aa2af085e62551e1`; neither tag nor snapshot was changed.
Created `round2/supervised` from clean `feature/supervised-F` at `2bebbd4`.
Task 0 integration completed in the required order with `--no-ff` merges,
preserving all entries and tests. Successive full-suite counts were 553, 554,
556, 559 and 559, all with zero failures. No scientific attempt has run yet.

## Task 1 attempt 1 — positive-SD protection and the prescribed rank panel

**Status before execution:** pre-registered; no panel fit has run. This design
and its implementation are committed before execution.

**Hypothesis:** excluding zero-variance β point masses from protection allows
dead programs to be removed while retaining credibly prognostic programs.
**Exact scientific change:** protection now requires |Eβ|≥1.96 SD(β) **and**
SD(β)>0 in both nullchecks. Protection remains opt-in. Default/off paths retain
bit-for-bit identity to `feature/supervised-F` in regression tests. The specified
convergence clarification is applied to assessment, not to fitting defaults.

**Arms and fixed data:** 800 cells, seeds 20261016:20261035, noise {.25,1},
K_init {5,7,10,15}, calibrated n_train400/n_validation120, p50 per modality,
β3=3.3. Single: selected model with protected pruning; α_F=0 reference with
signed point-Laplace L/F. Multimodal: selected model with/without protection,
and 10/1 framework reference. All other settings match the handoff.

**Pass/fail:** in every K/noise group, recovery≥16/20, mean false prognostic
programs≤that model's reference+.5, 20/20 converged, median retained K=3.
Single uses its ELBO flag; multimodal uses its output-change flag. Single's
last-iteration output criterion and intermediate pruning caps are also reported;
neither changes the adoption rule. The returned fit's convergence flag is the
criterion, with pre-pruning caps exposed separately. Errors/missing cells cannot
pass. No threshold, default, test bed or parameter tuning is permitted.

**Execution:** after this entry's commit, run
`Rscript results/benchmark_sim/run_k_selection_panel.R 6`. Each cell, including
errors/warnings, checkpoints independently to
`results/benchmark_sim/outputs/k_selection/round2/attempt1/cells/`; aggregates,
summary and adoption assessment are in that attempt folder. The existing Round 1
outputs are preserved. No real-data cohort is loaded. Conditional training-only
1-SE fitting is eligible only if single-modality simulation adoption passes.

**Verification before execution:** dead-boundary tests failed against old code,
then 23/23 protection tests and 2/2 panel tests passed. Full-suite result is
563 passed, 0 failed (fresh full suite before the pre-run commit). Results pending.

**Result, Task1-1:** execution failed before the grid or any fit was dispatched.
Exact error: `cannot change value of locked binding for 'grid'`, from `grid <<-`.
The function's superassignment searched into the graphics package rather than
creating a local grid. Zero scientific cells completed; no criterion was assessed.
Output: `outputs/k_selection/round2/attempt1/execution_status.md`.

## Task 1 attempt 2 — repair runner scope; unchanged scientific design

**Before execution:** pre-registered and committed before retrying.
**Hypothesis and exact change:** local configuration variables plus an explicit
`clusterExport(..., envir=environment())` prevent inherited locked bindings from
intercepting the entry point. Add an output-folder argument solely for isolated
entry tests/fresh attempts. No fitter, data, arm, metric or adoption rule changes.
**Arms and pass/fail:** exactly the 800 cells and every criterion of Task1-1.
A zero-worker entry-point test failed before the fix and now must reach worker
setup after writing the 800-row design. The previous failed launch is preserved.
**Execution after this entry's commit:** same command, fresh
`results/benchmark_sim/outputs/k_selection/round2/attempt2/`. No overwriting,
external cohorts or tuning. Scientific results pending.

Task1-2 verification before execution: entry tests 3/3; fresh full suite
**564 passed, 0 failed**. Execution waits until the added Task 7 housekeeping
finishes; no scientific cell has yet run in either attempt.
