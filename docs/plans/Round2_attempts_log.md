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

Task 7 housekeeping completed at this pre-run pause: fourteen process records
relocated, both renders successful with native Quarto, protected deliverables
unchanged, full suite 564/0. No scientific configuration or attempt design changed.
Task1-2 is still the next scientific execution.

**Result, Task1-2:** stopped the coordinator and six identified workers before
editing code (exit143). Saved 209 cell rows: 48 complete metric rows and 161
reporting errors, with zero malformed checkpoint files. The error was
`invalid argument type` from applying `!` to NULL when no pruning history exists.
It occurred after fitting; these are reporting errors, not established model
failures. Every raw checkpoint is retained; `results_partial.csv` and
`execution_status.md` in the attempt2 folder document the partial execution.
No adoption criterion is assessed from this incomplete/error-affected panel.

## Task 1 attempt 3 — repair absent-history reporting; unchanged science

**Before execution:** pre-registered and committed before running.
**Hypothesis/exact change:** coerce absent pruning-history flags to logical(0)
before negation in the cap counter, so unpruned fits report zero intermediate
caps instead of losing metrics. Apply the same two-line reporting fix to both
models. NA flags remain NA; no warnings are silently removed.
**Arms and pass/fail:** exactly Task1-1's 800 cells and every unchanged criterion.
No fitter, prior, threshold, test bed or stopping rule changes. A real-fit
integration test reproduces the prior reporting error before the fix.
**After this entry's commit:** same command, fresh
`results/benchmark_sim/outputs/k_selection/round2/attempt3/`. Task1-1 and Task1-2
outputs remain untouched. Scientific results pending.

Task1-3 verification before execution: real-cell reporting test RED→GREEN,
4/4 panel tests; fresh full suite **565 passed, 0 failed**.

Task1-3 interim single portion (all320 cells, no errors): protected recovery2/20
in each of eight groups; ELBO convergence20/20, false-active margins all pass,
medianK2 (noise.25)/1 (noise1) at every startingK. Recovery and rank fail, so
no adoption. Single's objective is already variance-corrected; no new approved
ladder rung is available. Issue needs Andrew; see the interim handoff and
single_summary_interim.csv. Multimodal results are still pending; no run was
changed or terminated on this assessment.


### Andrew's clarification while Task1-3 was running

Andrew approved option 1: Task5 may use the unchanged K7/no-pruning single
fallback. He also asked that supervised K selection remain an active objective,
recognizing the strong earlier alpha_F=0 result. The failed protected-pruning
panel does not invalidate that benchmark or the unpruned supervised simulation
result (17/20 and 19/20 recovery). Any viable prior remedy from Task2 can motivate
rechecking supervised K selection under the unchanged criteria; this approval
does not authorize changing a threshold, test bed, default or adoption rule.
The recorded 2026-07-13 real-data K4/K5 warm-start alternatives have K_eff=2 and
mean external C=.6270, versus K7's .6267/K_eff=2. The older K_eff=4 beta-boost
entry is explicitly superseded in DECISIONS.md. Retained rank and survival-active
rank must remain distinct in subsequent reports.

## Task 1 attempt 4 — observe the single-modality pruning transition

**Before execution:** register and commit before running this diagnosis. The
completed 320-cell single portion of Task1-3 has 2/20 recovery in every group,
although the earlier unpruned selected model recovered17/20 and19/20. Andrew
asked that supervised K selection continue while approving the unchanged K7
fallback for Task5. The multimodal part of Task1-3 is still running unchanged.

**Hypothesis/exact change:** instrumentation only. Record each posterior state
presented to the existing single pruning nullcheck and its exact returned gains.
This distinguishes loss of credible protection before deletion from loss during
a subsequent refit. No fitter source, fitting argument, moment initialization,
prior, threshold, objective or convergence rule is changed. The observer is a
process-local binding, restored on success or error. It consumes no RNG draws.
Saved real-data supervised fit has a nondegenerate shared prior; do not presume
that a point-mass prior explains every single-modality pruning loss.

**Arms/data:** no new statistical arm; exactly the selected protected single
configuration from Task1-3, first prescribed seed20261016, K7, noise.25 and1,
same calibrated n400/p50 data, original caps/tolerance. These two fixed case
traces diagnose a transition, not its frequency or a new adoption rule. No seed
substitution is allowed if a case does not exhibit the failure. No external data.

**Pass/fail/report:** record loading correlation/recovery, matched beta/SD,
credible protection, matched deletion gain, actual dropped program, shared
prior slab probability/scale and convergence at each pruning stage. Identify
the first loss only if observed. No remedy is adopted from a case study.
The unit test must show an identical fit and identical nullcheck gains, and
restoration of the original binding after both success and error. The existing
attempt folder must refuse overwrite. Full suite must pass before committing.

**Execution after commit:**
`Rscript results/benchmark_sim/run_single_pruning_trace.R`;
fresh `results/benchmark_sim/outputs/k_selection/round2/attempt4_trace/`.
No running panel or fitter source is edited.

Task1-4 pre-run verification: the state-capture assertion failed against the
non-observing stub, then passed with an identical full fit and exact nullcheck
gains. Original bindings restore after success and error; old attempt folders
refuse overwrite. Fresh full suite **567 passed, 0 failed**. Task1-3 is still
running, with 783/800 checkpoints and no recorded execution errors; its source
has not been edited. This diagnosis does not change that panel or its design.
