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

**Result, Task1-4:** both pre-specified cases completed without errors or
warnings, and every observed fit converged. At noise .25 the K = 7 fit recovered
the planted program (|cor| .895844,beta−.914905,SD.056769), protected it with
nullcheck gain−Inf, and deleted a different program3. The K = 6 refit lost recovery
(best |cor| .390941); the terminal K2 fit has |cor| .477163. At noise 1 recovery
survives (.958283→.715245 after the first refit→.734053 at terminal K1), so the
rank criterion fails there independently of recovery. Effective shared-prior
slab probability stays positive in both cases (.5 initially, near1 terminal).
A point-mass beta prior is not the cause observed in these cases. Source audit:
mean-only refits reset EL2/EF2 to squared means, beta moments to zero, Tau to
initial feature variances, and the persistent loading-prior state to empty.
Output: `outputs/k_selection/round2/attempt4_trace/stages.csv`; raw state objects
are ignored RDS files. This two-case diagnosis estimates no population frequency.

Andrew explicitly approved testing full-posterior pruning refits as an additional
opt-in remedy beyond the listed Task1 ladder. All criteria/defaults remain fixed.
Implementation/tests are being prepared in separate temporary copies while the
unchanged original panel finishes; no candidate scientific fit has run yet.

## Task 1 attempt 5 — full posterior state in single protected-pruning refits

**Before execution:** this entry and the exact runner/opt-in implementation are
committed before any candidate scientific fit. Andrew explicitly approved this
additional remedy after Task1-4 showed the protected program was lost during a
refit after a different program was deleted. This approval adds a remedy; it
changes no criterion, threshold, data, test bed or default.

**Hypothesis/one change:** preserving the retained posterior state prevents the
refit initialization from erasing the prognostic direction. Carry EL/EL2, EF/EF2,
beta/second moments, feature precisions, persistent F priors and the shared beta
prior into each refit. The L prior is fit from its pseudo-observations on every
coordinate update and has no separately held state. Keep the existing refit
warm-up 0, update order, EBNM calls, objective and nullcheck unchanged. The new
`prune_preserve_posterior=TRUE` path is opt-in; FALSE/default stays bit-identical.

**Arms/data:** all 160 protected single cells, seeds20261016:20261035, noise .25/1,
K_init5/7/10/15, unchanged n = 400 / p = 50,β₃ = 3.3 calibrated bed, selected signed
point-Laplace/shared-beta configuration. Reuse the 160 unchanged alpha_F=0
reference rows and retain the 160 old mean-only protected rows from Task1-3.
The new runner refuses to start until the original800-cell panel is complete,
unique and error-free; it records that source table's checksum. No external data.

**Adoption rule:** the same four conditions in all eight K/noise groups:
recovery>=16/20; mean false prognostic<=reference+.5; all 20 final fits converge
under the single ELBO rule; median K=3. Single last-iteration output flags and
all intermediate refit caps remain separate reported diagnostics as in Task1-3.
Missing/error cells cannot pass. No additional seed case, threshold, start,
prior or fitting tolerance may be selected from these results. If simulation
adopts, run the handoff's training-only 1-SE analysis with this same candidate.
The benchmark can change only after Claude review and Andrew approval.

**Software verification:** original API fails the new tests as unsupported.
Temporary implementation tests verify exact state subsetting, nonzero retained
beta in the first working predictor, retained F variance in beta precision,
FALSE/default identity, malformed-state errors and exact recursive state transfer
using a deliberately dead-column unit fixture. The scientific bed is untouched.
The adapter test verifies the flag is forwarded only when requested. Full suite
and base-branch identity tests must pass after the original panel stops and the
prepared copies are applied, before this attempt's commit.

**Execution after commit:**
`Rscript results/benchmark_sim/run_posterior_k_panel.R 6`;
fresh `outputs/k_selection/round2/attempt5_posterior/`. Six socket workers,
one serial fit per worker, scheduling chunks of1 (no scientific change).
Longleaf160-cell wrapper prepared; no jobs submitted. Results pending.


**Full result, Task1-3:** exited 0;800/800 unique cells,errors 0,warnings 0,93.23min.
Single protected recovery2/20 and median K2(.25)/1(1) at all four starting K;
false margins and all returned ELBO flags pass. Four intermediate Single stages
cap. MM protected recovered19/20;20/20(K = 5),18;20(K = 7),17;19(K = 10),18;20(K = 15).
Every MM false-positive margin passes. MedianK3 exceptK = 5, noise = 1=2. Final output
convergence20/20 exceptK = 15, noise = 1=19/20 (seed 20261033,finalK = 5,cap 500,
cor .934102). MM protected also has one intermediate cap. Reference MM has zero
final/one intermediate cap; ordinary selected MM one final/one intermediate.
Officialadoption.csv rejects Single(recovery,rank) and MM(convergence,rank).
No newK/protection is adopted; no 1-SE run. Full tables are in the attempt3 folder
and the interim handoff/progress chapter. This completes the prior pending MM
portion; no failed attempt or earlier interim entry was rewritten. Andrew's
approved Task1-5 remedy continues for Single; Task2 will address MM collapse.

Task1-5 final pre-run verification: **573 passed, 0 failed** after applying the
prepared code only once the original panel and workers had exited. All four
posterior-state tests and six panel/adapter tests pass; existing regression
checks confirm bit-for-bit base-branch/default identity. A standalone regression
invocation initially omitted the prediction module (one harness error); loading
the master suite's normal dependencies gives17/17. The full suite passes after
a final logging-only clarification and matching the adoption CSV's established
`adopted` column. No candidate scientific fit has run before this commit.


**Result, Task1-5 (after execution):**

### Task1-5 completed: full-state refits improve recovery but are not adopted

All 160 new fits completed with zero execution errors or warnings. The comparison
retains all 160 original reference and 160 mean-only protected rows. Full-state
recovery is 15–19/20 instead of 2/20, and every false-positive margin passes.
However, K_init=5/noise=.25 misses the recovery threshold (15/20); five final
fits hit the 300-iteration ELBO cap; high-noise median K is 1 at every starting K.
The candidate therefore fails recovery, convergence and rank. No 1-SE run is
eligible, no configuration is adopted, and the benchmark/defaults stay unchanged.

| K_init | Noise | Recovery /20 | False active | Reference + .5 | Final ELBO converged /20 | Median K | Intermediate caps |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 0.25 | 15 | 1.15 | 1.75 | 19 | 3 | 2 |
| 5 | 1 | 18 | 0.20 | 0.90 | 20 | 1 | 0 |
| 7 | 0.25 | 17 | 1.10 | 1.65 | 19 | 3 | 4 |
| 7 | 1 | 19 | 0.15 | 0.85 | 19 | 1 | 5 |
| 10 | 0.25 | 17 | 1.15 | 1.80 | 19 | 3 | 6 |
| 10 | 1 | 19 | 0.20 | 0.90 | 19 | 1 | 8 |
| 15 | 0.25 | 16 | 1.10 | 1.80 | 20 | 3 | 0 |
| 15 | 1 | 19 | 0.30 | 0.90 | 20 | 1 | 0 |

The five final caps are seed20261026 at noise .25 with K_init5/7/10 (final K3),
and seed20261017 at noise1 with K_init7/10 (final K2). Intermediate caps total25,
reported separately from returned-fit flags. Outputs:
`results/benchmark_sim/outputs/k_selection/round2/attempt5_posterior/`;
`adoption.csv` records `adopted=FALSE`. Original reference table MD5:
`d45c207936855fe274ad14209474bf34`.

This completes Andrew's additional full-state experiment. The listed Single
ladder has no remaining applicable remedy: its survival objective is already
variance-corrected, and the failure is not convergence alone. Further changes
require a new decision; the approved K7/no-pruning fallback remains available
for Task5. Continue Task2's separate multimodal diagnosis and prescribed ladder.
No external data were scored or used for this assessment.

Task1-5 result checkpoint: fresh full suite **573 passed, 0 failed**. All ten
benchmark manifest checks pass; annotated tag and peeled commit remain unchanged.

## Task 2 attempt 1 — trace the specified TCGA contrasts before choosing a remedy

**Before execution:** this entry, runner header, opt-in controls and tests are
committed before any diagnostic fit. Prepared tiny synthetic unit tests are
software verification and are not used to select a statistical configuration.

**Hypothesis:** refitting the shared point-Laplace beta prior produces a point
mass at zero; subsequent updates under that prior keep beta exactly zero. A
fixed starting prior separates this mechanism from the supervised initialization.
The Normal/warm-up and per-factor historical arms distinguish shared prior
estimation from its family. The no-start contrast isolates supervised initialization.

**Exact implementation change:** opt-in multimodal trace records the prior before
and after each inner sweep, posterior beta means/SDs, raw update x and A,
canonical scale, live loading columns, ELBO parts and the actual output changes.
Additional read-only diagnostics measure full molecular prediction changes
including the intercept, centered molecular changes, and centered Cox predictor
changes. They diagnose intercept allocation without changing any stopping rule.
Single's existing trace opt-in also records its actual sweep prior before/after,
raw coefficient x/s/A and posterior beta moments.
All recursive pruning stages remain in the trace. Explicit trace=FALSE and
beta_prior_fix=FALSE preserve the original object. The specified diagnostic
controls beta_prior_fix=TRUE and beta_prior="shared_normal" are new opt-ins;
no default, objective, stopping rule or threshold changes.

**Arms/data:** K7, seed42, unchanged 144 TCGA matched patients, expression3000
and methylation10000; selected, fixed starting shared Laplace(pi=.5,scale1),
shared Normal(0,1) with the same warmup10 then EB refit, historical per-factor
Normal without supervised start, selected without supervised start. Repeat all
five for the expression block alone with fit_cox_on_yf and its selected settings.
The single per-factor counterpart retains w_S=1 and uses per-program Normal
beta/no supervised start. Single remains unpruned as the handoff fallback;
multimodal retains ordinary unprotected pruning. Use each model's original cap
and convergence rule. The explicit historical reference is not presented as a
one-option contrast. Every other diagnostic contrast changes exactly one option.

**Input authorization:** Andrew approved extracting the existing cache's unchanged
$training block (including its original methylation imputation) to ignored
`data/round2_tcga_training.rds`. Frozen patient/feature IDs and order match.
No ICGC outcomes, predictions or metrics select these arms or any remedy.

**Pass/fail and report:** diagnose the hypothesis only when a point-mass refit
(effective slab probability=sum of weights on positive-scale components=0)
is followed by all beta means/SDs exactly zero in the next and subsequent
updates. Record the first occurrence and any contrary case without assuming
causality. For each cap, report which actual reconstruction/eta change exceeds
the unchanged tolerance, plus centered-eta change as a diagnostic only.
A remedy still requires the later fixed simulation and 5-fold training-CV
criteria; this diagnostic comparison adopts nothing. All warnings/errors are
checkpointed per arm. Comparison tables carry the relevant frozen TCGA benchmark
rows. Raw update pseudo-observations/precisions and canonical scales are distinct
columns, so unit changes remain visible.

**Execution after commit:**
`Rscript results/multimodal_real/run_round2_collapse_diagnosis.R all`.
Fresh output `results/multimodal_real/outputs/round2/diagnosis1/`; existing arm
checkpoints refuse overwrite. Six PSOCK workers, one serial fit each. Longleaf
array wrapper prepared (10 cells,8GB/cell,12h); no jobs submitted. Default identity
against feature/supervised-F and the full suite must pass before this commit.

Task2-1 software verification update (before any real-data fit): the first full
suite returned581 passed/3 failed. All three failures came from the new Single
trace trying to read g_beta in historical per-program fits, where no shared
prior exists. Restricted the additional prior records to shared-prior fits;
existing per-program traces remain unchanged. The affected supervision tests
and new diagnosis tests now pass22/22; the full suite is being rerun. This is a
software correction before the diagnostic attempt, not a scientific retry or
a change to its arms, criteria, inputs or defaults.

Task2-1 final pre-run verification: **584 passed, 0 failed**, including
base-branch/default identity checks, all6 multimodal trace tests and all5
diagnosis tests. Longleaf wrapper syntax passes. No diagnostic real-data fit
has run before this preregistration commit.


## Task3 bounded head-to-head design (prepared; no scientific execution yet)

Approved package/paper SHAs are recorded in DECISIONS. Use the paper's reported
TCGA+CPTAC parameters, not new local BO: its full search is documented as4–8h
on30 cores (50 initial+100 BO evaluations,30 starts,5 folds). Andrew explicitly
authorized this fallback. Parameters: k3,alpha .333935556177224,
lambda .349448737941329,nu .0556512168045738,ntop270;lambdaW=lambdaH=0,
ngene3000,within-sample ranks. Parameter binary SHA256
23ec3fe34cfd306cabfef64ee64edd25d60eddfd0adf756be5d6d3e977d842ec.

Use authors' unchanged raw loader, same ordered eligible TCGA_PAAD/CPTAC
patients/outcomes as D4 (check independently against the D4 project loader),
then authors' outer-union/zero-fill/package preprocess_data(). No second log
is applied in DeSurv; no D4 preprocessing change. Fit100 seeds, each ninit1,
tol1e-7,tol_init1e-5,maxit=imaxit4000,parallel_initFALSE; record errors,
warnings, convergence/nan flags for every seed. The paper keeps successful
non-nan fits for consensus, min_frequency=max(1,ceiling(.3*n_success)),ntop270.
Final fit uses consensusW0/H0/beta0,seedNULL,tol1e-7,maxit4000 exactly as paper04.
No external metric chooses a seed, parameter or orientation.

Refit selected Single on exact D4 preprocessing/patients with set.seed42 and
unchanged selected K7/no-pruning controls. D4 remains the frozen benchmark,
using cached risk and its documented negative-eta orientation. DeSurv and
selected Single use fixed positive-hazard Cox orientation. Score all5 held-out
cohorts once, with each method's pipeline; DeSurv projection is paper05's UNION
of top genes, X'Wbeta. Include missing-gene counts and patient/outcome alignment.
Report per-cohort and mean C, model convergence, and paired per-cohort differences
plus mean-of-cohorts stratified bootstrap B2000/seed1 versus D4 and selected
Single. Frozen benchmark rows come from snapshot, never a D4 refit. Undefined
bootstrap replicates fail loudly, not silently dropped. This is a comparison,
not benchmark replacement. Add a Longleaf full-BO wrapper; no job submission.


Task3 software preparation only (while the unchanged Task2 diagnosis runs):
pinned package installed in isolated ignored library, tiny software fit
converges in22 iterations without nan_flag, and9 focused adapter/driver tests
pass. They verify exact preprocessing/projection, full-precision reported
parameters, seed-fit identity to the author call, per-seed error checkpoints,
no overwrite, raw in-fold BO inputs, invalid-mode failure and the literal
paper03 one-SE postprocessing. A missing/changed cached author helper fails
its checksum check. The full suite and preregistration commit remain required
before any head-to-head scientific fit; none has started. Current preparatory
work does not choose a remedy for Task2 or score external cohorts.


## Task2-1 completed TCGA diagnosis (2026-10-15)

Commit `ccebdf8` preceded all ten fits; its full suite passed584/0. The process
completed with exit0, ten status checkpoints, no errors or warnings. Single
fits took100–102 seconds; multimodal fits took2384–2556 seconds each. All returned
K7, with no pruning refits. Fits and trace CSVs are in
`results/multimodal_real/outputs/round2/diagnosis1/`; real-data RDS files are
ignored. All new training results are shown beside frozen benchmark rows:

| Source/model | Arm | K | Active | Risk SD | Training C | Converged/iterations |
|---|---|---:|---:|---:|---:|---|
| round2/single | selected | 7 | 1 | 0.350339 | 0.654791 | TRUE/297 |
| round2/single | fixed | 7 | 0 | 0.309121 | 0.638645 | TRUE/297 |
| round2/single | normal | 7 | 3 | 0.396799 | 0.661636 | TRUE/298 |
| round2/single | framework | 7 | 3 | 0.456210 | 0.651632 | FALSE/300 |
| round2/single | no_start | 7 | 1 | 0.346748 | 0.654967 | FALSE/300 |
| round2/multimodal | selected | 7 | 0 | 0.000000 | NA | FALSE/500 |
| round2/multimodal | fixed | 7 | 0 | 0.059893 | 0.635837 | FALSE/500 |
| round2/multimodal | normal | 7 | 0 | 0.000000 | 0.629695 | FALSE/500 |
| round2/multimodal | framework | 7 | 1 | 0.255337 | 0.646367 | FALSE/500 |
| round2/multimodal | no_start | 7 | 0 | 0.000000 | NA | FALSE/500 |
| benchmark_snapshot/multimodal_framework | pruned | 7 | 1 | NA | 0.637768 | FALSE/NA |
| benchmark_snapshot/single_alpha_F0 | pruned_expression_only | NA | NA | NA | 0.593717 | NA/NA |
| benchmark_snapshot/expression_only_framework | pruned_expression_only | 7 | 1 | NA | 0.643559 | FALSE/NA |

The absorbing-prior hypothesis holds in both affected arms. In selected
multimodal, slab probability reaches exactly0 at outer99/inner1; all beta means
and SDs become exactly0 at outer99/inner2 and stay0 thereafter. Without the
supervised start, these events occur at outer54/inner3 and54/inner4. The prior
at extinction retains a nominal Laplace scale1, with zero slab weight; this is
not a numerical scale-zero diagnosis. Fixed Laplace never loses its .5 slab
weight, but has zero credible programs. Shared Normal avoids exact point-mass
absorption but ends with risk SD1.61e-7 and zero credible programs. Historical
per-factor Normal retains one credible program. None of Single's shared-prior
contrasts has a zero slab probability. No configuration is adopted.

Every multimodal arm hits the unchanged500-iteration output cap. Final raw
expression/methylation reconstruction changes are respectively:
selected .000446/.008929; fixed .001263/.015696; Normal .000456/.007998;
framework .001139/.006451; no-start .001139/.006375. Both exceed1e-4 in all
arms. Including intercepts or centering predictions does not remove molecular
movement (see `cap_diagnosis.csv`). Selected and no-start eta changes are0;
fixed/framework eta changes are .000424/.001163, Normal2.23e-10. Every arm's
last100 ELBO differences is positive: no observed ELBO oscillation supports
choosing damping at this point. This does not prove every coordinate is
monotone. No stopping criterion, iteration cap, prior default or benchmark changed.

The newly fitted framework training C .646367 differs from the frozen
framework's .637768. It is a current-code diagnostic refit, not a tagged-code
reproduction; no equality or superiority claim follows. This difference and
the framework's output cap remain explicit. No external score selects a remedy.

Next causally applicable rung: constrain the refitted shared Laplace slab
probability below at1/K (R2); R1 fixed-prior diagnostic already fails the
necessary full-training credible-program and convergence requirements, but
its required simulation/CV evaluation has not yet been performed. Those
unrun parts must not be described as completed. Continue the prescribed
ladder with a committed design and unchanged simulation/CV criteria.

Task3 scope update, Andrew's explicit answer (2026-10-15): **defer the real-data
head-to-head for review**. The prospective design above is preserved as a
proposal, not execution authorization. Installation and nine software tests
are preparation only. No DeSurv scientific fit or external score has run.
Andrew also explicitly approved omitting DeSurv from Task4's signed-Gaussian
panel, which will retain the five specified comparators.

Task2-1 result-checkpoint verification: fresh full suite **584 passed,0 failed**;
all ten benchmark SHA256 entries pass and the annotated/peeled tag IDs remain
unchanged. The result checkpoint adopts no configuration.

Task3 preparation/defer checkpoint: fresh full suite **593 passed,0 failed**,
including all nine author-adapter/driver tests. Both Longleaf scripts pass
bash syntax checks. No scientific DeSurv run was started.

## Task4 bounded rerun — fixed small-rank handling and full denominators

Before scientific execution: this design and all code/tests are committed.
The original calibrated bed, sample sizes, seed vector20261016:20261035,
noise .25/1, p50, beta3=3.3, and five original methods are unchanged. Andrew
explicitly approved omitting DeSurv from signed-Gaussian data; Task3's
separate real-data comparison is deferred. This bounded task is independent
of Task2's open remedy ladder and follows its completed diagnostic checkpoint.

Exact changes: preserve every method row even when another arm errors; use
ordinary Cox for an EBMF rank1 fit and explicitly mark lasso inapplicable;
for rank0, keep all rows in the20-dataset denominator, count recoveryFALSE,
C=NA. Check actual flash$n_factors before ldf(), which throws on valid rank0.
Record numbers of rank0/rank1/rank>=2 separately from failures. Cox-z selects
the two-step program; lasso uses that ranking and its own coefficient sign.
Lasso nonzeros are not posterior credible intervals, so their false-active
quantity staysNA. Method-specific warnings/errors stay explicit.

Arms: EBMF alone, EBMF->Cox, EBMF->lasso-Cox, PCA->Cox(K7), unchanged selected
canonical shared-point-Laplace supervised model(K7,no pruning,w_S1,warmup10).
There are40 independent datasets,200 method rows. Recovery and C definitions
are unchanged. Summaries keep n20 per method/noise, with denominators for
finite C/correlation explicitly labelled; errors count as not recovered and
cannot silently disappear. No new adoption rule or tuning occurs here.

Verification: six new tests were RED against the old runner, then GREEN in
an isolated draft. A first draft passed5/6; the remaining failure exposed
ldf() on rank0, corrected using the fitted object's actual retained count.
Tests cover rank0, rank1, isolated Cox and supervised errors, full denominator,
and distinguishing EBMF errors from zero factors. No package code changes.
Full-suite verification is required before committing/running.

Run: Rscript results/benchmark_sim/run_supervision_comparators.R 6.
Fresh output results/benchmark_sim/outputs/supervision_comparators/round2/;
40 per-dataset checkpoints, results.csv, summary.csv, ebmf_rank_counts.csv.
Fresh-folder guard prevents overwrites. Six serial PSOCK workers, one dataset
per task; no nested BLAS parallelism. All tables write before any final
nonzero exit for retained per-arm errors. No external cohorts or real data.

Task4 final pre-run verification: fresh full suite **599 passed,0 failed**,
including all six small-rank/error-isolation intent tests. The five-arm design
is committed before the40-dataset scientific rerun.

Task4 execution failure before any fit: committed dd11923, fresh599/0 suite.
The sandbox blocked PSOCK serverSocket port11504. Only design.csv exists;
zero datasets/fit results. Preserve the original round2 folder and its
execution_status.md. Retry the exact committed design/code using host socket
access, via source()+run_comparator_panel(6L,out_dir=.../round2_retry1).
The output folder and execution permission are the only changes. No arm,
seed, threshold, sample, prior or stopping criterion changes. This entry is
committed before retrying. No failed run is overwritten or discarded.

Task4 socket-retry pre-run verification: fresh **599 passed,0 failed**;
no scientific code changed from dd11923.

## Task2 attempt2 — complete fixed-shared-prior simulation/CV evaluation (R1)

**Before execution:** commit this entry, runner header, wrapper and tests before
any candidate simulation or training-CV fit. The ten-arm diagnosis is already
complete and separately committed589656a. This attempt changes exactly one
option from each selected model: beta_prior_fix=TRUE. No extra warmup, scale,
threshold, weight, cap or data change. Trace is disabled only as an observational
performance option whose numerical identity is already tested.

**Hypothesis:** preventing shared-prior refitting avoids exact absorbing slab
weight0, but may still fail credible-program/convergence or recovery criteria.
Full-TCGA fixed-prior diagnostics already fail the necessary credible-program
and output-convergence requirements. Complete the handoff's specified
simulation/CV evaluation; do not select or adopt R1 by ignoring those failures.

**Arms/test bed:** selected K7 Single (no pruning) and selected K7 multimodal
(ordinary unprotected pruning), fixed initial shared Laplace(pi=.5,scale1).
Original calibrated n400/validation120,p50 per modality,beta3=3.3,
seeds20261016:20261035,noise .25/1:80 new simulation fits. Compare false-active
counts to the80 unchanged saved K7 reference rows in Task1-3 results.csv:
alpha_F0 Single and10/1 multimodal. They are reused, never overwritten/refit.

**Training evaluation:** same144 TCGA patients,3000 genes,10000 imputed CpGs
from the approved unchanged cache, with IDs/order checked against frozen fit
metadata. Five event-stratified folds, seed42, shared across models/candidate
and reference. For Single, estimate column center/SD within training folds
and apply these to held-out patients; no validation values set the transform.
Multimodal retains its existing training-defined intercept/projection behavior.
Fit fixed-prior and each model's diagnostic-framework reference in each fold,
seed42:20 CV fits. References are exactly collapse_arm_controls(...,"framework"):
Single canonicalw_S1,per-programNormal,no supervised start; multimodal10/1
per-factorNormal,no supervised start. These are the already committed diagnosis
controls, not a new reference chosen after CV. The simulation Single reference
remains the handoff's alpha_F0 reference; training and simulation references are
explicitly labelled. No ICGC or other external cohort is loaded/scored.

**Pass/fail:** each model must recover>=16/20 at both noises, mean null
survival-active<=its specified simulation reference+.5, and all20 returned
fits converge under its own original rule. Intermediate pruning caps remain
separate. No median-K gate until K selection is adopted. The paired fivefold
mean deltaC must be>=minus one paired SE (sd of five differences/sqrt5).
An error/missing seed/fold or undefined C cannot pass by dropping rows. Full
TCGA credible/nonconstant/convergence requirements also remain necessary; R1
currently fails those. No external check or replacement proposal follows a
failure. Simulation validation C is descriptive, never used to select an arm.

**Outputs/run:** results/benchmark_sim/run_round2_fixed_prior.R6; fresh
results/benchmark_sim/outputs/round2/fixed_prior/.100 independent checkpoints,
results.csv,simulation_summary.csv,training_cv_summary.csv,
training_cv_with_benchmark.csv,folds.csv,simulation_references.csv. Frozen
benchmark rows are explicitly non-CV rows alongside new CV results; no frozen
training C is presented as a CV estimate. Every fit records errors/warnings,
iterations,convergence,pruning caps,risk SD and active-program count. RDS inputs/
fits remain ignored. Same full-prior diagnostic rows remain in diagnosis1.
Longleaf wrapper100 cells,8GB/cell,12h, serial fits; prepare/finalize modes
support independent cells, all100 required for aggregation. Never submit jobs.
Fresh-folder/MD5 guards prevent overwriting and a changed training cache.

**Software verification:** seven tests were RED0/7 before the new runner,
then GREEN7/7; the added invalid-mode guard passes (8/8 total). Tests verify
the exact80+20 design, balanced deterministic folds, training-only scaling,
one-option controls, no denominator deletion after failure, paired SE and
undefined-C handling, isolated exact error checkpoint and no overwrite.
Fresh full suite must pass before the preregistration commit and fits.


## Task4 completed five-comparator panel (2026-10-15)

Pre-run commits dd11923 and d353322; full suite599/0 before each. The first
sandbox launch failed before any dataset because a PSOCK socket could not
open; its design and execution_status.md remain in round2/. The unchanged
host-socket retry completed with exit0:40 datasets,200 method rows, zero
fitting errors. Outputs:
`results/benchmark_sim/outputs/supervision_comparators/round2_retry1/`.
DeSurv was omitted by Andrew's explicit decision; no arm or threshold changed.

| Method | Noise | Recovered/n | Validation C (finite n) | Null active/fit | Median K | ELBO-converged |
|---|---:|---:|---:|---:|---:|---|
| ebmf_alone | 0.25 | 0/20 | NA (0) | NA | 2 | not reported |
| ebmf_cox | 0.25 | 0/20 | 0.634182 (20) | 1.15 | 2 | not reported |
| ebmf_lasso_cox | 0.25 | 0/20 | 0.638527 (19) | NA | 2 | not reported |
| pca_cox | 0.25 | 0/20 | 0.687753 (20) | 4.25 | 7 | not reported |
| supervised | 0.25 | 17/20 | 0.763449 (20) | 1.1 | 7 | 19/20 |
| ebmf_alone | 1 | 0/20 | NA (0) | NA | 1.5 | not reported |
| ebmf_cox | 1 | 0/20 | 0.538504 (18) | 0.45 | 1.5 | not reported |
| ebmf_lasso_cox | 1 | 0/20 | 0.553546 (6) | NA | 1.5 | not reported |
| pca_cox | 1 | 0/20 | 0.596066 (20) | 3.3 | 7 | not reported |
| supervised | 1 | 19/20 | 0.765290 (19) | 0.25 | 7 | 19/20 |

EBMF rank counts: low noise0 zero-factor,0 one-factor,20 at least2; high
noise2 zero-factor,8 one-factor,10 at least2. No EBMF fitting error. Every
method/noise denominator is20, including constant-risk and inapplicable
cases. At high noise lasso is applicable to10 datasets and has finite C in6;
its undefined scores are retained. EBMF-alone has no fitted survival risk,
so its C is undefined by design. The summaries label every finite-C count.
Lasso's nonzero count is not Bayesian credibility; false-active isNA.

The supervised fit caps at300 iterations for seed20261026/noise .25 and
seed20261017/noise1; both recover the planted program. At seed20261018/noise1
it has no active program and undefined constant-risk C, counted as not
recovered. Consequently this panel does not establish universal convergence
or satisfy all project criteria. It is a comparator report, not adoption.
There are30 warning-bearing lasso rows: installed glmnet warns that a future
release will change the Cox tie default from Breslow to Efron. The exact
warning is retained per row; no package/default was changed to silence it.

This completes the bounded Task4 rerun, including rank/denominator reporting.
Task2's remedy ladder remains open and is the priority. No real-data or
external-cohort selection occurred in this simulation task.

Task2-2 software update before any scientific fit: the ninth test reproduced
CSV's all-NA integer column reloading as logical, which would reject valid
Longleaf checkpoints. Compare design-key values as character (includingNA),
still rejecting a changed seed. RED8/9, then GREEN9/9. No arm/data/criterion
change; this correction is part of the original pre-run implementation.

Task2-2 final pre-run verification: fresh full suite **608 passed,0 failed**,
including all nine new fixed-prior protocol tests and unchanged-default checks.
Longleaf wrapper syntax passes. No candidate simulation/CV fit ran before
this preregistration commit.

Task4 result-checkpoint verification: fresh full suite **608 passed,0 failed**.
The40-dataset retry exited0; all200 rows and40 cell checkpoints are preserved.
No model/default change follows the two supervised caps.

Task2-2 interim simulation result (all80 simulation cells complete; CV still
running, no adoption): errors0. Single low/high recovery17/20,19/20;
null-active1.00,.35 against reference1.15,.35 (all margins pass);
ELBO-converged19/20,20/20. Seed20261026/noise .25 caps at300.
Multimodal low/high recovery19/20,20/20; null-active .35,0 against
reference .75,.65 (all margins pass); output-converged20/20 at both noises.
Single therefore fails universal simulation convergence. Multimodal passes
the simulation gates but the full-TCGA fixed-prior diagnosis has zero credible
programs and caps, so this is not an adoption. Saved via the committed summary
function as outputs/round2/fixed_prior/simulation_summary_interim.csv.
All20 training-CV fit specifications remain unchanged and are completing.

Task2-2 data handling: folds.csv contains patient IDs and event indicators and
remains local/ignored, alongside controlled input/model files. Aggregate and
per-fit metrics contain no individual outcomes and may be committed. The
seed42 fold rule and a hash of the local fold manifest preserve the design
without committing a patient-level clinical manifest. This changes git
tracking only, not the running analysis or reference reuse.

Read-only interpretation audit during R1 CV: the frozen bed's survival eta
uses both observed modalities (including molecular residual variation).
Single's stored oracle is joint eta, while Single fits expression only;
its recovery target remains original raw F after per-gene standardization of
Y. Preserve all approved generators, thresholds and scores exactly. These
facts are recorded as framing limitations in the Claude handoff, not a reason
to change the bed or reinterpret a failed adoption criterion as passing.


Task2-2 framework CV collapse and Andrew's direction: completed multimodal
framework fold1 has beta means and second moments exactly0, no credible
program, risk SD0, C=NA, and hits500 iterations. Other completed framework
folds2/3 retain one prognostic program. Andrew explicitly requests continued
diagnosis of the constant prediction, not a tied-risk C=.5 substitution.
Keep the committed scorer and all five finite paired-fold requirements.
Compare this failed training fold with the successful full-TCGA framework
diagnosis and the other folds; do not tune on external data. The ongoing
R1 evaluation remains unchanged. No scoring amendment or new remedy is
adopted. This is an additional instability in the reference itself.

### Task2-3: unchanged framework fold1 trace (preregistered before fitting)

Hypothesis: the historical per-factor Normal model also loses survival signal
on a TCGA training fold because each variance MLE max(0,x²−1/A) can be zero.
Its variance is refitted every coordinate, so unlike a fixed shared point mass
this need not be absorbing. The saved fold1 has beta and second moments0;
its beta KL is0 from outer167 onward, while folds2/3 and full TCGA retain
signal. Read-only final-state pseudo-z magnitudes are .245,.503,.651,.910,
.019,.311,.524; all seven variance MLEs are0. Projection uncertainty
contributes at most .000381 of A at this saved state. This does not yet
explain the earlier trajectory or imply an external-data defect.

Exact change: enable observational trace only for R1 cell84, seed42/fold1,
same frozen TCGA matrices, per-factor Normal, no supervised initialization,
K7, ordinary pruning, original500 cap/output rule. Use the existing full-TCGA
framework trace as a descriptive comparison; do not refit or select on it.
Runner: results/multimodal_real/run_framework_fold_trace.R.
Output: results/multimodal_real/outputs/round2/framework_fold_trace/.
No external data or C=.5 substitution. Andrew requests continued collapse
diagnosis and leaves C=NA. No criterion, prior, default or test bed changes.

Pass/fail for diagnostic validity: all saved posterior fields and history
must be identical to original cell84. Report every inner-sweep pseudo-z and
zero/nonzero posterior sequence, including revival, persistent extinction,
and cap. A mismatch blocks interpretation. This diagnostic cannot adopt a
configuration and does not count as a successful remedy. Fresh outputs only.
Expected serial runtime about1hour,8GB; no HPC submission. Code/report
checkpoint and full R suite precede this diagnostic fit.

Task2-3 preregistration verification: new focused checks3/3 passed after
RED0/3; full suite **611 passed,0 failed** before the diagnostic commit.
The trace has not yet been fitted at this checkpoint.


### Auxiliary mathematical audit during Task2-3 (no real-data fit or remedy)

Canonicalization maps (L,F,beta) to (aL,F/a,a beta), factor by factor.
It transforms L/F priors but has no beta-prior argument; the shared prior
g_beta is not transformed by the caller. Preserving fitted products alone
therefore does not establish preservation of KL(q_beta||g_beta) under a
shared fixed prior. A deterministic two-program software calculation uses
scales2.581989/1.154701 and a shared N(0,1) prior: eta changes0, while the
correct Normal KL changes1.629438/1.209438 to3.627545/1.098930.
The current carried beta KL is not recomputed by canonicalization.

This is a prior/parameterization consistency question for review, distinct
from the historical per-factor Normal CV collapse; that update re-estimates
its own variance at each coordinate. It does not demonstrate that scale
handling causes the observed zero beta. In the saved selected/no-start
traces the per-sweep scales are already close to1 at prior extinction
(selected outer99 range .998530–1.000556; no-start outer54 .998964–1.003347),
although first-outer scales are much larger. No prior, canonicalization,
objective, default or candidate is changed on this evidence. Any proposed
shared-prior configuration needs this formulation checked in addition to
the numerical adoption criteria; product-invariance tests alone are
insufficient mathematical validation. No new real-data fit or score follows.


### Saved framework loading-pattern comparison (descriptive, no new fit)

All comparisons use existing fits and the same frozen feature order. For
benchmark program4 (the frozen fit's sole credible program), the largest
absolute loading correlations by block are:

| Saved fit | Expression r | Methylation r | Combined r | Credible programs | Training C | Fold-validation C |
|---|---:|---:|---:|---:|---:|---:|
| Frozen benchmark full TCGA | 1 | 1 | 1 | 1 | .637768 | NA |
| Current full-TCGA framework diagnosis | .915133 | .898948 | .895973 | 1 | .646367 | NA |
| Current framework CV fold1 | .787112 | .685382 | .752163 | 0 | NA | NA |

Each strongest match above is program4. Relative to the current full-TCGA
program4, fold1's correlations are .938139 (expression), .846790
(methylation), .933739 (combined), yet its beta is exactly0. Thus the
loading pattern remains recognizable; this is not removal of the program
by pruning (fold1 has K7 and no pruning refit). This comparison is a
description of saved fits, not a cross-validation estimate using a
full-training loading matrix, an external selection criterion, or proof of
benchmark reproducibility. No fit, score convention or model is changed.


### Task2-2 final fixed-prior evaluation: not adopted

The100-fit local evaluation exited0:80 simulations and20 TCGA CV fits,
with0 fitting errors and0 nonempty warning records. All100 design keys,
the unchanged training-cache MD5 and the local fold-manifest SHA256 were
verified. Output: results/benchmark_sim/outputs/round2/fixed_prior/.

| Model | Noise | Recovered | Null active / reference | Converged |
|---|---:|---:|---:|---:|
| Single | .25 | 17/20 | 1.00 / 1.15 | 19/20 |
| Single | 1 | 19/20 | .35 / .35 | 20/20 |
| Multimodal | .25 | 19/20 | .35 / .75 | 20/20 |
| Multimodal | 1 | 20/20 | 0 / .65 | 20/20 |

Recovery and false-positive margins pass throughout. Single's low-noise
seed20261026 hits300 iterations, so universal simulation convergence fails.
These are unchanged K7 beds; K selection has not been adopted.

| Source/model | Fixed-prior CV mean C | Framework CV mean C | Delta / paired SE | Converged candidate/reference | Frozen training C |
|---|---:|---:|---:|---:|---:|
| Round2 Single (5 pairs) | .608225 | .623985 | -.015761 / .006968 | 5/5, 5/5 | NA |
| Round2 multimodal (5 specified pairs) | NA | NA | NA / NA | 0/5, 0/5 | NA |
| Frozen multimodal benchmark, full TCGA | NA | NA | NA | NA | .637768 |
| Frozen alpha_F0 expression-only benchmark, full TCGA | NA | NA | NA | NA | .593717 |
| Frozen expression-only framework, full TCGA | NA | NA | NA | NA | .643559 |

Single fails the one-paired-SE CV rule. Multimodal's reference fold1 has
C=NA/risk SD0; no finite five-pair difference exists, so the rule cannot
pass. The candidate's own five C values are all finite, mean .599978,
but none of its five fits has a credible program, and all cap at500.
All five reference fits also cap at500; reference credible-program counts
by fold are0,1,1,0,0. Its finite fold C values2–5 are .541502,.586047,
.485437,.359223. No conditional four-fold mean is substituted.

The frozen training rows above are labelled non-CV references, not
comparisons against a refitted benchmark. Earlier full-TCGA fixed-prior
diagnosis likewise has0 credible programs and a500 cap. Thus R1 fails
necessary real-data criteria despite the multimodal simulation passes.
No fixed prior, K, default or benchmark change is adopted. No external
cohort is scored. The unchanged-fold trace continues, and Andrew approves
the separate prospective profiled-curvature experiment after exact trace
reproduction; its implementation is still private and no scientific
profiled-curvature fit has run.

Task2-2 final checkpoint verification: fresh full suite **611 passed,0 failed**.
All100 fit records are retained; clinical folds.csv and all RDS/data remain
ignored. Only aggregate/per-fit metrics and hash/package provenance are staged.


### Task2-3 result: exact reproduction; loss at outer167/inner1

The unchanged-fold trace exited0 in2142.164 seconds, with no warnings.
All13 checks in identity.csv are TRUE: every saved posterior/prior field
and the complete500-iteration history equal original R1 cell84 exactly.
No numerical fit changed when observation was enabled.

Of2500 inner sweeps, the first all-zero posterior occurs at sweep831,
outer167/inner1. The maximum absolute pseudo-z falls from1.368511,
1.262920 and1.128647 in outer166/inner3–5 to .966195 at the next sweep.
Programs3 and4 are the last two nondegenerate coefficients. All1670
sweeps from831 through2500 have zero means/SDs; there is no later revival.
Earlier revival is real: programs1/3/4/5/6 revive7/8/2/1/2 times,
respectively. The existing full-TCGA framework trace never loses all
coefficients, and its program4 is nondegenerate throughout2500 sweeps.
This supports the per-coordinate variance-boundary diagnosis, not shared
point-mass absorption, a scoring defect or removal by pruning.

The repeated fold retains K7, has no credible program, risk SD0/C=NA,
and reaches the unchanged500 cap. Its comparison_with_benchmark.csv retains
all frozen benchmark training rows. The frozen multimodal benchmark
training C remains .637768 (one credible program); this is not a CV row.

Outputs: results/multimodal_real/outputs/round2/framework_fold_trace/:
identity.csv, fold1_trace.csv, fold1_boundaries.csv,
full_training_boundaries.csv, comparison_with_benchmark.csv,
training_source.txt and warnings.txt. The fit RDS remains ignored.

Andrew's explicit approval of the separate profiled-curvature test was
conditional on this reproduction, and that prerequisite now passes.
No remedy has been fitted at this checkpoint. Source edits still wait for
the fresh full suite/result commit. The prospective test uses the existing
framework priors/initialization and changes only the working-intercept
profile flag; it will not combine the unresolved shared-prior scale question.

Task2-3 result-checkpoint verification: fresh full suite **611 passed,0 failed**.
All13 exact reproduction checks pass; only safe factor-level diagnostic CSVs
and hash/status records are committed, with the model RDS kept ignored.

### Task2-4: approved profiled Cox working intercept (preregistered before fitting)

**Authorization/prerequisite.** Andrew explicitly approves this additional
remedy outside the prior ladder. Task2-3 must first reproduce R1 cell84 exactly;
its13 posterior/history checks pass, including all500 history rows.
The collapse at outer167/inner1 coincides with every per-factor pseudo-z
falling below1, and occurs with the loading pattern retained. Read-only
saved-state calculations show58% constant-component curvature for the matched
prognostic program. The diagonal approximation penalizes a direction along
which the actual Cox partial likelihood is invariant.

**One change.** Set profile_cox_intercept=TRUE on each model's existing
framework controls, otherwise unchanged: per-program/per-factor Normal beta,
no supervised initialization, point-Laplace molecular priors, w_S1 in Single,
K7, Single unpruned, multimodal ordinary unprotected pruning. This does not
combine the separate shared-prior/canonicalization issue with a curvature
remedy. Unsupported shared-prior combinations fail explicitly.

Profile the arbitrary constant in the fixed quadratic, giving
Hc=diag(w)-ww'/sum(w), hc=u+Hc eta0. An exactly equivalent implementation
weighted-centers only the survival design under frozen w, including its
loading-uncertainty moments, and centers the old predictor in its working
response. Genomic residuals/precisions, raw prediction, canonicalization,
reported ELBOs/pruning scores and original convergence monitors/rules/caps stay
unchanged. The multimodal fixed-working monitor reports its actual profiled
quadratic; Single retains its original ELBO monitor and rule. Refreshing either
working approximation does not guarantee monotonicity of the original monitor.
Zero total curvature preserves the zero-information path; no tuned constant
or new prior, observed-data change, default or threshold is introduced.
The derivation and limitations are in
docs/plans/Round2_profiled_curvature_proposal.md.

**Simulation and training arms.**80 new simulations on the original beds:
20 seeds20261016:20261035, noise .25/1, both models, n400/validation120,
p50 per modality, beta3=3.3. Reuse the original80 K7 simulation reference
rows, with unchanged false-positive bounds. TCGA uses the approved144-patient,
3000-gene/10000-CpG imputed training cache and original event-stratified fivefold
assignments(seed42). Single scaling is estimated within each training fold.
Screening/imputation stays fixed: this is conditional CV, not a nested raw
pipeline. Fit10 new candidate CV fits and both full-TCGA candidates. Reuse
exactly10 unchanged R1 framework CV references after checking the complete100
source keys, cache/fold hashes, per-cell keys/errors and source checksums.
Total102 records,92 new fits; no external cohort for selection.

**Unchanged pass/fail.** Recovery>=16/20 at each noise; null-active mean
<=original model/noise reference+.5; every simulation fit converged under
its original rule (ELBO for Single, output change for multimodal). The
retained-rank gate is not adopted. Training CV requires all5 finite matched
pairs and mean delta>=-sd(delta)/sqrt5, without dropped/error/constant folds.
The existing multimodal framework fold1 remains C=NA, so that paired
criterion cannot pass unchanged; success on another criterion cannot conceal
it. Report a candidate's own mean only when all its5 folds are finite,
separately from the paired comparison and labelled finite counts.

Full-training no-collapse requires>=1 credible program (SD>0 and
|Ebeta|>=1.96SD) AND nonconstant risk. All caps remain failures against the
convergence goal. Every real-data table retains labelled frozen benchmark
rows. This is not automatic adoption: only all-criterion evidence can support
a proposal pending Claude review and Andrew approval, with the prescribed
one-time external report only for a final configuration selected without
external values. No external tuning, rule amendment, C=.5 substitution or
new reference is authorized by this experiment.

**Software and execution.** The new behavior first fails8 model/math and6
protocol checks on the old project. Private implementations pass8/8,
6/6 and9/9 existing R1 checks. Dense independently computed covariance
quadratics agree with actual precision, linear terms and the fixed-working
monitor; raw genomic residuals are preserved. Default/off objects match
the previous implementation; existing base-branch identity tests also remain
required. Before scientific execution, apply the opt-in code, register tests,
run the full suite, record its actual count and commit this entry/code.
Do not infer the full count from private tests.

Runner: results/benchmark_sim/run_round2_profiled_cox.R6.
Fresh output: results/benchmark_sim/outputs/round2/profiled/.
Longleaf wrapper: round2_profiled_cox.slurm,92 new serial tasks,8GB each,
12hours; prepare freezes the10 references, finalize aggregates102 records.
No jobs are submitted. Local execution uses six serial socket workers and
the recorded R1 BLAS environment; package/runtime provenance is preserved.
Patient-level folds.csv and all data/RDS stay ignored; only safe aggregate/
per-fit metrics, hashes and provenance may be committed. Any error/warning
is retained per fit, never silently retried or overwritten.

Task2-4 preregistration verification: focused23/23 (14 new plus9 existing)
and fresh full suite **625 passed,0 failed**. Default/off fits match
protected-base fitter definitions. The SLURM wrapper passes bash syntax
checking. No scientific profiled-curvature fit has run before this commit.


### Additional software invariant found during Task2-4 (no model change)

A deterministic zero-pseudo-observation helper call returns a refitted
point-mass molecular prior while returning posterior moments computed under
the previous nondegenerate prior. For F with Normal prior variance1, the
returned second moment is .1851852 and stored KL .4357921, while the new
prior is delta0. For the actual point-Laplace loading family (pi=.5,rate1),
the returned second moment is .04679413 and KL .2366599, again with a
returned delta0 prior. A shared-score example similarly returns mean
.3246915, second moment .4220161 and KL .1185925 with a new delta0 prior.

These q distributions are not absolutely continuous with respect to the
returned point prior; their KL to that returned prior is infinite. The
finite stored KL instead refers to the prior used before hyperparameter
refitting. This is a concrete prior/posterior consistency gap, predating
Round2, that product-invariance/default-identity tests do not check.
It affects interpretation of the reported ELBO/pruning criterion as the
variational objective of the returned state. It is not established as the
cause of the real-data collapse or Single oscillation. The completed full-TCGA
selected/framework states do not have a point-mass/positive-second-moment
support conflict on inspection; finite prior-lag discrepancies require
further observation. No benchmark result, default, helper or running script
is changed on this finding; record it for methodological review.

Single's completed simulation seed20261026/noise .25 shows oscillation in
both fixed-prior and profiled-framework fits. Over the last50 iterations,
relative ELBO increments are negative16 and17 times, respectively; ranges
are [-.000197505,.000128715] and[-.000180143,.000116265].
Final output changes .143724 and .105693 also fail the output tolerance.
This is not simply a near-convergence cap; the original ELBO rule remains
the criterion. No extra iterations, damping or rule change is applied.


### Pause checkpoint preparation — prior audit prepared, not executed

Andrew requested a pause after the running Task2-4 panel finishes and its
results are committed. No additional scientific fit will start before that
pause. The next diagnostic is prepared only, not an executed attempt and not
a model remedy; it must receive its prospective attempt-log entry and commit
before any execution after resumption.

`round2_prior_consistency_helpers.R` computes posterior slab probabilities
and E[|F|] from each existing point-Laplace quadratic (A,B), then evaluates
E_q[log g] under the prior used for the update and the returned refitted prior.
The difference yields the corresponding F-prior KL correction; a support
conflict is retained as infinity. The observer returns the original update
object unchanged and records aggregate iteration/factor diagnostics. It uses
base R probability functions and the existing truncated-Normal moment helper;
no dependency, fitted prior or posterior is replaced.

`run_round2_prior_consistency_trace.R` is prepared to repeat only the already
completed Single simulation seed20261026/noise .25 with identical controls.
It must reproduce the entire saved fit exactly. It records every warning,
partial loading trace on error, prior-density differences and a diagnostic
F-prior-adjusted monitor. That monitor corrects only the measured F-prior
contribution; it is not a corrected full objective or an adoption criterion.
The analytic/native slab-probability tolerance is a software check only.
Six focused tests cover numerical quadrature, units, incompatible support,
a fixed-q prior-density identity, the no-spike boundary and unchanged update
objects. These tests do not establish a remedy for any scientific failure.

Checkpoint software verification (2026-10-15): fresh
`Rscript tests/run_tests.R` completed with **631 passed,0 failed**.
This includes the six prior-consistency fixtures. The prepared observational
simulation rerun has not been executed; no new remedy was tested by those fixtures.


## Task2-4 completed: profiled Cox curvature is not adopted (2026-10-15)

The preregistered panel finished with exit0:102 records,92 new fits and10
unchanged cached framework CV references, zero execution errors and zero
warning-bearing records. All original controls, criteria, fixed simulation
beds, TCGA cache and folds were retained. Output directory:
`results/benchmark_sim/outputs/round2/profiled/`.

| Model | Noise | Recovered /20 | False active | Reference false active | Converged /20 |
|---|---:|---:|---:|---:|---:|
| Single | .25 | 18 | 1.35 | 1.15 | 19 |
| Single | 1 | 20 | .40 | .35 | 20 |
| Multimodal | .25 | 16 | .80 | .75 | 20 |
| Multimodal | 1 | 19 | .40 | .65 | 20 |

Every recovery and false-positive criterion passes. Single seed20261026 at
noise .25 reaches300 iterations; its oscillation is documented above. All40
multimodal simulations converge. These are the approved recovery/false-positive
criteria, not proof that protected K selection is resolved.

| Source/model | Mean candidate CV C | Reference CV C | Paired delta | Paired SE | CV converged | Frozen training C |
|---|---:|---:|---:|---:|---:|---:|
| Round2 Single, profiled | .624776 | .623985 | .000791 | .000791 | 5/5 | — |
| Round2 multimodal, profiled | .563678 | NA | NA | NA | 0/5 | — |
| Frozen multimodal framework | — | — | — | — | — | .637768 |
| Frozen single alpha_F=0, matching TCGA comparison | — | — | — | — | — | .593717 |
| Frozen expression-only framework | — | — | — | — | — | .643559 |

Single passes its five-pair one-SE comparison. Multimodal has five finite
candidate scores, but only four finite original-reference scores: the
predeclared paired comparison remains undefined and fails eligibility. The
candidate's own mean uses all five folds; no undefined reference is replaced
or dropped. All five profiled multimodal folds now have nonconstant risk,
including formerly collapsed fold1; their credible-program counts are
0,1,0,0,0. Every multimodal candidate CV fit reaches500 iterations, as did
all five cached reference fits. All five Single candidate and reference CV
fits converge. Frozen training C values are labelled historical training
quantities, not CV estimates or external scores.

| Source/model | Training C | Credible programs | Nonconstant risk | Converged | Iterations | K_final |
|---|---:|---:|---|---|---:|---:|
| Round2 Single, profiled | .651632 | 3 | yes | no | 300 | 7 |
| Round2 multimodal, profiled | .639172 | 1 | yes | no | 500 | 7 |
| Frozen multimodal framework | .637768 | 1 | — | no | — | 7 |
| Frozen single alpha_F=0, matching TCGA comparison | .593717 | — | — | — | — | — |
| Frozen expression-only framework | .643559 | 1 | — | no | — | 7 |

Both full-training candidates pass the defined no-collapse check and fail
convergence. No full-training factor is deleted. Eight of92 new fits cap:
one Single simulation, the full Single fit, five multimodal CV fits and the
full multimodal fit; another five caps belong to the unchanged cached
multimodal references. No intermediate pruning cap is reported. Detailed
per-fit C, activity counts, risk SD, rank, caps and warnings are in
`results.csv`, `training_cv_with_benchmark.csv` and
`training_full_with_benchmark.csv`; the latter two retain frozen rows.

**Decision: do not adopt.** Profiling removes exact constant-risk collapse
in this candidate's five training-CV folds, but does not deliver stable,
credibly prognostic programs throughout CV or universal convergence. The
original paired-reference issue remains open. No external cohort is scored,
no criterion/default is changed, and the benchmark remains the current best
state. Further Task2 work is paused by Andrew, not resolved or ladder-exhausted.
The prior/posterior KL and common-beta-prior scale questions remain separate
review items, not established explanations for these results.

Final checkpoint verification: full R suite **631 passed,0 failed** (including
six new prior-consistency fixtures); prepared observer syntax parses. The
observer itself has not been run and requires prospective registration and
commit before execution. All ten benchmark SHA256 checks pass after the
panel; tag object0fea14903a75606e495267f2f35712eac10d9fc6 still peels to
a352f17ca9cff523caa80756aa2af085e62551e1. The writeup.qmd/html/pdf files
have no diff against that tag. Patient-level folds.csv and RDS/data are
ignored and excluded from the commit. No push, protected-branch edit or
merge has occurred.


## Task5-QC1: full-probe TCGA methylation variability (prospective)

Andrew approved the revised order before implementation. This descriptive
attempt addresses the10/2 per-probe mean/SD histogram request, before any
feature-count model fits. It changes no model, filter, preprocessing,
imputation, criterion or default and uses no survival outcomes or external
scores. No full-probe mean/SD analysis has run at this preregistration.

Inputs: raw TCGA `tcga_meth_filtered_sex.rds`$meth (305492 probes x150 samples),
existing TCGA-only training cache (144 patients x10000 screened probes), and
the frozen joint fit's training_spec feature/patient manifest as the screened reference.
Patient matching is by cache IDs, preserving its144-patient order. The raw
universe is all305492 available TCGA probes, not a newly approved shared
TCGA/ICGC universe. No ICGC values, outcomes or scores are used. Input MD5
and a non-identifying patient-ID hash are retained; raw/cached inputs and the
snapshot must not change. Cross-reactive/SNP-filter provenance is currently
unverified; the run applies no additional QC.

For each probe on observed raw beta values: count observed/missing, mean,
sample variance/SD using n_observed-1, arcsine-scale SD after asin(2*beta-1),
and proportions <.1, >.9 and outside[.1,.9]. Probes with n<2 have undefined
SD; n=0 also has undefined mean/filter membership. They remain in inventory
and aggregate denominators. No missing value is filled. Constants stay zero.
The unchanged cached screened imputation is summarized separately, without
rerunning imputation.

Descriptive groups: all available; mean filter .1<=mean<=.9; proportion
filter retaining probes with neither extreme fraction>=.95; both; mean-only;
proportion-only; current screened set. The two extreme fractions are tested
separately, not combined. Unknown membership is reported. Screened probe IDs
must match the frozen benchmark exactly as sets. There is no new SD cutoff.

Figures: per-probe SD histograms on beta/arcsine scales,80 bins, separate axes
for each group/scale; mean-beta versus SD binned density. Histograms display
all/mean/proportion/screened groups; finite/missing SD counts remain explicit.
They do not estimate a survival benefit or choose features from outcomes.
Probe inventory RDS stays ignored/local. Commit only aggregate tables,
provenance and plots. Fresh output only:
`results/multimodal_real/outputs/feature_sensitivity/cpg_qc/`.

Software validity: known means/sample variances and transform moments match
independent calculations; block size does not change values; missing/constant
and boundary cases are correct; invalid range/identities fail; filter counts
and histogram bins conserve finite observations; the driver preserves input
hashes and refuses overwrite. All warnings and errors are retained. A failure
is reported, not silently retried or used to change a cutoff. This attempt
has no model-adoption pass rule; it cannot satisfy supervision/pruning goals.

Runner: `Rscript results/multimodal_real/summarize_cpg_variability.R run`.
Serial local run, minutes, several-GB input peak;10k-row summaries bound the
extra working memory. No jobs or parallel workers submitted. Preregister,
run the full suite, record the actual count and commit BEFORE execution.

Software preflight clarification before any descriptive run: the frozen
loading matrix has no row names. Its training_spec stores feature and patient
identities; the driver verifies those fields, not row positions. The
end-to-end fixture was tightened to this real archive schema (RED), then the
driver was corrected (GREEN). No result, test bed or statistical rule changes.

Task5-QC1 preregistration verification: focused11/11 pass after real-archive
manifest RED/GREEN check. Fresh full suite **642 passed,0 failed**. Frozen
feature IDs and patient order match the unchanged cache. No full-probe
variability result has been computed before this commit.


## Task5-QC1 result: full-probe methylation variability completed

Preregistered in dc69f6d; serial run completed with exit0, no errors and zero
warning records. All305492 available TCGA probes are described on the unchanged
144 training patients. Frozen screened probe IDs and training patient order
match the cache. No imputation, survival fit, filter/default adoption or
external scoring occurred. Outputs:
`results/multimodal_real/outputs/feature_sensitivity/cpg_qc/`.

| Descriptive group | Probes | Overlap with frozen screened10k |
|---|---:|---:|
| All available TCGA |305492|10000|
| Mean beta in[.1,.9] |182996|9990|
| Neither observed extreme fraction>=.95 |221222|10000|
| Both filters |182996|9990|
| Mean-only |0|0|
| Proportion-only |38226|10|
| Frozen/current screened reference |10000|10000|

Every probe has defined sample SD; no probe is constant. There are29414
missing raw entries, retained as missing in observed-value summaries. Of
these,2130 occur in the screened reference; the unchanged cached imputation
contains none. Median beta SD is .060008 across all probes, .179288 in raw
screened probes and .179186 in the imputed screened cache. Median arcsine SD
is .153757/.417231/.416941, respectively. Median probe mean beta is
.461005/.350031/.350125. These are data descriptions, not survival results.
The current10k screen is strongly enriched for variability; that does not
establish the best feature universe. The mean filter would discard10 of the
currently selected probes and is not adopted. Upstream cross-reactive/SNP
filtering remains unverified; this is an explicit provenance gap.

The requested SD histograms and mean-versus-SD binned plot were visually
checked. No SD threshold was selected. An independent check on21 probes
(including the most-missing probe) agrees with base-R mean/SD to maxima
4.44e-16/4.16e-17, and transformed SD to5.55e-17. Cached values at observed
raw entries agree exactly. Input hashes remain unchanged; zero warnings.
The full probe inventory RDS remains local/ignored; only aggregate tables,
plots and non-identifying hashes/provenance are committed.

The advisor progress PDF is updated to six pages with these results and the
SD figure; changed pages were visually checked. The original write-up is
unchanged. Fresh post-run full suite **642 passed,0 failed**. No new
performance, biological or supervision claim is made. Next: the fixed-K7
expression feature-count study and its training-only factor-agreement curve.


## Task5 expression preflight and agreement software (no feature fit)

The D4 training universe has17900 shared genes on273 patients. Its frozen
baseline gene order is reproduced exactly, and its standardized matrix is
exactly the corresponding subset of the all-gene matrix. Existing per-cohort
combined ranking reaches shared counts2064/5001/10000/17900 at screening
sizes3000/6607/11602/17900. No ranking rule or tie-breaking changes. The nominal
5000 target is crossed at5001 and will be plotted at its actual count.
Safe design tables are in
`results/benchmark_sim/outputs/feature_expansion/design_preflight/`.

The frozen D4 activity matrix has rank4 of7 and exact zero columns1,2,4.
Thus original Task5's inverse(L-transpose-L) continuation is undefined. A
cold-start primary curve at all four counts is proposed, preserving fixed K7,
alpha_F=0 benchmark settings, five training folds and all original criteria;
warm continuation would wait for a rank-deficient policy. This decision is
pending Andrew; no feature-count fit has started.

`code/feature_sensitivity.R` supplies only independent preparation/agreement
utilities. `feature_count_manifest()` consumes cohort-specific log-scale
training matrices and the existing rankings, reports the reached shared counts
and unavailable targets, and returns nested gene sets. It changes neither
preprocessing nor model fitting. `feature_cor()` calculates centered Pearson
correlations, preserving NA for constant columns. `feature_program_agreement()`
uses the existing greedy non-reused activity matching pattern, aligns loading
rows on common gene identities, compares proper normalized projection/risk
scores and reports beta times projection SD. This avoids raw coefficient-unit
comparisons. Matching matrices, unmatched slots and continuous best-match
margins stay visible; no match-quality threshold is selected.

Inputs to agreement must represent the same ordered patients. Named row-order
mismatches fail, and the future runner must retain the ordered identity manifest
for unnamed matrices. Constant factors remain unmatched, rather than forcing
correlation0 or a fabricated match. Signed correlations remain signed; frozen
biological labels are not automatically assigned. Existing base-R operations,
preprocessing and prediction avoid dependencies. Nine focused checks pass
from RED to GREEN, including scale/permutation/common-feature/constant and
patient-order cases. Full-suite count is recorded after completion below.

The matched expression-only multimodal framework itself retains K7 in the
frozen results. The pooled Single classified4 versus joint7 comparison also
changes training population and model settings; it does not isolate a
methylation effect. A retained K7 can be appropriate if its nonredundant
molecular/prognostic roles and stability justify it. Fewer factors alone are
not an adoption criterion for the real-data study.

Feature preparation checkpoint: fresh full suite **651 passed,0 failed**,
including nine new agreement/manifest checks. No feature-performance fit has
run. The original benchmark/write-up are unchanged; the initialization
decision remains pending.


## Task2-5: unchanged Single loading-prior consistency audit (prospective)

This is the bounded observational diagnostic already included in the approved
priority plan, executed while the independent feature-curve initialization
decision is pending. It is not a remedy, new test bed, parameter selection or
claim that debugging replaces the feature deliverables. The report/CpG tasks
are complete; no feature-performance fit is authorized by this entry.

Hypothesis: the F update returns moments under its incoming point-Laplace
prior, then refits and returns a different prior, while retaining KL against
the incoming prior. Quantify that discrepancy during one actual completed
simulation, beyond the existing synthetic counterexample. One observational
change only: a process-local wrapper records old/new priors and moments and
returns the original loading-update object unchanged.

Repeat ONLY Task2-4 cell021: seed20261026, noise .25,400 training/120 simulated
validation subjects,p50,K7; framework per-program Normal beta, no supervised
start, profile_cox_intercept TRUE, w_S1, no pruning, the original300-iteration
ELBO rule/cap. No weight, prior, initialization, bed, convergence criterion or
scorer is changed. No external or patient-level real-data input is read.
The saved original fit is
`results/benchmark_sim/outputs/round2/profiled/cells/cell_021.rds`.
Its original result is recoveryTRUE and convergenceFALSE at300; neither can
be reclassified by this diagnostic. The simulated validation outcomes are
not used to select or repair anything.

Per loading block record incoming/returned pi and rate, reconstructed slab
mass and absolute moment, native returned slab mass and second-moment sum,
old/new expected log prior, cached KL and the F-prior KL correction. Analytic
support conflicts and returned point priors with strictly positive native
second moments are counted separately; a tiny probability discrepancy alone
does not establish the native posterior support. The monitor adjustment is
only the measured F-prior contribution, not a corrected full ELBO.

Validity requires exact identity of the ENTIRE saved fit, plus every named
field and complete history, and max analytic/native slab-probability gap
<=1e-7 (existing software tolerance, not an adoption threshold). Preserve all
warnings, partial rows and errors; any mismatch blocks interpretation. Keep
all finite/nonfinite adjustments and support conflicts; no missing cases,
extra iterations, retry, damping or new prior update is allowed. Report the
original convergence failure unchanged. No adoption is possible in this audit.

Runner: `Rscript results/benchmark_sim/run_round2_prior_consistency_trace.R`.
Fresh output:
`results/benchmark_sim/outputs/round2/prior_consistency_seed26/`.
Local serial, minutes,1GB; no HPC job. Safe CSVs/provenance only may be
committed; fit RDS remains ignored. The already prepared observer is extended
only with native moment telemetry, tested RED/GREEN before this commit.
Run the full suite and record its actual count before committing this entry;
no audit fit runs until the commit exists.

Task2-5 preregistration verification: strengthened native-telemetry fixture
RED5/6 then GREEN6/6; fresh full suite **651 passed,0 failed**. No audit fit
has run before this commit.


Task2-5 result: exit0,13.544s, zero warnings; exact entire-fit identity and
18/18 component checks. Analytic/native probability error6.66e-16. Three
support conflicts independently confirmed by positive native EF2 with a
returned point prior, at iterations2/3/8; none in the last50. The F-prior-only
adjustment retains17 negative late increments; the fit remains unconverged
at300. No full-objective correction or model remedy is established. Safe
aggregate outputs are in `outputs/round2/prior_consistency_seed26/`; local fit
RDS is excluded. Fresh suite651/0. Andrew subsequently authorized cold-primary
feature sensitivity and requested HTML progress output; no feature fit yet.

## Task5-F1: fixed-K7 cold-primary expression feature sensitivity (prospective)

Andrew approved the revised advisor-first order and authorized judgment to
proceed after the prescribed frozen warm inverse proved singular. One
experimental axis: feature count. Four nested sets use the unchanged D4
combined mean/variance cohort screening and intersection:2064/5001/10000/17900
shared genes (per-cohort cutoffs3000/6607/11602/17900). The5000 target crosses
at5001; no arbitrary trim. Full common measured universe is17900. Same273
TCGA+CPTAC training patients, same existing log flags and per-platform z
standardization, no external loading. Existing constant columns are retained.

All24 fits use K7, cold SVD, alpha_F0, alpha/configured maximum iterations and
tolerance from globals.yml, point-exponential L/F, per-program Normal beta,
no canonicalization, no supervised start, no cohort indicator, no pruning.
This is sensitivity of the benchmark formulation, not a supervised remedy.
Warm continuation is deferred; no model default or adoption rule changes.
Five event-stratified folds from create_stratified_folds(status,5,seed42) are
identical across counts, plus one full fit per count; seed42 for each fit.
Selection/scaling are fixed before patient CV, so label it CONDITIONAL CV,
not a fully nested preprocessing assessment. CV sign orientation uses only
training outcomes under the existing CV convention; full fits retain the
benchmark sign correction. No feature count or weight is adopted from this
study, and no external cohort is rescored.

Record all24 rows: held-out C and Cox partial log-likelihood, training C for
full fits, risk SD, operational shrinkage categories using existing configured
thresholds, convergence/caps, iterations, warnings/errors, runtime, R allocation
high-water MB (not OS RSS), and constant-feature counts. Constant risk C=NA,
never substituted0.5; errors stay visible. Global ELBOs across different p
are not comparable. Operational activity is not posterior credibility.

After full fits, use the tested agreement utility against the fresh2064-gene
fit and the preceding feature count: non-reused activity matching, common-gene
signed loading correlations, normalized YF projection/risk correlations,
beta per projection SD and best-match margins; undefined/dead entries remain
NA/unmatched. Save full matrices locally and safe aggregate pairs/tables.
Compare the fresh2064 full fit directly with frozen D4; discrepancies must
be explained before biological claims. Label frozen benchmark training rows
alongside new C/logPL tables; the frozen external mean0.627 remains separate
historical evidence, not a new CV estimate. No external scores select settings.

Runner: `Rscript results/benchmark_sim/run_feature_expansion.R prepare`, then
`... all` (serial) or a single cell1-24. Prepare freezes ordered local patients,
genes/folds/controls and input hashes; that RDS stays ignored. Aggregate design
and per-fit CSVs are safe; no patient-level predictions or outcomes committed.
Expected local budget: up to a few hours, increasing with p; each cell is
checkpointed. Independent HPC cells require one CPU,4GB; no submission here.
Do not edit any running script or loaded helper. The runner and this design
must be committed, with a fresh passing full suite, before any model fit.

New scorer tests require actual risk direction and equality to the existing
partial-likelihood calculation, constant-risk C=NA and null logPL, and rejection
of nonfinite predictions with the intended error. RED9/12 then GREEN12/12.
The author-line cleanup and HTML/PDF regeneration are presentation changes,
not a second scientific change. Scientific defaults remain unchanged.

Pre-fit reference verification: archived D4 predates the sign correction, as
already documented in the2026-10-02 decision. Preserve its frozen risk=-eta
orientation, yielding training C0.6563073062 (raw stored eta C0.3436926938).
This is an existing convention, not orientation selected on held-out outcomes.
New fits use the current training-only orientation. No archived fit is edited.
Prepared-input checks pass, exact baseline/order identity and all24 design
cells verified; no model fit has run. Full suite **654 passed,0 failed**.
