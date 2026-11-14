# S-16 continuation checkpoint, paused at Andrew's request

**Historical pause record:** execution subsequently resumed and completed Steps
2 to 5 under the registered fallback. Use the [final handoff](2026-11-13-share-targeted-weight-handoff.md)
and [findings](../findings/2026-11-13-share-targeted-weight-findings.md) for the
completed evaluation. The pause and resume instructions below remain provenance.

Date: 2026-11-13. Author: Andrew Walther.

## Resume here

Andrew requested a pause because the rolling usage limit was about to reset
in approximately two hours. Do not resume until he asks. This pause is not
a statistical stop condition. Work remains in Steps 3 to 5.

Branch: `s16-share-targeted-weight`. Latest implementation commit:
`c294866b`. Amendment 1 is `5e51438b`, and the amended helper is
`564dfe66`. These commits were pushed to the authorized origin,
`https://github.com/andrew-walther/multiomicsGEP.git`.
Do not merge into main or run Step 6. The original prompt and amended plan
remain controlling. Shared status files, the progress book and prelim chapter
must remain untouched. Final proposed text belongs in the findings note.

## Completed and verified

- Amended exact curvature uses Σ_i E[l_ik²] Σ_j τ_j f̂_jk².
  The mean-τ expression remains only in the approximation.
- The pilot target is fixed and tracked by absolute loading correlation
  through recalibration, regardless of later activity.
- Historical checks use final states, without comparing the inside-sweep trace.
- Implemented opt-in final-state tracing and a 50-sweep geometric weight ramp.
  Defaults are unchanged. Exact pre-change unit-weight comparison passes.
- All 16 Step 2 diagnostic fits completed and their checkpoints were inspected.
  No ramp passes the full convergence gate. The ramp is not adopted.
- Only target 0.01 is carried into Step 3 under the registered raw-fit fallback.
  Its full-training diagnostic weight is 82.75367, convergence occurs in
  14 sweeps, and the tracked achieved share is 0.00873483 with loading
  correlation 0.999432 to the pilot. Targets 0.05, 0.1, 0.25 and 0.5
  reach the 500-sweep cap. See committed `diagnostic_summary.csv` and
  `diagnostics_tracked_programs.csv` for every result.
- Full suite: 879 counted checks pass, zero fail. New standalone assertions
  also run but are outside that count. One existing bounded-β optimizer
  fallback warning was surfaced. Derivation verification: 109/109 pass.
- Independent Sol code review cleared the raw 0.01 workflow without the ramp,
  including simulation, real CV and program-structure integration. Final
  scientific findings review is still required after results are assembled.
- Off-grid smoke seed 99 passed: pilot, two-step and target converged.
  Target 0.01 achieved 0.009967125 in one pass.

The historical amended exact audit has directional shares 0.683, 0.505 and
0.201 in the three checked S-9 cells. The largest directional closed-form
ratio error is 2.213-fold, so do not repeat the amendment's preliminary
mean-precision numbers as revised exact values. Report this discrepancy
without changing the registered threshold.

## Screen state at pause

The two-worker screen initially failed to create its localhost PSOCK socket
in the sandbox. Running the same frozen command with elevated local socket
access succeeded. Andrew's later pause request interrupted the parent and
both workers. Their exit was confirmed. No screen process should remain.

One completed atomic task checkpoint is preserved:
`data/share_targeted_ws_2026-11-13/screen/weak_s01.rds`.
It has schema 1, `complete = TRUE`, all `signal`, `null` and `paired` results,
no task error, elapsed time 27.067 seconds, and unchanged source fingerprints.
The interrupted in-progress tasks have no completed checkpoint. Do not
delete the completed task or count interrupted work as scientific failure.
No full-screen verdict has been calculated yet.

Resume the exact frozen command from the repository root, with permission
for local PSOCK sockets if the sandbox requires it:

```sh
Rscript data/share_targeted_ws_2026-11-13/screen_driver.R screen 2 \
  data/share_targeted_ws_2026-11-13/screen_settings.rds \
  > data/share_targeted_ws_2026-11-13/screen_stdout.log 2>&1
```

Settings contain `list(allowed_targets = 0.01, fit_control = list())`.
The cache will reuse the completed task. Do not edit source while running.
The frozen driver hash is
`c24c0b7b8e0348074dbc0d4bc5d6a8af275556d014cbc2cd9d17f4df6b98fa09`.
Worker source dependencies are hashed and must remain unchanged.
At most two CPU workers, one BLAS thread each.

Expected completion: 20 task checkpoints covering 30 primary datasets
(10 weak, 10 paired permuted-survival null, 10 low-variance), 10 secondary
paired-weight null arms, and both references per primary dataset.
For one target the aggregate table should contain 100 rows:
30 target rows, 10 secondary-null rows and 60 reference rows.
Check field completeness, all warning/error records and convergence counts.

## Remaining execution

1. Finish Step 3, review every denominator and compute the registered primary
   null gate. Commit and push the aggregate screen results. No target that
   fails this gate proceeds to real data.
2. Write a new ignored real settings RDS with only the passing targets and
   no ramp. `numeric(0)` targets runs the unit-weight and two-step references
   even if target 0.01 fails. Run:
   `Rscript results/benchmark_sim/run_share_targeted_weight.R real 2 <settings.rds>`.
   Same local-socket permission may be needed. Complete five folds and full
   fits, inspect eligibility and select by mean fold CV partial likelihood.
3. Run `structure 1 <settings.rds>` with that same real settings file.
   It produces full correlation matrices, nested-Cox Δℓ, semi-partial R²,
   the S-11 program rule and Classical, basal-like, iCAF and myCAF enrichment.
   Signature inputs already exist locally and the annotation check passed.
   Commit and push Steps 4 to 5. Step 7 remains optional after required work.
4. Rewrite the existing findings and review handoff, whose initial stop
   narrative is superseded by Amendment 1. Include all numeric gates,
   exact commands, judgments, limits, driver errors and draft shared text.
   Obtain a fresh independent scientific findings review before final commit.
   Push and stop. Do not run external scoring or merge.

All new fits save posterior moments, τ and final Cox weights under ignored
`data/`. Only aggregate CSVs are tracked. The real runner applies fold-specific
training centers and verifies any target's null gate against current source.

## Decisions and limitations to preserve in the final report

- At most three target-fit passes follow the pilot. Keep its loading vector
  fixed. Loss of activity is recorded without reapplying the pilot guard.
- Signal false activity follows S-9: exclude the matched program only when
  |r| ≥ 0.7. Under the null every active program is false.
- The primary false-active variance fraction is
  Var(sum of false-active risk components) / Var(total risk), averaged over
  replicates. Signed covariance allocation is reported separately. This
  definition was fixed before any null outputs. No false actives gives zero.
  An undefined fraction with false actives does not pass the gate.
- Null gate requires ten unique finite evaluated rows, mean false-active ≤ 1
  and mean variance fraction < 0.01. Convergence is reported separately.
- CV selection uses mean fold CVPL, not the historical summed quantity.
  Within 0.5 units, weaker supervision means smaller median fold weight,
  then smaller target. Complete finite CVPL and C are required in addition
  to the registered ≥4 converged folds and converged full fit. Disclose
  this complete-score requirement and tie definition.
- The ramp failed convergence. There is no successful divergent-case recovery
  fixture. Its schedule tests do not establish an effective numerical fix.
  The pre-change unit-weight fixture has zero final factors and is narrow.
- Two diagnostic drivers failed after all requested fits were saved because
  the driver was edited while R was reading it. Both failures are disclosed
  in the diagnostic report. Complete state fields and trace row counts were
  independently checked. Raw transient source was not snapshotted before
  editing. Later ramp/fallback provenance and frozen drivers were recorded.
- Multi-target cap aliases are explicit unfit provenance rows. They do not
  borrow another fit's final prediction and do not affect this single-target
  fallback.

## Local evidence and reports

Aggregate outputs are under
`results/benchmark_sim/outputs/share_targeted_ws_2026-11-13/`.
Private checkpoints, frozen drivers and logs are under
`data/share_targeted_ws_2026-11-13/`.
The full-suite log is `tests_final.log` there.
Diagnostic summaries rebuild with
`Rscript results/benchmark_sim/summarize_share_targeted_diagnostics.R`.

Preserved execution notes under that data directory's `execution_notes/`:
`s16-code-review.md`, `s16-diagnostic-report.md`,
`s16-real-runner-report.md`, `s16-structure-module-report.md`, and
`s16-execution-ledger.md`. These supply the detailed code review,
test evidence, driver-error record and interfaces for final documentation.
