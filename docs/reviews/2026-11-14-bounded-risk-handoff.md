# Review handoff: fixed-prior survival supervision

**Complete.** The registered diagnostic and independent reporting review are finished.

Date: 2026-11-14. Author: Andrew Walther.

Addressed to the reviewing Claude Code session.

## What to check first

Read the [S-17 registration](../plans/2026-11-14-bounded-risk-supervision-plan.md)
and [execution prompt](../prompts/2026-11-14-codex-bounded-risk-prompt.md), then
recompute the real-data gate from the saved final passes. Every Normal scale
has a failed required target. Steps 2 and 3 are therefore not reached. There
is no new CV, structure, null, recovery or external-scoring result to adopt.

The [findings](../findings/2026-11-14-bounded-risk-findings.md) distinguish the
limited success at target 0.05 from the complete scale gate. Normal σ_b = 0.5
and 1 achieve stable shares of 0.042223 and 0.062140 at weights 392.274 and
372.044, with tracked loading correlations of 0.997342 and 0.997177 against
their own fixed-prior pilots. These individual arms cannot proceed under the
registered all-target gate.

## Commits and execution scope

- `e0c5697b`: opt-in canonicalization trace, 50-sweep stable stopping, coefficient
  prior override, retained recalibration passes, separate reference controls,
  resumable screen and regression tests. Independently reviewed before fitting.
- `1421716a`: requested usage pause, 60 saved fits, partial-case inventory and
  resume instructions.
- `19f3cf29`: complete Step 1 evidence, 105 retained fits, final exports and
  1,611 successful saved-state checks.
- Final documentation commit: the branch commit containing this completed
  handoff, available in the branch history.

The branch is `s17-bounded-risk-supervision`, created from main at `6d12a0dd`.
The implementation and completed Step 1 evidence are pushed. Every new scientific fit uses its
frozen numerical source. Later changes are tests and saved-state audits, outside
those numerical source hashes. No merge or external scoring is authorized here.
No protected status, roadmap, decision, progress-book or prelim source is edited.

## Gate verdicts

All eight real unit-weight pilots converge stably and meet their four-SD
coefficient limits. No family-by-scale pair passes all targets.

| Family | σ_b | Coefficient limit | Final target 0.05 | Final target 0.25 | Final target 0.5 | Pair |
|---|---|---|---|---|---|---|
| Normal | 0.1 | 0.4 | not stable, β max 0.610 | not stable, 0.556 | not stable, 60.358 | fail |
| Normal | 0.25 | 1 | not stable, 1.556 | not stable, 1.598 | not stable, 24.977 | fail |
| Normal | 0.5 | 2 | stable, 0.421 | not stable, 4.133 | not stable, 8.844 | fail |
| Normal | 1 | 4 | stable, 0.434 | not stable, 8.818 | not stable, 16.320 | fail |
| Laplace | 0.1 | 0.4 | not stable, 21.319 | not stable, 11.043 | stable, 17.349 | fail |
| Laplace | 0.25 | 1 | not stable, 15.194 | not stable, 14.939 | not stable, 23.525 | fail |
| Laplace | 0.5 | 2 | not stable, 14.969 | not stable, 18.341 | not stable, 17.523 | fail |
| Laplace | 1 | 4 | stable, 0.438 | not stable, 19.598 | not stable, 20.141 | fail |

There are 36 actual real fits and 32 final arm records. Four targets require
two calibration passes. The final Laplace σ_b = 0.1, target-0.5 fit reaches the
100,000 weight cap. It is stably converged but fails its coefficient limit.
Its achieved tracked share is 0.444248, but loading correlation to the original
pilot is only 0.029019. No real guard, warning or task error is present.

Only Step 1 is implemented for S-17. The conditional stop leaves S-17 CV
selection and the revised null-gate orchestration unimplemented. Historical
S-16 selection, variance-null gates and unfit capped aliases must not be used
as S-17 substitutes.

Step 2 eligibility and selection are untested, not failed. Step 3's revised
false-active coefficient null check and both recovery comparisons are also
untested. The S-16 risk-variance null rule is not used. Step 4 is not run.

## Synthetic convergence diagnostics

All 24 S-9 strong-signal fits are complete at the registered weights 40, 100
and 1,000. Every fit reaches 500 sweeps without stable convergence. Three
Laplace fits meet the final native ELBO criterion but fail the risk window.
No task error or warning is recorded in these fits. The planted standardized
coefficient is about 6.7, far from every prior scale, so these are convergence
diagnostics under a misspecified prior rather than recovery evidence.

The weak-cell seed-1 grid is complete: 45 retained fits, 32 final arm states.
All eight pilots are stable. Four of 24 final targets are stable, all at nominal
share 0.05. These are Normal σ_b = 0.1 and 1, and Laplace σ_b = 0.5 and 1.
Normal 0.1 achieves share 0.006221, while Laplace 0.5 achieves 0.218469 after
three converged passes. Both miss the factor-of-two odds tolerance. The scale-1
arms achieve about 0.0511. Every final weak target 0.25 and 0.5 reaches 500
sweeps. The synthetic results do not enter the real-data pair gate.

## Stable convergence and statistical interpretation

The new grid driver delegates to `s16_fit()` and `s16_calibrate()` instead of
altering the historical global-folder diagnostic CLI. `s16_diagnostic_fit()`
remains unchanged and unused by this study. This implementation choice avoids
mixing the historical namespace with the recalibrated eight-prior case grid.
Check it against the prompt's allowance for a new file when extending the old
workflow is awkward. No analytical choice depends on this file layout.


A positive `stable_convergence_window` is opt-in and equals 50 here. It keeps
iterating until the current native ELBO tolerance and the registered risk window
both pass, or the fit reaches 500 sweeps. The native ELBO flag is retained
separately. This avoids accepting or rejecting a short native run solely because
its stability window was not observed. A fresh independent review endorsed this
operationalization before the first fit.

Each of the last 50 L2 relative risk-vector changes must be below 0.001. The
1% maximum-risk condition compares the first sweep's starting max absolute η
with the final sweep's ending value. Intermediate maximum deviation is exported
but does not enter this gate. The four-SD coefficient check uses final means,
while peak coefficients are also exported. These choices were fixed before
outcomes, without changing the grid or thresholds.

The Normal and equal-variance Laplace priors are fixed shrinkage distributions,
not hard coefficient bounds. Canonicalization rescales q(β), but the fixed g_β
is unchanged. It preserves risk and reconstruction, while generally changing
the coefficient KL. The native objective is calculated before that rescaling.
The historical invariance comment must not be treated as a proof for these
fixed-prior arms. The results test the registered algorithm, not an exact
standardized-coordinate CAVI derivation or a general impossibility of bounded
survival supervision.

Calibration tracks the pilot loading by maximum absolute Pearson correlation,
including inactive programs. The initial no-signal guard alone can return weight
one. Guarded targets reuse that pilot and do not establish that the requested
share was attained. A nonconverged pilot is preserved, with dependent target
fits explicitly unavailable. Unexpected errors are saved and stop the case.

## Provenance and verification

The implementation and extended reference-routing tests pass the full suite:
879 counted checks, zero failures, plus the new standalone assertions. The
existing bounded-coefficient optimizer profile-likelihood fallback warning is
reported in the test log. It is not a scientific fit warning.

The first suite run found a degenerate noise-only fixture in the new window
test. Its coefficients collapsed to zero. A planted-signal fixture now verifies
nonzero coefficients, rejection of an unstable 50-sweep window, and stopping at
a later stable window. The historical S-16 default-fit fixture still passes.
No update algebra changed, so the conditional derivation-verifier requirement
was not triggered.

A fresh Sol code review cleared Step 1 before fitting. It checked the fixed
prior families and warm starts, exact moments and Cox weights, original pilot
tracking, stable-window implementation, local mock isolation and unchanged
baseline behavior. Its private report is
`data/bounded_risk_2026-11-14/execution_notes/code-review.md`.

The source snapshot is archived privately. Every case saves the same numerical
source/input hashes and package versions. The prepared training matrix and
patient-level fit states remain in ignored `data/`. At most two serial family
processes run at once, with one BLAS thread each. No live numerical source is
edited while a partition is running.

The complete grid has 24 complete cases, 105 retained fitted passes, 88 final
arm states and 663 program-share rows. All 1,611 saved-state audit checks pass.
There are no guards, unavailable targets, recorded fit warnings or fit errors.
The sole capped fit is described above. Retained-fit wall seconds sum to
8,275.282 across concurrent workers. The findings give per-cell totals and
partition timings, with interrupted work treated separately.

At Andrew's usage pause, 60 fitted passes were preserved. They were reused on
resume, and 45 further retained passes completed. An export I/O error later
partially rewrote one CSV, while its full checkpoint trace remained intact.
One unchanged retry restored it and succeeded. The repaired trace matched its
checkpoint at tolerance 10⁻¹², then the complete audit passed. No model was
rerun for this repair. Both original and retry logs are preserved privately.

The fresh final suite also passes 879/0. Its log is
`execution_notes/full-suite-resumed-final.log`.

## Final independent review

A fresh independent final review passes with no blocking issue or factual
correction. Before reading the reports, it reconstructed the gates and
all 105 trace windows. It then checked all 242,047 trace rows against the
checkpoints, all 37,963 sweep-prior records, and all 81 pilot-loading matches.
Direct checkpoint calculations confirmed the reported moderate real shares,
the capped Laplace match, and the weak target misses.

All 481 finite program shares agree with its direct calculations. The remaining
182 rows have zero loading directions or loading second-moment denominators.
They are explicitly undefined and retained. All 60 pause fit records match the
final records, and all 24 case source/input signatures agree with current files.

The reviewer inspected the saved audit and suite logs, rather than rerunning
those complete commands, and fit no model. Its private report is
`data/bounded_risk_2026-11-14/execution_notes/findings-independent-review-final.md`.
The verdict clears reporting of the registered diagnostic. Scientific adoption,
any new registration and any later external scoring remain separate decisions.

## Open questions and results ready for Chapter 4

The real convergence table and moderate-share loading-tracking result can be
reported as diagnostics with their scope stated. The synthetic strong-signal
screen is deliberately prior-misspecified and cannot validate recovery. No
predictive improvement, full program comparison, null safety or adoption claim
is ready for Chapter 4.

Questions for Andrew after review are whether a separately registered study
should evaluate the individually stable moderate-share arms, and whether a
standardized-coordinate constraint or derivation should precede another
strong-supervision grid. Neither question authorizes further work on this branch.

The finite 500-sweep traces do not establish unboundedness. Predictive and
program-structure quantities for S-17 remain unknown because their registered
entry gate fails. The fixed-D4 preprocessing is historical and no new fully
nested pipeline estimate is produced.

The findings include proposed shared-status text and a labelled draft paragraph
for the progress book and Chapter 4. Andrew's scientific adoption decision and
any later external scoring remain separate from this diagnostic report.
