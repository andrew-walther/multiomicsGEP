# Fixed-prior survival supervision findings

**Complete.** The registered diagnostic and independent reporting review are finished.

Date: 2026-11-14. Author: Andrew Walther. Registration: S-17.

## Result

The registered Step 1 grid is complete: 24 cases and 105 retained fitted
passes across real data and the two seed-1 synthetic diagnostics.

No tested Normal prior scale passes the registered real-data convergence screen.
All four unit-weight Normal pilots converge stably, but every scale has a failed
higher-target fit. Normal priors with σ_b = 0.5 and 1 support stable target-0.05
fits within the coefficient limit. Neither supports the required targets 0.25 and
0.5. The Laplace comparison also has no passing scale.

The registration requires a scale to pass every target before training CV or the
ten-seed simulation and null screen. Steps 2 and 3 are therefore not run. No new
arm is selected. The prespecified full program-structure analysis and predictive evaluation
under these priors are not run. The working model is unchanged. External scoring is not run.

[The registration](../plans/2026-11-14-bounded-risk-supervision-plan.md) and
[execution prompt](../prompts/2026-11-14-codex-bounded-risk-prompt.md) fix the
families, scales, targets and stopping rule. The reviewing notes at the end of the
[S-16 handoff](../reviews/2026-11-13-share-targeted-weight-handoff.md) replace
S-16's risk-variance null gate. S-17 does not use that historical gate.

## Step 1: real-data convergence

Training uses 273 TCGA and CPTAC patients, 139 deaths and the prepared 2,218-gene
D4 matrix. The model has K = 5, point-Laplace L and F, pooled survival without
strata and the A1 survival block weight. Every same-prior pilot uses weight one.
The target loading vector is selected at that pilot and tracked through at most
three warm-started recalibrations. There is no ramp. Every fit has a maximum of
500 sweeps.

The survival predictor is η = (YF)β. The monitored criterion before
canonicalization is expression log likelihood plus w_S times the Cox term,
minus loading and score KL terms and w_S times the coefficient KL. Weighting
both coefficient terms leaves the conditional β update untempered at fixed F,
while increasing survival curvature in the F update. The target share is a
working quadratic precision ratio along the selected loading direction. It is
not a fraction of total Cox information or a measured change in orientation.

Normal coefficients use a fixed shared N(0, σ_b²) prior. Laplace coefficients use
fixed scale σ_b/√2, giving the same variance. Fixing these distributions imposes
shrinkage. They have unbounded support and are not hard bounds on coefficients.

The table gives final target passes. Target cells list **sweeps, stable
convergence, maximum absolute coefficient mean**. All pilots are stable, and their
cells give sweeps and maximum absolute coefficient mean. A scale passes only
when its pilot and every final target pass are stable and their maximum
coefficient means are at most 4σ_b.

| Prior | σ_b | 4σ_b limit | Pilot | Target 0.05 | Target 0.25 | Target 0.5 | Scale gate |
|---|---|---|---|---|---|---|---|---|
| Normal | 0.1 | 0.4 | 143, 0.219 | 500, no, 0.610 | 500, no, 0.556 | 500, no, 60.358 | fail |
| Normal | 0.25 | 1 | 66, 0.350 | 500, no, 1.556 | 500, no, 1.598 | 500, no, 24.977 | fail |
| Normal | 0.5 | 2 | 67, 0.383 | 64, yes, 0.421 | 500, no, 4.133 | 500, no, 8.844 | fail |
| Normal | 1 | 4 | 67, 0.392 | 64, yes, 0.434 | 500, no, 8.818 | 500, no, 16.320 | fail |
| Laplace | 0.1 | 0.4 | 145, 0.281 | 500, no, 21.319 | 500, no, 11.043 | 268, yes, 17.349 | fail |
| Laplace | 0.25 | 1 | 98, 0.358 | 500, no, 15.194 | 500, no, 14.939 | 500, no, 23.525 | fail |
| Laplace | 0.5 | 2 | 66, 0.373 | 500, no, 14.969 | 500, no, 18.341 | 500, no, 17.523 | fail |
| Laplace | 1 | 4 | 66, 0.384 | 65, yes, 0.438 | 500, no, 19.598 | 500, no, 20.141 | fail |

Of the 24 final target states, four are stable and three also meet the
coefficient limit. All eight pilots meet both conditions. No complete scale
passes all targets.

The Laplace σ_b = 0.1, target-0.5 arm uses two calibration fits. Both converge
stably, at 330 and 268 sweeps, but the final coefficient maximum is 17.349 against
a limit of 0.4. Target 0.05 also uses two fits for Normal σ_b = 0.25 and Laplace σ_b =
0.25 and 0.5. Those initial passes converge, but recalibration fails. The 32 final
real arm records therefore represent 36 actual fits. No real pilot triggers the
no-signal guard. The final Laplace σ_b = 0.1, target-0.5 pass reaches the 100,000
weight cap. Its achieved tracked share is 0.444248, with loading correlation
0.029019 to the pilot, so the tracked program identity is weak. No real fit emits
a warning or task error.

For Normal σ_b = 0.1, target 0.05 and target 0.25 meet the risk-stability
thresholds at the sweep cap, but fail the current ELBO criterion. Risk stability
alone is not counted as convergence. At target 0.5 the same prior ends with
maximum absolute β of 60.358 and maximum absolute η of 19.606. In its final sweep,
the program with that coefficient has β = 36.046 before canonicalization and a
score SD of 1.674, giving β = 60.358 after scaling. This is a finite-run diagnostic,
not proof that the objective is unbounded.

[Detailed fit records](../../results/benchmark_sim/outputs/bounded_risk_2026-11-14/screen_fits.csv)
retain the native ELBO flag separately from stability, all calibration passes,
weights, warnings and elapsed seconds. The per-program share table reports every
program, including inactive programs. A fixed continuous prior does not switch
off survival curvature through an estimated zero prior variance. Its activity
counts therefore have a different interpretation from the working model.

Three of the eight final target-0.05 arms meet both the stability and coefficient
conditions. They give a limited positive result. Fixed priors permit
moderate directional shares that failed under the S-16 estimated prior. Their
tracked loading remains close to the same-prior unit-weight pilot:

| Prior | σ_b | Final weight | Achieved directional share | Absolute loading correlation with pilot |
|---|---|---|---|---|
| Normal | 0.5 | 392.274 | 0.042223 | 0.997342 |
| Normal | 1 | 372.044 | 0.062140 | 0.997177 |
| Laplace | 1 | 391.387 | 0.044017 | 0.996892 |

These correlations describe the single direction already tracked for
calibration. They do not replace the unrun all-program comparison, annotation
or CV analysis in Step 2. Target attainment is approximate: the inherited rule
accepts an achieved directional odds ratio within a factor of two of the target
odds, or stops after three converged passes. Its tolerance is not an additional
convergence gate.

## Step 1: synthetic convergence diagnostics

The strong-signal S-9 diagnostic uses p = 2,000, 1% molecular variance,
seed 1 and K = 7. Every family and prior scale is fitted at fixed weights
40, 100 and 1,000. All 24 fits reach 500 sweeps without stable convergence.
Three Laplace fits meet the native ELBO criterion at the final sweep but fail
the risk-stability window. No fit records a task error or warning.

The table summarizes the three fixed weights for each prior. Coefficient ranges
are final maxima across those three fits. The risk column is the largest final
maximum absolute η. The planted standardized coefficient is about 6.7, far
from every tested prior scale. These are convergence diagnostics under a
misspecified prior, not recovery or predictive comparisons.

| Prior | σ_b | Native ELBO criterion met | Range of final maximum absolute β | Largest final maximum absolute η |
|---|---|---|---|---|
| Normal | 0.1 | 0/3 | 1.144 to 112.202 | 109.651 |
| Normal | 0.25 | 0/3 | 32.872 to 35.245 | 104.142 |
| Normal | 0.5 | 0/3 | 33.862 to 45.657 | 126.401 |
| Normal | 1 | 0/3 | 31.695 to 40.777 | 107.320 |
| Laplace | 0.1 | 0/3 | 31.332 to 39.201 | 106.395 |
| Laplace | 0.25 | 2/3 | 31.445 to 33.148 | 100.605 |
| Laplace | 0.5 | 1/3 | 31.364 to 33.409 | 93.038 |
| Laplace | 1 | 0/3 | 30.831 to 31.620 | 99.688 |

The weak-cell seed-1 diagnostic uses the registered S-13 generator and K = 7.
Its eight pilots all converge stably. Four of the 24 final target states are
stable, all at nominal target 0.05. Every final target-0.25 and target-0.5 fit
reaches 500 sweeps. The table gives final-pass sweeps and stable convergence.
No guard, cap, warning or task error is recorded in this cell.

| Prior | σ_b | Retained fits | Pilot sweeps | Target 0.05: sweeps, stable | Target 0.25 | Target 0.5 |
|---|---|---|---|---|---|---|
| Normal | 0.1 | 9 | 97 | 107, yes | 500, no | 500, no |
| Normal | 0.25 | 5 | 105 | 500, no | 500, no | 500, no |
| Normal | 0.5 | 4 | 105 | 500, no | 500, no | 500, no |
| Normal | 1 | 6 | 105 | 53, yes | 500, no | 500, no |
| Laplace | 0.1 | 4 | 196 | 500, no | 500, no | 500, no |
| Laplace | 0.25 | 5 | 107 | 500, no | 500, no | 500, no |
| Laplace | 0.5 | 6 | 106 | 208, yes | 500, no | 500, no |
| Laplace | 1 | 6 | 105 | 53, yes | 500, no | 500, no |

The four stable target-0.05 arms all use three calibration passes. Stability
does not establish that the target share was attained:

| Prior | σ_b | Achieved directional share | Loading correlation with pilot | Final maximum absolute β |
|---|---|---|---|---|
| Normal | 0.1 | 0.006221 | 0.202400 | 0.687 |
| Normal | 1 | 0.051130 | 0.999821 | 0.660 |
| Laplace | 0.5 | 0.218469 | 0.999934 | 4.412 |
| Laplace | 1 | 0.051186 | 0.999817 | 0.658 |

Normal σ_b = 0.1 and Laplace σ_b = 0.5 remain outside the factor-of-two
odds tolerance after three converged passes. The former also has weak loading
identity relative to its original pilot. The Normal and Laplace σ_b = 1 arms
are within tolerance. The four-SD acceptance screen selects priors from the real
data only. It is not applied as a new gate to these synthetic diagnostics.
No recovery, validation concordance or null-safety conclusion is drawn here.



## Steps 2 and 3: not reached

No Normal σ_b passes Step 1, so training CV, program-structure analysis, the
weak/null ten-seed screen, paired-weight nulls and the S-14 recovery comparison
are not run. Their criteria are not failed results. They are untested. No
adoption candidate is established, including for the individually stable
Normal target-0.05 arms.

## Implementation and verification

The opt-in fitter trace records each YF score's raw SD, applied canonical scale,
and coefficient means before and after scaling. For a live program,
β_after = c_k β_before. A collapsed score uses scale one, so raw SD and applied
scale are recorded separately. Training matrices and posterior moments enter
this block. Its output is one aggregate row per program per sweep in
`history$s16`. It adds no diagnostic field when tracing is off.

The fitter's `stable_convergence_window` defaults to zero. A positive value
requires tracing and ELBO stopping. S-17 sets it to 50. The stopping block first
checks the current native ELBO criterion, then all 50 relative L2 changes in η
and the endpoint change in maximum absolute η. The window's initial maximum,
floored at machine epsilon, is the latter denominator. Iteration continues until
both criteria hold or sweep 500. `history$native_elbo_converged` records the
current native criterion separately. Default fitting and update formulas do not
change.

The shared `s16_fit()` helper accepts `prior_beta` through `fit_control`, with
Normal unchanged as its default. Each same-prior pilot and warm start receives
the fixed shared prior. `s16_calibrate(..., retain_pass_fits = TRUE)` retains
every fitted pass and trace. Its default is false. The dataset and real-data
helpers also accept separate two-step reference controls. An estimated-prior
working-model reference still requires its own unit-weight fit. A fixed-prior
pilot is not that reference.

`s17_controls(family, sigma)` returns the named opt-in fitter arguments.
`ashr::normalmix(1, 0, sigma)` creates the Normal prior, while
`ebnm::laplacemix(1, 0, sigma/sqrt(2))` creates the equal-variance Laplace prior.
These existing packages supply prior objects that the fitter accepts. No new
dependency is introduced. `yaml` reads project ranks and tolerance.
`s17_screen(family, stage)` prepares the registered training data and runs the
family's serial case grid. The separate real and simulation stages enforce the
registered ordering and permit two independent single-threaded processes.

`s17_case()` takes a dataset, family, scale, configuration and source signature.
It returns and atomically saves a private case object keyed by dataset, family
and scale. The object contains fits, calibration records, warnings, elapsed
seconds, second moments and final Cox weights. Complete cases and completed
arms are reused only with a matching signature. A nonconverged pilot remains
reported, with targets explicitly unavailable. Unexpected errors are saved and
rethrown. A no-signal guard reuses the pilot at weight one and does not establish
the requested stronger share.

`s17_stability(fit, sigma)` independently evaluates the saved trace, returning
a one-row data frame of window, convergence and coefficient-limit results.
Repeated sweep-level quantities are deduplicated across programs. `s17_export()`
writes these fit rows, per-program shares and traces as aggregate CSVs. It also
reports final and peak coefficients and risk, plus the largest within-window
risk deviation. Only the endpoint deviation enters the registered 1% gate.
`s17_summarize()` requires all 24 complete, matching case checkpoints before
combining exports. Its real-data pair gate uses only the pilot and final pass
of all three targets, with the final four-SD coefficient limit. The synthetic
results do not select priors for later steps.

The separate `tests/audit_bounded_risk_screen.R` script verifies saved states
without fitting a model. `audit_inputs()` reconstructs the three training
inputs. `audit_equal()` checks shapes and finite numerical agreement to relative
tolerance 10⁻¹⁰. This audit tolerance does not alter an acceptance threshold.
`audit_track()` returns original-pilot matches and achieved shares for every
available fit, using absolute Pearson loading correlation from base R's
`stats::cor()`. A missing match remains undefined with a reason.
`audit_screen()` checks the complete grid, source/input hashes, package versions,
fixed priors, moments, traces and recomputed Cox weights. It writes aggregate
check and tracking tables and exits nonzero for an audit failure. A recorded
nonconverged fit is a scientific outcome and can pass this consistency audit.

Fixing a Normal or Laplace scale imposes shrinkage, not a hard coefficient bound.
Canonicalization preserves LF′ and η, but rescaling q(β) while holding its prior
fixed generally changes the coefficient KL. The native objective is calculated
before rescaling. These fits test the registered algorithm and do not establish
an exact CAVI derivation under a fixed prior expressed throughout in standardized
coordinates. The added trace makes this limitation observable. The four-SD rule
is an empirical acceptance threshold.

## Reproduction and evidence

Run from the repository root with the original prepared training input and
checkpoint directory available. The source/input signatures and package
versions must match. Each family process is serial. Run each pair serially or
concurrently, with at most two processes:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
Rscript results/benchmark_sim/run_bounded_risk_screen.R normal real
Rscript results/benchmark_sim/run_bounded_risk_screen.R laplace real
# Complete both real partitions before starting these:
Rscript results/benchmark_sim/run_bounded_risk_screen.R normal sim
Rscript results/benchmark_sim/run_bounded_risk_screen.R laplace sim
# After both simulation partitions finish:
Rscript results/benchmark_sim/run_bounded_risk_screen.R summarize
Rscript tests/audit_bounded_risk_screen.R
Rscript tests/run_tests.R
```

Complete cases and completed arms within partial cases are reused. A stale
signature or unexpected case error stops execution. No external input is scored.
Patient-level matrices, risk vectors and checkpoints remain under ignored
`data/bounded_risk_2026-11-14/`. Only aggregate outputs are committed.

The output directory is
`results/benchmark_sim/outputs/bounded_risk_2026-11-14/`:

| File or pattern | Contents |
|---|---|
| `screen_fits.csv` | Every calibration pass, final-state flag, native and stable convergence, coefficient limit, weight, cap, warning and elapsed seconds |
| `screen_real_pair_gate.csv` | The eight family-by-scale real-data gates |
| `screen_program_shares.csv` | Activity, standardized coefficients and exact directional and per-gene shares for every program and fitted pass |
| `screen_loading_tracking.csv` | Original pilot matching and achieved shares, including nonconverged states, with explicit undefined-match reasons |
| `screen_completeness_audit.csv` and `_summary.csv` | Saved-state checks and their counts |
| `*_trace.csv` | Per-sweep objective, risk, coefficient, loading norm and canonicalization diagnostics |
| `*_fits.csv` and `*_shares.csv` | Case-specific exports of the fit and program tables |
| `pause_inventory.csv` and `pause_saved_fits.csv` | Historical 60-fit inventory at the requested usage pause, not the final denominator |

The numerical implementation is `e0c5697b`, the requested pause is `1421716a`,
and completed Step 1 evidence is `19f3cf29` on `s17-bounded-risk-supervision`.
No merge is performed.

The completed grid has **24 cases, 105 retained fitted passes and 88 final
arm states**. The 17 extra fitted passes are recalibrations, four on real data
and 13 in the weak cell. All 16 unit-weight pilots are stable. Across all final
arm states, 24 are stable and 64 are not. The real-data gates alone determine
whether the study proceeds. All eight gates fail.

The complete saved-state audit passes **1,611 checks with zero failures**.
It verifies all case and arm records, unchanged source/input hashes and package
versions, priors at every sweep, trace dimensions, posterior moments and final
Cox weights. The program table has 663 rows across all retained passes. There
are no guards, unavailable targets, recorded fit warnings or fit errors. There
is one capped fit, the real Laplace σ_b = 0.1 target-0.5 final pass.

The final full suite passes **879 counted checks, zero failures**, plus the
new standalone assertions. Its one previously recorded bounded-β
profile-likelihood fallback warning is a test diagnostic, not a warning from a
study fit. The default S-16 fixture passes exactly. The first suite attempt
found a collapsing noise-only stability fixture. A planted, nonzero-β fixture
now checks both rejection of an unstable window and later stable stopping.
No update algebra changed, so the conditional algebra-verifier requirement was
not triggered.

The final suite log is `full-suite-resumed-final.log`. Summary and audit logs
are `summarize-final-retry1.log` and `audit-final.log`, under the private
`data/bounded_risk_2026-11-14/execution_notes/` directory.

A fresh independent final review passes the completed diagnostic and reports
with no blocking issue or factual correction. It independently reconstructed
all 105 trace windows and native ELBO flags, checked 242,047 exported trace rows
against checkpoints, and verified the fixed prior over 37,963 saved sweeps.
It directly recalculated all finite directional shares and all 81 pilot matches,
without using the audit's matching helper as proof. Of 663 program-share rows,
481 are finite and 182 are undefined because a loading direction or its loading
second-moment denominator is zero. Those rows remain explicitly undefined.
They are not dropped or recoded as zero survival share.

The review also verified all 60 pause records against the final retained records
and all 24 current source/input signatures. It inspected the saved 1,611/0 audit
and 879/0 suite logs, but did not rerun those complete commands or fit a model.
The private review is
`execution_notes/findings-independent-review-final.md` under the same ignored
study directory. Its verdict approves reporting of the registered diagnostic,
not a model, predictive claim, adoption or external scoring.

### Timing and the requested pause

Elapsed times below are sums of the retained fits' recorded wall-clock seconds.
Concurrent workers overlap. These sums are not elapsed study duration and do
not include discarded in-flight work at the requested pause.

| Data | Prior | Retained fits | Stable fitted passes | Final arm states | Stable final states | Retained-fit seconds |
|---|---|---|---|---|---|---|
| Real | Normal | 17 | 7 | 16 | 6 | 1009.438 |
| Real | Laplace | 19 | 9 | 16 | 6 | 1131.662 |
| Strong S-9 | Normal | 12 | 0 | 12 | 0 | 1313.321 |
| Strong S-9 | Laplace | 12 | 0 | 12 | 0 | 1309.244 |
| Weak S-13 | Normal | 24 | 14 | 16 | 6 | 1965.601 |
| Weak S-13 | Laplace | 21 | 11 | 16 | 6 | 1546.016 |

The retained-fit total is 8,275.282 seconds. The two simulation partition logs
also preserve command wall times, including export and discarded interrupted
work in the original attempts:

| Partition | Before requested pause, seconds | Resumed successful command, seconds |
|---|---|---|
| Normal simulation | 1,384.57 | 1,957.36 |
| Laplace simulation | 1,253.29 | 1,672.67 |

At the usage-limit pause, 60 retained fitted passes were saved in 15 complete
and two partial cases. The two workers were deliberately interrupted. On
Andrew's authorization, the 60 saved passes were reused and 45 further retained
passes completed. Interrupted unsaved work is excluded from the final fit
count, rather than counted as a scientific failure. The historical pause
inventory is preserved separately from the final denominator.

The first combined export failed with a filesystem I/O error while rewriting
one generated trace CSV. The disk reported about 100 GiB free. The complete
trace remained in its private checkpoint. One unchanged export retry succeeded,
and the repaired CSV matched that checkpoint at tolerance 10⁻¹². The subsequent
full audit passed. No model was rerun and no numerical source changed for this
recovery. The original failure log is retained as `summarize-final.log`.


## Decisions and limits

The study uses a new grid orchestrator around the existing `s16_fit()` and
`s16_calibrate()` helpers. The historical diagnostic script initializes its
fixed output folders when sourced and fits isolated weights without
recalibration. Extending that top-level workflow would mix the historical
S-16 namespace with the new case grid. The S-17 wrapper instead supplies
opt-in prior controls, tracing and retained calibration passes through the
shared fitting path. `s16_diagnostic_fit()` and its historical CLI remain
unchanged and are not used for S-17. This file-layout choice uses the prompt's
allowance for a new file where extending the old workflow is awkward.


The 50-sweep stability window must be observed before a fit can pass. An opt-in
stopping control continues after an early ELBO threshold until both the current
ELBO criterion and the registered window hold, or the fit reaches 500 sweeps.
This operationalization was independently reviewed and recorded before fitting.
It prevents an unobserved window from becoming either a pass or an automatic
failure. Default fitting behavior is unchanged.

The 1% maximum-risk check compares the beginning and end of the window. The
largest intermediate deviation is reported separately and does not tighten the
registered gate. The four-SD coefficient threshold applies to final means, while
peak coefficients are also reported. A no-signal guard is a weight-one fallback,
not an achieved stronger target.

The prepared D4 screen and cohort scaling are fixed historical preprocessing.
There is no nested preprocessing or new CV performance estimate here. The
500-sweep limit cannot distinguish every slow fit from a fit without a stable
fixed point. The synthetic strong-signal coefficient is about 6.7, outside all
registered fixed-prior scales. That cell assesses convergence, not recovery under
a correctly specified prior.

## Text for the shared status files

S-17 tested fixed shared Normal and equal-variance Laplace priors on survival
coefficients at four scales, with same-prior pilots and share-targeted weights.
No Normal scale passed the complete real-data stable-convergence and coefficient
screen. Moderate target-0.05 fits converged for Normal σ_b = 0.5 and 1, but the
required stronger targets failed. The registration therefore stopped before
training CV, program characterization and the ten-seed null/recovery study.
The unit-weight working model remains unchanged. No new external score was
computed and no fixed-prior arm is adopted.

## Draft paragraph for the progress book and Chapter 4

**Draft, for review.** We tested whether fixing the prior variance of the
standardized survival coefficients would permit stronger supervision of the
molecular loadings. Priors with standard deviations of 0.1, 0.25, 0.5 and 1 were
compared at unit survival weight and at weights calibrated to directional
curvature shares of 0.05, 0.25 and 0.5. The unit-weight Normal pilots converged
stably. A target share of 0.05 also converged under the two weaker Normal priors,
with standard deviations of 0.5 and 1. Neither prior supported stable fits at
both stronger targets, and no prior scale passed the complete registered
screen. We therefore did not evaluate these models by training cross-validation
or the registered null and recovery checks. The experiment does not establish
a new supervision weight or show that predictive or molecular-program structure
improved. Fixing a prior scale supplies shrinkage, while canonicalization still
rescales the coefficient posterior, so it does not itself impose a hard bound
on the standardized coefficients.
