# S-16: share-targeted survival weight, stopped at the checkpoint audit

Date: 2026-11-13. Author: Andrew Walther.

Registration: [share-targeted survival weight plan](../plans/2026-11-13-share-targeted-survival-weight-plan.md).
Execution instructions: [S-16 prompt](../prompts/2026-11-13-codex-share-targeted-weight-prompt.md).
Branch: `s16-share-targeted-weight`, from `main` at `7b5ea2c8`.

## Summary

Execution stopped during Step 1. The saved S-13 and S-15 states do not contain
the noise precisions needed for the registered exact precision-share checks.
The plan permits replacing missing E[L²] with E[L]² as a labelled bound. It
does not specify a replacement for missing τ. No new study fit was run, no survival
weight was selected, and no fitter update was changed. Steps 2 to 5 remain
unassessed. External scoring was not run.

The implemented `share_targeted_ws()` calculates the registered directional
ratio and its inverse at a supplied converged training state. Its behavioral
tests use known numerical answers. The audit also calculates the available
S-9 shares and the permitted K = 5 bound. These are historical diagnostics,
not evidence that recalibration, large-weight convergence or the null gate
has passed.

## Step 1: saved-state availability

The inventory reads the local RDS files without fitting a model. Historical
checkpoints are left in their existing locations. New tracked outputs contain
aggregate quantities only.

| Source | Files | States available | Limitation |
|---|---:|---:|---|
| S-9 grid | 450 | 450 full fits | Final moments and precision trace available, but trace inputs were not saved at the time of each update |
| S-13 | 840 | 839 reduced states | No τ, E[L²] or full trace. Saved `L_norm` is Σ E[L]² |
| S-15, historical rank 7 | 216 | 36 reduced states | No τ or E[L²]. The 180 inner fits have no saved state |
| S-15, working rank 5 | 216 | 36 reduced states | Same omissions as rank 7 |
| K = 5 S-11 | 144 | 144 reduced states, 24 with E[L] | τ available. No E[L²] or precision trace. Full-data fits permit the registered E[L]² bound |

The missing S-13 state belongs to the previously reported fitting error, not
a failure in this audit. The source save calls independently confirm the
field omissions in `run_yfb_supervision_strength.R`,
`run_yfb_supervised_real_refit.R` and `run_yfb_supervised_real.R`.
The search found no linked full-state copies of the S-13 or S-15 fits.
It did not establish that an unlinked copy cannot exist elsewhere.

The fitted τ cannot be reconstructed exactly from these reduced states.
Its update uses both loading second moments, which are also absent in
S-13 and S-15. Using τ = 1, fitting τ again, or substituting a different
study's state would introduce an unregistered approximation or a new fit.
Those substitutions were not made.

### Precision-share calculation and interpretation

For each program the helper evaluates the Cox working weights at the
supplied training predictor η = Y E[F] E[β]. It computes

\[
r_k = \frac{w_S E[\beta_k^2]\sum_i W_{ii}(Y\hat f_k)_i^2}
{\bar\tau\sum_i E[l_{ik}^2]},\qquad s_k=\frac{r_k}{1+r_k}.
\]

The denominator uses mean τ, as registered. With heterogeneous gene
precisions, actual expression curvature along the unit loading direction
instead uses Σ_j τ_j f̂_jk². Mean τ is a proxy for that weighted average,
so the registered directional quantity is not an exact Hessian share in
general. The code retains the registered definition. This distinction
requires review before interpreting the target as an actual curvature share.
The per-gene calculation uses the individual τ_j. The closed forms use standardized coefficient b_k,
program PVE and mean working weight, and never set the returned weight.
Replacing E[L²] with E[L]² reduces the denominator. The resulting ratio
and share are upper bounds conditional on the saved remaining quantities.
The inverse weight is consequently a lower bound, rather than an exact
calibration of the target.

The historical `shareA_median` trace uses W from the start of a sweep and
moments from each program's sequential update. Later program updates,
the τ update and canonicalization precede the saved final state. A
final-state recomputation therefore measures a different state. A discrepancy
cannot by itself establish that the precision formula is wrong, and an
exact reproduction would require unsaved update-state quantities.

### Available numerical checks

The diagnostic runner checks seed 1 at p = 2,000 and unit weight in the
three S-9 variance cells. It also checks the three active programs in
the K = 5 S-11 full-data checkpoint. These four historical fits converged.
They are a subset audit, not the required comparison across every cell.

| S-9 planted variance | Matched program | Registered r_k | Registered s_k | Per-gene median share | Recorded trace share | Difference |
|---|---:|---:|---:|---:|---:|---:|
| 0.25% | 3 | 2.47117 | 0.711913 | 0.21332190 | 0.21330614 | +1.58 × 10⁻⁵ |
| 1% | 1 | 1.47598 | 0.596120 | 0.05018175 | 0.05018299 | −1.24 × 10⁻⁶ |
| 5% | 1 | 0.62983 | 0.386439 | 0.00869508 | 0.00869559 | −5.10 × 10⁻⁷ |

The per-gene shares match the trace when rounded to three decimal places
in these examples. The registration does not specify rounding precision.
Four-decimal equality fails in the 0.25% example, so no stronger equality
claim is made. Maximum symmetric fold error for the closed-form ratios
is 1.135 for the directional quantity and 1.649 for the per-gene quantity.
Both are within factor 2 in these three examples.

| K = 5 S-11 program | Registered directional share, upper bound | Per-gene median share, upper bound | Closed directional ratio fold error | Closed per-gene ratio fold error |
|---|---:|---:|---:|---:|
| 2 | 0.000147213 | 7.027 × 10⁻⁷ | 1.048 | 1.225 |
| 4 | 0.000095243 | 1.164 × 10⁻⁶ | 1.041 | 1.289 |
| 5 | 0.000170209 | 2.478 × 10⁻⁶ | 1.039 | 1.218 |

These bounds are conditional on the registered mean-τ denominator. They
cannot replace exact S-15 values or verify S-15 trace agreement. No
numerical shares are imputed for S-13 or S-15.

Aggregate sources under
`results/benchmark_sim/outputs/share_targeted_ws_2026-11-13/` are
`checkpoint_inventory.csv` and `saved_share_diagnostics.csv`.
`audit_run.csv` records 6.571 elapsed seconds and zero new fits for the
recorded audit run. A rebuild updates elapsed time.
Rebuild from the repository root without running a fit:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript results/benchmark_sim/audit_share_targeted_ws.R
```

The runner loads the existing generator, reconstructs the three synthetic
training inputs, reads the saved fits and prepared real training matrix,
then writes only per-program diagnostics and inventory counts. It does not
load external outcomes. Correlations with empty S-9 loading columns produce
three zero-SD warnings. They remain visible. Nonfinite candidate correlations
are assigned zero for matching, following the existing scorer.

`audit_one()` takes a converged historical fit and selected program indices.
It joins the helper's diagnostics to available trace shares and computes
symmetric approximation errors as max(exact/approximate, approximate/exact).
Its output has one aggregate row per selected program. The runner's
`inventory_checkpoints()` reads each directory's RDS files and counts field
availability. Missing directories fail loudly. It returns one count row per
source without extracting patient-level values.

### Weight function: purpose, logic and use

`code/share_targeted_ws.R` takes a converged canonical fit, its numeric
training matrix with n patients and p genes, survival times, event indicators,
the target share and the current positive weight. The posterior matrices
are n × K for L and p × K for F, with K coefficient moments and p noise
precisions. No held-out outcome enters the calculation.

The function validates the inputs, calculates the per-program diagnostics,
and chooses the survival-active program with the largest standardized
absolute coefficient. Activity follows the existing 1.96-SD rule with a
positive posterior SD and a live loading column. It then returns
`w_current * (target_share / (1 - target_share)) / r_k`, capped at 100,000.
With no active program it returns 1 and flags the guard. The return value
is a list containing the weight, chosen program, per-program data frame,
ratio and guard, cap and missing-second-moment flags.

The calculation holds the posterior state fixed. It does not fit a model
or guarantee that a subsequent fit will achieve the target. No recalibration
runner was implemented before the stop. The helper applies its no-signal
guard whenever called. A future runner must preserve the distinction between
the registered unit-weight pilot guard and loss of activity during later
passes. This later-pass behavior needs an explicit decision before execution.

The implementation uses base R matrix operations and the existing fitter's
`calc_cox_taylor_yf()` for Cox weights. It introduces no package dependency.
Tests isolate the weight algebra with known working weights and exercise
selection, inversion, cap, guard and factor rescaling. The audit uses the
actual Cox function on historical states.

Run the behavioral checks and full suite from the repository root:

```sh
Rscript tests/test_share_targeted_ws.R
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript tests/run_tests.R
```

Final verification passed: 879 counted checks and zero failures, plus the
new `stopifnot()` behavioral assertions, which are not included in that
count. The focused script and actual-Cox integration check passed. The
suite emitted its existing bounded-β optimizer profile-likelihood fallback
warning. `git diff --check` also passed. No fitter source changed.

## Steps 2 to 5: not run

| Registered check | Verdict | Evidence |
|---|---|---|
| Step 1 trace matches to rounding | Not established | Update-state inputs are absent from final checkpoints |
| Step 1 closed forms within factor 2 in every required cell | Incomplete | Exact S-13/S-15 quantities unavailable |
| Step 1 recalibration, at most 3 passes | Not run | Zero new study fits |
| Step 2 convergence within 500 sweeps without excessive post-ramp decrease | Not assessed | Zero diagnostic or stabilization fits |
| Step 3 primary null mean false-active ≤ 1 and false-active risk variance < 1% | Not assessed | Zero S-16 screen fits |
| Step 4 at least 4/5 converged folds and converged full fit | Not assessed | Zero S-16 CV/full fits |
| Step 4 likelihood selection and weaker-supervision tie rule | Not applied | No eligible S-16 arms evaluated |
| Step 5 program structure and annotation | Not run | No eligible S-16 full fits |

The Step 2 fallback to convergent targets was not invoked. It applies after
diagnosing and testing a stabilization change, which did not occur here.
Step 7 was not started. Step 6 remains reserved for review and Andrew's
approval.

## Stop decision and limits

The execution prompt says to stop and write the handoff if a registered
requirement cannot be met. The repository instructions also require reporting
a blocker without silently substituting a workaround. Execution stopped
because the exact historical comparison cannot be completed with the saved
states under the no-new-fits requirement. This is a provenance limitation,
not a scientific rejection of share-targeted supervision.

The available calculations were completed to leave the review with concrete
evidence. Their scope is recorded alongside the outputs. No target, cap,
selection rule or pass threshold was changed. No shared status file,
progress-book chapter, week plan or prelim chapter was edited. The root
README documents only the new utility and how to run its audit.

## Independent review

A fresh Sol agent reviewed the code, tests, audit runner, findings and handoff
against the registration. It independently inspected representative actual
S-13, S-15 and K = 5 S-11 checkpoint fields and checked the aggregate tables.
It found no blocking defect in the registered weight inversion or the
stopped deliverable. Its approval was conditional on passing final tests
and completing commit bookkeeping. It confirmed the missing-state blocker
and identified the mean-τ curvature limitation documented above.
The final tests passed and the Step 1 commit is `6531b160`.

The reviewer ran the focused behavioral test but did not independently
rerun the full suite, every checkpoint inventory read or every numeric audit
calculation. The review supports preserving the partial work for further
review. It does not approve scientific adoption or completion of S-16.

## Branch delivery

Step 1 is committed locally as `6531b160`. Findings, handoff and their
indexes are preserved in a separate documentation commit on the same branch.
The requested push was rejected by automatic approval review. A retry after
the GitHub connector confirmed matching authenticated-user and repository-owner
IDs, with push permission, was also rejected. The reviewer requires trusted
user authorization naming the external destination and payload. No indirect
push or other upload was attempted. The branch remains local pending that
approval, and no merge was performed.

## Text for the shared status files

Proposed status text, pending review: S-16 stopped at the historical
checkpoint audit. The share-target weight calculator is implemented and
tested. Exact comparisons for S-13 and S-15 require noise precisions that
their reduced checkpoints did not save. No new study fits or external scores
were produced, and the working model remains unchanged. Review must decide
whether to authorize reconstruction fits or revise the historical check.

**Draft paragraph for the progress book and Chapter 4, pending review.**
We registered a study that sets survival weight through a target directional
share of loading precision. An audit found that the historical S-13 and
S-15 checkpoints lack the noise precisions required for the planned exact
comparison. We therefore stopped before fitting the proposed procedure.
Its convergence, null false activity and predictive performance have not
been evaluated.
