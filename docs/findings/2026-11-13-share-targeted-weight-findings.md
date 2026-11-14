# S-16: survival weight set by a target directional precision share

Date: 2026-11-13. Author: Andrew Walther.

Registration: [plan and Amendment 1](../plans/2026-11-13-share-targeted-survival-weight-plan.md).
Execution: [original prompt](../prompts/2026-11-13-codex-share-targeted-weight-prompt.md).
Branch: `s16-share-targeted-weight`. No external cohort was scored and no model was adopted.

## Summary

S-16 completed through Steps 2 to 5 under the registered fallback. No new weight is
selected for adoption. Target 0.01 failed the primary null gate: mean false-active
count was 0.2, but the mean false-active risk variance fraction was 0.19343, above
the 0.01 limit. No share-targeted arm entered real-data CV. The unit-weight and
two-step references both converged in all five folds and the full fit. Unit weight
remains the only registered selection candidate, with mean CV C 0.6604 versus
0.6563 for two-step, paired difference 0.00415 with SE 0.00415. Their loading
correlations are all above 0.9999935. No external scoring was run.

Amendment 1 resolved the original missing-checkpoint blocker. The weight calculator now
uses the gene-precision-weighted expression curvature and tracks the program selected
at the unit-weight pilot. The 50-sweep geometric ramp did not meet the registered
convergence gate. Under the plan's fallback, only target share 0.01 was carried into
the simulation screen. Targets 0.05, 0.1, 0.25 and 0.5 were dropped because their raw
real-data diagnostic fits reached 500 sweeps.

## Fixed design and evidence

Every S-16 fit uses at most 500 sweeps, with `cfg$cavi$max_iter <- 500L`, canonical
YFB projection, point-Laplace L and F, per-program Normal β and the survival-block
weight. Real training data have K = 5, 273 TCGA and CPTAC patients, 139 deaths and
2,218 D4 genes. Training cohorts are pooled without Cox strata. Each fold supplies
its own training center, applied to both training and held-out projections.

The two-step reference fits the same-prior, same-rank canonical model with w_S = 0
for the molecular loadings, followed by an exact-times Cox model on its live scores.
This operationalizes the unsupervised reference through canonical w = 0 CAVI plus
Cox rather than a separate EBMF software fit. The fixed choices specify w = 0
loadings, while Step 5 calls the reference EBMF. This implementation interpretation
is explicit, and no numerical equivalence to another EBMF backend is claimed. No alternative rank or prior is selected by this study.

Tracked outputs are aggregate CSVs under
`results/benchmark_sim/outputs/share_targeted_ws_2026-11-13/`, abbreviated below as
`OUT`. Full fits, posterior moments, final Cox working weights, fold predictions,
centers, frozen drivers and logs remain under ignored
`data/share_targeted_ws_2026-11-13/`, abbreviated as `PRIVATE`.

## Step 1: final-state share and pilot tracking

For program k, the exact registered ratio is

\[
r_k(w_S)=\frac{w_S E[\beta_k^2]\sum_i W_{ii}(Y\hat f_k)_i^2}
{\{\sum_i E[l_{ik}^2]\}\{\sum_j\tau_j\hat f_{jk}^2\}},
\qquad s_k=\frac{r_k}{1+r_k}.
\]

Here f̂_k is the unit vector in the mean loading direction, W contains the Cox
working weights evaluated at the final fitted predictor, and τ_j is gene j's noise
precision. The calculation is exact for these saved working-precision ingredients.
It is not a direct measurement of the Cox likelihood's complete information matrix
or of a change in loading orientation. Mean τ is used only in the closed-form
approximation.

The target program at the w_S = 1 pilot is the survival-active program with the
largest standardized absolute coefficient. `share_targeted_ws()` returns
w_S = w_current × [s*/(1−s*)] / r_k(w_current), capped at 100,000. If the pilot has
no active program, it returns 1 and flags the no-signal guard. Subsequent passes
match loadings to that same pilot vector by largest absolute correlation, whether
or not the matched program remains active. Loss of activity is recorded. There
are at most three target-fit passes, warm-started with all posterior moments and
persistent priors. Failure to reach the ratio tolerance after three converged
passes does not remove an arm from eligibility.

### Historical check under Amendment 1

| Source, program | Exact directional share | Per-gene median share | Directional approximation fold error |
|---|---:|---:|---:|
| S-9, 0.25% variance, seed 1, matched P3 | 0.683248 | 0.213322 | 1.273 |
| S-9, 1% variance, seed 1, matched P1 | 0.504894 | 0.050182 | 1.348 |
| S-9, 5% variance, seed 1, matched P1 | 0.200542 | 0.008695 | 2.213 |
| K = 5 S-11 P2, bound | 0.000122247 | 7.027 × 10⁻⁷ | 1.149 |
| K = 5 S-11 P4, bound | 0.000076192 | 1.164 × 10⁻⁶ | 1.301 |
| K = 5 S-11 P5, bound | 0.000153867 | 2.478 × 10⁻⁶ | 1.065 |

Source: `OUT/saved_share_diagnostics.csv`. K = 5 checkpoints lack E[L²], so their
E[L]² substitution gives upper bounds on the ratio and share, conditional on the
other saved ingredients. These bounds do not set new-study weights. S-13 and S-15
historical comparisons are waived. Every new S-16 fit retains τ, E[L²], E[F²], E[β²]
and final Cox weights.

The factor-two approximation check fails in the checked S-9 5% cell, with a
2.213-fold directional error. It is below four. The weight calculator uses the
exact expression, so execution continued without changing any threshold. The
preliminary numbers in Amendment 1 were based on the earlier mean-precision
calculation and are superseded by this amended exact audit. Final-state shares are
not compared with the old inside-sweep trace.

Rebuild the historical audit without fitting:

```sh
Rscript results/benchmark_sim/audit_share_targeted_ws.R
```

## Step 2: high-weight convergence diagnostics

### Raw fits and the ramp candidate

| Dataset | Raw weight | Raw sweeps, converged | Ramp sweeps, converged | Largest raw relative decrease | Largest post-ramp decrease |
|---|---:|---|---|---:|---:|
| S-9, p = 2,000, 1%, seed 1 | 10 | 237, yes | 241, yes | 0 | 0 |
| Same | 40 | 500, no | 500, no | 0.002123 | 0 |
| Same | 100 | 500, no | 500, no | 0.000247 | 0 |
| Same | 1,000 | 500, no | 500, no | 0 | 0 |
| Full real training, target 0.25 | 2,730.871 | 500, no | 500, no | 0 | 0 |
| Full real training, target 0.5 | 8,192.613 | 500, no | 500, no | 0 | 0 |

Source: `OUT/diagnostic_summary.csv` and per-case `*_trace.csv`. The raw S-9 w = 40
and 100 fits have early decreases larger than tolerance 10⁻⁵. This is consistent
with overshoot from a local quadratic approximation that is not a global lower
bound. The traces do not isolate Taylor error from empirical-prior updates.
The real high-weight fits also fail with monotone objectives, so avoiding decreases
is not sufficient.

Large coefficients and risk scores continue to change near the sweep cap. For
example, the raw real target 0.5 reaches max |η| = 94.89 and max |Eβ| = 25.91. The
finite 500-sweep record supports persistent drift and possible separation, but it
does not establish that the parameter sequence is mathematically unbounded.

The smallest tested change is a geometric continuation of the survival-block weight
from 1 to w* over R = 50 sweeps:

\[
w_t=\exp\{(t-1)\log(w^*)/(R-1)\},\quad t<R,
\qquad w_t=w^*,\quad t\ge R.
\]

At sweep t the objective remains
E[log p(Y|L,F,τ)] + w_t E[log PL] − KL_L − KL_F − w_t KL_β.
The coefficient update retains its A1 untempered form because the survival term
and coefficient KL have the same multiplier. The loading update and objective
use w_t. This is a continuation scheme, not a new lower-bound derivation or a
monotonicity guarantee. Objectives during the ramp are different and are not
compared for convergence. The first permitted stopping comparison is sweep 51,
when both objectives use w*.

The ramp removes the tested early decreases, but only 1/4 simulation cases and
0/2 real target cases converge. It therefore fails the registered joint convergence
and decrease criterion. No ramp length was tuned after these results. The ramp
remains opt-in and is not used in Steps 3 to 5. Its tests verify schedule, exact
unit-weight default compatibility and stopping restrictions. There is no successful
divergent-case recovery fixture, so that part of the proposed-fix test criterion
is not met.

### Raw-fit fallback

The new full-training unit-weight pilot converges in 52 sweeps. Its selected
program is P2. The lower target diagnostics are:

| Target | Pilot-derived weight | Sweeps | Converged | Tracked final share | Loading correlation to pilot |
|---|---:|---:|---|---:|---:|
| 0.01 | 82.75367 | 14 | yes | 0.00873483 | 0.999432 |
| 0.05 | 431.1902 | 500 | no | 0.947203 | 0.879260 |
| 0.1 | 910.2903 | 500 | no | 8.53 × 10⁻¹⁰ | 0.523009 |
| 0.25 | 2,730.871 | 500 | no | 2.29 × 10⁻⁸ | 0.537359 |
| 0.5 | 8,192.613 | 500 | no | 4.57 × 10⁻⁸ | 0.108225 |

Tracked shares for nonconverged states are diagnostic values, not successful
calibrations. At targets 0.1 to 0.5 the closest pilot match is P3. Reporting the
large share at the original P2 column would incorrectly substitute a different
program. Only target 0.01 enters the screen. This is the registered no-fix fallback.

Convergence throughout this report means the fitter's relative-ELBO stopping rule.
It does not imply that the predictor itself has stabilized to 10⁻⁴. The raw S-9
w = 10 fit has final relative η change 0.00194, and the real target 0.01 has 0.00234.
The complete traces and maximum coefficients are retained rather than hidden by
the convergence flag.

### Reproduction and driver failures

From the repository root, the original fit modes are:

```sh
Rscript results/benchmark_sim/run_share_targeted_diagnostics.R
Rscript results/benchmark_sim/run_share_targeted_diagnostics.R real
Rscript results/benchmark_sim/run_share_targeted_diagnostics.R ramp_sim
Rscript results/benchmark_sim/run_share_targeted_diagnostics.R ramp_real
Rscript results/benchmark_sim/run_share_targeted_diagnostics.R fallback_real
Rscript results/benchmark_sim/summarize_share_targeted_diagnostics.R
```

Fit modes protect existing checkpoints and stop if a destination already exists.
The last command rebuilds aggregate summaries from the saved 16 fits without
refitting. The diagnostic fit times sum to 1,285.988 seconds.

Both synthetic driver processes exited nonzero after writing all four requested
fit checkpoints because the live driver was edited while R was reading it. One
ended with a parse error and one attempted to rerun an existing real pilot.
No checkpoint was overwritten. Saved fields and all eight synthetic trace row
counts were checked independently. The errors are driver failures, distinct from
model nonconvergence. The transient raw traced source was not snapshotted before
editing. Later ramp/fallback hashes were retained, and the screen uses a frozen
driver with source fingerprints. A sandbox localhost-socket error occurred before
any screen fit and was resolved with local socket permission. The user-requested
usage pause interrupted unsaved in-progress tasks. Completed atomic tasks were
preserved and reused.

## Step 3: simulation screen and null gate

The screen completed 20 task checkpoints with 100 aggregate arm records: 30 primary
target records, 10 secondary paired-weight null records and 60 reference records.
There were 104 numerical fit calls: 60 pilots and two-step factorizations, 34 target
calibration passes and 10 secondary null fits. Six guarded target arms reused their
unit-weight pilots. All 100 arm checkpoint records retain the required finite
moments, τ and final working weights. This is not a count of 100 distinct fits.
No arm error or warning was recorded. There are 99 converged arm records and one
nonconverged capped null arm.

| Cell | Arm | Recovered / 10 | Mean false-active | Mean FA risk fraction | Converged / 10 | Mean validation C | Median slope |
|---|---|---:|---:|---:|---:|---:|---:|
| Weak | Unit pilot | 10 | 1.0 | 0.0256 | 10 | 0.6717 | 1.047 |
| Weak | Target 0.01 | 10 | 0.8 | 0.0198 | 10 | 0.6719 | 1.054 |
| Weak | Two-step | 10 | 0.5 | 0.0122 | 10 | 0.6682 | 1.039 |
| Low variance | Unit pilot | 5 | 3.0 | 0.5132 | 10 | 0.7995 | 1.579 |
| Low variance | Target 0.01 | 0 | 3.4 | 0.9826 | 10 | 0.7177 | 2.095 |
| Low variance | Two-step | 0 | 3.8 | 0.9463 | 10 | 0.6437 | 1.186 |
| Permuted weak | Unit pilot | 0 | 0.5 | 0.3877 | 10 | 0.4613 | −1.114 |
| Permuted weak | Own-pilot target 0.01 | 0 | 0.2 | 0.1934 | 9 | 0.4606 | −1.482 |
| Permuted weak | Two-step | 0 | 0.4 | 0.2256 | 10 | 0.4461 | −1.581 |
| Permuted weak | Paired signal weight | 0 | 0.4 | 0.2855 | 10 | 0.4602 | −1.115 |

Sources: `OUT/screen_metrics.csv`, `screen_null_gate.csv` and `screen_task_status.csv`.
Recovery includes all ten replicates. Under null survival there is no prognostic
truth, so all active programs are false. Mean FA risk fraction is the registered
operational nonnegative fraction defined below, not a percentage of active programs.

| Target-0.01 cell | Final weight range | Achieved directional share range | Mean per-gene share | Target passes | Guarded datasets | Largest standardized FA coefficient |
|---|---:|---:|---:|---:|---:|---:|
| Weak | 0.726 to 1.438 | 0.009719 to 0.010329 | 0.000561 | 10 | 0 | 0.1473 |
| Low variance | 0.004901 to 1.356 | 0.005497 to 0.014319 | 0.002224 | 18 | 0 | 1.3049 |
| Own-pilot null | 0.07664 to 100,000 | 0.000001878 to 0.010084 | 0.000750 | 6 | 6 | 28.5048 |

The null share summaries use only the four datasets with a defined pilot target.
The six guarded cases have no target direction, so their achieved share is undefined,
not zero. All ten weak and low-variance targets met the factor-two ratio tolerance.
Three of four defined null targets did so. Null targets lost activity in two
replicates, which were retained as Amendment 1 requires. Many simulation weights
are below 1, particularly in the low-variance cell, and were fitted as registered.

**Primary null verdict: fail.** Mean false-active 0.2 satisfies ≤1. The mean risk
fraction 0.193428 fails <0.01. Seed 2 converges with one false active carrying
fraction 0.934282. Seed 10 reaches weight 100,000 and 500 sweeps after two passes,
with one false active carrying approximately all of the risk variance. It remains
in the registered denominator. Excluding it only as a descriptive check leaves
fraction 0.103809 across nine datasets, still above 0.01. This does not replace
the registered ten-replicate calculation. The secondary paired-weight null also
fails the variance threshold, with mean false-active 0.4 and fraction 0.285512.
This gate describes ten replicates and does not establish a general false-positive
guarantee.

The target preserves weak-cell recovery but loses the low-variance recovery seen
at unit weight. Mean loading correlation of the best truth candidate with its
best two-step loading is 0.999979 for weak pilot and target, 0.474353 for the
low-variance pilot, and 0.494452 for the low-variance target. The two-step
self-comparisons are 1. The supplementary `screen_loading_agreement.csv` supplies
all 100 arm comparisons and flags 25 candidate records below the recovery correlation
threshold. These remain descriptive candidates, not recovered programs. Its 30
previously reported target comparisons agree with the raw metrics to 4.45 × 10⁻¹⁶.

The original fresh-screen commands, from the repository root, are below. They
require an empty private screen cache. Preserve historical caches rather than
overwriting them. Do not rerun fits merely to refresh a table.

```sh
Rscript -e 'dir.create("data/share_targeted_ws_2026-11-13", recursive=TRUE, showWarnings=FALSE); saveRDS(list(allowed_targets=.01,fit_control=list()), "data/share_targeted_ws_2026-11-13/screen_settings.rds")'
Rscript results/benchmark_sim/run_share_targeted_weight.R screen 2 data/share_targeted_ws_2026-11-13/screen_settings.rds
```

The observed screen used the frozen driver saved under PRIVATE and the source at
`c294866b`, unchanged through `3c18d39f`. Reusing those original task caches requires
that source state. A later change only fixes empty-target labels for real-mode
reporting, and intentionally makes the broad source fingerprint differ. It does
not change the screen algorithm or invalidate the historical results. Source
fingerprints matched at screen completion. The resumed batch took approximately
four minutes with two workers. Saved task times sum to 556.564 seconds, including
the pre-pause cached task.

The permuted-survival experiment permutes training time and event together.
Validation outcomes retain the planted weak signal, so its validation C and
calibration slope are not an independent validation null distribution. The
simulation preprocessing separately z-scores training and validation matrices,
following S-13. This is the inherited simulation protocol, not a prospective
training-only scaling procedure.

### Rebuild all-arm loading agreement without fitting

The following audit rebuilds `screen_loading_agreement.csv` from the saved states.
It regenerates synthetic truth, computes correlations and checks agreement with
the existing target rows. It does not refit any model. A convenience copy is saved
under `PRIVATE/execution_notes/s16-screen-loading-agreement.R`.

```r
# ============================================================
# Script: s16-screen-loading-agreement.R
# Purpose: Compare saved screen loadings with their two-step references.
# Author: Andrew Walther
# Created: 2026-11-13
# Dependencies: saved S-16 checkpoints and S-13/S-14 simulation generators
# ============================================================
Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
source('results/benchmark_sim/run_yfb_survival_weight_p2000.R'); s9_load()
source('results/benchmark_sim/run_yfb_supervision_strength.R')
source('results/benchmark_sim/run_yfb_near_threshold.R')
#' Compare the best truth candidate with its best two-step loading
#' @param fit,reference Saved fitted states.
#' @param truth Planted gene-direction vector.
#' @return One row with candidate indices and absolute correlations.
#' @examples
#' # loading_agreement(fit, reference, truth)
loading_agreement <- function(fit, reference, truth) {
  r <- abs(as.numeric(suppressWarnings(cor(truth, fit$EF))))
  if (!any(is.finite(r))) return(data.frame(truth_matchcor=NA_real_, bestprogramindex=NA_integer_, besttwostepindex=NA_integer_, correlation=NA_real_, candidate_below_0.7=NA))
  k <- which.max(replace(r,!is.finite(r),-Inf))
  rr <- abs(as.numeric(suppressWarnings(cor(fit$EF[,k], reference$EF))))
  j <- if (any(is.finite(rr))) which.max(replace(rr,!is.finite(rr),-Inf)) else NA_integer_
  data.frame(truth_matchcor=r[k], bestprogramindex=k, besttwostepindex=j,
    correlation=if(is.na(j)) NA_real_ else rr[j], candidate_below_0.7=r[k]<.7)
}
rows <- list()
for (cell in c('weak','lowvariance')) for (seed in 1:10) {
  ck <- readRDS(sprintf('data/share_targeted_ws_2026-11-13/screen/%s_s%02d.rds',cell,seed))$result
  truth <- if(cell=='weak') s13_data(seed,'weak',FALSE)$truth else s14_data(seed,.0005,'strong')$truth
  sets <- if(cell=='weak') c('signal','null') else 'lowvariance'
  for (nm in sets) {
    set <- ck[[nm]]; label <- if(nm=='signal') 'weak' else if(nm=='null') 'weak_null' else 'lowvariance'
    ref <- set$twostep$checkpoint$fit
    arms <- list(pilot=set$pilot$checkpoint$fit, twostep=ref,
                 primary=set$targets[['0.01']]$checkpoint$fit)
    if(nm=='null') arms$paired_signal_weight <- ck$paired[['0.01']]$checkpoint$fit
    for (kind in names(arms)) {
      r <- loading_agreement(arms[[kind]],ref,truth)
      rows[[length(rows)+1L]] <- cbind(data.frame(cell=label,seed=seed,
        arm=if(kind %in% c('primary','paired_signal_weight')) 's0.01' else kind,
        analysis=kind,target=if(kind %in% c('primary','paired_signal_weight')) .01 else NA_real_),r)
    }
  }
}
out <- do.call(rbind,rows)
raw <- read.csv('results/benchmark_sim/outputs/share_targeted_ws_2026-11-13/screen_metrics.csv')
join <- merge(out,raw,by=c('cell','seed','arm','analysis','target'))
valid <- is.finite(join$matched_twostep_correlation)
stopifnot(nrow(out)==100L,nrow(join)==100L,all(abs(join$correlation[valid]-join$matched_twostep_correlation[valid])<1e-10))
write.csv(out,'results/benchmark_sim/outputs/share_targeted_ws_2026-11-13/screen_loading_agreement.csv',row.names=FALSE)
cat('rows',nrow(out),'existing comparisons',sum(valid),'max difference',if(any(valid)) max(abs(join$correlation[valid]-join$matched_twostep_correlation[valid])) else NA_real_,'\n')
print(aggregate(correlation~cell+analysis,out,mean))
cat('below .7',sum(out$candidate_below_0.7,na.rm=TRUE),'\n')
```

## Step 4: training cross-validation

No share target passed the primary null gate. The real evaluation therefore
contains only the unit-weight and two-step references. Both use the same five
S-11 folds and fold-specific training centers. The twelve full/fold fit records
have no errors or warnings, correct centers and complete final-state ingredients.
Every fit converged, with at most 152 sweeps.
CV is conditional on the fixed prepared D4 feature matrix. Feature selection and
per-cohort standardization were not fully redone inside folds, so this does not
evaluate the entire prospective preprocessing pipeline.

| Arm | Converged folds | Full converged | Mean fold CVPL | Mean CV C | SE C |
|---|---:|---|---:|---:|---:|
| Unit weight | 5/5 | yes | −154.697688 | 0.660413 | 0.027184 |
| Two-step reference | 5/5 | yes | −154.713830 | 0.656264 | 0.028765 |

Sources: `real_fold_metrics.csv`, `real_summary.csv`, `real_selection.csv` and
`real_paired_c.csv`. The paired C difference, unit weight minus two-step, is
0.00414886 with SE 0.00415221. The mean CVPL difference is 0.01614185.
Two-step is a comparison reference rather than a candidate in the registered
unit-weight-versus-target selection rule. Unit weight is the sole surviving
candidate. Its retention is not evidence of a predictive improvement from
share-targeted supervision. Both references satisfy the literal convergence
criterion and the implementation's additional complete-score check.

The full unit-weight pilot's selected P2 has exact directional share 0.0001220463
and per-gene median share 7.01582 × 10⁻⁷. The selected-program fold shares range
from 0.00007644 to 0.00016121. All two-step loading-supervision shares are zero
because w = 0. These values use saved second moments, not the historical bound.

```sh
Rscript -e 'saveRDS(list(allowed_targets=numeric(0),fit_control=list()), "data/share_targeted_ws_2026-11-13/real_settings.rds")'
Rscript results/benchmark_sim/run_share_targeted_weight.R real 2 data/share_targeted_ws_2026-11-13/real_settings.rds
```

The observed run used one PSOCK worker while the full regression suite occupied
the second CPU slot. It completed in approximately 2 minutes 37 seconds. Its
frozen real driver and current source fingerprints match commit `883505dd`.
The empty-target label bug was caught before this run, reproduced in a failing
test and corrected. No phantom target rows or target fits were produced.

## Step 5: program structure and annotation

All five one-to-one unit-weight loading matches with the same-rank two-step
reference have absolute correlation 0.999993593 to 0.999999663. The full K × K
matrices are in `structure_cor_w1.csv` and `structure_cor_twostep.csv`, with
one-to-one matches and fold correspondence in the adjacent match/folds tables.
No share-targeted full fit was eligible for this step.

| Unit-weight program | Active | Main program | Active fold matches | Δ log partial likelihood | Semi-partial expression R² |
|---|---|---|---:|---:|---:|
| P1 | no | no | 0/5 | 0.28745 | 0.151153 |
| P2 | yes | yes | 5/5 | 9.52637 | 0.114596 |
| P3 | no | no | 2/5 | 1.50866 | 0.069406 |
| P4 | yes | no | 3/5 | 3.43704 | 0.053510 |
| P5 | yes | no | 2/5 | 4.53023 | 0.040484 |

The nested-Cox likelihood differences refit all K coefficients and then refit
without each program, with `timefix = FALSE`. They do not merely set a coefficient
to zero. Semi-partial R² is the registered increment in reconstruction sum of
squares explained conditional on the other fitted programs, divided by total
centered expression sum of squares.

Both arms have active P2, P4 and P5. Only P2 satisfies the unit-weight main-program
rule. The two-step reference has main P2 and P5. P5 has five active fold matches
under two-step and two under unit weight. This is an activity/stability difference
under different coefficient estimation, despite nearly identical molecular
loadings. It is not evidence that unit-weight supervision discovered a new
molecular direction. No fold match was unknown.

Annotation used fixed Moffitt Classical/basal-like and human Elyada iCAF/myCAF
sets from the existing local resources. Measured overlaps are respectively 21,
18, 25 and 15 genes. All four signatures were evaluated in all five programs of
each arm, producing 40 rows. Enrichment uses absolute loading ranks, so it
summarizes loading magnitude rather than subtype activation or a signed pole.
For unit weight, P3 has Classical NES 1.851, adjusted p 4.47 × 10⁻⁵, P2 has
basal-like NES 1.554, adjusted p 0.0481, and P5 has myCAF NES 1.651, adjusted p
0.00563. Two-step has the same descriptive patterns. Adjusted p-values are from
the existing per-program four-signature analysis, not a new multiplicity procedure
across all programs and both arms. The full table includes every signature,
including nonsignificant results.

```sh
Rscript results/benchmark_sim/run_share_targeted_weight.R structure 1 data/share_targeted_ws_2026-11-13/real_settings.rds
```

`structure_status.csv` reports no errors for either arm. Annotation and overlap
outputs are `structure_annotation.csv` and `structure_signature_overlap.csv`.

## Implementation and verification

`code/share_targeted_ws.R` accepts a converged canonical training fit, its n × p
training matrix, outcomes, desired share and optional fixed pilot loading vector.
It validates dimensions and moments, calculates per-program final-state diagnostics,
selects or tracks the target and returns a weight plus cap, guard, match and
activity flags. Its exact ratio uses individual gene precisions. Invalid states
fail explicitly. It uses base R and the fitter's existing Cox working-weight helper.

`trace_s16 = TRUE` adds per-sweep, per-program objective, η-change, coefficient,
loading-norm and final-state share records. `survival_weight_ramp = 0` is off by
default. A positive value enables the documented geometric schedule only for the
supported pooled canonical A1 setting. Unsupported model controls are rejected.
The molecular and survival model defaults remain unchanged.

`run_share_targeted_weight.R` performs dataset-level pilot fitting, up to three
posterior-preserving recalibrations, exact-times two-step scoring, primary and
paired-weight nulls, fold-centered real CV and checkpoint validation. Successful
arms survive failures in other arms. It saves complete states and warns explicitly.
The structure module computes nested-Cox likelihood differences, reconstruction
variance increments, one-to-one loading matches and the S-11 program rule. Its
annotation uses the existing fgsea helper and fixed local tumor/fibroblast gene
sets. No new dependency was introduced.

The tests check hand-calculated heterogeneous-precision shares and inversion,
activity and guard behavior, sign/reordered fixed-pilot tracking, loss of activity,
cap transitions, warm starts, failure preservation, fold centering, CV scores,
selection ties, structure quantities and synthetic signature membership.
The pre-change w_S = 1 fixture checks exact full-object equality, but its tiny
example has zero final factors and supplies narrow compatibility coverage.

```sh
Rscript tests/test_share_targeted_ws.R
Rscript tests/test_share_targeted_diagnostics.R
Rscript tests/test_share_targeted_runner.R
Rscript tests/test_share_targeted_structure.R
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 Rscript tests/run_tests.R
Rscript derivations/SSBMF-YFB-model-derivations/verify_derivation.R
```

The full suite passed 879 counted checks with zero failures. New `stopifnot()`
assertions also passed but are not included in that count. The existing bounded-β
optimizer profile-likelihood fallback warning was surfaced. The derivation
verifier passed 109/109. Code review cleared the raw 0.01 workflow after fixes to
reference routing, final-state reporting, partial-failure retention, dependencies
and source validation. A separate fresh independent statistical review passed
the final reporting of the fallback evaluation with no blocking issue. It checked
the amended plan, code, provenance, both reports and aggregate tables. It
independently recomputed the null gate, CV means, SEs, paired differences and
loading-agreement summaries from CSVs. It did not rerun fits, the full test suite,
the derivation verifier or all R audit calculations. The review approves reporting,
not adoption, merging or external scoring. Its fixed-D4 CV clarification is
incorporated above. The full review is preserved under PRIVATE/execution_notes/.

## Judgment calls, deviations and limits

- The Step 5 term EBMF is operationalized using the fixed-choice w = 0 canonical
  factorization plus Cox, not a separately fitted EBMF backend. No equivalence to
  a separate software implementation was tested.
- Three passes means three target fits after the unit-weight pilot. No target is
  removed just for missing the ratio tolerance after three converged passes.
- The primary false-active variance fraction is Var(sum of false-active risk
  components) / Var(total risk), averaged across replicates. The S-11 signed
  covariance allocation is also reported. This operational definition was fixed
  before null results. An undefined fraction with false actives cannot pass.
- Under permuted survival, all active programs are false. Under signal, S-9's
  |r| ≥ 0.7 match is required before excluding a program from false activity.
  Earlier S-14 scoring always excluded the best match and is not used here.
- CVPL uses the mean fold contribution, so the 0.5-unit tie window is on that
  scale. Fold-varying supervision is ordered by median final fold weight, then
  target. These definitions were fixed before CV results.
- The implementation requires complete finite CVPL and C in all five folds in
  addition to the registered convergence rule. This is stricter than the literal
  eligibility requirement. Its practical effect is reported with the results.
- Cap duplicates are explicit unfit provenance records. General multi-target cap
  alias resolution is limited, and this single-target fallback does not test it.
- No empirical finding licenses a new model, an altered null threshold, an
  alternative rank or external model selection. Scientific adoption waits for
  Andrew after review. Optional Step 7 was not run. No target passed the primary null gate, so no
  wider share-target selection or rank surface was added to this failed candidate
  evaluation. Rank K = 5 remains fixed by the prior adopted configuration.

## Text for the shared status files

**Proposed status text, pending review and Andrew's adoption decision.** S-16
completed under the registered no-fix fallback. The 50-sweep ramp did not make
large-weight fits converge within 500 sweeps. Only target share 0.01 entered the
screen, and it failed the primary null risk-variance gate. No share-targeted arm
advanced to real data. The unit-weight working model remains unchanged. Its new
five-fold reference C is 0.6604 versus 0.6563 for two-step, paired difference
0.00415 with SE 0.00415. The molecular loading correlations exceed 0.9999935.
No external scoring was run.

**Draft paragraph for the progress book and Chapter 4.** We evaluated survival
weights defined by a target directional share of the loading update's working
precision. A geometric ramp removed early objective decreases but did not resolve
large-weight nonconvergence within 500 sweeps. The sole target allowed by the
registered fallback, 0.01, retained recovery in the weak-survival cell but reduced
low-variance recovery from 5/10 at unit weight to 0/10. Its primary permuted-survival
null had mean false-active count 0.2, but false-active programs accounted for a
mean 19.34% of fitted risk variance, above the 1% limit. No new target proceeded
to real-data evaluation. The unit-weight and two-step reference loadings remained
nearly identical, with a small training-CV concordance difference of 0.00415 and
SE 0.00415. These results do not support adopting share-targeted weighting.
