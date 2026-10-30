# Normalized-projection supervision: precision and numerical diagnosis

The approved specification is the October 30 realignment plan. This document
records diagnoses during execution; the saved gate tables and final execution
report determine the outcome. The external validation cohorts were not opened
for model choice or scored by these experiments.

## Precision balance

`run_yfb_precision_audit.R` loads only the TCGA/CPTAC training cohorts and the
saved single-log D4 fit. It aligns the saved genes, divides the nonnegative raw
expression by cohort-specific SD, and centers the projection. Saved activity
moments and loading moments are canonically rescaled, and saved coefficients
are expressed in normalized-projection units. This is a conditional precision
audit at a transformed saved state, not a refit under the new preprocessing.
Saved residual precisions are retained; it does not claim a new noise fit.

For gene j and factor k the reconstruction precision is τ_j Σ_i E[L_ik²],
and survival precision is γ E[β_k²] s_k² Σ_i w_i Yc_ij². Direct per-gene sums
are checked against the matrix identity. The median share is monotone in γ;
the γ for a target share follows from the median precision odds, with direct
root solving for the simulation grid. Dead columns have no defined relative
precision and are excluded from the revised calibration.

Outputs are under `results/benchmark_sim/outputs/yfb_supervision_2026-10-30/`:
`precision_audit.csv` and `sign_variance_audit.csv`. The patient-level cache
`training_only.rds` is explicitly git-ignored and must never be staged.

## Comparison with DeSurv

The installed package is version 1.0.1, revision
`cc9d4c2486b3125704ab71eea61f35e040a40080`. Inspection of the R wrappers and
that exact revision's `src/functions.cpp` shows a loss with reconstruction
SSE/(2‖X‖²), survival contribution 2 log PL/N_event, penalties on standardized
coefficients and normalized factor norms. Its loading step matches gradient
norms, clamps multiplicative ratios, and backtracks against its actual loss.
This differs from multiplying the unnormalized Bayesian survival likelihood
by feature count divided by event count. Its existence does not prove the
specified variational updates converge at large γ.

[Exact source](https://github.com/rashidlab/DeSurv/blob/cc9d4c2486b3125704ab71eea61f35e040a40080/src/functions.cpp).

## Registered execution failures and corrections

- The sandbox initially disallowed a local socket; no fit ran in that launch.
- The first registered run included dead columns in the calibration median.
  Their zero reconstruction precision made the calibration undefined. Saved
  checkpoints and task errors remain; missing arms are explicit failures.
- The second run excludes dead columns, but strong messages expose numerical
  disagreement between the compiled posterior and ebnm. A saved pseudo-message
  shows loss of the slab variance in ebnm; a second message at very large
  precision yields a negative KL due to numerical error.
- A comparison defect allowed the nominal unsupervised reference to use a
  survival-score start. The corrected reference disables that start.
- The third registration rescales normal-means problems, checks both moments,
  reduces precision damping, and saves coherent states and traces on failure.
  Smaller damping alone did not repair the reproduced case; all tested values
  failed by a posterior, KL or prior-optimizer diagnostic.

Every registration precedes its corresponding run in the attempts log.
Do not combine these attempts or exclude their failed records from a denominator.

## Mathematical limits of the proposed objective

The diagonal Cox working curvature is an approximation. Recomputing the loading
normalizer changes the projection between sweeps, while the update treats that
normalizer as fixed. The damped F message also differs from the full-weight β
message. Therefore the displayed objective is monitored rather than advertised
as a monotone exact ELBO. Decreases and failure sweeps are retained. The code
does not label a numerical failure as convergence or substitute a silently
corrected posterior when the two numerical calculations disagree.

## Stable posterior repair

Independent quadrature identifies the source of the discrepancy, not just its
symptom. Canceling N(x;0,s²) gives a stable spike/slab normalizer and the exact
KL identity in the decision log. A new isolated kernel evaluates this same
posterior; the shared historical kernel is unchanged. The saved weak-message
fixture and moderate messages agree with independent integration. A warned
profile-likelihood fallback solves a separately reproduced constrained-prior
optimizer failure without changing the parameter bounds.

The fourth registered run evaluates this repair. Earlier normal-means errors
remain in their original attempts. Diagnostic replay of the previously failing
seed now converges at every tested damping value but still produces several
false survival-active factors. This separates numerical repair from the
unresolved statistical gate.

## Precision trace convention

The original aggregate sweep medians included dead factors, whose projection
scale is an arbitrary floor. They must not be interpreted as a meaningful
loading balance. Updated traces exclude dead columns and record factor-wise
shares, coefficients and norms. Derived final-state matched-factor and
live-factor shares are recomputed from the saved fit and regenerated paired
training data. These are explicitly final-state evaluations, not a replacement
for the original sweep traces. No recovery or gate metric depends on this
reporting correction, and no original checkpoint is overwritten.


## Subsequent boundary and variance-storage corrections

Fixed-seed inspection identified exactly zero activities with nonzero survival
loadings, and tiny positive activity moments skipped by an arbitrary RMS
threshold. Canonicalize every positive RMS; normalize a surviving loading when
its activity is exactly zero and carry its prior scale. Tests preserve genomic
reconstruction, projection moments, coefficient prior and KL. This addresses
an unidentified scale without changing the statistical model.

The next run then exposed lost coefficient uncertainty in large-mean second
moments. Stop that invalid driver and retain all completed and missing records.
The final registered run stores variances separately and evaluates centered
Normal-mixture posteriors/KL by conjugate and mixed-support formulas. It also
computes expected reconstruction and risk uncertainty as nonnegative variance
sums. This is the same posterior model with a stable numerical representation.
The variance-storage audit distinguishes raw subtraction failures from the
separately carried variances. No variance-loss state is accepted as a valid
gate fit.

## Shared-loading simulation coordinates

A further algebra audit found a simulation design error introduced in this
execution: raw expression used F, prepared reconstruction used F/SD_train,
but generated risk used raw-expression F and the recovery target was F×SD.
Thus attempts 1–6 cannot establish recovery of the shared reconstruction/Cox
loading. Their numerical failures and repairs remain informative, but their
negative gate tables are not decisive scientific tests of that claim.

Attempt 7 generates risk from the unit-normalized prepared reconstruction
loading F/SD_train. Validation retains that same training direction, while its
features receive their own cohort SD division. The new test checks both
rank-one reconstruction and Cox direction algebraically. A deterministic
comparison confirms that the legacy generator still reproduces old matrices,
outcomes and recovery targets exactly. Raw and prepared variance shares and
oracle C are saved separately. The coefficient is now per unit normalized
prepared projection; across-attempt discrimination is not directly comparable.
The complete corrected grid retains the originally registered gate.
