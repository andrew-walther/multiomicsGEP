# Feature-group YFB fitter and registered simulation

Author: Andrew Walther

Date: 2026-11-12

The [informative-scenario study](#2026-11-12-informative-scenario-study)
passes the expression-only gate in one cell. Both joint priors pass there,
so the registered rule supports neither a benefit nor a disadvantage from
separate priors. Andrew separately authorized the matched-TCGA comparison
as **not validated by simulation**. Earlier runs are retained below as
historical evidence, including their original stop decisions.

## First run, before the empty-block fix

The registered 80-fit simulation is complete. **Every cell and arm fails the
acceptance gate.** The grouped model has zero recoveries, 38 explicit
refitted-prior support failures and two converged fits. All 40 pooled fits
converge, but their best cell recovers only 4/10. The rule leaves grouped as
the default preference because pooled passes no cell. It permits neither
prior to proceed to step 4 because grouped fails every cell. The descriptive
matched-TCGA comparison did not run. These results do not establish that
separate modality priors improve recovery or prognosis.

The full suite passes **844 tests with zero failures**, including 14 fitter
checks and 11 simulation checks. It emits the existing bounded-β optimizer
profile-likelihood fallback warning. The study itself emitted **zero warning
messages**. All returned fits converged. The 38 errored grouped attempts remain
failed and nonconverged in the denominator, rather than being treated as
iteration-cap fits or silently removed.

The feature-group fitter was implemented in `6f31a0d7` after tests in
`8a3b5517`. The initial runner was committed in `ae5da859`. The finish prompt
resolved the earlier recovery-target ambiguity. This continuation committed
its failing tests in `d96df715`, the tested runner and smoke aggregates in
`2c672530`, and the full simulation aggregates in `e065da5c`. Main was pulled
with `--ff-only`. No worktree or branch change was used.

NULL and one-group fits still preserve EL, EF, EBeta and the ELBO history
bit for bit on pooled and stratified-centered fixtures. The continuation does
not change the fitter, its update modules or the historical multimodal fitter.
The earlier derivation verification passed 88 checks. No algebra changed here,
and that verification was not rerun.

Run the focused regression checks from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript tests/test_feature_group_simulation.R
```

## Test implementation and interpretation

**Purpose.** The tests specify the feature-group extension while preserving the
existing single-modality numerical path. The frozen outputs in
`tests/fixtures/yfb-feature-group-reference.rds` were generated from the unchanged
fitter at `19f910bd`. They contain synthetic data results only.

**Logic.** `fg_fixture()` builds a 48-patient, 24-feature signed example with
two programs and censored survival. Its second version adds two strata and
stratum-specific feature means. `fg_fit()` resets the fitting seed and applies
the working point-Laplace L/F and per-program Normal β configuration for five
sweeps. Tests require `identical()` for EL, EF, EBeta and the ELBO history under
the omitted argument, explicit NULL and one group.

The noncanonical F-update check uses 40 fixed pseudo-observations whose two
groups differ strongly in sparsity. It compares each block's posterior means,
second moments and fitted prior with a separate `ebnm()` call. A pooled prior
cannot satisfy both reference results.

`fg_observe()` temporarily wraps the real canonical F update, KL evaluation
and rescaling functions. It delegates their calculations unchanged, records
their inputs and outputs, and restores the original bindings on exit. The
canonical tests require one update per group and program, the correct input
columns, a separate prior refit on each block's pseudo-observations, and
rescaling of every live group's prior. The implemented result contract is
`fit$prior_F[[group]][[program]]` for explicit multiple groups. At least one
prior must remain live so the rescaling check cannot pass vacuously.

`fg_loading_kl()` recomputes each block's loading KL under the returned group
prior. The sweep's pseudo-likelihood and incoming prior define q. If the prior
refit changes g, the check adds E_q[log(g_old/g_new)] to the original KL.
Numerical integration supplies E_q|F| for the Laplace slab, and the atom at zero
is treated separately. An incompatible returned point mass produces infinite
KL. A separate test integrates q log(q/g) directly and checks the same answer,
including the canonical change of units F → F/c and rate → c rate.

The objective check combines this independently computed F KL with the captured
L and β KL contributions and recomputes the canonical objective from the
returned posterior moments. Prediction is checked against the sum of the two
new-data block projections, with unit canonical normalization asserted.

**Design choices.** Frozen numerical outputs detect changes even when NULL and
one-group fits change together. The observer checks real update calculations
rather than substituting artificial returns. The loading KL check uses the
returned refitted prior, because the existing helper's cached KL refers to its
incoming prior. These can differ. The grouped path corrects that accounting
without changing the required single-group baseline.

**Inputs and outputs.** Fixtures return a patient-by-feature matrix, survival
vectors and fitter arguments. The observer returns the fit and captured update
records. The KL helper returns one scalar per block. The standalone test file
prints named assertions and returns a nonzero process status for any failure.

**Limits.** The objective check independently verifies F accounting but uses
the existing canonical objective helper for the likelihood terms and captured
production KL values for L and β. It runs after one and four sweeps. Additional
analytic checks cover Normal and point-exponential KL corrections. Grouped
restart and pruning are explicitly rejected. Centered grouped prediction is
not separately tested. Bitwise frozen fixtures depend on R, BLAS and package versions. Do not
regenerate them with a fitter that already implements feature groups.

**R details.** `ebnm()` supplies the independent block-prior reference.
Second moments are reconstructed as posterior mean squared plus posterior SD
squared, because its default posterior output has no `second_moment` field.
Base R `integrate()` provides the independent loading-prior density calculation.
Observer bindings use `match.call()` so named helper arguments are supported.

## Fitter implementation

`code/feature_groups_yfb.R` contains four helpers. They add no package dependency.

`yf_feature_groups(feature_group, p)` validates one nonmissing, nonempty label
per feature and returns column indices in first-appearance order. NULL and
single-group labels return NULL, which selects the unchanged arithmetic path.
Factor and numeric labels are converted to character. This keeps the existing
model intact when there is only one block and avoids changing positional callers.

`yf_group_ebnm(x, s, groups, family)` solves each group's Gaussian
pseudo-observation problem with `ebnm()`. It returns combined posterior means,
SDs and second moments, each group's full fitted prior and each negative KL.
Both `update_F_k()` and `update_F_surv_YFB_k()` use it for explicit multiple
groups. Their likelihood precisions and signals are unchanged. The former is
the reconstruction-only update. The latter also retains its existing
noncanonical survival terms. Groups may occupy nonadjacent matrix columns.

`yf_update_F_groups()` implements the canonical grouped F step. For each group
and program it initializes a prior from that group's reconstruction
pseudo-observations, then calls the existing sequential loading sweep on that
block's columns. The sweep receives the full current Z_k and V^Z_k, including
all other groups. Updated joint moments pass to the next group. The helper
assembles the loading moments, A/B diagnostics and summed loading KL. Its
dead-block threshold is the existing 10⁻⁸ likelihood-precision threshold.
A dead block is zeroed and removed from the joint projection before continuing.

`yf_group_loading_kl()` evaluates the sweep posterior under the refitted group
prior. Normal priors need only second moments. Exponential slabs use the
posterior nonzero probability and mean. Laplace slabs use the existing analytic
absolute-moment helper. The density correction is
KL(q||g_new) = KL(q||g_old) + E_q[log(g_old/g_new)]. An incompatible prior raises
an explicit error. No clipping or replacement prior hides that failure.

`fit_ssbmf_yfb()` accepts the appended `feature_group` argument. Shared L,
feature-specific τ, Cox scores and coefficient updates remain in the existing
loop. Grouped priors persist across sweeps. `canonicalize_yf_factors()` scales
every group's prior along with its F column. Canonical Normal variances divide
by c², and exponential/Laplace rates multiply by c. The fitter returns
`prior_F[[group]][[program]]` and the column labels for multiple groups.
Canonical priors use the existing family/pi/rate/variance representation.
Noncanonical priors are the `ebnm()` fitted-g objects.

The grouped canonical objective contains the corrected F KL. **The historical
single-group path still caches KL under the incoming sweep prior rather than
the refitted prior.** Preserving bitwise single-group compatibility requires
leaving it unchanged here. Grouping can therefore change the monitored
objective and stopping time through both separate priors and this accounting
correction. Absolute objectives should not be compared across these paths as
though only the number of priors changed.

Multiple groups reject pruning, posterior restarts, frozen F and fixed-q EM
prior updates. Any explicit feature_group rejects cohort offsets. These
combinations require additional state and validation and are unnecessary for
the fixed study. Unsupported combinations fail before fitting. Existing
single-group options otherwise retain their historical behavior.

Inputs remain an n-by-p numeric matrix, survival vectors and existing fitter
controls. Outputs remain the original posterior moments and history, with
group-prior metadata added only for multiple groups. Prediction needs no new
formula: concatenate new blocks in training order and use
`predict_ssbmf_yfb()` with the stored normalization and, when requested,
training projection center. Training-only preprocessing is the caller's
responsibility.

## Registered simulation implementation

Andrew resolved the recovery target in the
[finish prompt](../prompts/2026-11-12-codex-d-finish-prompt.md). The experiment
uses the standardized reconstruction loading g = D⁻¹F[,3], stacked in fitted
column order and normalized to unit Euclidean norm. Both survival and recovery
now use g. Raw-F and DF correlations are not study outputs.

The four cells cross 50:50 or 50:250 expression:methylation features with
methylation active fractions 0.15 or 0.03. Expression retains 0.15. Each seed
has grouped and pooled F-prior arms, giving 80 fits over seeds 1–10. The
working configuration is K = 5, point-Laplace L and F, per-program Normal β,
w_S = 1, canonicalization, pooled cohorts without strata and no pruning.
Stopping uses max_iter = 300 and tol = 10⁻⁵ from the global configuration.
The generator remains at 120 training and 120 validation patients, noise scale
1 and truncation at zero. These differ from the stronger-signal S-9 cells.

`simulate_multimodal_yfb_data()` adds the optional
`methylation_active_fraction`. NULL uses the existing loading-generator default
and preserves the serialized pre-change output exactly. A supplied fraction
changes methylation alone. The existing generator selects ceiling(p × fraction)
nonzero entries per program, with at least one. Thus the nominal 15% and 3%
fractions realize as 16% and 4% at 50 features, and 15.2% and 3.2% at 250.
The new argument rejects values outside (0,1]. No package dependency changed.

`fgsim_prepare()` takes the generator output and seed. It returns training and
validation matrices, their training-only scaling statistics, g, a, both risk
vectors and redrawn outcomes. It uses `fgsim_scale()` to standardize each
block, subtracts the standardized training centers for both projections and
sets a = SD(original training η) / SD(X_c g). Survival uses the generator's
Weibull shape 1.5 and scale 0.01. The fixed survival seed is seed + 12000.
Training draws precede validation draws. Each cohort uses the generator's own
censoring calibration, with target 0.30. Validation never estimates centers,
feature SDs, g or a. A nonfinite or nonpositive multiplier stops the run.

`fgsim_assert_inputs()` checks identical matrices, survival vectors and common
signed-SVD initial L/F matrices before and after both arms. The matrices are
computed once per cell and seed. `fgsim_fit()` passes those inputs and the
fixed controls to the existing fitter. Its local child environment wraps
canonicalization only to observe the final pooled priors, which the historical
fitter does not return. The wrapper delegates all calculations unchanged and
leaves global bindings unchanged. A two-sweep regression test requires the
entire pooled fit to be identical to a direct call. `fgsim_assert_priors()`
requires expression and methylation entries, each with K valid fitted
point-Laplace priors. Empty placeholders fail.

`fgsim_recovery()` matches the column maximizing |cor(g,F_k)|. It requires a
correlation of at least 0.7, posterior survival activity and β_k cor(g,F_k) > 0.
Activity uses |E[β_k]| ≥ 1.96 SD(β_k), with positive SD and VBeta when available.
The second-moment fallback applies only to fit objects without VBeta. A wrong-sign
match is not recovered and remains in the false-active count. The returned
program table contains signed correlations and coefficient uncertainty. A
separate block table gives each modality's signed matched-program correlation.
Dead columns have undefined correlations.

`fgsim_precision()` retains Σ_j τ_mj E[F²_mjk] for each modality and program.
`fgsim_one()` adds feature shares, live-program flags and methylation dominance
indicators. Its prior table records point-mass weight 1 − π, Laplace scale
1/rate and the true nonzero fraction. Here π in the fitter is the slab
probability. A point-mass prior has undefined slab scale. Pooled-prior parameters
are repeated by modality for comparison and explicitly marked `prior_scope =
pooled`. Full moments, inputs and generated data remain under ignored `data/`.
Only aggregate summaries leave that directory.

`fgsim_one()` runs the two arms sequentially, captures every warning occurrence
and retains errors as failed rows. `fgsim_run()` writes the aggregate tables
after each seed pair and refuses to overwrite study directories. It uses one
CPU worker with BLAS thread counts set to one by the launch command. Its
`--smoke` mode uses seed 1 in each cell, giving eight fits. Smoke rows cannot
pass the ten-seed gate.

`fgsim_gate()` reports recovery, mean false activity, convergence and mean
validation C with its defined count. Failed seeds remain in the denominator.
The registered thresholds are recovery ≥ 8/10, mean false activity ≤ 1 and
convergence ≥ 9/10. The implementation conservatively fails any cell containing
an error because false activity is unknown for that fit. This explicit
no-error condition is additional to the three numerical thresholds.
`fgsim_select()` applies the registered preference for grouped priors unless
pooled passes a cell that grouped fails. If grouped fails all cells, step 4
stops regardless of pooled results. No weights or prior-support changes are
introduced.

The test commit is `d96df715`. Before implementation, the two retained checks
passed and nine new checks failed. Tests freeze the NULL generator behavior,
verify exact active counts, reconstruct the Weibull and censoring draws,
check the signed recovery boundary and VBeta preference, and reject altered
arm inputs. The pooled observer and gate-selection checks use real fitter
output and explicit synthetic decision tables, respectively.

### Smoke check

The eight-fit smoke check completed. Each grouped fit stopped with
“Nonfinite feature_group F KL: the refitted prior is incompatible with its
sweep posterior.” All four pooled fits converged, and none recovered the
planted direction. There were no warning messages. These are smoke outcomes,
not the registered ten-seed gate. The error was retained without changing
settings or prior support. The registered ten-seed results follow below.

## Simulation results and adoption rule

The [aggregate output directory](../../results/benchmark_sim/outputs/feature_group_simulation_2026-11-12/)
contains `metrics.csv`, `summary.csv`, `selection.csv`, `programs.csv`,
`blocks.csv`, `priors.csv`, `precision.csv`, `settings.txt` and source hashes.
Each cell has ten attempts per arm, including every error. C is oriented toward
greater hazard without validation-outcome sign changes. Undefined C values
are excluded only from its descriptive mean, with the defined count shown.

| Expression:methylation | Nominal methylation nonzero | Arm | Recovered / 10 | Mean false active | Converged / 10 | Errors | Mean validation C (defined n) | Gate |
|---|---:|---|---:|---:|---:|---:|---:|---|
| 50:50 | 15% | grouped | 0 | unknown | 1 | 9 | 0.492 (1) | fail |
| 50:50 | 15% | pooled | 2 | 0.8 | 10 | 0 | 0.585 (10) | fail |
| 50:50 | 3% | grouped | 0 | unknown | 0 | 10 | undefined (0) | fail |
| 50:50 | 3% | pooled | 4 | 0.9 | 10 | 0 | 0.658 (10) | fail |
| 50:250 | 15% | grouped | 0 | unknown | 0 | 10 | undefined (0) | fail |
| 50:250 | 15% | pooled | 0 | 0.8 | 10 | 0 | 0.544 (10) | fail |
| 50:250 | 3% | grouped | 0 | unknown | 1 | 9 | 0.599 (1) | fail |
| 50:250 | 3% | pooled | 0 | 0.8 | 10 | 0 | 0.547 (10) | fail |

The gate requires recovery ≥ 8/10, mean false activity ≤ 1 and convergence
≥ 9/10. Pooled passes the latter two in every cell but fails recovery in all
four. Grouped fails recovery and convergence in every cell, independently of
the unknown false-activity means. Consequently, the conservative no-error
condition in the implementation has no effect on any of these verdicts.

Every grouped error is the same explicit nonfinite loading-KL error. The two
returned grouped fits are seed 4 at 50:50 and 15% methylation, and seed 3 at
50:250 and 3% methylation. They converge in 24 and 35 sweeps, respectively.
Their matched absolute loading correlations are 0.118 and 0.135. Each has one
false-active program. All pooled attempts return converged fits. No completed
fit reaches the iteration cap without convergence. Warning counts are zero
in every row, so there is no study `warnings.csv` to populate.

### Block recovery and modality precision

`blocks.csv` reports signed block correlations for each returned fit. In the
two returned grouped fits, expression/methylation correlations are
0.126/0.109 at 50:50 and 0.049/−0.255 at 50:250. These are not recovered
shared programs. Their methylation precision shares in the matched program
are 0.262 and 0.650, below feature shares 0.5 and 0.833.

For pooled priors, the matched program's mean methylation precision share is
0.495 at 50:50 with 15% methylation, 0.518 at 50:50 with 3%, 0.839 at 50:250
with 15%, and 0.785 at 50:250 with 3%. Methylation exceeds its feature share
in 7/10, 4/10, 7/10 and 5/10 returned pooled fits, respectively. Dominance
is defined per fit, so it can occur even when the mean share is below the
feature share. `precision.csv` contains both modality shares for every
program, with live-program flags. No weights were introduced.

`priors.csv` gives each returned fit's point-mass weight and Laplace scale
for every modality and program, alongside its true nonzero fraction. The
38 failed grouped attempts have no final fitted prior to report. The two
returned grouped fits cannot support a general comparison of learned
sparsity. Both grouped and pooled fits include dead or nearly uninformative
programs, so averaging scale parameters across all programs is not an
estimate of the planted loading distribution.

### Interpretation and limits

The registered fixed implementation does not pass its simulation criteria.
The pooled comparator recovers some prognostic directions at 50:50 but none
at 50:250. The many grouped errors also prevent a meaningful paired
performance comparison. A mean C based on one returned grouped fit is not a
ten-seed performance estimate and should not be contrasted with the pooled
mean as evidence of improvement.

The grouped failure has a specific mathematical cause. The sequential
loading sweep holds its incoming prior fixed. Refitting the prior can remove
the continuous slab while its sweep posterior retains positive slab mass.
KL(q || g_new) is then infinite. An independent reproduction found incoming
slab probability 0.9999983 and rate 68.24035, followed by a pure point-mass
prior, while q retained positive second moments from 0.0003867 to 0.0010562.
The fitter's explicit failure preserves this incompatibility. Clamping the
prior, recomputing q or changing the update order would require a new design
decision and was not attempted.

The pooled path retains its historical incoming-prior KL accounting for
bitwise compatibility, whereas grouped checks the returned prior. Thus this
experiment evaluates the two fixed implementations, including that stopping
and failure-handling difference. It cannot attribute all differences to
modality-specific shrinkage alone, and their absolute monitored objectives
must not be compared as evidence favoring either prior.

The generator truncates molecular data at zero. D⁻¹F describes its underlying
linear reconstruction component rather than an exact decomposition of the
truncated matrix. The regenerated hazard nonetheless uses this exact g, so
recovery and survival have one registered target. Training/validation
censoring calibration follows the generator separately in each cohort.
Neither that calibration nor this low-signal truncated scenario reproduces
all aspects of the S-9 design.

## Step 4: matched TCGA comparison

**Not run.** Grouped fails the gate in every simulation cell, triggering the
registered stop before step 4. Pooled passes no cell either. Task C's prepared
matched-patient data do not override this rule. No matched-TCGA model, CV
comparison, joint-only program list or real-data adoption was produced.

## 2026-11-12 rerun: empty blocks and recovery diagnostics

The [approved registration](../prompts/2026-11-12-codex-d-rerun-prompt.md)
keeps the same four cells, seeds, survival draws and working configuration.
It changes the grouped empty-block handling and adds expression-only fits
and true-risk C. The [new aggregate directory](../../results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_rerun/)
and corresponding ignored `data/feature_group_simulation_2026-11-12_rerun/`
private directory preserve the first run's directories unchanged.

**All 120 fits converge, with no errors and no study warnings.** The empty-block
crash is resolved. Every arm and cell still fails the recovery threshold.
Grouped is the default selected prior because pooled passes no cell, but the
selected arm passes no cell either, so step 4 did not run. This rerun does not
adopt a multimodal model or a weighting scheme.

### Empty-block implementation and tests

`yf_update_F_groups()` still obtains a sequential block posterior under its
incoming prior. If the refitted prior is a pure point mass at zero, it now
subtracts that block's just-returned mean contribution from the returned
E[Z_k] and its just-returned variance contribution from the returned V[Z_k].
It then sets both loading moments to exactly zero before the next group
updates. The refitted point-mass prior is stored. Since zero second moment
makes q exactly delta zero as well, that block's KL is zero. Non-point-mass
refits still call the existing KL calculation and retain its explicit
nonfinite-KL error. This is the approved posterior reset, rather than
clamping a prior probability or hiding an incompatible KL.

The helper's inputs and sweep calculations are unchanged. Its return list
adds final EZ/VZ and named reset-event counts for the groups. The fitter
accumulates those counts in `fit$empty_blocks`, a modality-by-program integer
matrix, and returns it only for canonical grouped fits. Every point-mass
refit counts, including repeated resets of an already empty block. The total
is a count of update events, not distinct removed programs. Low-information
blocks handled by the pre-existing precision threshold do not increment this
new refit counter. No VF carrier exists on this kernel path. Zero first and
second moments also imply zero loading variance.

The first new real-update regression uses a weak but nonzero block whose
sweep returns positive mean and variance before the prior collapses. It
checks that the next group receives only the remaining block's projection,
that final EZ = YF, and that VZ equals the projection variance reconstructed
from the returned moments. It independently checks that only the live
block contributes loading KL. The second fixture fits a program present in
one modality and absent in the other. It requires exact zero moments, a
point-mass prior, a positive counter, nonzero companion loadings and a
monitored objective recomputed with independent live-block loading KL.

Tests were committed first in `041a3e5d`. The new grouped checks initially
failed, including the original incompatible-KL error, while all 14 existing
fitter checks passed. After implementation, all 16 fitter checks pass,
including the frozen pooled and stratified-centered NULL/one-group outputs.
The full suite passes **849 checks with zero failures**. Its only warning is
the existing bounded-β optimizer profile-likelihood fallback. The derivation
verification passes **88/88**. The implementation and twelve-fit smoke
aggregates are in `a0bdb81d`. All smoke fits converged without errors or
warnings. Production changes are confined to the grouped helper, conditional
grouped metadata in the fitter and the simulation runner. The historical
multimodal fitter and single-group updates are unchanged.

### Expression-only and true-risk diagnostics

`fgsim_arm_inputs()` takes the prepared joint data and K, computes the same
signed joint SVD used in the first run, and returns the grouped, pooled and
expression-only input lists. The two joint lists are identical. Expression
uses the expression columns of training and validation data, the matching
rows of that joint F start and the same initial L. Its recovery truth is the
expression part of g renormalized to unit norm. Its survival outcomes and
true η remain those drawn from the joint model. It does not regenerate
expression-only outcomes or use a new SVD.

`fgsim_assert_expression_inputs()` checks the feature subsets, normalized
truth, training centers, group labels, both starting matrices and shared
survival/risk fields. The existing joint-arm identity assertion runs before
and after fitting, and the new expression assertion does likewise. Regression
checks mutate each field in turn and require failure.

`fgsim_true_risk()` evaluates `frozen_reverse_cindex()` on the validation η
used to generate outcomes, with no outcome-based sign change. This value is
computed before any fit, retained even if a fit fails, and repeated for all
three arms in the seed/cell. It is an oracle discrimination reference, not a
strict upper bound on an empirical C in a finite validation sample.

`fgsim_one()` fits the three arms sequentially, using each arm's own feature
matrix and recovery target. Its metrics add signed expression-block loading
correlation for the matched program, empty-block events and true-risk C.
Joint programs are matched on the full g before reporting the expression
block correlation. Expression-only programs are matched on its normalized
expression part. Its precision and prior tables use only the modalities
actually fitted. The `p_expression` and `p_methylation` metadata identify the
source simulation cell, so expression-only rows still identify that cell's
methylation count even though only 50 expression features are fitted.
`empty_blocks.csv` gives every grouped fit's counts by modality and program.
Non-grouped arms report zero events in metrics because they do not execute
the grouped reset. These zeros do not assert that their priors stay nonempty.

`fgsim_gate()` now includes total reset events with their defined-fit count
and mean true-risk C with its defined count. The numerical acceptance
thresholds are unchanged. `fgsim_select()` follows the rerun's explicit
selected-arm rule: grouped is preferred unless pooled passes a cell that
grouped fails, and step 4 requires the selected arm to pass at least one cell.
This supersedes the earlier additional all-grouped-fail stop if pooled passes
some cell. The distinction has no effect on this rerun because neither joint
arm passes any cell. The runner uses one CPU worker and one BLAS thread.
Only aggregate CSVs/settings/hashes are committed. Fits stay under ignored
`data/`. The full rerun aggregates are committed in `d526d066`.

### Rerun gates and discrimination

| Features, expression:methylation | Methylation nonzero | Arm | Recovery /10 | Mean false active | Converged /10 | Errors | Empty-block events | Validation C (defined n) | True-risk C (defined n) | Gate |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---|
| 50:50 | 15% | grouped | 1 | 1.2 | 10 | 0 | 258 | 0.577 (10) | 0.773 (10) | fail |
| 50:50 | 15% | pooled | 2 | 0.8 | 10 | 0 | 0 | 0.585 (10) | 0.773 (10) | fail |
| 50:50 | 15% | expression only | 1 | 0.6 | 10 | 0 | 0 | 0.539 (10) | 0.773 (10) | fail |
| 50:50 | 3% | grouped | 5 | 0.9 | 10 | 0 | 572 | 0.685 (10) | 0.790 (10) | fail |
| 50:50 | 3% | pooled | 4 | 0.9 | 10 | 0 | 0 | 0.658 (10) | 0.790 (10) | fail |
| 50:50 | 3% | expression only | 5 | 0.3 | 10 | 0 | 0 | 0.617 (10) | 0.790 (10) | fail |
| 50:250 | 15% | grouped | 0 | 0.9 | 10 | 0 | 427 | 0.533 (10) | 0.787 (10) | fail |
| 50:250 | 15% | pooled | 0 | 0.8 | 10 | 0 | 0 | 0.544 (10) | 0.787 (10) | fail |
| 50:250 | 15% | expression only | 0 | 0.2 | 10 | 0 | 0 | 0.513 (9) | 0.787 (10) | fail |
| 50:250 | 3% | grouped | 0 | 0.9 | 10 | 0 | 153 | 0.559 (10) | 0.781 (10) | fail |
| 50:250 | 3% | pooled | 0 | 0.8 | 10 | 0 | 0 | 0.547 (10) | 0.781 (10) | fail |
| 50:250 | 3% | expression only | 2 | 0.3 | 10 | 0 | 0 | 0.576 (8) | 0.781 (10) | fail |

Every cell has ten attempted fits per arm. The gate requires recovery ≥ 8/10,
mean false activity ≤ 1 and convergence ≥ 9/10. All arms pass convergence.
Grouped at 50:50 and 15% methylation also fails false activity, with mean 1.2.
Every arm fails recovery. The conservative no-error condition carried from
the original runner has no effect here because no fit errors.

The three undefined validation C values are expression-only seed 3 at
50:250 with 15% methylation, and seeds 7 and 10 at 50:250 with 3%. Their
returned risk scores are constant. They remain in all gate denominators and
are excluded only from the descriptive C means, whose defined counts are
shown. True-risk C is defined for every seed/cell and identical across arms.
There are no unconverged fits or study warning messages.

All 40 pooled fits reproduce the first run exactly in recovery, false
activity, matched program/correlation, convergence, iteration count,
validation C and monitored objective. This verifies that the grouped fix
and new diagnostics did not change the pooled comparator.

### Expression recovery and modality precision

Mean absolute expression-block loading correlation for the matched program
is reported below, with the defined count in parentheses. A zero or constant
expression loading block has undefined correlation. These values are excluded
only from this descriptive mean. Matching remains arm-specific as described
above. Signed values for every seed are in `metrics.csv` and `blocks.csv`.

| Features, expression:methylation | Methylation nonzero | Grouped | Pooled | Expression only |
|---|---:|---:|---:|---:|
| 50:50 | 15% | 0.231 (10) | 0.237 (10) | 0.241 (10) |
| 50:50 | 3% | 0.553 (9) | 0.541 (10) | 0.521 (10) |
| 50:250 | 15% | 0.214 (9) | 0.168 (10) | 0.179 (9) |
| 50:250 | 3% | 0.259 (9) | 0.258 (10) | 0.377 (8) |

The matched program's methylation L-precision share and the number of seeds
where it exceeds methylation's feature share are:

| Features, expression:methylation | Methylation nonzero | Grouped mean share | Grouped dominance /10 | Pooled mean share | Pooled dominance /10 |
|---|---:|---:|---:|---:|---:|
| 50:50 | 15% | 0.500 | 7 | 0.495 | 7 |
| 50:50 | 3% | 0.508 | 5 | 0.518 | 4 |
| 50:250 | 15% | 0.545 | 3 | 0.839 | 7 |
| 50:250 | 3% | 0.667 | 2 | 0.785 | 5 |

All matched-program shares above are defined. The full table contains both
modality shares for every live program. Feature-share thresholds are 0.5
and 5/6. Expression-only has precision share 1 for its sole modality when
defined. The grouped reset totals are 258, 572, 427 and 153 across the four
cells in this table, summing to **1,410 events**. They split by expression and
methylation as 105/153, 305/267, 336/91 and 139/14, respectively. Per-program
counts are retained in `empty_blocks.csv`. The fitted-prior table contains
all returned group/program priors, point-mass weights, Laplace scales and
true nonzero fractions, as in the first run.

### Registered reading of problem 2 and step 4

Expression-only fails the 8/10 recovery criterion in every cell, with
1/10 and 5/10 at 50:50 and 0/10 and 2/10 at 50:250. There is no cell where
expression-only passes recovery while both joint arms fail. Under the
registered descriptive rule, this scenario is too weak at the current sample
size for the gate to separate the arms. It does not establish that adding
methylation causes the recovery failure. The proposed next decision is to
choose a scenario in which the expression-only arm can pass, before using
the gate to assess an added-modality benefit. No scenario or weights were
changed in this task.

The true-risk C means of 0.773–0.790 show that the known joint direction
carries prognostic information. They do not establish that its loadings are
identifiable from 120 training patients or that this optimizer can recover
them. Expression-only also omits the true methylation contribution to the
hazard and uses matching rows of the joint SVD start, as registered. Its
failure cannot distinguish weak molecular recovery signal from fitting or
initialization limitations. Grouped recovery improving from 4/10 pooled to
5/10 in one cell is below the gate and does not establish a benefit.

**Step 4 was not run.** The registered rule selects grouped by default, but
it passes no cell. There are no matched-TCGA rerun fits, CV results or
joint-only program claims. The returned-prior KL in grouped fits and the
historical incoming-prior KL in pooled fits still differ. The reset makes
q and g compatible in an empty grouped block, but does not remove that
accounting difference for live blocks. Absolute objectives remain unsuitable
for a grouped-versus-pooled model comparison.

### Rerun independent review

A read-only statistical/code review found no blocking defect before the full
run. It verified returned-sweep projection removal, zero empty-block KL,
conditional counts, preservation of the single-group path, expression inputs,
true-risk scoring and the selected-arm rule. Its requests to distinguish
true-risk C from a strict finite-sample bound and to qualify the descriptive
weak-scenario interpretation are addressed above.

A fresh final reviewer independently audited the code and all 120 aggregate
attempts before reading these findings. The review verified all gate and
diagnostic tables, reconciled the 400 modality/program event rows with the
1,410 total resets, confirmed identical true-risk C across each seed/cell
triplet, checked all 12 source hashes and confirmed exact reproduction of
all 40 original pooled results. The suite and derivation logs confirmed
849/0 and 88/88. The subsequent findings review requested defined counts
for expression-block correlation means. Those counts and the handling of
zero/constant blocks are now explicit. The final verdict was **PASS for the
approved rerun and stop decision, with no blocking code, statistical or
reporting defects**. No scientific adoption is supported.

## 2026-11-12 informative-scenario study

The [next-steps registration](../prompts/2026-11-12-codex-d-next-prompt.md)
changes the generator scenario to `both_informative` while retaining 120
training patients, 120 validation patients, the four feature-count/sparsity
cells, seeds 1–10 and the three paired arms. Its outcome does not change the
[matched-TCGA comparison](2026-11-12-multimodal-tcga-comparison-findings.md),
which remains **not validated by simulation**.
Earlier study directories remain unchanged. The
[new aggregate directory](../../results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_informative/)
contains the 120-fit results, including descriptive program-2 recovery and
per-cell separation verdicts. Synthetic fits remain under ignored `data/`.

### Two-direction simulation and tests

`fgsim_prepare()` now accepts the registered target programs. Its default
retains the earlier program-3 calculation. For `programs = 1:2`, it forms each
g_k by stacking D⁻¹F[,k] across modalities in fitted order and normalizing
that column to unit norm. It constructs
η = a X_c(g_1 β_1 + g_2 β_2), using the generator's β_1 = 0.70 and β_2 = 0.35.
The multiplier a preserves the original generator's training η SD. Validation
uses the same directions, multiplier, training centers and feature SDs.
Event and censoring draws retain the generator settings and seed + 12000.
Program 3 is not prognostic. The generator itself and the fitter are unchanged.

`fgsim_arm_inputs()` and `fgsim_assert_expression_inputs()` handle two-column
g in this mode. Each expression subcolumn is renormalized separately. All
arms receive the same outcomes and true joint risk. Expression-only still
uses the matching rows of the common joint SVD start and the same initial L,
rather than its own expression SVD. These choices preserve the registered
paired-start diagnostic.

`fgsim_recovery_two()` applies the existing maximum absolute loading-correlation
match independently to programs 1 and 2. Recovery requires |cor| ≥ 0.7,
posterior survival activity and E[β_k] cor(g,F_k) > 0. Program 1 alone controls
the gate. Program 2 is descriptive. During execution, Andrew confirmed the
S-9 false-activity convention: survival-active fitted columns are exempt when
they have a threshold-qualified loading match to either truth, even if the
coefficient sign is wrong. A wrong sign still fails recovery. Exempt indices
are unique, so choosing the same fitted column for both directions cannot
subtract it twice. The return includes both matches/correlations, recovery
flags, per-modality correlations and a shared-best-column indicator.

The tests reconstruct both unit directions and the exact Weibull/censoring
draws, check the risk-SD calibration and training-only validation transform,
and require paired starts/outcomes and normalized expression truth. Scoring
checks distinguish loading matches from correct-sign recovery, qualify both
true survival programs for false-activity accounting, retain unmatched active
programs, reject subthreshold exemptions and avoid double exemption when
best columns coincide. Tests were committed first in `ba50fbb4`, with the
confirmed false-activity checks in `dab028bf`.

`fgsim_separation()` reads the existing gate table. A cell is eligible only
when expression-only passes the complete gate. In an eligible cell, grouped
passing while pooled fails is labelled helped, pooled passing while grouped
fails is labelled hurt, and matching gate outcomes are reported directly.
Ineligible cells are labelled not testable, regardless of their joint-arm
results. Regression checks cover an eligible benefit and a pooled advantage
that remains uninterpretable until the expression gate passes.

The new `--informative` command writes new directories, keeping the historical
low-variance mode available. Source and configuration hashes are saved.
Program-2 counts, both matching fields, true-risk C, reset counts and precision
for both candidate matches are retained. The prior table's true nonzero
fraction is a modality-level reference: all three true programs have the
same fixed nonzero count within a modality. It is not a matched support
precision metric. Twelve smoke fits completed before the full run. The
runners and smoke evidence are committed in `c2bde939`. The final full suite
passes **856 tests with zero failures**, with the existing bounded-β fallback
warning. No fitter or update kernel changed in this continuation.

### Informative-scenario results

| Features, expression:methylation | Methylation nonzero | Arm | Program 1 recovery /10 | Program 2 recovery /10 | Mean false active | Converged /10 | Empty events | Validation C (n) | True-risk C (n) | Gate |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---|
| 50:50 | 15% | grouped | 9 | 5 | 0.9 | 10 | 285 | 0.785 (10) | 0.817 (10) | pass |
| 50:50 | 15% | pooled | 9 | 5 | 0.9 | 10 | 0 | 0.783 (10) | 0.817 (10) | pass |
| 50:50 | 15% | expression only | 5 | 2 | 0.9 | 10 | 0 | 0.708 (10) | 0.817 (10) | fail |
| 50:50 | 3% | grouped | 10 | 6 | 0.4 | 10 | 673 | 0.797 (10) | 0.827 (10) | pass |
| 50:50 | 3% | pooled | 9 | 7 | 0.1 | 10 | 0 | 0.789 (10) | 0.827 (10) | pass |
| 50:50 | 3% | expression only | 9 | 4 | 0.2 | 10 | 0 | 0.746 (10) | 0.827 (10) | pass |
| 50:250 | 15% | grouped | 4 | 0 | 1.6 | 10 | 309 | 0.733 (10) | 0.806 (10) | fail |
| 50:250 | 15% | pooled | 4 | 0 | 1.7 | 10 | 0 | 0.716 (10) | 0.806 (10) | fail |
| 50:250 | 15% | expression only | 7 | 2 | 0.4 | 10 | 0 | 0.706 (9) | 0.806 (10) | fail |
| 50:250 | 3% | grouped | 5 | 3 | 0.9 | 10 | 465 | 0.706 (10) | 0.809 (10) | fail |
| 50:250 | 3% | pooled | 4 | 3 | 1.0 | 10 | 0 | 0.676 (10) | 0.809 (10) | fail |
| 50:250 | 3% | expression only | 5 | 4 | 0.3 | 10 | 0 | 0.693 (10) | 0.809 (10) | fail |

All 120 fits converge with **zero errors and zero study warnings**. All ten
seeds remain in every gate denominator. The expression-only validation C at
50:250, 15% methylation, seed 3 is undefined because its score is constant.
It is excluded only from the descriptive validation C mean, whose defined
count is shown. True-risk C is defined in every seed/cell and identical
across its three arms.

Grouped fits record **1,732 empty-block reset events**, including repeated
resets of an already empty modality/program block. The totals are 285, 673,
309 and 465 in the displayed cell order. Event counts by modality/program
are in `empty_blocks.csv`. They are not counts of unique programs removed.

Independent best-column matching selects the same fitted column for both
truth directions in 33 fits: 21 expression-only, five grouped and seven
pooled. None has |cor| ≥ 0.7 for both directions. Thus these are coinciding
best candidates, not successful recovery of two distinct programs by one
column. Program-2 recovery is lower than program 1 in most cells and is
reported without changing the primary gate.

For the program-1 matched column, methylation's mean L-precision share and
its dominance count relative to its feature share are:

| Features, expression:methylation | Methylation nonzero | Grouped mean share (n) | Grouped dominance /10 | Pooled mean share (n) | Pooled dominance /10 |
|---|---:|---:|---:|---:|---:|
| 50:50 | 15% | 0.470 (10) | 5 | 0.482 (10) | 5 |
| 50:50 | 3% | 0.424 (10) | 4 | 0.514 (10) | 6 |
| 50:250 | 15% | 0.524 (10) | 1 | 0.669 (10) | 1 |
| 50:250 | 3% | 0.600 (10) | 1 | 0.680 (10) | 1 |

These are defined in all joint fits. The feature-share thresholds are 0.5
and 5/6. `precision.csv` retains all live-program shares and both truth-match
flags. Per-modality loading correlations for each matched truth are in
`blocks.csv`, and the final fitted priors are in `priors.csv`. No modality
weights were introduced.

### Registered separation verdict

| Features, expression:methylation | Methylation nonzero | Expression gate | Grouped gate | Pooled gate | Registered comparison |
|---|---:|---|---|---|---|
| 50:50 | 15% | fail | pass | pass | not testable |
| 50:50 | 3% | pass | pass | pass | both priors pass |
| 50:250 | 15% | fail | fail | fail | not testable |
| 50:250 | 3% | fail | fail | fail | not testable |

The informative scenario passes the expression-only gate in one cell,
50:50 with 3% methylation nonzero loadings. It can therefore evaluate the
priors in that cell under the registered rule. Grouped recovers 10/10 and
pooled 9/10, with false activity 0.4 and 0.1, respectively. Both pass. The
rule supports neither a helped nor a hurt verdict. Their mean validation C
values, 0.797 and 0.789, are descriptive and do not replace that rule.

Expression-only fails the gate in the other three cells. Joint recovery
at 50:50 with 15% methylation is reported, but that cell is ineligible for
the registered grouped-versus-pooled comparison. Both joint arms fail at
50:250. This study does not establish a general advantage from separate
priors, nor does it establish that the sample size cannot test the priors
in every cell: one cell is eligible. It does not alter Part A's descriptive
status or adopt a new model.

The true-risk C is an oracle reference rather than a strict finite-sample
upper bound. The scenario still truncates molecular measurements at zero.
The standardized directions define the regenerated hazard exactly, but
represent the underlying reconstruction component rather than an exact
linear decomposition of the truncated observations. Expression-only omits
the true methylation contribution to risk and uses the joint SVD start.
These design limits remain relevant when interpreting an expression-gate
failure. Grouped uses returned-prior loading KL and pooled retains historical
incoming-prior accounting, so absolute objectives must not be compared as
though only the prior grouping changed.

### Informative-scenario independent review

A fresh reviewer read the plan, code and tests before the findings
conclusions. It independently recomputed recovery/sign checks, both best
matches, threshold-qualified false activity, candidate collisions, modality
precision, validation/true-risk C, convergence and iteration counts from all
120 private checkpoints. It verified grouped truth normalization, training
risk-SD calibration and validation centers, rebuilt all 12 gate rows, and
reconciled all 600 program rows and 400 grouped event rows. The review found
no blocking code or aggregate defect. The aggregate evidence is committed in
`413481e0`. The subsequent prose review verified all displayed gate,
precision, collision and separation values and returned **PASS**, with no
blocking code, statistical or reporting defect. The TCGA findings receive
a separate review after their full run.

## Independent review

A fresh reviewer inspected the tests and existing update code without editing
files or running fits. The initial review identified that the proposed objective
test used the incoming prior's KL, potentially missing disagreement with the
returned refitted prior. It also identified positional argument capture,
potentially vacuous rescaling assertions and an incorrect assumption about an
`ebnm()` output field. The tests were revised. Follow-up review found the
Laplace KL correction and canonical rate transformation mathematically
consistent and approved the tests as a pre-implementation checkpoint. The
coverage limits above remain explicit. That checkpoint preceded the final
review of the fitter, scientific results and completed findings.

A separate independent statistical/code reviewer found no blocking fitter
defect in the working configuration. The review verified common EZ/VZ
propagation, group reconstruction inputs, the returned-prior KL correction
and prior scaling. Its requests for a multi-sweep objective test, analytic
Normal/exponential KL tests and explicit disclosure of the single-group
objective difference were addressed. Review of the proposed simulation then
identified the reconstruction/prognostic target ambiguity described above.
That earlier pause was resolved by the finish prompt. The completed study
uses the shared standardized target described above.

### Continuation review

The smoke-check reviewer independently reproduced the grouped support
incompatibility and found no blocking runner defect. The review confirmed
that completing the registered attempts with explicit failures was valid
without changing the design. It requested explicit disclosure that the
implementation conservatively fails cells with fit errors. That condition
and its lack of effect on the observed verdicts are documented above.

A fresh final reviewer first read the requirements and code without the
findings conclusions, then independently audited all 80 unique attempts and
all 42 saved fits. The reviewer recomputed recovery and validation C,
regenerated truth and survival from the private inputs, checked exact paired
inputs in the two returned grouped pairs, verified all 12 source hashes and
all 420 prior-table rows, and reran the 11 focused regression checks.
The subsequent findings review checked every gate and diagnostic table.
The final verdict was **PASS for the registered diagnostic execution and
stop decision, with no blocking code, statistical or reporting errors**.
No scientific adoption was supported. Two minor wording requests were
addressed: the title now describes the completed simulation, and the earlier
review is explicitly identified as a historical checkpoint.

## Text for the shared status files

Task D's informative-scenario study is complete after the empty-block fix.
All 120 fixed fits converge without errors or study warnings. The two-direction
hazard uses unit standardized loading directions for programs 1 and 2,
β = (0.70, 0.35), and preserves the original training risk SD. Program 1
controls recovery. Program 2 is descriptive. Andrew confirmed that
threshold-qualified loading matches to either true program are exempt from
false activity even with a wrong-sign coefficient, while recovery still fails.

The expression-only gate passes in one cell, 50:50 features with 3%
methylation nonzero loadings. Grouped recovers 10/10 and pooled 9/10 there,
with mean false activity 0.4 and 0.1. Both pass, so the registered rule supports
neither a benefit nor a disadvantage from separate priors. The other three
cells are ineligible for the prior comparison because expression-only fails
the gate. Both joint arms pass at 50:50 with 15% methylation, but neither
passes at 50:250. There is no general separate-prior adoption. Part A's
matched-TCGA analysis remains descriptive and **not validated by simulation**,
regardless of these results. The full suite passes 856 tests with zero failures.
