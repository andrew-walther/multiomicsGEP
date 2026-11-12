# Feature-group YFB fitter and registered simulation

Author: Andrew Walther

Date: 2026-11-12

## Current result

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

Task D's registered grouped-versus-pooled comparison is complete at K = 5,
with 50:50 and 50:250 expression:methylation features and equal or unequal
modality sparsity. Survival was regenerated from the unit standardized
reconstruction direction, preserving the generator's training risk SD.
Every cell and arm fails the acceptance gate. Grouped has zero recoveries,
38 explicit refitted-prior support errors and two converged fits. Pooled
converges in all 40 attempts, with recovery 2/10 and 4/10 at 50:50 and 0/10
in both 50:250 cells. The rule retains grouped as its default preference but
authorizes neither prior for step 4. The matched-TCGA comparison did not run.
There were no study warnings. The full suite passes 844 checks with zero
failures and its existing bounded-β fallback warning. No modality weights,
prior-support fix or model adoption was introduced. The next decision is
whether to register a prior-update strategy that preserves compatible
support, then reevaluate the grouped implementation. These results do not
establish a benefit from separate modality priors.
