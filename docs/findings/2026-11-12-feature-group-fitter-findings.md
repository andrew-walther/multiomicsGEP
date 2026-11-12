# Feature-group YFB fitter and simulation preparation

Author: Andrew Walther

Date: 2026-11-12

## Current result

The feature-group fitter is implemented. NULL and one-group fits preserve EL,
EF, EBeta and the ELBO history bit for bit on pooled and stratified-centered
fixtures. Multiple groups fit separate loading priors with a shared L and joint
Cox projection. The full working-checkout suite passes **836 tests with zero
failures**, including 14 fitter checks and three simulation-helper checks.
All **88 derivation checks pass**. The full suite emits the existing bounded-β
optimizer profile-likelihood fallback warning.

Tests were written first and committed as `8a3b5517`. Fitter edits began only
after Part 1 finished and pushed its aggregate results and findings through
`cbe90e44`. The simulation log reached 104/104 and no Part 1 fitting or scoring
process remained. The repository stayed on main without a worktree.

The focused test command is:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript tests/test_feature_groups_yfb.R
```

The focused fitter command now passes 14/14. Before implementation, three
prerequisite checks passed and seven feature-group checks failed as expected.
Master-suite registration was deferred during Part 1 and enabled afterward.

The simulations have not run. Independent review identified a recovery-target
ambiguity created by the fixed z-scored input. Andrew has been asked which
target should define recovery and false activity. This is a software result,
not acceptance of the multimodal model.

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

## Simulation and matched-data evaluation

Neither simulation ratio has run. Recovery ≥ 8/10, mean false activity ≤ 1 and
convergence ≥ 9/10 remain unassessed. Modality precision shares remain unassessed.
The five-fold matched-data comparison has not run because the fitter and its
simulation acceptance gate are pending. No modality weights or model choices
have changed.

The prepared sequential runner is
`results/benchmark_sim/run_feature_group_simulation.R`. It uses generator
defaults for `low_variance_prognostic_factor`, seeds 1–10, K = 5 from the
working-rank config, and 50 expression features against 50 or 250 methylation
features. Training means and SDs scale validation without refitting. Helpers
return per-cell recovery, false activity, convergence and precision shares
Σ_j τ_mj E[F²_mjk]. Undefined dead-program precision shares remain NA.
Warnings are retained individually. Failed cells stay in the requested
denominator and cannot pass the gate. Synthetic fit checkpoints stay under
ignored data/, with aggregate-only CSVs in the output directory. Existing
study directories are refused. BLAS thread limits must be set in the launching
shell before numerical libraries load. Direct execution is currently blocked
with an explicit recovery-target decision message.

The provisional runner compares reconstruction truth F/s to the fitted
loadings, following the earlier raw-input comparator's correlation threshold
of 0.7 plus survival activity. This target is not yet approved. If
X = (Y − μ)D⁻¹, the underlying reconstruction loading is D⁻¹F, while the generator's
survival score satisfies YFβ = X(DF)β + μFβ. The Cox constant does not affect
ranking, so DF is the standardized prognostic direction. The generator truncates
observations at zero by default. D⁻¹F therefore describes its underlying linear
reconstruction component, not an exact decomposition of the truncated data.
The survival projection identity remains exact. The earlier raw-input
comparison did not distinguish these targets. Reporting raw-F correlations
alongside F/s does not resolve this issue. The recovery and false-active
definition must be settled before executing the study. An unequal-SD
regression test will verify the selected target and projection identity.

## Independent review

A fresh reviewer inspected the tests and existing update code without editing
files or running fits. The initial review identified that the proposed objective
test used the incoming prior's KL, potentially missing disagreement with the
returned refitted prior. It also identified positional argument capture,
potentially vacuous rescaling assertions and an incorrect assumption about an
`ebnm()` output field. The tests were revised. Follow-up review found the
Laplace KL correction and canonical rate transformation mathematically
consistent and approved the tests as a pre-implementation checkpoint. The
coverage limits above remain explicit. This is not the required final review
of the fitter, scientific results and completed findings.

A separate independent statistical/code reviewer found no blocking fitter
defect in the working configuration. The review verified common EZ/VZ
propagation, group reconstruction inputs, the returned-prior KL correction
and prior scaling. Its requests for a multi-sweep objective test, analytic
Normal/exponential KL tests and explicit disclosure of the single-group
objective difference were addressed. Review of the proposed simulation then
identified the reconstruction/prognostic target ambiguity described above.
The simulation remains paused for Andrew's decision. Final scientific review
cannot occur before the study runs.

## Text for the shared status files

Task D's grouped fitter is implemented after Part 1's completion and output
commit. NULL and one-group calls reproduce the frozen pooled and
stratified-centered moments and ELBO histories exactly. The full suite passes
836 checks and the derivation passes 88 checks. The grouped path corrects
loading KL to its refitted group priors, while the single-group path retains
historical accounting for compatibility. Simulation and matched-data
evaluation have not started. Simulation recovery needs a decision between
the standardized reconstruction loading D⁻¹F and prognostic direction DF.
