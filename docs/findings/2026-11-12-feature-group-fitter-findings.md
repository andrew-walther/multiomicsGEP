# Feature-group YFB fitter: tests before implementation

Author: Andrew Walther

Date: 2026-11-12

## Current result

Implementation is waiting for Part 1 of the D/E follow-up prompt to finish and
commit its outputs. This session has not edited the fitter or any file it
sources. At the last check, the prescribed Part 1 simulation and training logs
were absent from `data/`, and no corresponding Rscript process was visible.
Absence of a process does not establish completion. Part 1 has uncommitted code
changes in the shared checkout, including its starting-state instrumentation in
`code/fit_ssbmf_yfb.R`. Those changes belong to Part 1.

The focused test command is:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript tests/test_feature_groups_yfb.R
```

It exits with status 1: **3 checks pass and 7 fail**. Both unchanged-fitter
baselines and the independent KL calculation pass. The feature-group tests fail
because the new argument is not implemented. The malformed-input test rejects
the current generic unused-argument error because it requires explicit validation.
This is a test checkpoint, not a completed fitter or scientific result.

The new file is deliberately not yet registered in `tests/run_tests.R` so that
its expected failures do not interrupt Part 1's required regression run. Add it
when implementing the fitter. The full suite has not been run for this checkpoint.

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
rescaling of every live group's prior. The proposed result contract is
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
incoming prior. These can differ. Implementation must resolve that accounting
without changing the required single-group baseline.

**Inputs and outputs.** Fixtures return a patient-by-feature matrix, survival
vectors and fitter arguments. The observer returns the fit and captured update
records. The KL helper returns one scalar per block. The standalone test file
prints named assertions and returns a nonzero process status for any failure.

**Limits.** The objective check independently verifies F accounting but uses
the existing canonical objective helper for the likelihood terms and captured
production KL values for L and β. Grouped restart, pruning and centered
prediction need checks or explicit rejection once supported options are
settled. Bitwise frozen fixtures depend on R, BLAS and package versions. Do not
regenerate them with a fitter that already implements feature groups.

**R details.** `ebnm()` supplies the independent block-prior reference.
Second moments are reconstructed as posterior mean squared plus posterior SD
squared, because its default posterior output has no `second_moment` field.
Base R `integrate()` provides the independent loading-prior density calculation.
Observer bindings use `match.call()` so named helper arguments are supported.

## Simulation and matched-data evaluation

Neither simulation ratio has run. Recovery ≥ 8/10, mean false activity ≤ 1 and
convergence ≥ 9/10 remain unassessed. Modality precision shares remain unassessed.
The five-fold matched-data comparison has not run because the fitter and its
simulation acceptance gate are pending. No modality weights or model choices
have changed.

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

## Text for the shared status files

Task D has reviewed tests and frozen pooled and stratified-centered baselines.
The standalone file passes three prerequisite checks and has seven expected
failures before feature groups are implemented. Fitter edits wait for Part 1's
completed runs and committed outputs. Simulation and matched-data evaluation
have not started. The grouped objective must account for the returned F priors
rather than silently reuse KL values computed under a different prior.
