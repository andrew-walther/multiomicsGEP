# Normalized-projection supervision evidence

These outputs support the October 30 execution report and progress-book chapter.
They do not replace the immutable historical benchmark.

- `precision_audit.csv`: conditional training-state precision audit, not a refit.
- `sign_variance_audit.csv`: aggregate coefficient variance check.
- `attempt1`–`attempt6`: retained numerical/design diagnostics. Their raw-risk
  simulation direction differs from the prepared reconstruction loading.
- `attempt7`: corrected shared-loading simulation; controlling scientific gate.
- `attempt7/simulation_units.csv`: raw/prepared signal variance shares, training
  predictor SD and synthetic oracle discrimination. Normalized coefficient
  units differ from the older simulation, so across-attempt C is not comparable.
- Each attempt's raw RDS checkpoints retain fits, original metrics, warnings
  and failures. They are not overwritten by the summarizer.
- `metrics.csv` adds derived diagnostics and preserves raw metrics alongside
  conservative failed-fit penalties. Missing records stay in the denominator.
- `summary.csv`, `gate_cells.csv`, `gate.csv` and `gate_criteria.csv` contain the
  complete-denominator decision. `raw_gate_diagnostic.csv` is descriptive only.
- In `counts.csv`, `raw_checkpoint_errors` excludes missing files; `raw_errors`
  and `errors` include missing records. `errors` also includes invalid stored
  coefficient variance. `nonconverged` includes all conservative failures.
  Error, variance-loss and cap counts can overlap. Undefined synthetic C is
  counted even when the available-record mean is finite. `mean_C` retains
  the all-finite-state descriptive mean; `mean_valid_C` and `valid_C_n`
  restrict prediction summaries to converged fits with valid uncertainty.
- `all_zero_returned_states` and `nonfinite_returned_eta` concern returned
  states only. `unavailable_factor_states` counts absent non-lasso states.
  Inactive factors and survival-only loading columns are separate descriptors.
- Attempts 1 and 5 are intentionally incomplete. Never fill their missing
  arms using a later engine. Attempt 5's stop reason is saved in its directory.
- `attempt7_diagnostics`: deterministic fixed-seed replays of saved registered
  fits, not additional gate observations. Loading means are compared directly.
- Root `diagnostic_*` files are historical attempt-4 replays. Do not confuse
  these with the final corrected-design replays.
- `provenance.csv`, `checkpoint_hashes.csv`, package/R versions, test summaries
  and compressed logs identify implementation and verification evidence.

All simulation checkpoints and failing pseudo-messages are synthetic.
`training_only.rds` contains patient-level training inputs, is git-ignored,
and MUST NOT be committed or copied into presentation resources.

Reproduce only using the registered engine and simulation revision. The
attempts log records each registration before execution. Do not pool repeated
attempts, choose a weight from synthetic validation C, or change the gate.
