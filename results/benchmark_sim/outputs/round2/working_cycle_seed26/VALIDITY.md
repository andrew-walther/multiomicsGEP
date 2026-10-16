# Invalid objective diagnostic — preserve original outputs

The one Task2-6 replay completed after preregistration88df7d3 in42.345 seconds,
with no warnings, exact entire-fit/history identity and complete late-cycle
density coverage. Native/ledger moment disagreement was at most2.22e-16.
The original fit remains unconverged at300; no scientific criterion changed.

Post-run review found an observer-only reconstruction error: its expected
squared residual omitted cross-products between different programs. The
original one-program objective fixture could not detect this. A new independent
two-program fixture failed by.0064, then passed after the observer was corrected
to use (Y−EL EFᵀ)² + EL2 EF2ᵀ − EL²(EF²)ᵀ.

The original telemetry, iteration/transition summaries and hashes are retained
unchanged. They are NOT valid complete-objective or ascent evidence. Do not
use their posterior-step decrease counts to diagnose the fitter. Exact fit
identity is still valid evidence that observation did not change fitting.

The corrected observer passes8 focused tests. No corrected scientific replay
has run at this checkpoint. Task2-6's no-retry rule applies: a further replay
needs Andrew's decision and its own prospective committed entry. This is an
observer defect, not a demonstrated new defect in either model update.
