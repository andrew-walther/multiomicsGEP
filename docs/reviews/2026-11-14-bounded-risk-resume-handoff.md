# Resume handoff: fixed-prior convergence screen

Date: 2026-11-14. Author: Andrew Walther.

## Paused state

Andrew requested a usage-limit pause. Both Step 1 simulation workers were
interrupted deliberately and are stopped. Their exit status 1 is an interruption,
not a new scientific failure. No worker or delegated task remains active.

Branch: `s17-bounded-risk-supervision`, from main `6d12a0dd`.
Implementation: `e0c5697b`, reviewed, tested and pushed. The pause commit follows it.
The numerical source is unchanged. Do not edit it while resuming cached fits.

The real-data screen is complete: 36 actual fits, all eight pilots stable, no
Normal or Laplace scale passes all required targets. Normal scales 0.5 and 1
have individually stable target-0.05 arms, but fail stronger targets. Steps 2
and 3 must not run under this registration. No external scoring or merge.

Saved state is recorded in
`results/benchmark_sim/outputs/bounded_risk_2026-11-14/pause_inventory.csv` and
`pause_saved_fits.csv`. It contains 60 saved fit records:

- All eight real cases complete, 36 fits.
- Seven S-9 cases complete. The Laplace scale-1 case retains completed weights
  40 and 100. Weight 1,000 was interrupted. Total saved S-9 fits: 23.
- The weak Normal scale-0.1 case retains its completed pilot. Its target arm
  was interrupted. No other weak case is started.

There are 15 complete cases, two partial cases and seven absent cases out of
24. Every present case has matching numerical source/input hashes. Every saved
fit has finite τ, E[L²], E[F²], E[β²] and final Cox weights, with the expected
trace dimensions. The full saved-state audit remains to be run after completion.
The interrupted in-flight fits were not saved and will be rerun. Do not count
those interruptions as scientific nonconvergence or task errors.

Private checkpoints remain under `data/bounded_risk_2026-11-14/screen/`.
The frozen source archive and logs are under `execution_notes/` beside it.
These are local and ignored by Git. Only aggregate traces and tables are pushed.

## Resume commands

Read the original [prompt](../prompts/2026-11-14-codex-bounded-risk-prompt.md)
and [registration](../plans/2026-11-14-bounded-risk-supervision-plan.md).
Run the following two partitions, at most two concurrently, each with one
BLAS thread. They reuse complete cases and saved arms from partial cases.
Redirect to new log filenames so the interruption logs remain intact.

```sh
/usr/bin/time -p env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript results/benchmark_sim/run_bounded_risk_screen.R normal sim \
  > data/bounded_risk_2026-11-14/execution_notes/normal-sim-resumed.log 2>&1
/usr/bin/time -p env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript results/benchmark_sim/run_bounded_risk_screen.R laplace sim \
  > data/bounded_risk_2026-11-14/execution_notes/laplace-sim-resumed.log 2>&1
```

After both exit successfully:

```sh
Rscript results/benchmark_sim/run_bounded_risk_screen.R summarize
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript tests/audit_bounded_risk_screen.R
```

Check all 24 cases, logical arm and actual fit denominators, guards, cap use,
recalibration passes, stable convergence and matching source hashes. Report
synthetic convergence only. The S-9 strong-signal prior is deliberately
misspecified, and those fits do not determine the real-data pair gate.

Then complete the [findings draft](../findings/2026-11-14-bounded-risk-findings.md)
and [review handoff draft](2026-11-14-bounded-risk-handoff.md), including synthetic
tables, case/fit timing, the complete audit and reproduction ledger. Commit and
push Step 1 results. Get a fresh final independent findings review, address its
issues, finalize the documents and indexes, push, and stop. No follow-up study.

## Verification and reviews already performed

The latest full suite passed 879 counted checks with zero failures, plus new
standalone assertions. The existing bounded-beta optimizer profile-likelihood
fallback warning is recorded. The default S-16 fit fixture passes exactly.
No update algebra changed. The final suite log is
`data/bounded_risk_2026-11-14/execution_notes/full-suite-reviewed.log`.

Fresh code review cleared Step 1 before fitting. An independent interim
findings reviewer reconstructed all 36 real trace summaries and gates before
reading the drafts. Its verdict has no blocking real-data reporting issue.
Final synthetic, loading-tracking, audit and provenance review remains pending.
Private reports are `code-review.md` and `findings-review-interim.md` under
`execution_notes/`.

Keep the caveat that a fixed prior is shrinkage, not a hard coefficient bound.
Canonicalization changes q(β) while the fixed prior stays unchanged, so an exact
standardized-coordinate CAVI derivation is not established. The achieved
moderate shares and single tracked loading correlations do not substitute for
the unrun full program-structure or predictive evaluation.

Do not edit protected shared status files. Proposed text is already in the
findings draft. Preserve the requested pause history in the final reports.
