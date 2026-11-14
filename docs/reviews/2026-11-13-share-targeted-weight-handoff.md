# S-16 review handoff: missing historical precision states

Date: 2026-11-13. Author: Andrew Walther.
To: the reviewing Claude Code session.

Read the [registration](../plans/2026-11-13-share-targeted-survival-weight-plan.md),
the [execution prompt](../prompts/2026-11-13-codex-share-targeted-weight-prompt.md)
and the [findings](../findings/2026-11-13-share-targeted-weight-findings.md).

## Check first

Execution stopped during Step 1, before any S-16 study fit. S-13 and S-15
reduced checkpoints do not save τ. The registered E[L]² substitution for
missing E[L²] does not supply τ. Verify the actual field lists and source
save calls before considering any new study. A fresh independent review
confirmed these omissions and agreed that the no-new-fits historical check
is blocked.

Also check the definition of directional share. The registered mean-τ
denominator is implemented literally. With heterogeneous gene precisions,
actual expression curvature in direction f̂ instead uses Σ_j τ_j f̂_j².
No formula has been changed on this branch. This distinction needs a
decision before calling the registered share exact curvature.

The intra-sweep trace and final saved moments describe different states.
Do not interpret their disagreement as an algebra failure without the
unsaved working weights and intermediate moments.

## Branch and commits

The branch is `s16-share-targeted-weight`, based on `main` at `7b5ea2c8`.
It must remain unmerged pending review. The final handoff commit is the
branch tip and contains this file and its findings note.

| Step | Commit | Coverage |
|---|---|---|
| Step 1 | `6531b160` | Weight helper, tests, no-fit inventory and available numeric audit, README usage |
| Steps 2 to 5 | None | Stopped before study fits |
| Step 6 | None | Explicitly excluded pending review and approval |
| Step 7 | None | Not started |

**Delivery blocker:** automatic approval review rejected `git push -u origin
s16-share-targeted-weight` twice. The second attempt followed a read-only
GitHub ownership check. The authenticated user and repository owner both
have login `andrew-walther` and ID 69438120, and repository metadata report
push and admin permissions. Automatic review still requires trusted user
authorization establishing the destination for this payload. No alternate
upload was attempted. The commits are local until Andrew explicitly approves
pushing this code, tests, documentation and aggregate diagnostics to
`https://github.com/andrew-walther/multiomicsGEP.git` on
`s16-share-targeted-weight`.

## Registered checks and verdicts

| Check | Verdict | Number or evidence |
|---|---|---|
| Step 1 historical state available | Blocked | S-13: 839/840 reduced states, 0 with τ. Each S-15 rank: 36/216 states, 0 with τ |
| Step 1 per-gene trace agreement to rounding | Not fully established | Three seed-1 S-9 examples agree at 3 decimals. Maximum absolute difference 1.5761 × 10⁻⁵. One fails at 4 decimals |
| Step 1 both closed forms within factor 2 in every required cell | Incomplete | Three checked S-9 cells have maximum directional fold error 1.135 and per-gene error 1.649. K = 5 bounds also within factor 2. Missing S-13/S-15 values prevent every-cell claim |
| Step 1 weight inversion | Software checks implemented | Hand-calculated directional ratio and per-gene share, target inversion, guard, cap, active-program selection, rescaling and malformed inputs |
| Step 1 recalibration | Not run | Zero S-16 recalibration passes |
| Step 2 convergence and post-ramp monotonicity | Not assessed | Zero S-16 diagnostic or stabilization fits |
| Step 3 null mean FA ≤ 1 and FA risk variance < 1% | Not assessed | Zero S-16 screen datasets fitted |
| Step 4 ≥ 4/5 converged folds and converged full fit | Not assessed | Zero S-16 training-CV fits |
| Step 4 selection by CV partial likelihood, weaker supervision within 0.5 units | Not applied | No S-16 candidate selected |
| Step 5 structure and annotation | Not run | No eligible S-16 full fits |

The recovery, false-active and convergence gates remain required. The
software tests and historical arithmetic checks do not satisfy them.

Final software verification: 879 counted checks passed, zero failed. The
new `stopifnot()` behavioral assertions also passed within the suite but
are not included in that count. Standalone and actual-Cox integration
checks passed. The existing bounded-β optimizer profile-likelihood fallback
warning was surfaced. `git diff --check` passed.

## Reproduce and inspect

From the repository root:

```sh
Rscript tests/test_share_targeted_ws.R
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript tests/run_tests.R
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript results/benchmark_sim/audit_share_targeted_ws.R
```

The last command needs the existing local RDS caches and
`data/yfb_supervised_real_2026-11-04.rds`. It uses training data only and
writes aggregate outputs to
`results/benchmark_sim/outputs/share_targeted_ws_2026-11-13/`:

- `checkpoint_inventory.csv`: field counts from all 1,866 checkpoints in five directories.
- `saved_share_diagnostics.csv`: three matched S-9 programs and three active K = 5 S-11 programs.
- `audit_run.csv`: elapsed seconds and zero new fits for the audit itself.

The full test suite runs its existing small model checks. These are software
regressions, not new S-16 study observations. Three zero-SD warnings in the
historical audit come from empty loading columns when matching the planted
program. They are surfaced and nonfinite correlations receive zero as in
the existing scorer.

No fitter file was edited, so no update algebra changed and the derivation
verification requirement was not triggered. No opt-in stabilization control
was added and no default changed. Exact-times Cox calls for new study fits
were not needed because none were run.

## Judgment calls and unresolved decisions

1. Stopped at missing S-13/S-15 τ rather than estimating τ, substituting S-11,
   or refitting historical models. The prompt requires reporting impossible
   requirements and stopping, and Step 1 forbids new fits.
2. Calculated a bounded subset of available checks for the handoff: S-9
   seed 1, p = 2,000, unit weight at each variance, and the K = 5 S-11
   full-data unit-weight bound. The subset is explicitly incomplete.
3. Reported three-decimal trace agreement without redefining the registered
   rounding criterion. The registration does not specify decimal precision.
4. Kept the registered mean-τ formula despite its difference from the
   direction-weighted precision. Adoption and formula revision remain for
   Andrew and the reviewing session.
5. No recalibration behavior was chosen for disappearance of pilot activity
   after a later pass. The helper guards on any no-active input. A future
   runner needs a decision before using it beyond the pilot.
6. The helper restricts its supported input to the pooled unweighted canonical
   path. Training centers must be applied consistently before fitting and
   calculating weights. It rejects identifiable controls outside this scope
   instead of silently using the wrong survival design.
7. Used the requested dated document names and November 13 registration date.
   Shared status and manuscript edits are deferred exactly as requested.

No cap, target, selection rule, null gate, program rule or sweep limit was
changed. No study was rerun and no follow-up was started.

Questions for Andrew after review:

- Is a registration amendment authorizing historical reconstruction fits
  acceptable, or should the unavailable S-13/S-15 comparison be waived?
- Should the directional expression precision use Σ_j τ_j f̂_j², or should
  the mean-τ definition remain an explicitly approximate targeting metric?
- Should future trace verification save the exact update-state ingredients
  or compare consistently timed final-state diagnostics?
- What should happen if an active pilot loses all activity during recalibration?

Uncertain quantities are the missing exact S-13/S-15 shares and trace inputs,
not values silently estimated in the output. The K = 5 values are bounds.
The inventory search cannot exclude unlinked checkpoint copies outside the
documented paths. The full every-cell approximation verdict is unknown.

## What is ready for Chapter 4

The missing-state finding and the statement that the procedure has not yet
been evaluated are ready for review. The findings include a labelled draft
paragraph. Available historical shares may be cited only with their subset,
state-timing and bound qualifications.

No claim about a selected target, improved convergence, null calibration,
new real-data programs, training-CV improvement or external performance is
ready. Do not update the chapter or shared status files until this review.

## Independent review of this partial deliverable

A fresh Sol agent reviewed the code and both documents against the plan,
checked representative actual checkpoint field lists, compared the tables
with their CSV sources, and ran the focused behavioral script. It found no
blocking code or reporting defect. Approval was conditional on the final
full-suite result, exact cap/guard assertions and completed commit records.
It agreed that continued execution is blocked. Its mathematical concern
about heterogeneous τ is included above.
The final tests and exact cap/guard assertions subsequently passed, and
the Step 1 commit record is complete.

The independent review did not rerun the entire inventory, numeric audit
or full suite. This handoff preserves the separate review of scientific
adoption and any authorization to amend the registration.
