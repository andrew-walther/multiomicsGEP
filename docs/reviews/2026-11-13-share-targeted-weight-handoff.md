# S-16 review handoff: completed fallback evaluation

Date: 2026-11-13. Author: Andrew Walther.
To: the reviewing Claude Code session.

Read the [amended registration](../plans/2026-11-13-share-targeted-survival-weight-plan.md),
[execution prompt](../prompts/2026-11-13-codex-share-targeted-weight-prompt.md) and
[findings](../findings/2026-11-13-share-targeted-weight-findings.md). This replaces
the pre-amendment checkpoint stop handoff.

## Check first

The 50-sweep ramp fails the registered convergence gate. Only raw target 0.01
proceeds under the fallback, and it fails the primary null gate: mean false-active
count 0.2 satisfies ≤1, but mean false-active risk variance fraction 0.193428
exceeds 0.01. No target enters real-data CV. Unit weight is the sole remaining
registered candidate, not a winner among evaluated target shares.

Check null seeds 2 and 10, including the capped nonconverged seed 10 and the
complete denominator of ten. Descriptively removing seed 10 still gives fraction
0.103809 across nine. Keep the registered ten-case failure.

Every real reference fit converges. All five loading correlations exceed
0.99999359. Unit weight has active P2, P4 and P5, with only P2 passing the main
program rule. Two-step also has main P5 because of coefficient activity across
folds. That difference does not establish a new molecular direction.

No model is adopted, no shared status or chapter source is changed, no external
cohort is scored, and no merge into main is authorized.

## Branch and commits

Branch: `s16-share-targeted-weight`, based on main at `7b5ea2c8`. Implementation
and result commits below are pushed to the user-authorized origin
`https://github.com/andrew-walther/multiomicsGEP.git`. Final documents are at the
tip of this branch. Do not merge.

| Work | Commit | Contents |
|---|---|---|
| Original Step 1 | `6531b160`, `22ea81ae` | Original helper and stop audit, superseded |
| Amendment 1 | `5e51438b` | Waiver, gene precision, final states and pilot tracking |
| Amended Step 1 | `564dfe66` | Exact helper, tracking tests, revised audit |
| Step 2 and study implementation | `c294866b` | Sixteen diagnostics, ramp, fallback, tested runners |
| Requested pause | `20a0eef4` | Validated resume checkpoint |
| Step 3 | `3c18d39f` | Complete screen and failed primary null gate |
| Steps 4 to 5 | `883505dd` | Reference CV and structure, empty-target reporting fix |
| Final evidence and review | This handoff's commit at branch tip | Final reports, all-arm loading audit, indexes and review record |

## Verdicts against the registration

| Criterion | Verdict | Number and evidence |
|---|---|---|
| Historical S-13/S-15 comparison | Waived | Amendment 1 item 1 |
| Exact gene-weighted directional ratio | Implemented | Tested heterogeneous τ and Σ E[L²] Σ τ f̂² |
| Final-state ingredients retained | Pass | All new arm records retain τ, E[L²], E[F²], E[β²] and Cox W |
| Fixed pilot tracking | Implemented | Absolute loading correlation, loss of activity retained |
| Closed-form factor-two check | Fails one checked cell | S-9 5% error 2.21308, below factor 4. Exact expression sets weight |
| Step 2 full convergence/decrease gate | Fail | Raw and ramp simulation 1/4 converge. Ramp real 0/2 converge |
| Divergent-case recovery regression | Not satisfied | Ramp schedule tests pass, but no successful convergence fix |
| Raw-fit fallback | Only 0.01 passes | 14 sweeps, no decrease, achieved 0.00873483, weight 82.75367 |
| Primary null mean false-active ≤1 | Pass | 0.2 over ten own-pilot nulls |
| Primary null risk fraction <0.01 | Fail | 0.193428 over ten, nine converge, one capped and nonconverged |
| Secondary paired-weight null | Descriptive failure | Mean FA 0.4, risk fraction 0.285512, ten converge |
| Simulation recovery | Descriptive | Weak 10/10 all arms. Low variance unit 5/10, target and two-step 0/10 |
| Real convergence eligibility | Pass for references | Both 5/5 folds plus full, twelve records, no scoring errors |
| Real selection | Unit weight retained | Sole remaining candidate. Mean CVPL −154.697688 |
| Structure and annotation | Complete for references | Ten program rows, full matrices and matches, forty annotation rows |
| Software checks | Pass | 879 counted checks, zero failures, additional standalone assertions pass |
| Derivation verifier | Pass | 109/109 |
| Step 6 | Not run | Explicitly excluded |
| Step 7 | Not run | Optional wider surface omitted after no target passed null |

Convergence is the existing relative-ELBO flag, not proof of predictor stability
or finite-parameter convergence under separation. The raw S-9 w = 10 and real
target 0.01 stop with η relative changes about 0.00194 and 0.00234. Do not claim
unbounded separation from a finite trace or a stabilized risk score from this flag.

## Evidence and reproduction

Aggregate directory: `results/benchmark_sim/outputs/share_targeted_ws_2026-11-13/`.

| Files | Check |
|---|---|
| `checkpoint_inventory.csv`, `saved_share_diagnostics.csv` | Historical waiver, exact amended audit, labelled K = 5 bounds |
| `diagnostic_summary.csv`, `diagnostics_tracked_programs.csv`, `*_trace.csv` | Sixteen fits, objective changes, drift and tracked programs |
| `screen_metrics.csv`, `screen_null_gate.csv`, `screen_task_status.csv` | All 100 arm records, ten-case gate, nonconverged capped seed |
| `screen_loading_agreement.csv` | One hundred comparisons. Thirty prior target values agree within 4.45 × 10⁻¹⁶ |
| `real_fold_metrics.csv`, `real_summary.csv`, `real_selection.csv`, `real_paired_c.csv` | Twelve fit records, mean CVPL and paired C |
| `structure_programs.csv`, `structure_cor_*.csv`, `structure_match_*.csv`, `structure_folds_*.csv` | Nested-Cox Δℓ, semi-partial R², complete matches and program rule |
| `structure_annotation.csv`, `structure_signature_overlap.csv`, `structure_status.csv` | All four signatures per program and no hidden annotation errors |

The findings give exact commands for all stages and embed the supplementary
no-fit loading-audit script. Private fits, posterior moments, centers, predictions,
frozen drivers and logs stay under ignored `data/share_targeted_ws_2026-11-13/`.
Execution reports and the full independent code review are in its `execution_notes/`.
Final regression output is `tests_resumed_final.log` there.

The screen source is `c294866b`, unchanged through `3c18d39f`. Frozen driver SHA-256:
`c24c0b7b8e0348074dbc0d4bc5d6a8af275556d014cbc2cd9d17f4df6b98fa09`.
Every source fingerprint matched when it finished. The later empty-target real
label fix changes the source hash, so original screen caches intentionally require
the recorded source state for exact cached replay. Fresh current-code tasks remain
possible. No screen fitting method changed.

Real frozen driver SHA-256:
`1e404057abedc70dc07383ad40c6b2b85a259292105c4535aa6864828ad2c2a0`.
Its source matches `883505dd`, and all six real checkpoint fingerprints match.
At most two CPU workers ran, each with one BLAS thread. Real used one worker
while the full suite occupied the other CPU slot.

## Every implementation interpretation to preserve

1. Three passes means three target fits after the unit pilot. Later activity loss
   never resets the selected pilot or invokes the no-signal guard again.
2. The primary false-active fraction is Var(sum of FA components)/Var(total risk),
   averaged over replicates. Signed covariance allocation is also reported. This
   operational definition was fixed before null results. Six guarded nulls have
   no defined target direction, not zero achieved share.
3. S-9 requires |r| ≥0.7 before excluding the matched program from false activity.
   Under null survival all active programs are false. Earlier S-14 counts that
   always excluded the best match were not used.
4. Fixed choices specify w = 0 loadings plus Cox, while Step 5 calls the reference
   EBMF. The implementation uses same-prior canonical CAVI plus Cox. No separate
   EBMF backend was fitted or established numerically equivalent.
5. CVPL is the mean fold contribution, so the 0.5 tie window uses that scale.
   Weaker fold-varying supervision means smaller median fold weight, then target.
   No target survives here, so this ordering does not affect selection.
6. Five finite CVPL and C values are required in addition to the convergence
   criterion. This is stricter than literal eligibility but excludes neither arm.
7. Real CV uses the fixed prepared D4 matrix with fold training centers. It does
   not redo all feature selection and per-cohort standardization inside folds.
8. Simulation train and validation are separately z-scored as in S-13. Nulls
   permute training outcomes only. Validation retains planted survival and is not
   an independent null distribution or a primary-gate input.
9. Multi-target capped aliases are explicit unfit provenance rows, with limited
   numerical reuse semantics. No duplicates occur here. The one capped null fit
   remains in the registered result and gate denominator.
10. The exact compatibility fixture has zero final factors and narrow coverage.
    No successful divergent-case ramp fixture exists. Schedule tests and 879
    counted checks do not establish an effective convergence fix.
11. Only fallback target 0.01 was screened. Targets 0.05 to 0.5 have no simulation
    or CV verdict. Failed large-weight fits do not show insensitivity of stable
    high-share real programs.
12. Optional Step 7 was omitted after the sole target failed null. No rank,
    threshold, target or program rule was changed after seeing results.

## Failures and pause history

Two synthetic diagnostic drivers exited nonzero after saving all four requested
fits because the live top-level script was edited while R read it. One parse
failed, and the other tried to rerun an existing real pilot. All eight states and
trace row counts were checked, with no overwrite. The transient raw traced fitter
was not separately snapshotted, limiting exact forensic source replay. Later
ramp/fallback hashes and frozen drivers were saved. Do not call these clean batches.

The initial screen sandbox could not create its local PSOCK socket and ran no
fit. Local socket permission resolved it. Andrew's usage pause stopped all workers.
The completed weak-seed-1 task was validated and reused, while interrupted unsaved
work was not counted as a scientific failure. The resumed screen exited successfully.

Before real execution, a test reproduced `paste0("s", numeric(0))` inventing an `s`
arm. The corrected helper returns only reference labels for an empty target list.
The test failed before and passed after. This is reporting, not a model change.

## What is ready for Chapter 4

The reviewed convergence limitation, failed null gate, low-variance recovery loss
and reference comparison are candidates for reporting with their qualifications.
The findings provide a labelled draft paragraph and proposed status text.

No new share-targeted model, stronger-supervision program, external result or rank
choice is ready for adoption. Enrichment is descriptive for absolute loading
magnitudes, not signed activation, subtype assignment or survival causation.

## Open questions and uncertain quantities

- Does a future study need a different numerical method before higher target
  shares? The current finite traces do not separate unbounded separation from a
  slow optimization limitation.
- Should future convergence also require risk-score stability? The rule was not
  changed after results here.
- Does Andrew want a separate EBMF backend in a new comparator registration?
  Numerical equivalence was not tested in this study.
- Should future null calibration complement relative risk variance with an
  absolute spurious-effect measure? A modest coefficient can dominate a small
  null risk variance. The current threshold remains unchanged and failed.

Missing S-13/S-15 historical shares are waived. All new-state gate numbers are
available and checked. Exact posterior recovery and unboundedness remain unproven.
The raw transient source has the provenance limitation above. No missing numerical
result was silently imputed.

## Independent review

A fresh code review cleared the raw 0.01 workflow after concrete fixes. A separate
fresh Sol scientific review passed the final fallback reporting with no blocking
statistical or reporting issue. It read the amended registration, code, provenance
and both documents, and independently recomputed the null gate, real means, SEs,
paired differences and loading summaries from CSVs. The fixed-D4 CV limitation
it requested is now stated in both reports.

The reviewer did not refit models, rerun the full suite, rerun the derivation
verifier or execute every R audit. The private report is
`data/share_targeted_ws_2026-11-13/execution_notes/s16-findings-independent-review.md`.
Its verdict approves reporting, not scientific adoption, merge or external scoring.
Andrew's adoption decision remains separate and precedes manuscript, status or
external-scoring changes.

## Reviewing session notes (2026-11-14)

The reviewing Claude Code session checked the tables against the aggregate CSVs and
reran the four S-16 test scripts, which pass. The fitter changes (`survival_weight_ramp`,
`trace_s16`) are opt-in, and the w_S = 1 baseline test confirms the default fit is
unchanged. Andrew approved merging the branch into `main` as a recorded result. Nothing
is adopted. Three points change how the results should be read.

1. **The registered null gate was misspecified (registration error, not an execution
   error).** Under permuted survival every active program is false, so any active program
   carries all of the risk variance, and a 1% threshold amounts to requiring no false
   actives in any dataset. The references fail the same gate by more. Over the same ten
   null datasets: target 0.01 had mean false-active 0.2 and risk fraction 0.19, w_S = 1 had
   0.5 and 0.39, and the two-step fit had 0.4 and 0.23. The failure is not evidence against
   the 0.01 target. The real warning is null seed 10. Its w_S = 1 pilot had two false
   actives, so the no-signal guard did not fire, the weight reached the 10⁵ cap and the fit
   diverged (standardized |β| 28.5). A guard based on pilot activity is not enough.
2. **In simulation the 0.01 target was weaker supervision.** The simulation cells sit near
   directional share 0.5 at w_S = 1, so target 0.01 gave w_S ≈ 0.16 in the low-variance
   cell. Its loss of recovery (0/10 against 5/10) repeats that supervision matters there.
   It is not a cost of share targeting.
3. **Large-weight non-convergence looks like an unbounded objective, not a stalled
   optimizer.** In the real target-0.5 fit the ELBO rose every sweep while max |η| grew
   13, 28, 41, 59 and 85 at sweeps 10, 50, 100, 200 and 400, and max |β| reached 26.
   Canonicalization fixes each program score at unit SD, so β is a standardized log
   hazard ratio, and its prior variance is re-estimated each sweep, so nothing bounds it.
   With 2,218 genes and about 220 training patients, a gene combination can nearly order
   the survival times. This is consistent with overfitting when genes outnumber patients
   (a finite trace cannot prove unboundedness). It explains why the ramp did not help.
   The follow-up is S-17: a fixed prior on the standardized coefficients, the Bayesian
   counterpart of DeSurv's fixed penalty on β.

A replacement for the draft chapter paragraph in the findings note was given to Andrew
in the review session. Codex's draft should not be used as written.
