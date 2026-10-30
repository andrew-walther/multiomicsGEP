# Execution report: normalized-projection survival supervision

Date: October 30, 2026. This report covers the authorized execution prompt only.
The section beginning “After Codex: review prompt for Claude” was not executed.
The independent review verdict remains for the subsequent review session.

## Outcome

The cleanup audit and preprocessing records are complete. An opt-in normalized
YF fitter, independent numerical checks, registered paired simulation attempts
and an advisor-facing chapter have been implemented. The scientific outcome
and exact failure counts below are populated from the completed saved results.
Implementation completion is not a claim that the supervised method works.

The corrected shared-loading grid completed all **360/360** fits, with
**0 raw numerical-error checkpoints**,
**0 lost-variance records**, and
**153 nonconverged fits**. **Gate S failed:** none of the
three share-target arms qualifies in the required two large-feature cells.
The original criteria are unchanged. Earlier grids are retained as numerical
and design diagnostics because their risk and reconstruction directions used
different units.

| attempt | saved | missing | task_level_errors | raw_checkpoint_errors | variance_loss_records | sweep_cap_records | nonconverged | undefined_C |
|---|---|---|---|---|---|---|---|---|
| attempt1 | 80 | 280 | 56 | 4 | 0 | 1 | 285 | 284 |
| attempt2 | 360 | 0 | 0 | 149 | 0 | 12 | 161 | 149 |
| attempt3 | 360 | 0 | 0 | 155 | 0 | 6 | 161 | 0 |
| attempt4 | 360 | 0 | 0 | 0 | 53 | 115 | 140 | 0 |
| attempt5 | 90 | 270 | NA | 0 | 4 | 22 | 293 | 270 |
| attempt6 | 360 | 0 | 0 | 0 | 0 | 129 | 129 | 0 |
| attempt7 | 360 | 0 | 0 | 0 | 0 | 153 | 153 | 0 |

Missing files, fitting exceptions, variance loss and caps are distinct and can
overlap in the conservative nonconvergence total. Task-level status is
unavailable for the deliberately interrupted fifth driver; its stop reason
and complete missing-record denominator are saved. Undefined C is counted
explicitly. Partial-state C does not establish a successful fit.

| attempt | all_zero_returned_states | nonfinite_returned_eta | inactive_factor_total | orphan_factor_total | unavailable_factor_states | warning_records | objective_decreases |
|---|---|---|---|---|---|---|---|
| attempt1 | 0 | 0 | 261 | 0 | 228 | 4 | 479 |
| attempt2 | 0 | 0 | 431 | 45 | 149 | 60 | 1325 |
| attempt3 | 0 | 0 | 882 | 302 | 0 | 60 | 2951 |
| attempt4 | 0 | 0 | 781 | 267 | 0 | 124 | 13708 |
| attempt5 | 0 | 0 | 238 | 82 | 225 | 31 | 3318 |
| attempt6 | 0 | 0 | 844 | 305 | 0 | 112 | 16600 |
| attempt7 | 0 | 0 | 823 | 293 | 0 | 145 | 16862 |

Warning-record counts are cumulative within datasets in the original driver;
they are not independent warning-event counts. Zero/nonfinite counts concern
returned states; unavailable states are explicit. Objective decreases concern
the monitored approximation and cannot be described as monotone exact ascent.

| variance | arm | recovered | baseline_recovered | recovery_gain | mean_false | converged | qualifies |
|---|---|---|---|---|---|---|---|
| 0.0025 | share0.05 | 0 | 0 | 0 | 6.8 | 2 | FALSE |
| 0.01 | share0.05 | 0 | 0 | 0 | 6.8 | 1 | FALSE |
| 0.05 | share0.05 | 0 | 1 | -1 | 7 | 1 | FALSE |
| 0.0025 | share0.2 | 0 | 0 | 0 | 7 | 0 | FALSE |
| 0.01 | share0.2 | 0 | 0 | 0 | 7 | 0 | FALSE |
| 0.05 | share0.2 | 0 | 1 | -1 | 6.8 | 1 | FALSE |
| 0.0025 | share0.5 | 0 | 0 | 0 | 7 | 0 | FALSE |
| 0.01 | share0.5 | 0 | 0 | 0 | 7 | 0 | FALSE |
| 0.05 | share0.5 | 0 | 1 | -1 | 7 | 0 | FALSE |

The exact cell and arm verdicts are in `attempt7/gate_cells.csv` and
`attempt7/gate.csv`. All fit-level failures and derived penalties are in the
attempt-specific `metrics.csv`; raw checkpoints are unchanged.

The fixed diagnostic replays match registered loading means exactly. A loading
match can coexist with a sweep cap, unmatched activity and severe distortion
of reconstruction. Feature-specific residual precisions become extremely
small on affected genes; all posterior calculations can remain finite while
joint fitting is unstable. The fixed cases are diagnostic, not extra gate
observations.

| variance | recovered | false_active | converged | decreases | largest_mean_difference | first_decrease |
|---|---|---|---|---|---|---|
| 0.0025 | FALSE | 6 | TRUE | 1 | 0 | 2 |
| 0.01 | FALSE | 4 | FALSE | 17 | 0 | 2 |
| 0.05 | TRUE | 4 | FALSE | 84 | 0 | 2 |

| variance | expected_reconstruction_loss | min_tau | max_loading_energy_fraction |
|---|---|---|---|
| 0.0025 | 0.103143039382764 | 0.00554563962473137 | 0.0148814777294625 |
| 0.01 | 939.106205910707 | 4.16326340616329e-08 | 0.999074191437772 |
| 0.05 | 453803162951241 | 8.7483777770172e-20 | 0.999999999999997 |


The supervision branch remains unmerged. The existing benchmark defaults and
the immutable benchmark snapshot/tag are unchanged. No new external cohort
was scored. Real-data training CV and component ablations are conditional on
passing Gate S, rather than deliverables to run after a failed gate.

Source and saved-evidence checkpoint: **28a1837**. The registered fitter is
**fe7064b**, corrected design **ad331c8**, and replay guard **fe97c04**.

## Plan coverage and evidence

| Item | Disposition | Implementation checkpoint and evidence |
|---|---|---|
| C1 scalar KL normalizer | Inherited complete; retained | 31db892; exact scalar/vector KL checks in the full suite |
| C2 single logarithm | Inherited complete; retained | 6a0b96f; globals and previously saved loader checks |
| C3 benchmark recomputation | Inherited complete; not repeated | Saved `kl_fix_rerun_2026-10-30` outputs and chapter cleanup table |
| C4 comparator orientation | Inherited complete; retained | b85d390; saved orientation rescore, not a new external evaluation |
| C5 test portability | Inherited complete; retained | 1703bc6; versioned fixtures and full suite |
| C6 hazard direction and β variance | Done | 4feec24; hazard regression, exact fixed-prior variance check, saved aggregate variance audit |
| C7 local branch inventory | Done, refreshed | `docs/progress_book/data/2026-10-30-branch-inventory.csv`; nothing deleted |
| C8 Linux-only all-zero case | Deferred as specified | Tracker retains the platform risk; new fitter reports degeneracy explicitly |
| P preprocessing decisions | Recorded | 12ebee8; DECISIONS and audit §5 dispositions; chapter decision table |
| S0 precision and DeSurv audit | Done | Conditional analytical/direct-sum audit; `precision_audit.csv`; exact installed DeSurv revision inspected |
| S1 opt-in fitter and tests | Implemented with documented numerical refinements | 12ebee8, 2be72e5, 06a14bc, 1d644e5, fe3ca0c, dae6eef, fe7064b, ad331c8; focused fixtures and full suite |
| S2 registered decision test | Complete; corrected scientific gate failed | Seven separately retained attempts, registrations and complete-denominator summaries |
| S3 real training CV | Conditional; see final gate | No external selection or unauthorized real-data model choice |
| D chapter and living documents | Updated | October 30 chapter, issue tracker, decisions, roadmap, project status and document indexes |
| Review handoff | Supplied | This report and the plan's handoff section; independent review verdict deliberately unfilled |

The original cleanup work was already on main when this execution began. The
remaining audit was committed separately, fast-forwarded to main, and the
supervision branch was created from that tested state. No worktree was used.

## Registrations, failures and deviations

The attempts log is the controlling chronological record. Original checkpoints
are retained under `results/benchmark_sim/outputs/yfb_supervision_2026-10-30/`.
`provenance.csv` identifies fitting revisions, including the intermediate
live-factor calibration change that was registered before execution but did
not receive a separate intermediate commit before launch.

1. **Initial socket launch:** the sandbox prevented opening the local server
   socket. No fit began. The identical command was launched with local-socket
   permission, without switching to forked workers.
2. **Attempt S-1:** the calibration median included dead factors. Their zero
   reconstruction precision made relative precision undefined; many arms
   could not begin. Missing records are counted explicitly, not omitted.
3. **Attempt S-2:** calibrate only live factors and solve the median-share
   equation in log γ. The corrected calibration exposed numerical posterior
   and constrained-prior optimization failures.
4. **Attempt S-3:** remove the Cox-score start from the unsupervised reference,
   rescale normal-means problems, reduce F-message damping, check both loading
   moments, and return coherent failure states with completed sweep traces.
   Smaller damping alone did not repair the reproduced example.
5. **Attempt S-4:** independent integration identified catastrophic cancellation
   in the exponential posterior. Retain ebnm's fitted prior but evaluate the
   same posterior and KL after canceling the common Gaussian term. Use a
   separate experimental kernel, leaving the historical kernel unchanged.
   A saved constrained-prior optimizer failure is handled by a warned
   profile-likelihood fallback over the same parameter space.

6. **Attempt S-5:** fix tiny-positive activity RMS normalization, remove the
   positive projection-norm floor and identify the free scale of survival-only
   loadings. Defer reference coefficient fitting until after reconstruction.
   This run was stopped after a completed checkpoint showed zero stored
   coefficient variance under a nondegenerate prior; workers finished current
   checkpoints and every unstarted record remains missing explicitly.
7. **Attempt S-6:** carry variances separately, evaluate centered Normal
   mixtures and their KL stably at the fitted prior, and use nonnegative
   variance algebra for reconstruction and prediction. The prior families,
   bounds, weight grid and scientific gate remain unchanged.

8. **Attempt S-7:** an algebra audit found that the simulation's raw Cox
   direction differed from the prepared reconstruction direction. Attempts
   1–6 remain numerical/design diagnostics, not valid tests of shared-loading
   recovery. Generate Cox risk from the normalized prepared reconstruction
   loading F/SD_train and retain that same direction in validation. Preserve
   old generated data exactly; save raw/prepared signal shares and oracle C.
   Register the identical gate and run the full six-arm design again.

The numerical corrections deliberately refine the specification's literal
“KL from the ebnm return value” requirement: direct returned numerical moments
were demonstrably wrong on the saved fixture. Later strong-message fixtures found the same numerical issue in returned
Normal/point-normal moments and marginal likelihoods. The final engine uses
stable analytical posteriors and KL at the same fitted priors for all supported
families. Exponential calculations are verified by independent quadrature;
centered Normal mixtures by conjugate formulas and mixed-support KL. No posterior interpolation or post-hoc prior clipping
is substituted for that identity.

Additional explicit refinements:

- Estimate the bounded point-normal prior by directly optimizing its marginal
  likelihood, then evaluate the posterior at that same g. The requested joint
  slab-proportion bounds are not supplied as a public ebnm fitting argument.
  The bounded model can be fitted, so no inverse-gamma Normal family is adopted.

- Use the exact loading norm for live factors, with a floor only at zero,
  instead of an additive epsilon that breaks exact rescaling invariance.
- Perform activity RMS canonicalization at sweep boundaries, carrying prior
  scales and both moments, so the projection scale stays fixed during the
  gene-wise sweep. β and its prior need no rescaling.
- Fit the loading prior on the vector message and hold it fixed during the
  sequential gene updates. Recompute its KL atomically with the returned q.
- Calibrate a fixed γ separately for each paired dataset from its unsupervised
  reference. The grid denotes common precision-share targets, not one common
  numerical γ across all seeds. For the reduced-damping attempts, γ retains
  the original calibration denominator while actual F damping is smaller;
  attained shares must be reported instead of calling them exact targets.
- The calibrated A1 bed is extended by setting the raw centered signal's
  variance share and regenerating expression-only outcomes. The final correction makes the prepared loading
  shared between reconstruction and risk; the coefficient now has normalized
  prepared-projection units, so C cannot be compared across attempts as a
  numerical-repair effect. Feature offsets
  preserve Gaussian noise while making reconstruction inputs nonnegative.
  This is not the signed A1 input or the original two-modality survival DGP.
- The existing canonical comparator uses the same prepared simulated matrix;
  it does not reproduce the historical real-data preprocessing as another
  simulation arm.
- Interpret the false-active and convergence requirements in each qualifying
  variance cell, rather than averaging successful cells with unsuccessful
  ones. This is the conservative gate interpretation registered in advance.
- Older aggregate trace medians included dead-factor scale floors. They are
  not meaningful loading-balance estimates. Derived final-state shares are
  recomputed for matched and live factors without changing gate metrics or
  overwriting original checkpoints. New diagnostic replays retain factor-wise
  traces and exclude dead columns.

No gate threshold was relaxed. Repeated attempts are not pooled. No external
validation outcome or validation cohort was used to choose a setting.

## Numerical and statistical diagnoses

The saved pseudo-message has a huge pseudo-observation and a large standard
error but a relatively small linear message B=A x. The true posterior is close
to its nondegenerate prior. Direct subtraction of Gaussian terms loses its
variance. Integrating the prior-tilted kernel independently confirms the
variance and provides a stable calculation of its KL. The earlier compiled
approximation also has smaller cancellation error, which is why it is not
used as an unquestioned oracle.

With slab probability π and rate λ, define

I = ∫₀∞ λ exp((B−λ)θ − Aθ²/2) dθ,

N = 1 − π + πI.

The exact posterior KL is

KL(q ∥ g) = B E[θ] − A E[θ²]/2 − log N.

The exponential kernel evaluates these terms without the large common normal density,
using tail expansions when necessary and a Gaussian entropy form for a
well-separated positive slab. The prior remains the ebnm marginal-likelihood
estimate. For centered Normal mixtures the exact posterior KL is the Bernoulli
spike/slab KL plus posterior slab probability times the conditional Normal KL.
`VBeta`, `VL`, and `VF` retain variance even when the second moment rounds to
the squared mean; the predictor variance is computed from nonnegative terms. The independent weak and moderate quadrature fixtures test actual
moments and KL, rather than output dimensions.

The constrained β optimizer can fail at a boundary even with finite inputs.
Its fallback profiles over the same bounded slab probability and log SD,
checks the likelihood and emits a captured warning. It is not the
inverse-gamma Normal fallback because the requested bounded point-normal
family can be fitted.

Scientific limitations remain distinct from these numerical repairs. The
projection normalizer is frozen in the working update and recomputed between
sweeps. Precision damping affects the F message while β retains the full Cox
weight. The monitored objective is therefore a diagonal-Taylor approximation,
not a guaranteed monotone bound. Every decrease is reported. A larger γ also
increases coefficient precision; a nonzero slab does not by itself control
false survival activity. The gate evaluates those issues jointly.

## Teach-back: fitter and mathematical helpers

**Purpose.** The experimental engine tests whether survival can change gene
loadings under η = Zβ with normalized Z = centered(Y) F / ‖F_k‖, while retaining
a nonnegative reconstruction Y ≈ LF′. It does not replace the benchmark fitter.

**Inputs and outputs.** Input is a named list containing one n×p nonnegative
matrix, n follow-up times, n event indicators, factor count, optional cohort
and risk-stratum labels, and fitting controls. Output includes n×K activity
moments, list-wrapped p×K loading moments, K coefficient moments, p residual
precisions, fitted priors, n×K projection moments, explicit activity/loading/coefficient variances, n risk values, role labels,
per-sweep diagnostics, convergence status and explicit failure information.

**Logic.** Scale each cohort's features by SD without subtracting means. Build
centered projection data separately. Initialize reconstruction directions and,
for supervised fits, positive/negative Cox-score directions; warm up β. At
each sweep, freeze loading normalizers and log-space Breslow working weights.
Fit the bounded common β prior, then update L, F gene by gene, and β factor by
factor. Carry the within-factor gene coupling through the sequential projection
refresh. Canonicalize activities and inverse-rescale loadings and priors;
update residual precisions from expected squared error. Evaluate the approximate
objective and output changes. Classify the fitted factors without deletion
or reduced-rank refitting.

**Design choices.** A list interface preserves the modality boundary, but only
one block is implemented. Nonnegative data match point-exponential priors;
centering only the projection removes baseline drift. A point-normal slab
bound prevents an absorbing point-mass coefficient prior. Damping acts before
posterior calculation, so q and its KL stay synchronized. Stable numerical
calculations repair the same posterior rather than silently replacing it.

**Gotchas.** Coefficient uncertainty is `VBeta`, not a subtraction of second
moments at extreme means. `VL` and list-wrapped `VF` serve the same purpose
for activities and loadings. Exactly zero activities with surviving loadings
have an unidentified loading scale; the canonicalization fixes that scale
without altering prediction and reports the columns. Numerical precision
floors remain in zero-information update denominators. Zero cohort SD, invalid survival outcomes and nonfinite inputs
raise errors. A collapsed whole factorization, sweep cap or numerical failure
cannot be silent convergence. Numerical failures return the last coherent
state and a message; these partial fits are not successful models. Dead
columns have no meaningful normalized precision share. Fixed-column K roles
reuse the established absolute coefficient threshold from globals; the
registered gate instead uses |E[β]|/SD(β) ≥ 1.96. The role labels are
descriptive and do not validate posterior survival activity. The proposed scaling
uses each validation cohort's distribution and is not a frozen transform for
one future patient. No multimodal extension or external prediction pipeline
has been adopted here.

**R-specific notes.** `ebnm` fits priors and evaluates normal means; `ashr`
constructs point-normal mixture objects; base matrix operations implement
moment algebra; `yaml` supplies existing project constants; Rcpp provides the
sequential gene loop. The isolated exponential kernel uses the same fitted
prior with stable analytical moments and KL. Existing log-space risk sets are
reused with the specified curvature floor. No new package dependency was added.

| Function | Responsibility and key input/output | Important limit |
|---|---|---|
| `fit_yfb_supervised` | Molecular/survival inputs → moments, priors, roles, histories and flags | M=1 only; approximate objective, no adoption claim |
| `yfb_ebnm` | Pseudo-observations and precision → posterior moments, fitted prior, positive KL | Median-SE rescaling; stable exponential evaluation at the same g |
| `yfb_beta_prior` | K coefficient messages → bounded point-normal mixture | Warned same-model profile fallback; bounds are not false-positive calibration |
| `yfb_cox` | Times, events, η and optional strata → score, floored curvature and linear message | Breslow working approximation, not a global Hessian bound |
| `yfb_projection` | Centered Y and loading moments → normalized Z, its variance and scales | Exact scale invariance for every live column; zero columns are reported separately |
| `yfb_normal_posterior` | Centered Normal mixture messages → mean, separate variance and exact mixed-support KL | Retains variance at large means |
| `yfb_predictor_variance` | Projection/coefficient means and variances → nonnegative risk variance | Avoids large second-moment subtraction |
| `yfb_canonicalize` | Moments and priors → equivalent canonical state | Handles tiny activities and survival-only loading scales |
| `yfb_scale` | Nonnegative X and cohort labels → SD-divided matrix | Constant features stop rather than being dropped silently |
| `yfb_stable_exp` | A, B, π and rate → exact exponential posterior moments and coordinate KL | Same posterior model; stable evaluation replaces unreliable returned moments |
| `yfb_stable_compile` | Compile/load the isolated numerical kernel once | Build failures stop; local SDK flags are restored |

## Teach-back: evaluation and documentation

**Purpose and inputs.** The precision audit uses training data and a saved fit
without fitting or scoring an external cohort. The paired simulation uses
fixed seeds, feature counts and variance cells. The summarizer consumes the
registered design and all checkpoints, rather than only successful records.
The chapter consumes saved CSV/RDS evidence through R chunks.

**Logic and outputs.** The audit checks the analytical precision identity by
direct per-gene sums. Simulation generation solves for the planted raw variance
share, offsets expression to be nonnegative, verifies the shared
prepared loading direction, and regenerates expression-only
survival, and prepares each cohort consistently. Each socket task saves every
arm's moments, trace, warnings and metric row. The summarizer creates the
complete expected denominator, inserts explicit missing rows and evaluates the
unchanged gate. Diagnostic replays use fixed saved γ values and are never
added to the gate denominator. Document tables are rendered from saved files.

**Choices and limits.** Socket workers avoid the earlier macOS fork failure.
BLAS thread controls and memory notes prevent unintentional HPC oversubscription.
Synthetic validation C is descriptive and never reoriented using validation
outcomes. The chapter reports C only for converged fits with valid uncertainty and
shows that denominator; all-finite-state partial-fit C remains descriptive in
the CSV. Undefined C is counted explicitly; a conditional mean is not allowed
to conceal its denominator. Lasso is a predictive ceiling, with no factor-count
false-active gate. Saved warning strings in earlier checkpoints accumulated
within each dataset, so a count of records carrying warnings is not a count
of independent warning events.

**R-specific notes.** Base integration supplies an independent numerical
oracle; `parallel::makePSOCKcluster` supplies workers; `glmnet` performs the
training-CV lasso. Saved synthetic checkpoints are separate from the ignored
patient-level training cache. `knitr` renders tables; Quarto renders the book.
All mathematical helpers have documented interfaces and behavioral regression
checks. The base-R diagnostic plot is checked by rendering and visual inspection.
The geometry script reports expected residual energy, precision tails, loading
concentration, projection correlations and end-of-sweep jumps. Those jumps
locate observed events; they do not identify a first within-sweep block loss.

| Evaluation function | Input/output and purpose | Gotcha |
|---|---|---|
| `yfb_decision_data` | Seed, p and variance share → paired synthetic cohorts and scaled truth | Raw variance target, prepared shared loading, fixed validation direction; legacy mode preserves earlier diagnostics |
| `yfb_decision_score` | Fit and truth → recovery, activity, C and convergence metrics | Failed fits remain failures; no validation sign optimization |
| `yfb_calibrate_gamma` | Live reconstruction/survival precisions → fixed γ | Dead columns excluded; calibrated once, not updated within a fit |
| `yfb_decision_task` | One design row → arm checkpoints and task status | Existing checkpoints are reused rather than overwritten |
| `yfb_precision_share_final` | Saved fit plus regenerated training data → per-factor final-state shares | Distinct from original sweep trace; dead columns are NA |
| `yfb_apply_failure_penalty` | Raw records → conservative gate metrics with originals preserved | Non-converged/invalid factor fits cannot count as recovered |
| `yfb_gate_cells` | Full cell summaries → exact cell/arm gate verdict | Requires paired baseline; all boundaries checked in qualifying cells |
| `summarize_yfb_attempt` | Design/checkpoints → metrics, summaries, counts and gate CSVs | Missing and errored rows remain in the denominator |

## Verification and development failures

`Rscript tests/run_tests.R` completed with **742 passed, 0 failed** and a
clean process exit after the final fitting and shared-loading corrections.
The final suite emits the intentional saved-prior optimizer-fallback warning.
The original cleanup sign fixture produced a Cox convergence warning under
perfectly ordered outcomes; the supervision branch uses finite simulated
outcomes and retains the same direction assertion. No real fit is changed
by that fixture refinement.
The complete output is saved as `final_test_suite.log.gz`; `test_results.csv`
records the historical cleanup and final counts. Targeted red/green checks
also verify explicit posterior variances, nonzero projection uncertainty,
canonicalization at zero/tiny activity scales, strict unsupervised starts,
and shared simulation coordinates. The final added check integrates the
conditional posterior for correlated genes independently and verifies both
sequential mean updates and projection-variance refresh. An exact comparison with the historical
generator confirms unchanged legacy matrices, outcomes and targets.

`quarto render docs/progress_book` completed with a clean exit. The full log
is saved as `book_render.log.gz`. The rendered chapter shows the final 742/0
suite, failed gate and saved-result tables. Its figure resource is byte-
identical to the saved diagnostic PNG in both source and rendered resources.
The PNG was inspected visually with common axes across each row. The first
Cairo request failed because the local X11 library was unavailable; the
existing base-R Quartz device produced the checked image without a new dependency.

Rendered main-text checks find no reviewer names, paths, branch names or
placeholders. The previous chapter source, execution prompt, frozen snapshot
and benchmark tag are unchanged. The source/results checkpoint is pushed;
final documentation is committed separately. The chapter is rendered locally
on the supervision branch; this execution does not publish the website.
`closeout_checks.csv` records these checks. The independent review verdict
remains reserved, and the supervision branch remains unmerged.

Development red/green checks intentionally reproduced missing implementation,
the residual-variance matrix-expression error, undefined dead-prior messages,
contaminated unsupervised initialization, a constrained-prior optimizer failure,
and numerical variance loss. A test-isolation bug temporarily caused unrelated
KL failures by overwriting global pi; the test is local and a regression checks
that pi is intact. An intermediate isolation edit had a syntax error, caught
before a clean full-suite rerun. A later manual gate-check wrapper used the
wrong test-state variable name after all four assertions passed; correcting
the wrapper produced a clean 4/4 exit without changing implementation or tests. An audit also corrected a summary implementation error: failed partial fits
receive the registered penalty, and the unintended zero-error restriction
was removed because the gate permits 90% convergence. Raw counts and raw gate
diagnostics are retained. Earlier test/driver parse errors after editing
an actively executing script are reported separately from model failures.
The frozen completed drivers for attempts 2–4, 6 and 7 were resumed from all saved
checkpoints and exited cleanly without refitting. Attempts 1 and 5 are not
filled in with a later engine; their missing records remain explicit.

Quarto's first preview completed all R chunks, then the sandbox blocked its
Sass cache database. Rendering with normal cache permission succeeded. The
final render uses completed saved results. The previous meeting chapter's
source is unchanged. No website publication is performed by this execution.

## Decisions needed next

1. Decide the inferential meaning of γ: a joint generalized-Bayes learning
   rate, or a loading-estimation penalty with separately calibrated survival
   uncertainty. These are different statistical methods.
2. Decide whether to derive a fixed normalized projection objective and an
   acceptance rule, including its normalization derivatives and the role of
   freely estimated feature-noise variances, before another
   registered fitting study.
3. Decide whether survival-only loading directions are admissible expression
   programs, and what identification assumptions support a uniquely
   prognostic-factor recovery claim.
4. Decide whether the planned η=Lβ fallback with fold-in prediction is the
   next registered comparison. No fallback or new default was adopted.
5. Complete the independent review before any supervision merge. Multimodal
   work, real-data model selection and external validation remain downstream.

The standalone problem statement presents the exact failure and distinguishes
these model decisions from the repaired numerical calculations. No remaining
numerical patch is supported by the current failing fixtures. Another search
over weights, prior bounds or screens would not resolve the uncertainty or
objective questions without a new stated model and registration.
