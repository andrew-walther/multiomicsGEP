# Advisor priorities and revised Round2 execution plan

> **For agentic workers:** Use `superpowers:executing-plans` to implement the approved tasks in this document. Work directly in the existing checkout; do not create a worktree or delegate. Andrew approved the earlier sequence on2026-10-15. The restart priorities below supersede that sequence for planning; new scientific attempts still require prospective registration and any decision required by the original handoff. Individual task checkboxes record implementation; planned analyses are not completed results.

**Goal:** Produce clear advisor-facing evidence on survival supervision, feature-count sensitivity, parsimony and multimodal value while retaining the approved interim write-up and benchmark.

**Architecture:** Separate benchmark sensitivity from experimental method development. On restart, review the saved instability and update equations, then test one justified supervision remedy under the original gates. Resume the remaining fixed-K feature curves after that review; keep rank, features, filters and likelihood changes in distinct experiments.

**Tech stack:** Existing R fitters, loaders, test suite and plotting tools; Quarto with the working external Pandoc configuration; existing Longleaf SLURM conventions for large fits.

**Specification:** Andrew's resumed-work messages and screenshots on 2026-10-15; `Meeting_Notes_Followup_10_2_26.md` (meeting conclusions and Workstreams A/B/C); `2026-10-16-round2-handoff-plan.md`; `2026-10-16-round2-return-note.md`; `paper/writeup/_process/implementation_plan.md` (R1–R9 revision dependencies).

## Start here when work resumes

**State:** The bounded mathematical/evidence review and approved corrected
observer are delivered below for Andrew's review. Stop before any new replay, remedy
implementation or fitting panel. The project objective remains unmet; the
former broad goal is not resumed by this review.
Last scientific/report checkpoint: `0c0c717` on `round2/supervised`; no fitting
jobs are running. The current write-up is a bounded revision of the interim
account; the approved interim tag and snapshot remain unchanged. Do not
repeat the completed expression panel, histograms or25k diagnostic by default.

**North star:** Survival must inform F in a defensible, convergent model,
with a parsimonious set of molecular/prognostic programs, preserved external
performance and a demonstrated contribution from methylation. More fits or
a higher training C are not substitutes for these requirements.

### Recommended goal hierarchy

Use a hierarchy in the plan, with **one bounded active execution goal at a
time**. Keep the north star as the project objective; do not use the entire
remaining research program as the next uninterrupted autonomous goal.

| Level | Goal | Completion evidence |
|---|---|---|
| Project objective | Valid, convergent survival-supervised Single and expression+methylation models with justified programs, preserved performance and interpretable multimodal contribution. | Every original scientific criterion satisfied; proposed configuration reviewed and approved. This objective is not complete. |
| Bounded review delivered | Explain the remaining supervision/optimization failure and specify one justified prospective next experiment. | Equation/evidence review, remedy-ladder ledger, exact observer design with unchanged gates, and decision request are in the current handoff. Stop and regroup before implementation/scientific fitting. |
| Subsequent goal | Evaluate the agreed single-change remedy. | Behavior/algebra/default-identity tests plus the registered simulation/training results; explicit pass/fail by criterion. A failed remedy can complete this bounded evaluation goal without completing the project objective. |
| Later milestones | Complete the remaining feature curves; justify rank/program roles and multimodal value; finalize the reviewed write-up. | The tables, stability checks, paired comparisons and review gates specified below, with unavailable or failed results explicit. |

**Suggested next-goal text:**

> Review the saved supervised-model failures and the F/β update equations.
> Distinguish repaired Cox arithmetic from unresolved working-step, prior,
> scale and convergence issues. Produce a concise evidence-based diagnosis,
> a ledger of tried/inapplicable/untried remedies, and one prospective next
> experiment under the original criteria. Preserve the benchmark and all
> failures. Do not run a new scientific fit. Stop with the review package
> and regroup with Andrew before implementation or a fitting panel.

The existing broad goal remains paused. Do not mark it achieved to switch
focus, or silently resume it while editing these documents. On resumption,
use this smaller scope explicitly; if the goal interface supports editing
its objective, use the bounded text above and retain the project objective
here. Otherwise retain the paused project goal and execute a separately
agreed bounded work block. No goal object is created or changed by this plan.

### Which document controls what?

| Document | Role |
|---|---|
| This plan | Current priorities, next bounded goal, sequencing and restart checkpoints. |
| `2026-10-16-round2-return-note.md` | Concise current status, result tables, task/commit index and review questions. |
| `2026-10-16-round2-handoff-plan.md` | Original criteria, remedy ladders and per-task rules, with later explicit Andrew decisions summarized in the current handoff. |
| `2026-10-16-round2-attempts-log.md` | Append-only exact prospective designs and results; do not rewrite failed attempts. |
| `2026-10-15-round2-history-note.md` | Verbatim former handoff; chronological explanations and historical checkpoints, not restart instructions. |
| `Meeting_Notes_Followup_10_2_26.md` | Advisor motivations, original workstreams and writing requirements; historical branch/prompt sequence superseded for Round2. |
| `2026-10-15-m6-m7-handoff-plan.md` and comparator inventory | Historical Round1 specifications/survey; do not re-execute their branch instructions. |

Direct instructions from Andrew take precedence. Historical plans do not
reactivate completed jobs or override the current pause, branch restriction,
benchmark preservation or unchanged scientific criteria. The dated decision
record remains the audit trail; no dated entry is rewritten by consolidation.

### Restart priorities and decision checkpoints

| Order | Work | Concrete deliverable and gate |
|---:|---|---|
| 1 | Review the current evidence and update equations | One short diagnosis linking the25k extreme-predictor trajectory, Cox working approximation, prior/objective consistency and F-survival balance. Separate established defects from hypotheses; identify the next applicable untried remedy. Regroup before another fitting panel. |
| 2 | Restore reliable survival-supervised fitting | One prospectively registered remedy at a time, first checked algebraically and on small existing fixtures/saved training cases, then judged by the unchanged simulation and training-CV criteria. No replacement unless every required gate passes. |
| 3 | Complete feature-count sensitivity on a stable formulation | Resume methylation-only, then joint one-modality-at-a-time curves through the full eligible universe; report prediction, loading agreement and resources. Preserve count-only and filter comparisons separately. |
| 4 | Justify retained programs and multimodal value | Assess molecular/prognostic roles, redundancy and stability; then matched predictive increment or reproducible additional biology. K7 is acceptable if justified; forcing deletion is not the objective. |
| 5 | Complete the reviewed write-up | Keep verified progress and failures current; integrate final model/biology claims only after their validation and review gates. Preserve the approved interim benchmark throughout. |

The write-up remains the main deliverable at every checkpoint. This order
prioritizes the unresolved scientific bottleneck while preserving the
completed advisor-facing evidence. It does not prescribe a new likelihood,
weight, convergence rule or unreviewed remedy.

### First resumed work block: bounded review, then regroup

- [x] Read the latest checkpoint in `2026-10-16-round2-return-note.md`, the M2/M3/A1
  and prior-consistency entries in `2026-10-16-round2-attempts-log.md`, and the updated
  progress report before implementing or fitting. Verify branch, clean state,
  benchmark hashes and that no prior job remains live.
- [x] Compare the actual F/β updates and priors with the derivations, using
  `code/cox_working_logspace.R`, `code/profile_cox_intercept.R`,
  `code/multimodal_yfb_updates.R`, `code/multimodal_yfb_helpers.R`,
  `code/fit_multimodal_yfb.R` and `code/fit_cox_on_yf.R` as the implementation
  references. Read the multimodal derivation and corrected Single update
  notes first. Review the working approximation versus the intended objective,
  the constant-risk direction, scale transformations of posterior moments
  and priors, and prior/posterior support consistency. Do not assume a new
  kernel alone repairs the optimization algorithm.
- [x] Use saved evidence first: M3 `case.csv` and `terminal_history.csv`,
  M2's ignored failure state, the exact prior-consistency trace, and
  `results/multimodal_real/outputs/round2/loading_balance/`. Small terminal conditional survival
  fractions are descriptive; they do not prove that a weight is needed or
  diagnose the full trajectory. A falling Cox likelihood alone does not prove
  a faulty joint update. Any new execution requires its own committed attempt.
- [x] Produce a short findings table: observation, supported cause, remaining
  hypothesis, minimal proposed check/remedy and the existing ladder rung or
  approval it uses. Explicitly distinguish prior absorption to β=0,
  nonconvergence with nonzero risk, and repaired arithmetic underflow.
- [x] Regroup with Andrew at that deliverable, before a scientific fitting
  panel. Prefer review of a specific derivation or falsifiable next experiment
  over another broad sweep. A fresh session may use this plan and the handoff;
  neither a new session nor clicking Resume changes scientific authorization.

### Bounded review delivered — 2026-10-16

**Execution update:** Andrew approved the observer. Task2-6 reproduced the
complete saved fit, but an observer reconstruction cross-term error invalidates
its full-objective interpretation. The correction is tested; original outputs
are retained with a validity warning. Stop the scientific diagnostic under
its no-retry rule pending a separate replay decision/registration. Continue
the authorized October16 chapter and GitHub Pages publication. See the current
handoff and append-only attempt log; this is not a new model failure or remedy.

**Superseding corrected result:** Andrew approved a separate Task2-7 replay.
It passes exact fit/density-coverage checks. Posterior fixed-prior steps show
no late losses; F-prior refits and working refreshes do. Joint L/beta
comparisons retain unresolved losses. The next recommendation is a reviewed
fixed-q loading-prior EM step with coherent KL, not an implemented remedy or
chosen weight. Read Task2-7/VALIDITY.md for numbers and reference-prior limits;
the full suite passes694/0. Update/publish the meeting chapter and regroup.
The meeting chapter/report is now refreshed with this result. GitHub Pages
built commitb5a413a on the publication branch; source/model branches remain
local and protected branches untouched. Publication and this bounded observer
are complete; regroup before a new remedy. No project-level adoption gate passed.

The [current handoff review package](2026-10-16-round2-return-note.md#bounded-supervision-review--2026-10-16)
contains the equation/evidence table, complete Task2 remedy ledger and one
prospective experiment. It supersedes the earlier generic request to choose
another remedy immediately; no new scientific result or configuration is adopted.

**Diagnosis:** The within-program F residual correction is present. Fixed
quadratic coordinate algebra, prior/objective consistency, and ascent after
refreshing Cox working terms are separate questions. Shared-prior absorption,
per-factor Normal variance-boundary collapse and nonconvergence with nonzero
risk have different evidence. The F-prior audit establishes a support/KL gap,
but its partial monitor correction does not explain late oscillation.
Canonicalization preserves products but does not by itself preserve a carried
shared beta prior. Stable Cox arithmetic repairs M2's underflow, while M3's
expanding predictor and changing reconstruction remain unresolved. Small
terminal survival fractions exclude priors and do not select a mixing weight.

**Ledger:** R1 and the additional approved profiled-intercept remedy failed
adoption; log-space arithmetic is a demonstrated numerical repair only.
Exact R2 (pi floor1/K) and R3 (per-factor Normal with supervised start only)
remain untried, as do longer warm-up, damping, joint/collapsed F–beta and
last-resort survival weighting. The staged-converged-start rung is currently
inapplicable to the saved capped framework fits. No ladder is exhausted.

**Recommended next experiment:** An isolated observer of the complete
working cycle in one unchanged Single seed20261026/noise.25 simulation replay,
compared exactly with saved Task2-4 cell021. This changes observation only,
not a model update. Independently account for current L/F/beta distributions
and priors, coordinate steps, hyperparameter refits, scale transformations,
tau changes and anchored working refresh. Unlike the completed F-only audit,
it must locate loss within a frozen quadratic versus after refresh, or report
why support/coverage prevents localization. Full posterior distributions are
required; moments alone do not determine entropy. Exact controls, arms,
software tests, decision branches and stopping rule are in the handoff.

**Andrew's decision:** Approve or decline this diagnostic extension before
implementation/execution. It requires its own prospective attempt-log entry,
full tests and commit before replay. Execute once, then regroup; no remedy
or fitting panel follows automatically. A positive diagnostic does not waive
the original recovery, false-positive, convergence, rank-when-applicable,
credibility, paired-CV or final external gates. The existing undefined
reference fold remains undefined and prevents the paired-CV gate passing
as written; changing it requires a separate decision. External data remain
unavailable for remedy selection. R2/R3 are still options for the distinct
absorption issue, subject to a coherent shared-prior parameterization.

**This review's scope:** Documentation only, based on code, derivations and
saved aggregate evidence. No remedy, new fit, thresholds/defaults, benchmark,
test bed or external scoring change. Scientific/report checkpoint0c0c717 and
planning basef8c895b are retained. Completion of this bounded review is not
completion of the north star. Verification details are recorded in the handoff.

### Requirements before and after the next remedy

- [x] Reconcile the completed attempt history with the original Task2 ladder;
  list tried, inapplicable and untried rungs with evidence. Do not call the
  entire ladder exhausted. In particular, a staged start requires a genuinely
  converged starting fit; the saved capped framework is not such a start.
- [ ] Follow the existing order: applicable prior remedies, staged fitting,
  damping if oscillation is supported, derived joint/collapsed F–β update,
  then last-resort training-CV survival weight. No new damping value, weight,
  parameter grid or outside-ladder change is chosen in this plan.
- [ ] Before implementing a chosen opt-in change, register its hypothesis,
  one exact change, fixed arms and pass/fail rule; commit before execution.
  Add behavior/algebra tests and protected-base default-identity checks.
- [ ] Keep the original gates: recovery≥16/20 at both noise levels; null-active
  mean≤reference+0.5; every fit converges under its model's rule; median rank3
  when evaluating rank adoption; nonconstant real-data risk and≥1 credible
  coefficient with positive SD; training-CV noninferiority within one paired
  SE. Undefined reference C remains undefined. External reporting occurs only
  at the permitted final gate, with lower95% ΔC limit>−.02; it cannot select
  the remedy. No benchmark replacement without review and Andrew's approval.
- [ ] End each bounded attempt with a result table and a progress decision:
  resolved, a specific failed/inapplicable rung, or needs Andrew. Apply the
  handoff's six-attempt interim-report limit and the existing stop rules.
  Do not start another long panel merely because a process completed normally.

### Deferred work to resume only after that review

- Methylation-only counts remain10k/25k/50k/100k/all305352; joint CpG counts
  remain the separately specified10k/30k/100k/all. Do not silently harmonize
  their test beds. Reuse exact input reconstruction, fixed K7/no deletion,
  original conditional-CV caveats and common-feature matching. Any changed
  arithmetic/source path needs a new prospective design and fingerprints;
  the old driver's hash guard must not be bypassed or its design overwritten.
- Evaluate the supplied mean[.1,.9] and95%-extreme-proportion filters as
  separate frozen comparisons, without selecting an SD cutoff. Report probe
  eligibility/QC gaps, unchanged imputation and the full feasible endpoint.
- Compare expression-only, methylation-only and joint models on identical
  donors/folds. Report credible status, modality-specific reconstruction and
  matched loadings/activities as well as C. Nonzero dense loadings do not
  demonstrate feature selection; comparable C does not establish stable biology.
- Return to rank selection after the fitting/objective issue is resolved or
  reaches its formal stopping point; retain the failed pruning results.
  Characterize justified K7 roles without inventing a desired deletion count.
- DeSurv head-to-head stays deferred by Andrew; optional Task6 robustness
  cleanup is lower priority. Tasks0/4/7 are complete and should not be repeated.

## Current status: what is and is not solved

| Advisor question | Verified evidence | Status |
|---|---|---|
| Does survival inform F? | The opt-in canonical Single update uses A_gen+w_S A_surv and B_gen+w_S B_surv, w_S=1; this is actual supervision of F. The evaluated benchmark still has alpha_F=0. | Implemented experimentally; a stable replacement meeting all criteria is **not established**. |
| Is supervised recovery better? | The completed comparator panel recovers17/20 and19/20 planted low-variance programs. The latest profiled panel gives Single18/20 and20/20 and multimodal16/20 and19/20. | Recovery/false-positive evidence is positive; some convergence requirements fail. Simulation success alone does not validate the complete method. |
| Has external performance been preserved? | The earlier selected supervised Single arm had mean external C=.620725 vs frozen D4 .626699; delta=-.005974,95% interval[-.031708,.020662]. | It **fails** the predeclared lower-limit>- .02 requirement. No Round2 external scoring has occurred. |
| Do more genes work? | Completed fixed-K7 alpha_F0 curve at2064/5001/10000/17900 shared genes; all24 fits converge. Mean conditional CV C=.6084/.5773/.6088/.6068. Loading/activity/risk agreement is reported. | Full-range expression sensitivity complete; no validated gene-selection or noninferiority claim, and program agreement changes. |
| Have we supplied the requested methylation histogram? | Task5-QC1 covers305492 TCGA probes on144 patients, observed means/SDs, SD histograms and filter-overlap counts. Exact model-input reconstruction covers305352 eligible probes. | Completed descriptively; no new filter or SD cutoff adopted. |
| Does real multimodal pruning select a smaller K? | The latest full supervised-framework candidate remains K7; the frozen framework also removes none. Protected/full-state pruning fails adoption criteria in simulation. | **Unresolved**. K7 may appropriately represent expression, methylation and prognostic structure. Its adequacy needs role, stability and redundancy evidence; deletion itself is not the objective. |
| Does the new multimodal method work? | Algebra/default-identity tests, replicated simulations, isolated failures and exact failed-fold tracing are available. Profiling produces nonconstant risk in all five candidate CV folds. | Partly validated; all five CV fits and full training still cap, four CV folds lack credible programs. Not a fully validated final method. |
| Does methylation add value? | The write-up reports methylation-dominant structure and shared immune loading, without demonstrated predictive improvement; the discussed methylation components are not prognostic. | Added molecular structure is descriptive evidence. Better survival explanation/prediction and robust additional biological value remain unproven. |
| Are priors alone balancing all terms adequately? | Scale freedom, beta gating, shared-prior absorption and constant-direction working curvature have been investigated. Sparse loading priors do not automatically balance total likelihood contributions as p grows. | **Unresolved**. No universal balancing remedy is adopted. The prior/KL consistency question needs separate review. |

Engineering issues actually addressed include positive-SD protection (dead point-mass programs are not protected), opt-in retention of full posterior state during pruning refits, comparator handling of rank0/rank1 and isolated failures, exact saved-fold reproduction, and default/off identity checks. These are useful corrections and diagnoses, not substitutes for the unresolved scientific goals above. Latest scientific/report checkpoint0c0c717 has full suite684/684; counts are software evidence, not a claim of scientific correctness.

## Global constraints

- Continue on `round2/supervised`; no pushes, protected-branch edits or merges.
- Preserve tag `benchmark-interim-writeup-2026-10-14` and the snapshot. Verify its10 SHA256 entries before and after work; keep benchmark rows in every new real-data result table.
- Preserve the approved interim in its immutable tag/snapshot. Andrew authorized bounded updates to the working write-up with verified sensitivity/diagnostic findings; final model adoption and biological claims retain their review gates.
- Register and commit each scientific attempt before running. One analytical change per attempt; fixed original criteria/caps/defaults. Never edit a running script or its loaded helpers.
- Training data/simulation only for selection. Keep external outcomes and scores out of feature, model and weight selection. The later Round2 rule of one final external report supersedes the earlier Workstream B request to score ICGC at each feature count.
- The feature-count study is sensitivity analysis, not automatic adoption. Bayesian shrinkage's ability to tolerate many features is a hypothesis to test.
- Run `Rscript tests/run_tests.R` before each commit; stage explicit own files; run `git branch --show-current` first. New code requires tests, README/CLAUDE updates and teach-back.
- Separate requested features per cohort from effective shared features. Use actual available counts, not assumed20k/300k labels. Do not silently replace imputation, preprocessing, probe QC, feature ranking or thresholds.
- No new SD cutoff is chosen from the histogram. Mean/proportion filters can be described using the meeting's specified thresholds; adopting a new filter requires prospective design and approval.
- Large-feature fits may need Longleaf. Prepare memory/runtime estimates and wrappers, but do not submit jobs. Report local feasibility and explicit unrun cells.

## Review focus

1. A fixed K study must keep K=7 and deletion off; shrunken/inactive factors are classified, not silently deleted. Varying K is a later experiment.
2. Feature sets must be nested and restricted to training-derived eligibility/ranking. Report exact counts, IDs/hashes, missingness, QC and preprocessing provenance.
3. Compare program loadings only on common features and patient activities on the same training patients; report undefined correlations for constant factors and ambiguous matches.
4. More CpGs change computation and total reconstruction information. Do not label a larger ELBO or nonzero risk as better prediction, convergence or survival-informed biology.
5. Preserve partial failures: caps, rank-deficient warm starts, memory limits, invalid scores, unavailable raw/QC information and cold-versus-warm discrepancies remain visible.

## Earlier approved deliverables and completion checklist

The restart order above takes precedence for the next work block. These
checklists retain the original scope and document completed versus open work.

### 1. Advisor-facing account of progress since 10/2 — first

**Files:** Create `docs/findings/2026-10-16-advisor-progress-since-10-02.qmd` and its rendered PDF; update the Round2 handoff and progress chapter. Do not add process records beside the write-up deliverable.

- [x] Build the report from saved evidence only, with no new fits. Show the benchmark, the supervised candidate's failed external non-inferiority result, planted recovery, K-selection failures, the limited feature experiment and the current multimodal failure.
- [x] Label each advisor request complete, partial or open. Include the actual main result tables and distinguish historical results from proposed work.
- [x] Explain alpha_F=0 precisely: survival fits coefficients and may affect the returned stopping point, but does not guide F; it is close to factorization-then-survival in this important sense. Experimental canonical fitting has a different F update and is not alpha_F=0.
- [x] Audit the seven joint programs by expression/methylation contribution, credible survival association, overlap of roles and redundancy. Seven retained programs may be appropriate; do not make deletion itself a success criterion. Frozen roles do not establish roles for an unadopted newer fit. Distinguish Single four classified programs from deletion to K4.
- [x] Map available new evidence to write-up revision markers R1–R9. Do not claim a replacement model or add unvalidated biology.
- [x] Render the addendum successfully with the established Quarto/Pandoc workaround and verify the main write-up still matches the benchmark.

**Completion:** A readable advisor-facing PDF and a short list of findings eligible for review/integration. This gives a visible deliverable before another long fitting panel.

### 2. Full-probe methylation descriptive analysis — next, no model selection

**Files:** Create `results/multimodal_real/summarize_cpg_variability.R`, `tests/test_cpg_variability.R` and fresh `results/multimodal_real/outputs/feature_sensitivity/cpg_qc/` aggregate outputs. Read `code/load_multiomics_data.R` and existing distribution checker first; reuse their established scales and data conventions.

- [x] Inventory the available TCGA probes, patient set, raw beta scale, missingness, existing QC and imputation provenance. Reconcile these with the screened10k baseline; do not overwrite it. Inspect shared feature identifiers without external outcome-based selection.
- [x] For each CpG across TCGA training patients compute observed count, mean beta, sample SD/variance, and proportions below .1 and above .9. Also report transformed SD for asin(2*beta-1). Keep pre-/post-imputation statistics distinctly labelled when both are available.
- [x] Summarize counts under the mean rule (.1<=mean<=.9) and the meeting's95% extreme-proportion rule, their overlap and overlap with current10k. Show distributions before and after these descriptive filters.
- [x] Produce the requested **histogram of per-probe SD**, on beta and arcsine scales, and a mean-versus-SD plot. Keep the feature-level inventory local if necessary; commit only approved aggregate outputs and plots, never patient data.
- [x] Report whether cross-reactive/SNP probe filtering is verified upstream. Unknown QC is an explicit gap, not permission to claim probes are clean or silently change the universe.
- [x] Test known means/sample variances, missing values, boundary means=.1/.9, raw-beta range, constant probes, transform-specific SD and inventory/count consistency using software fixtures before the TCGA descriptive run.

**Completion:** Full-probe inventory counts, filter-overlap table and SD histograms. No survival outcomes, C-index or newly chosen cutoff enters this step.

### 3. Expression feature-count curve at fixed K=7

**Files:** Create the specified Task5 runner `results/benchmark_sim/run_feature_expansion.R`, continuation helper/test files after reading existing warm-start interfaces, and fresh aggregate outputs under `results/benchmark_sim/outputs/feature_expansion/`. Write a prospective attempt entry before fitting.

**Approved sensitivity baseline:** The approved alpha_F=0 D4 configuration, not an unadopted supervised remedy. This is benchmark sensitivity and cannot be described as resolving supervision. Andrew approved starting this labelled sensitivity study before the supervision remedies settle.

- [x] Freeze current model settings, patient set, preprocessing, ranking and training folds. Publish a nested feature-set manifest and actual counts before fitting. Reconcile Task5's current2064→5000→10000→all-shared targets with per-cohort screening/intersection; stop for a decision if the existing ranking cannot produce those sets without a new rule. Never describe5000 requested per cohort as5000 shared genes.
- [x] Fit the increasing nested sets with **K=7 and pruning off**. Include the current screened baseline and the maximum eligible shared set, reporting when a nominal target exceeds availability. Keep model, folds, initialization protocol and convergence rule fixed.
- [x] Resolve the prescribed warm continuation explicitly: frozen activity rank4/7 makes its inverse undefined. Andrew authorized cold-primary judgment; all counts use cold starts. Warm continuation and its comparison test are deferred, not implemented or passed.
- [x] At every step report5-fold training-CV C and partial log-likelihood, convergence/caps, runtime/memory, active/retained program classifications and loading shrinkage summaries.
- [x] Match programs by training-activity correlation; report the matching matrix, loading correlations on common genes, projection-score/risk agreement and coefficients only in aligned units. Flag ambiguous or constant factors rather than forcing an interpretation.
- [x] Plot C and factor agreement versus **actual** feature count. Do not choose a winner using external C; preserve the specified final external-report gate.

**Completion:** One expression feature-count table and performance/agreement curves through the full feasible universe, with failures explicit. No K sweep or model-weight change is combined with this experiment.

### 4. Methylation-only curve, then joint curves one modality at a time

**Files:** Create the specified Task5 runner `results/multimodal_real/run_feature_expansion_multimodal.R`, focused continuation/protocol tests and large-feature SLURM wrapper; fresh results under `results/multimodal_real/outputs/feature_sensitivity/`.

- [x] Freeze the approved historical framework settings and a fixed K=7/no-deletion sensitivity arm. This is a labelled fallback experiment; its existing convergence limitations remain part of the result.
- [ ] First evaluate methylation alone over the original meeting checkpoints10k/25k/50k/100k/all eligible, subject to actual availability and a prospectively frozen eligibility/filter rule. Do not introduce the mean filter into an otherwise count-only experiment without a distinct predeclared comparison.
- [ ] Then vary expression in joint fits with CpGs fixed at10k; vary CpGs over Task5's10k/30k/100k/all with expression fixed. Matched expression-only/methylation-only/joint comparisons must use the same TCGA patients, folds, expression count, K and controls.
- [ ] Report the same C, partial-likelihood, convergence, feature/program retention, matching and resource measures as the expression study. Record expression/methylation reconstruction separately; raw objectives across different p are not a predictive comparison.
- [ ] Prepare Longleaf arrays for infeasible local cells and leave them explicitly unrun. Do not silently downsample, shrink the maximum universe or alter convergence caps.

**Completion:** Per-modality feature curves followed by comparable joint results. This follows the advisor's sliding feature-count experiment; it does not substitute another K-selection sweep.

### 5. Supervision and pruning remedies — bounded and distinct

**Files:** Existing `run_round2_prior_consistency_trace.R`, attempt log and method-specific derivations/tests; any actual remedy gets its own prospective entry and commit.

- [x] Run the already prepared observational audit only after logging and committing its exact design. It must reproduce the complete saved fit and does not repair the model. Report its finding briefly; do not let diagnosis become an open-ended substitute for the feature deliverables.
- [x] At fixed settings measure per-loading reconstruction/survival precision and signal contributions, coefficient gating and scale changes; report total objective/removal-score components separately. A raw global likelihood ratio is not the per-coordinate F-update balance.
- [ ] Use the approved remedy ladders and original adoption criteria. Do not add inverse-p weights, estimated mixing parameters, SD cutoffs, longer caps or a new reference without the required decision and prospective design.
- [ ] Evaluate whether adding noise features in the originally specified large-p simulation changes pruning/recovery, separately from true additional signal. This addresses the meeting's B5 hypothesis and does not redefine the original adoption panel.
- [ ] Keep real-data K7/no-deletion results and simulated rank-selection results separate. Do not force deletion merely because fewer factors are preferred.

**Completion:** Each bounded attempt yields a pass/fail table or a precise decision request. Supervision is solved only when the original simulation, convergence, CV and external requirements are all met; no current option does so.

### 6. Establish the multimodal contribution and revise the larger write-up

- [ ] Report two separate questions: predictive increment over matched expression-only, and additional reproducible molecular/biological structure. A matched C is not evidence of predictive improvement.
- [ ] First summarize the existing benchmark biology honestly: methylation-dominant and shared immune structure, with no demonstrated prognostic contribution of those specific components. Compare stability across feature sets and loading support, not only the number of factors.
- [ ] For a final stable configuration, predeclare the permitted matched comparison and paired uncertainty. Score external cohorts only at the original final-report gate; never optimize features or remedies against ICGC.
- [ ] Plan CpG annotation, enrichment and composition/batch checks under the final-program characterization/review process. Novel biological interpretation requires those checks; methylation dominance alone is not value-add validation.
- [ ] Revise the write-up's affected R1–R9 passages after Claude's review and Andrew's approval of the findings/configuration. Keep the benchmark record intact and retain partial/negative findings.

**Completion:** A defensible predictive or biological contribution, or an explicit negative result; a reviewed manuscript revision that distinguishes these outcomes.

## Approved sequencing decision (2026-10-15)

The original Round2 Task5 says Tasks1/2 must settle the configurations first.
The earlier10/2 Workstream B explicitly allowed feature sensitivity with the
current method while supervision work continued. Andrew's latest request
puts these experiments back near the front. Andrew approved using the
frozen benchmark/fallback for a **labelled fixed-K sensitivity study now**,
while the unadopted supervised method remains under investigation. It does
not waive any adoption criterion or claim that alpha_F=0 solves the goal.

Andrew approved this revised order in the pending-question reply before implementation.
Descriptive evidence gathering and this plan do not constitute a new
scientific attempt. Keep the same goal/branch; clicking Resume in the app
reactivates the paused goal. A new chat is optional, using this plan and the
Round2 handoff, and is not necessary to change priorities.

## Execution checkpoints

After each deliverable, update a short advisor-request checklist with links
to actual tables/figures and one of: complete, partial, unresolved, or needs
decision. Report what new evidence was produced, the interpretation it
supports, what failed and the next bounded action. Do not use elapsed time,
fit counts or test counts alone as evidence of scientific progress.


## Expression preflight: available universe and continuation decision

The existing D4 preprocessing yields17900 common genes over273 training
patients. Existing per-cohort combined ranks give the nested shared counts:
2064 (top3000/cohort),5001 (top6607/cohort),10000 (top11602/cohort), and
17900 (all). The nominal5000 target is crossed at5001 because both cohort
rank lists are intersected; no tie-breaking/ranking rule is changed. The
baseline gene order equals frozen D4's; its standardized matrix is exactly
the corresponding subset of the full17900-gene matrix. No model is fitted
in this preflight and no external cohort is loaded.

The frozen D4 activity matrix has rank4 of7; columns1,2,4 are exactly zero.
Task5's prescribed continuation formula therefore has singular L-transpose-L.
Carrying zero columns into new loading rows also risks keeping them dead,
which is consequential when asking whether more features reveal new programs.

**Approved decision:** Andrew authorized using judgment and prioritizing progress; use the same cold SVD start protocol at each of the
four fixed-K7 feature counts as the primary sensitivity curve. Keep the
frozen alpha_F=0 settings, cohort preprocessing, convergence rule and original
training-CV fold protocol. Evaluate full training plus five folds per count
(24 new fits); report C, partial log-likelihood, shrinkage, program roles,
agreement, caps and resource use with benchmark rows. Do not score external
cohorts or choose a weight/feature count from them. Defer warm continuation
until a rank-deficient policy is reviewed. This changes Task5's continuation
requirement with Andrew’s approval before fitting; no adoption rule changes. The prospective performance-panel design must still be committed before fitting.


### Completed expression checkpoint and interpretation

Task5-F1 now completes the expression curve through17900 common genes:
24/24 converged; CV means .6084/.5773/.6088/.6068, full retention4/4/4/3.
Agreement changes materially; no new count/default is adopted. Dense nonzero
loading means do not establish gene-level signal/noise selection. Andrew
emphasizes that full inputs can reduce arbitrary filtering even without a C
gain if performance/stability hold up; this remains a motivation, not a changed
adoption gate or a demonstrated noninferiority result. See the updated HTML/PDF
advisor report and safe cold_primary result tables. Methylation/joint curves
remain open; original shared-probe IDs number305352 versus305492 raw TCGA.

### Historical full-CpG preparation checkpoint (superseded by terminal checkpoint below)

Task5-P1 (da9bce1/aff0d03) reconstructs305352 eligible probes on144 patients
with exact frozen10k values/ranking and no warnings. Both supplied descriptive
masks are available:182934 mean-rule,221142 proportion-rule probes. The
30-cell methylation-only count design through all probes is committed at
aff0d03. Its first10k full fit capped at500 with zero credible coefficients,
training C=.612671 and nonconstant risk. The25k full cell stopped at the F-update moment-validity guard after28s.
No usable fit/score. The approved unchanged-settings diagnostic now
identifies global exponential underflow in Cox risk sets at outer6; valid
projection moments, but rapidly expanding predictors. At that earlier checkpoint no remedy had run. Since then the opt-in stable
arithmetic case completed at500 without convergence; provisional CpG-count
correlations are now available. CV and larger counts remain unrun.

Andrew reiterates the intended benefit: broader inputs can reduce arbitrary
screening if prediction and meaningful molecular structure hold up, allowing
Bayesian weighting/shrinkage rather than enforcing a small top-variance set.
Distinguish successful input reconstruction from demonstrated model-level
performance or validated CpG selection. Evaluate both supplied filters
separately; retain current assay QC and probe-availability rules and introduce
no unspecified SD cutoff or new default. Report agreement/roles descriptively;
no new correlation threshold becomes an adoption criterion.

### Historical reporting authorization (terminal results now incorporated)

Andrew requests applicable sensitivity findings in the larger write-up as
well as the progress report. Add completed expression performance/loading
agreement and partial CpG results/diagnosis with explicit benchmark context.
The initial update had only a saved10k methylation fit. The terminal update
now includes the saved25k stable-arithmetic fit and provisional common-CpG
agreement; the separate joint-fit audit is not substituted for that comparison. This is a
bounded manuscript update, not approval of a new model or altered criterion.
The approved interim remains preserved in the immutable tag/snapshot.


## Approved pause checkpoint: 2026-10-15

Andrew requests pausing after the current25k fit, saved loading comparison
and report/commit. Those bounded runs now terminate: arithmetic crash removed,
scientific fit still capped500 with extreme predictors; provisional common-CpG
agreement recorded. Do not start another fitting experiment before regrouping.
The current progress report/manuscript include the results; benchmark unchanged.
Review working-step stability and prior/objective consistency before deciding
the next prospective remedy. Larger count/CV/filter/joint panels remain open.


### Terminal evidence for restart priorities

The25k case completes at500 with no execution errors but fails convergence:
training C=.983854, riskSD39086, logPL=-40177, final relative risk/reconstruction
changes=.00589/.0622 against tolerance.0001. It does not establish successful
supervision or a useful high-dimensional fit. The original10k also caps.
Common-CpG loading r=.945/.922/.925/.118/.038/-.087/.923; risk r=.245.
Both are capped and numerical arithmetic differs: provisional agreement only.
The numerical underflow failure is repaired opt-in; imputation is not the
identified cause. A stable fully supervised replacement is still unproven.
The original selected supervised Single external noninferiority check failed;
the paused goal's older state description is stale on that point. This plan
and the latest handoff, not that description, govern interpretation.


## Approved implementation and meeting-document consolidation — 2026-10-16

Andrew instructed proceeding with the proposed prior-consistency repair and
updating all three deliverables. Its fixed-posterior point-Laplace prior and
KL correction are now opt-in;12 focused checks and full suite706/0 pass.
No scientific fit has run. The existing tail-moment approximation has a
quantified limitation, so exact global-objective claims remain unsupported.
Use 2026-10-16-round2-fixed-q-loading-prior-proposal.md for the derivation, code
walkthrough, supported controls, exact one-case prospective arms, original
criteria and stopping rule. Before a scientific run, its fixed driver and
final execution registration must be committed. This does not resume the
broad paused goal or authorize an adoption/external panel.

Meeting documents now distinguish strong5,001→10,000 gene program/risk
agreement from less uniform agreement against2,064 or all17,900 genes.
The completed Single training/CV curve is not the still-open matched TCGA
framework expression curve. CpG/filter/joint curves, large-p noise studies,
rank and multimodal/biology validation remain open. October16 is a concise
meeting chapter; internal registration/approval/observer bookkeeping remains
in this plan and the handoff. The long October15 methods log is preserved
in the book, clearly labelled as a methods-development record. Refresh all HTML/PDF outputs
and publish only rendered gh-pages; source branch stays local.

Andrew clarifies the deliverable: a concise advisor update since10/2 plus
the full method/results/biology write-up, supporting manuscript drafting and
journal choice. The new 2026-10-16-advisor-meeting-summary HTML/PDF is the
brief meeting document; the longer report retains detailed evidence.
Journal_Targeting_10_16_26.md under paper/writeup/_process/ records a
provisional publisher-scope shortlist and remaining evidence requirements.
Retain the methods-development chapter in the book with an explicit purpose
note, as Andrew subsequently requested; do not remove it from navigation.

Final publication verification: GitHub Pages reports rendered commit
`c12c3aa36e7f675291111f260fe77c70244f9e80` built successfully. October16 is
concise and separates adjacent-count program agreement from whole-range
agreement; chapter11 remains with its explicit methods-record purpose.
Meeting handout2 pages, detailed report13 pages, full write-up43 pages;
all HTML/PDF outputs render. Local book checks find108 references with no
missing target. Full suite706/0 and focused12/12 under both R/C++ paths;
all10 benchmark hashes match before/after. No scientific fit or external
scoring in this block; source/model branch remains local.


## Current authorized order — advisor preparation, audit, method continuation

Andrew requests a structured overview and a brief route through the write-up
for the upcoming meeting, informed by DeSurv's narrative but supported by this
project's evidence. First render the updated full write-up/meeting summary and
refresh the concise book. Review code/test bloat before methodological work.
The audit distinguishes accumulated branch diff from uncommitted changes;
no evidence/test deletion or history rewrite follows merely from line count.
Then prepare the exact one-case fixed-q proposal's driver, tests and prospective
attempt entry; commit before any execution. Keep its original controls,
one-change comparison, density limitations and stop rule. No broad panel,
new external evaluation, default change or protected-branch/source push.


### Terminal one-case result and next priority

Task2-9 completes once after02a658a. Fixed-q local gains/density spot checks
pass; native ELBO rule stops at10 iterations, but output movement remains
.01915 versus existing1e-4 diagnostic. Both reference/candidate recover the
planted program with1 false-active program; no single-case performance claim
or adoption. The next priority is stopping/objective consistency, before an
expanded panel. Register any new observation/analytical change separately;
no automatic replay, higher cap, stopping-rule substitution or external data.

Final book publication: gh-pages commitb5fe3c322e5fbfc9ed6486e3f33e6c9bcbdcc063
built successfully and live content was verified. Only rendered book content
was pushed; source/model branch remains local. Full tests710/0, benchmark10/10;
working scientific recommendation remains review before another registered case.


## Meeting discussion guide and write-up order

The October16 book chapter now opens with the meeting brief's four discussion
questions: model interpretation, validation/feature sensitivity, biological
value, and the supported manuscript claim/evidence still needed. Existing
result tables remain the progress update; the opening provides a route through
the full write-up without duplicating its detail. Andrew requested restoring
the table of contents before the summary; the PDF now does so, keeping the
short Background/Aim, Method, Results and Conclusions paragraphs together.
This is a documentation-only follow-up; no fit, default, criterion or scientific
recommendation changes. The next methodological priority remains review of
stopping/objective consistency before separately registered execution.

Verification: full suite710/0; benchmark manifest10/10 before/after;
108 local book references resolve. The write-up HTML/PDF render, and the
contents-first PDF/structured summary layout was visually checked (44 pages).
Rendered book publication commit27bae873bfc42436baa76f383cc3956809c47ea6;
source/model branch remains local.


### Meeting chapter executive summary

Andrew requests replacing the October16 opening discussion-guide table with
an executive summary. The chapter now summarizes the write-up's model,
external prediction evidence and provisional biology, then gene/CpG study
progress and the numerical/prior-consistency findings. It retains the
alpha_F0 benchmark limitation and unresolved output stability; no new fit,
criterion or adoption decision follows. The existing detailed result sections
remain below the overview. This supersedes the earlier discussion-guide
layout, not the scientific priorities.

Verification: full suite710/0; benchmark hashes10/10 before/after;
108 local book references resolve. Rendered publication is gh-pages b27b73e.
