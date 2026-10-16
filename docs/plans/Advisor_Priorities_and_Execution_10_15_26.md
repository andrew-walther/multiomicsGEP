# Advisor priorities and revised Round2 execution plan

> **For agentic workers:** Use `superpowers:executing-plans` to implement the approved tasks in this document. Work directly in the existing checkout; do not create a worktree or delegate. Andrew approved this order on2026-10-15. Individual task checkboxes record implementation; planned analyses are not completed results.

**Goal:** Produce clear advisor-facing evidence on survival supervision, feature-count sensitivity, parsimony and multimodal value while retaining the approved interim write-up and benchmark.

**Architecture:** Separate benchmark sensitivity from experimental method development. First report the evidence already available, then describe all eligible CpGs and vary nested feature sets at fixed K=7, one modality at a time. Revisit pruning and supervision remedies as distinct experiments so changes in rank, features and likelihood formulation are not confounded.

**Tech stack:** Existing R fitters, loaders, test suite and plotting tools; Quarto with the working external Pandoc configuration; existing Longleaf SLURM conventions for large fits.

**Specification:** Andrew's resumed-work messages and screenshots on 2026-10-15; `Meeting_Notes_Followup_10_2_26.md` (meeting conclusions and Workstreams A/B/C); `Codex_Handoff_Round2_10_16_26.md`; `Codex_to_Claude_Round2.md`; `paper/writeup/_process/implementation_plan.md` (R1–R9 revision dependencies).

## Current status: what is and is not solved

| Advisor question | Verified evidence | Status |
|---|---|---|
| Does survival inform F? | The opt-in canonical Single update uses A_gen+w_S A_surv and B_gen+w_S B_surv, w_S=1; this is actual supervision of F. The evaluated benchmark still has alpha_F=0. | Implemented experimentally; a stable replacement meeting all criteria is **not established**. |
| Is supervised recovery better? | The completed comparator panel recovers17/20 and19/20 planted low-variance programs. The latest profiled panel gives Single18/20 and20/20 and multimodal16/20 and19/20. | Recovery/false-positive evidence is positive; some convergence requirements fail. Simulation success alone does not validate the complete method. |
| Has external performance been preserved? | The earlier selected supervised Single arm had mean external C=.620725 vs frozen D4 .626699; delta=-.005974,95% interval[-.031708,.020662]. | It **fails** the predeclared lower-limit>- .02 requirement. No Round2 external scoring has occurred. |
| Do more genes work? | Round1 requested1000/2000/3000/5000 features per cohort, yielding662/1365/2064/3600 shared genes. Mean CV C=.650200/.648933/.651199/.655413; capped folds1/0/1/2. All four full fits converge. | Limited evidence exists, with a documented nonzero runner exit after saved tables. No full-range experiment, 10k/all-gene result or completed factor-agreement curve. |
| Have we supplied the requested methylation histogram? | Existing plots show pooled values and per-probe skewness for screened10k CpGs, beta and arcsine scales. | **Completed in Task5-QC1:** all305492 TCGA probes on144 patients, observed-value means/SDs and requested histograms; no filter adopted. |
| Does real multimodal pruning select a smaller K? | The latest full supervised-framework candidate remains K7; the frozen framework also removes none. Protected/full-state pruning fails adoption criteria in simulation. | **Unresolved**. K7 may appropriately represent expression, methylation and prognostic structure. Its adequacy needs role, stability and redundancy evidence; deletion itself is not the objective. |
| Does the new multimodal method work? | Algebra/default-identity tests, replicated simulations, isolated failures and exact failed-fold tracing are available. Profiling produces nonconstant risk in all five candidate CV folds. | Partly validated; all five CV fits and full training still cap, four CV folds lack credible programs. Not a fully validated final method. |
| Does methylation add value? | The write-up reports methylation-dominant structure and shared immune loading, without demonstrated predictive improvement; the discussed methylation components are not prognostic. | Added molecular structure is descriptive evidence. Better survival explanation/prediction and robust additional biological value remain unproven. |
| Are priors alone balancing all terms adequately? | Scale freedom, beta gating, shared-prior absorption and constant-direction working curvature have been investigated. Sparse loading priors do not automatically balance total likelihood contributions as p grows. | **Unresolved**. No universal balancing remedy is adopted. The prior/KL consistency question needs separate review. |

Engineering issues actually addressed include positive-SD protection (dead point-mass programs are not protected), opt-in retention of full posterior state during pruning refits, comparator handling of rank0/rank1 and isolated failures, exact saved-fold reproduction, and default/off identity checks. These are useful corrections and diagnoses, not substitutes for the unresolved scientific goals above. Current full suite631/631 at f0858ec; counts are software evidence, not a claim of scientific correctness.

## Global constraints

- Continue on `round2/supervised`; no pushes, protected-branch edits or merges.
- Preserve tag `benchmark-interim-writeup-2026-10-14` and the snapshot. Verify its10 SHA256 entries before and after work; keep benchmark rows in every new real-data result table.
- Preserve the approved `paper/writeup/writeup.qmd`, HTML and PDF until a reviewed revision is ready. The initial unified write-up is already delivered; this work supplies evidence for its revision.
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

## Approved order and deliverables

### 1. Advisor-facing account of progress since 10/2 — first

**Files:** Create `docs/reports/Advisor_Progress_since_10_2_26.qmd` and its rendered PDF; update the Round2 handoff and progress chapter. Do not add process records beside the write-up deliverable.

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

- [ ] Freeze the approved historical framework settings and a fixed K=7/no-deletion sensitivity arm. This is a labelled fallback experiment; its existing convergence limitations remain part of the result.
- [ ] First evaluate methylation alone over the original meeting checkpoints10k/25k/50k/100k/all eligible, subject to actual availability and a prospectively frozen eligibility/filter rule. Do not introduce the mean filter into an otherwise count-only experiment without a distinct predeclared comparison.
- [ ] Then vary expression in joint fits with CpGs fixed at10k; vary CpGs over Task5's10k/30k/100k/all with expression fixed. Matched expression-only/methylation-only/joint comparisons must use the same TCGA patients, folds, expression count, K and controls.
- [ ] Report the same C, partial-likelihood, convergence, feature/program retention, matching and resource measures as the expression study. Record expression/methylation reconstruction separately; raw objectives across different p are not a predictive comparison.
- [ ] Prepare Longleaf arrays for infeasible local cells and leave them explicitly unrun. Do not silently downsample, shrink the maximum universe or alter convergence caps.

**Completion:** Per-modality feature curves followed by comparable joint results. This follows the advisor's sliding feature-count experiment; it does not substitute another K-selection sweep.

### 5. Supervision and pruning remedies — bounded and distinct

**Files:** Existing `run_round2_prior_consistency_trace.R`, attempt log and method-specific derivations/tests; any actual remedy gets its own prospective entry and commit.

- [x] Run the already prepared observational audit only after logging and committing its exact design. It must reproduce the complete saved fit and does not repair the model. Report its finding briefly; do not let diagnosis become an open-ended substitute for the feature deliverables.
- [ ] At fixed settings measure per-loading reconstruction/survival precision and signal contributions, coefficient gating and scale changes; report total objective/removal-score components separately. A raw global likelihood ratio is not the per-coordinate F-update balance.
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

### Full-CpG preparation checkpoint

Task5-P1 (da9bce1/aff0d03) reconstructs305352 eligible probes on144 patients
with exact frozen10k values/ranking and no warnings. Both supplied descriptive
masks are available:182934 mean-rule,221142 proportion-rule probes. The
30-cell methylation-only count design through all probes is committed at
aff0d03. Its first10k full fit capped at500 with zero credible coefficients,
training C=.612671 and nonconstant risk. The25k full cell stopped at the F-update moment-validity guard after28s.
No usable fit/score. The approved unchanged-settings diagnostic now
identifies global exponential underflow in Cox risk sets at outer6; valid
projection moments, but rapidly expanding predictors. No remedy has been
run. CV and larger counts remain unrun; CpG count correlations are unavailable.

Andrew reiterates the intended benefit: broader inputs can reduce arbitrary
screening if prediction and meaningful molecular structure hold up, allowing
Bayesian weighting/shrinkage rather than enforcing a small top-variance set.
Distinguish successful input reconstruction from demonstrated model-level
performance or validated CpG selection. Evaluate both supplied filters
separately; retain current assay QC and probe-availability rules and introduce
no unspecified SD cutoff or new default. Report agreement/roles descriptively;
no new correlation threshold becomes an adoption criterion.

### Reporting update authorized

Andrew requests applicable sensitivity findings in the larger write-up as
well as the progress report. Add completed expression performance/loading
agreement and partial CpG results/diagnosis with explicit benchmark context.
Only10k methylation has a saved full fit, so count agreement is unavailable;
do not substitute saved joint-fit agreement for that comparison. This is a
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
