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

- [ ] Freeze current model settings, patient set, preprocessing, ranking and training folds. Publish a nested feature-set manifest and actual counts before fitting. Reconcile Task5's current2064→5000→10000→all-shared targets with per-cohort screening/intersection; stop for a decision if the existing ranking cannot produce those sets without a new rule. Never describe5000 requested per cohort as5000 shared genes.
- [ ] Fit the increasing nested sets with **K=7 and pruning off**. Include the current screened baseline and the maximum eligible shared set, reporting when a nominal target exceeds availability. Keep model, folds, initialization protocol and convergence rule fixed.
- [ ] Use Task5's continuation protocol; initialize new loading rows from current activities by least squares, preserving existing rows, then let all rows update. Test default cold-fit identity and the prescribed small-example objective comparison. Report rather than hide any warm-start/cold-start failure. Cold-fit the largest set as the prescribed comparison.
- [ ] At every step report5-fold training-CV C and partial log-likelihood, convergence/caps, runtime/memory, active/retained program classifications and loading shrinkage summaries.
- [ ] Match programs by training-activity correlation; report the matching matrix, loading correlations on common genes, projection-score/risk agreement and coefficients only in aligned units. Flag ambiguous or constant factors rather than forcing an interpretation.
- [ ] Plot C and factor agreement versus **actual** feature count. Do not choose a winner using external C; preserve the specified final external-report gate.

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

- [ ] Run the already prepared observational audit only after logging and committing its exact design. It must reproduce the complete saved fit and does not repair the model. Report its finding briefly; do not let diagnosis become an open-ended substitute for the feature deliverables.
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
