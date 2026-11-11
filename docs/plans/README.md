# Plans

Plans, design specs and proposals, written before the work. Results of plans
that ran are in [`../findings/`](../findings/README.md); meeting follow-ups that
set the work are in [`../meetings/`](../meetings/README.md).

Current priorities are maintained in the [roadmap](../../ROADMAP.md).

| Date | Document | Purpose | Status |
|---|---|---|---|
| 2026-11-11 | [results-rebuild-design](2026-11-11-results-rebuild-design.md) | Unified result dependencies, model provenance, current releases, proposed layout and checkpoint retention, built on 0b | design for Andrew's review; no implementation or moves |
| 2026-11-11 | [results-rebuild-script-dispositions](2026-11-11-results-rebuild-script-dispositions.csv) | Proposed current, historical or retired role for all 204 inventoried scripts | proposed mapping accompanying G; not migration authorization |
| 2026-03-25 | [fit-modular-design-spec](2026-03-25-fit-modular-design-spec.md) | Design for `fit_modular.R` (factor-wise CAVI) and its companion walkthrough | complete |
| 2026-03-25 | [fit-modular-plan](2026-03-25-fit-modular-plan.md) | Task-by-task implementation plan for `fit_modular.R`, its walkthrough and the factor-wise simulation report | complete |
| 2026-04-29 | [beta-zero-fix-design-spec](2026-04-29-beta-zero-fix-design-spec.md) | Design for fixing β = 0: changes inside the LB model, or the Cox-on-YF reformulation | complete (the reformulation became YFB, `fit_cox_on_yf.R`) |
| 2026-05-20 | [cohort-lmm-review/](2026-05-20-cohort-lmm-review/) | Literature review of Bayesian LMMs, executive summary and model comparison proposing cohort indicator columns in L (`.qmd` set with its `references.bib`) | complete (implemented and evaluated in [cohort-lmm-benchmark-report](../findings/2026-05-22-cohort-lmm-benchmark-report.qmd)) |
| 2026-05-25 | [merged-cohort-benchmark-plan](2026-05-25-merged-cohort-benchmark-plan.md) | Six-configuration merged TCGA+CPTAC benchmark with K chosen by cross-validation | complete |
| 2026-05-25 | [preprocessing-comparison-k-floor-plan](2026-05-25-preprocessing-comparison-k-floor-plan.md) | Extends the benchmark to 18 preprocessing × model configurations with K ≥ 3 | complete |
| 2026-05-27 | [desurv-aligned-preprocessing-plan](2026-05-27-desurv-aligned-preprocessing-plan.md) | Align gene selection with DeSurv (mean+variance rank, per-cohort top 3000, intersection) and compare D1–D4 | complete (produced D4, the recommended configuration) |
| 2026-06-14 | [multicohort-sim-proposal](2026-06-14-multicohort-sim-proposal.qmd) | Design of the shared versus study-specific multi-cohort simulation, with first results | complete |
| 2026-07-12 | [joint-k-alpha-bayesopt-plan](2026-07-12-joint-k-alpha-bayesopt-plan.md) | Joint Bayesian-optimization tuning of K and α, so the K = 7 vs DeSurv k = 3 comparison uses a comparable search | complete (run as Step 3 of the K-parsimony plan) |
| 2026-07-13 | [ssbmf-k-parsimony-followup-plan](2026-07-13-ssbmf-k-parsimony-followup-plan.md) | Six steps to test whether K = 7 can shrink (multistart, deflation or joint BO, validation, enrichment) | complete |
| 2026-07-14 | [pathway-enrichment-plan](2026-07-14-pathway-enrichment-plan.md) | fgsea/ORA enrichment of the survival-active D4 programs, with subtype and DeSurv concordance | complete (see the 2026-07-15 and 2026-08-21 enrichment reports) |
| 2026-07-14 | [ssbmf-progress-consolidation-and-remaining-work-plan](2026-07-14-ssbmf-progress-consolidation-and-remaining-work-plan.md) | Summary of work since 6/18 and remaining Items 1–6 | complete for Items 1–5; Item 6 (HTML visual summary) not started |
| 2026-07-15 | [ssbmf-item2-executive-summary-plan](2026-07-15-ssbmf-item2-executive-summary-plan.md) | Spec for the long progress report and the short status summary of 7/15 | complete |
| 2026-08-13 | [ssbmf-factor-classification-k-selection-plan](2026-08-13-ssbmf-factor-classification-k-selection-plan.md) | Three-way factor classification and validation of ARD-based K selection | complete |
| 2026-08-20 | [yfb-dual-source-f-experiments-plan](2026-08-20-yfb-dual-source-f-experiments-plan.md) | Whether survival can inform F (α_F > 0) via a non-degenerate β estimator or a β warm-up | complete (null result; α_F = 0 kept) |
| 2026-09-04 | [review-findings-fix-plan](2026-09-04-review-findings-fix-plan.md) | Ordered fixes (Breslow ties, orientation, pooled β, CV log-likelihood) after the 9/4 chapter review | complete |
| 2026-10-01 | [prelim-proposal-plan](2026-10-01-prelim-proposal-plan.md) | Plan for the prelim Chapter 4: sync to bios-dissertation, content constraints, outline, figures | active (prelim due 11/16) |
| 2026-10-15 | [cox-logspace-arithmetic-plan](2026-10-15-cox-logspace-arithmetic-plan.md) | Mathematical note on numerically stable Cox risk-set sums (`cox_arithmetic = "logspace"`) | complete (implemented opt-in) |
| 2026-10-15 | [round2-profiled-curvature-proposal](2026-10-15-round2-profiled-curvature-proposal.md) | Opt-in profiled Cox intercept for the working curvature | complete (evaluated; not adopted) |
| 2026-10-16 | [round2-handoff-plan](2026-10-16-round2-handoff-plan.md) | Specification for the 10/15–16 diagnostic block: criteria, remedy order, per-task rules | complete (merged 2026-10-26) |
| 2026-10-16 | [round2-return-note](2026-10-16-round2-return-note.md) | Status of the 10/15–16 diagnostic work: results, task and commit index, open questions | complete; the reference summary for that work |
| 2026-10-16 | [round2-attempts-log](2026-10-16-round2-attempts-log.md) | Append-only record of each diagnostic study: design registered before it runs, then its outcome | historical work closed at merge; new October 30 supervision attempts appended |
| 2026-10-16 | [round2-fixed-q-loading-prior-proposal](2026-10-16-round2-fixed-q-loading-prior-proposal.md) | Opt-in fixed-posterior loading prior with consistent KL | complete (one case run; not adopted) |
| 2026-10-16 | [round2-closeout-and-docs-plan](2026-10-16-round2-closeout-and-docs-plan.md) | Freeze `round2/supervised` under a tag without merging, write one findings report, port vetted code, restructure the canonical documents | partly done: tagged `archive/round2-supervised-2026-10` and merged into main 2026-10-26 (Task 3 replaced by the merge); findings report folded into the 10/26 progress-book chapter; Tasks 4–5 after the meeting |
| 2026-10-16 | [preprocessing-alignment-and-unified-model-plan](2026-10-16-preprocessing-alignment-and-unified-model-plan.md) | One documented preprocessing pipeline and one model implementation for single-modality and multimodal results | not started (Task A1 becomes rebuild Task 13; A2–A4 replaced by the rebuild plan) |
| 2026-10-16 | [feature-selection-experiment-plan](2026-10-16-feature-selection-experiment-plan.md) | Unsupervised versus survival-screened gene selection at matched gene counts | training CV completed October 26, external selection deferred |
| 2026-10-30 | [Unified-Architecture-Alignment](2026-10-30-Unified-Architecture-Alignment.md) | Gemini diagnostic and solution proposal: consolidated diagnostics and a proposed unified architecture for the single-modality and multimodal models | proposal (unchecked external advice, not decisions) |
| 2026-10-30 | [realignment-plan](2026-10-30-realignment-plan.md) | Cleanup of open review items, then survival supervision of F under the YFB specification, the 10/30 chapter and the issue tracker | executed; corrected supervision gate failed; independent review pending; supersedes the 10/26 parallel-work plan |
| 2026-10-30 | [yfb-targeted-diagnostics-plan](2026-10-30-yfb-targeted-diagnostics-plan.md) | Matched noise-feedback, guarded-update and coefficient-uncertainty screens | complete; twelve fits and coefficient screen; no scientific adoption |
| 2026-11-02 | [package-and-paper-pathway-plan](2026-11-02-package-and-paper-pathway-plan.md) | When and how work moves from this repository into the SSBMF package and the SSBMF-paper manuscript: roles, triggers per stage, what moves and what stays | active |
| 2026-11-01 | [gemini-multiomicsGEP-diagnostic-and-tdd-guide](2026-11-01-gemini-multiomicsGEP-diagnostic-and-tdd-guide.md) | Gemini diagnostic guide for the YFB supervision and convergence blockers: fractional genomics likelihood with matched prior scaling, anchoring the residual precisions, a consistent normalized objective, too many survival-active factors, and a proposed test gate suite | proposal (unchecked external advice, not decisions) |
| 2026-11-03 | [survival-weight-p2000-plan](2026-11-03-survival-weight-p2000-plan.md) | Attempt S-9: the survival-weighted YFB fit on the S-7 grid at p = 50 and 2,000 over a range of w_S, with Gate S unchanged | complete; Gate S passes at w_S = 1 |
| 2026-11-04 | [supervised-real-data-plan](2026-11-04-supervised-real-data-plan.md) | Attempt S-11: the w_S = 1 fit on the real training data (S3), feature-set-size curve, survival-screening comparison and null check | complete; see findings |
| 2026-11-05 | [supervised-external-scoring-plan](2026-11-05-supervised-external-scoring-plan.md) | Attempt S-12: one external scoring of the fixed w_S = 1 fit and the two-step comparison | complete; see findings |
| 2026-11-05 | [supervision-strength-plan](2026-11-05-supervision-strength-plan.md) | Attempt S-13 design: survival-block-weighted objective (untempered β), weight chosen by CV partial log-likelihood, weak-signal test bed | complete 2026-11-05; see findings |
| 2026-11-05 | [near-threshold-plan](2026-11-05-near-threshold-plan.md) | Attempt S-14: supervised fit against the matched two-step fit for planted programs near the detection threshold, strong and realistic survival | complete 2026-11-05; see findings |
| 2026-11-06 | [supervised-real-refit-plan](2026-11-06-supervised-real-refit-plan.md) | Attempt S-15: real-data A1 fit with w_S chosen by five-fold CV partial log-likelihood, nested-CV estimate, adoption rule and conditional external scoring | complete 2026-11-06; w* = 1, not adopted |
| 2026-11-06 | [kl-fix-pruning-rerun-plan](2026-11-06-kl-fix-pruning-rerun-plan.md) | Attempt KL-1 (issue #6 item 1): October 15 single-modality K_init = 7 cells rerun after the KL fix, with pruning off, ordinary and protected | complete 2026-11-06; see findings |
| 2026-11-06 | [Flat β pruning rerun](2026-11-06-flat-beta-pruning-rerun-plan.md) | Repeat KL-1 after fixing the flat β likelihood | complete |
| 2026-11-05 | [chapter4-revamp-plan](2026-11-05-chapter4-revamp-plan.md) | Revamp prelim Chapter 4 around the supervised w_S = 1 YFB fit: outline, results checklist (tables, figures, biology), affected documents | complete 2026-11-06: chapter assembled and revised after review (`be75ade`, `502cd93`, `0560150`); Andrew's review of the assembled prelim pending (due 11/16) |
| 2026-11-05 | [chapter4-revamp-prompt](2026-11-05-chapter4-revamp-prompt.md) | Session prompt for the Chapter 4 revamp | used 2026-11-06 |
| 2026-11-09 | [week-plan](2026-11-09-week-plan.md) | Week of 11/9: L/F prior (point-exponential vs point-Laplace), β prior, methylation alone, M-modality fitter, LB vs YFB with fold-in; training-CV redo of choices that used external cohorts; results-rebuild design; inventory of open plans | active; nothing run yet |

## Superseded and historical plans

| Date | Document | Purpose | Status |
|---|---|---|---|
| 2026-04-16 | [ssmf-desurv-sim-benchmark-plan](2026-04-16-ssmf-desurv-sim-benchmark-plan.md) | Planned SSBMF-versus-DeSurv benchmark: weighted ELBO, ARD over K, simulation then PDAC | superseded by the `results/benchmark_sim/` pipeline; DeSurv head-to-head deferred |
| 2026-04-27 | [improvement-plan](2026-04-27-improvement-plan.md) | Options for fixing all-zero β on merged TCGA+CPTAC: harmonized feature selection, batch correction, prior changes | superseded by per-platform z-standardization and YFB (DECISIONS 2026-05-06) |
| 2026-05-08 | [yfb-fit-improvements-plan](2026-05-08-yfb-fit-improvements-plan.md) | Five-phase plan to make YFB risk orientation, preprocessing, K cross-validation and instrumentation consistent | superseded by the 2026-05-25/27 preprocessing plans and later K-selection work |
| 2026-07-08 | [ssbmf-post-lab-meeting-action-plan](2026-07-08-ssbmf-post-lab-meeting-action-plan.md) | Six-phase triage of the 6/18 lab feedback (objective normalization, joint vs two-step, K parsimony, stratified hazard, enrichment, manuscript) | superseded by the 2026-07-14 consolidation plan |
| 2026-07-15 | [ssbmf-open-items-followup-plan](2026-07-15-ssbmf-open-items-followup-plan.md) | Prioritized backlog from the 7/15 report (concordance CIs, factor stability, point-normal collapse, DeSurv, HPC) | superseded (only item 1a done; later plans re-prioritized the rest) |
| 2026-10-01 | [working-plan](2026-10-01-working-plan.md) | Working plan for multimodal YFB, real data and open single-modality items | superseded by [the 10/2 meeting follow-up](../meetings/2026-10-02-meeting-notes-followup.md) |
| 2026-10-15 | [advisor-priorities-and-execution-plan](2026-10-15-advisor-priorities-and-execution-plan.md) | Priorities and goal hierarchy after the 10/15 regrouping (progress report, CpG inventory, feature-count curves) | superseded by the 2026-10-16 plans and the 2026-10-26 parallel-work plan |
| 2026-10-15 | [m6-m7-handoff-plan](2026-10-15-m6-m7-handoff-plan.md) | Handoff for the robustness and comparator milestones | superseded |
| 2026-10-15 | [m6-m7-return-note](2026-10-15-m6-m7-return-note.md) | Results returned from the robustness and comparator milestones | superseded (summary in [m6-m7-summary](../findings/2026-10-15-m6-m7-summary.md)) |
| 2026-10-15 | [round2-history-note](2026-10-15-round2-history-note.md) | Long status log of the 10/15 diagnostic work, kept verbatim | superseded by the round2 return note |
| 2026-10-15 | [comparator-inventory](2026-10-15-comparator-inventory.md) | Inventory of comparator methods and their availability | superseded |
| 2026-10-22 | [unified-method-rebuild-plan](2026-10-22-unified-method-rebuild-plan.md) | Gated rebuild of the method from a written specification as the `ssbmf` R package | historical reference |
| 2026-10-26 | [parallel-work-sessions-plan](2026-10-26-parallel-work-sessions-plan.md) | Coordinates work on the current method (session A) with the rebuild (session B); round2 merge and the 10/26 meeting deliverables | superseded by the October 30 realignment |
| 2026-11-06 | [real-data-comparators-plan](2026-11-06-real-data-comparators-plan.md) | Fixed S-11/S-12 references and four training-only prognostic comparators | registered |
| 2026-11-06 | [desurv-style-analyses-plan](2026-11-06-desurv-style-analyses-plan.md) | Fixed risk groups, program classification, rank comparisons and mixed-signal simulation | registered |
