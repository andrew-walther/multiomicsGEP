# Findings

Results, reports and data audits from work that ran. Most are Quarto reports
(`.qmd` with rendered `.html`/`.pdf`); the link goes to the source. The plan
behind each one is in [`../plans/`](../plans/README.md). Shorter advisor-facing
status updates are in [`../progress_report/`](../progress_report/README.md).

| Date | Document | Purpose | Status |
|---|---|---|---|
| 2026-03-31 | [pdac-data-audit/](2026-03-31-pdac-data-audit/PDAC_data_audit.qmd) | The seven PDAC expression cohorts (platform, survival, censoring, file layout) and whether the loader handles them | complete (corrected 2026-10-16: CPTAC is RNA-seq) |
| 2026-04-29 | [ssbmf-summary-report](2026-04-29-ssbmf-summary-report.qmd) | Simulation validation and DeSurv-matched PDAC benchmark of the original LB model | superseded by 2026-05-05 |
| 2026-05-05 | [ssbmf-summary-report](2026-05-05-ssbmf-summary-report.qmd) | LB versus YFB re-benchmark after the first round of fixes | superseded by 2026-05-06 |
| 2026-05-06 | [ssbmf-summary-report](2026-05-06-ssbmf-summary-report.qmd) | Synthesis of four fix rounds (per-platform standardization, multistart, K-CV, LB vs YFB) with external validation | superseded by the 2026-05-25 and 2026-05-27 reports |
| 2026-05-22 | [cohort-lmm-benchmark-report](2026-05-22-cohort-lmm-benchmark-report.qmd) | Cohort-indicator extension on synthetic and real PDAC data | complete (indicator not adopted) |
| 2026-05-25 | [merged-benchmark-report](2026-05-25-merged-benchmark-report.qmd) | 18-configuration merged-cohort benchmark of preprocessing × model × cohort indicator | superseded (its M5/K = 3 recommendation was replaced by D4/K = 7 on 2026-05-27) |
| 2026-05-27 | [desurv-alignment-report](2026-05-27-desurv-alignment-report.qmd) | How DeSurv-style per-cohort gene selection changes external C for the candidate configurations | complete (its 0.636 later corrected to 0.627) |
| 2026-05-27 | [desurv-factor-diagnostics-report](2026-05-27-desurv-factor-diagnostics-report.qmd) | Convergence, survival associations, Kaplan–Meier curves and top genes per configuration | superseded in part: its adverse/protective labels were reversed by the 2026-06-16 direction correction |
| 2026-07-12 | [joint-vs-twostep-sweep-report](2026-07-12-joint-vs-twostep-sweep-report.qmd) | Simulation sweep over survival effect size: joint YFB/LB versus PCA→Cox and EBMF→Cox | complete |
| 2026-07-15 | [pathway-enrichment-report](2026-07-15-pathway-enrichment-report.qmd) | Gene-set enrichment, subtype concordance and DeSurv overlap for the two survival-active programs | superseded by 2026-08-21 |
| 2026-07-15 | [ssbmf-progress-report](2026-07-15-ssbmf-progress-report.qmd) | Long progress report since 6/18: objective fixes, simulations, K parsimony, external validation, biology | complete |
| 2026-08-21 | [pathway-enrichment-report](2026-08-21-pathway-enrichment-report.qmd) | Enrichment for all four retained K = 7 programs, including the genomics-only programs | complete (current reference) |
| 2026-08-21 | [session-progress-findings](2026-08-21-session-progress-findings.md) | Phase checklist of the 8/21 meeting-prep work (model chapter, factor classification, K recovery, EBMF comparisons) with results | complete for Phases 0–4; Phase 5 never scoped |
| 2026-10-15 | [m6-m7-summary](2026-10-15-m6-m7-summary.md) | Robustness (jittered starts, Normal loading priors, feature counts) and comparator results, with branch commits | complete (no configuration adopted) |
| 2026-10-16 | [advisor-progress-since-10-02](2026-10-16-advisor-progress-since-10-02.qmd) | Five-page advisor report of the evidence since 10/2 and the revised execution order (`.qmd`, `.html`, `.pdf`) | complete |
| 2026-10-16 | [advisor-meeting-summary](2026-10-16-advisor-meeting-summary.qmd) | Two-page brief for the 10/16 meeting: model, validation, biology, manuscript (`.qmd`, `.html`, `.pdf`) | complete |
| 2026-10-26 | [single-log-check](2026-10-26-single-log-check.md) | Effect of the double log₂(x + 1) on three RNA-seq cohorts: refit of the benchmark on single-logged data, training CV and one external scoring | complete (no material change; recommendation for Andrew) |
| 2026-10-26 | [survival-screen-cv](2026-10-26-survival-screen-cv.md) | Unsupervised versus survival-screened gene sets at 251/500/1,000 genes: training CV and program agreement with the benchmark fit | complete for training CV; external scoring waits for Andrew's review |
| 2026-10-26 | [preprocessing-audit](2026-10-26-preprocessing-audit.md) | Step-by-step audit of single-modality and multimodal preprocessing: the double log, twelve discrepancies, the sign floor of nonnegative factors on z-scored data, proposed alignments | complete; alignments await Andrew's decision |
| 2026-10-27 | [operational-roadblocks](2026-10-27-operational-roadblocks.md) | Ranked roadblocks to an operational supervised model, with evidence, rebuild status and remedies; what is resolved; the 10/17 outside review's claims checked | complete |
| 2026-10-30 | [current-method-a1](2026-10-30-current-method-a1.md) | Whether a larger survival weight w_S lets survival guide the loadings F in the current single-modality fitter: simulation grid and TCGA + CPTAC 5-fold CV | complete (no: w_S = 1 is best; survival's share reaches O(1) only where the fit fails); merged 2026-10-30, Gate A1 recommendation (drop w_S > 1) awaits Andrew |
| 2026-10-30 | [yfb-supervision-diagnosis](2026-10-30-yfb-supervision-diagnosis.md) | Precision balance, DeSurv loss audit and registered numerical diagnoses for normalized YF supervision | complete; corrected shared-loading gate failed |
| 2026-10-30 | [yfb-supervision-problem-statement](2026-10-30-yfb-supervision-problem-statement.md) | Standalone statistical failure, numerical repairs, exact gate evidence and next model decisions | unresolved method question; no external selection |

| 2026-10-30 | [yfb-targeted-diagnostics-findings](2026-10-30-yfb-targeted-diagnostics-findings.md) | Twelve matched noise/guarded-update fits and exact versus mean-field Gaussian coefficient diagnostics | complete; partial stabilization; no screening success |
