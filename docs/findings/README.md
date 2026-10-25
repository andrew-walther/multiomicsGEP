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
