# docs/

Project documents, sorted by what kind of document they are. Each type folder
has a README that lists its documents with a one-line purpose and a status.

| Folder | What goes there |
|---|---|
| [`plans/`](plans/README.md) | Plans, design specs and proposals, written before the work |
| [`findings/`](findings/README.md) | Results, reports and data audits from work that ran |
| [`reviews/`](reviews/README.md) | Reports from independent reviewers, and external advice as received |
| [`prompts/`](prompts/README.md) | Paste-in prompts, briefs and session handoffs |
| [`meetings/`](meetings/README.md) | Meeting notes, follow-ups and summaries |
| [`notes/`](notes/README.md) | Derivations and technical notes |
| [`progress_book/`](https://andrew-walther.github.io/multiomicsGEP/) | Progress notebook, one chapter per advisor meeting (Quarto book published to GitHub Pages; chapter names are the meeting date) |
| [`progress_report/`](progress_report/README.md) | Earlier advisor status updates (April–August 2026), before the notebook |
| `lb_walkthroughs/` | Walkthroughs of each LB update (`update_beta`, `update_L`, `update_F`, `update_tau`) and the LB fitting loop (`fit_modular`), plus [`update_L_fix.md`](lb_walkthroughs/update_L_fix.md), the debugging guide for the L update. Undated: each documents a code file |

`make -C docs` re-renders the walkthroughs and the data audit
(`findings/2026-03-31-pdac-data-audit/`).

## Naming

`YYYY-MM-DD-topic-type.ext`, all hyphens. The date is the day the document was
written: from its old file name, or failing that its first commit. The type word
matches the folder (`-plan`, `-spec`, `-report`, `-findings`, `-review`,
`-prompt`, `-notes`); meeting records end in the meeting's own name. A Quarto
document keeps its rendered `.html`/`.pdf` under the same stem. A group of
documents that belong together sits in a dated subfolder.

When you add a document, give it a row in its folder's README.

## Status words

- **active**: still in use, or still guides open work.
- **complete**: the work is done. The document is the record.
- **superseded**: a later document replaces it (named in the row).
- **not started**: approved or drafted, waiting on something.

## Renamed on 2026-10-25

Before 2026-10-25 documents sat in `plans/` (current) and `plans/archive/`
(superseded plans, meeting follow-ups and design specs, moved there on
2026-10-16 from `docs/plans/` and `docs/superpowers/`), `reports/`, `notes/` and
the top of `docs/`, with dates in several formats at the end of the name. The
`archive/` folder is gone: superseded documents stay in their type folder with
status **superseded**.

A few files still cite the old paths and were deliberately left alone:

- `paper/prelim/project3-ssbmf.qmd` (a source comment). Committing it fires the
  post-commit hook that syncs Chapter 4 into bios-dissertation.
- `paper/writeup/_process/`: drafting and review records, including
  `source_inventory.csv`, that point at evidence pinned to a research commit.
- Generated output: the progress book's `_book/` and the rendered `.html` of
  each report. They update when re-rendered.
- Files in the bios-dissertation repository (`ROADMAP.md`,
  `prelim/references.bib`, literature-review notes) that cite
  `docs/reports/cohort_lmm_review/references.bib` and two plans.
- Paths on the `round2/supervised` branch, which was not merged.

Use the table to find the new name for an old path. A row for a `.qmd` covers its
`.html`/`.pdf`.

| Old path | New path |
|---|---|
| `docs/plans/2026-10-16_feature_selection_experiment_plan.md` | `docs/plans/2026-10-16-feature-selection-experiment-plan.md` |
| `docs/plans/2026-10-16_preprocessing_alignment_and_unified_model_plan.md` | `docs/plans/2026-10-16-preprocessing-alignment-and-unified-model-plan.md` |
| `docs/plans/2026-10-16_round2_closeout_and_docs_plan.md` | `docs/plans/2026-10-16-round2-closeout-and-docs-plan.md` |
| `docs/plans/2026-10-22_unified_method_rebuild_plan.md` | `docs/plans/2026-10-22-unified-method-rebuild-plan.md` |
| `docs/plans/Meeting_Notes_Followup_10_2_26.md` | `docs/meetings/2026-10-02-meeting-notes-followup.md` |
| `docs/plans/Prelim_Proposal_Plan_10_1_26.md` | `docs/plans/2026-10-01-prelim-proposal-plan.md` |
| `docs/plans/archive/2026-03-25-fit-modular-design.md` | `docs/plans/2026-03-25-fit-modular-design-spec.md` |
| `docs/plans/archive/2026-03-25-fit-modular.md` | `docs/plans/2026-03-25-fit-modular-plan.md` |
| `docs/plans/archive/2026-05-25-merged-cohort-benchmark.md` | `docs/plans/2026-05-25-merged-cohort-benchmark-plan.md` |
| `docs/plans/archive/2026-05-25-preprocessing-comparison-k-floor.md` | `docs/plans/2026-05-25-preprocessing-comparison-k-floor-plan.md` |
| `docs/plans/archive/2026-05-27-desurv-aligned-preprocessing.md` | `docs/plans/2026-05-27-desurv-aligned-preprocessing-plan.md` |
| `docs/plans/archive/Improvement_Plan_4_27_26.md` | `docs/plans/2026-04-27-improvement-plan.md` |
| `docs/plans/archive/Meeting_Notes_Followup_4_27_26.md` | `docs/meetings/2026-04-27-meeting-notes-followup.md` |
| `docs/plans/archive/Meeting_Notes_Followup_8_21_26.md` | `docs/meetings/2026-08-21-meeting-notes-followup.md` |
| `docs/plans/archive/Meeting_Notes_Followup_8_27_26.md` | `docs/meetings/2026-08-27-meeting-notes-followup.md` |
| `docs/plans/archive/Meeting_Notes_Followup_8_3_26.md` | `docs/meetings/2026-08-03-meeting-notes-followup.md` |
| `docs/plans/archive/SSMF_DeSurv_Sim_Benchmark_4_16_26.md` | `docs/plans/2026-04-16-ssmf-desurv-sim-benchmark-plan.md` |
| `docs/plans/archive/SSMF_Meeting_Synthesis_4_30_26.md` | `docs/meetings/2026-04-30-ssmf-meeting-synthesis.md` |
| `docs/plans/archive/Working_Plan_10_1_26.md` | `docs/plans/2026-10-01-working-plan.md` |
| `docs/plans/archive/YFB-fit-improvements-05_08_26.md` | `docs/plans/2026-05-08-yfb-fit-improvements-plan.md` |
| `docs/plans/archive/beta_zero_fix_design.md` | `docs/plans/2026-04-29-beta-zero-fix-design-spec.md` |
| `docs/plans/archive/cross_modality_diagnostics_prompt_4_28_26.md` | `docs/prompts/2026-04-28-cross-modality-diagnostics-prompt.md` |
| `docs/plans/archive/joint_k_alpha_bayesopt_plan_07_12_2026.md` | `docs/plans/2026-07-12-joint-k-alpha-bayesopt-plan.md` |
| `docs/plans/archive/meeting_2026-05-29_multi_cohort_simulation.md` | `docs/meetings/2026-05-29-multi-cohort-simulation-meeting.md` |
| `docs/plans/archive/meeting_2026_08_21_session_progress.md` | `docs/findings/2026-08-21-session-progress-findings.md` |
| `docs/plans/archive/pathway_enrichment_plan.md` | `docs/plans/2026-07-14-pathway-enrichment-plan.md` |
| `docs/plans/archive/pathway_enrichment_prompt.md` | `docs/prompts/2026-06-16-pathway-enrichment-prompt.md` |
| `docs/plans/archive/rashid_lab_meeting_notes_06_18_2026.md` | `docs/meetings/2026-06-18-rashid-lab-meeting-notes.md` |
| `docs/plans/archive/rashidlab_meeting_feedback_04_09_26.md` | `docs/meetings/2026-04-09-rashidlab-meeting-feedback.md` |
| `docs/plans/archive/review_findings_fix_plan_09_04_2026.md` | `docs/plans/2026-09-04-review-findings-fix-plan.md` |
| `docs/plans/archive/ssbmf_factor_classification_k_selection_08_13_2026.md` | `docs/plans/2026-08-13-ssbmf-factor-classification-k-selection-plan.md` |
| `docs/plans/archive/ssbmf_item2_executive_summary_plan_07_15_2026.md` | `docs/plans/2026-07-15-ssbmf-item2-executive-summary-plan.md` |
| `docs/plans/archive/ssbmf_k_parsimony_followup_plan_07_13_2026.md` | `docs/plans/2026-07-13-ssbmf-k-parsimony-followup-plan.md` |
| `docs/plans/archive/ssbmf_open_items_followup_plan_07_15_2026.md` | `docs/plans/2026-07-15-ssbmf-open-items-followup-plan.md` |
| `docs/plans/archive/ssbmf_post_lab_meeting_action_plan_07_08_2026.md` | `docs/plans/2026-07-08-ssbmf-post-lab-meeting-action-plan.md` |
| `docs/plans/archive/ssbmf_progress_consolidation_and_remaining_work_07_14_2026.md` | `docs/plans/2026-07-14-ssbmf-progress-consolidation-and-remaining-work-plan.md` |
| `docs/plans/archive/yfb_dual_source_F_experiments_08_20_2026.md` | `docs/plans/2026-08-20-yfb-dual-source-f-experiments-plan.md` |
| `docs/notes/YFB_derivation_05_08_26.qmd` | `docs/notes/2026-05-08-yfb-derivation-notes.qmd` |
| `docs/notes/YFB_mathematical_consistency_audit_05_08_26.qmd` | `docs/notes/2026-05-08-yfb-mathematical-consistency-audit-notes.qmd` |
| `docs/notes/YFB_supplemental_proofs_05_08_26.qmd` | `docs/notes/2026-05-08-yfb-supplemental-proofs-notes.qmd` |
| `docs/notes/YFB_update_derivations_detailed_05_08_26.qmd` | `docs/notes/2026-05-08-yfb-update-derivations-detailed-notes.qmd` |
| `docs/notes/sign_orientation_issue_05_08_26.qmd` | `docs/notes/2026-05-08-sign-orientation-issue-notes.qmd` |
| `docs/notes/gemini/Gemini_review_10_17_26.md` | `docs/reviews/2026-10-17-gemini-review.md` |
| `docs/reports/pathway_enrichment_overview.qmd` | `docs/notes/2026-05-30-pathway-enrichment-overview-notes.qmd` |
| `docs/reports/cohort_lmm_benchmark_report_5_22_26.qmd` | `docs/findings/2026-05-22-cohort-lmm-benchmark-report.qmd` |
| `docs/reports/desurv_alignment_report_05_27_26.qmd` | `docs/findings/2026-05-27-desurv-alignment-report.qmd` |
| `docs/reports/desurv_factor_diagnostics_05_27_26.qmd` | `docs/findings/2026-05-27-desurv-factor-diagnostics-report.qmd` |
| `docs/reports/joint_vs_twostep_sweep_07_12_2026.qmd` | `docs/findings/2026-07-12-joint-vs-twostep-sweep-report.qmd` |
| `docs/reports/merged_benchmark_report_5_25_26.qmd` | `docs/findings/2026-05-25-merged-benchmark-report.qmd` |
| `docs/reports/multicohort_sim_proposal_06_14_26.qmd` | `docs/plans/2026-06-14-multicohort-sim-proposal.qmd` |
| `docs/reports/pathway_enrichment_report_07_15_26.qmd` | `docs/findings/2026-07-15-pathway-enrichment-report.qmd` |
| `docs/reports/pathway_enrichment_report_08_21_26.qmd` | `docs/findings/2026-08-21-pathway-enrichment-report.qmd` |
| `docs/reports/ssbmf_progress_report_07_15_26.qmd` | `docs/findings/2026-07-15-ssbmf-progress-report.qmd` |
| `docs/reports/ssbmf_summary_report_04_29_26.qmd` | `docs/findings/2026-04-29-ssbmf-summary-report.qmd` |
| `docs/reports/ssbmf_summary_report_05_05_26.qmd` | `docs/findings/2026-05-05-ssbmf-summary-report.qmd` |
| `docs/reports/ssbmf_summary_report_05_06_26.qmd` | `docs/findings/2026-05-06-ssbmf-summary-report.qmd` |
| `docs/reports/cohort_lmm_review/` | `docs/plans/2026-05-20-cohort-lmm-review/` |
| `docs/progress_report/SSBMF_Status_Update_4_29_26.qmd` | `docs/progress_report/2026-04-29-ssbmf-status-update.qmd` |
| `docs/progress_report/SSBMF_Status_Update_5_07_26.qmd` | `docs/progress_report/2026-05-07-ssbmf-status-update.qmd` |
| `docs/progress_report/SSBMF_Status_Update_5_27_26.qmd` | `docs/progress_report/2026-05-27-ssbmf-status-update.qmd` |
| `docs/progress_report/SSBMF_Status_Update_5_28_26.qmd` | `docs/progress_report/2026-05-28-ssbmf-status-update.qmd` |
| `docs/progress_report/SSBMF_Status_Update_07_15_26.qmd` | `docs/progress_report/2026-07-15-ssbmf-status-update.qmd` |
| `docs/progress_report/SSBMF_Status_Update_08_03_26.qmd` | `docs/progress_report/2026-08-03-ssbmf-status-update.qmd` |
| `docs/reviews/2026-09-04_progress_notebook_review.md` | `docs/reviews/2026-09-04-progress-notebook-review.md` |
| `docs/reviews/2026-10-02_chapter_review.md` | `docs/reviews/2026-10-02-chapter-review.md` |
| `docs/reviews/2026-10-02_prelim_proposal_review.md` | `docs/reviews/2026-10-02-prelim-proposal-review.md` |
| `docs/PDAC_data_audit.qmd` | `docs/findings/2026-03-31-pdac-data-audit/PDAC_data_audit.qmd` |
