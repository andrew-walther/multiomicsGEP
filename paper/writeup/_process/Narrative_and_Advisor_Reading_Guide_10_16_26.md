# Narrative and advisor discussion guide

Current evidence through October16,2026. This note supports the next advisor
meeting; the scientific results remain those in the write-up and frozen tables.

## Narrative to develop

Lead with the scientific question: bulk tumors contain tumor-state and
microenvironmental programs whose molecular prominence and prognostic roles
need not agree. Explain the model as a way to represent those roles, predict
risk in new cohorts and investigate matched modalities. Then provide evidence
in this order: actual benchmark prediction; identity and clinical associations
of programs; feature-range stability; what methylation contributes; unresolved
requirements for survival-guided discovery. A list of update/debugging efforts
is supporting history, not the manuscript's central argument.

The present Single benchmark supports cross-cohort prediction and interpretable
program structure, but alpha_F0 prevents attributing its programs to survival-
guided loading discovery. Bayesian shrinkage/uncertainty and shared multimodal
activities are methodological opportunities; calibrated uncertainty, valid
rank selection and added multimodal biology still need evidence. Projection
of new patients is useful but is also available in DeSurv; it is not unique.
Avoid a superiority claim against DeSurv before the deferred head-to-head.

## Useful lessons from the reference materials

Amber's paper connects a biological problem to one architectural choice and
then asks separately whether true programs are recovered, biological roles
make sense and findings transfer to held-out cohorts. Her defense deck makes
that separation visible before technical detail. The thesis chapter supplies
a fuller problem-to-method-to-evidence argument. We can use that organization
without copying wording, numerical claims, tumor-stroma interpretations or
claims of successful supervision into this project.

For our biology, emphasize what the evidence shows: adverse basal-like and
protective classical tumor-associated programs, plus stromal and immune
structure that need not be prognostic on its own. The joint fit additionally
contains methylation-dominant structure, whose cause and added value are
unvalidated. Biological usefulness need not equal a higher C-index, but needs
annotation, reproducibility and tests of composition/technical explanations.

## Advisor route and decisions

1. Read the new structured summary; agree which model is evaluated versus intended.
2. Review the model table and the projection/prediction rule, leaving update
   derivations for questions rather than reading them sequentially.
3. Review external prediction with the actual independently tuned comparator;
   keep C differences, uncertainty, calibration and program recovery distinct.
4. Review the four Single programs and joint roles; identify the biological
   question that a multimodal analysis could answer beyond prediction.
5. Agree the central claim, essential methodological gates and manuscript
   scope. Use the existing journal-shortlist note after agreeing that claim.

## Sources inspected, not new analyses

- DeSurv paper source: /tmp/round2_DeSurv_paper_0a870a7/paper/paper.Rmd,
  02_introduction_REVISED.Rmd, 04_results_REVISED.Rmd and 05_discussion_REVISED.Rmd;
  author repository https://github.com/rashidlab/DeSurv-paper, pinned local
  commit0a870a78673296a3802aeddfee0514e4a8b0f9ae.
- Package README/API: https://github.com/rashidlab/DeSurv,
  local commitcc9d4c2486b3125704ab71eea61f35e040a40080; no fit/install rerun.
- Amber thesis: bios-dissertation/templates/source/dissertation-resources-amber/
  dissertation 1.pdf, Chapter2 (printed pp22–42; PDF pp54–74), with the
  dissertation.Rmd wrapper. Referenced chapter source files are absent from
  that local copy, so the rendered chapter is authoritative for this review.
- Defense: defense.qmd in the same directory, DeSurv portion and speaker notes.
- Project claims: writeup.qmd, source_inventory.csv, frozen benchmark tables,
  and completed feature-count tables. The official related-method citation
  is retained as Young et al., submitted; no publication status is inferred
  from a repository README's journal label.

No private biological result is copied, no external outcome is used for
selection, and no DeSurv experiment is restarted by this narrative review.
