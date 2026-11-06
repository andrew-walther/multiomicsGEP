# Real-data comparators

Registration for [Prompt 6](../prompts/2026-11-06-codex-handoff-prompts.md).
Commit this registration with its index row and attempts-log entry before any
fit or external evaluation. Run Prompt 7 after the Prompt 6 results commit.

1. **Reconstruct the reference scores.** Use the S-11 D4 full-training w1 and
   w0 checkpoints and the saved TCGA+CPTAC training input. Do not refit either
   factorization. Reproduce the S-12 Cox regression for w0 and the S-12
   projection and preprocessing for both arms. Verify the two scores
   reproduce S-12 aggregate results on the full eligible samples.
2. **Fit comparators in TCGA+CPTAC.** Fit a univariable Cox model on the
   available PurIST probabilities. Audit stage, grade, age and margin
   mappings without inspecting external survival outcomes. Use harmonizable
   fields available in both training cohorts and all five external cohorts,
   with complete-case training and scoring. Report unavailable fields,
   coverage and reasons for exclusions. If no common field can be harmonized,
   record the clinical arm as unavailable. Fit D4 elastic-net Cox using the
   S-11 training folds, with alpha in {0, .25, .5, .75, 1} and lambda chosen by
   minimum training CV partial-likelihood deviance. Use the larger penalty
   for ties within alpha and the larger alpha for ties across alpha. Verify
   fold alignment and record the selected settings and denominators.
3. **Fit actual DeSurv.** Reuse the pinned author implementation,
   preprocessing adapter and training-only tuning workflow in
   `results/benchmark_sim/run_desurv_h2h.R`. Use
   `run_desurv_comparison.R` and `code/preprocess_desurv.R` for the SSBMF
   preprocessing context. Verify training patient identity and order before
   tuning and report any incompatibility. Freeze genes, model settings and
   score definitions before accessing external outcomes.
4. **Evaluate the five S-12 external cohorts.** Use Cox-oriented risks and
   `frozen_reverse_cindex()` for all concordance calculations. Use 2,000
   percentile bootstrap replicates with seed 1. Report per-cohort C and
   calibration slopes on the original score scale. Estimate paired
   supervised-minus-comparator C on each comparator's common patients.
   Average the within-cohort paired differences equally, resampling linked
   donors across assays via `code/concordance_ci.R`. Report omitted strata.
   Estimate pooled HR per pooled external SD using cohort-stratified Cox.
   Verify the donor linkage and report the eligible, scored, paired and
   bootstrap denominators and any warnings or undefined estimates.
5. **Publish aggregate outputs.** Add the runner, aggregate CSVs, findings
   and findings index row. Include exact regeneration commands and a
   source table for every reported result. Run regression checks and the
   project test suite before the results commit. Keep patient-level inputs
   and scores local and excluded from Git. Leave the prelim, writeup and
   progress book to the reviewer.
