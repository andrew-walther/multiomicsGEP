# Attempt S-12: one external scoring of the fixed survival-weighted YFB fit

Registered 2026-11-05, before any external cohort was loaded for this model. Andrew reviewed the
S-11 training results and asked for this scoring the same day. It is the external step named in
the S-11 registration ([plan](2026-11-04-supervised-real-data-plan.md)) and in DECISIONS.md
(2026-11-04): every choice is already fixed by training cross-validation, and nothing here is
chosen from the external results.

## Fixed models (no refit)

- **Supervised fit:** the S-11 all-patient fit `fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1)` on
  the 2,218 D4 genes, saved as `fits/size_D4_N2218_w1_f0.rds` in the S-11 output folder. Risk
  score: `predict_ssbmf_yfb(Y_ext, E[F], E[β], EF_norms)`.
- **Two-step comparison:** the S-11 all-patient w_S = 0 fit on the same genes
  (`fits/size_D4_N2218_w0_f0.rds`), followed by a Cox regression of training survival on the
  live program scores Z = Y E[F], exactly as the `w0cox` arm of the S-11 summarizer. Risk score:
  Y_ext E[F] b.
- The historical frozen-YFB external C values (`results/benchmark_snapshot_2026-10-14/`) are
  shown beside these, labelled as the earlier pipeline (double log, different gene handling).
  They are not a paired comparison.

## External data and preprocessing

The five cohorts in `cfg$pdac$external_cohorts`: Dijk, Moffitt_GEO_array, PACA_AU_array,
PACA_AU_seq, Puleo_array. Each is loaded with `load_pdac_raw()` under the current single-log
flags (all false), restricted to the D4 genes, and z-scored per gene within the cohort, the same
transform the training cohorts received. A D4 gene absent from a cohort is set to 0 (its cohort
mean after z-scoring), which equals dropping that gene's row of E[F]. The number of missing genes
per cohort is reported. No external outcome enters preprocessing.

## Reported

Orientation is fixed by the fits (positive risk means higher hazard). No score is flipped.

1. **Primary:** C per cohort for both models, with a 95% percentile bootstrap interval
   (`bootstrap_concordance_ci`, B = 2,000, seed 1, `flip = FALSE`).
2. Paired difference supervised minus two-step per cohort (`bootstrap_concordance_diff_ci`,
   `flip_a = flip_b = FALSE`), and the cohort-stratified pooled difference
   (`bootstrap_concordance_diff_ci_stratified`).
3. Calibration slope per cohort: the coefficient of a Cox regression of survival on the risk
   score. A slope of 1 means the fitted risk differences have the right size. This checks the
   S-9 calibration limitation on real data.
4. Per-cohort n, events and missing D4 genes.

PACA_AU_array and PACA_AU_seq contain overlapping donors. They are reported separately. The
pooled difference is also reported without PACA_AU_seq, as a sensitivity check.

## Software

`results/benchmark_sim/score_yfb_supervised_external.R`, writing aggregate CSVs to
`results/benchmark_sim/outputs/yfb_supervised_external_2026-11-05/`. No patient-level output is
written to the repository. No fitter or update code changes.
