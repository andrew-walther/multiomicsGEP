# Attempt S-11: the survival-weighted YFB fit on the real training data

Registered 2026-11-04, before any fit below was run. This follows S-9, which passed Gate S at
w_S = 1 ([findings](../findings/2026-11-04-survival-weight-p2000-findings.md)), and the
2026-11-04 decision that real-data work follows a working simulation result. It covers step S3
of the [realignment plan](2026-10-30-realignment-plan.md), the feature-set-size curve and the
survival-screening comparison listed in the
[S-9 plan](2026-11-03-survival-weight-p2000-plan.md), and the null calibration check of the
work order's Part B step 3.2. The attempts log carries the matching S-11 entry.

## Fixed model

`fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1)`: point-Laplace L and F, per-program Normal β,
K = 7 (`cfg$benchmark$k_merged_yfb_desurv`), at most 300 sweeps, tolerance 10⁻⁵, α from
`config/globals.yml`, SVD start, no supervised start, no sign correction. This is the S-9 call.
The comparison arm is the same fitter at w_S = 0 (survival does not enter the loadings; β is
still fitted), labelled "unsupervised loadings". It is the same-pipeline counterpart of the
historical α_F = 0 results, which stay on record unchanged.

*Correction (2026-11-04, after the run):* at w_S = 0 the canonical fitter also drops the Cox
term from the β update, so β = 0 and the arm has no risk score (C undefined). The summarizer
therefore adds a two-step arm, `w0cox`, not in the original registration: Cox regression of
survival on the w_S = 0 training program scores Z = Y E[F], scored on the held-out fold. It is
the unsupervised-loadings comparison the registration intended.

## Data

- **Training patients:** TCGA_PAAD and CPTAC (`cfg$pdac$training_cohorts`). No external cohort
  is loaded by any step here.
- **Preprocessing:** the current single-log configuration (every `platform_log_transform` flag
  false; the loaders return log2-scale values once), then a z-score of every gene within each
  cohort, on all genes shared by the two cohorts. This is the centered input S-9 used.
- **Gene ranking for the size curve:** the existing combined mean/variance rank within each
  cohort (`feature_count_manifest()` in `code/feature_sensitivity.R`), which uses no outcome.
  A gene enters at the larger of its two cohort ranks, so the sets are nested. The smallest set
  is the D4 size (per-cohort top 3,000, `cfg$preprocessing$top_n_genes_desurv`).
- **Folds:** five folds stratified by cohort × event, `set.seed(2026)`, as in
  `run_a1_real_ws_cv.R`. Gene ranking and z-scoring precede the folds and use no outcome, so
  CV C is conditional on them. The Cox screen (below) is the only outcome-based step, and it runs
  inside the folds.
- The prepared matrix is saved under `data/` (git-ignored). No patient-level file is committed.

## Part 1. Step S3 and the feature-set-size curve

**Sizes:** the D4 set, then about 3,000, 5,000, 8,000, 12,000 shared genes and all shared genes
(actual counts reported; ties can move a target). **Arms:** w_S = 1 and w_S = 0. Per size and
arm: five fold fits and one fit on all training patients.

Per fit: held-out C with the orientation fixed at fit time (positive β means higher hazard; no
flipping), convergence and sweeps, live and survival-active programs (|E β_k| ≥ 1.96 SD), final
survival share of the loading precision, max |E β|, reconstruction loss and runtime. Per size and
arm: mean fold C ± SE, pooled out-of-fold C, CV partial log-likelihood, and folds converged.

**S3 deliverables** come from the D4 size at w_S = 1: the programs of the all-patient fit, which
are survival-active and with which sign, their top-loading genes, and the CV C, reported next
to w_S = 0.

**Program agreement across sizes:** for each size, the all-patient fit is matched to the
all-genes fit by patient-level program scores (`feature_program_agreement()`), and the
risk-score correlation is reported.

**Preferred size (registered rule).** The smallest size whose mean fold C is within one SE of
the all-genes set's mean fold C and whose all-patient fit has a match with |r| ≥ 0.7 for every
survival-active program of the all-genes fit. If the all-genes set is not the best size, the
one-SE rule against the best size is also reported. An elbow plot shows mean fold C ± SE against
size for both arms.

## Part 2. Survival screening before the joint fit

The October 16 design (`2026-10-16-feature-selection-experiment-plan.md`) with the fixed model:

- U(N): the nested unsupervised mean/variance screen inside the D4 set, N ∈ {250, 500, 1,000}
  (actual counts from the entry ranks).
- S(N): the D4 genes ranked by the Wald p-value of `coxph(Surv ~ y + strata(cohort))` on
  fold-training patients only, top N matched to U(N)'s actual count. The all-patient fit uses
  all training patients for the screen.
- D4 itself (from Part 1) is the reference.
- Arms: w_S = 1 and w_S = 0, so the screening effect is reported both with and without survival
  in the loadings, on the same pipeline.

Reported per N: mean fold C ± SE for U and S, the paired S(N) − U(N) difference over folds with
its SE, convergence, survival-active programs, and the overlap of S(N) with U(N). The question is
whether a survival screen adds anything once survival already shapes the loadings. The
October 26 α_F = 0 results (S above U by 0.038, −0.016 and 0.014) stay on record beside these,
labelled as the historical double-log, unsupervised-loadings model. No arm is adopted from this
part without Andrew's review.

## Part 3. Null calibration check (simulation)

The S-7 datasets (seeds 1–10, p ∈ {50, 2,000}, variance 0.25%, 1%, 5%) with the training
survival outcome permuted across patients (`set.seed(seed + 9000)`), so survival is independent
of expression. Fit at w_S = 1 on the z-scored input. Every survival-active program is false.
**Criterion (Part B step 3.2):** mean false-active at most 1 per p × variance cell. Convergence
is also reported.

## External cohorts

None here. After Andrew has reviewed Parts 1 and 2, and the size and screening arm are fixed
by training CV, the fixed model is scored once on the external cohorts (DECISIONS.md
2026-11-04). That scoring gets its own registered entry.

## Software and compute

`results/benchmark_sim/run_yfb_supervised_real.R` with modes `prepare`, `size`, `screen`,
`null` and `summary`. PSOCK workers, one BLAS thread each, one checkpoint per fit under
`results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04/`. Fits on all shared genes take
the longest; they are scheduled first. No fitter or update code changes.
