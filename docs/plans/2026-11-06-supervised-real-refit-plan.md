# Attempt S-15: real-data refit with the survival weight chosen by cross-validation

Registered 2026-11-06, before any fit. Andrew asked for this refit after reviewing S-13 and S-14
(2026-11-06). It follows the 2026-11-05 decision entry: nothing from S-13 is adopted by default,
and the candidate for real data is the survival-block weight (A1) with w_S chosen per dataset by
five-fold CV partial log-likelihood.

## Why

- S-11 fixed w_S = 1 on the real training data. At that weight survival supplies 10⁻⁷ to 10⁻⁶ of
  the loading precision, so the loadings equal the unsupervised ones.
- In S-13 the CV rule chose w_S between 1 and 10, rising with signal strength, and never a weight
  at which fits stop converging.
- S-14 showed supervision helps most where a prognostic program is near the detection threshold.

This attempt asks what weight the same rule chooses on the real data, and whether the result
predicts at least as well as w_S = 1.

## Data and model

- **Data:** the S-11 prepared training matrix (`data/yfb_supervised_real_2026-11-04.rds`, not in
  Git). It holds TCGA_PAAD and CPTAC (n = 273), single log and per-cohort z-scores, restricted to
  the 2,218 D4 genes. The outer folds are the S-11 five folds, stratified by cohort × event. No
  external cohort is loaded.
- **Model:** `fit_ssbmf_yfb(canonicalize = TRUE, survival_block_weight = TRUE, w_S = w)`, with
  the S-11 settings otherwise (K = 7, point-Laplace L and F, Normal β, `config/globals.yml` CAVI
  settings). At w = 1 this is identical to the S-11 fit.
- **Weights:** {1, 2, 5, 10, 20, 50}, as in S-13.

## Procedure

1. **Weight curve and final weight.** For every weight, fit the five S-11 folds and the full
   training set. The CV partial log-likelihood (Verweij–van Houwelingen) is summed over folds.
   The final weight w* maximizes it among weights whose five fold fits ran without error. Ties
   go to the smaller weight.
2. **Honest estimate (nested CV).** Within each outer training set, run the same rule on inner
   five-fold splits, stratified by cohort × event with `set.seed(2026 + outer fold)`. Then fit
   the outer training set at that fold's chosen weight and score the outer test fold. The mean
   of the five outer C values is the estimate for the procedure. Step 1's C at w* is not, since
   w* was chosen on those folds.
3. **Reproducibility check.** The w = 1 fold fits must reproduce the S-11 w_S = 1 fold C values
   (mean 0.663), since A1 equals A0 at w = 1.

## Reported

- The weight curve: CV partial log-likelihood, CV C, convergence, and the median survival share
  of the loading precision for live and active programs.
- w*, the five per-outer-fold chosen weights, and the nested-CV mean C with its SE. These are
  shown beside the S-11 w_S = 1 CV C (0.663, SE 0.024) and the S-11 two-step CV C (0.652).
- For the full fit at w*:
  - program roles (prognostic, survival-only direction, genomics-only, inactive; PVE cut 1%);
  - each active program's standardized coefficient (log hazard ratio per training SD of its
    score);
  - the share of the training risk-score variance carried by each program;
  - the best |r| between each program's loadings and the S-11 unsupervised (w_S = 0) loadings.

## Decision rule (registered)

- **Adopt** w* as the real-data configuration if all three hold:
  1. w* > 1;
  2. the full fit at w* runs without error and reaches convergence or the iteration cap with no
     ELBO decrease larger than 10⁻⁶ relative;
  3. the nested-CV mean C is at least the S-11 w_S = 1 CV C minus one SE (0.663 − 0.024 = 0.639).
- **If w* = 1:** the S-11/S-12 model stands, and no new external scoring is done.
- **If adopted:** run one external scoring of the fixed full fit at w*, with the S-12 procedure
  and script (`results/benchmark_sim/score_yfb_supervised_external.R`, given the new fit path).
  The S-11 w_S = 1 fit and the two-step fit are the comparisons. This happens only after the
  training results above are committed. Nothing is chosen from the external results.
- **False actives:** following Andrew's 2026-11-06 guidance and the S-13/S-14 contribution
  analysis, false activity is reported as contribution, not as a count threshold. That means the
  standardized coefficient and the risk-variance share of every active program other than the
  dominant one.

## Software

- **Runner:** `results/benchmark_sim/run_yfb_supervised_real_refit.R`. It runs 216 fits: 6 outer
  levels (all patients and five outer folds) × 6 inner levels × 6 weights. It uses PSOCK workers
  with one BLAS thread each, and writes one checkpoint per fit.
- **Output:** aggregate CSVs in `results/benchmark_sim/outputs/yfb_supervised_refit_2026-11-06/`.
  Checkpoints (`.rds`) are git-ignored.
- No fitter or update code changes. Training patients only.
