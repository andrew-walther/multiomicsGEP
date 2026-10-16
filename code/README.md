# code/

R implementation of the survival-supervised Bayesian matrix factorization model.
Scripts use paths relative to the repository root, so source them from there.
Hyperparameters are read from `config/globals.yml`.

## Model fitting

| File | Contents |
|------|----------|
| `fit_cox_on_yf.R` | `fit_cox_on_yf()`: CAVI fit of the YFB model (η = (YF)β), the recommended configuration. Optional cohort offsets, stratified baseline hazard, cohort-specific β and ELBO-based program pruning |
| `fit_modular.R` | `fit_supervised_mf_modular()`: CAVI fit of the LB model (η = Lβ) |
| `fit_modular_multistart.R` | Multiple-initialization wrapper for the LB fit |
| `deflation_init.R`, `warmstart_from_fit.R` | Initialization: deflation start, warm start of a smaller-K refit |

Each source file ends in a stand-alone runner that stops when no data are set,
after the functions are defined. Load one with
`tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) NULL)`.

## CAVI updates

One file per variational factor; the fitters source these.

| File | Update |
|------|--------|
| `update_L.R`, `update_L_surv_YFB.R` | q(L), patient scores (LB and YFB forms) |
| `update_F.R`, `update_F_surv_YFB.R` | q(F), feature loadings (LB and YFB forms) |
| `update_beta.R` | q(β), survival coefficients, via a second-order Taylor approximation of the Cox partial likelihood |
| `update_beta_cohort.R` | Cohort-specific β |
| `update_F_cohort.R` | Normal conjugate update for cohort offset columns |
| `update_tau.R` | q(τ), feature-specific noise precision |

## Prediction and evaluation

| File | Contents |
|------|----------|
| `predict_cox_on_yf.R` | YFB risk scores for new patients by projection, Y_new F β |
| `predict.R` | `predict_supervised_mf()`: LB hold-out prediction |
| `compute_elbo.R` | Full ELBO, including EBNM and survival terms |
| `compute_bic.R` | Joint log-likelihood and BIC across K |
| `compute_cv_loglik.R` | Cross-validated held-out survival and genomics log-likelihood |
| `concordance_ci.R` | Bootstrap confidence intervals for Harrell's C |

## Model selection

| File | Contents |
|------|----------|
| `select_K.R` | `auto_prune_K()`, `select_K_cv()`: number of programs |
| `select_alpha_cv.R` | Cross-validated supervision weight α |
| `select_k_alpha_bo.R` | Joint (K, α) search by Bayesian optimization |

## Data preparation

| File | Contents |
|------|----------|
| `preprocess_desurv.R` | DeSurv-style gene selection and per-platform standardization across cohorts |
| `feature_selection.R` | `cox_feature_selection()`: univariate Cox screening |
| `train_test_split.R` | `stratified_split()`: event-stratified train/test split |
| `load_multiomics_data.R` | `build_multiomics_cohorts()`: matched TCGA/ICGC expression, methylation and survival |

## Multimodal YFB (expression + methylation, in progress)

| File | Contents |
|------|----------|
| `fit_multimodal_yfb.R` | Joint fitter with per-modality loading priors and optional pruning |
| `multimodal_yfb_updates.R`, `multimodal_yfb_helpers.R` | Updates and helper calculations |
| `multimodal_yfb_sweep.cpp` | Compiled loading sweep (falls back to R if compilation fails) |
| `preprocess_multimodal_yfb.R`, `predict_multimodal_yfb.R` | Input checks and held-out scoring |
| `simulate_multimodal_yfb.R`, `run_multimodal_yfb_simulation.R` | Simulation generator and study runner |

## Interpretation

| File | Contents |
|------|----------|
| `pathway_enrichment.R` | Gene-set enrichment (fgsea, ORA) of fitted programs and subtype concordance |

## Diagnostic and sensitivity-study utilities

Opt-in helpers for specific diagnostic studies; the default fits do not use them.
Each study's design and outcome is recorded in `docs/plans/`.

`cox_working_logspace.R`, `profile_cox_intercept.R`, `loading_prior_fixed_q.R`,
`loading_balance_audit.R`, `multimodal_moment_audit.R`,
`multimodal_program_agreement.R`, `feature_sensitivity.R`,
`cpg_feature_sensitivity.R`, `multimodal_feature_sensitivity.R`

## Reference and archived code

- `Supervised_Bayesian_MF_V2.R`: earlier single-file implementation, kept for comparison. Do not extend.
- `legacy/`: original V1 implementation and early scripts.
