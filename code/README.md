# code/

R implementation of the survival-supervised Bayesian matrix factorization model.
Scripts use paths relative to the repository root, so source them from there.
Hyperparameters are read from `config/globals.yml`.

## Model fitting

| File | Contents |
|------|----------|
| `fit_ssbmf_yfb.R` | `fit_ssbmf_yfb()`: CAVI fit of the YFB model (η = (YF)β), the recommended configuration. Optional cohort offsets, stratified baseline hazard, cohort-specific β, ELBO-based program pruning, and projection-only centering for nonnegative input (`center_projection`) |
| `fit_ssbmf_lb.R` | `fit_ssbmf_lb()`: CAVI fit of the LB model (η = Lβ) |
| `fit_ssbmf_multistart.R` | Multiple-initialization wrapper for the LB fit |
| `fit_cox_on_yf.R`, `fit_modular.R`, `fit_modular_multistart.R`, `predict_cox_on_yf.R`, `predict.R` | Aliases kept for historical scripts: each sources the renamed file and binds the old function name (`fit_cox_on_yf()`, `fit_supervised_mf_modular()`, the multistart wrappers, `predict_cox_on_yf()` and `predict_supervised_mf()`). New code should not use them |
| `deflation_init.R`, `warmstart_from_fit.R` | Initialization: deflation start, warm start of a smaller-K refit |

`source()` a fitter to load its functions. Each fitter also ends in a runner
block that executes only under `Rscript code/<file>.R`.

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
| `predict_ssbmf_yfb.R` | YFB risk scores for new patients by projection, Y_new F β |
| `predict_ssbmf_lb.R` | `predict_ssbmf_lb()`: LB hold-out prediction |
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

## Experimental supervision

`fit_yfb_supervised.R` provides the opt-in YFB fitter used in the October 30
supervision evaluation. See the [execution report](../docs/reviews/2026-10-30-codex-execution-report.md)
for the evaluated configurations and unresolved criteria.

## Model families and shared code

| Family | Fit and prediction | Dedicated updates |
|---|---|---|
| LB | `fit_ssbmf_lb.R`, `fit_ssbmf_multistart.R`, `predict_ssbmf_lb.R` | `update_L.R` and `update_F.R`. YFB also sources `update_L.R` for the reconstruction residual helper. |
| Historical YFB | `fit_ssbmf_yfb.R`, `predict_ssbmf_yfb.R` | `update_L_surv_YFB.R`, `update_F_surv_YFB.R`, `update_beta_cohort.R` |
| Experimental YFB | `fit_yfb_supervised.R` | `yfb_exponential_stable.R`, `yfb_diagnostic_controls.R` and the coupled loading sweep |
| Multimodal YFB | `fit_multimodal_yfb.R`, `predict_multimodal_yfb.R` | `preprocess_multimodal_yfb.R`, `multimodal_yfb_helpers.R`, `multimodal_yfb_updates.R`, `multimodal_yfb_sweep.cpp` |

LB and historical YFB source `update_beta.R`, `update_tau.R`,
`compute_elbo.R`, `update_F_cohort.R` and `deflation_init.R`. The experimental
and multimodal YFB paths share parts of the Cox and loading update machinery.
The historical YFB canonical path also calls the multimodal loading update.
The V1 and V2 monoliths are reference implementations.

Start a physical split by extracting and testing the shared reconstruction
residual helper. Align the Cox calculations only after comparing tie and
stratum handling, extreme predictors, and benchmark parity. Keep the current
entry paths callable during any directory migration.

## Dependencies

The core single-modality fit uses `survival` and `ebnm`. Analysis runners also
use `yaml` to read `config/globals.yml`. Rcpp and a C++ compiler provide the
compiled canonical loading sweep. If unavailable, the fitter warns and uses
the base R loop. Set `options(multimodal_yfb.engine = "R")` to select it explicitly.

| Use | Packages |
|---|---|
| Multimodal model and methylation loading | `Rcpp` with a C++ compiler, Bioconductor `impute` |
| Comparators and tuning | `flashier`, `glmnet`, `rBayesianOptimization`, `GPfit` |
| Pathway analysis and figures | `fgsea`, `msigdbr`, `ggplot2`, `pheatmap`, Bioconductor `clusterProfiler` |
| Optional quantile normalization | Bioconductor `preprocessCore` |

Install Bioconductor packages with `BiocManager::install()`.

## Stratified YFB cross-validation

For canonical YFB with projection centering, pass `strata_id`,
`canonicalize = TRUE`, and `center_projection = TRUE` to
`select_K_cv(model = "YFB", ...)` or `cv_survival_loglik()`. Each fold uses its
training means and within-stratum risk comparisons. Each held-out stratum
must be represented in that fold's training data. To score a previously
unseen stratum on its own, call `predict_ssbmf_yfb()` with a fixed vector
such as `projection_center = rep(0, ncol(Y_ext))`.
See the [review handoff](../docs/reviews/2026-11-06-stratified-yfb-cv-handoff.md)
for regression checks and the separate canonical scaling diagnostic.

## Experimental normalized-projection interface

The experimental `fit_yfb_supervised()` takes a named list containing one
nonnegative patient-by-gene matrix plus survival times and events. It scales
genes by within-cohort SD and returns posterior moments, factor roles,
per-sweep diagnostics and explicit failure flags. After sourcing its file,
the interface is `fit_yfb_supervised(list(expression = X), time, event, K = 7, gamma = 10)`.
This example documents an experimental interface, not an adopted configuration.

The returned `VBeta`, `VL` and `VF` fields store variances explicitly.
Use `sqrt(fit$VBeta)` for coefficient SD to retain small variances when means
are large. Second moments remain available for the update formulas.

The opt-in controls `tau_init`, `tau_fixed = TRUE` and `guard = TRUE` isolate
noise feedback and enable guarded sweeps. The guard uses a fixed Cox curvature
bound and restores all moments, priors and KL terms on rejection.
`backtrack_stalled` is a failure, not convergence. These controls preserve
default numerical results and are not adopted scientific defaults.
`yfb_joint_beta.R` is an exact fixed-Gaussian-design reference, not a joint
coefficient update integrated into this experimental fitter or full-Cox
coverage validation. The canonical fitter's separate joint Gaussian option
is documented in its derivation.

The [targeted findings](../docs/findings/2026-10-30-yfb-targeted-diagnostics-findings.md)
report partial stabilization with fixed noise and no complete screening success.
The [results guide](../results/README.md#historical-supervision-diagnostics)
retains the registered commands and checkpoint restrictions.
