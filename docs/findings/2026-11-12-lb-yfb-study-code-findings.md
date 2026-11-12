# Matched LB and YFB study code

Date: 2026-11-12. Author: Andrew Walther.

## Scope and result

The study code is implemented. No scientific simulation grid, training CV or external scoring was run in this session. Task E's scientific comparison remains pending the priors chosen by tasks A and B. The software checks use small synthetic matrices only. This note accompanies Prompt 2 of the [execution instructions](../prompts/2026-11-12-codex-c-d-e-prompts.md).

The new runner is `results/benchmark_sim/run_lb_yfb_study.R`. It uses `code/lb_yfb_study.R`, with an opt-in survival-block weight and explicit coefficient initialization added to `code/fit_ssbmf_lb.R`. Historical LB defaults remain available when `w_S = NULL` and `EBeta_init = NULL`. No modular update formula or YFB fitter was changed for E.

## Matched fitting and prediction

The runner reads one raw training matrix and one saved five-fold assignment. Each training complement determines the feature means and SDs within cohort. The resulting matrix is constructed once and supplied to both models. Z-score scaling subtracts the training-cohort mean and divides by its SD. SD scaling divides without subtracting a mean. A zero or invalid training SD causes an error rather than an unreported feature exclusion. Cohort labels balance or identify the saved folds and define evaluation pairs. They do not enter the pooled fits as Cox strata, molecular offsets or separate coefficients.

One SVD supplies the same L and F initialization. Point-exponential starts use absolute SVD entries. The start is canonicalized once to unit SD of the projected program scores. Both models receive exactly the same zero-variance second moments, zero β means and variances, and feature precisions. YFB's existing `posterior_init` interface prevents it from canonicalizing the start again. LB's `EBeta_init` prevents its historical Cox warm start. The check rejects different matrices, folds or starting states, and also compares LB's actual recorded starting state with the supplied state.

The predictors retain Cox orientation, with larger values indicating greater hazard. They are never flipped to increase observed concordance. Each held-out LB fold is predicted by `predict_ssbmf_lb_foldin()` under the fitted loading priors and untempered genomics precision. YFB uses `predict_ssbmf_yfb()` with the training projection center and factor normalization. A training predictive score for LB also uses fold-in. The requested correlation of fitted linear predictors instead uses fitted Lβ for LB and the fitted YFB projection predictor. These are different quantities and are kept separate.

Primary concordance uses `survival::concordancefit()` with cohort labels as scoring strata, pooling comparable pairs within cohorts by pair count. Pooled concordance is reported alongside it. Exact observed event times are preserved with `timefix = FALSE`. Constant scores and no comparable pairs yield undefined concordance, recorded as `NA`. The five-fold mean and paired YFB-minus-LB difference are undefined if any needed fold is undefined. The summary reports the number of defined folds and pairs instead of silently dropping them. The paired SE is the SD of the five fold differences divided by √5.

## Separating survival weight from genomics tempering

For the opt-in LB path, let G be the expected Gaussian reconstruction log likelihood and S the second-order expected Cox log partial likelihood. The monitored objective is

\[
G+w_S\{S-\mathrm{KL}_\beta\}-\mathrm{KL}_L-\mathrm{KL}_F.
\]

This matches the block-weight convention of YFB with `survival_block_weight = TRUE`. The L update has A = A_gen + w_S A_surv and B = B_gen + w_S B_surv. It reuses the existing update through α = 1/2 and divisors 1/2 and 1/(2w_S). The F update uses untempered genomics. The β likelihood and β KL have the same multiplier, so β's coordinate update has weight one. The current-posterior Cox expansion is recomputed after each sweep before monitoring the objective. The stored fold-in α is zero, which gives the same genomics precision for new patients.

The opt-in path accepts positive finite w_S and rejects an α schedule, normalized A/B, β boosting or β burn-in. Those controls have different meanings in the historical algorithm and are outside the matched study. The two models retain their existing update order and posterior approximations. Their ELBO values are not a between-model performance metric. No method or prior is adopted by adding this code.

## Function walkthrough and interfaces

All molecular matrices have patients in rows and features in columns. K is the number of fitted columns, not the number of nonzero or survival-active programs.

| Function or block | Purpose and logic | Inputs and outputs |
|---|---|---|
| `lb_yfb_study_assert_identical()` | Checks each shared design component with exact identity and names the mismatch. | Two design lists containing Y, fold and posterior state. Returns invisible TRUE or errors. |
| `lb_yfb_study_scale()` | Learns cohort-specific means and SDs from training rows and applies them to every row in that cohort. | Matrix, cohort vector, training indices, scaling name. Returns a matrix with the same ordering and dimensions. |
| `lb_yfb_study_start()` | Uses one SVD and one canonicalization to construct the common posterior start. It forms squared moments after rescaling to make floating-point rounding identical. | Training matrix, K, molecular prior and centering option. Returns posterior moments, precisions and empty persistent-prior fields. |
| `lb_yfb_study_concordance()` | Counts comparable pairs with fixed Cox orientation, within cohorts and pooled. | Aligned scores, outcomes and cohort labels. Returns a one-row data frame. |
| `lb_yfb_study_cor()` | Computes Pearson correlation only when both vectors vary. | Two score vectors. Returns a scalar or NA. |
| `lb_yfb_study_pair()` | Fits both models from the common state, predicts training patients, and computes predictor and program correlations. | Shared matrix, outcomes, cohort labels, priors and fitting settings. Returns fits, start, predictive scores and aggregate tables. |
| `lb_yfb_study_cv()` | Validates the full outcome vectors before subsetting, then scales and fits each complement once. Held-out scoring uses fold-in for LB and projection for YFB. | Raw matrix, aligned outcomes/cohorts, five fold labels and model settings. Returns fold metrics, paired summaries and private paired fits/scores. |
| LB `w_S` and `EBeta_init` blocks | Separate likelihood weighting from α tempering, record the actual start, and monitor the current-posterior objective. | Existing fitter inputs plus optional positive weight and coefficient vector. Return fields add weight, initial state and coordinate negative KL values in the opt-in path. |
| Runner option and design helpers | Require explicit prior/scaling choices and enumerate the registered simulation design before fitting. | Command-line strings. Return validated options and 104 dataset-weight combinations. |
| Runner planted-data and loading-truth helpers | Reproduce the November 6 generator and transform molecular truth by inverse feature SD. | Seed/amplitude, or loading vector and training matrix. Return synthetic data or transformed truth. |
| Runner cell capture and simulation loop | Save errors as explicit failed model rows and private checkpoints, continue through the registered design, then report a terminal error after aggregate tables are saved. | Design rows and settings. Return aggregate metrics plus private fits and validation scores. |
| Runner write and main functions | Refuse reused output/checkpoint directories, record settings/source hashes, dispatch the requested mode, and separate aggregate tables from private objects. | Validated options. Write CSVs/settings under the requested output and patient-level objects under ignored `data/`. |

Base R supplies SVD, scaling and correlations. `ebnm` remains the posterior and empirical-prior solver used by the fitters and LB fold-in. `survival` supplies Cox calculations and comparable-pair concordance. `yaml` reads existing configuration values for rank and iteration controls. No new dependency was introduced. All execution is serial with one BLAS/OpenMP thread, avoiding interference with the concurrent studies.

## Registered runs and limits

Simulation mode contains 104 paired cells, or 208 model fits. It includes the three S-9 shares 0.25%, 1% and 5% and the S-14 strong-survival share 0.05%, each for seeds 1 to 10 and w_S = 1 and 10. It also includes the November 6 matched generator at amplitudes 0.5 and 1, seeds 1 to 6 and both weights. The latter has four true factors. Both fitted models use the configured seven columns. The S-9 generator historically SD-scales its training and validation cohorts separately before returning them. That preparation is preserved and documented. Any additional study scaling uses the returned training matrix for both datasets, with molecular loading truth divided by the additional feature SD.

S-9/S-14 recovery uses molecular loading correlation, coefficient activity and the registered sign rule. An unmatched or dead factor does not receive a recovery success. The November 6 comparison retains its continuous, one-to-one matched mean absolute score correlation for the three true survival programs. It has no newly registered binary recovery gate, so its `recovery_success` is NA. Its false-active count exempts only estimated programs matched to true survival programs at absolute correlation at least 0.7. These continuous matched-score results are descriptive. The binary carry-forward comparison uses the registered S-9/S-14 cells. A survival-active coefficient needs positive posterior SD and |E[β]| at least 1.96 SD, so a zero mean and zero SD do not count as activity.

Each failed simulation cell is retained with its error message and a failed status for both arms. If an error interrupts a paired cell after one fit, the whole paired cell is marked failed rather than presenting a partial pair as complete. Successful cell checkpoints retain fits, validation predictions and outcomes under ignored `data/`. Aggregate simulation tables contain no patient scores or matrices. Training mode checkpoints every completed or failed fold under ignored `data/`. It writes CV aggregate tables before starting the full-data pair. A later fold error produces explicit failed rows and a terminal error after CV evidence is saved. A full-data error saves its message and failure rows while preserving completed CV evidence. Full-data program and predictor correlations are descriptive. Top-gene and DeSurv-overlap interpretation still requires the scientific fits and review.

## Commands after the priors are fixed

Run from the repository root. Set these three shell variables to the registered choices from A and B. For example, zscore and point_laplace describe A's current PL setup. The β choice must come from B.

```sh
E_SCALING=zscore
E_PRIOR_LF=point_laplace
E_PRIOR_BETA=normal  # replace with B's registered choice before scientific fitting

Rscript results/benchmark_sim/run_lb_yfb_study.R \
  --mode=simulation --scaling="$E_SCALING" \
  --prior_LF="$E_PRIOR_LF" --prior_beta="$E_PRIOR_BETA" \
  --out=results/benchmark_sim/outputs/lb_yfb_simulation_2026-11-12

Rscript results/benchmark_sim/run_lb_yfb_study.R \
  --mode=training --scaling="$E_SCALING" \
  --prior_LF="$E_PRIOR_LF" --prior_beta="$E_PRIOR_BETA" \
  --input=data/lf_prior_real_2026-11-12.rds \
  --out=results/benchmark_sim/outputs/lb_yfb_training_2026-11-12
```

The training input is the existing ignored preparation from A, containing raw D4 X, outcomes, cohort labels and the saved S-11 fold labels. If it is absent, prepare it through A's authorized workflow first. This runner does not load external cohorts or choose priors. For a later SD-scaled point-exponential choice, use `E_SCALING=sd` and `E_PRIOR_LF=point_exponential`. That selects centered YFB projection while keeping reconstruction input uncentered. Use fresh output directory names for later runs. Existing outputs and private checkpoint directories are never overwritten.

A tiny reproducible software smoke run is:

```sh
Rscript results/benchmark_sim/run_lb_yfb_study.R \
  --mode=smoke --scaling=zscore \
  --prior_LF=point_laplace --prior_beta=normal \
  --out=/tmp/lb_yfb_study_smoke_new
```

## Verification

Nine focused study checks and three runner checks pass. A fresh independent review approved the final patch, including the later-fold failure check. They cover exact-input rejection, training-only scaling, pair-count-weighted concordance with unequal cohort sizes, fixed anti-concordant direction, a hand-computed weighted LB sweep with nonzero initial β, the final-state objective, invalid controls and outcome-vector lengths, exact shared initialization, prediction formulas and all five held-out folds. An injected failure in fold 3 checks that all five checkpoint files are retained, earlier fits survive, failure rows name the error and incomplete summaries stay undefined. Runner checks cover settings/design, failure retention and inverse-SD molecular truth. A mutation check removed survival weighting from a temporary copy of the L update. The hand-computed sweep test failed with a maximum L difference of 39.25, while the other eight checks passed. The repository source is unchanged by this mutation check.

Two smoke runs used a 40 by 20 synthetic matrix, two fitted programs, five folds and a two-iteration cap. PL with z-score scaling and PE with SD scaling both returned their aggregate tables. Each correctly reported that the fitted models reached the deliberate two-iteration cap. LB fold-in converged in both full-data smoke fits. These are software checks, not evidence that either model passes a scientific convergence or recovery gate.

The final full suite, `Rscript tests/run_tests.R`, passed 816 tests with zero failures. It emitted the existing bounded β optimizer profile-likelihood fallback warning. Fresh independent final review approved the complete patch after the per-fold checkpoint correction. Final PL and PE smoke runs also passed on the reviewed code, with the expected two-iteration fit caps reported and fold-in converged. Aggregate smoke evidence is under `/tmp/lb_yfb_study_final_pl_20261112/` and `/tmp/lb_yfb_study_final_pe_20261112/`; their private checkpoints are under ignored `data/`. No scientific fits were run.

## Text for the shared status files

Task D steps 1 and 2 have been audited separately. Its KL edge cases are corrected, but its M = 1 equivalence gate remains open and needs a fitter-design decision. Task E's matched-input, five-fold fold-in/projection CV code and configurable runner are implemented and tested on synthetic data. LB has an opt-in untempered-genomics survival-block weight and explicit β initialization to support the matched comparison. Historical defaults are unchanged. No scientific E grid, training CV or external fit was run. Run the simulation and training commands only after A and B fix the priors, then independently review the results before entering them in the research documents. The registered recovery, false-active and convergence gates remain scientific requirements.
