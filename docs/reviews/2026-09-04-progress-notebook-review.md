# Review of the September 4 progress-notebook chapter

**Reviewed target:** `main` at commit `a650491`  
**Chapter:** `docs/progress_book/chapters/2026-09-04.qmd`  
**Scope:** correctness, statistical methodology, and traceability; not prose or presentation order.

## Executive recommendation

Do not revise only the chapter. Correct the model and scoring issues first, rerun the affected cohort-beta comparison and uncertainty analysis, then update the chapter and `DECISIONS.md` from regenerated outputs. The highest-priority corrections are the cohort-beta external-prediction fallback, the Cox tied-time likelihood, and the estimand/bootstrap mismatch.

## Findings that require correction

### 1. Cohort-beta external prediction mixes incompatible sign orientations

**Locations:** `code/fit_cox_on_yf.R`, `compute_pooled_beta()` call near the returned fit; `code/update_beta_cohort.R`.

After Phase C may globally negate `EBeta`, the returned `EBeta_pooled` is computed using post-flip `rowMeans(EBeta)` as an initializer but pre-flip final `w`, `z`, and `ZF`. `compute_pooled_beta()` performs only one Gauss--Seidel sweep, so this is not merely an arbitrary sign convention.

On the cached D4 cohort-beta fit, recomputing the same fit with Phase C disabled leaves `EF` and the cohort beta matrix equal up to sign, but changes the pooled beta materially:

| Program | Phase C disabled | Cached Phase-C result |
|---|---:|---:|
| 3 | -0.0124904 | -0.0098190 |
| 6 | approximately 0 | -0.0293623 |
| 7 | 0.0403865 | 0.0203912 |

This can alter external risk-score rankings. Rerun the cohort-beta benchmark and all bootstrap outputs after fixing it.

**Recommended fix:** preserve a coherent pre- or post-Phase-C orientation for every input to the pooled update, and iterate the shared-beta fallback to convergence if it is intended to be a shared-beta refit. Add a regression test asserting that a global sign change in the fitted cohort-beta matrix changes `EBeta_pooled` only by a global sign.

### 2. Cox likelihood is wrong when event times are tied

**Location:** `code/fit_cox_on_yf.R`, `calc_cox_taylor_yf()`.

The code computes one risk-set denominator per sorted row. For tied events, Breslow requires the same denominator for all events at that time. The current result also changes if rows with the same event time are reordered.

A three-observation check gave log partial likelihood -1.864706 from the implementation versus -2.102889 from `coxph(..., ties = "breslow")`; permuting the two tied rows changed the implementation to -1.244592. The actual training data contain 9 TCGA and 12 CPTAC event rows whose times are shared with another row.

This affects fitting, `compute_joint_ll_bic()`, and `cv_survival_loglik()`.

**Recommended fix:** implement tied-time Breslow score, diagonal Hessian, and log partial likelihood consistently; compare each against `survival::coxph(..., ties = "breslow")` and add a row-permutation regression test.

### 3. The pooled bootstrap answers a different question from the reported mean external C-index

**Locations:** `results/benchmark_sim/run_cohort_beta_bootstrap_ci.R`, `code/concordance_ci.R`, chapter section 3.

The chapter reports the equally weighted mean of five within-cohort C-indices. The pooled bootstrap concatenates all 616 patients, includes between-cohort pairs, and gives patients in larger cohorts more weight. These are different estimands: for baseline minus cohort-beta the reported mean difference is 0.0135, while the pooled concordance difference is 0.0221.

Patient-only resampling also treats the five external cohorts as fixed. It cannot quantify uncertainty for transfer to a new cohort and may be optimistic if cohort effects matter.

**Recommended fix:** bootstrap paired within-cohort C-index differences and average them using the same weights as the headline metric. Report this as conditional on the five cohorts. If the scientific target is a new external cohort, use a cohort-aware analysis and state the severe small-sample limitation with only five cohorts.

### 4. The chapter overstates a non-significant result as equivalence

**Location:** `docs/progress_book/chapters/2026-09-04.qmd`, summary and section 3.

The supplied pooled CI for baseline minus `beta_cohort_id` is [-0.0093, 0.0505]. It includes a meaningful loss for cohort beta. A confidence interval that includes zero does not demonstrate equivalence, and the CI itself must be regenerated after findings 1--3 are addressed.

**Recommended wording after rerun:** “The point estimate favored the baseline; this analysis did not detect a difference, while the interval remained compatible with a meaningful loss.”

### 5. External C-index scoring selects its direction using external outcomes

**Locations:** `code/concordance_ci.R`; benchmark `oriented_cindex()` helpers.

The code evaluates `max(C, 1-C)` or flips each risk score after inspecting the same outcome data being scored. This makes the discrimination statistic invariant to a global sign but is not performance of a prospectively oriented risk predictor.

**Recommended fix:** establish the risk direction from training data only and retain it across every external cohort and bootstrap replicate. State clearly if a sign-invariant association metric is intentionally being reported instead.

### 6. Phase C also replaces posterior second moments incorrectly

**Location:** `code/fit_cox_on_yf.R`, Phase C sign flip.

When Phase C flips beta, it sets `EBeta2 <- EBeta^2`. A sign flip preserves the existing second moment E[beta^2]; replacing it with the squared posterior mean removes posterior variance. This affects shared and cohort-beta fits.

**Recommended fix:** negate only `EBeta`; leave `EBeta2` unchanged. Add a test that verifies posterior second moments are invariant to a sign correction.

## Chapter/output discrepancies

1. The chapter says ELBO, BIC, and in-sample log-likelihood prefer K=3. The CSV shows ELBO and BIC prefer K=3, but joint log-likelihood is highest at K=15 (-750479.9093 versus -768147.9308 at K=3). The plot script also uses the ELBO-optimal K marker for the log-likelihood panel.

2. The chapter says retained-factor counts stay flat around K=11 and K=13. The threshold output instead gives 3 active factors at each K with `beta_thresh=0.001`; at 0.01, K=11 has 3 and K=13 has 1. The external-C dip is not a threshold artifact, but this stated support is wrong.

3. The under-specified-K simulation statement is too broad. In the hybrid scenario, the joint minus EBMF mean C-index difference is +0.1129 at K=3 and +0.1571 at K=4, compared with +0.0952 at K=6. It is not flat or negative for every K below six.

4. The six hardcoded coefficient values in the cohort-heterogeneity table match the cached fit and `DECISIONS.md`. The prose 0.543 to 0.617 correctly rounds the source values 0.5431 and 0.6172. The remaining reviewed live tables trace to their CSVs.

## Claims that should be narrowed

### Partial pooling

`update_beta_cohort_k()` correctly computes cohort-restricted A/B terms and makes one vectorized `ebnm()` call per factor. With `prior_family = "normal"`, however, `ebnm()` defaults to a prior with mode zero. This shares an estimated prior scale across cohorts; it does not estimate a nonzero common coefficient toward which cohort coefficients shrink. Describe it as zero-centered shared shrinkage unless a hierarchical common-mean model is added.

### Unsupervised factor discovery

At `alpha_F=0`, the L/F/Tau update equations have no survival input; this was verified directly from `update_L_surv_YFB.R`, `update_F_surv_YFB.R`, and `compute_R_k()`. The stronger claim that the returned factorization is mathematically identical is not airtight because the full, survival-containing ELBO controls convergence in `fit_cox_on_yf.R`.

Suggested wording: “For fixed initialization and iteration count, the genomics updates are unsupervised; survival can affect which iterate is returned through the stopping rule.”

### Simulation evidence for cohort-specific beta

The simulation zeros each study-specific loading outside its owning cohort, then applies a shared beta to that factor. It demonstrates attribution where the other cohort has no factor exposure; it does not test recovery of genuinely different slopes for the same observed factor in both cohorts. Also disclose the documented unstable K=6 seed rather than reporting only robust medians.

### Gene-set overlap

The 2,064-gene intersection is the correct background universe for the two fitted gene lists. The hypergeometric p-values should remain descriptive/nominal: both factor sets were learned from the same expression matrix, the gene-selection events are correlated, and the report selects best matches among 40 EBMF factors.

## Recommended implementation sequence

1. Correct Cox ties and Phase-C second moments, with exact unit tests.
2. Repair the cohort-beta pooled prediction path and its sign-invariance test.
3. Decide and implement the external-performance estimand and orientation convention.
4. Rerun the cohort-beta benchmark, uncertainty analysis, and any K-sweep likelihoods affected by the tie correction.
5. Regenerate the figures/tables and revise the chapter plus `DECISIONS.md` from the new outputs.
6. Run the full project test suite, then add targeted tests for the new failure modes.

## Update after commits `e1fac72` through `6196b59`

The later commits add a strata-only arm, held-out survival likelihoods, a pooled-versus-single-cohort subanalysis, and a top-two factor-match diagnostic. They do not correct Findings 1--6 above. In particular, the new likelihoods retain the tied-time implementation and the cohort-beta CV path still uses the sign-inconsistent `EBeta_pooled` fallback.

### 7. New chapter claim: held-out survival likelihood does not give the same ranking as external C

**Locations:** `docs/progress_book/chapters/2026-09-04.qmd`, section 3; `results/benchmark_sim/outputs/cohort_beta_comparison/cohort_beta_heldout_survival_ll.csv`.

The chapter says held-out survival likelihood “reaches the same ranking” as external validation. It does not. External C ranks: `joint_yfb` (0.6267), `strata_only` (0.6263), `joint_yfb_beta_c` (0.6132), `joint_yfb_all_c` (0.6060), then `joint_yfb_cohort_L` (0.5431). Held-out log partial likelihood per event ranks: `strata_only` (-3.220025), `joint_yfb` (-3.221213), `joint_yfb_all_c` (-3.221584), `joint_yfb_beta_c` (-3.232926), then `joint_yfb_cohort_L` (-3.248388).

The criteria agree only that `joint_yfb_cohort_L` is worst. This disagreement is useful to report, rather than evidence that the methods agree. Regenerate the likelihoods after the tie correction and the stratified test scoring correction below.

### 8. New stratified-CV likelihood scores pooled risk sets on held-out rows

**Locations:** `code/compute_cv_loglik.R`, held-out `calc_cox_taylor_yf()` call; `results/benchmark_sim/run_cohort_beta_supplementary.R`.

The new CV support correctly subsets `strata_id` when fitting each training fold. It then calls `calc_cox_taylor_yf(pred$risk_scores, time[test_idx], status[test_idx])` without `strata = strata_id[test_idx]` when scoring the fold. A fit trained under stratified Cox risk sets is therefore evaluated under a pooled test-fold risk set. This affects `strata_only` and `joint_yfb_all_c`, central to the new claimed ranking.

**Recommended fix:** pass held-out strata to the scoring call and add a test that the stratified test-fold likelihood equals the sum of cohort-specific likelihoods. State that the stratified likelihood evaluates within-cohort ranking rather than cross-cohort survival-time comparisons.

### 9. New cohort-beta CV analysis discards known cohort membership at scoring time

**Locations:** `code/compute_cv_loglik.R`, cohort-beta path; `results/benchmark_sim/run_cohort_beta_supplementary.R`.

For `beta_cohort_id`, every held-out fold contains patients from the two known training cohorts, but the implementation scores with `EBeta_pooled`, as if the fold were an unseen cohort. That is a defensible new-cohort target, but it is not ordinary patient-level CV within the training cohorts, where known held-out labels can select the relevant beta column. It also depends on Finding 1's flawed fallback.

**Recommended fix:** decide which target is intended and label it. For within-cohort CV, score with `fit$EBeta` and held-out `beta_cohort_id`; for an unseen-cohort target, use the repaired pooled fallback. Reporting both would answer distinct questions.

### 10. The pooling subanalysis does not isolate pooling

**Locations:** `results/benchmark_sim/run_training_set_subanalysis.R`; chapter section 3.

The pooled fit has n=273 and 2,064 genes, whereas the single-cohort fits have n=129/144 and 3,000 cohort-selected genes. The script discloses this, but “Pooling clearly helps” is too causal: both sample pooling and feature selection change, and the outcome-informed external orientation from Finding 5 remains.

**Recommended wording:** “The pooled configuration outperformed these single-cohort configurations; this supplementary comparison changes both training sample size and feature selection, so it does not isolate the effect of pooling.” A gene-set-matched sensitivity analysis would isolate it.

### 11. The factor-merge result is a single-replicate exploratory diagnostic

**Locations:** `results/multi_cohort_sim/run_top2_match_diagnostic.R`; chapter section 4.

The diagnostic examines one saved example fit per scenario/K setting, then uses the post-hoc rule `second_cor > 0.5`. “7 of 27 versus 0 of 114” is not a result aggregated over the 15 simulated seeds. Call it a first-seed diagnostic consistent with merging under under-specification, or calculate it across all seeds and report the distribution before making the stronger claim.

### What the recent updates did correctly

- `cv_survival_loglik()` now subsets cohort-related vectors to training folds; its NULL path is tested for identity and the new passthrough tests prevent the prior dimension mismatch.
- The single-cohort script explicitly states that its gene sets are not matched.
- Adding a `strata_only` arm is an appropriate decomposition of the cohort-aware extensions.

These additions are worth retaining once the likelihood and scoring issues are corrected.
