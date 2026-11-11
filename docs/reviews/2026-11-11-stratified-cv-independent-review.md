# Independent review of stratified YFB cross-validation

Author: Andrew Walther
Date: 2026-11-11
Reviewed commits: `68e364662225cc88423047bdba4fa4628686c527` and `dbb411322b1cafb268aa1fd0471ac7186596307d`
Checkout tested: `85f1c81cbf90250453ebd6c488e788e19a7fe637`

## Verdict

Approved for the reviewed implementation scope. No blocking or nonblocking code defect was found in the two commits jointly. The item 0 completion condition in the November 9 week plan is met by this review: the independent note exists, the focused regression checks passed, and no defect requires correction. This approves the CV integration and cohort/event fold option. It does not establish the scientific performance of a model or select a prior, rank or scaling rule.

The review read the project instructions, the week plan's pre-task-A review and cohort/event requirements, the November 6 handoff, both diffs, the current implementations, their tests, the predictor and fitter centering code, and a current historical sweep caller. No patient data were read and no analysis study was run. Fits used only small synthetic regression fixtures.

## Implementation checks

| Requirement | Review evidence |
|---|---|
| Training centers only | Both helpers pass only training rows to the fitter. With strata, the fitter computes `rowsum(Y, stratum) / table(stratum)` on those rows. The predictors receive that returned matrix directly. Neither helper estimates a test-fold mean. Independent instrumentation checked the exact center matrix and risk vector in all four held-out folds of each helper. |
| Fold-specific strata | `select_K_cv()` extracts the full labels from `...`, removes them from the unsliced arguments, and subsets them separately for fitting and scoring. `cv_survival_loglik()` subsets its explicit `strata_id` argument on both sides. Both reject missing or wrong-length labels. |
| Within-stratum concordance | Both use `concordancefit(..., strata = ..., reverse = TRUE)` when strata are present. This pools comparable pairs within strata. KCV-T20 checks every held-out fold against the formula interface with `strata()`. CVLL-T20 independently checks the training diagnostic. |
| Held-out partial likelihood | The likelihood helper passes held-out labels and held-out scores to `calc_cox_taylor_yf()`. Its risk sets contain only the held-out members of each stratum. CVLL-T19 distinguishes this from pooled scoring. An independent tied-event fixture also matched a direct Breslow risk-set sum. |
| Missing training stratum | Centered prediction matches labels against the training center's row names and errors before metric computation if a label is absent. Both CV helpers propagated the expected unknown-label error for a singleton stratum. No substitute test mean or global center was used. |
| Unstratified compatibility | On an independent synthetic fixture, both current helpers returned objects identical to their pre-`68e3646` implementations with no strata, centering or fold cohort option. The fold generator's explicit-NULL regression also passed. |
| Cohort/event balancing | The opt-in generator shuffles original indices within cohort/event cells and assigns them round-robin, carrying the offset across cohorts within each event status. Thus cell counts and pooled event-status counts differ by at most one between folds. Original patient indices are preserved even for singleton cells. |
| Training cohort coverage | The generator explicitly checks each training complement for all supplied cohort labels and rejects an allocation that loses a cohort. CF1 and CF6 verify coverage, and CF5 verifies the sparse-cohort failure. |
| Argument separation | `fold_cohort_id` controls assignment only. It does not silently set genomic cohort offsets, Cox strata or coefficient cohorts. CF integration tests inspect actual training rows in rank CV for both models, alpha CV, survival-likelihood CV and genomic bi-CV. |

## Focused test results

These results were obtained afresh for this review. The previously reported 794-test full suite was not used as a substitute.

| Command | Result |
|---|---|
| `Rscript tests/test_select_K_cv.R` | 20 passed, 0 failed, exit 0. This script lacks a final count print, so the 20 passing test records were counted in its output. |
| `Rscript tests/test_compute_cv_loglik.R` | 20 passed, 0 failed, exit 0. |
| `Rscript tests/test_cohort_folds.R` | 11 passed, 0 failed, exit 0. |
| Projection-centering wrapper below | 11 passed, 0 failed, exit 0. |

All 62 focused repository tests passed. Their captured outputs contained no warnings or failing test records. The wrapper supplies the dependencies that the projection test expects from the full runner and makes assertion failures determine its exit status:

```sh
Rscript -e 'source("tests/test_helpers.R"); library(ebnm); library(survival); source("code/fit_ssbmf_yfb.R"); source("code/predict_ssbmf_yfb.R"); source("code/simulate_multimodal_yfb.R"); source("tests/test_projection_centering.R"); report_results("focused projection"); if (.test_env$failed) quit(status=1)'
```

Additional independent checks ran through `Rscript -` and exited 0. They used seed 725, 40 rows, eight columns, two equally sized strata, alternating event status, four folds and at most three fitting iterations:

1. Intercept every prediction call from both centered helpers. Identify rows against the synthetic input, recompute centers from the complementary training rows, and compare each held-out center matrix exactly. Compare scores numerically with `(Y_test - center[strata_test, ]) %*% sweep(EF, 2, EF_norms, "/") %*% EBeta`. All eight held-out calls passed. The first draft of this review-only assertion retained matrix-derived names in the oracle vector. Converting the oracle to an unnamed vector corrected that assertion, without changing project code.
2. Replace one stratum label with a singleton label. Call each centered helper without cohort-balanced folds. Both stop with `projection_strata contains level(s) not seen in training`.
3. Load both helper implementations with `git show 68e3646^:code/<file>` into a separate R environment. Compare complete returned objects with the current helpers on the same synthetic data, using unstratified defaults, four folds, ranks 1 and 2 for rank CV, and rank 1 for likelihood CV. Both comparisons were identical.
4. Check a six-row, two-stratum fixture with tied event times against the independent sum `sum_i delta_i * (eta_i - log(sum_{j: stratum_j = stratum_i, time_j >= time_i} exp(eta_j))))`. It matched the stratified Breslow helper within `1e-12`.

## Limits and follow-up

The numerical canonicalization diagnostic in the November 6 handoff remains an open scientific decision. A nearly constant within-stratum projection can produce a large loading rescale and small genomic precision, so the separate dead-update condition is not scale-invariant. The reviewed commits do not change that rule. Any replacement still needs its registered recovery, false-active and convergence evaluation.

The new fold option guarantees training coverage of the labels actually supplied as `fold_cohort_id`. It does not infer coverage of a different `strata_id` vector. For task A, pass the cohort labels explicitly to the fold option and the intended Cox labels explicitly to `strata_id`. The planned two-cohort use satisfies that contract. Sparse allocations can fail deliberately, even when another fold design might exist.

The historical `run_k_init_sweep.R` caller still omits the new fold option and therefore retains event-only assignment. This is consistent with the opt-in compatibility contract. A new task-A caller must request cohort/event folds explicitly. Historical and study-specific runner migration was outside this review.

Concordance and Cox partial likelihood are invariant to a constant score shift within a stratum. Passing either metric alone would therefore be insufficient evidence for correct centers. The direct prediction checks above and CVLL-T20 address that limitation. The full suite and algebra-verification script were not rerun, because this review changes no update equations or implementation files.
