# Independent review of Task F training CV

Date: 2026-11-12. Author: Andrew Walther.

**Verdict: the saved results reproduce and support D4, per-cohort z-scores and K = 5 under the registered rules. Full compliance with the week plan remains incomplete because the required additional fold seeds were not run.** K = 5 is a parsimonious choice on one split, not an established optimum or a demonstrated improvement over K = 7. No arithmetic error changes the three saved choices.

## Scope and independence

The review followed the [Task F prompt](../prompts/2026-11-12-codex-task-f-review-prompt.md). The registration, historical D1–D5 definitions, runner, aggregate outputs and private inputs were examined before reading the Task F findings, decision, result entry or status section. The initial verdict was recorded in `/tmp/task-f-review-2026-11-12/precomparison-verdict.txt` before that comparison. A separate statistical reviewer independently checked the code and aggregate calculations without reading those conclusions.

The checkout was `main` at `1bf8e6f60f0ef188a00b9d35cdf652e11056759f`. Existing unrelated work included uncommitted Task B additions to the Task F runner and methylation work. Those changes were not edited. Runner line references below refer to the **committed file at that revision**, before those additions. Saved grid provenance is `d721bb409ee1287e1fa124c6a9236487a4430402`. There were no committed changes between that provenance and the reviewed revision in `code/` or `run_lf_prior_comparison.R`. During review, concurrent work committed the Task B additions and an arm-isolation guard in `fcdf5f15`. That guard is acknowledged below, separately from the historical execution code.

Only this note and the reviews index are review deliverables. Patient-level files and scratch reruns remain outside Git. No fitter, selection rule, registration or model decision was changed. The review used at most two simultaneous R processes, with one BLAS/OpenMP thread per process. It did not run the full test suite because no production code changed.

## Choices and preprocessing checks

The [May 27 definitions](../findings/2026-05-27-desurv-alignment-report.qmd), lines 73–90, and `run_desurv_comparison.R`, lines 94–137, confirm that D1 and D3 use LB, D2 and D4 use YFB, and D5 adds a cohort indicator to D4. Reducing these to D2 and D4 is justified for the settled YFB model without a cohort indicator. It is the explicit narrowing in the November 12 registration, rather than a literal rerun of all five historical model configurations. D2 versus D4 changes normalization as well as gene selection, so its difference cannot be attributed to the gene-ranking rule alone.

`f_prepare()`, lines 44–66, calls the same D2 preprocessing configuration as the historical runner: the common gene universe, platform-specific single log, per-cohort z-scores, merged quantile normalization, then the 2,000 highest-variance genes. A fresh reconstruction from raw TCGA and CPTAC inputs reproduced the saved D2 gene order and matrix exactly, with maximum absolute difference zero. D4's raw combined-rank selection and intersection were also rebuilt. They reproduced the saved 2,218-gene matrix exactly. D2 and D4 share 258 genes.

The screens and transformations use expression and cohort labels without survival outcomes. Computing them on all 273 training patients is explicitly authorized. These are CV results conditional on preprocessing fitted to the full training cohort, rather than an evaluation of a pipeline that refits preprocessing within each fold. No external cohort enters these choices.

All three scaling inputs match the registration. The D4 z-score matrix agrees with independently standardized raw expression to within 4.5 × 10⁻¹⁵. X/sd and within-patient ranks followed by per-cohort z-scores reproduce exactly. All input matrices are finite. Only D4_XSD enables `center_projection` in `f_task()`, lines 103–108. Task A's separate positive-rank point-exponential arm remains descriptive and is not replaced by the signed rank-z arm.

The centered fits use training-fold means internally. The runner scores raw `Y %*% EF %*% EBeta` rather than subtracting the stored center. With one pooled center and no strata, this adds the same constant to every score, which cancels from concordance and Cox partial likelihood. It does not change these metrics. This argument would need reconsideration for stratum-specific centers or predictions requiring absolute score levels.

### Reference pairing

Task F and Task A have identical outcomes, cohort order, fold vectors and D4 reference inputs. All five saved Task A test-index vectors equal the indices selected by Task F's fold vector. The raw loader's ordered times and event indicators also agree with the saved inputs. The cache does not retain patient identifiers as an additional identity check.

`f_fit()`, lines 81–88, and `lfp_fit()` in `run_lf_prior_comparison.R`, lines 82–93, use the same PL call at K = 7: point-Laplace L/F, per-program Normal β, canonicalization, w_S = 1, alpha = 0.5, 300 sweeps, tolerance 10⁻⁵, no strata or pruning, and no centered projection. Both task functions set seed 42. Reuse is therefore appropriate. The PL reference includes two capped, nonconverged folds, which must remain visible.

### Independent calculations from fits.csv

The saved file contains 90 rows, covering 15 configurations with folds 0 through 5. Keys are unique, every CV configuration has five finite concordances, and no row records a fit error. Private checkpoints contain no nonempty fit warnings for these configurations. Successful return is distinct from convergence.

| Comparison at K = 7 | Mean within-cohort C | Difference from D4 z-score | Paired SE | Registered choice |
|---|---:|---:|---:|---|
| D4 z-score reference | 0.657842 | 0 | 0 | Keep |
| D2 preprocessing | 0.646338 | −0.011503 | 0.013921 | Keep D4 |
| D4 X/sd | 0.638107 | −0.019734 | 0.006578 | Keep z-score |
| D4 rank then z-score | 0.641171 | −0.016671 | 0.009011 | Keep z-score |

Each alternative must exceed the reference by more than one SE of the paired fold difference. None does. Within-cohort concordance uses pair-count pooling through `survival::concordance(... + strata(cohort), reverse = TRUE)`, with higher scores fixed as higher risk. Pooled and cohort-specific concordance and pooled held-out partial log-likelihood are also saved. ELBO is not used to compare inputs.

For rank selection, the registered calculation is

`threshold = max_K mean(C_K) − sd(C_best)/sqrt(5)`.

The best mean occurs at K = 10: 0.665899487, with fold SE 0.026942489. The threshold is 0.638956999. K = 5 through 12 qualify, and the smallest is **5**, with mean 0.654620082. K = 4 has mean 0.620338514 and does not qualify. Lines 218–221 implement the rule correctly.

The rank rule uses the marginal SE at the best K, not the SE of a paired difference. It includes fold variation shared across ranks. For K = 5 minus K = 10, the paired difference is −0.011279405 with SE 0.009805347. A different rule requiring each rank to lie within one paired SE of K = 10 would select K = 10, not K = 5. That is a sensitivity calculation, not grounds to change the registered rule after seeing results. K = 5 minus K = 7 is −0.003221642 with paired SE 0.006506817.

K = 5 converges in all five CV folds and in the full fit. K = 10 converges in only two folds and its full fit does not converge. Only K = 2 through 5 converge in all folds. Restricting descriptively to those ranks would still select K = 5. The registered rule did not prescribe that restriction. A single split and correlated training folds do not establish rank stability or provide a formal equivalence test.

Every rank retains K columns with nonzero loadings. At K = 5, the mean survival-active count across folds is 2.6 and the full-fit count is three. The number of retained columns and the survival-active count are therefore distinct from rank selection.

### Task A rerun at K = 5

The rerun is the registered real-data comparison, not a new simulation study. It compares point-exponential L/F on X/sd with centered projection against the PL z-score rank-sweep fits at K = 5, with the same folds and remaining options. PE has mean C 0.519192047 versus PL 0.654620082, a paired difference of −0.135428036 with SE 0.021221023. PE converges in zero of five folds, versus five of five for PL, and its maximum absolute cohort correlation is 0.998081. Task A's real-data criterion fails and PL remains the supported choice.

## Fresh fit reruns

Two fits were run sequentially from scratch using the saved inputs, current fitter and seed 42, with a new checkpoint directory under `/tmp/task-f-review-2026-11-12/`. The registered checkpoints were not overwritten or reused as fit results. The harness called `f_task()` and compared every saved metric except elapsed time, test indices, EF, EBeta, EBeta2 and held-out scores. All comparisons passed. Maximum absolute differences in the saved states and risk scores were zero.

| Fit, fold 1 | Sweeps | Converged | Within-cohort C | Pooled C | TCGA C | CPTAC C | Held-out pooled log-PL |
|---|---:|---|---:|---:|---:|---:|---:|
| D4 z-score, K = 5 | 53 | Yes | 0.675545 | 0.676505 | 0.668142 | 0.684492 | −154.9576 |
| D2 preprocessing, K = 7 | 56 | Yes | 0.702179 | 0.696576 | 0.712389 | 0.689840 | −154.5105 |

The first fit has five live and three survival-active programs. The second has seven live and two survival-active programs. Neither fresh fit emitted a warning. Aggregate reconstruction and fit-check logs are in `check.log` and `rerun.log` in the scratch directory. The private RDS files there contain patient-level information and must not be committed.

## Findings requiring follow-up

1. **[P2] The required repeat splits are missing.** The week plan, lines 132–135, requires four additional fold seeds when a paired difference is within one SE. D2 versus D4 meets that condition, as does K = 5 versus K = 7. The registration does not explicitly waive the requirement. The findings disclose its omission, but Task F is incomplete against that part of the plan. Repeats could change a choice. Preserve the single-split result and distinguish Andrew's adoption decision from evidence of stability.

2. **[P2, addressed during review] Summary inputs were not isolated by study arm.** `run_training_cv_redo.R`, lines 188–206 at the reviewed revision, reads every RDS in the directory, groups by input and K without `part`, and merges full fits by input and K alone. The Task B extension writes three additional β-prior arms as D4_Z/K5 to that directory. Without a guard, running the general summary after those files exist would mix priors, inflate the fold count and duplicate merged rank rows. This does not invalidate the saved 90-row Task F table. Concurrent commit `fcdf5f15` added an explicit exclusion of `b_rerun_` parts at current line 201, addressing that contamination route. Exact design-completeness checks remain a separate issue below.

3. **[P2] Missing checkpoints can yield decisions from fewer folds.** Lines 188–215 discover existing files rather than join to a required design. A missing CV checkpoint can disappear, and paired merges retain only the remaining shared folds. The summary reports counts but does not stop before selection. Require folds 1 through 5 exactly once, finite required metrics, and explicit error checks. This is a reproducibility defect rather than an observed omission in the reviewed outputs. Convergence should be reported separately, preserving the registered rule unless a new rule is prospectively approved.

4. **[P3] The explanation for bypassing select_K_cv() is partly inaccurate.** The registration and findings say the helper returns only the chosen K. `code/select_K.R`, lines 471–478, also returns the CV table and fold results. The valid reason for custom scoring is that its concordance grouping follows the model's `strata_id`, whereas Task F requires pooled fitting and within-cohort evaluation. Correct the return-value claim without changing these results.

The workflow appropriately keeps patient-level outputs out of Git and fixes score orientation. No security issue was identified in this bounded local review. Runtime was modest for the requested fits. Maintainability concerns are concentrated in summary completeness and arm isolation, rather than the numerical fitting call.

## Comparison with the reported findings

After recording the independent verdict, the [Task F findings](../findings/2026-11-12-training-cv-redo-findings.md), Task F decision, F result entry and status section were read. Their numerical choices, fold means, paired SEs, convergence counts and Task A rerun agree with the independent calculations. The findings explicitly disclose both the missing repeat seeds and capped high-rank fits.

The description of K = 5 through 12 as flat should be qualified as a descriptive range under the registered one-SE rule. It does not establish equivalence across those ranks. The paired-rule sensitivity above makes that limitation concrete. The K = 5 versus K = 7 comparison supports choosing the smaller rank under the rule, but does not show superior prediction. No claim about external performance at K = 5 follows from Task F. Existing K = 7 results must keep their model-state labels.

## Proposed changes, not applied

- Record Task F as reproducible under its single-split rules, with repeat-seed evaluation still pending. Complete the four additional seeds or explicitly record a decision to defer that requirement.
- Retain D4, z-scores and the adopted K = 5 as the current choices while describing the limits above. Keep historical K = 7 results separate from new K = 5 work.
- Correct the `select_K_cv()` return-value description, qualify the rank plateau language, and state that the D2/D4 contrast compares whole preprocessing configurations.
- Retain the concurrently added Task B arm-exclusion guard and add exact design-completeness checks before the general summary is reused. Preserve warning fields in future aggregate summaries so warnings remain reviewable without private files.
- Include K = 6 in the findings' convergence limitation. It has four converged folds, so the affected range begins at 6, rather than 7.
- Consolidate the Task F decision's obsolete sentence saying that config was unchanged with its later statement that `k_working_yfb: 5` was added. The key is present, and the historical K = 7 key is preserved.

No edits to DECISIONS.md, PROJECT_STATUS.qmd, the week plan or attempts log were made in this review.
