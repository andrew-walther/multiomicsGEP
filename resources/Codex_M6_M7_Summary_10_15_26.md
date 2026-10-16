# Supervised model robustness and comparator assessment — 2026-10-15

The five prompts were executed in order E, A, B, C, D, on separate branches.
Nothing was pushed or merged. Each branch starts from `feature/supervised-F` because E did not adopt a configuration change. The branches contain independent changes; this report gathers their findings without merging them.

| Prompt | Branch | Commit | Full-suite passed / failed | Outcome |
|---|---|---|---|---|
| E | `m5/k-selection` | `ac572f6` | 553 / 0 | Protected pruning implemented and tested; adoption assessment blocked |
| A | `m7/desurv-h2h` | `54b3df8` | 547 / 0 | Stopped: no installed or vendored DeSurv fitter |
| B | `m7/sim-comparators` | `5ebf47f` | 548 / 0 | Low-noise comparison complete; high-noise comparison partial |
| C | `m7/multimodal-comparators` | `33a716a` | 549 / 0 | Rank audit and four-arm real-data comparison complete; model failures reported |
| D | `m6/robustness` | This report's commit | 550 / 0 | All tables saved; trailing process error (exit 1) |

Counts differ because these are independent branches, not cumulative changes.

## E: configuration remains unchanged

The implemented protection excludes a coefficient from the pruning nullcheck when
|E[β]| ≥ 1.96√(E[β²] − E[β]²), including the literal zero-variance boundary 0 ≥ 0.
An opt-in single-modality output stopping rule also compares successive LF′ and YFβ;
default fitting remains unchanged. Tests compare default and explicit-false fits against the base branch and exercise protection when removal would increase the objective.

The adoption rule requires every fit to converge under an ELBO stopping rule. The multimodal fitter uses output convergence and does not expose the required full ELBO convergence history. Adding a new multimodal stopping rule exceeds the permitted edits. The panel was stopped before a completed simulation batch. This is an **unassessed convergence requirement**, not empirical evidence that protected pruning failed. No other adoption conditions were established. No tuning was attempted and the conditional real-data 1-SE step was not run.

Controlling configuration for A–D: K_init=7; selected canonical single-modality model with supervised initialization, shared point-Laplace β, warm-up 10, point-Laplace L/F, no pruning; selected multimodal model with ordinary ELBO pruning and no protection. The experimental E options have not been adopted.

## A: DeSurv could not be fitted

Neither DeSurv package name was installed, and no usable fitter was found in the repository, local code, or data folders. The authors' [software page](https://naimurashid.github.io/software/) and [replication repository](https://github.com/rashidlab/DeSurv-paper) point to `remotes::install_github("rashidlab/DeSurv")` and a Bayesian optimization procedure. A usable implementation, its dependencies and documented training-only selection procedure are needed before this comparison can proceed. Compatibility with the centered D4 preprocessing also needs to be established. No replacement implementation, fitting, external scoring or bootstrap was performed. The local citation is unpublished; a PNAS publication claim was not independently established.

## B: planted simulation comparators

Same calibrated data: n_train=400, β₃=3.3, 20 seeds per noise level. Recovery requires |cor|≥0.7 with the correct sign. Two-step programs are chosen by largest absolute ordinary-Cox z; the lasso fit supplies its coefficient for the sign check. Lasso nonzero counts are not Wald significance counts.

| Method | Noise | Datasets completed/planned | Recovered | C: mean over defined / count | Mean false active | Median retained K |
|---|---:|---:|---:|---:|---:|---:|
| EBMF alone | .25 | 20/20 | 0 | — | — | 2 |
| EBMF → Cox | .25 | 20/20 | 0 | .6342 / 20 | 1.15 | 2 |
| EBMF → lasso-Cox | .25 | 20/20 | 0 | .6385 / 19 | — | 2 |
| PCA → Cox | .25 | 20/20 | 0 | .6878 / 20 | 4.25 | 7 |
| Selected supervised | .25 | 20/20 | 17 | .7634 / 20 | 1.10 | 7 |
| EBMF alone | 1 | 10/20 | 0 | — | — | 2 |
| EBMF → Cox | 1 | 10/20 | 0 | .5283 / 10 | .60 | 2 |
| EBMF → lasso-Cox | 1 | 10/20 | 0 | .5535 / 6 | — | 2 |
| PCA → Cox | 1 | 10/20 | 0 | .5785 / 10 | 3.00 | 7 |
| Selected supervised | 1 | 10/20 | 9 | .7933 / 9 | .20 | 7 |

At noise 1, eight one-factor EBMF results could not enter glmnet (requires at least two columns), and two zero-factor fits failed. No padding or tuning was attempted. The whole-dataset error handler discarded earlier arm results and skipped the later independent PCA and supervised arms for those ten datasets. Thus the high-noise table is conditional on ten successful datasets, not a full or unbiased comparison. The low-noise supervised arm has one nonconverged fit; the successful high-noise subset has none. Undefined C values are retained as NA in full-denominator means; the table explicitly gives defined-only means and denominators. Worker warnings were not retained as a complete audit. DeSurv was unavailable.

## C: multimodal real-data comparison

The earlier starting-rank ceilings 6, 7 and 8 were passed correctly (`results/multimodal_real/run_multimodal_real_fit.R:144`, `code/run_multimodal_yfb_simulation.R:54`, `greedy_Kmax=K`). All retained six EBMF factors, with identical loadings, Cox coefficients and training risks (maximum difference zero). Greedy fitting saturated at six; the rank argument was not ignored. Reporting scripts now show the actual retained rank instead of NA.

New fixed K_init=7 fits use TCGA training and ICGC validation, with frozen orientation. Primary validation contains 50 donors; the all-donor set contains 67.

| Model | Retained K | Active programs | Converged | Primary C | All-donor C |
|---|---:|---|---|---:|---:|
| Selected multimodal | 7 | 0 | No, cap 500 | NA | NA |
| Pruned 10/1 framework | 7 | 1 (program 4) | No, cap 500 | .7006 | .6647 |
| Stacked EBMF → Cox | 6 | 0 (Wald) | Not exposed | .4237 | .4182 |
| Selected expression-only | 7 | 1 (program 2) | Yes | .5882 | .6163 |

The selected multimodal fit returned β first and second moments exactly zero and constant risk. The pre-fixed scoring helper reports NA for constant scores; the older multimodal helper would assign .5 for tied predictions. Either convention indicates no discrimination. Both joint fits reached the 500-iteration cap. Flashier emitted a maximum-iteration warning; its wrapper does not expose convergence, so that field is NA rather than an inferred success.

Before expression scoring, review detected that scaling primary and all-donor ICGC separately would give shared donors different predictions. The new helper scales all ICGC expression donors once and subsets by ID. A regression test and the saved four-arm scoring audit show maximum shared-donor risk difference zero. No core fitter was changed and no model was retuned against ICGC.

## D: selected single-modality robustness

All 40 simulated datasets completed their default, 20 jittered and Normal-loading-prior fits: **880 fits**. Jitter perturbs both L and F by independent Gaussian noise with SD .1 times the corresponding column SD, after reconstructing the original supervised SVD start. The selected start is the largest final monitored ELBO among the 20 jittered starts, without looking at recovery.

| Arm | Noise | Recovered / 20 | Mean false active | Nonconverged / 20 | Recovery differs from default |
|---|---:|---:|---:|---:|---:|
| Default point-Laplace L/F | .25 | 17 | 1.10 | 1 | — |
| Best of 20 jittered | .25 | 17 | 1.15 | 0 | 0 |
| Normal L/F | .25 | 15 | .95 | 0 | — |
| Default point-Laplace L/F | 1 | 19 | .25 | 1 | — |
| Best of 20 jittered | 1 | 17 | .45 | 0 | 2 losses, 0 gains |
| Normal L/F | 1 | 17 | .40 | 0 | — |

Across all 800 jittered fits, 28 reached the iteration cap; all 40 chosen best-ELBO starts converged. Two of 40 default fits reached the cap; all 40 Normal-loading fits converged. No model warning strings were recorded. The process launched before a diagnostic-only NA-warning guard was corrected; a spurious aggregate warning about retained warnings must not be interpreted as a model warning. High-noise validation-C means remain NA where constant risks make C undefined. No alternate start or prior is adopted.

The original training CV script was rerun unchanged. Old files are preserved as `*_prefix.csv`, refreshed files as `*_postfix.csv`.

| β prior (L/F point-Laplace in both) | Old mean C | Refreshed mean C | Old capped folds | Refreshed capped folds |
|---|---:|---:|---:|---:|
| Point-Laplace | .651322 | .651199 | 5/5 | 1/5 |
| Normal | .652394 | .653153 | 5/5 | 1/5 |

The paired point-Laplace-minus-Normal mean changes from −.001072 (SE .001515) to −.001954 (SE .003055). CV partial log-likelihood is now −779.067 vs −779.594. This remains an inconclusive comparison; the point-Laplace β choice stands. Fold 4 reaches 300 iterations under both priors.

| Top features per cohort | Intersecting genes | Mean fold C | SE | Mean active / fold | Capped folds / 5 | Full-training active | Full fit converged |
|---|---:|---:|---:|---:|---:|---:|---|
| 1000 | 662 | 0.650200 | 0.017200 | 1.2 | 1 | 3 | Yes |
| 2000 | 1365 | 0.648933 | 0.022598 | 2.4 | 0 | 3 | Yes |
| 3000 | 2064 | 0.651199 | 0.023631 | 2.2 | 1 | 2 | Yes |
| 5000 | 3600 | 0.655413 | 0.027386 | 2.0 | 2 | 3 | Yes |

The observed means differ by at most .0065, with SEs .017–.027. These results do not justify changing the selected 3,000-per-cohort setting. All four full-training fits converged.


Feature selection and standardization precede the folds, matching the existing CV protocol. Absolute C may therefore be optimistic; these are training-only sensitivity results, not a nested feature-selection validation. No external cohorts were loaded by D.


**Execution anomaly:** After writing and printing all simulation, CV and feature
summary tables, the D Rscript process exited 1 with `unexpected 'else' in
"nteger(a[1])else"`. A diagnostic warning-guard edit was made while Rscript
was still reading its source; this likely shifted its read position. This was
an avoidable execution mistake. The fit functions were already loaded before
that edit, which changed no fitting or scoring logic. No experiment was rerun
or tuned. The saved 880 distinct simulation fits, 20 finite-C feature folds,
four converged full-training fits and all summary arithmetic were independently
verified. The current on-disk runner parses successfully and passes the 550-test
suite. The process must nevertheless be reported as a **nonzero exit**, not a
clean successful run. `execution_status.md` preserves the exact error and scope.

## Implementation walkthrough and reproducibility

**E options and runners.** Protection modifies candidate eligibility in existing greedy nullchecks, passing its flag through refits; it does not change the gain calculation. The single-modality opt-in stopping criterion uses maximum relative output changes. Inputs and outputs retain existing fitter shapes (L: n×K, F: p×K, β: length K). Tests verify default identity and mathematically specified boundary behavior. The fixed panel and conditional 1-SE runner document their design, but their presence does not imply completed experiments. The zero-point-mass boundary and unavailable multimodal ELBO convergence are the main limitations.

**B runner and scoring test.** The runner loads existing calibrated generators, dispatches datasets through an R `parallel` socket cluster, fits flashier EBMF, survival Cox, glmnet lasso and base-R PCA, then applies the existing recovery and frozen-C conventions. Inputs are seed/noise pairs; outputs are per-arm CSV rows and an explicit failure table. The new test checks actual score direction/recovery behavior. Whole-dataset failure scope and differing valid-C denominators limit interpretation.

**C runner, scaling helper and reporting fixes.** The rank audit compares retained factors and predictions rather than trusting requested ceilings. The comparison runner uses fixed per-arm local checkpoints, fits four models once, and scores the two ICGC subsets with shared preprocessing. `standardize_icgc_expression()` accepts all-donor expression plus primary IDs and returns consistently scaled full/subset matrices; its regression test requires equal shared-donor predictions. Existing sweep reporting now uses actual EBMF factor count. Inputs are existing local expression/methylation and outcome loaders; outputs are CSV comparisons, warning/rank/scoring audits and ignored local fit RDS files. No data or RDS files were committed.

**D functions.** `make_jitter_start(tb,start)` reconstructs signed SVD, replaces program 7 by Yᵀu₀ and its residual least-squares L column, then adds reproducible column-scaled noise; it returns EL_init n×7 and EF_init p×7. `best_elbo_start(rows)` rejects nonfinite objectives and chooses the largest ELBO, with first-row ties. `load_robustness_code()` reuses the original generator/scorer rather than copying statistical logic. `robustness_dataset(seed,noise)` fits 22 prescribed arms, writes completed starts immediately and returns explicit errors. `run_supervised_robustness(cores)` dispatches 40 socket tasks, aggregates results, preserves and reruns the old CV script, then evaluates four feature counts with fixed folds and full-training active counts. CSVs expose convergence and undefined metrics. Existing survival/ebnm/yaml and base/parallel functions avoid new dependencies. Tests verify ELBO-only selection, repeatable jitter magnitude and direction, and exact unjittered custom/default fit equivalence. Reruns deliberately refuse to overwrite preserved CV backups; use a fresh output location/archive the old outputs explicitly before a new run.

Run commands on their respective branches:

```sh
Rscript tests/run_tests.R
# B
Rscript results/benchmark_sim/run_supervision_comparators.R
# C
Rscript results/multimodal_real/run_supervised_multimodal_comparison.R
# D
Rscript results/benchmark_sim/run_supervised_robustness.R 6
```

These are local socket experiments. Worker count must respect allocated CPUs and memory if subsequently run under SLURM. They were not submitted to Longleaf. Test counts confirm code checks, not scientific adoption, convergence or successful completion of blocked comparisons.
