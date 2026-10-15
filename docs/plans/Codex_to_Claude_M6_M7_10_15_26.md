# Handoff to Claude: supervised robustness and comparators, 2026-10-15

## Completion status

The authorized five-prompt pass is finished, including stop-and-document actions
where the requested work could not proceed. **The full scientific work is not
complete:** E's adoption assessment and A's DeSurv comparison are blocked, B is
partial, and D produced verified tables but its process exited 1 after completion.
C's experiments finished, with substantive model failures rather than a favorable
validation result. Do not describe all five prompts as successfully completed.

Read `CLAUDE.md`, then `docs/reports/Codex_M6_M7_Summary_10_15_26.md` for the full
tables, implementation explanations and caveats. The original specification is
`docs/plans/Codex_Handoff_M6_M7_10_15_26.md`.

## Repository state and branch map

Work ran directly in the checkout, in order E → A → B → C → D. Each branch is
independently based on `feature/supervised-F`; these changes are **not cumulative**.
Nothing was pushed or merged, and neither `feature/supervised-F` nor `main` was
modified. Current checkout is `m6/robustness` at `2729b57`. It was clean after that
commit; this handoff is a later documentation-only addition, not part of that commit.

| Prompt | Branch | Commit | Full suite passed / failed |
|---|---|---|---|
| E | `m5/k-selection` | `ac572f6` | 553 / 0 |
| A | `m7/desurv-h2h` | `54b3df8` | 547 / 0 |
| B | `m7/sim-comparators` | `5ebf47f` | 548 / 0 |
| C | `m7/multimodal-comparators` | `33a716a` | 549 / 0 |
| D | `m6/robustness` | `2729b57` | 550 / 0 |

The consolidated report is committed on D. Other branches' tracked runners and
CSV files are not present in the D checkout. Inspect them without switching via
`git show <branch>:<path>`. Ignored local fit RDS files may remain across checkouts;
their presence does not mean their corresponding code was merged. Real data and
fit RDS files were not committed. Any integration requires review and separate
user authorization; do not automatically push, merge or overwrite base branches.

## Controlling configuration: unchanged

E did not adopt protection or a new starting rank. All subsequent experiments used
K_init=7 and the original handoff settings:

- Single modality: canonicalization, w_S=1, supervised initialization, shared
  point-Laplace β, beta_warmup=10, point-Laplace L/F, no pruning.
- Multimodal: intercept, ebnm prior updates, point-Laplace loadings in both
  modalities, ordinary pruning, plugin survival ELBO, supervised initialization,
  shared point-Laplace β, warm-up 10; no protected pruning.

The historically reported D4 model is a separate α_F=0 reference. Do not confuse
its external C≈.627 with the selected supervised models' new results.

## E: implementation exists; adoption is unassessed

Protected nullchecks and an opt-in single-modality output stopping rule were
implemented. Default fits were tested against the base branch. The protection
uses the literal inequality |Eβ| ≥ 1.96√(Eβ²−(Eβ)²), so a point mass at zero also
satisfies 0≥0; this boundary was tested and reported, not silently changed.

The fixed adoption rule requires every fit to converge under the ELBO rule. The
multimodal fitter stops on output changes and lacks the required full ELBO
convergence history. Adding a multimodal ELBO stopping rule exceeded the allowed
edits. The panel was stopped before a completed batch. **No empirical adoption
failure, passing K panel or selected 1-SE rank is claimed.** The conditional real
1-SE analysis was not run. Do not tune or adopt these options from unit tests.

Files on E: `results/benchmark_sim/run_k_selection_panel.R`,
`run_k_selection_real_1se.R`; outputs `results/benchmark_sim/outputs/k_selection/`.
Next decision: clarify/authorize how multimodal ELBO convergence is to be assessed
before resuming the fixed adoption experiment.

## A: no usable DeSurv implementation

Neither DeSurv package name was installed; no usable fitter was found in the repo,
local code or data folders. Authors' software/replication materials point to
`remotes::install_github("rashidlab/DeSurv")` and Bayesian optimization. No substitute
was written or installed. No head-to-head fits, five-cohort scoring or bootstrap
were performed. Compatibility with centered D4 inputs remains to be established.
The local bibliography's PNAS publication claim was not independently verified.

Audit on A: `results/benchmark_sim/outputs/desurv_h2h/BLOCKED.md`.
Next requirement: obtain a usable implementation and its documented training-only
rank/penalty selection procedure.

## B: complete low-noise results, conditional high-noise subset

Calibrated n_train=400, β3=3.3, seeds 20261016:20261035, noise .25 and 1.
At noise .25, selected supervised recovery was 17/20, mean C .7634; all EBMF,
EBMF→Cox/lasso and PCA comparators recovered 0/20. Ordinary EBMF→Cox mean C .6342,
PCA→Cox .6878; lasso's defined-only mean .6385 uses 19/20 scores.

At noise 1, eight one-factor EBMF results could not enter glmnet (needs ≥2 columns)
and two zero-factor fits failed. The whole-dataset error handler **lost earlier
arm rows and skipped later PCA/supervised arms** on those ten datasets. Remaining
results are conditional on 10/20 datasets: supervised recovery 9/10 versus 0/10
comparators. They cannot establish the full high-noise comparison. Undefined C
values and defined-only denominators are explicit in the report/CSV. Worker
warnings were not completely retained. No DeSurv arm was available.

Runner on B: `results/benchmark_sim/run_supervision_comparators.R`.
Outputs: `results/benchmark_sim/outputs/supervision_comparators/`.
Before a new run, agree on zero/one-factor handling and preserving independent arm
results; do not silently pad factors or reinterpret the conditional subset.

## C: rank saturation explained; selected real-data multimodal fit collapsed

The cached K ceilings 6/7/8 were passed correctly and all retained the identical
six EBMF factors, explaining identical C. Actual retained rank reporting was fixed
in result scripts. No core fitter was edited.

| TCGA → ICGC model | Retained K | Active | Converged | Primary C (n50) | All C (n67) |
|---|---:|---:|---|---:|---:|
| Selected multimodal | 7 | 0 | No, cap500 | NA | NA |
| Pruned 10/1 framework | 7 | 1 | No, cap500 | .7006 | .6647 |
| Stacked EBMF → Cox | 6 | 0 Wald | Not exposed | .4237 | .4182 |
| Selected expression-only | 7 | 1 | Yes | .5882 | .6163 |

Selected multimodal β first/second moments were exactly zero: constant risk,
no discrimination. The frozen helper returns NA; the older Harrell-tie helper
would return .5. Both joint fits capped. Flashier emitted a max-iteration warning;
its wrapper does not expose convergence. Do not claim successful validation.

Review caught separate scaling of primary/all ICGC expression subsets; this was
fixed **before expression scoring** by scaling all donors once then subsetting IDs.
Regression tests and all four arms' shared-donor audits show zero risk difference.
No model was retuned against ICGC.

Runner on C: `results/multimodal_real/run_supervised_multimodal_comparison.R`.
Outputs: `results/multimodal_real/outputs/supervised/` (comparison, rank, warnings,
scoring-consistency CSVs; local ignored per-arm fit checkpoints).
Next scientific priority: investigate collapse/convergence using training data or
simulation; do not choose remedies from the ICGC validation outcomes.

## D: robustness evidence; process error is retained

All 880 simulation fits are saved (40 datasets × default, 20 jittered, Normal L/F).
Highest final monitored ELBO chooses among 20 jitters, without using recovery.

| Arm | Recovery .25 / 1 (each /20) | Mean false-active .25 / 1 |
|---|---|---|
| Default | 17 / 19 | 1.10 / .25 |
| Best of20 | 17 / 17 | 1.15 / .45 |
| Normal L/F | 15 / 17 | .95 / .40 |

Best-ELBO selection loses two high-noise recoveries and gains none. All selected
starts converged, but 28/800 jitters and 2/40 defaults capped; all Normal L/F fits
converged. No alternate initialization/loading prior was adopted.

The unchanged training CV runner compares **β priors**, keeping L/F point-Laplace.
Refreshed mean C: point-Laplace .651199 vs Normal .653153; paired difference
−.001954 (SE .003055), still inconclusive. One capped fold per prior, previously
five. Prefix backups exactly match baseline; both old/new tables are retained.

| Features/cohort | Intersected p | Mean CV C | SE | Full active | Capped folds |
|---|---:|---:|---:|---:|---:|
| 1000 | 662 | .650200 | .017200 | 3 | 1/5 |
| 2000 | 1365 | .648933 | .022598 | 3 | 0/5 |
| 3000 | 2064 | .651199 | .023631 | 2 | 1/5 |
| 5000 | 3600 | .655413 | .027386 | 3 | 2/5 |

All four full-training fits converge. No feature-count change is supported.
Screening/scaling precede folds, matching existing code: absolute CV may be
optimistic. D does not load or score external cohorts.

**Nonzero exit:** After saving and printing every final table, Rscript exited 1
with `unexpected 'else' in "nteger(a[1])else"`. Codex edited a diagnostic warning
guard while Rscript was running, likely shifting its source read position. This
was an avoidable execution mistake. No fit/metric logic changed; no rerun or tuning
occurred. All 880 unique fits, 20 finite feature folds, four full fits and summary
arithmetic were independently verified. Current source parses and a fresh full
suite passes550/0. Report as verified saved results with an exit1, not a clean run.
The aggregate warning about saved warnings was spurious: none of the 880 rows
contains a model warning string.

Runner: `results/benchmark_sim/run_supervised_robustness.R`.
Outputs: `results/benchmark_sim/outputs/supervised_robustness/`, including
`execution_status.md`, all starts, best starts, CV prefix/postfix/comparison and
feature tables. The original `outputs/supervised_train_cv/` CSVs were refreshed.
Do not casually rerun: the runner protects existing prefix backups from overwrite.

## Review and remaining work

Read-only code/result review found no unresolved Important/Critical findings in
the recorded results, after correcting C's scaling problem. Unit-test success
does not establish model convergence, adoption or successful blocked comparisons.
Documentation was updated on each independent branch; D's consolidated report,
DECISIONS, PROJECT_STATUS and progress chapter gather cross-branch findings.

Before further implementation, discuss with Andrew: the E convergence requirement,
DeSurv installation/procedure, B low-rank failure handling, and investigation of C's
collapse/caps. A clean end-to-end D rerun would establish clean process completion,
but is not necessary to recompute the already independently verified summaries.
Do not automatically restart experiments, change methods, integrate branches or
push. This handoff is informational; obtain a concrete next task from Andrew.
