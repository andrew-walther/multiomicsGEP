# multiomicsGEP

Survival-supervised Bayesian matrix factorization for prognostic molecular
programs in pancreatic cancer.

The model represents a patient-by-feature molecular matrix as Y ≈ LFᵀ and
links the programs to survival through a Cox proportional hazards model.
L contains patient activities, F contains feature loadings, β contains survival
coefficients, and τ contains feature-specific noise precisions. Empirical Bayes
priors learn sparsity and shrinkage through coordinate ascent variational inference.

Two predictors are available: **YFB**, η = (YF)β, scores new patients by direct
projection, and **LB**, η = Lβ, uses latent activities. Expression is the current
single-modality application. A joint expression and methylation model is in development.

The [one-modality equivalence audit](docs/findings/2026-11-12-multimodal-m1-equivalence-findings.md)
records the remaining differences between the fitters and the KL regression checks.
Run `Rscript tests/test_multimodal_m1_equivalence.R` for the acceptance test,
which currently fails because the multimodal interface lacks required options.

The [modality-weight study](docs/findings/2026-11-12-modality-weights-findings.md)
derives reconstruction-only power weights for the grouped YFB model. Verify the
weighted equations with
`Rscript derivations/SSBMF-YFB-model-derivations/verify_derivation.R`.

## Getting started

Developed with R 4.5. Install the core dependencies, then run code from the
repository root:

```r
install.packages(c("survival", "ebnm", "yaml"))
```

Rcpp and a C++ compiler accelerate the loading sweep. A base R fallback is
available. See the [code guide](code/README.md#dependencies) for analysis dependencies.
Patient data are not distributed with the repository. Real-data runners read
cohorts from `PDAC_DATA_ROOT`.

### Fit and predict

This small synthetic example uses signed point-Laplace molecular priors and
a Normal prior on survival coefficients. It illustrates fitting and prediction.

```r
source("code/fit_ssbmf_yfb.R")
source("code/predict_ssbmf_yfb.R")

set.seed(1)
n <- 150; p <- 200
L <- matrix(rexp(n * 3), n, 3)
F <- matrix(rexp(p * 3) * rbinom(p * 3, 1, 0.2), p, 3)
Y <- L %*% t(F) + matrix(rnorm(n * p), n, p)
time <- rexp(n, exp(drop(L %*% c(0.8, -0.8, 0))))
status <- rbinom(n, 1, 0.7)

train <- 1:100
Y_train <- scale(Y[train, ])
Y_test <- scale(Y[-train, ], center = attr(Y_train, "scaled:center"),
                scale = attr(Y_train, "scaled:scale"))
fit <- fit_ssbmf_yfb(
  Y_train, time[train], status[train], K = 5,
  canonicalize = TRUE, w_S = 1,
  prior_LF = "point_laplace", prior_beta = "normal", verbose = FALSE
)

fit$history$converged  # check before interpreting a fit
head(fit$EF)           # p × K posterior mean loadings
fit$EBeta              # K posterior mean survival coefficients
risk <- predict_ssbmf_yfb(Y_test, fit$EF, fit$EBeta, fit$EF_norms)$risk_scores
```

`risk` has one Cox score per held-out patient. Higher scores indicate shorter
survival. The example applies training means and SDs to the held-out matrix.
The [code guide](code/README.md) describes other priors, prediction and stratified CV.

For pooled-cohort CV, use `fold_cohort_id = cohort` in the CV helpers to balance
cohort and event jointly. See the [fold example](code/README.md#cohort-and-event-folds).

Run the regression suite with:

```sh
Rscript tests/run_tests.R
```

A [small local results pilot](results/pilot/README.md) demonstrates synthetic
result rebuilding, versioned artifacts and freshness checks. It does not rebuild
the historical analyses.
Early benchmark tables and figures are preserved in the
[history archive](results/history/early-benchmarks/README.md).

The [matched LB/YFB runner](results/benchmark_sim/run_lb_yfb_study.R) compares
LB fold-in with YFB projection using identical inputs and five training folds.
Its [code findings and run commands](docs/findings/2026-11-12-lb-yfb-study-code-findings.md)
describe the explicit prior/scaling arguments, synthetic smoke mode and private checkpoints.
For the registered K = 5 comparison, pass `--K=5 --w_S=1` with
`--scaling=zscore --prior_LF=point_laplace --prior_beta=normal`. The runner
records rank and survival weight in `settings.txt`, and reports held-out
concordance and comparable-pair counts for each cohort.

Score the saved full-data pair with
`Rscript results/benchmark_sim/score_lb_yfb_external.R` after training finishes.
The scorer uses the dated private checkpoint and writes aggregate external
concordance tables, including donor-aware pooled intervals. Existing output
directories are never overwritten.

After both study modes finish, rebuild their aggregate summaries with
`Rscript results/benchmark_sim/summarize_lb_yfb_study.R --date=2026-11-12`.
It reads the saved checkpoints, applies the registered simulation gate and
writes program-activity, gene-overlap and score-correlation tables without
refitting either model.
The [K = 5 study findings](docs/findings/2026-11-12-lb-yfb-study-findings.md)
report the registered gate, gene comparisons, external intervals and numerical limits.



The [methylation-only study runner](results/multimodal_real/run_methylation_single_modality.R)
compares five CpG counts with expression on matched TCGA patients at the fixed
historical rank of seven. Andrew stopped this run to prioritize rank five.
Completed checkpoints are preserved. See the [study note](docs/findings/2026-11-12-methylation-single-modality-findings.md)
for the partial results, private checkpoints, commands and K = 5 handoff.

## Selected findings

The working rank-5 expression fit uses 2,218 genes from 273 TCGA and CPTAC
patients. Its mean external concordance is 0.630 across five dataset series,
compared with 0.631 for the rank-matched two-step fit. The donor-bootstrap
interval for the paired difference is −0.010 to 0.008. On real data,
survival supervision barely changes the loadings at the selected weight.
The training-only weight search selects w_S = 1, with nested CV concordance
0.658 and fold SE 0.027. Of its 216 fits, 27 inner fits reach the iteration
cap. All selected outer fits and the full-data fit converge.

The [rank-5 exhibit ledger](docs/findings/2026-11-12-k5-exhibits-findings.md)
lists the aggregate tables, figures, interpretation changes and rebuild
commands. Set `EXHIBIT_K=5 OMP_NUM_THREADS=1` for those commands to write
separate rank-5 outputs. With `EXHIBIT_K` unset, historical paths and ranks
are preserved. To summarize saved weight-search checkpoints, run
`EXHIBIT_K=5 OMP_NUM_THREADS=1 Rscript results/benchmark_sim/run_yfb_supervised_real_refit.R summary`.

In simulations with 2,000 genes, the planted prognostic program explained
0.05%–0.2% of centered molecular variance before gene standardization. Strong
survival effects gave mean oracle concordance of 0.936–0.958, using the true
generating risk score. The two-step fit recovered the program in 0/10 datasets
at 0.05%, 8/10 at 0.1% and 10/10 at 0.2%. These counts locate the empirical
recovery transition in this generator. Recovery required loading correlation
|r| ≥ 0.7, survival activity and the correct coefficient sign. At 0.05%, the
supervised fit at unit survival weight increased mean validation concordance
from 0.644 to 0.800 and recovered the program in 5/10 datasets. It still
averaged 2.5 false-active programs and failed the registered recovery criterion.
See the [simulation findings](docs/findings/2026-11-05-near-threshold-findings.md)
for the full comparison and limits.

## Repository and documentation

| Path | Contents |
|---|---|
| [code/](code/README.md) | Fitters, prediction, preprocessing and model selection |
| [config/globals.yml](config/globals.yml) | Hyperparameters and simulation settings |
| [tests/](tests/) | Regression tests |
| [results/](results/README.md) | Analysis runners, saved outputs and rebuild commands |
| [derivations/](derivations/README.md) | Full model derivations and algebra checks |
| [docs/](docs/README.md) | Plans, findings, reviews and progress-book source |
| [paper/](paper/writeup/README.md) | Writeup versions and manuscript material |
| [presentation/](presentation/) | Meeting slides |
| [longleaf_setup/](longleaf_setup/) | UNC Longleaf setup and SLURM scripts |

Documentation:

- [Progress notebook](https://andrew-walther.github.io/multiomicsGEP/)
- [November 13 writeup](paper/writeup/2026-11-13-writeup/README.md), with rank-5 exhibits and the descriptive joint evaluation
- [Chapter 4 numerical sources](paper/prelim/notes/chapter4-number-sources-2026-11-06.md)
- [Project status](PROJECT_STATUS.qmd)
- [Decisions](DECISIONS.md)
- [Roadmap](ROADMAP.md)

## Feature-group fitter tests

`fit_ssbmf_yfb()` accepts `feature_group`, one label per column of a concatenated
input matrix. Multiple groups fit separate F priors while sharing L, the joint
survival projection and feature-specific noise precisions. NULL and a single
group preserve the earlier numerical path. Standardize each block using its
training rows before concatenation. Grouped fits return
`prior_F[[group]][[program]]`; prediction uses the same concatenated column order.
Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript tests/test_feature_groups_yfb.R`
for the focused tests, also registered in the full suite. The
[findings](docs/findings/2026-11-12-feature-group-fitter-findings.md) explain
supported options and KL accounting. The registered simulation compares grouped
and pooled F priors, plus an expression-only diagnostic, at 50:50 and 50:250
expression:methylation features, with
15% or 3% nonzero methylation loadings. It regenerates survival from the
standardized reconstruction direction and preserves the original training risk SD.
A refitted point-mass prior empties that modality-program block and contributes
zero loading KL. Grouped fits return a modality-by-program `empty_blocks` count
matrix, including repeated resets of an already empty block. Other nonfinite KL
errors remain explicit. The rerun also reports validation C for the true risk.
Run the twelve-fit smoke check with
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript results/benchmark_sim/run_feature_group_simulation.R --smoke`.
Omit `--smoke` for the fixed 120-fit rerun. Each command refuses to overwrite
its dated output directory. Fits stay under ignored `data/`, while aggregate
results are written under `results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_rerun`.
The original 80-fit output directories are preserved.
Use `Rscript tests/test_feature_group_simulation.R` for its focused regression checks.

## Matched-TCGA multimodal comparison

This descriptive comparison is **not validated by simulation**. It requires
Task C's existing private preparation and reuses its 144 patients and five
folds without modifying it. Run
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript results/multimodal_real/run_multimodal_tcga_comparison.R --smoke`
for one fold, then omit `--smoke` for the full cohort and all five folds.
The full run reuses the smoke fold only when source and input hashes agree.
It retains capped fits. Aggregate tables carry the validation label, and
patient fits/scores remain under ignored `data/`.

For the separate two-direction `both_informative` simulation, run
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript results/benchmark_sim/run_feature_group_simulation.R --informative --smoke`.
Omit `--smoke` for all 120 fits in the new `_informative` directories. Program 1
controls the recovery gate, and program 2 is descriptive. Grouped and pooled
priors are compared only in cells where expression alone passes the gate.
The simulation does not change the TCGA comparison's descriptive status.

## License and author

[MIT](LICENSE). Patient data remain under their own terms of use.
Author: Andrew Walther.

## Training-only β-prior rerun at K = 5

With the prepared local training inputs, run
`OMP_NUM_THREADS=1 Rscript results/benchmark_sim/run_training_cv_redo.R brerun 3 5`.
This fits three shared β priors on the full training data and five saved folds,
then compares them with the existing per-program Normal rank-sweep checkpoints.
Patient-level checkpoints remain under ignored `data/`. The aggregate
`b_rerun_paired.csv` and the [β-prior findings](docs/findings/2026-11-12-beta-prior-comparison-findings.md)
record paired concordance differences and the registered one-SE CV criterion.
Run `Rscript tests/test_training_cv_brerun.R` for its synthetic regression checks.

For the matched-TCGA expression and methylation comparison at K = 5, use
`Rscript results/multimodal_real/run_methylation_single_modality.R rerun --K=5`,
then `Rscript results/multimodal_real/run_methylation_single_modality.R sweep --K=5`.
The `K05/` directories keep its private checkpoints and aggregate results separate
from the stopped K = 7 study. The [methylation findings](docs/findings/2026-11-12-methylation-single-modality-findings.md)
explain the revised timing screen, scope and rebuild commands.

For the approved K = 5 expression convergence follow-up, run
`OMP_NUM_THREADS=1 Rscript results/multimodal_real/run_expression_iteration_extension_k5.R run`.
The eight capped K = 5 expression fits restart with a 1,000-sweep cap in the private
`K05/iteration1000/` directory. Their first 300 ELBO values must reproduce before use.
The new K = 5 aggregate tables preserve the original 300-sweep results.
Then run `Rscript results/multimodal_real/characterize_methylation_factors_k5.R` to
describe saved K = 5 factors, CpG annotations, D4 expression enrichment and local
ABSOLUTE purity associations. It reads the approved 450K annotation package from
the ignored `K05/R-library/` directory and writes aggregate CSVs under
`K05/factor_characterization/`. No external scores or model adoption are involved.
The K = 5 synthetic checks are `tests/test_expression_iteration_extension_k5.R`
and `tests/test_methylation_factor_characterization_k5.R`.
