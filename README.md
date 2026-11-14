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

Run the regression suite with:

```sh
Rscript tests/run_tests.R
```

The experimental [share-target weight calculator](code/share_targeted_ws.R)
inverts the gene-precision-weighted directional curvature ratio at a converged
training state. Recalibration tracks the loading selected at the unit-weight
pilot even if its activity changes. Run `Rscript tests/test_share_targeted_ws.R` for its behavioral checks
and `Rscript results/benchmark_sim/audit_share_targeted_ws.R` for the saved-state
audit. The audit needs the existing local checkpoints and prepared training
matrix. See the [S-16 findings](docs/findings/2026-11-13-share-targeted-weight-findings.md)
for Amendment 1, the registered convergence/null gates and interpretation limits.
The study uses `run_share_targeted_diagnostics.R` for raw-weight and optional
50-sweep ramp diagnostics, and `summarize_share_targeted_diagnostics.R` to
rebuild their aggregate summaries. `run_share_targeted_weight.R` provides
`smoke`, `screen`, `real` and `structure` modes with at most two workers.
Fit checkpoints and settings belong under ignored `data/`. Only targets
allowed by the diagnostic and null gates proceed to real-data evaluation.
The completed S-16 fallback did not pass its primary null gate. Unit weight
remains the working model, with no new external scoring or model adoption.

The registered [fixed-prior supervision study](docs/plans/2026-11-14-bounded-risk-supervision-plan.md)
tests fixed Normal and Laplace priors on the survival coefficients. The opt-in
`stable_convergence_window = 50L` control requires the existing ELBO criterion
and the registered 50-sweep risk-score stability check. Its default is zero.
The screen reuses the share-target weight calculator and retains every
calibration pass. Run the two prior families with at most two concurrent R
processes and one BLAS thread each:

```sh
Rscript results/benchmark_sim/run_bounded_risk_screen.R normal real
Rscript results/benchmark_sim/run_bounded_risk_screen.R laplace real
# After both real-data partitions finish:
Rscript results/benchmark_sim/run_bounded_risk_screen.R normal sim
Rscript results/benchmark_sim/run_bounded_risk_screen.R laplace sim
Rscript results/benchmark_sim/run_bounded_risk_screen.R summarize
```

This requires the existing prepared training data under ignored `data/`.
The study is paused at Andrew's request after 60 saved fits. Follow the
[resume handoff](docs/reviews/2026-11-14-bounded-risk-resume-handoff.md).
Patient-level checkpoints stay there. Aggregate traces and gate tables go to
`results/benchmark_sim/outputs/bounded_risk_2026-11-14/`.

## Selected findings

The working rank-5 expression fit uses 2,218 genes from 273 TCGA and CPTAC
patients. Its mean external concordance is 0.630 across five dataset series,
compared with 0.631 for the rank-matched two-step fit. The donor-bootstrap
interval for the paired difference is −0.010 to 0.008. On real data,
survival supervision barely changes the loadings at the selected weight.
See the [rank-5 findings](docs/findings/2026-11-12-k5-exhibits-findings.md) for details.

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
- [Latest writeup](paper/writeup/2026-11-13-writeup/README.md)
- [Project status](PROJECT_STATUS.qmd)
- [Decisions](DECISIONS.md)
- [Roadmap](ROADMAP.md)

## License and author

[MIT](LICENSE). Patient data remain under their own terms of use.
Author: Andrew Walther.
