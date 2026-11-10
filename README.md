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

## Selected findings

The reported expression fit uses 2,218 genes from 273 TCGA and CPTAC patients.
Its mean external concordance is 0.620 across five dataset series, compared
with 0.613 for the rank-matched two-step fit. The donor-bootstrap interval
for the paired difference includes zero. In a simulation with weak molecular
signal and strong survival effects, supervision improved prediction near the
detection threshold. It barely changes real-data loadings at the selected weight. See the [status report](PROJECT_STATUS.qmd)
for the findings and their limits. The earlier benchmark and its estimates
remain documented as historical results.

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

Read the [progress notebook](https://andrew-walther.github.io/multiomicsGEP/),
[project status](PROJECT_STATUS.qmd), [decisions](DECISIONS.md) and
[roadmap](ROADMAP.md) for the research record and continuing work.
The [prelim Chapter 4 source ledger](paper/prelim/notes/chapter4-number-sources-2026-11-06.md)
includes its exhibit rebuild commands.

## License and author

[MIT](LICENSE). Patient data remain under their own terms of use.
Author: Andrew Walther.
