# multiomicsGEP

**Survival-Supervised Bayesian Matrix Factorization for Prognostic Molecular Programs in Pancreatic Cancer**

## Overview

This project develops a supervised Bayesian matrix factorization model that
decomposes a high-dimensional molecular matrix (gene expression, and in an
extension DNA methylation) while jointly modelling patient survival with a Cox
proportional hazards model.

Unsupervised factorizations such as PCA or NMF recover the structure that explains
the most molecular variance, which need not be the structure that predicts
outcome. Adding a survival likelihood to the factorization targets latent
programs that are both molecularly coherent and prognostic.

```
Y (n×p) = L (n×K) F' (K×p) + E          matrix factorization (molecular data)
h(t_i)  = h₀(t_i) exp(η_i)              Cox proportional hazards (survival)
```

- **L**: patient scores on each latent program
- **F**: feature (gene) loadings that define each program
- **β**: survival coefficient of each program
- **τ**: feature-specific noise precision

Two linear predictors are implemented:

| Model | Linear predictor | Notes |
|-------|------------------|-------|
| LB  | η = Lβ     | Patient scores are learned jointly with survival (`code/fit_modular.R`) |
| YFB | η = (YF)β  | Predictor computed from the observed data and the loadings, so new patients are scored by projection without refitting L (`code/fit_cox_on_yf.R`) |

Inference uses coordinate ascent variational inference (CAVI). Each coordinate
update is an empirical Bayes normal means (EBNM) problem solved with the
[`ebnm`](https://cran.r-project.org/package=ebnm) package, so sparsity in the
loadings and shrinkage of the survival coefficients are learned from the data
through empirical Bayes priors rather than tuned penalties. The non-conjugate Cox likelihood is handled
with a second-order Taylor approximation, and the β update accounts for
uncertainty in the program scores (E[l²] rather than l̄²).

## Current results

The recommended single-modality configuration is YFB with per-platform
z-standardized expression and DeSurv-style gene selection (genes ranked by mean
and variance within each platform; 2,064 genes shared across cohorts). Trained on
273 patients (TCGA-PAAD and CPTAC-PDAC) with K = 7 initial programs, the model
retains two prognostic programs, one associated with worse and one with better
survival, and reaches a mean concordance of 0.627 across five held-out PDAC
cohorts (n = 616).

**Multimodal extension (in progress).** A joint YFB model for matched gene
expression and DNA methylation (TCGA training, ICGC validation) is implemented in
`code/fit_multimodal_yfb.R` with prediction, simulation and K-selection
utilities. Its derivation is in
[`derivations/multimodal_YFB/`](derivations/multimodal_YFB/multimodal_YFB_derivation.pdf).
It is not yet the production model.

Development history, benchmark tables and the reasoning behind each modelling
choice are in [`PROJECT_STATUS.qmd`](PROJECT_STATUS.qmd) and
[`DECISIONS.md`](DECISIONS.md). Advisor-meeting updates are collected in the
[progress notebook](https://andrew-walther.github.io/multiomicsGEP/).

## Getting started

### Requirements

R (developed with R 4.5) and

```r
install.packages(c("survival", "ebnm", "yaml"))
```

Scripts use paths relative to the repository root; run them from there.

### Minimal example

```r
# The script ends with a stand-alone runner that stops when no data are set;
# the function definitions are already loaded by then.
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) NULL)
source("code/predict_cox_on_yf.R")

set.seed(1)
n <- 150; p <- 200
L <- matrix(rexp(n * 3), n, 3)
F <- matrix(rexp(p * 3) * rbinom(p * 3, 1, 0.2), p, 3)
Y <- scale(L %*% t(F) + matrix(rnorm(n * p), n, p))
time   <- rexp(n, exp(drop(L %*% c(0.8, -0.8, 0))))
status <- rbinom(n, 1, 0.7)

train <- 1:100
fit <- fit_cox_on_yf(Y[train, ], time[train], status[train], K = 5, verbose = FALSE)

fit$EBeta   # posterior mean survival coefficients (K); unused programs shrink to 0
fit$EF      # p × K posterior mean loadings

# Score held-out patients by projection
risk <- predict_cox_on_yf(Y[-train, ], fit$EF, fit$EBeta, fit$EF_norms)$risk_scores
```

### Benchmarks

```bash
Rscript results/benchmark_sim/run_YFB_benchmark.R   # YFB (recommended)
Rscript results/benchmark_sim/run_LB_benchmark.R    # LB
```

Outputs are written to `results/benchmark_sim/outputs/`. Hyperparameters and
simulation settings are read from `config/globals.yml`.

### Tests

```bash
Rscript tests/run_tests.R
```

### Data

Every cohort is bulk tumour gene expression, measured by RNA-seq or microarray.
Patient data are not distributed with this repository. Scripts read the PDAC
cohorts from the directory in `PDAC_DATA_ROOT`. Sample sizes and censoring are
for the patients used in the analysis.

| Cohort | Expression platform | n | Censored | Role |
|--------|---------------------|---|----------|------|
| TCGA-PAAD | RNA-seq | 144 | 48% | Training |
| CPTAC-PDAC | RNA-seq | 129 | 50% | Training |
| Dijk | RNA-seq | 90 | 10% | External validation |
| Moffitt | Microarray | 123 | 33% | External validation |
| PACA-AU (array) | Microarray | 63 | 40% | External validation |
| PACA-AU (RNA-seq) | RNA-seq | 52 | 40% | External validation |
| Puleo | Microarray | 288 | 37% | External validation |

The multimodal analysis uses matched TCGA and ICGC expression and methylation,
assembled by `code/load_multiomics_data.R`.

## Repository structure

```
multiomicsGEP/
├── code/            model implementation (see code/README.md)
├── config/          globals.yml: hyperparameters and simulation settings
├── tests/           test suite (run_tests.R)
├── results/         benchmark and real-data analysis runners and outputs
├── derivations/     variational update derivations
├── docs/            update walkthroughs, progress notebook, plans and reports
├── paper/           write-up (see paper/writeup/README.md) and manuscript sources
├── presentation/    lab-meeting slides
├── longleaf_setup/  SLURM scripts for the UNC Longleaf cluster
├── PROJECT_STATUS.qmd   full project documentation and development log
├── DECISIONS.md         dated modelling and analysis decisions
└── ROADMAP.md           priorities and completed milestones
```

## Documentation

- [Corrected CAVI derivations](derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf)
- [Multimodal YFB derivation](derivations/multimodal_YFB/multimodal_YFB_derivation.pdf)
- [Project write-up](paper/writeup/writeup.pdf): methods, evaluation, results and full derivations
- [Progress notebook](https://andrew-walther.github.io/multiomicsGEP/)

## License

MIT (see [`LICENSE`](LICENSE)). The PDAC, TCGA and ICGC data are not included and
remain under their own terms of use.

## Author

Andrew Walther
