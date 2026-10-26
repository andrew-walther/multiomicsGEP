# ssbmf

Survival-supervised Bayesian matrix factorization for one or more molecular
data blocks measured on the same patients.

For blocks $Y_m$ (samples by features),

$$Y_m = \mathbf 1\mu_m^\top + L F_m^\top + E_m,\qquad (E_m)_{ij}\sim N(0, 1/\tau_{mj}),$$

with patient activities $L$ shared across blocks and loadings $F_m$ specific
to each block. Survival follows a Cox model (optionally stratified, e.g. by
cohort) with linear predictor $\eta = \sum_m \tilde Y_m F_m \beta$, where
$\tilde Y_m$ is the block centered at its training means, so a new patient is
scored in closed form from their molecular data. Priors on $L$, $F_m$ and
$\beta$ are fitted by empirical Bayes, and the fit maximizes one stated
variational objective: every accepted step increases it, and the fit reports
any step that did not.

The model, the objective and every update are derived in
`derivations/unified_ssbmf/ssbmf_derivation.pdf` of the analysis repository.

## Installation

From the analysis repository root (the package is not yet on GitHub on its own):

```sh
R CMD INSTALL ssbmf
```

## Quickstart

```r
library(ssbmf)
sim  <- ssbmf_simulate(n = 120, p = c(expression = 200, methylation = 400), K = 3, seed = 1)
dd   <- ssbmf_data(sim$X, sim$time, sim$status)   # features x samples by default
fit  <- ssbmf_fit(dd, K = 5)                       # K is an upper bound
coef(fit)                                          # beta per SD of program score, with intervals
ssbmf_classify(fit)                                # survival-active / reconstruction-only / inactive
risk <- predict(fit, newdata = sim$validation$X, type = "risk")   # log hazard, higher = worse
```

`ssbmf_cv()` chooses K (and, if asked, the survival weight) by training
cross-validated concordance with the one-standard-error rule and refits on
all training data. `ssbmf_top_features()`, `ssbmf_loadings()` and
`ssbmf_activities()` return the programs on a canonical scale.

## Status

Development version (0.0.0.9000). Settings marked provisional in the
derivation are awaiting review; see `docs/findings/2026-10-27-rebuild-first-results.md`
in the analysis repository.
