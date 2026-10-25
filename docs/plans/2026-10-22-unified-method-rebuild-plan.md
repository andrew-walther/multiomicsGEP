# Unified SSBMF method and R package: first-principles rebuild plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Andrew approves every **Gate** before work past it starts. Tasks marked **[unattended-safe]** may run in an overnight session once the gate before them is approved; every other task needs Andrew in the loop.

**Goal:** An installable R package, `ssbmf`, that fits one survival-supervised Bayesian matrix factorization model to one or more data blocks (expression only, or expression + methylation) and returns programs, loadings, survival coefficients with uncertainty, and risk predictions for new patients. Its fitting algorithm must provably improve a stated objective at every step, converge in simulation and on the PDAC cohorts, recover planted programs (including a low-variance prognostic one), and predict survival at least as well as the current benchmark.

**Architecture:** Start from a written model specification, derive every update from it, and implement the derivation as an R package in `ssbmf/` inside this repository, following the layout of DeSurv (`rashidlab/DeSurv`): a data constructor, a fitter returning an S3 object with `predict()`/`coef()` methods, a cross-validation helper, feature filtering, top-feature extraction, testthat tests and vignettes. Single-modality is the one-block case of the same `ssbmf_fit()`. This repository becomes the analysis repository (as `rashidlab/DeSurv-paper` is for DeSurv): its runners call the package, and the existing fitters stay untouched as references and to reproduce the frozen benchmark. Once stable, the package moves to its own repository.

**Tech stack:** R (≥ 4.5); package imports `ebnm` (1.1.42), `survival` (3.8.6), `stats`; suggests `flashier` (1.0.7; reference fits in tests), `testthat` (≥ 3.0.0), `knitr`, `rmarkdown`; development with `roxygen2`, `devtools`; Quarto for derivations; Longleaf SLURM for simulation panels. `Rcpp` only if Task 10 timings require compiled inner loops.

**Spec:** the *Requirements* and *Package interface* sections below (from the 10/16 advisor meeting and Andrew's 10/22 feedback), plus the derivation produced by Tasks 1–3 (`derivations/unified_ssbmf/`), which becomes the controlling specification for Tasks 4–8.

**Relation to earlier plans.** This plan sets the order of method work from 2026-10-22.
- `2026-10-16-round2-closeout-and-docs-plan.md`: Task 1 (freeze the round2 branch under a tag) is a prerequisite here (Task 0). The documentation restructure can proceed independently.
- `2026-10-16-preprocessing-alignment-and-unified-model-plan.md`: Task A1 becomes Task 13 here. Tasks A2–A4 are replaced by Tasks 4–8 and 15–17.
- `2026-10-16-feature-selection-experiment-plan.md`: unchanged in design, deferred to Part 7 (run with the package).
- Round2 remedies not yet tried (lower bound on the β inclusion probability, survival weighting, block weighting, different scale fixing) are candidates in Tasks 2, 3 and 11–12, not separate experiments.

---

## Requirements

**What the method must do**
1. Accept one or more data blocks for the same training patients, with survival (time, status), and optionally a cohort label used as a Cox stratum.
2. Learn K latent programs: patient activities L (shared across blocks) and block-specific feature loadings F_m.
3. Score new patients in closed form from their molecular data alone (projection, η = Σ_m Ỹ_m F_m β), with no refitting.
4. Let survival inform the loadings F_m, so that a prognostic program that explains little variance can be found.
5. Learn as much as possible from the data (priors fitted by empirical Bayes) and avoid tuning parameters. A weight between the survival and reconstruction terms is considered only after testing whether the model can set that balance itself (Part 4).
6. Given expression alone, behave as the single-modality model: there is only one implementation.
7. Report posterior uncertainty for the survival coefficients β, not only point estimates.
8. Classify every fitted program by role (survival-active, reconstruction-only, inactive) on one fitted model, without deleting and refitting.
9. Be usable by someone outside the project: installable package, documented functions, a quickstart that runs in seconds on simulated data, and vignettes for one- and two-block use.

**What the algorithm must do**
- Every coordinate update must not decrease one stated objective (the ELBO, or a stated lower bound on it). This is tested, not assumed.
- Prior refits, rescalings and Cox approximations are part of that objective's bookkeeping; any step outside it is removed or proven not to decrease the objective.
- Numerics are stable at realistic dimensions (2,000–18,000 genes, 10,000–300,000 CpGs): log-space risk-set sums and no unbounded growth of η with the number of features.

**Operational criteria (proposed; confirm at Gate 0).** A configuration is operational when all of these hold:

| # | Criterion | Proposed threshold |
|---|-----------|--------------------|
| O1 | Monotone objective | No coordinate step lowers the objective by more than 1e-8 × \|objective\| in any test or panel fit |
| O2 | Convergence | Relative objective change < 1e-6 within 1,000 sweeps in ≥ 95% of simulation fits and in every real-data fit at the selected K |
| O3 | High-variance program recovery (simulation) | Best-match \|r\| ≥ 0.9 in ≥ 90% of fits at or above the true rank |
| O4 | Low-variance prognostic recovery (simulation, lower noise) | Matched (\|r\| ≥ 0.7), survival-active, correct sign in ≥ 80% of datasets, and more often than the two-step comparator |
| O5 | Expression-only real data | Mean external C over the five cohorts no lower than 0.627 − 0.01, risk direction fixed in training |
| O6 | Multimodal real data | Joint and expression-only fits both converge; ICGC C reported with a paired, donor-level bootstrap interval for the difference |
| O7 | Coefficient uncertainty | 90% posterior intervals for β cover the true value in 80–95% of simulation fits |
| O8 | Package | `R CMD check` with no errors or warnings; quickstart and vignettes run |

## Package interface (modeled on DeSurv and flashier)

Exported functions and methods. Names and return fields are fixed here so the tasks agree; arguments may gain defaults during the derivation gates.

| Export | Role | DeSurv analogue |
|--------|------|-----------------|
| `ssbmf_data(X, time, status, strata = NULL, orientation = c("features_by_samples", "samples_by_features"))` | Validate and package inputs. `X` is a matrix (one block) or a named list of matrices with matching sample names; returns class `ssbmf_data` | `desurv_data()` |
| `ssbmf_control(...)` | All algorithm settings with documented defaults (`max_sweeps`, `tol`, prior families, β prior, `weight_rule`, `init`, `seed`, `verbose`) | arguments of `desurv_fit()` |
| `ssbmf_fit(data, K, control = ssbmf_control())` | Fit; returns class `ssbmf_fit` | `desurv_fit()` |
| `predict(fit, newdata = NULL, type = c("risk", "scores"))` | `"risk"`: η (log hazard; higher = worse) for `newdata` (an `ssbmf_data` or matrix/list aligned by feature name) or the training data; `"scores"`: per-program projection scores Z | `predict.desurv_fit()` |
| `coef(fit)` | data.frame of β: `program`, `mean`, `sd`, `lower`, `upper`, `lfsr` | `coef.desurv_fit()` |
| `ssbmf_loadings(fit, block = NULL)`, `ssbmf_activities(fit)` | Posterior-mean F_m (features × K) and L (samples × K), canonically scaled for reporting | fields of the fit |
| `ssbmf_classify(fit)` | Program roles with `pve` by block, β summary | — |
| `ssbmf_top_features(fit, n = 50, block = NULL)` | Top-loading features per program | `desurv_get_top_genes()` |
| `ssbmf_cv(data, K_grid, nfolds = 5, control = ssbmf_control(), seed = 1)` | Training cross-validated C over `K_grid` (and weight settings if Part 4 adopts any), one-standard-error rule, refit on all training data; returns `list(cv = data.frame, best, fit)` | `desurv_cv()` |
| `ssbmf_filter_features(X, method = "mean_variance", n_top)` and `ssbmf_preprocess(raw, spec)` | Feature screening and the approved preprocessing (Task 15) | `gene_filter()`, `preprocess_data()` |
| `ssbmf_simulate(...)` | Simulated data for examples, vignettes and tests | — |
| `print()`/`summary()` for `ssbmf_data`, `ssbmf_fit` | Dimensions, convergence, roles, β table | `print.desurv_data()` |

**Fitted object (`ssbmf_fit`) fields:** `L` (list: `mean`, `mean2`), `F` (named list by block, each `mean`, `mean2`), `beta` (`mean`, `mean2`), `mu`, `tau` (by block), `center` (training feature means by block), `priors` (fitted g for L, F_m, β), `objective` (trace by sweep and its components), `step_decreases`, `converged`, `sweeps`, `roles`, `features` (names by block), `control`, `call`, `package_version`.

**Quickstart (target for the README and `?ssbmf_fit`):**
```r
library(ssbmf)
sim  <- ssbmf_simulate(n = 120, p = c(expression = 200, methylation = 400), K = 3, seed = 1)
dd   <- ssbmf_data(sim$X, sim$time, sim$status)
fit  <- ssbmf_fit(dd, K = 5)
coef(fit)                                   # beta with intervals
ssbmf_classify(fit)                         # program roles
risk <- predict(fit, newdata = sim$validation$X, type = "risk")
```

## Global constraints

- Do not modify `results/benchmark_snapshot_2026-10-14/`, the tag `benchmark-interim-writeup-2026-10-14`, the existing `code/` fitters or `tests/run_tests.R`; it keeps passing unchanged.
- The package never sources files from `code/`. Analysis runners load it with `pkgload::load_all("ssbmf")` locally or `library(ssbmf)` after `R CMD INSTALL ssbmf` (Longleaf).
- Package defaults live in `ssbmf_control()` and are documented there. Project analysis settings (cohorts, K grids, thresholds for runners) live in a new `ssbmf:` section of `config/globals.yml` and are passed in by runners.
- Package code reads no data files, writes nothing, and contains no patient data. Never commit patient data or fitted `.rds` objects derived from it; real-data fits go to git-ignored `results/ssbmf/outputs/`.
- External cohorts never choose K, weights, priors, preprocessing or stopping rules. Each candidate configuration is scored on external data once, after it is frozen on training data.
- Risk is η = log hazard (higher = worse) by construction; no direction is chosen on evaluation data (never `oriented_cindex()`).
- Roxygen documentation for every function (exported functions with `@examples` that run in under 5 s), math notation in comments, no assistant-tool names in code or deliverables. Author in DESCRIPTION: Andrew Walther (aut, cre); co-authors decided at Gate 0.
- Attempt budget: if the same failure survives two fix attempts, stop and report. Each Part ends with keep / revise / drop, recorded in `DECISIONS.md`.

## Review focus

1. **Many correlated features.** At 25k+ CpGs the old fit's risk scores reached SD ≈ 39,000. Expected: η stays bounded as the number of features grows. Test in Task 7 (`fit is scale-stable as p grows`).
2. **All β near zero early in fitting.** The old fits fell into an absorbing β = 0 state. Expected: a program's β can become nonzero again. Test in Task 6 (`beta update leaves the zero state`).
3. **Tied event times and strata with few events.** Expected: the Cox term matches `survival::coxph` (Breslow) to 1e-8. Test in Task 4.
4. **Constant or all-missing features, or mismatched sample names across blocks.** Expected: `ssbmf_data()` stops with an informative error. Test in Task 4.
5. **New data with features in a different order, extra features or missing features, or the other orientation.** Expected: prediction aligns by feature name, ignores extras, errors on missing names. Test in Task 8.

---

## File structure

| Path | Responsibility |
|------|----------------|
| `derivations/unified_ssbmf/ssbmf_derivation.qmd` | Model specification, objective, every update, Cox treatment, scale, K, weights (Tasks 1–3, 11) |
| `derivations/unified_ssbmf/verify_derivation.R` | Numerical checks of the derivation: bounds, stationarity of each update by finite differences |
| `ssbmf/DESCRIPTION`, `ssbmf/NAMESPACE` (roxygen), `ssbmf/LICENSE` (MIT), `ssbmf/.Rbuildignore`, `ssbmf/README.Rmd` | Package metadata and front page |
| `ssbmf/R/data.R` | `ssbmf_data()`, validation, orientation, print |
| `ssbmf/R/cox.R` | Log-space risk-set sums, partial log-likelihood, the chosen Cox surrogate (internal) |
| `ssbmf/R/priors.R` | `ebnm` wrapper returning posterior moments and the KL against the prior actually used (internal) |
| `ssbmf/R/updates.R`, `ssbmf/R/objective.R` | Coordinate updates; the full objective and its components (internal) |
| `ssbmf/R/init.R`, `ssbmf/R/control.R`, `ssbmf/R/fit.R` | Initialization; `ssbmf_control()`; `ssbmf_fit()` |
| `ssbmf/R/methods.R` | `predict`, `coef`, `print`, `summary`, `ssbmf_loadings`, `ssbmf_activities`, `ssbmf_top_features` |
| `ssbmf/R/classify.R`, `ssbmf/R/cv.R` | `ssbmf_classify()`; `ssbmf_cv()` |
| `ssbmf/R/simulate.R`, `ssbmf/R/preprocess.R` | `ssbmf_simulate()`; `ssbmf_filter_features()`, `ssbmf_preprocess()` |
| `ssbmf/tests/testthat/test-*.R` | Package tests (testthat edition 3) |
| `ssbmf/vignettes/quickstart.Rmd`, `ssbmf/vignettes/multiblock.Rmd` | One-block and two-block walkthroughs on simulated data |
| `results/ssbmf/` | Panel and real-data runners (analysis repository side); compact committed summaries; ignored fits |

---

## Part 0: Set-up and decisions

### Task 0: Freeze, branch, scaffold, and confirm decisions (Gate 0)

**Files:** Create the package skeleton (`ssbmf/DESCRIPTION`, `NAMESPACE`, `LICENSE`, `.Rbuildignore`, `R/ssbmf-package.R`, `tests/testthat.R`); modify `docs/plans/README.md` (add this plan; mark the 10/16 ordering superseded); modify `DECISIONS.md` (Gate 0 answers); modify `.gitignore` if build artifacts appear (`ssbmf/*.tar.gz`, `ssbmf.Rcheck/`).

- [ ] **Step 1:** Carry out Task 1 of the round2 close-out plan (tag `archive/round2-supervised-2026-10`; verify `git show <tag>:docs/plans/Round2_attempts_log.md` resolves).
- [ ] **Step 2:** `git switch main && git switch -c rebuild/unified-ssbmf`.
- [ ] **Step 3:** Install development tools (`testthat`, `devtools` or `pkgload` + `rcmdcheck`) after Andrew approves D8.
- [ ] **Step 4:** Scaffold the package; `R CMD build ssbmf && R CMD check ssbmf_*.tar.gz` passes on the empty package. Commit.
- [ ] **Step 5:** Andrew answers the decisions below; record them in `DECISIONS.md`. **Gate 0.**

**Decisions for Gate 0**

| ID | Question | Recommendation |
|----|----------|----------------|
| D1 | Fresh implementation or repair the old fitters? | Fresh, as the `ssbmf` package; old code kept as reference. The old fitters carry ~30 option flags each from successive patches; a derivation-first rebuild is easier to verify than auditing every flag combination |
| D2 | Operational criteria O1–O8 and thresholds | As tabled; adjust the numbers now, not after results are seen |
| D3 | Survival predictor | Keep the projection form η = Σ_m Ỹ_m F_m β (closed-form scoring of new patients) |
| D4 | Cohorts in training | Cox strata by cohort only; no cohort columns in L |
| D5 | Patients missing a block at prediction time | Out of scope for version 1 |
| D6 | Comparators | Two-step (flashier then Cox) in all panels; DeSurv on real data only after O1–O5 hold, with renewed authorization |
| D7 | Compute | Unit tests local; panels on Longleaf (SLURM array, one fit per task; memory confirmed in Task 10) |
| D8 | New development dependencies | `testthat` (≥ 3.0.0) and `devtools`/`pkgload` + `rcmdcheck`, as in DeSurv; the old repository tests keep their own harness |
| D9 | Package name and authorship | `ssbmf`; DESCRIPTION authors per the manuscript list (Walther, Young, Liu, Rashid) or Andrew alone until release |
| D10 | Input orientation | Accept both via `orientation`; default `"features_by_samples"` to match DeSurv and Bioconductor convention; internally samples × features |
| D11 | Where the package lives | `ssbmf/` in this repository until Gate 8, then its own repository `rashidlab/ssbmf`; this repository becomes the analysis/paper repository |

---

## Part 1: Model specification and derivations (Andrew in the loop)

### Task 1: Write the model specification (Gate 1a)

**Files:** Create `derivations/unified_ssbmf/ssbmf_derivation.qmd` §1–§2.

- [ ] **Step 1: Generative model.** For blocks m = 1…M:
  Y_m = 1 μ_mᵀ + L F_mᵀ + E_m, (E_m)_ij ~ N(0, 1/τ_mj); L_ik ~ g_L,k; (F_m)_jk ~ g_m,k; β_k ~ g_β; Cox partial likelihood with η_i = Σ_k β_k Z_ik, Z_ik = Σ_m (Ỹ_m F_m)_ik, Ỹ_m = Y_m − 1 ȳ_mᵀ (training means, fixed). Include w_m (block weights) and w_S (survival weight) as explicit symbols, both 1 by default, so Part 4 can change them without re-deriving.
- [ ] **Step 2: Choices this section settles, each with a stated reason.**
  - Feature intercept μ_m in every block (removes the reason the old single-modality model needed z-scoring; lets nonnegative loadings model log-scale data).
  - Loading prior families: point-exponential on L; for F_m, point-exponential (nonnegative, interpretable programs) vs point-Laplace (signed). State what each implies for below-average values. Confirmed empirically in Task 10.
  - Identifiability: no rescaling inside the fitting loop (as in flashier; empirical-Bayes priors absorb scale); a canonical scaling applied only to reported output (`ssbmf_loadings()`, `ssbmf_activities()`, `coef()`), shown to leave L F_mᵀ and η unchanged.
  - Meaning of K: K_init is an upper bound; point-mass priors may zero a program during fitting; roles are assigned after fitting by posterior quantities (β: local false sign rate or interval excluding 0; reconstruction: proportion of variance explained).
- [ ] **Step 3: M = 1 reduction.** Show the M = 1 case is the single-modality model; list every difference from the frozen D4 configuration (intercept, no z-scoring, no in-loop canonicalization, β prior). These are why O5 is a non-inferiority check, not an exact reproduction.
- [ ] **Step 4:** Render the PDF; Andrew reviews. **Gate 1a.**

### Task 2: Objective and treatment of the Cox term (Gate 1b)

**Files:** `ssbmf_derivation.qmd` §3; `verify_derivation.R` (Cox checks).

The round2 audit found that refreshing the Cox working approximation lowered the objective in 18 of 50 refreshes. A second-order Taylor expansion at the current η̄ is not a lower bound, so coordinate ascent on it is not guaranteed to increase anything. This task chooses a treatment that makes the objective well defined.

- [ ] **Step 1:** Write the full objective: E_q log p(Y | ·) over blocks, the Cox term, and KL(q ‖ g) for L, F_m, β (τ, μ point-estimated). Show where E_q[η] and Var_q[η] enter, given that η is a product of β and F.
- [ ] **Step 2:** Derive three treatments of E_q log PL and state for each whether coordinate ascent is guaranteed monotone: (a) the current working quadratic, as reference; (b) a fixed-curvature quadratic lower bound on each risk-set log-sum-exp (Böhning-type: the Hessian of log Σ exp is ≼ ½ I), giving a minorize–maximize scheme; (c) a Jensen/moment bound E log Σ_j exp η_j ≤ log Σ_j E exp η_j under a Gaussian approximation to q(η_j).
- [ ] **Step 3:** In `verify_derivation.R`, check numerically on random states (n = 50, ties included) that (b) and (c) lie below the Monte Carlo value of E_q log PL; report the gap.
- [ ] **Step 4:** Recommend one treatment with the reason; Andrew decides. **Gate 1b.**

### Task 3: Coordinate updates (Gate 1c)

**Files:** `ssbmf_derivation.qmd` §4–§6; `verify_derivation.R` (stationarity checks).

- [ ] **Step 1:** Derive each update as the maximizer of the Gate 1b objective with everything else fixed: μ_m and τ_mj (closed form; τ floored); q(L_·k) via `ebnm` on pseudo-observations pooled over blocks; q(F_m,jk) feature by feature with the coupling term through η (the term the old α_F update dropped); q(β_k); each prior refit with the KL recomputed against the refitted prior.
- [ ] **Step 2:** Absorbing β = 0: show which prior choices make it absorbing (per-factor Normal with σ̂² = max(0, x² − s²); shared point-Laplace with π̂ → 0) and specify a default that is not (candidates: π_β ≥ 1/K; Normal with a weak inverse-gamma hyperprior on σ²_β). Derive the chosen update.
- [ ] **Step 3:** Initialization: SVD-based start; optional survival-score start for one column (F_·K ∝ max(0, ±Ỹᵀu⁰), u⁰ the null-model Cox score residuals), since the survival contribution to F is proportional to the current β.
- [ ] **Step 4:** In `verify_derivation.R`, check by central finite differences that the objective's gradient in each updated block of variational parameters is ≤ 1e-6 after its update.
- [ ] **Step 5:** Numbered pseudocode for the full algorithm. Andrew reviews. **Gate 1c.** After this gate, Tasks 4–10 are unattended-safe.

---

## Part 2: Package implementation (TDD) — [unattended-safe after Gate 1c]

Each task: write failing tests in `ssbmf/tests/testthat/test-<name>.R`; run `Rscript -e 'devtools::test("ssbmf", filter = "<name>")'` (fail); implement with roxygen docs; run again (pass); run `Rscript -e 'devtools::document("ssbmf"); devtools::test("ssbmf")'` (all pass) and `Rscript tests/run_tests.R` (old suite unchanged); commit.

### Task 4: Data constructor and Cox primitives

**Files:** Create `ssbmf/R/data.R`, `ssbmf/R/cox.R`, `ssbmf/tests/testthat/test-data.R`, `test-cox.R`.

**Interfaces — Produces:**
- `ssbmf_data(X, time, status, strata = NULL, orientation)` → class `ssbmf_data`: `Y` (named list of samples × features matrices; a single matrix becomes block `"block1"`), `time`, `status`, `strata` (factor or NULL), `n`, `p` (named integer), `samples`, `features` (named list).
- Internal: `.logsumexp_risk(eta, time, strata)` → numeric vector of log Σ_{j∈R_i} exp η_j (Breslow); `.cox_loglik(eta, time, status, strata = NULL)` → scalar; `.cox_surrogate(...)` with the signature fixed by Gate 1b.

- [ ] **Step 1: Failing tests.**
```r
test_that("cox loglik matches coxph (Breslow, ties, strata)", {
  set.seed(1); n <- 60; x <- rnorm(n); eta <- 0.7 * x
  time <- round(rexp(n), 1); status <- rbinom(n, 1, 0.7); s <- rep(1:2, 30)
  ref <- survival::coxph(survival::Surv(time, status) ~ x + strata(s), ties = "breslow",
           init = 0.7, control = survival::coxph.control(iter.max = 0))$loglik[1]
  expect_equal(ssbmf:::.cox_loglik(eta, time, status, s), ref, tolerance = 1e-8)
})
test_that("log-space risk sums are finite for |eta| up to 1e4", {
  expect_true(all(is.finite(ssbmf:::.logsumexp_risk(c(1e4, -1e4, 0, 5e3), 4:1, NULL))))
})
test_that("constant feature is rejected", {
  X <- matrix(rnorm(40), 4, 10, dimnames = list(paste0("g", 1:4), paste0("s", 1:10)))
  X[2, ] <- 3
  expect_error(ssbmf_data(X, rexp(10), rep(1, 10)), "constant")
})
test_that("block sample names must match", {
  A <- matrix(rnorm(40), 4, 10, dimnames = list(paste0("g", 1:4), paste0("s", 1:10)))
  B <- matrix(rnorm(60), 6, 10, dimnames = list(paste0("c", 1:6), paste0("t", 1:10)))
  expect_error(ssbmf_data(list(expression = A, methylation = B), rexp(10), rep(1, 10)), "sample")
})
test_that("matrix input becomes one block; orientation is honored", {
  X <- matrix(rnorm(40), 4, 10, dimnames = list(paste0("g", 1:4), paste0("s", 1:10)))
  dd <- ssbmf_data(X, rexp(10), rep(1, 10))
  expect_identical(names(dd$Y), "block1"); expect_identical(dim(dd$Y$block1), c(10L, 4L))
  dd2 <- ssbmf_data(t(X), rexp(10), rep(1, 10), orientation = "samples_by_features")
  expect_identical(dd$Y, dd2$Y)
})
```
- [ ] **Steps 2–5:** run (fail), implement, run (pass), commit `Add ssbmf_data() and log-space Cox primitives`.

### Task 5: Prior fitting with exact KL

**Files:** Create `ssbmf/R/priors.R`, `ssbmf/tests/testthat/test-priors.R`.

**Interfaces — Produces (internal):** `.ebnm_fit(x, s, family, g_fixed = NULL, pi_min = 0)` → list(`mean`, `mean2`, `kl` (summed KL(q ‖ g) against the returned `g`), `g`, `loglik`). `family` ∈ {"point_exponential", "point_laplace", "normal", "point_normal"}.

- [ ] **Step 1: Failing tests.** `kl` equals the closed-form Normal–Normal KL for `family = "normal"` (tolerance 1e-8); `kl ≥ 0` and `mean2 ≥ mean^2` for every family on random inputs; EBMF identity `loglik == sum(E_q log N(x; θ, s²)) - kl` (1e-6); with `pi_min = 1/7` and all x ≈ 0 the fitted slab weight is ≥ 1/7.
- [ ] **Steps 2–5:** run, implement (call `ebnm::ebnm()`; `pi_min` per the Task 3 derivation), run, commit.

### Task 6: Updates and the objective

**Files:** Create `ssbmf/R/updates.R`, `ssbmf/R/objective.R`, `ssbmf/tests/testthat/test-updates.R`.

**Interfaces — Produces (internal):** state list `st` (`mu`, `EL`, `EL2`, `EF`, `EF2` by block, `Ebeta`, `Ebeta2`, `tau`, `gL`, `gF`, `gbeta`, `kl`); `.update_mu(st, d)`, `.update_tau(st, d)`, `.update_L(st, d, k)`, `.update_F(st, d, m, k)`, `.update_beta(st, d, k)` each returning `st`; `.objective(st, d, w)` → list(`total`, `recon` (by block), `surv`, `kl_L`, `kl_F`, `kl_beta`), `w` = list(`block`, `surv`).

- [ ] **Step 1: Failing tests.**
  - Monotonicity: for 20 random states (n = 40; blocks of 30 and 60 features; K = 3), each update applied alone lowers `.objective()$total` by at most 1e-8 × |total|.
  - Stationarity: after `.update_F(st, d, m, k)`, the finite-difference gradient in the updated coordinates is ≤ 1e-5.
  - With `w$surv = 0`, the F and L updates equal the unsupervised EBMF updates (hand-computed one-feature case).
  - `beta update leaves the zero state`: from Ebeta = Ebeta2 = 0 on data with planted β = 1.5, after 5 sweeps |Ebeta_k| > 0.5 for the planted program.
- [ ] **Steps 2–5:** run, implement per the Gate 1c pseudocode, run, commit.

### Task 7: `ssbmf_control()` and `ssbmf_fit()`

**Files:** Create `ssbmf/R/init.R`, `ssbmf/R/control.R`, `ssbmf/R/fit.R`, `ssbmf/tests/testthat/test-fit.R`, `ssbmf/tests/testthat/helper-sim.R` (`sim_small(n, p, K, beta, seed)`: nonnegative programs, Weibull survival, features × samples; used by Tasks 7–8 until `ssbmf_simulate()` exists, then replaced by it in Task 9).

**Interfaces — Consumes:** Tasks 4–6. **Produces:** `ssbmf_control(...)`; `ssbmf_fit(data, K, control)` → `ssbmf_fit` with the fields listed under *Package interface*.

- [ ] **Step 1: Failing tests.**
  - On `sim_small(n = 120, p = c(expression = 50, methylation = 50), K = 3, beta = c(0.7, 0.35, 0), seed = 1)`: `fit$step_decreases == 0` and `fit$converged`.
  - With `control$w_surv = 0`, fitted L Fᵀ agrees with a `flashier` backfit from the same initialization and priors, intercept as a fixed all-ones factor: relative Frobenius difference < 1e-3 (if the conventions differ, record the measured difference and agree the tolerance at review). Skipped when flashier is not installed.
  - `fit is scale-stable as p grows`: same latent structure with a correlated block of 500 vs 5,000 features; sd of the fitted η differs by less than a factor of 2.
  - Matrix input and one-element list input give identical fits.
  - A non-`ssbmf_data` input errors with a message pointing to `ssbmf_data()`.
- [ ] **Steps 2–5:** run, implement, run, commit.

### Task 8: Methods, roles, top features

**Files:** Create `ssbmf/R/methods.R`, `ssbmf/R/classify.R`, `ssbmf/tests/testthat/test-methods.R`.

**Interfaces — Produces:** `predict.ssbmf_fit`, `coef.ssbmf_fit`, `print`/`summary` methods, `ssbmf_loadings()`, `ssbmf_activities()`, `ssbmf_top_features()`, `ssbmf_classify()` (data.frame `program`, `role` ∈ {"survival_active", "reconstruction_only", "inactive"}, `pve_<block>`, `beta_mean`, `beta_sd`, `lfsr`).

- [ ] **Step 1: Failing tests.**
  - Permuting features of `newdata`, or adding extra features, leaves `predict()` unchanged; a missing feature errors naming it; `newdata` in the other orientation gives the same result when `orientation` is passed.
  - `predict(fit)` (no `newdata`) reproduces the fitted η̄ (1e-10).
  - Direction: with one planted adverse program (β > 0), `frozen_reverse_cindex()`-equivalent C of `predict(fit, val, "risk")` on validation data > 0.6, computed with `survival::concordance(Surv(time, status) ~ risk, reverse = TRUE)`.
  - `coef(fit)` has one row per program, `lower ≤ mean ≤ upper`.
  - `ssbmf_classify()` labels exactly the planted prognostic program survival-active.
  - The canonical scaling in `ssbmf_loadings()`/`ssbmf_activities()`/`coef()` reproduces `predict()` exactly.
- [ ] **Steps 2–5:** run, implement, run, commit.

---

## Part 3: Simulation verification at w = 1 — [unattended-safe after Gate 1c]

### Task 9: Simulators and registered simulation designs

**Files:** Create `ssbmf/R/simulate.R`, `ssbmf/tests/testthat/test-simulate.R`, `results/ssbmf/sim_design.md`.

The existing simulations (50 + 50 features; 1,000 genes) converged where real data did not, so they do not reproduce the real-data difficulty.

**Interfaces — Produces:** `ssbmf_simulate(n, p, K, low_variance_prognostic = FALSE, feature_correlation = c("none", "block"), noise = c("lower", "higher"), n_validation = n, seed)` → list(`X` (features × samples, by block), `time`, `status`, `truth` (`L`, `F`, `beta`), `validation`). Semi-synthetic data from real expression stays in the analysis repository: `results/ssbmf/simulate_semisynthetic.R` (`simulate_semisynthetic(Y_real, K_true, beta, seed)`), because it needs patient data.

- [ ] **Step 1:** Tests: dimensions and orientation; censoring 0.30 ± 0.05; planted low-variance program explains < 2% of variance; same seed gives identical output. Replace `sim_small()` in Tasks 7–8 tests with `ssbmf_simulate()` and confirm they still pass.
- [ ] **Step 2:** Write `sim_design.md`: realistic (2,000 expression + 10,000 correlated methylation features) and semi-synthetic designs × noise × K_init ∈ {K_true, K_true + 2, 2 K_true} × 20 datasets; arms: `ssbmf_fit` at w = 1 for each F prior family from Task 1, and the two-step comparator; metrics O1–O4, O7; runtime and memory. Commit before any run.
- [ ] **Step 3:** Implement, test, commit.

### Task 10: Run the w = 1 panel (Gate 3)

**Files:** Create `results/ssbmf/run_sim_panel.R`, `results/ssbmf/run_sim_panel.slurm`, `results/ssbmf/summarize_sim_panel.R`; outputs `results/ssbmf/outputs/sim_panel/` (fits ignored; summary CSVs committed).

- [ ] **Step 1:** Install the package on Longleaf (`R CMD INSTALL ssbmf`); smoke-run 2 datasets per cell; record runtime and memory (D7; decide whether `Rcpp` inner loops are needed).
- [ ] **Step 2:** Submit the array; summarize against O1–O4, O7 per cell.
- [ ] **Step 3:** Report which criteria hold at w = 1. **Gate 3:** if O1–O3 fail, return to Tasks 2–3 (derivation), not to parameter changes.

---

## Part 4: Balancing survival against reconstruction, and block against block (Andrew in the loop)

### Task 11: Can the model set the balance itself? (Gate 4a)

**Files:** `ssbmf_derivation.qmd` §7; `results/ssbmf/balance_diagnostics.R`.

Under the Gaussian model each of the n × p_m entries counts as an independent observation, so reconstruction outweighs a survival term with n observations, and a block with more features outweighs one with fewer. If features are strongly correlated, that counting overstates the information.

- [ ] **Step 1: Measure.** On Task 10 fits and the real training data: survival share of the F-update precision; effective number of independent features per block (e.g., participation ratio of residual-covariance eigenvalues).
- [ ] **Step 2: Model-based options (no tuning parameter).** Derive the effective-feature weight w_m = p_eff,m / p_m from training residuals; examine whether a hyperprior on w_S is identifiable. A power-likelihood weight cannot generally be learned by maximizing the marginal likelihood; review calibration approaches (Holmes & Walker 2017, *Biometrika*; Bissiri, Holmes & Walker 2016, *JRSS-B*; Grünwald's SafeBayes) and record which apply. Verify each citation before relying on it.
- [ ] **Step 3: Fallback.** Training-CV choice of w_S on a short grid with the one-standard-error rule inside `ssbmf_cv()` (as DeSurv chooses α), stated as the last resort.
- [ ] **Step 4:** Andrew and advisor choose which options to test. **Gate 4a.**

### Task 12: Test the chosen balance options (Gate 4b)

**Files:** Modify `ssbmf/R/control.R`, `ssbmf/R/fit.R` (`control$weight_rule`); `ssbmf/tests/testthat/test-weights.R`; reuse the Task 10 panel.

- [ ] **Step 1:** Tests: `weight_rule = "none"` reproduces the w = 1 fit exactly; the effective-feature rule gives w_m ≈ 1 for independent features and < 1 for block-correlated features; O1 still holds.
- [ ] **Step 2:** Implement; rerun the panel for each chosen rule.
- [ ] **Step 3:** Compare on O3–O4 and against the two-step arm; keep/drop in `DECISIONS.md`. **Gate 4b.**

---

## Part 5: Data transformation and preprocessing (can start in parallel with Part 1)

### Task 13: Preprocessing specification (documentation)

**Files:** Create `docs/methods/preprocessing.md` (Task A1 of the 10/16 preprocessing plan, extended).

- [ ] **Step 1:** For each data type and cohort: source file, transformation, missing values, screening (rule, count, which patients it sees), scaling/centering, cross-cohort alignment, where implemented. Include the confirmed double log₂(x + 1) on TCGA, Dijk and PACA-AU RNA-seq (`load_pdac_raw()` then `preprocess_merged_cohorts()` with `cfg$pdac$platform_log_transform`); CPTAC and microarray cohorts are logged once.
- [ ] **Step 2:** Methylation: β-values vs M-values; imputation; probe filters (from the CpG inventory). Expression: nonnegative log scale with intercept vs per-cohort standardization.
- [ ] **Step 3:** Validation accounting: PACA-AU has 47 donors in both series plus one repeated RNA-seq donor (568 distinct people in 616 records); per-cohort C stays valid; pooled analyses use one record per donor.
- [ ] **Step 4:** Proposed pipeline, one choice per step, identical across blocks except stated exceptions; mark which steps belong in the package (`ssbmf_preprocess()`: transformations, screening, alignment) and which stay in the analysis repository (cohort loading, file paths). Andrew reviews.

### Task 14: Effect of the double log on the frozen benchmark (Gate 5a)

**Files:** Create `results/benchmark_sim/run_single_log_check.R`; outputs to a new folder (snapshot untouched).

- [ ] **Step 1:** Rerun D4 with the existing `fit_cox_on_yf()` (same seed, K = 7, α_F = 0) on single-logged data: training 5-fold CV C; overlap of the selected genes with the frozen 2,064; external C once per cohort with a donor-level paired bootstrap against the frozen risk scores.
- [ ] **Step 2:** Report the change; if material, update the writeup numbers and `PROJECT_STATUS.qmd` with Andrew's approval. **Gate 5a.**

### Task 15: Preprocessing in the package

**Files:** Create `ssbmf/R/preprocess.R`, `ssbmf/tests/testthat/test-preprocess.R`.

**Interfaces — Produces:** `ssbmf_filter_features(X, method = "mean_variance", n_top, orientation)` → character vector of retained features; `ssbmf_preprocess(raw, spec)` → list(`data` (`ssbmf_data` for training), `apply` (function applying the frozen training choices to a new cohort)). `spec` is a list documented in `?ssbmf_preprocess`; the analysis repository builds it from the approved Task 13 table in `config/globals.yml`.

- [ ] **Step 1:** Tests: log applied exactly once per platform flag (raw input → output equals log₂(x + 1) once); screening sees training samples only; `apply()` on a validation cohort uses training feature order and errors on missing features; `ssbmf_filter_features()` on a matrix with known means/variances returns the expected top features.
- [ ] **Steps 2–5:** run, implement, run, commit.

---

## Part 6: Real data (Andrew in the loop)

### Task 16: Expression only, TCGA + CPTAC (Gate 6a)

**Files:** Create `results/ssbmf/run_real_single.R`.

- [ ] **Step 1:** Training only: `ssbmf_cv()` over K_init ∈ {5, 7, 10, 15} with the approved configuration; report convergence (O2), roles, coefficient intervals.
- [ ] **Step 2:** Freeze; score the five external cohorts once; per-cohort C, mean C (O5), donor-aware pooled interval, comparison with D4 and the two-step comparator. **Gate 6a.**

### Task 17: Expression + methylation, TCGA → ICGC (Gate 6b)

**Files:** Create `results/ssbmf/run_real_multimodal.R`.

- [ ] **Step 1:** Training (TCGA, n = 144): joint and expression-only fits through the same `ssbmf_fit()`, K_init by `ssbmf_cv()`; convergence, roles, block shares per program.
- [ ] **Step 2:** Freeze; score ICGC primary (n = 50) once; paired donor-level bootstrap for joint vs expression-only (O6). **Gate 6b:** decide whether the new model replaces the benchmark (snapshot rule: at least as good on every goal, better on one). Update `PROJECT_STATUS.qmd`, `ROADMAP.md`, `DECISIONS.md`; add a progress-book chapter.

---

## Part 7: Experiments unlocked by an operational model

Start only after Gate 6a (expression experiments) or 6b (CpG and joint experiments); all use the package.
- Gene-count (2k → 17.9k), CpG-count (10k → 305k) and joint curves: where performance levels off, and how stable the programs are.
- Feature screening: unsupervised mean/variance vs + survival screen at matched sizes (`2026-10-16-feature-selection-experiment-plan.md`).
- Biology of the new programs (enrichment, CpG annotation to genes); DeSurv head-to-head if authorized.

---

## Part 8: Package release readiness (Gate 8)

### Task 18: Documentation, vignettes, check

**Files:** Create `ssbmf/README.Rmd` (rendered `README.md`), `ssbmf/vignettes/quickstart.Rmd`, `ssbmf/vignettes/multiblock.Rmd`; modify `ssbmf/DESCRIPTION` (version 0.1.0, `URL`, `BugReports`, `VignetteBuilder: knitr`); modify the project `README.md` quickstart to use the package.

- [ ] **Step 1:** README: what the model is (one paragraph with the model equations), installation (`remotes::install_github(...)` once the repository exists; `R CMD INSTALL ssbmf` before that), the quickstart above, outputs, citation.
- [ ] **Step 2:** Vignettes on simulated data only: one block (fit, roles, top features, prediction on validation data, `ssbmf_cv()` with a small grid) and two blocks (block shares per program, prediction).
- [ ] **Step 3:** `R CMD build ssbmf && R CMD check --as-cran ssbmf_0.1.0.tar.gz`: 0 errors, 0 warnings; notes listed and justified (O8).
- [ ] **Step 4:** Andrew decides whether to move the package to its own repository (D11), add a Zenodo DOI, and when to make it public. **Gate 8.**

---

## Open questions (for discussion before Gate 0)

1. **Timeline against the prelim.** Proposal due 11/16, exam 12/7. Parts 1–3 are the minimum for the proposal to describe a working supervised method; if they are not done by early November, the chapter describes the method and the benchmark and states supervision of F as the proposed aim.
2. **The paper's central claim** depends on Part 4: if survival supervision recovers low-variance programs that the two-step method misses, that is the claim; if not, the claim is a unified multi-block Bayesian factor–Cox model with closed-form scoring and coefficient uncertainty, distributed as a package.
3. **Runtime for the full CpG universe** (305k features) may need compiled inner loops (`Rcpp`, as DeSurv uses); decide after Task 10 timings.
4. **Package scope for version 1:** cross-validation and preprocessing helpers are included (as in DeSurv); Bayesian optimization over hyperparameters (DeSurv's `desurv_bayesopt()`) is not, because the method aims to need no tuning beyond K.
