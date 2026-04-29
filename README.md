# multiomicsGEP

**Supervised Bayesian Matrix Factorization for Joint Genomics and Survival Modelling**

---

## Overview

This project implements a **Supervised Bayesian Matrix Factorization** model that jointly decomposes high-dimensional genomics data — such as gene expression or DNA methylation matrices — while simultaneously modelling patient survival outcomes via a Cox proportional hazards model.

The core insight is that standard unsupervised factorization (e.g. PCA, NMF) finds latent structure that explains genomic variance, but has no reason to recover factors that are *clinically meaningful*. By supervising the factorization with survival data, this model discovers **Gene Expression Programs (GEPs)** that are both genomically coherent and prognostically relevant.

**The model in one equation:**

```
Y (n×p) = L (n×K) × F' (K×p) + E        ← genomics
h(t_i)  = h₀(t_i) exp(Lᵢ · β)           ← survival (Cox PH)
```

- **L** — patient loading matrix: where each patient sits in the latent space
- **F** — factor weight matrix: which genes/features define each program
- **β** — survival coefficients: how much each program contributes to prognosis
- **τ** — feature-specific noise precision

Inference is performed via **Coordinate Ascent Variational Inference (CAVI)**, where each update is solved as an **Empirical Bayes Normal Means (EBNM)** problem with point-normal priors — naturally promoting sparsity in both F and β.

---

## Repository Structure

```
multiomicsGEP/
│
├── README.md                          ← You are here
├── PROJECT_STATUS.qmd/.pdf           ← Full project documentation & session log
├── DECISIONS.md                       ← Architectural/analytical decision log
├── ROADMAP.md                         ← Prioritised next steps & completed items
├── CLAUDE.md                          ← Claude Code entry point (thin; defers to above)
│
├── config/
│   └── globals.yml                    ← Single source of truth for hyperparameters
│                                         (K_max, alpha_grid, tol, synthetic DGP params)
│
├── code/
│   ├── fit_modular.R                  ← ✅ Canonical CAVI loop (factor-wise Gauss-Seidel)
│   ├── update_beta.R                  ← Modular β update (scalar EBNM, Cox survival)
│   ├── update_L.R                     ← Modular L update (vector EBNM, dual-source)
│   ├── update_F.R                     ← Modular F update (vector EBNM, pure genomics)
│   ├── update_tau.R                   ← Modular τ update (closed-form MLE)
│   ├── compute_elbo.R                 ← Full ELBO: genomics + survival + KL divergences
│   ├── preprocess_desurv.R            ← DeSurv-aligned preprocessing (log2, top-2000, rank)
│   ├── select_alpha_cv.R              ← Alpha mixing CV selection via 1-SE rule
│   ├── predict.R                      ← Hold-out prediction (SVD pseudoinverse projection)
│   ├── train_test_split.R             ← Stratified 80/20 split preserving event rate
│   ├── feature_selection.R            ← Univariate Cox gene filtering (train-only)
│   ├── select_K.R                     ← K selection: auto_prune_K() + select_K_cv() stub
│   ├── Supervised_Bayesian_MF_V2.R    ← Monolithic reference implementation (V2, reference only)
│   ├── SupervisedMF_Context.md        ← AI/developer quick-reference for the code
│   └── legacy/                        ← Archived files (V1, early scripts)
│       ├── Supervised_Bayesian_MF.R   ← V1 original (archived, known issues)
│       ├── execute_update_beta.R      ← Early demo (superseded by demos/)
│       └── multiomicsGEP_code.Rmd     ← Early exploratory notebook
│
├── docs/                              ← Companion documentation (PDF + HTML)
│   ├── Makefile                       ← `make all` renders .qmd → .pdf + .html via Quarto
│   ├── fit_modular.qmd/.pdf/.html     ← fit_modular.R walkthrough (full CAVI loop)
│   ├── update_beta.qmd/.pdf/.html     ← β update: derivation, code, tests, demos
│   ├── update_L.qmd/.pdf/.html        ← L update: derivation, code, tests, demos
│   ├── update_F.qmd/.pdf/.html        ← F update: derivation, code, tests, demos
│   ├── update_tau.qmd/.pdf/.html      ← τ update: derivation, code, tests, demos
│   ├── PDAC_data_audit.qmd/.pdf/.html ← Audit of available PDAC cohorts & data quality
│   ├── SSMF_DeSurv_Sim_Benchmark.md  ← DeSurv benchmark design notes
│   └── update_L_fix.md               ← Debugging guide: A_surv/A_gen imbalance → β=0 fix
│
├── tests/                             ← 171 tests (run: Rscript tests/run_tests.R)
│   ├── run_tests.R                    ← Master test runner
│   ├── test_helpers.R                 ← Lightweight assertion framework (no testthat)
│   ├── test_update_beta.R             ← 24 tests for update_beta.R
│   ├── test_update_L.R                ← 28 tests for update_L.R
│   ├── test_update_F.R                ← 26 tests for update_F.R
│   ├── test_update_tau.R              ← 27 tests for update_tau.R
│   ├── test_predict.R                 ← 19 tests for predict.R + train_test_split.R
│   ├── test_elbo.R                    ← 15 tests for compute_elbo.R
│   ├── test_preprocess_desurv.R       ← Tests for preprocess_desurv.R
│   ├── test_select_alpha_cv.R         ← Tests for select_alpha_cv.R
│   └── test_real_data_loading.R       ← Real-data pipeline tests (auto-skip if data absent)
│
├── demos/                             ← Interactive demonstrations (5 per module)
│   ├── demo_update_beta.R
│   ├── demo_update_L.R
│   ├── demo_update_F.R
│   └── demo_update_tau.R
│
├── results/
│   ├── benchmark_sim/                 ← ✅ DeSurv benchmark (current, canonical)
│   │   ├── run_ssbmf_benchmark.R      ← Entry point: synthetic + PDAC cross-cohort
│   │   ├── ssbmf_summary_report.qmd/.pdf/.html  ← 24-page benchmark report
│   │   └── outputs/
│   │       ├── synthetic/
│   │       │   ├── point_normal/      ← tables/ + figures/ (11 figure types)
│   │       │   └── point_laplace/     ← tables/ + figures/
│   │       ├── real_data/
│   │       │   ├── tcga_only/{prior}/    ← primary single-cohort results
│   │       │   ├── cptac_only/{prior}/
│   │       │   └── merged/
│   │       │       ├── v2_{prior}/        ← v2 preprocessing (QN, intersect-first)
│   │       │       └── v2_lambda{X}_{prior}/ ← λ-sweep results
│   │       ├── diagnostic_heatmaps/       ← Phase 1 loading heatmaps (all modes × priors)
│   │       ├── ebmf_diagnostic/           ← EBMF Cox survival check (confirms data has signal)
│   │       └── ebmf_warmstart/            ← Warm-start experiments (β-only + full CAVI)
│   ├── modular_sim_factor/            ← Prior-family × K comparison (legacy, superseded)
│   │   ├── run_factor_modular_simulation.R
│   │   ├── run_prior_k_comparison.R
│   │   ├── synthetic/
│   │   └── PDAC/
│   ├── full_sim/                      ← V2 monolithic simulation (legacy)
│   │   ├── run_simulation.R
│   │   └── simulation_report.qmd/.pdf
│   ├── modular_sim_block/             ← Block-wise modular (deprecated)
│   ├── figures/                       ← Legacy per-cohort figures (modular_sim_factor era)
│   └── tables/                        ← Legacy per-cohort tables (modular_sim_factor era)
│
├── derivations/
│   ├── MF_UpdateDerivations/
│   │   ├── MF_Derivations_UpdateAlgo_REVISED.pdf  ← ✅ Corrected derivations (21 pages)
│   │   ├── MF_Derivations_UpdateAlgo_REVISED.tex
│   │   ├── MF_V2_Companion.pdf        ← ✅ Math ↔ code companion (17 pages)
│   │   ├── MF_V2_Companion.tex
│   │   └── MF_Derivations_UpdateAlgo_*.pdf  ← Historical drafts (contain errors R1–R8)
│   ├── qB/                            ← q(β) derivation (11 pages)
│   ├── qL/                            ← q(L) derivation (vector EBNM, dual-source)
│   ├── qF/                            ← q(F) derivation (τ cancellation property)
│   ├── qTau/                          ← q(τ) derivation (variance correction)
│   ├── EBMF/                          ← Empirical Bayes MF background theory
│   └── SurvivalMF/                    ← Survival + MF background notes
│
├── presentation/                      ← Lab meeting slide decks
│   └── walther_lab_meeting_04_09_2026/
│
├── longleaf_setup/                    ← UNC Longleaf HPC SLURM scripts
│   ├── README.md
│   ├── install_packages.R
│   └── run_*.sl                       ← SLURM job scripts
│
└── paper/
    ├── multiomicsGEP_manuscript.qmd   ← Manuscript draft (in progress)
    └── abstract_example.md
```

---

## Quickstart

### Prerequisites

```r
install.packages(c("survival", "ebnm"))
```

### Run the DeSurv Benchmark (current entry point)

```bash
# Synthetic validation + PDAC cross-cohort benchmark (both priors)
Rscript results/benchmark_sim/run_ssbmf_benchmark.R
```

This runs: (1) synthetic benchmark at n=300, p=1000, K_true=5 for both `point_normal` and `point_laplace` priors; (2) PDAC cross-cohort benchmark for TCGA-only and CPTAC-only training modes. Outputs go to `results/benchmark_sim/outputs/`. Render the summary report with:

```bash
quarto render results/benchmark_sim/ssbmf_summary_report.qmd
```

### Run the Real PDAC Analysis

The runner supports 7 PDAC cohorts via environment variable control:

```bash
# Single dataset
DATA_MODE=real DATASET_NAME=TCGA_PAAD \
  Rscript results/modular_sim_factor/run_factor_modular_simulation.R

# All 7 cohorts + pooled RNA-seq (with batch correction) + hold-out evaluation
DATA_MODE=real RUN_ALL=TRUE HOLDOUT_EVAL=TRUE \
  Rscript results/modular_sim_factor/run_factor_modular_simulation.R

# Longleaf HPC (override data path)
export PDAC_DATA_ROOT=/proj/rashidlab/data/PDAC
DATA_MODE=real RUN_ALL=TRUE \
  Rscript results/modular_sim_factor/run_factor_modular_simulation.R
```

**Advanced options** (all accept environment variable overrides):

| Argument | Default | Options | Description |
|----------|---------|---------|-------------|
| `prior_family` | `"point_normal"` | `"point_laplace"`, `"normal_scale_mixture"` | EBNM prior for L, F, β |
| `n_init` | `1` | any integer | Number of random restarts (best ELBO kept) |
| `init_method` | `"svd"` | `"random"` | Initialization strategy |
| `batch_correct` | `TRUE` | `FALSE` | limma batch correction for pooled data |
| `holdout_eval` | `FALSE` | `TRUE` | 80/20 stratified hold-out prediction |
| `feature_selection` | `"variance"` | `"cox"` | Gene selection method |
| `k_select` | `"fixed"` | `"auto_prune"` | K selection strategy |

```bash
# Example: 3 random inits with point_laplace prior + hold-out evaluation
DATA_MODE=real DATASET_NAME=CPTAC \
  PRIOR_FAMILY=point_laplace N_INIT=3 INIT_METHOD=random HOLDOUT_EVAL=TRUE \
  Rscript results/modular_sim_factor/run_factor_modular_simulation.R

# Example: K auto-pruning on Puleo array (fit K=10, prune to active)
DATA_MODE=real DATASET_NAME=Puleo_array K_SELECT=auto_prune \
  Rscript results/modular_sim_factor/run_factor_modular_simulation.R
```

**Available datasets:**

| Dataset | Platform | n | Censoring |
|---------|----------|---|-----------|
| TCGA_PAAD | RNA-seq | 144 | 48% |
| CPTAC | Proteomics | 129 | 50% |
| Dijk | RNA-seq | 90 | 10% |
| Moffitt_GEO_array | Microarray | 123 | 33% |
| PACA_AU_array | Microarray | 63 | 40% |
| PACA_AU_seq | RNA-seq | 52 | 40% |
| Puleo_array | Microarray | 288 | 37% |

**Note:** PDAC data files are stored locally (not in git). The default path is
`~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/UNC Dissertation (Liu)/PDAC_data`.
Override with `PDAC_DATA_ROOT` (e.g. `export PDAC_DATA_ROOT=/proj/rashidlab/data/PDAC` on Longleaf).

### Apply to Real Data

The recommended entry point for real data is `fit_supervised_mf_modular()` in `code/fit_modular.R`:

```r
source("code/fit_modular.R")   # also sources update_L/F/beta/tau.R automatically

res <- fit_supervised_mf_modular(
  Y      = your_matrix,    # numeric matrix: n patients × p genes (pre-normalised, column-centred)
  time   = your_time,      # numeric vector: survival/censoring time
  status = your_status,    # integer vector: 1 = event, 0 = censored
  K      = 5,              # number of latent factors (select via cross-validation)
  max_iter = 300,
  tol      = 1e-3,
  verbose  = TRUE
)

# Access results
res$EL     # n×K posterior mean patient loadings
res$EF     # p×K posterior mean factor weights (GEP signatures)
res$EBeta  # K posterior mean survival coefficients
res$EBeta2 # K posterior second moments (uncertainty)
res$Tau    # p noise precision per feature
res$history$rmse        # RMSE per iteration
res$history$elbo_proxy  # genomics ELBO per iteration
```

---

## What to Expect

On the simulated benchmark:

| Metric | Expected |
|--------|----------|
| Reconstruction RMSE | Converges near **1.0** (true noise SD) |
| ELBO proxy | **Non-decreasing** across iterations |
| β sign recovery | **(+, −, +, −, 0)** — matches ground truth |
| C-index (modular) | ~**0.86** on held-out tertiles |
| Factor 5 (β=0) | Correctly shrunk toward zero |

---

## Mathematical Background

The model is derived and documented in two companion documents (both in `derivations/MF_UpdateDerivations/`):

| Document | Description |
|----------|-------------|
| **`MF_Derivations_UpdateAlgo_REVISED.pdf`** | Full corrected CAVI derivations with all 8 errata (R1–R8) from the original working document resolved |
| **`MF_V2_Companion.pdf`** | 17-page companion that walks through every equation and maps it to the exact R variable in V2.R |

The key mathematical concepts are:

- **Mean-field CAVI:** The variational posterior factorises as q(L,F,β) = q(L)q(F)q(β), enabling coordinate-wise updates.
- **EBNM sub-problems:** Each coordinate update has the form: given pseudo-observations x with noise s, estimate a sparse prior g and return posterior moments. Solved via the `ebnm` R package.
- **Cox Taylor expansion:** The non-conjugate Cox likelihood is linearised via a 2nd-order Taylor expansion into a weighted Gaussian form, enabling EBNM updates for L and β.
- **Error-in-variables correction:** The β precision uses E[l²] (not l̄²), preventing survival coefficients from overfitting to uncertain loadings.
- **Gauss-Seidel updates:** Within each CAVI iteration, updates to earlier factors are immediately used when computing later ones — equivalent to block coordinate descent with sequential incorporation of new information.

---

## Version History

| Version | File | Status | Notes |
|---------|------|--------|-------|
| V1 | `code/legacy/Supervised_Bayesian_MF.R` | Archived | Original implementation; 6 known algorithmic issues |
| V2 | `code/Supervised_Bayesian_MF_V2.R` | Reference | Monolithic; all V1 issues corrected (A1–A6); kept for comparison |
| Modular | `code/fit_modular.R` + `update_*.R` | ✅ **Current** | Factor-wise Gauss-Seidel CAVI; tested (124/124); recommended for all new work |

**V2 improvements over V1:**

| ID | Fix |
|----|-----|
| A1 | True Gauss-Seidel CAVI (z_no_k from current EL/EBeta inside k-loop) |
| A2 | Orthogonalisation behind `orthogonalize=FALSE` flag (default off) |
| A3 | Numerical floors `pmax(..., 1e-10)` on all EBNM precision inputs |
| A4 | Dual convergence: both ΔL and Δβ must fall below tolerance |
| A5 | ELBO proxy (genomics log-likelihood) tracked per iteration |
| A6 | `refresh_taylor` flag for optional per-factor Taylor recomputation |

---

## Project Status

The model is fully implemented, tested, benchmarked against DeSurv, and documented in a 24-page report. See [`PROJECT_STATUS.qmd`](PROJECT_STATUS.qmd) for the complete session log.

**Completed:**
- Modular CAVI implementation (171/171 tests passing)
- DeSurv-aligned preprocessing pipeline (`code/preprocess_desurv.R`)
- Alpha CV selection via 1-SE rule (`code/select_alpha_cv.R`)
- SVD pseudoinverse prediction fix (`code/predict.R`)
- Synthetic validation: supervised C-index 0.79 > PCA 0.76 at n=300, p=1000
- PDAC cross-cohort benchmark: TCGA-only median external C-index 0.60 (5 cohorts),
  competitive with DeSurv's reported 0.60–0.65 range
- Prior sensitivity: `point_normal` vs `point_laplace` compared; `point_normal` recommended
- v2 preprocessing for merged cohort: intersect-first → log₂ → quantile normalization →
  top-2000 by merged variance → rank transform; fixes 838-gene selection bug
- Lambda sweep (λ ∈ {1, 5, 10, 20} × 3 priors): all β=0; λ≥5 collapses EL matrix — λ tuning ruled out
- EBMF diagnostic: 5/20 unsupervised factors are Cox-significant (C-index up to 0.63);
  confirms survival signal exists in merged data — failure is a **model problem, not data**
- EBMF warm-start experiments: β-only (EL fixed) → 6/20 factors active ✓; full CAVI → β
  collapses in 23 iters ✗. Root cause localised to `update_L_k()` A\_surv/A\_gen imbalance
- 24-page benchmark report (`results/benchmark_sim/ssbmf_summary_report.pdf`)

**Highest-priority next step:**
- Fix `update_L_k()` A\_surv/A\_gen scale imbalance so survival signal survives the joint CAVI.
  Full debugging plan in `docs/update_L_fix.md`. Four candidate fixes in priority order:
  (1) reorder β before L in inner loop, (2) β-only burn-in, (3) ridge Cox warm-start,
  (4) normalise A\_surv/A\_gen to comparable scales.

---

## Author

Andrew Walther — April 2026
