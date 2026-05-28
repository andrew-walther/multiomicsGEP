# multiomicsGEP

**Supervised Bayesian Matrix Factorization for Joint Genomics and Survival Modelling**

---

## Overview

This project implements a **Supervised Bayesian Matrix Factorization** model that jointly decomposes high-dimensional genomics data — such as gene expression or DNA methylation matrices — while simultaneously modelling patient survival outcomes via a Cox proportional hazards model.

The core insight is that standard unsupervised factorization (e.g. PCA, NMF) finds latent structure that explains genomic variance, but has no reason to recover factors that are *clinically meaningful*. By supervising the factorization with survival data, this model discovers **Gene Expression Programs (GEPs)** that are both genomically coherent and prognostically relevant.

**The model:**

```
Y (n×p) = L (n×K) × F' (K×p) + E        ← matrix factorization (genomics)
h(t_i)  = h₀(t_i) exp(ηᵢ)               ← Cox proportional hazards (survival)
```

- **L** — patient loading matrix: coordinates of each patient in the latent factor space
- **F** — factor weight matrix: gene loadings that define each program
- **β** — survival coefficients: prognostic weight of each factor
- **τ** — feature-specific noise precision

Two parameterizations of the linear predictor η are implemented and benchmarked:

| Model | Linear predictor | Key property |
|-------|-----------------|--------------|
| LB    | η = Lβ          | Factor scores learned jointly with survival; `code/fit_modular.R` |
| YFB   | η = (YF)β       | Predictor computed directly from observed expression; eliminates train/test mismatch in projection; `code/fit_cox_on_yf.R` |

Inference is performed via **Coordinate Ascent Variational Inference (CAVI)**, where each variational update is solved as an **Empirical Bayes Normal Means (EBNM)** problem — promoting sparsity in both F and β through point-normal or normal priors.

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
│   ├── benchmark_sim/                 ← ✅ Formal benchmark pipeline (canonical)
│   │   ├── run_LB_benchmark.R         ← LB model runner (η = Lβ, alpha CV, external validation)
│   │   ├── run_YFB_benchmark.R        ← YFB model runner (η = (YF)β, Cox-on-YF)
│   │   ├── run_phase1_diagnostics.R   ← Loading heatmaps
│   │   ├── archive/                   ← 9 retired scripts (see archive/README.md)
│   │   └── outputs/
│   │       ├── LB_benchmark/          ← LB model outputs
│   │       └── YFB_benchmark/         ← YFB model outputs
│   ├── figures/                       ← Active per-cohort figure outputs
│   ├── tables/                        ← Active per-cohort table outputs
│   └── legacy/                        ← Retired simulation generations
│       ├── full_sim/                  ← V2 monolithic simulation
│       ├── modular_sim_block/         ← Block-wise modular (deprecated)
│       ├── modular_sim_factor/        ← Factor-wise exploratory runner + PDAC/synthetic reports
│       ├── figures/                   ← Legacy figure outputs (full_sim, modular_sim, synthetic)
│       └── tables/                    ← Legacy table outputs
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

### Run the Formal Benchmark (current entry point)

```bash
# LB model benchmark (η = Lβ; alpha mixing CV-selected per training set)
Rscript results/benchmark_sim/run_LB_benchmark.R

# YFB model benchmark (η = (YF)β; Cox-on-YF reformulation)
Rscript results/benchmark_sim/run_YFB_benchmark.R
```

Outputs go to `results/benchmark_sim/outputs/LB_benchmark/` and `outputs/YFB_benchmark/` respectively. Benchmark reports are versioned by date in `docs/reports/` — the most recent is `ssbmf_summary_report_05_05_26.pdf`.

### Run the Real PDAC Analysis

For real PDAC data, use `fit_supervised_mf_modular()` in `code/fit_modular.R` directly, or use the archived exploratory runner `results/legacy/modular_sim_factor/run_factor_modular_simulation.R`. The active benchmark runners (`run_LB_benchmark.R`, `run_YFB_benchmark.R`) are the canonical entry points for formal evaluation.

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

The model is fully implemented, tested, and benchmarked. Two model variants are evaluated:
**Cluster A (LB)** uses η = Lβ; **Cluster B (YFB)** uses η = (YF)β (Cox-on-YF reformulation).
See [`PROJECT_STATUS.qmd`](PROJECT_STATUS.qmd) for the complete session log.

**Current state (2026-05-28):** 246/246 tests passing. Recommended configuration: YFB with
DeSurv-aligned gene selection (combined mean+variance rank, top-3000 per cohort, 2064 genes),
per-platform z-standardization, K=7 — mean external C-index 0.636 across 5 held-out PDAC
cohorts. Two active prognostic programs identified (one adverse, one protective); signal
replicates across RNA-seq, microarray, and proteomics platforms.

**Next:** pathway enrichment on active factor gene weights; merge to main. See `ROADMAP.md`.

---

## Author

Andrew Walther — May 2026
