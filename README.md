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
│   ├── benchmark_sim/                 ← ✅ Formal benchmark pipeline (canonical)
│   │   ├── run_LB_benchmark.R         ← Cluster A runner (Lβ linear predictor)
│   │   ├── run_YFB_benchmark.R        ← Cluster B runner (YFβ / Cox-on-YF)
│   │   ├── run_phase1_diagnostics.R   ← Loading heatmaps
│   │   ├── archive/                   ← 9 retired scripts (see archive/README.md)
│   │   └── outputs/
│   │       ├── LB_benchmark/          ← Cluster A outputs
│   │       └── YFB_benchmark/         ← Cluster B outputs
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
# Cluster A benchmark (Lβ linear predictor, alpha CV, external validation)
Rscript results/benchmark_sim/run_LB_benchmark.R

# Cluster B benchmark (Cox-on-YF / YFβ reformulation)
Rscript results/benchmark_sim/run_YFB_benchmark.R
```

Outputs go to `results/benchmark_sim/outputs/LB_benchmark/` and `outputs/YFB_benchmark/` respectively. Benchmark reports are versioned by date and live in `docs/reports/` (e.g., `docs/reports/ssbmf_summary_report_04_29_26.pdf` for the archived DeSurv benchmark).

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

**Completed:**
- Modular CAVI implementation (171/171 tests passing)
- Cox-on-YF reformulation (Cluster B): `code/fit_cox_on_yf.R`, `code/predict_cox_on_yf.R`
- DeSurv-aligned preprocessing pipeline (`code/preprocess_desurv.R`)
- Alpha CV selection via 1-SE rule (`code/select_alpha_cv.R`)
- Phase C sign correction for both LB (`code/fit_modular.R`) and YFB (`code/fit_cox_on_yf.R`):
  post-convergence training concordance check; flips EBeta if C_train < 0.5
- SVD pseudoinverse prediction fix (`code/predict.R`)
- Synthetic validation: LB supervised C-index 0.79 > PCA 0.76; YFB C=0.906 (K_eff=4)
- PDAC cross-cohort benchmark (Phase C fixed): LB external C=0.55–0.67 across all training
  modes (tcga_only, cptac_only, merged); YFB external C=0.55–0.63 (single-cohort, normal prior)
- Prior sensitivity: `point_normal` vs `point_laplace` vs `normal` compared across both models;
  `normal` prior needed for YFB on real PDAC (spike-and-slab collapses beta when signal is weak)
- v2 preprocessing for merged cohort: intersect-first → log₂ → quantile normalization →
  top-2000 by merged variance → rank transform
- DeSurv benchmark report (`docs/reports/ssbmf_summary_report_04_29_26.pdf`) and full
  Phase A–C re-benchmark report (`docs/reports/ssbmf_summary_report_05_05_26.pdf`)

**Current priorities:** K selection via CV (`code/select_K.R` stub); YFB merged beta collapse
(K_eff=0 for YFB on mixed RNA-seq + proteomics); prior comparison follow-up. See `ROADMAP.md`.

---

## Author

Andrew Walther — May 2026
