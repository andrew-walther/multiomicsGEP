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
├── README.md                         ← You are here
├── PROJECT_STATUS.md                 ← Full project documentation & session log
│
├── code/
│   ├── fit_modular.R                 ← ✅ Canonical CAVI loop (factor-wise, calls update_*.R)
│   ├── Supervised_Bayesian_MF_V2.R   ← ✅ Monolithic reference implementation (V2)
│   ├── update_beta.R                 ← Modular β update (scalar EBNM, Cox survival)
│   ├── update_L.R                    ← Modular L update (vector EBNM, dual-source)
│   ├── update_F.R                    ← Modular F update (vector EBNM, pure genomics)
│   ├── update_tau.R                  ← Modular τ update (closed-form MLE)
│   ├── SupervisedMF_Context.md       ← AI/developer quick-reference for the code
│   └── legacy/                       ← Archived files (V1, early scripts)
│       ├── Supervised_Bayesian_MF.R  ← V1 original (reference only, known issues)
│       ├── execute_update_beta.R     ← Early demo (superseded by demos/)
│       └── multiomicsGEP_code.Rmd    ← Early exploratory notebook
│
├── docs/                              ← Companion documentation (PDF + HTML)
│   ├── Makefile                       ← `make all` renders .qmd → .pdf + .html via Quarto
│   ├── fit_modular.qmd/.pdf/.html     ← ✅ fit_modular.R walkthrough (full CAVI loop)
│   ├── update_beta.qmd/.pdf/.html     ← β update: code walkthrough, tests, demos
│   ├── update_L.qmd/.pdf/.html        ← L update: code walkthrough, tests, demos
│   ├── update_F.qmd/.pdf/.html        ← F update: code walkthrough, tests, demos
│   └── update_tau.qmd/.pdf/.html      ← τ update: code walkthrough, tests, demos
│
├── tests/                             ← 105 tests (run: Rscript tests/run_tests.R)
│   ├── run_tests.R                    ← Master test runner
│   └── test_update_*.R                ← Per-module test suites
│
├── demos/                             ← Interactive demonstrations (5 per module)
│   └── demo_update_*.R                ← Run: Rscript demos/demo_update_*.R
│
├── results/                           ← Simulation outputs (grouped by implementation)
│   ├── modular_sim_factor/            ← ✅ Factor-wise CAVI (canonical)
│   │   ├── run_factor_modular_simulation.R
│   │   └── factor_modular_sim_report.qmd/.pdf/.html
│   ├── full_sim/                      ← V2 monolithic simulation
│   │   ├── run_simulation.R
│   │   └── simulation_report.qmd/.pdf
│   ├── modular_sim_block/             ← Block-wise modular (deprecated)
│   │   ├── run_modular_simulation.R
│   │   └── modular_sim_report.qmd/.pdf/.html
│   ├── figures/
│   │   ├── modular_sim_factor/        ← 8 figure pairs (factor-wise)
│   │   ├── full_sim/                  ← 8 figure pairs (V2 monolithic)
│   │   └── modular_sim/               ← 8 figure pairs (block-wise, deprecated)
│   └── tables/
│       ├── modular_sim_factor/        ← 11 CSV tables (factor-wise)
│       ├── full_sim/                  ← 11 CSV tables (V2 monolithic)
│       └── modular_sim/               ← 11 CSV tables (block-wise, deprecated)
│
├── derivations/
│   ├── MF_UpdateDerivations/
│   │   ├── MF_Derivations_UpdateAlgo_REVISED.pdf  ← ✅ Corrected derivations
│   │   ├── MF_Derivations_UpdateAlgo_REVISED.tex  ←    LaTeX source
│   │   ├── MF_V2_Companion.pdf                    ← ✅ Math ↔ code companion doc
│   │   ├── MF_V2_Companion.tex                    ←    LaTeX source
│   │   └── MF_Derivations_UpdateAlgo_*.pdf        ←    Historical derivation drafts
│   ├── EBMF/
│   │   └── EBMF_Derivations*.pdf     ← Empirical Bayes Matrix Factorization theory
│   └── SurvivalMF/
│       └── *.pdf                     ← Survival + MF background notes
│
└── paper/
    └── multiomicsGEP_manuscript.qmd  ← Manuscript draft (in progress)
```

---

## Quickstart

### Prerequisites

```r
install.packages(c("survival", "ebnm"))
```

### Run the Simulation Benchmark

**V2 implementation (monolithic):**
```r
source("code/Supervised_Bayesian_MF_V2.R")
```

**Modular implementation (canonical, recommended):**
```r
Rscript results/modular_sim_factor/run_factor_modular_simulation.R
```

Both scripts use the same parameters (n=250, p=1000, K=5, seed=42) and produce equivalent results. The modular script uses `fit_supervised_mf_modular()` from `code/fit_modular.R`, which implements factor-wise Gauss-Seidel CAVI via the four standalone update modules. Results are written to `results/tables/modular_sim_factor/` and `results/figures/modular_sim_factor/`; the rendered report is at `results/modular_sim_factor/factor_modular_sim_report.qmd`.

The V2 run generates a synthetic dataset and prints:

- Factor summary table (β estimates, log-rank p-values, sparsity, PVE)
- Proportional hazards test (cox.zph)
- C-index: supervised latent space vs. top-5 PCA
- Estimated vs. true β coefficients
- Top 5 features per GEP

It also renders a diagnostic dashboard: RMSE trace, ELBO proxy, GEP heatmap, Kaplan-Meier curves, and signal recovery plot.

### Apply to Real Data

The recommended entry point for real data is `fit_supervised_mf_modular()` in `code/fit_modular.R`:

```r
source("code/fit_modular.R")   # also sources update_L/F/beta/tau.R automatically

res <- fit_supervised_mf_modular(
  Y      = your_matrix,    # numeric matrix: n patients × p genes (pre-normalised)
  time   = your_time,      # numeric vector: survival/censoring time
  status = your_status,    # integer vector: 1 = event, 0 = censored
  K      = 5,              # number of latent factors (select via cross-validation)
  max_iter = 100,
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
| Modular | `code/fit_modular.R` + `update_*.R` | ✅ **Current** | Factor-wise Gauss-Seidel CAVI; tested (105/105); recommended for all new work |

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

The model derivation and V2 implementation are complete and validated on simulated data. See [`PROJECT_STATUS.md`](PROJECT_STATUS.md) for the full status including derivation review notes, session log, and prioritised next steps.

**Immediate next steps:**
- Real data application (TCGA or similar)
- Cross-validated selection of K (number of factors)
- Full ELBO tracking (currently genomics-only proxy)

---

## Author

Andrew Walther — March 2026
