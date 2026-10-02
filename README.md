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

## Multimodal YFB (in progress)

An extension of YFB to matched gene expression and DNA methylation (TCGA/ICGC),
developed on the `codex/multimodal-yfb` branch.
[PDF](derivations/multimodal_YFB/multimodal_YFB_derivation.pdf) and
[editable Quarto source](derivations/multimodal_YFB/multimodal_YFB_derivation.qmd).
The nine-page document derives linear joint YFB with separate loading priors,
the diagonal Cox working approximation, and L/F/β/Tau updates. The linked
[R implementation specification](derivations/multimodal_YFB/multimodal_YFB_implementation.qmd)
maps the derivation to the isolated implementation. The runnable matched-block
fitter, frozen predictor, simulation generator, and K-selection runner are in
`code/fit_multimodal_yfb.R`, `code/predict_multimodal_yfb.R`,
`code/simulate_multimodal_yfb.R`, and `code/run_multimodal_yfb_simulation.R`.
The multimodal module is not the production single-modality PDAC model; its
replicated K-selection study is still in progress.

Options added on 2026-10-01 (defaults keep the original behaviour):

- **Loading priors per modality.** `control$prior_F` sets point-exponential,
  point-Laplace or Normal, with the signed priors fit by `ebnm`.
- **Factor pruning, which together selects K.**
  - a per-feature intercept (`intercept`);
  - marginal-likelihood prior fits (`prior_update = "ebnm"`);
  - an ELBO nullcheck that removes factors (`prune`), with the Cox partial
    log-likelihood as the survival term.
- **Faster sweep.** The per-feature loading sweep is compiled C++
  (`code/multimodal_yfb_sweep.cpp`), with an automatic fallback to R.

Matched TCGA (training) and ICGC (validation) expression + methylation cohorts are
built by `code/load_multiomics_data.R`, which reads a local, git-ignored `data/`
folder. The real-data runner is `results/multimodal_real/run_multimodal_real_fit.R`;
the simulation comparison of fitting variants is
`results/benchmark_sim/run_multimodal_pruning_comparison.R`.
The original extensive derivation is preserved as technical notes in the same folder.
Run the base-R mathematical checks with
`Rscript derivations/multimodal_YFB/verify_derivation.R`.
Rebuild the PDF with
`quarto render derivations/multimodal_YFB/multimodal_YFB_derivation.qmd --to pdf`;
rendering also executes the checks and stops on assertion failure.

---

## Repository Structure

```
multiomicsGEP/
├── README.md  LICENSE (MIT)
├── CLAUDE.md                ← agent instructions (AGENTS.md is a symlink to it)
├── PROJECT_STATUS.qmd/.pdf  ← full project documentation and development log
├── DECISIONS.md             ← dated architectural and analytical decisions
├── ROADMAP.md               ← prioritized next steps and completed items
├── config/globals.yml       ← single source of truth for hyperparameters
├── code/
│   ├── fit_modular.R        ← canonical CAVI loop for LB (η = Lβ)
│   ├── fit_cox_on_yf.R      ← YFB (η = (YF)β); predict_cox_on_yf.R for hold-out scoring
│   ├── update_*.R           ← modular CAVI updates (β, L, F, τ, and cohort/YFB variants)
│   ├── compute_elbo.R  compute_bic.R  compute_cv_loglik.R   ← ELBO, BIC, held-out log-likelihood
│   ├── select_K.R  select_alpha_cv.R  select_k_alpha_bo.R   ← K and alpha selection
│   ├── preprocess_desurv.R  feature_selection.R  train_test_split.R  predict.R
│   ├── pathway_enrichment.R  concordance_ci.R
│   ├── *multimodal_yfb*.R   ← multimodal YFB fitter, updates, prediction, simulation (in progress)
│   ├── Supervised_Bayesian_MF_V2.R   ← monolithic V2 reference (do not extend)
│   ├── SupervisedMF_Context.md       ← math ↔ code quick reference
│   └── legacy/              ← V1 and early scripts (archived)
├── tests/                   ← run_tests.R plus one test file per module
├── results/
│   ├── benchmark_sim/       ← formal benchmark pipeline: LB/YFB runners, K and alpha CV,
│   │                          external validation, DeSurv and EBMF comparisons
│   ├── figures/  tables/    ← per-cohort outputs
│   └── legacy/              ← retired simulation generations
├── derivations/             ← corrected CAVI derivations (MF_UpdateDerivations/), per-update
│                              derivations (qB, qL, qF, qTau), multimodal_YFB/
├── docs/
│   ├── progress_book/       ← Quarto book, one chapter per advisor meeting
│   │                          (https://andrew-walther.github.io/multiomicsGEP/)
│   ├── plans/               ← working plans (current: Working_Plan_10_1_26.md)
│   ├── reports/             ← dated analysis reports
│   └── *.qmd/.pdf/.html     ← per-update walkthroughs, PDAC data audit
├── presentation/            ← lab-meeting decks (latest: walther_lab_meeting_08_27_2026/)
├── longleaf_setup/          ← UNC Longleaf SLURM scripts
├── paper/                   ← manuscript draft (multiomicsGEP_manuscript.qmd)
└── data/                    ← local TCGA/ICGC copies; git-ignored, never committed
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

Outputs go to `results/benchmark_sim/outputs/`. Reports are versioned by date in `docs/reports/` — see `desurv_alignment_report_05_27_26.pdf` for the most recent external validation results.

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
| Modular | `code/fit_modular.R` + `update_*.R` | ✅ **Current** | Factor-wise Gauss-Seidel CAVI; tested (`Rscript tests/run_tests.R`); recommended for all new work |

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

The single-modality model is implemented, tested, and benchmarked. The recommended
configuration uses the YFB linear predictor (η = (YF)β) with DeSurv-aligned gene
selection (genes ranked jointly by mean expression and variance within each platform,
top 3,000 per cohort, about 2,064 after intersection). Trained on 273 PDAC patients
(TCGA + CPTAC), it identifies two active prognostic programs, one associated with worse
survival and one with better survival, with mean external concordance 0.627 across five
independent PDAC cohorts (RNA-seq, microarray, proteomics).

Test suite (2026-10-01, `codex/multimodal-yfb`): 509 passing.

**Current work** (prioritized in `docs/plans/Working_Plan_10_1_26.md`; see `ROADMAP.md`):
the multimodal YFB extension, matched TCGA/ICGC preprocessing, and whether survival
supervision can be made to inform the factor matrix F (the `alpha_F` question).
Full history: [`PROJECT_STATUS.qmd`](PROJECT_STATUS.qmd); decision rationale: `DECISIONS.md`.

---

## License

MIT (see `LICENSE`). PDAC and TCGA/ICGC data are not included and keep their own terms.

---

## Author

Andrew Walther — May 2026

## Chapter 4 proposal

The canonical chapter is `paper/prelim/project3-ssbmf.qmd`. Its installed
post-commit hook renders and commits the generated Chapter 4 copy in
bios-dissertation; rerun that repository's `prelim/tools/build_prelim.py`
to update the full prelim. Continuing work is described as methodological
aims and expected results without a fixed-date timeline. A substantive
revision will follow the advisor-meeting plan and verified additional results.
