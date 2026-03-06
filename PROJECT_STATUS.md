# multiomicsGEP — Project Status

**Supervised Bayesian Matrix Factorization for Joint Genomics + Survival Modelling**

Author: Andrew Walther
Last updated: March 2026

---

## Project Overview

This project implements a **Supervised Bayesian Matrix Factorization** model
that jointly decomposes high-dimensional genomics data (gene expression, DNA
methylation, etc.) and models patient survival outcomes.  The core idea is to
find a shared latent space that simultaneously explains genomic variation *and*
predicts clinical outcomes.

**The model:**
- **Genomics:** Y = L F' + E, where L (n×K) are patient loadings, F (p×K) are
  biological factor weights, and E is Gaussian noise with feature-specific
  precision τ_j.
- **Survival:** Cox proportional hazards h(t_i) = h₀(t_i) exp(Σ_k l_{ik} β_k),
  linking the same latent loadings L to time-to-event data via survival
  coefficients β.

**Inference:** Coordinate Ascent Variational Inference (CAVI), where each
parameter update reduces to an Empirical Bayes Normal Means (EBNM) problem
with point-normal priors.

---

## Repository Map

```
multiomicsGEP/
├── PROJECT_STATUS.md                ← This file (project documentation)
├── README.md                        ← GitHub front page
├── code/
│   ├── SupervisedMF_Context.md      ← AI/collaborator quick-reference for V2 code
│   ├── Supervised_Bayesian_MF.R     ← V1: Original CAVI implementation (reference only)
│   ├── Supervised_Bayesian_MF_V2.R  ← V2: Corrected implementation (all fixes) ← USE THIS
│   └── multiomicsGEP_code.Rmd       ← Earlier exploratory R Markdown notebook
├── derivations/
│   ├── EBMF/
│   │   ├── EBMF_Derivations.pdf           ← Empirical Bayes MF theory notes
│   │   └── EBMF_Derivations_Latex.pdf     ← LaTeX-compiled version
│   ├── MF_UpdateDerivations/
│   │   ├── MF_Derivations_UpdateAlgo_1_15_26.pdf  ← Derivation v1 (Jan 15)
│   │   ├── MF_Derivations_UpdateAlgo_1_27_26.pdf  ← Derivation v2 (Jan 27)
│   │   ├── MF_Derivations_UpdateAlgo_1_29_26.pdf  ← Derivation v3 (Jan 29)
│   │   ├── MF_Derivations_UpdateAlgo_2_5_26.pdf   ← Derivation v4 (Feb 5)
│   │   ├── MF_Derivations_UpdateAlgo_2_12_26.pdf  ← Derivation v5 (Feb 12) — REVIEWED, HAS ERRORS
│   │   ├── MF_Derivations_UpdateAlgo_REVISED.tex/.pdf  ← Corrected derivations (R1-R8)
│   │   │                                                  + full step-by-step algebra (21 pages)
│   │   └── MF_V2_Companion.tex/.pdf        ← V2 companion: math ↔ code mapping (17 pages)
│   └── SurvivalMF/
│       ├── SupervisedMF_Likelihood_Estimation_Derivations.pdf
│       └── SurvivalMF_Derivations_Yusha.pdf
├── results/                         ← Simulation outputs (added Session 2)
│   ├── simulation_report.md/.pdf    ← 12-section report + rendered PDF (23 pages)
│   ├── run_simulation.R             ← Standalone script to reproduce all results
│   ├── simulation_console_output.txt← Verbatim run log
│   ├── figures/                     ← 8 figure pairs (PDF + PNG)
│   │   ├── fig1_rmse_trace.*        ← RMSE convergence over iterations
│   │   ├── fig2_elbo_proxy.*        ← ELBO proxy (genomics log-likelihood)
│   │   ├── fig3_beta_comparison.*   ← Estimated vs true survival coefficients
│   │   ├── fig4_gep_heatmap.*       ← Top-feature heatmap across GEPs
│   │   ├── fig5_kaplan_meier.*      ← KM curves by GEP loading tertile
│   │   ├── fig6_signal_recovery.*   ← Reconstructed Y vs true signal
│   │   ├── fig7_loading_correlations.* ← Permuted L̂ vs L_true correlation matrix
│   │   └── fig8_tau_distribution.*  ← Estimated noise precision per feature
│   └── *.csv                        ← 11 numeric summary tables (beta, C-index,
│                                       convergence history, factor summary, GEP top
│                                       features, loading correlations, PH test)
└── paper/
    └── multiomicsGEP_manuscript.qmd  ← Manuscript draft (Quarto)
```

---

## Mathematical Framework

### Model Components

| Symbol | Dimension | Meaning |
|--------|-----------|---------|
| Y      | n × p     | Observed genomics data matrix |
| L      | n × K     | Patient loading matrix (shared latent space) |
| F      | p × K     | Biological factor weight matrix |
| β      | K × 1     | Survival coefficient vector |
| τ_j    | p-vector  | Feature-specific noise precision |
| (t_i, δ_i) | n-vectors | Survival time and event indicator |

### Inference Strategy

1. **Mean-field factorisation:** q(L, F, β) = q(L) q(F) q(β)
2. **CAVI loop:** For each factor k = 1, ..., K:
   - Update q(l_k) via EBNM(x = B_L/A_L, s = 1/√A_L)
   - Update q(f_k) via EBNM(x = B_F/A_F, s = 1/√A_F)
   - Update q(β_k) via EBNM(x = B_β/A_β, s = 1/√A_β)
3. **Update τ_j** from expected squared residuals
4. **Cox Taylor expansion** linearises the non-conjugate survival likelihood
   into a weighted Gaussian form with working response z_i and weight W_{ii}

### Key Mathematical Concepts

- **Posterior means vs second moments:** EL = E_q[l], EL2 = E_q[l²] =
  Var_q(l) + mean². This distinction is critical for the variance inflation
  term in the τ update and the error-in-variables correction in the β update.
- **EBNM sub-problems:** Each coordinate update reduces to: given
  pseudo-observations x with noise s, estimate a prior g and posterior q.
- **Error-in-variables:** Using E[l²] (not l̄²) in the β precision prevents
  overfitting to uncertain loadings.

---

## Derivation Review Summary

A thorough review of the derivation document (MF_Derivations_UpdateAlgo_2_12_26.pdf)
identified **8 errors** that were corrected in `MF_Derivations_UpdateAlgo_REVISED.tex`:

| ID | Location | Error | Correction |
|----|----------|-------|------------|
| R1 | Sec. 1A | Wrong dimension subscripts on Y, L, F, E | Fixed: Y_{n×p} = L_{n×K} F'_{K×p} + E_{n×p} |
| R2 | Eq. 22 | z^{-k}_i defined with f_{jk'} instead of β_{k'} | Fixed: uses β_{k'} (survival coefficient) |
| R3 | Eqs. 32-35 | f²_{jk}, β²_k where E_q[·] second moments needed | Added explicit E_q[·] notation |
| R4 | Eqs. 47-48 | EBNM inputs subscripted i (sample) | Fixed to j (feature) for F update |
| R5 | Eq. 62 | Missing exponent: β_k should be β²_k | Added squared exponent |
| R6 | Sec. 6 | β update as generic WLS; no EBNM form | Derived full A_k, B_k for EBNM |
| R7 | Eq. 83 | Sign error: log τ + τR̄² | Fixed to log τ − τR̄² |
| R8 | Eq. 85 | l²_{jk} wrong index (j should be i) | Fixed to l²_{ik} |

### Document Relationships

- **2_12_26.pdf** → the reviewed document (contains errors R1-R8)
- **REVISED.tex** → corrected version with all 8 fixes + errata table
- **MF_V2_Companion.tex** → comprehensive companion mapping math to V2 code

---

## Code Versions

### V1: `Supervised_Bayesian_MF.R`

The original implementation. Functionally correct but with several
known issues:

- Orthogonalisation always ON; resets EL2/EF2 to squared means (drops
  posterior variance every 10 iterations)
- z_no_k computed from start-of-iteration η only (not true Gauss-Seidel)
- No numerical floors on EBNM precision inputs
- Convergence checks only ΔL (not Δβ)
- Only RMSE tracked (no ELBO monitoring)

### V2: `Supervised_Bayesian_MF_V2.R`

The corrected implementation with 6 algorithmic improvements:

| ID | Change | Description |
|----|--------|-------------|
| A1 | Gauss-Seidel CAVI | z_no_k recomputed from current EL, EBeta inside k-loop |
| A2 | Orthogonalisation flag | Behind `orthogonalize=FALSE` (default OFF) |
| A3 | Precision floors | `pmax(..., 1e-10)` on A_L, A_F, A_Beta |
| A4 | Dual convergence | Both mean\|ΔL\| and mean\|Δβ\| must be < tol |
| A5 | ELBO proxy | Genomics log-likelihood tracked per iteration |
| A6 | Taylor refresh flag | `refresh_taylor=FALSE` for optional per-k recomputation |

**New function signature:**
```r
fit_supervised_mf(Y, time, status, K=5, max_iter=100, tol=1e-5,
                  orthogonalize=FALSE, refresh_taylor=FALSE, verbose=TRUE)
# Returns: list(L, F, Beta, Beta2, Tau, history)
```

---

## How to Run

### Prerequisites

R packages required:
```r
install.packages(c("survival", "ebnm"))
```

### Running the V2 Simulation

```r
source("code/Supervised_Bayesian_MF_V2.R")
```

This will:
1. Generate synthetic data (n=250 patients, p=1000 features, K=5 factors)
2. Run the CAVI algorithm with default settings
3. Print: factor summary table, PH test, C-index comparison, estimated vs true β
4. Render diagnostic plots: RMSE trace, ELBO proxy, GEP heatmap, KM curves,
   signal recovery scatter

### Expected Output

- **RMSE** converges near 1.0 (the true noise SD)
- **ELBO proxy** is non-decreasing
- **Estimated β** recovers sign pattern (+, −, +, −, 0)
- **Supervised C-index** exceeds Top-5 PCA C-index
- **Censoring rate** ~30-40%

### Compiling LaTeX Documents

```bash
cd derivations/MF_UpdateDerivations/
pdflatex MF_Derivations_UpdateAlgo_REVISED.tex
pdflatex MF_V2_Companion.tex
# Run twice for cross-references:
pdflatex MF_V2_Companion.tex
```

---

## Current Status & Next Steps

### Completed

- [x] Derivation review: 8 errors (R1-R8) identified and corrected
- [x] REVISED.tex/PDF: Corrected derivation document with errata table (21 pages)
- [x] REVISED.tex expanded: Full step-by-step algebra added throughout (Taylor,
  L/F/β/τ updates) with new `derivbox` gray-box environment for intermediate steps
- [x] V2.R: Corrected R implementation with all 6 algorithmic fixes (A1-A6)
- [x] V2_Companion.tex/PDF: Comprehensive math-to-code companion document (17 pages)
- [x] PROJECT_STATUS.md + SupervisedMF_Context.md: Project documentation
- [x] Simulation validation: V2 runs end-to-end on synthetic data (n=250, p=1000, K=5)
- [x] results/: Full simulation outputs — 8 figures, 11 CSV tables, console log
- [x] results/simulation_report.md/.pdf: 12-section simulation report (23 pages)

### Potential Next Steps

- [ ] **Real data application:** Apply V2 to actual genomics + survival datasets
  (TCGA, GEO, etc.) using the `DATA_MODE = "real"` pathway in V2.R
- [ ] **Select K:** Implement cross-validated C-index or held-out log-likelihood
  to choose the number of factors K objectively (currently set to 5)
- [ ] **Additional prior families:** Test `prior_family = "point_laplace"` or
  `"normal_scale_mixture"` in EBNM calls for different sparsity structures
- [ ] **Full ELBO:** Track the complete ELBO (survival Taylor term + all KL
  divergences) instead of the genomics-only proxy currently used
- [ ] **Scalability:** Profile and optimise for large p (e.g., p > 10,000)
  — focus on the n×p `Var_Term` and `R_k` matrix operations
- [ ] **Manuscript:** Complete `paper/multiomicsGEP_manuscript.qmd`; reference
  V2 results and corrected derivations
- [ ] **Bibliography:** Add `\bibliography{refs}` to REVISED.tex to resolve
  the `\citep{wang2022}` undefined-reference warning in the compiled PDF

---

## Session Log

### Session 1 (March 4–5, 2026)

**Derivation review and code correction session.**

1. Reviewed `MF_Derivations_UpdateAlgo_2_12_26.pdf` against the code
2. Identified 8 derivation errors (R1–R8)
3. Created `MF_Derivations_UpdateAlgo_REVISED.tex` with all corrections + errata table
4. Analysed V1 code against corrected derivations; identified 6 improvements
5. Created `Supervised_Bayesian_MF_V2.R` with all algorithmic fixes (A1–A6)
6. Created `MF_V2_Companion.tex` (math ↔ code companion, 17 pages)
7. Created `PROJECT_STATUS.md` and `SupervisedMF_Context.md`

**Key decisions:**
- Orthogonalisation defaults to OFF in V2 (Point-Normal EBNM promotes
  sparsity/distinctness without rotation)
- Taylor refresh defaults to OFF (standard IRLS-within-VI; more stable)
- V1 preserved as-is for reference; V2 is the recommended implementation

---

### Session 2 (March 5, 2026)

**Simulation validation, derivation expansion, and documentation session.**

1. **Expanded `REVISED.tex`** from ~650 lines to 1,199 lines with full
   intermediate algebra throughout every update section:
   - Added new `derivbox` tcolorbox (gray) for algebra steps, distinct from
     `correctionbox` (red) and `ebnmbox` (blue)
   - **Taylor (Sec 3):** Full 15-step completing-the-square derivation; derives
     working response $z_i = \hat{\eta}_i + u_i/W_{ii}$ from first principles
   - **L update (Sec 4):** Explicit residual substitution, squared-term expansion,
     $\mathbb{E}_q[\cdot]$ application, and A/B identification for genomics and
     survival terms separately then combined
   - **F update (Sec 5):** Explicit justification why F absent from survival term;
     full quadratic grouping to $A_{jk}$, $B_{jk}$ with R4 subscript fix shown
   - **β update (Sec 6):** $g(\beta)$ log-expectation derivation; variance-mean
     decomposition $\mathbb{E}[l^2] = \mathrm{Var}(l) + \bar{l}^2$ shown explicitly;
     cross-term expansion; mean-field independence argument; error-in-variables
     correction $A_k = \sum_i W_{ii}\overline{l^2_{ik}}$ (R5, R6)
   - **τ update (Sec 7):** Full $\bar{R}^2_{ij}$ expansion under mean-field;
     explicit $\partial\mathcal{F}/\partial\tau = 0$ optimisation step (R7, R8)
   - Fixed R2 correction box layout (changed `\paragraph` to `\subsubsection*`
     to prevent run-in heading from clipping the tcolorbox off the page edge)
   - Recompiled REVISED.pdf: **21 pages** (was 12)

2. **Ran V2.R simulation** (n=250, p=1000, K=5); saved all outputs to `results/`:
   - 8 figure pairs (PDF + PNG), 11 CSV tables, console log, run script
   - Key results: RMSE = 0.9978, C-index PCA = 0.827 vs Supervised = 0.828,
     PH test p = 0.258 (no violation), convergence in 45 iterations

3. **Wrote `results/simulation_report.md`** — 12-section comprehensive report
   with all figures embedded, factor summary tables, GEP top-feature tables,
   discussion, and application roadmap for real multi-omics data

4. **Rendered `results/simulation_report.pdf`** — 23 pages via pandoc + xelatex

**Deliverables committed:** `5da14f3`, `e424d14` (pushed to `origin/main`)
