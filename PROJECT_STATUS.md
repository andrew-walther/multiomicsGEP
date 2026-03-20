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
│   ├── update_beta.R                ← Modular β update functions
│   ├── update_L.R                   ← Modular q(L) update functions (vector EBNM, dual-source)
│   ├── update_F.R                   ← Modular q(F) update functions (pure genomics, τ cancellation)
│   ├── update_tau.R                 ← Modular q(τ) update functions (closed-form MLE)
│   └── multiomicsGEP_code.Rmd       ← Earlier exploratory R Markdown notebook
├── derivations/
│   ├── qB/                              ← Self-contained β derivation
│   │   ├── qBeta_update_derivation.tex  ← Full step-by-step q(β) derivation (11 pages)
│   │   └── qBeta_update_derivation.pdf  ← Compiled PDF
│   ├── qL/                              ← Self-contained q(L) derivation
│   │   ├── qL_update_derivation.tex     ← Vector EBNM, dual-source (genomics+survival)
│   │   └── qL_update_derivation.pdf     ← Compiled PDF
│   ├── qF/                              ← Self-contained q(F) derivation
│   │   ├── qF_update_derivation.tex     ← τ cancellation property, pure genomics
│   │   └── qF_update_derivation.pdf     ← Compiled PDF
│   ├── qTau/                            ← Self-contained q(τ) derivation
│   │   ├── qTau_update_derivation.tex   ← Closed-form MLE, variance correction
│   │   └── qTau_update_derivation.pdf   ← Compiled PDF
│   ├── EBMF/
│   │   ├── EBMF_Derivations.pdf           ← Empirical Bayes MF theory notes
│   │   └── EBMF_Derivations_Latex.pdf     ← LaTeX-compiled version
│   ├── MF_UpdateDerivations/
│   │   ├── MF_Derivations_UpdateAlgo_1_15_26.pdf  ← Derivation v1 (Jan 15)
│   │   ├── MF_Derivations_UpdateAlgo_1_27_26.pdf  ← Derivation v2 (Jan 27)
│   │   ├── MF_Derivations_UpdateAlgo_1_29_26.pdf  ← Derivation v3 (Jan 29)
│   │   ├── MF_Derivations_UpdateAlgo_2_5_26.pdf   ← Derivation v4 (Feb 5)
│   │   ├── MF_Derivations_UpdateAlgo_2_12_26.pdf  ← Derivation v5 (Feb 12) — REVIEWED, HAS ERRORS
│   │   ├── MF_Derivations_UpdateAlgo_3_6_26.pdf   ← March 6 draft (Section 6A: EBMF focus)
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
├── tests/                           ← Test infrastructure (105/105 passing)
│   ├── test_helpers.R               ← Lightweight assertion framework (assert_near, run_test, etc.)
│   ├── test_update_beta.R           ← 24 tests for update_beta.R (9 groups, TDD)
│   ├── test_update_L.R              ← 28 tests for update_L.R (9 groups)
│   ├── test_update_F.R              ← 26 tests for update_F.R (9 groups)
│   ├── test_update_tau.R            ← 27 tests for update_tau.R (9 groups)
│   └── run_tests.R                  ← Master test runner — 105/105 tests passing
├── demos/                           ← Interactive demonstrations (5 per module)
│   ├── demo_update_beta.R           ← 5 demos for update_beta.R
│   ├── demo_update_L.R              ← 5 demos for update_L.R
│   ├── demo_update_F.R              ← 5 demos for update_F.R
│   └── demo_update_tau.R            ← 5 demos for update_tau.R
├── docs/                            ← Companion documentation (MD + PDF + HTML)
│   ├── Makefile                     ← `make all` renders .md → .pdf + .html
│   ├── update_beta.md/.pdf/.html    ← Companion doc for update_beta.R
│   ├── update_L.md/.pdf/.html       ← Companion doc for update_L.R
│   ├── update_F.md/.pdf/.html       ← Companion doc for update_F.R
│   └── update_tau.md/.pdf/.html     ← Companion doc for update_tau.R
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

# q(β) derivation (Session 3):
cd derivations/qB/
pdflatex qBeta_update_derivation.tex
pdflatex qBeta_update_derivation.tex
pdflatex qBeta_update_derivation.tex
```

### Running the Test Suite

```r
Rscript tests/run_tests.R
# Expected: 105/105 tests passed (24 β + 28 L + 26 F + 27 τ)
```

### Running the Modular Update Demos

```r
Rscript demos/demo_update_beta.R   # β update: 5 demos
Rscript demos/demo_update_L.R      # q(L) update: 5 demos
Rscript demos/demo_update_F.R      # q(F) update: 5 demos — see τ cancellation
Rscript demos/demo_update_tau.R    # q(τ) update: 5 demos — see variance correction
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
- [x] derivations/qB/qBeta_update_derivation.tex/.pdf: Self-contained q(β) derivation (11 pages)
- [x] code/update_beta.R: Modular β update (`compute_z_no_k`, `update_beta_k`, `update_beta_all`)
- [x] tests/: Test infrastructure + 24 tests for β update (9 groups, all passing)
- [x] demos/demo_update_beta.R: 5 interactive demos for β update (anatomy, recovery, shrinkage, error-in-variables, K=5)
- [x] code/update_L.R: Modular q(L) update (`compute_R_k`, `update_L_k`, `update_L_all`) — vector EBNM, dual-source
- [x] code/update_F.R: Modular q(F) update (`update_F_k`, `update_F_all`) — pure genomics, τ cancellation
- [x] code/update_tau.R: Modular q(τ) update (`compute_var_term`, `compute_expected_residual_sq`, `update_tau`) — closed-form MLE
- [x] tests/test_update_L.R: 28 tests for update_L.R (9 groups, all passing)
- [x] tests/test_update_F.R: 26 tests for update_F.R (9 groups, τ cancellation verified)
- [x] tests/test_update_tau.R: 27 tests for update_tau.R (9 groups, variance correction verified)
- [x] tests/run_tests.R: Updated master runner — 105/105 tests passing
- [x] demos/demo_update_L.R: 5 demos (anatomy, signal recovery, genomics vs combined, error-in-variables, K=5)
- [x] demos/demo_update_F.R: 5 demos (τ cancellation, sparse recovery, differential shrinkage, sparsity, K=5)
- [x] demos/demo_update_tau.R: 5 demos (anatomy, known noise, variance correction, heteroscedastic, ELBO proxy)
- [x] derivations/qL/qL_update_derivation.tex/.pdf: Self-contained q(L) derivation
- [x] derivations/qF/qF_update_derivation.tex/.pdf: Self-contained q(F) derivation
- [x] derivations/qTau/qTau_update_derivation.tex/.pdf: Self-contained q(τ) derivation
- [x] docs/: Companion documentation for all 4 modular update scripts (.md + .pdf + .html)
- [x] roxygen2 @export/@family/@seealso tags added to all 11 functions
- [x] Targeted inline comments added to all 4 R scripts

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

---

### Session 3 (March 10, 2026)

**q(β) derivation, modular implementation, and TDD test suite.**

1. **Verified q(β) derivation** across REVISED.tex (Sec. 6) and March 6 PDF (Sec. 6A):
   - Both agree on final result: A_k = Σᵢ Wᵢᵢ·E_q[l²ᵢₖ], B_k = Σᵢ Wᵢᵢ·z^{-k}ᵢ·l̄ᵢₖ
   - Confirmed March 6 PDF Eq. 68 typo (missing β² exponent) is intermediate-step
     only; already flagged as R5 in REVISED.tex; final Eqs. 92–93 are correct
   - Confirmed V2.R omits 1/σ² prior term intentionally (EBNM handles it empirically)
   - Proved z_no_k reuse is mathematically valid (does not depend on l_{ik} or β_k)
   - Confirmed error-in-variables direction: larger EL2 → larger A_k → smaller x_k
     (more shrinkage), and s_k = 1/√A_k is SMALLER (not larger) when A_k increases

2. **Installed TDD and GSD skills** (Matt Pocock TDD skill, GSD v1.22.4):
   - TDD enforces red-green-refactor discipline
   - GSD provides discuss→plan→execute→verify workflow structure

3. **Created test infrastructure** (`tests/` directory):
   - `tests/test_helpers.R`: Lightweight assertion framework (assert_near, assert_true,
     assert_equal, assert_length, assert_finite, assert_positive, run_test)
   - `tests/run_tests.R`: Master runner; exits with status 1 if any tests fail

4. **TDD Red phase**: Wrote 24 failing tests (`tests/test_update_beta.R`) first:
   - T1: Mathematical identities (A_k, B_k, x_k, s_k, second moment)
   - T2: WLS limit (when EL2 = EL², no posterior variance → reduces to WLS)
   - T3: K=1 signal recovery
   - T4: Multi-factor K=5 recovery (signs correct, zero shrunk)
   - T5: Null factor shrinkage (β≈0 when z_no_k is pure noise)
   - T6: Error-in-variables (higher EL2 → more shrinkage)
   - T7: Numerical stability (w=0, extreme weights, n=1, floor test)
   - T8: Gauss-Seidel ordering (updated β₁ propagates to β₂ computation)
   - T9: V2.R consistency (modular function matches V2.R lines 356–361 exactly)

5. **TDD Green phase**: Implemented `code/update_beta.R`:
   - `compute_z_no_k(z, EL, EBeta, k)` — partial working response
   - `update_beta_k(w, z_no_k, EL_k, EL2_k, prior_family, A_floor)` — single-factor update
   - `update_beta_all(w, z, EL, EL2, EBeta, prior_family, A_floor)` — full Gauss-Seidel loop

6. **Fixed 2 test logic errors** during green phase:
   - T4.1: `assert_equal(length(res$details), K)` failed due to `5L ≠ 5` (int vs
     numeric); fixed to `assert_length(res$details, K)`
   - T6.1: Test incorrectly asserted `s_high > s_low`; since s_k = 1/√A_k, larger
     A_k → smaller s_k; fixed assertion to `s_high < s_low`

7. **24/24 tests pass** (`Rscript tests/run_tests.R`)

8. **Created `derivations/qB/qBeta_update_derivation.tex`** (11-page LaTeX document):
   - Full step-by-step derivation with same macros as REVISED.tex
   - Same tcolorbox styles (derivbox/correctionbox/ebnmbox) + new verifybox (green)
     for self-consistency checks and codebox (orange) for math→code mapping table
   - Sections: Setup, Objective, Taylor, Partial Working Response, E_q[L] expansion,
     Coordinate-Ascent, EBNM Problem, Verification Checks, Code Mapping, Algorithm
   - Compiled cleanly: `pdflatex` 3× → 11 pages, no undefined references

**Key decisions:**
- Plain R test framework (no testthat/DESCRIPTION) — no external dependencies
- `update_beta_k` is decoupled from CAVI loop: takes pre-computed vectors
- Returns diagnostic fields (A, B, x, s) to enable white-box test assertions
- Gauss-Seidel propagation in `update_beta_all`: `EBeta_curr[k]` updated immediately

---

### Session 4 (March 11, 2026)

**Demo for β update + to-do setup for L/F/τ modules.**

1. Created `demos/demo_update_beta.R` — 5 interactive scenarios demonstrating the β
   update: anatomy of A_k/B_k, signal recovery K=1, sparsity/shrinkage, error-in-variables,
   multi-factor K=5.  Format: sep/sec/val helpers, pure narrative output (no PASS/FAIL).
2. Updated `PROJECT_STATUS.md` to add `demos/` to the repo map.
3. Updated `code/SupervisedMF_Context.md` with expanded L/F/τ next-step details.
4. Added explicit per-module to-dos for Sessions 5 in this status document.

**Deliverables committed:** `ad80c4e` (pushed to `origin/main`)

---

### Session 5 (March 12, 2026)

**Full implementation of modular updates for q(L), q(F), and q(τ).**

1. **`code/update_L.R`** — Three functions:
   - `compute_R_k(Y, EL, EF, k)`: Partial residual matrix (n×p), shared with F update
   - `update_L_k(...)`: Vector EBNM — A_L[i] and B_L[i] are n-vectors due to
     patient-specific Cox weights W_{ii}; two additive sources (genomics + survival)
   - `update_L_all(...)`: Gauss-Seidel loop over K factors

2. **`code/update_F.R`** — Two functions:
   - `update_F_k(...)`: Pure genomics EBNM — A_F and B_F are p-vectors;
     key property: τ_j **cancels** in x_j = B_F[j]/A_F[j] but does NOT cancel in s_j
   - `update_F_all(...)`: Gauss-Seidel loop over K factors

3. **`code/update_tau.R`** — Three functions:
   - `compute_var_term(EL, EL2, EF, EF2)`: Variance correction (EL2⊗EF2ᵀ − EL²⊗EF²ᵀ)
   - `compute_expected_residual_sq(...)`: R2_bar = naive residual² + Var_Term
   - `update_tau(...)`: Closed-form MLE τ̂_j = n/colSums(R2_bar) + ELBO proxy

4. **Tests** (9 groups each, all passing):
   - `tests/test_update_L.R`: 28 tests (T9 verifies bit-for-bit match with V2.R lines 310–321)
   - `tests/test_update_F.R`: 26 tests (T1.4/T1.5 verify τ cancellation numerically)
   - `tests/test_update_tau.R`: 27 tests (T4 verifies variance correction direction)
   - Updated `tests/run_tests.R`: 105/105 total tests passing

5. **Demos** (5 scenarios each):
   - `demos/demo_update_L.R`: anatomy, signal recovery, genomics vs combined, error-in-variables, K=5
   - `demos/demo_update_F.R`: τ cancellation, sparse recovery, differential shrinkage, sparsity, K=5
   - `demos/demo_update_tau.R`: anatomy, known noise recovery, correction comparison, heteroscedastic, ELBO proxy

6. **LaTeX derivations** (compiled to PDF):
   - `derivations/qL/qL_update_derivation.tex/.pdf`: dual-source vector EBNM
   - `derivations/qF/qF_update_derivation.tex/.pdf`: τ cancellation proof
   - `derivations/qTau/qTau_update_derivation.tex/.pdf`: variance correction + closed-form MLE

**Key design decisions:**
- `compute_R_k` is shared between update_L.R and update_F.R (F sources it from update_L.R)
- τ does not appear in the survival term — F and τ updates are pure genomics
- `update_tau` has NO per-k loop and NO ebnm dependency (only base R)
- ELBO proxy computed inside `update_tau` as a convergence monitor

---

### Session 6 (March 20, 2026)

**Comprehensive documentation: companion docs, roxygen2, and inline comments.**

1. **Created `docs/` directory** with companion documentation for all four modular update scripts:
   - `docs/update_beta.md` — Overview, math background, function reference, test explanations (24 tests), demo explanations (5 demos)
   - `docs/update_L.md` — Dual-source precision/signal, compute_R_k shared dependency, test explanations (28 tests), demo explanations (5 demos)
   - `docs/update_F.md` — τ cancellation property, pure genomics, test explanations (26 tests), demo explanations (5 demos)
   - `docs/update_tau.md` — Closed-form MLE, variance correction, ELBO proxy, test explanations (27 tests), demo explanations (5 demos)

2. **Created `docs/Makefile`** for rendering:
   - `make all` renders all 4 `.md` files to both PDF (via xelatex) and HTML
   - `make clean` removes generated outputs
   - Uses Helvetica/Menlo fonts for Unicode Greek letter support

3. **Added roxygen2 package-readiness tags** to all 11 functions across 4 R scripts:
   - `@export` on every function
   - `@family` tags: `beta_update`, `L_update`, `F_update`, `tau_update`
   - `@seealso` cross-references between modules (e.g., compute_R_k → update_F_k)

4. **Added targeted inline comments** (~3–5 per file) where "why" context was missing:
   - `update_beta.R`: Floor trigger conditions, return list purpose, Gauss-Seidel mutable copy rationale
   - `update_L.R`: A_gen scalar vs A_surv n-vector explanation, inline z_no_k rationale (module independence), dual-source diagnostics
   - `update_F.R`: EL (not EL_curr) passed to R_k because L already updated, sum_EL2_k scalar broadcast
   - `update_tau.R`: Denominator floor logic, no-ebnm dependency note

5. **Verification:** 105/105 tests still passing; all 8 rendered outputs (4 PDF + 4 HTML) generated successfully
