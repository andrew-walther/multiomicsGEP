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
├── code/
│   ├── Supervised_Bayesian_MF.R     ← V1: Original CAVI implementation
│   ├── Supervised_Bayesian_MF_V2.R  ← V2: Corrected implementation (all fixes)
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
│   │   ├── MF_Derivations_UpdateAlgo_2_12_26.pdf  ← Derivation v5 (Feb 12) — REVIEWED
│   │   ├── MF_Derivations_UpdateAlgo_REVISED.tex   ← Corrected derivations (R1-R8)
│   │   └── MF_V2_Companion.tex             ← V2 companion: math ↔ code mapping
│   └── SurvivalMF/
│       ├── SupervisedMF_Likelihood_Estimation_Derivations.pdf
│       └── SurvivalMF_Derivations_Yusha.pdf
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
- [x] REVISED.tex: Corrected derivation document with errata table
- [x] V2.R: Corrected R implementation with all 6 algorithmic fixes (A1-A6)
- [x] V2_Companion.tex: Comprehensive math-to-code companion document
- [x] PROJECT_STATUS.md: Project documentation (this file)
- [x] Simulation validation: V2 runs end-to-end on synthetic data

### Potential Next Steps

- [ ] **Real data application:** Apply V2 to actual genomics + survival datasets
  (TCGA, GEO, etc.) using the `DATA_MODE = "real"` pathway
- [ ] **Additional prior families:** Test `prior_family = "point_laplace"` or
  other EBNM families for different sparsity structures
- [ ] **Cross-validation:** Implement held-out log-likelihood or C-index CV
  for selecting K (number of factors)
- [ ] **Convergence diagnostics:** Track full ELBO (including survival and KL
  terms) instead of genomics-only proxy
- [ ] **Scalability:** Profile and optimise for large p (e.g., p > 10,000)
- [ ] **Manuscript:** Complete `paper/multiomicsGEP_manuscript.qmd`

---

## Session Log

### Session 1 (March 4-5, 2026)

**Derivation review and code correction session.**

1. Reviewed MF_Derivations_UpdateAlgo_2_12_26.pdf against the code
2. Identified 8 derivation errors (R1-R8)
3. Created MF_Derivations_UpdateAlgo_REVISED.tex with all corrections
4. Analysed V1 code against corrected derivations; identified 6 improvements
5. Created Supervised_Bayesian_MF_V2.R with all fixes (A1-A6)
6. Created MF_V2_Companion.tex (math-to-code companion document)
7. Created PROJECT_STATUS.md (this file)

**Key decisions:**
- Orthogonalisation defaults to OFF in V2 (Point-Normal EBNM promotes
  sparsity/distinctness without rotation)
- Taylor refresh defaults to OFF (standard IRLS-within-VI; more stable)
- V1 preserved as-is for reference; V2 is the recommended implementation
