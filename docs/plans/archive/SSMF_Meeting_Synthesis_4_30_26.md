# SSMF Meeting Synthesis: Supervised Survival Bayesian Matrix Factorization

**Date**: [Meeting Date]  
**Participants**: Andrew, Rashid  
**Key Context**: Addressing merged data cohort effects, initialization strategy, and factor validation methodology

---

## 🎯 Executive Summary

The meeting clarified that **previous "promising" validation results originated from TCGA-only models, not merged data**. The core issue—survival factors being heavily penalized to zero in merged datasets—remains unresolved. The session established a revised validation strategy prioritizing **factor projection onto new datasets** over direct coefficient interpretation, and laid out a detailed initialization and optimization roadmap based on Amber's DeSurv methodology.

---

## 📋 Current Status & Key Findings

### Pre-processing & Heatmaps
- ✅ Pre-processing reordered; successfully generated 2,000 gene set
- ❌ Gene selection did not improve survival-related factor identification
- ✅ Factor heatmaps generated for each prior (Normal/PN/PL) × training set (TCGA/CPTAC/Merged)

### Cohort-Associated Patterns & Validation Issue
- Unsupervised EBMF identified 20 factors; 5 showed overall survival association
- **Critical Finding**: "Good" external validation C-index values came from **TCGA-only model**, not merged data
- Merged cohort analysis: survival factors penalized to zero (fundamental blocker)
- Cohort-associated batch effects remain apparent in merged data

### Matrix Notation Clarification
- **DeSurv notation**: Y = WH (W = factor matrix, H = subject loadings)
- **SSMF notation**: Y = LF (L = subject loadings matrix, F = factor matrix)
- *All following discussion uses L/F notation for SSMF context*

---

## 🔑 Core Technical Insights

### 1. Validation Philosophy Shift
**Old Approach** ❌  
- Interpret individual regression coefficients (β) from joint model
- Problem: instability + factor correlation makes direct interpretation meaningless

**New Approach** ✅  
- Project learned **factor matrix (F)** onto new, independent datasets
- Validate via risk score: $\eta = Y_{new} F^T B$ (or equivalent projection)
- Assess whether learned factors are truly survival-associated in held-out data
- **Key insight**: Coefficients themselves are less important than demonstrating learned factors correlate with survival in external cohorts

### 2. Beta Prior as Ridge Penalty
- Prior on β acts as **ridge penalty** that activates/inactivates factors
- Even factors "turned off" in survival model may retain biological relevance
- Factor activation is **not** the same as biological relevance

### 3. K (Number of Factors) Dependency
- **Betas are highly dependent on K**; cannot be treated as a nuisance parameter
- Two complementary strategies:

| Strategy | Approach | When to Use |
|----------|----------|------------|
| **Sparsity via Prior** | Fix K large (e.g., K=20); let β prior penalize inactive factors | Computational efficiency; assume sparse set of truly active factors |
| **Cross-Validation** | Loop over range of K values; select "true" K via CV | Principled factor selection; unknown true dimensionality |

### 4. Scale & Initialization Dominance Problem
- NMF term can overwhelm survival term if scales not balanced
- Survival factors get penalized to zero when NMF dominance is high
- **Solution**: Strict normalization constraints during initialization

---

## 🔧 Implementation Roadmap

### Phase 1: Data Harmonization (Foundation)

#### Rank & Quantile Normalization
```
Goal: Harmonize cross-platform datasets before model fitting

For each training dataset (TCGA, CPTAC, Merged):
  1. Apply rank normalization (critical for cross-platform data)
  2. Optionally apply quantile normalization for additional harmonization
  3. Document normalization parameters for reproducibility
  4. Store normalized datasets separately to avoid reprocessing
```

**Why**: Cross-platform data reside on drastically different scales. Without normalization, scale disparities skew the model and push survival factors to zero.

---

### Phase 2: Matrix Initialization & Normalization (Core)

#### Amber's DeSurv Initialization Sequence
1. **Initialize β = 0**
   - Start with no survival signal
   - Let block coordinate updates incorporate survival gradually

2. **Normalize Data Matrix (Y)**
   - Divide each column by its sum: $Y_{norm} = Y / \text{colSums}(Y)$
   - Bounds all values to [0, 1] range
   - Prevents NMF from dominating survival term

3. **Normalize Factor Matrices**
   - One of {L, F} must have weights in each column sum to 1
   - Maintain non-negativity constraints on both L and F
   - *Clarification needed*: Confirm which matrix (L or F) gets normalized and verify exact constraint formulation from Amber's supplement

4. **Update Order** (Block Coordinate Scheme)
   - Cycle: F → L → B (repeat until convergence)
   - Fix other matrices while updating one
   - Prevents local minima in non-convex optimization

---

### Phase 3: Multi-Initialization Optimization (Robustness)

#### Random Restart Strategy
```
For each K value:
  For i = 1 to N_init (30-100 random initializations):
    - Re-sample random starting points for L, F
    - Run block coordinate updates F → L → B until convergence
    - Record:
      * Final objective function value
      * Final C-index (external validation set)
      * Model parameters (L, F, B)
  
  Select best model by:
    - Option A: Highest objective function
    - Option B: Highest C-index (clinical validation)
    - *Decision needed*: Objective function vs. C-index for selection criterion
```

**Why**: Single optimization run prone to local minima. 30-100 restarts dramatically increases chance of finding global optimum.

---

### Phase 4: Factor Projection Validation (Assessment)

#### Implementation Steps
1. **Learn factor matrix F** on training data (merged TCGA+CPTAC)
2. **Project F onto external cohorts**:
   - TCGA-only test set: $Y_{test} F^T$
   - CPTAC-only test set: $Y_{test} F^T$
   - Other external cohorts if available
3. **Compute risk scores** using projected factors:
   - $\eta = Y_{new} F^T B$
   - Use B from training phase
4. **Validate** via C-index on each held-out cohort
   - Assess whether learned factors generalize
   - Compare vs. baseline (unsupervised EBMF, Cox regression, etc.)

**Why**: External validation on projections shows whether factors capture genuine survival signal vs. cohort-specific artifacts.

---

## 📝 Pending Technical Clarifications

These items require input from Rashid's DeSurv supplement or manuscript before full implementation:

1. **Matrix Normalization Detail**
   - Which matrix (L or F) gets column-wise normalization?
   - Exact constraint formulation: sum-to-1 on each column only, or additional bounds?
   - How do constraints interact with non-negativity?

2. **Objective Function Selection**
   - When selecting best model from N_init restarts:
     - Use **objective function** (likelihood + penalties)?
     - Use **C-index** (external validation metric)?
     - Ensemble/voting scheme across both?

3. **Normalization Scope**
   - Normalize all three datasets (TCGA, CPTAC, Merged) separately?
   - Or normalize after merging?
   - How to handle batch effects introduced by preprocessing?

4. **K Selection Strategy Priority**
   - Should we prioritize fixed-K with sparsity prior, or CV-based K selection?
   - Computational budget vs. principled selection tradeoff?

---

## ✅ Action Items

### Andrew (Immediate)
- [ ] **Review Amber's DeSurv manuscript & supplement** (especially initialization & normalization sections)
- [ ] **Implement rank normalization pipeline**
  - [ ] Develop function(s) for rank/quantile normalization
  - [ ] Apply to TCGA, CPTAC, and Merged datasets
  - [ ] Document normalization parameters and preserve normalized datasets
- [ ] **Implement matrix normalization & initialization**
  - [ ] Build initialization function following Amber's DeSurv sequence
  - [ ] Verify column-wise normalization constraints on L/F
  - [ ] Add validation checks to ensure constraints maintained during updates
- [ ] **Implement multi-initialization wrapper**
  - [ ] Loop over N_init = 30-100 random starts
  - [ ] Run block coordinate updates (F → L → B) per initialization
  - [ ] Store all models; select best by objective function (or clarify selection criterion with Rashid)
  - [ ] Document convergence behavior and objective function landscapes
- [ ] **Recheck & verify recent changes**
  - [ ] Confirm update loop reordering effects
  - [ ] Document which changes improved vs. worsened survival signal
  - [ ] Create reproducible script for full pipeline end-to-end
- [ ] **Plan & design factor projection validation**
  - [ ] Outline how F will be applied to external cohorts
  - [ ] Specify risk score computation: $\eta = Y_{new} F^T B$
  - [ ] Design C-index comparison framework (learned F vs. baselines)
- [ ] **Clarify pending technical items** with Rashid before proceeding

### Rashid
- [ ] Share latest Amber DeSurv manuscript & supplement (if not already received)
- [ ] Clarify matrix normalization specifics (which matrix, exact constraints)
- [ ] Advise on objective function vs. C-index selection criterion for multi-start optimization

---

## 📚 Reference Materials

### Key Concepts
- **Block Coordinate Updating**: Cycle through updating one matrix at a time while fixing others; improves convergence vs. simultaneous updates
- **Ridge Penalty via Prior**: Putting a shrinkage prior on β allows sparsity without explicit variable selection
- **Factor Projection**: Applying learned F to new data (Y_new F^T) validates generalizability without re-fitting
- **Data Harmonization**: Rank/quantile normalization removes scale disparities in cross-platform genomic data

### Files to Reference
- Amber's DeSurv manuscript & supplement (initialization, normalization, update schemes)
- Current SSMF codebase in `~/GithubProjects/multiomicsGEP/`
- PDAC data directory: `~/OneDrive/UNC Dissertation (Liu)/PDAC_data/`

---

## 🚀 Next Immediate Steps (Priority Order)

1. **Acquire & review Amber's DeSurv materials** → resolve technical clarifications
2. **Build & test rank normalization** → apply to all three training sets
3. **Implement matrix initialization & constraints** → ensure NMF/survival balance
4. **Develop multi-start optimization loop** → N_init = 30-100 with full pipeline
5. **Create factor projection validation** → external cohort testing framework
6. **Document & git commit** → clean, logical commits with clear messages

---

## 💡 Key Takeaways for Code Implementation

✅ **Do:**
- Use block coordinate updates (F → L → B cycle)
- Initialize β = 0; build in survival signal gradually
- Normalize Y and one of {L, F} during initialization
- Run 30-100 random initializations per K
- Validate via factor projection on held-out cohorts
- Document normalization parameters for reproducibility

❌ **Don't:**
- Interpret individual coefficient values directly (due to multicollinearity)
- Rely on single optimization run (prone to local minima)
- Skip cross-platform normalization (scales matter critically)
- Assume learned factors are survival-related without external validation
- Treat K as fixed without exploring sensitivity

---

**Last Updated**: [Meeting Date]  
**Next Review**: After implementing Phases 1–2 and clarifying technical details with Rashid
