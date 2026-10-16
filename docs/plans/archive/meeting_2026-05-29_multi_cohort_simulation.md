# Meeting Notes: Multi-Cohort Simulation & Model Validation
**Date:** May 29, 2026  
**Attendees:** Andrew, Rashid (Naim), Yusha

---

## Fixed vs. Random Effects

- Rashid prefers handling cohort effects jointly — unclear what information is lost by separating
- Subject-level random effects are for mixed models; the relevant question here is at the study level
- Decision between fixed vs. random effects depends on:
  - **Number of units (studies):** With only 2 cohorts, fixed effects are preferred over random
  - **Degree of heterogeneity:** Random effects could be justified if heterogeneity across studies is high
- Risk of adding cohort covariates: factors compete with shared effects (e.g., basal/classical subtypes) for explaining signal — a factor gets "used up" explaining cohort differences rather than biology

---

## Batch Effects in the Common + Study-Specific Model Structure

- The common + study-specific model structure should automatically capture batch effects:
  - Study-level batch effects are relegated to non-common (study-specific) components
  - Individual-level effects specific to a study are also captured in study-specific components
- Implication: should not need to manually specify study-specific L — batch effects cannot fall into the common portion by model construction

---

## Task: Validate Most Updated Model on Synthetic Data

### Background
- Prior synthetic data validation was done only in a **single-cohort** setting
- Next step: simulate data under a **heterogeneous multi-cohort setting** with both shared and study-specific factors

### Simulation Design Using EBMF

- Apply EBMF to merged multi-cohort genomics matrix to obtain shared and study-specific factors as ground truth
- From the EBMF output (L and F matrices):
  - **Shared factors:** simulate non-zero regression coefficients (Beta ≠ 0) → affect survival
  - **Cohort-specific factors:** simulate Beta = 0 → do not affect survival
- Simulate genomics data from a Gaussian model using the EBMF-derived L and F
- Simulate survival outcomes using shared factors only (study-specific effects should not affect patient survival)
- Run the current model implementation **without** adding dummy cohort indicator variables — test whether the model natively identifies shared vs. study-specific structure
  - If shared/specific factors are already correctly identified, dummy columns may be unnecessary

### Controlling Signal-to-Noise: Study-Specific vs. Shared Effects

- The **F matrix** controls the effect of each individual gene; adjusting the relative magnitude of F values in study-specific vs. shared columns controls the relative strength of each effect type
- Three simulation scenarios:

| Scenario | Description | Interpretation |
|---|---|---|
| **All shared** | No variation across studies; common factors explain all variability; no study-specific effects | Studies are effectively from the same cohort |
| **Nothing shared** | No overlap between studies; all explanatory power in study-specific effects | No shared structure across cohorts |
| **Hybrid** | Some shared structure, some study-specific structure | Cohorts have both similarities and differences |

- Note: inclusion of the survival outcome variable leads to **"anchoring"** — a factor connected to outcome is more likely to be pulled into the shared component

### Simple 2-Cohort Simulation Setup

- Simulate with 2 cohorts
- Use loadings (L) and factors (F) from an EBMF model fit to real data
- Some factors specific to cohort 1, some specific to cohort 2, some shared across both
- Simulate genomics matrix X from a Gaussian model given L, F
- Share the simulation design for L / F / Beta with Rashid and Yusha **before coding**

---

## Benchmarking Against EBMF

- Comparing to EBMF in the paper is a strong way to illustrate the value-add of the Bayesian model
- EBMF on multi-cohort data provides a natural reference for identifying shared vs. study-specific factors

---

## Next Steps

- [ ] **Andrew:** Design simulation mechanism for multi-cohort data (L / F / Beta specification across the 3 scenarios) and share with Rashid and Yusha before implementing
- [ ] **Andrew:** Apply current Bayesian model to simulated multi-cohort data; assess ability to recover shared and study-specific factors and their survival effects
- [ ] **Andrew:** Share updates and questions in chat before next meeting
- [ ] **All:** Next regular meeting — **Monday, June 15, 2026 at 2:00 PM**
