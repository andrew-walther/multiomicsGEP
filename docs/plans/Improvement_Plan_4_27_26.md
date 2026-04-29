# Plan for Improving Multi-omics Integration in DeSurv Benchmark

## Problem Statement
The current "DeSurv Benchmark" simulation, when run in `merged` mode (TCGA RNA-seq + CPTAC Proteomics), consistently results in all survival coefficients (`EBeta`) being shrunk to zero by the Automatic Relevance Determination (ARD) mechanism. This is because the preprocessing strategy (independent top-variable gene selection per cohort followed by strict intersection) amplifies technical batch effects, making them the dominant signal in the combined dataset, which the SSBMF model correctly identifies as non-prognostic for survival.

## Goal
To mitigate the "all B=0" issue by addressing batch effects and enhancing prognostic signals, allowing the SSBMF model to identify relevant factors from merged multi-omics datasets.

## Proposed Strategies & Efforts

### 1. Refine Data Preprocessing and Feature Selection

*   **Option 1: Harmonized Feature Selection:**
    *   **Description:** Instead of independent top-gene selection per cohort, explore methods that select features that are consistently variable *across* cohorts or are known biologically relevant features.
    *   **Effort:** Modify `preprocess_desurv.R` to implement a unified feature selection approach (e.g., selecting top variable genes from a pre-aligned combined dataset or using a union of highly variable genes).
    *   **Impact:** Aims to produce a more biologically coherent gene set for the model.
*   **Option 2: Pre-integration Batch Effect Correction:**
    *   **Description:** Apply dedicated batch effect correction algorithms (e.g., ComBat, RUVseq, limma's `removeBatchEffect`) to individual datasets *before* merging.
    *   **Effort:** Introduce a new step/function within `preprocess_desurv.R` to perform batch correction. Requires careful evaluation of whether prognostic signals are retained.
    *   **Impact:** Directly reduces technical variation, potentially unmasking biological signals.
*   **Option 3: Pathway- or Module-Level Integration:**
    *   **Description:** Instead of gene-level integration, aggregate gene/protein expression to functional units (pathways, gene sets). These higher-level features are often more robust to platform noise.
    *   **Effort:** Develop a new preprocessing module to transform gene-level data into pathway-level data before model fitting.
    *   **Impact:** Shifts focus to more stable biological modules, potentially reducing noise from individual gene/protein measurements.

### 2. Model-Based Batch Effect Handling

*   **Option 1: Explicit Batch Covariate in Model:**
    *   **Description:** Concatenate datasets and explicitly include a batch covariate (e.g., a binary indicator for TCGA vs. CPTAC) in the SSBMF model's formulation.
    *   **Effort:** Requires significant modifications to `fit_modular.R` and its update functions to incorporate batch information into the CAVI updates.
    *   **Impact:** Allows the model to explicitly account for and factor out batch effects without removing the data itself.

### 3. Adjust Model Hyperparameters & Priors

*   **Option 1: Experiment with `lambda`:**
    *   **Description:** The `lambda` parameter in `fit_supervised_mf_modular` scales the survival likelihood. Increasing its value might give more weight to survival signals.
    *   **Effort:** Easily adjustable in `run_ssbmf_benchmark.R`.
    *   **Impact:** Could encourage the model to find prognostic factors even in noisy data, but risks overfitting.
*   **Option 2: Explore Different `prior_beta`:**
    *   **Description:** Test alternative priors for `beta` (e.g., `point_laplace`) within the ARD mechanism to see if they offer different shrinkage properties compared to `point_normal`.
    *   **Effort:** Easily adjustable in `run_ssbmf_benchmark.R`.
    *   **Impact:** Might change the sensitivity of ARD to weak prognostic signals.

## Immediate Next Steps / Experiments

For your meeting, consider proposing these as immediate next steps due to their relatively lower implementation complexity for initial exploration:

1.  **Modify `run_ssbmf_benchmark.R` to increase `lambda` for the "merged" mode:**
    *   Run the "merged" benchmark with `lambda` values like `2.0`, `5.0`, and `10.0` (or higher) to observe the effect on `EBeta` and C-index.
2.  **Modify `run_ssbmf_benchmark.R` to test `prior_beta = "point_laplace"` for the "merged" mode:**
    *   Compare results (`EBeta` values, C-index) with the `point_normal` prior for the merged dataset.
3.  **Investigate Harmonized Feature Selection (Conceptual Discussion):**
    *   Discuss the feasibility and potential approaches for implementing "Harmonized Feature Selection" in `preprocess_desurv.R`, such as taking a union of top variable genes, or selecting top genes from a combined, raw dataset *before* specific cohort-level preprocessing.

## Key Considerations / Challenges

*   **Preserving Biological Signal:** Any batch correction or integration method must be carefully validated to ensure that it removes technical variation without inadvertently eliminating true biological signals related to prognosis.
*   **Computational Cost:** More sophisticated integration methods or model extensions can significantly increase computational time.
*   **Interpretability:** While finding prognostic factors is the goal, maintaining interpretability of these factors and the overall model remains important.
*   **Ground Truth:** For real data, the "true" prognostic factors are unknown, making objective evaluation of integration strategies challenging. Reliance on external validation and comparison to published benchmarks is key.

This plan provides a structured approach to tackle the multi-omics integration challenge and leverage the SSBMF model's capabilities more effectively.