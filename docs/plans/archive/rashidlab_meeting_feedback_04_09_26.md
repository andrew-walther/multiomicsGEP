# Lab Meeting Feedback — April 9, 2026
## Supervised Bayesian Matrix Factorization (PDAC Project)

---

## Summary

Andrew presented preliminary work on supervised Bayesian matrix factorization for pancreatic ductal adenocarcinoma (PDAC) datasets. The meeting surfaced methodological feedback across model specification, validation strategy, and presentation clarity. Rashid and the group raised concerns about the scaling of survival vs. genomics components, C-index calculation correctness, prior/hyperparameter documentation, and the need for a clearer validation pipeline using merged cohort training.

---

## Some Key Notes/Feedback from Andrew's Presentation

- **Slide 5 — Motivation/Benefits:** Run both simulation and real-world evaluation on the top three points from slide 5. The goal is to actively demonstrate the benefits of this approach, not just assert them.

- **Survival Component — Use F^T Y instead of L:** The survival model should use F^T Y rather than L in the survival component. (Note: Amber is sending relevant papers on F vs. L usage in similar models.)

- **Slide 6 — Prior and Hyperparameter Specification:** Explicitly specify the priors for L, F, and B, and document how hyperparameters are set. Also clarify how the baseline hazard is represented — **parametrically or non-parametrically** — as this is currently unclear.

- **Slide 7 — Scaling of Genomics vs. Survival Components:** Think carefully about the scaling difference between L_genomics and L_cox. The gradients may be dominated by the component with the larger expectation, which could distort learning. This likely needs explicit balancing or re-scaling.

- **Slide 9 — Convergence Criteria:** Specify the convergence criterion being used for vector/matrix quantities (e.g., L-hat). The current presentation does not make this explicit.

- **Slide 9 — Define q:** The variational distribution q is referenced but not defined. Add a clear definition of q and its role in the ELBO/variational terms.

- **Cross-Dataset Normalization:** Normalize across all datasets before splitting into training and validation. This should precede any training/validation comparisons to avoid data leakage or scale artifacts.

- **Slide 11, Figure 3 — Factor Interpretation:** The difference between the first three factors and the last two factors in Figure 3 may be illustrating the scaling imbalance noted above (Slide 7 point). Investigate whether this factor split reflects the genomics vs. survival component dominance issue.

---

## Andrew's Action Items (from Next Steps)

1. **Baseline hazard handling:** Review code and model details to confirm correct handling of the baseline hazard — clarify whether it is parametric or non-parametric, and make this explicit in both the code and presentation.

2. **Model component documentation:** Define and document all model components clearly — priors, Q, P, expectations — especially in the context of the ELBO and variational terms. This should appear in both the presentation and the codebase.

3. **C-index < 0.5 investigation:** Investigate and resolve the issue of C-index values less than 0.5 in model validation. Check the direction of risk scores and the computation details — a flipped risk score direction is a common source of this issue.

4. **Merged cohort training:** Train the model on merged CPTAC + TCGA datasets and use the remaining datasets as fixed validation sets. Recompute all validation metrics under this setup for comparison with current results.

5. **Survival vs. genomics tuning:** Consider adding a tuning parameter to control the relative contribution of the survival component versus the genomics component. Also review the impact of prior choice and factor selection (K) on model performance.

---

## Incoming Resources

- **Amber → Andrew:** Amber will send papers on the use of F (features matrix) vs. L (loadings matrix) in the survival component of similar models. Review these before updating the model formulation.

---

## Broader Context (AI Tool Usage — Applicable to All)

Rashid led a separate discussion on best practices for using AI tools (e.g., Claude) in research:

- Sketch plans before prompting; don't let AI drive the design.
- Review generated code line by line before committing.
- Understand the underlying algorithm — you must be able to explain and defend it.
- Manage token usage efficiently.
- Maintain ownership of the research; AI is a collaborator, not the author.

---

## Slide-by-Slide Carry-Over Guidance (for New Deck)

For each slide from the old deck that might carry content forward, the table below summarizes the verdict and what needs to change before reuse. Slides not listed had no specific feedback and can be ported or dropped at discretion.

| Old Slide | Content | Verdict | Required Changes Before Reuse |
|-----------|---------|---------|-------------------------------|
| 5 | Motivation / claimed benefits of approach | **Revise** | Add simulation + real-world results that actively *demonstrate* each of the top 3 points rather than asserting them. |
| 6 | Model specification | **Revise** | (1) Add explicit prior distributions for L, F, and B with hyperparameter values/ranges. (2) State clearly whether baseline hazard is parametric or non-parametric. (3) Update survival component formula to use F^T Y instead of L. |
| 7 | Genomics vs. survival component loss / update equations | **Revise** | Add explicit discussion of scaling difference between L_genomics and L_cox. Document how gradient domination is addressed (e.g., re-scaling, tuning weight). |
| 9 | Algorithm / ELBO / variational updates | **Revise** | (1) Define q (variational distribution) explicitly. (2) Specify convergence criterion for vector/matrix quantities (e.g., L-hat). (3) Define all terms: Q, P, expectations. |
| 11 | Figure 3 — factor visualization | **Revise or replace** | Rerun after resolving the scaling issue (Slide 7). The first-three vs. last-two factor split likely reflects component imbalance and should be re-examined once scaling is corrected. Recomputed metrics should use the merged CPTAC + TCGA training setup. |
| Any validation/results slides | C-index and dataset comparisons | **Hold — recompute first** | Do not carry over existing C-index values. (1) Investigate and fix C-index < 0.5 (check risk score direction). (2) Retrain on merged CPTAC + TCGA; recompute metrics with remaining datasets as fixed validation sets. (3) Normalize across all datasets before any split. |

### New Slides Needed in the Next Deck

The following content was flagged as missing or underdeveloped and should be built fresh:

- **Simulation study slide:** Demonstrate model benefits (from slide 5 motivation points) on synthetic data.
- **Prior/hyperparameter specification slide:** Dedicated slide or appendix showing the full generative model with all priors, hyperparameters, and the baseline hazard formulation.
- **Scaling / component weighting slide:** Explain the L_genomics vs. L_cox scaling problem and the proposed solution (tuning parameter, re-scaling strategy, or both).
- **Updated validation results slide:** Post-retraining C-index and other metrics under the merged cohort training setup, with correct risk score direction confirmed.

---

*Notes synthesized from lab meeting on April 9, 2026. Original notes compiled by meeting recorder; Andrew-specific items extracted and organized for presentation revision use.*
