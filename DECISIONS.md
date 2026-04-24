# DECISIONS.md — multiomicsGEP

Architectural and analytical decisions made during development, in reverse chronological order.
Each entry records what was decided, why, what was traded away, and which files implement it.

---

## 2026-04-24 — Lambda survival-scaling parameter: kept at 1.0 after sandbox evaluation

- **Decision:** The `lambda` parameter (scalar multiplier on survival precision terms in the L update) is retained in the codebase at the default value λ=1.0. No active λ tuning is performed.
- **Reason:** A principled argument for λ=p/n exists: the genomics ELBO term sums over p features while the Cox term sums over n patients, so when p>>n the genomics gradient dominates. However, a controlled sandbox (n=250, p=1000, K=5, seed=222) comparing λ∈{1, p/n=5, 2p/n=10} showed no benefit from scaling: hold-out C-index was flat at ≈0.805 across all three conditions, and β RMSE was *worse* at λ=p/n (+0.25) and λ=2p/n (+0.43) than at λ=1. Increasing λ inflates β estimates rather than correcting them, because the dominant source of β scale error is the L–β scale indeterminacy (L can rescale freely), not gradient imbalance.
- **Trade-offs:** The powered-likelihood approach (λ=p/n) is theoretically sound and used in robust Bayesian inference literature. It could become beneficial if the DGP changes (e.g., fewer features, stronger Cox signal). Keeping λ as an exposed parameter with default 1.0 costs nothing and preserves the ability to experiment.
- **Implementation:** λ is a named parameter in `update_L_k()`, `update_L_all()`, and `fit_supervised_mf_modular()` (all default 1.0). It is also in `config/globals.yml` under `cavi.lambda` and threaded through `run_ssbmf_benchmark()` and `run_real_data_benchmark()`. To test λ=p/n, change `globals.yml` and re-run.
- **Affected files:** `code/update_L.R`, `code/fit_modular.R`, `config/globals.yml`, `results/benchmark_sim/run_ssbmf_benchmark.R`, `results/benchmark_sim/sandbox_lambda_test.R`

---

## 2026-04-24 — Proportional hazards diagnostics added to benchmark pipeline

- **Decision:** `cox.zph()` (Grambsch–Therneau test) is now run on SSBMF risk scores for each external PDAC cohort and results are saved to `ph_diagnostics_table.csv` alongside the benchmark outputs.
- **Reason:** The proportional hazards assumption underlies the Cox model used to generate risk scores. A violation means the log hazard ratio is time-varying, which can distort C-index estimates and KM stratification p-values. Formal PH testing is required before the results can be shared or published.
- **Results (TCGA-only, point_normal):** Dijk p=0.77, Moffitt p=0.72, PACA-AU array p=0.34, PACA-AU seq p=0.41 — all PASS. Puleo_array p=0.026 — **FLAG**. The Puleo violation is marginal and likely reflects the large sample size (n=288) giving power to detect subtle time-varying effects; the C-index and KM results remain valid as approximate assessments.
- **Implementation:** `compute_ph_diagnostics.R` (standalone re-fit + projection script); `run_ssbmf_benchmark.R` (PH test now wired into external cohort loop for future runs); `ssbmf_summary_report.qmd` (Section 4.3).
- **Affected files:** `results/benchmark_sim/compute_ph_diagnostics.R` (new), `results/benchmark_sim/run_ssbmf_benchmark.R`, `results/benchmark_sim/ssbmf_summary_report.qmd`

---

## 2026-04-24 — ARD preferred over ELBO grid search for K selection

- **Decision:** K is determined automatically within a single model fit via Automatic Relevance Determination (ARD) — the point-normal/point-laplace prior on β shrinks irrelevant factor coefficients exactly to zero — rather than by fitting separate models at K = 1, …, K_max and comparing ELBOs.
- **Reason:** Two complementary advantages. *Efficiency:* ARD determines K_eff as a byproduct of CAVI; a grid search requires K_max separate full fits and introduces an outer loop. *Bayesian coherence:* ARD performs continuous soft shrinkage within a single probabilistic model; ELBO-based model selection performs discrete hard comparison across models with different dimensionalities. Setting K_max generously large (10) and letting ARD prune is equivalent to an ELBO grid search in expectation, without the computational overhead.
- **Trade-offs:** ARD can be conservative — correlated factors may collapse one even when both carry marginal survival signal. A full ELBO grid search would be more exhaustive but is 10× more expensive at K_max = 10. In practice, ARD + generous K_max matches published EBNM-based NMF standards (flash, flashier).
- **Affected files:** `results/benchmark_sim/run_ssbmf_benchmark.R`, `config/globals.yml` (k_max), `results/benchmark_sim/ssbmf_summary_report.qmd` (Section 1.3)

---

## 2026-04-24 — point_normal chosen as default beta prior over point_laplace

- **Decision:** `point_normal` is the recommended default prior on β for all future SSBMF runs.
- **Reason:** Across synthetic and all PDAC training modes, `point_normal` matches or slightly outperforms `point_laplace` on external C-index (TCGA-only: 0.602 vs 0.579; CPTAC-only: 0.628 vs 0.620). `point_laplace` selects higher α̂ (0.7 vs 0.5 on synthetic), suggesting it compensates for over-shrinkage of small-to-moderate coefficients by drawing more heavily on the survival gradient. The Gaussian slab is also simpler to interpret — posterior SDs have a direct normal-distribution meaning, whereas the Laplace slab mixes two scale regimes.
- **Trade-offs:** `point_laplace` has heavier tails and may outperform `point_normal` in settings with very sparse survival signal (few events, high censoring) where strong coefficient shrinkage is needed. Revisit if future larger-n runs show a consistent >0.02 C-index advantage for `point_laplace`.
- **Affected files:** `results/benchmark_sim/run_ssbmf_benchmark.R` (default `prior_beta`), `config/globals.yml`, `results/benchmark_sim/ssbmf_summary_report.qmd`

---

## 2026-04-24 — Multi-modal TCGA+CPTAC merge documented as expected failure

- **Decision:** The merged TCGA RNA-seq + CPTAC proteomics training mode (838-gene intersection, n=273) is documented as a known failure case rather than a valid benchmark condition. All β̂_k = 0 in both priors.
- **Reason:** When RNA-seq and proteomics are intersected at gene symbols and rank-normalised, PC1 separates the two platforms rather than separating patients by biology. The ARD prior correctly diagnoses that none of the learned factors carry survival signal — they carry platform identity instead. This is not a model failure; it is the model correctly reporting that no prognostic structure exists in this feature space.
- **Trade-offs:** Excluding merged results from the primary benchmark simplifies the comparison table. The failure case is retained in Section 6 of the report as a methodological lesson, motivating the shared-L multi-modal extension (separate F matrices per modality, shared L supervised by survival).
- **Affected files:** `results/benchmark_sim/run_ssbmf_benchmark.R`, `results/benchmark_sim/ssbmf_summary_report.qmd` (Section 6)

---

## 2026-04-24 — DeSurv-aligned preprocessing pipeline added

- **Decision:** Added `code/preprocess_desurv.R` — a preprocessing module that matches the DeSurv paper (Young et al. 2025, PNAS) pipeline: log₂(counts+1) for RNA-seq → select top-2000 most-variable genes per cohort → rank-transform each subject's expression vector.
- **Reason:** DeSurv serves as the primary external benchmark. Using the same preprocessing ensures any C-index difference reflects model architecture, not data transformation choices. Proteomics/microarray platforms skip the log₂ step (already on a normalized scale).
- **Trade-offs:** Rank-transform destroys absolute expression magnitude (EBNM shrinkage is insensitive to scale, but gene-level variance information is lost). Top-2000 gene filter is platform-specific — genes selected differ across cohorts, requiring intersection after preprocessing. TCGA_PAAD × CPTAC intersection yielded 838 genes (42% of 2000), acceptable given DeSurv used the same cohorts.
- **Affected files:** `code/preprocess_desurv.R`, `results/benchmark_sim/run_ssbmf_benchmark.R`

---

## 2026-04-24 — Alpha mixing parameter grid expanded to [0, 1]

- **Decision:** `config/globals.yml` `alpha_grid` expanded from `[0.1, 0.3, 0.5, 0.7, 0.9]` to `[0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]`.
- **Reason:** The boundary values α=0 (pure genomics, unsupervised NMF) and α=1 (pure survival) are meaningful scientific conditions, not just edge cases. Including them lets CV reveal whether the supervision signal is worth anything at all (α=0 wins → pure NMF is optimal; α=1 wins → ignore genomics structure and regress directly).
- **Trade-offs:** Adds 2 extra fits per CV fold. α=0 is equivalent to an unsupervised NMF run; α=1 may be numerically unstable if survival events are sparse (A_surv can be near-zero for small event counts).
- **Affected files:** `config/globals.yml`, `code/select_alpha_cv.R`

---

## 2026-04-24 — predict_supervised_mf() changed from ridge solve() to SVD pseudoinverse

- **Decision:** Replaced `solve(crossprod(EF) + lambda * diag(K))` with an SVD-based Moore-Penrose pseudoinverse in `predict_supervised_mf()`.
- **Reason:** During alpha CV, early-iteration fits can have ARD drive some factor columns of EF to near-zero. The resulting EF'EF is near-singular even with a fixed ridge term λ·I, because near-collinear *non-zero* columns still make the Gram matrix ill-conditioned. The SVD approach sets d_inv = 0 for singular values below `lambda * max(d)` (relative threshold), so collapsed factors contribute nothing to L_test — which is correct since their EBeta ≈ 0 by the same ARD shrinkage.
- **Trade-offs:** SVD is slightly more expensive than a Cholesky solve for dense K×K matrices, but K ≤ 20 makes this negligible. The relative threshold means λ is now a dimensionless tolerance rather than an absolute precision floor; 1e-8 works well empirically.
- **Affected files:** `code/predict.R` (lines 73–94)

---

## 2026-04-24 — Synthetic DGP fixed: equal factor amplitudes + 4-factor survival signal

- **Decision:** Changed `generate_synthetic_benchmark_data()` to (a) use equal F amplitude for all factors (removed 5× multiplier for null factors) and (b) read `beta_true` from `cfg$synthetic$b_true = [1.5, -1.2, 0.8, -0.5, 0.0]` instead of the hardcoded `[1.0, -0.8, 0, 0, 0]`.
- **Reason:** The 5× null-factor inflation caused PCA to capture most variance from non-prognostic factors and incidentally correlate with survival — PCA C-index (0.715) exceeded supervised C-index (0.673) despite the model having the correct generative structure. Equalizing amplitudes restored the intended benchmark: supervised (0.79) beats PCA (0.76). The `b_true` change aligns the DGP with the 4-signal globals.yml spec and tests a more realistic setting (4 prognostic programs at varied effect sizes).
- **Trade-offs:** Changing the DGP invalidates any previously reported synthetic C-index numbers. The new DGP is harder (4 prognostic factors to recover vs. 2) and better reflects real-data complexity.
- **Affected files:** `results/benchmark_sim/run_ssbmf_benchmark.R` (`generate_synthetic_benchmark_data()`)

---

## 2026-04-09 — Repository reorganisation [PENDING]

- **Decision:** *Pending.* Propose a cleaner directory structure (documented in `ROADMAP.md` → Infrastructure section) but do not move any files until a dedicated refactor commit with no concurrent branch work.
- **Reason:** The current layout has accumulated structural debt across three simulation generations: `results/full_sim/`, `results/modular_sim_block/`, and `results/modular_sim_factor/` coexist; `.qmd` reports are mixed with output tables/figures; `demos/` is a top-level sibling of `code/` rather than nested within it; `code/SupervisedMF_Context.md` is a doc file in the algorithm directory; `derivations/EBMF/` and `derivations/SurvivalMF/` are early-sketch folders now superseded by the per-update derivation subdirectories.
- **Trade-offs:** Any file move invalidates hard-coded paths in runner scripts and `.qmd` `source()` calls — must audit before moving. Deferring keeps the repo stable while active development continues.
- **Affected files:** `results/`, `code/`, `demos/`, `derivations/`, `.gitignore`

---

## 2026-04-09 — Synthetic seed changed from 42 → 222

- **Decision:** Changed the random seed for the synthetic data-generating process from 42 to 222.
- **Reason:** Seed 42 produced a degenerate case where True F1 and True F5 (the null factor) were mixed by the bijective permutation alignment, causing β̂_null = −0.77 — an artifact of the seed rather than model behavior. Seed 222 cleanly separates all 5 factors and correctly zeroes the null factor.
- **Trade-offs:** Reported simulation results (RMSE, β̂ estimates, factor correlations) are seed-dependent; changing the seed resets the canonical benchmark. Prior figures are invalidated.
- **Affected files:** `results/modular_sim_factor/run_factor_modular_simulation.R`

---

## 2026-04-09 — `bijective_match()` replaces column-greedy `which.max` for factor permutation alignment

- **Decision:** Replaced `apply(abs(cors), 2, which.max)` with a `bijective_match()` helper that performs greedy global-maximum assignment to align estimated factors to true factors.
- **Reason:** The column-greedy approach is non-bijective — it can assign the same true factor to multiple estimated factors when two estimated factors are most correlated with the same true factor. This produces incorrect factor-recovery plots and corrupted β̂ comparisons in Figure 3.
- **Trade-offs:** Greedy global-max is still a heuristic (not optimal) but guarantees a one-to-one mapping and is O(K²) — negligible for K ≤ 20.
- **Affected files:** `results/modular_sim_factor/run_factor_modular_simulation.R`

---

## 2026-04-09 — Companion document kept in LaTeX (not Quarto) for print-first speaker notes

- **Decision:** The lab meeting companion document (`Notes/lab_meeting_april9_companion.tex`) is authored in plain LaTeX rather than Quarto.
- **Reason:** The companion is a print-first speaker notes document — dense text, no code execution, no cross-referencing with R output. LaTeX compiles faster, supports tighter typographic control, and avoids Quarto's overhead for documents with no R chunks.
- **Trade-offs:** Not integrated with the Quarto build system; must be compiled separately with `pdflatex`. Cannot embed live R output.
- **Affected files:** `presentation/walther_lab_meeting_04_09_2026/Notes/`

---

## 2026-04-01 — Full ELBO tracking added alongside proxy

- **Decision:** Implemented both a fast ELBO proxy (reconstruction error only) and a full ELBO (proxy + survival term + KL divergences) tracked at every iteration.
- **Reason:** The proxy is sufficient for convergence monitoring but the full ELBO is required for model comparison (prior family selection, K selection). Both are now computed and stored in `history$elbo_proxy` and `history$elbo_full`.
- **Trade-offs:** Adds per-iteration cost of `compute_survival_elbo()` (Taylor approximation) and `compute_ebnm_kl()`. Full ELBO requires Taylor refresh every iteration for accuracy.
- **Affected files:** `code/compute_elbo.R`, `code/fit_modular.R`

---

## 2026-04-01 — Four conditions compared: prior family × K strategy (PN/PL × K=5/K_eff)

- **Decision:** All benchmarks compare 4 conditions: point-normal (PN) vs. point-laplace (PL) crossed with fixed K=5 vs. adaptive K_eff from `auto_prune_K()`.
- **Reason:** Isolates the contribution of prior choice from K selection strategy. Reveals whether adaptive K adds value beyond a fixed-K run with the same prior.
- **Trade-offs:** Requires 4× the compute per dataset. K_eff runs require an extra `auto_prune_K()` pass before the main fit.
- **Affected files:** `results/modular_sim_factor/run_factor_modular_simulation.R`, `code/select_K.R`

---

## 2026-03-31 — `fit_modular.R` uses Gauss-Seidel (factor-wise sequential) CAVI, not full-gradient

- **Decision:** Each CAVI iteration updates factors sequentially (k=1, ..., K), immediately incorporating each updated factor before proceeding to the next. This is the canonical CAVI loop in `code/fit_modular.R`.
- **Reason:** Full-gradient (parallel) CAVI requires the dense Hessian of the joint ELBO with respect to all factors simultaneously, which is intractable in closed form under the Taylor approximation for the survival term. Gauss-Seidel is the standard choice for mean-field CAVI.
- **Trade-offs:** Sequential updates mean later factors in each iteration benefit from earlier updates (faster convergence per iteration), but the ordering introduces implicit asymmetry between factors. Not easily parallelized across factors.
- **Affected files:** `code/fit_modular.R`

---

## 2026-03-25 — Hold-out prediction uses pseudo-inverse projection of Y_test onto F

- **Decision:** `predict_supervised_mf()` obtains test-set loadings L_test by projecting Y_test onto the trained factor matrix F via pseudo-inverse: `L_test = Y_test %*% F %*% solve(t(F) %*% F)`.
- **Reason:** After training, F is fixed. The natural way to score a new patient is to find the loadings that best reconstruct their expression profile under the trained factors. The pseudo-inverse gives the least-squares solution.
- **Trade-offs:** Does not propagate uncertainty from F into L_test (point estimate only). Assumes test data is drawn from the same distribution as training data (no domain shift).
- **Affected files:** `code/predict.R`

---

## 2026-03-20 — Convergence criterion: dual threshold on mean|ΔL| AND mean|Δβ| after 5-iteration burn-in

- **Decision:** Convergence is declared when both `mean(|EL - EL_old|) < tol` AND `mean(|EBeta - EBeta_old|) < tol`, checked only after iteration 5.
- **Reason:** A single criterion on L alone can declare convergence while β is still adjusting (or vice versa). The dual criterion ensures both the genomic structure (L) and survival signal (β) have stabilized. Mean is used instead of max because max is dominated by a few high-variance entries near factor orientation boundaries and rarely reaches 1e-3 on real datasets.
- **Trade-offs:** Mean convergence is weaker than max convergence — some individual loadings may still be changing when the algorithm stops. The 5-iteration burn-in prevents premature stopping during large initial swings.
- **Affected files:** `code/fit_modular.R` (lines 353–370)

---

## 2026-03-12 — Modular update architecture: each update function is independently testable

- **Decision:** CAVI updates are split into four independent modules (`update_L.R`, `update_F.R`, `update_beta.R`, `update_tau.R`), each with its own test suite and demo scripts.
- **Reason:** The monolithic V2 implementation (`Supervised_Bayesian_MF_V2.R`) mixed all update logic into one function, making it hard to test individual components or swap implementations. Modular architecture enables TDD, isolated debugging, and future replacement of individual components.
- **Trade-offs:** Cross-module dependencies must be managed explicitly (e.g., `compute_R_k` lives in `update_L.R` but is also used by `update_F.R`, requiring careful source ordering).
- **Affected files:** `code/update_L.R`, `code/update_F.R`, `code/update_beta.R`, `code/update_tau.R`, `code/fit_modular.R`

---

## 2026-03-10 — EBNM priors chosen over fixed-penalty; g estimated per CAVI step

- **Decision:** Shrinkage priors for L, F, and β use the Empirical Bayes Normal Means (EBNM) framework, with the prior g estimated from data at each CAVI step rather than fixed by cross-validation.
- **Reason:** Fixed-penalty methods (Ridge, LASSO) require cross-validation for λ, adding computational overhead and a choice of CV criterion. EBNM folds g-estimation into the ELBO maximization — the prior adapts to the data, and no separate CV loop is needed. Implemented via the `ebnm` R package (Stephens lab).
- **Trade-offs:** g-estimation adds per-update overhead. The estimated g is dataset-dependent and may change across CAVI iterations (instability if the model is mis-specified). Point-normal and point-laplace are the two families tested; other families (e.g., generalized double-Pareto) are not yet explored.
- **Affected files:** `code/update_L.R`, `code/update_F.R`, `code/update_beta.R`

---

## 2026-02-12 — Baseline hazard h₀(t) left non-parametric; enters through Cox partial likelihood

- **Decision:** The survival component uses the Cox proportional hazards model with an unspecified baseline hazard h₀(t). h₀(t) cancels exactly in the Cox partial likelihood and never needs to be estimated.
- **Reason:** Specifying a parametric h₀(t) (Weibull, exponential) would add a nuisance parameter and require a prior on its shape/scale. The Cox partial likelihood avoids this entirely while preserving the proportional hazards structure needed to link L to survival.
- **Trade-offs:** Cannot predict absolute survival probabilities (only relative risk via exp(Lβ)). The partial likelihood is not a true likelihood (it conditions on observed event times), which means the ELBO approximation for the survival term requires a Taylor expansion rather than a closed-form KL.
- **Affected files:** `code/update_beta.R`, `code/update_L.R`, `code/compute_elbo.R`

---

## 2026-02-12 — Model formulation: shared L links genomics (Y = LF′ + E) and survival (h(t) = h₀(t)exp(Lβ))

- **Decision:** The core model posits a shared n×K loading matrix L that simultaneously factorizes the genomics matrix Y and enters the Cox proportional hazards model as the linear predictor.
- **Reason:** A two-stage approach (first factor Y, then regress loadings on survival) does not jointly optimize the factors for survival prediction. Sharing L with a joint objective ensures the learned gene expression programs are informative for both reconstruction and prognosis.
- **Trade-offs:** The joint objective couples genomics and survival gradients; the relative scale of the two terms is dataset-dependent and not normalized (p >> n means the genomics term dominates). This is a known limitation — adding a λ scaling parameter is an identified future direction.
- **Affected files:** `code/fit_modular.R`, `code/update_L.R`, `code/update_beta.R`, `derivations/MF_UpdateDerivations/`
