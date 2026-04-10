# DECISIONS.md — multiomicsGEP

Architectural and analytical decisions made during development, in reverse chronological order.
Each entry records what was decided, why, what was traded away, and which files implement it.

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
