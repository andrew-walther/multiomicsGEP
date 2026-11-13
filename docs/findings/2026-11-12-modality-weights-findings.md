# Reconstruction weights for expression and methylation

Date: 2026-11-12. Author: Andrew Walther.

## Scope and progress

The registered study is [the modality-weights prompt](../prompts/2026-11-12-codex-d-modality-weights-prompt.md), revised in `9c4cac23`.
The aim is recovery of joint and modality-exclusive programs with a large CpG set.
Prediction is secondary. Work proceeds on main, with at most two CPU workers.
Task C inputs are read-only. Patient fits and scores remain under ignored `data/`.

## Derivation and independent review

The weighted objective multiplies each entire Gaussian reconstruction log-likelihood
by its modality weight. All priors and KL terms retain weight one. L and the
reconstruction part of F use effective precision w_m τ_mj. Survival updates and
Z = Σ_m Y_m F_m are unchanged. The positive weight cancels from the featurewise
τ point estimate, n / Σ_i E[R²_mij]. The objective includes the weighted log τ term.
Common-program rescaling preserves the objective when priors transform with the
posterior. Data-block rescaling changes the relative survival projections.

The updated derivation verification passed 109/109 checks. A fresh independent
statistical review ran the same checks and found no blocking mathematical errors.
A stale sentence was corrected following review. The tests include brute-force
coordinate coefficients, numerical τ maximization, canonicalization and a negative
control for block-scale dependence of Z.

For the revised zero-survival scaling check, corresponding starts and uncanonicalized
reconstruction coordinates are required. Projection normalization and fresh SVD
starts themselves depend on block scale. Positive rescaling before training-fold
z-scoring is tested separately for the supervised fit.

## Fitter implementation and verification

`code/modality_weights_yfb.R` validates named positive powers and computes both
training-only rules. The leading eigenvalue uses the patient Gram matrix, avoiding
a CpG-by-CpG allocation. Constant training features fail explicitly.
`fit_ssbmf_yfb()` sends w_m τ_mj to reconstruction updates while retaining the
unweighted τ estimate and projection. The objective weights both Gaussian terms.
NULL and all-one weights bypass extra arithmetic.

Andrew approved `F_prior_by_group = TRUE` as a separate option on November 12.
FALSE pools the F prior while retaining feature-group labels for the weights.
The all-one pooled path reproduces the ungrouped posterior and history bit for bit.
Weighted pruning, posterior restarts and burn-in/warm starts are rejected because
this study does not implement or evaluate them. Defaults remain unchanged.

The focused run passed 25/25 tests, including the prior feature-group checks.
The tests failed before the new API was implemented. The weighted objective was
checked using numerical integration of loading KL under the returned prior, not
only its cached KL. A separately recomputed objective also matched. Full suite:
865 passed, zero failed. The existing bounded-β profile-likelihood fallback warning
was retained. The revised zero-survival check uses a controlled reconstruction
coordinate sweep with corresponding transformed moments, separately from the
supervised raw-unit preprocessing test.

## Text for the shared status files

Pending the registered simulation, TCGA decision and external scoring.
