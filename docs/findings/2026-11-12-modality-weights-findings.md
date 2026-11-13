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

## Text for the shared status files

Pending the registered simulation, TCGA decision and external scoring.
