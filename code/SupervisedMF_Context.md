# Math ↔ code reference

A compact map from the model's notation to the R implementation. For the file
layout see [`README.md`](README.md); for complete derivations see
`derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf` (LB) and
`derivations/cox_on_YF/` (YFB).

## Model

```
Y (n×p) = L (n×K) F' (K×p) + E,   E_ij ~ N(0, 1/τ_j)
h(t_i)  = h₀(t_i) exp(η_i)
```

| Model | Linear predictor | Fitter | Prediction for new patients |
|-------|------------------|--------|-----------------------------|
| LB  | η_i = l_i β | `fit_supervised_mf_modular()` in `fit_modular.R` | `predict_supervised_mf()` in `predict.R` (projects Y_new onto F, then applies β) |
| YFB | η_i = (y_i F) β | `fit_cox_on_yf()` in `fit_cox_on_yf.R` (recommended) | `predict_cox_on_yf()` in `predict_cox_on_yf.R`: Y_new F β, the same formula used in training |

Inference is mean-field CAVI, q(L, F, β) = q(L) q(F) q(β), with τ updated in
closed form. Each L, F and β coordinate update is an empirical Bayes normal
means (EBNM) problem: given pseudo-observations x with standard errors s, `ebnm`
estimates the prior and returns posterior moments. The Cox partial likelihood is
replaced in each outer iteration by its second-order Taylor expansion around the
current η, giving per-patient weights W_ii and a working response z_i. Updates
within an iteration are Gauss-Seidel: each factor's update uses the factors
already updated in that iteration.

## Variables

| R variable | Symbol | Dimension | Meaning |
|------------|--------|-----------|---------|
| `EL`, `EL2` | E[l_ik], E[l_ik²] | n × K | Posterior first and second moments of patient scores |
| `EF`, `EF2` | E[f_jk], E[f_jk²] | p × K | Posterior moments of feature loadings |
| `EBeta`, `EBeta2` | E[β_k], E[β_k²] | K | Posterior moments of survival coefficients |
| `Tau` | τ_j | p | Feature-specific noise precision |
| `w` | W_ii | n | Cox working weights (negative diagonal Hessian) |
| `z` | z_i | n | Cox working response, η̂_i + u_i / W_ii |
| `R_k` | R⁻ᵏ | n × p | Residual with every factor except k removed |
| `z_no_k` | z⁻ᵏ_i | n | Working response with every factor except k removed |
| `ZF` | Z = YF | n × K | YFB program scores (observed projections) |
| `A_*`, `B_*` | A, B | | EBNM precision and signal: x = B/A, s = 1/√A |

E[l²] ≠ E[l]² whenever the posterior has variance. Using second moments in the
precisions is what propagates uncertainty between updates.

## LB updates (`update_L.R`, `update_F.R`, `update_beta.R`, `update_tau.R`)

| Update | Precision A | Signal B | R |
|--------|-------------|----------|---|
| q(l_ik) | Σ_j τ_j E[f_jk²] + W_ii E[β_k²] | Σ_j τ_j R⁻ᵏ_ij E[f_jk] + W_ii z⁻ᵏ_i E[β_k] | `update_L_k()` |
| q(f_jk) | τ_j Σ_i E[l_ik²] | τ_j Σ_i R⁻ᵏ_ij E[l_ik] | `update_F_k()` |
| q(β_k)  | Σ_i W_ii E[l_ik²] | Σ_i W_ii z⁻ᵏ_i E[l_ik] | `update_beta_k()` |

- In q(f_jk), τ_j cancels from x = B/A but not from s.
- In q(β_k), the precision uses E[l²] rather than E[l]², so uncertain scores
  shrink β (an errors-in-variables correction).
- q(τ): τ_j = n / Σ_i E[(y_ij − l_i f_j')²], where the expected squared residual
  adds Σ_k (E[l²]E[f²] − E[l]²E[f]²) to the squared mean residual
  (`update_tau()`, closed form, no EBNM).

## YFB updates (`update_L_surv_YFB.R`, `update_F_surv_YFB.R`)

L no longer enters the survival model, and β multiplies the observed scores
Z = YF.

| Update | Precision A | Signal B | R |
|--------|-------------|----------|---|
| q(l_ik) | Σ_j τ_j E[f_jk²] | Σ_j τ_j R⁻ᵏ_ij E[f_jk] | `update_L_surv_YFB_k()` |
| q(f_jk), α_F = 0 | τ_j Σ_i E[l_ik²] | τ_j Σ_i R⁻ᵏ_ij E[l_ik] | `update_F_surv_YFB_k()` |
| q(β_k)  | Σ_i W_ii Z_ik² | Σ_i W_ii z⁻ᵏ_i Z_ik | `update_beta_k()` with Z in place of E[L] |

The F update can also take a survival term weighted by α_F, which adds
α_F E[β_k²] Σ_i W_ii y_ij² to A and α_F E[β_k] (Y'(W z⁻ᵏ))_j to B, with the
genomics terms scaled by (1 − α_F). The default α_F = 0 learns F from the
molecular data alone because α_F > 0 caused a feedback instability (see
`DECISIONS.md`, 2026-04-30). q(τ) is shared with LB.

## Fitter outputs

Both fitters return `EL`, `EL2`, `EF`, `EF2`, `EBeta`, `EBeta2`, `Tau` and
`history` (iteration count, convergence flag, objective traces). `fit_cox_on_yf()`
also returns `EF_norms`, which `predict_cox_on_yf()` needs to score new patients
on the training scale.

## Numerical safeguards

- EBNM precisions are floored (`A_floor = 1e-10`) so that s = 1/√A stays finite.
- The τ denominator is floored at n × 10⁻⁸.
- The Cox Taylor weights W_ii are floored at 10⁻⁶.
