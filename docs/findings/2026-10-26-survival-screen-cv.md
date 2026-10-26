# Survival-associated gene screening: training cross-validation

Date: 2026-10-26. Design: `docs/plans/2026-10-16-feature-selection-experiment-plan.md` (Tasks
E1–E4) with the defaults approved on 2026-10-26. Runner:
`results/feature_selection/run_feature_selection_comparison.R` (modes `manifest`, `cv`, `summary`,
`agreement`); helpers and their tests in the same folder. Outputs:
`results/feature_selection/outputs/` (fits git-ignored). No external cohort was loaded.

## Question

Does screening genes for association with survival let the expression-only model reach the same
prediction and the same programs with fewer genes than the current unsupervised mean/variance
screen?

## Design

- Data: the 273 training patients (TCGA-PAAD and CPTAC), with the benchmark's preprocessing
  (historical double log, per-cohort z-scores computed once on all training patients before the
  folds, so CV is conditional on that scaling). The single-log check found this preprocessing
  choice does not change external performance.
- Arms, all at K = 7, α_F = 0, seed 42, the same five folds:
  - **U2064**: the benchmark's 2,064 genes (combined mean/variance rank, top 3,000 per cohort,
    intersection).
  - **U(N)**: the same unsupervised screen with a smaller per-cohort count, giving about N shared
    genes (nested in U2064). Sizes 251, 500 and 1,000 (ties give 251 for N = 250).
  - **S(N)**: the 2,064 benchmark genes ranked by the Wald p-value of
    `coxph(Surv(time, status) ~ gene + strata(cohort))` on the fold's training patients only, top
    N, with N equal to U(N)'s count so each pair is matched in size.
- Risk direction from each fit's training patients; held-out Harrell's C. 42 fits (7 arm-size
  combinations × 5 folds + full training); every fit converged; no errors. The per-gene Cox fit
  failed to converge for one gene in fold 2 (assigned p = 1, ranked last, counted in
  `manifest.csv`) and for none elsewhere.
- A fold whose held-out risk score is numerically constant (SD < 10⁻⁶, which happens when the
  fit has no survival-active program) has C undefined. In fold 5 every unsupervised arm,
  including U2064, ended with no survival-active program and a risk SD of about 10⁻⁸, yet the
  concordance function returned C ≈ 0.70–0.75 from rounding noise. Those values are excluded
  below; the as-computed table is in `cv_summary_as_computed.csv`. No S(N) fit had a constant
  risk score.

## Results

**Cross-validated C (constant-risk folds excluded):**

| Genes | U(N) mean C (SE) | S(N) mean C (SE) | Paired S − U (SE), folds |
|---|---|---|---|
| 251 | 0.557 (0.044), 4 folds | 0.628 (0.037), 5 folds | +0.038 (0.032), 4 |
| 500 | 0.618 (0.031), 4 folds | 0.630 (0.035), 5 folds | −0.016 (0.015), 4 |
| 1,000 | 0.619 (0.027), 4 folds | 0.654 (0.029), 5 folds | +0.014 (0.008), 4 |
| 2,064 (benchmark genes) | 0.579 (0.046), 4 folds | — | — |

Mean number of survival-active programs per fold fit: S(N) 1.6, 2.2, 2.8; U(N) 1.4, 1.8, 1.6;
U2064 1.4.

**Programs (full-training fits compared with the U2064 fit).** The U2064 fit has two
survival-active programs (its programs 3 and 7, the adverse and protective programs reported
earlier). Every S(N) fit recovers both: activity correlation 0.87–0.97 and loading correlation on
shared genes 0.93–0.99, with coefficients of the same sign and similar size; risk scores correlate
with U2064 at 0.96 (251 genes), 0.97 (500) and 0.99 (1,000). U(500) and U(1,000) also recover both
(r ≥ 0.98). U(251) does not: its counterpart of program 3 has a zero coefficient and its
counterpart of program 7 is weakly matched (r = −0.27), and its risk score correlates with U2064
at only 0.42.

## Interpretation

At matched size, the survival screen gives equal or higher CV concordance than the unsupervised
screen at 251 and 1,000 genes and similar concordance at 500. The differences are within about two
standard errors with five folds of about 55 patients, so the evidence is suggestive, not
conclusive. Two results are clearer than the C differences. First, a survival screen of 251 genes
keeps both benchmark survival programs and nearly the same risk score, while an unsupervised
screen of the same size loses them. Second, no survival-screened fit collapsed to zero
survival-active programs, while every unsupervised arm did so in one fold. The screen makes the
current model's survival signal more stable, which is consistent with the finding that survival
has almost no influence on the loadings in the current fitter: the screen supplies, through
feature choice, the supervision that the fitting algorithm does not.

Limitations: training CV only and conditional on scaling done before the folds; one seed; the
double-log preprocessing; the screen ranks only genes already in the 2,064 benchmark set.

## Next

External scoring of the frozen full-training fits (Task E5) waits for Andrew's review of these
results. If the rebuilt method lets survival inform the loadings, the comparison should be
repeated with it, since the benefit of the screen may then shrink.
