# Attempt S-11 findings: the survival-weighted YFB fit on the real training data

Written 2026-11-04. Registration: [`2026-11-04-supervised-real-data-plan.md`](../plans/2026-11-04-supervised-real-data-plan.md).
Outputs: `results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04/` (`summary.csv`,
`fits.csv`, `size_preferred.csv`, `size_agreement.csv`, `size_elbow.png`, `programs_D4.csv`,
`screen_paired.csv`, `null_cells.csv`). All 144 real-data fits and 60 null fits ran without
error. Training patients only: 273 TCGA_PAAD and CPTAC patients, 139 events, single log,
per-cohort z-scores, 17,900 shared genes. No external cohort was used.

## Summary

- **Survival barely shapes the real-data loadings at w_S = 1.** Survival's share of the loading
  precision is about 10⁻⁸ to 10⁻⁶ in the size fits (down to 10⁻¹¹ in the screening fits),
  against 0.01 to 0.2 in the simulation. Two factors combine: the real coefficients are small
  (|E β| ≤ 0.6 per SD of the program score) and the share scales with E[β²], and the real
  programs are high-variance programs whose genomics precision is large, while the simulated
  program carries 0.25% to 5% of the variance. The w_S = 1 and w_S = 0 loadings agree at
  |r| = 1.000 for every program. The w_S = 1 programs have the same top genes as the unsupervised (w_S = 0) fit.
- **CV C at the D4 genes (2,218): 0.663 (SE 0.024)** for the w_S = 1 fit, against 0.652 (SE
  0.026) for the two-step fit (w_S = 0 loadings, then Cox). The paired fold difference is
  +0.011 (SE 0.005). Since the loadings agree, the gain comes from the β update, not from
  supervision of the loadings.
- **Three survival-active programs at D4** (|E β| ≥ 1.96 SD): a hypoxia and glycolysis versus
  normal-pancreas axis (protective at its normal pole), a ribosomal and proliferation versus
  classical-epithelial axis (adverse at its proliferative pole), and a third adverse program.
  Gene lists are below.
- **Feature-set size: the D4 set is preferred.** By the registered rule (the smallest set within
  one SE of the all-genes C whose fit matches every survival-active all-genes program at
  |r| ≥ 0.7), the preferred set is D4 for both arms. For w_S = 1, C falls from 0.663 at 2,218
  genes to 0.643 to 0.648 at 8,000 genes and above. For the two-step fit, C is flat at 0.650 to
  0.653.
- **Survival screening adds no clear gain.** Screened sets S(N) beat unsupervised sets U(N) by
  0.003 to 0.012 with w_S = 1 and by 0.001 to 0.021 for the two-step fit, each within about one
  paired SE. No screened set beats the w_S = 1 fit on all 2,218 D4 genes.
- **Null check passes.** On the S-7 datasets with permuted survival, mean false-active counts
  at w_S = 1 are 0.1 to 0.5 per cell (criterion ≤ 1), and all 60 fits converge.

## Correction to the registration: the w_S = 0 arm

With `canonicalize = TRUE`, w_S also multiplies the Cox term in the β update, so at w_S = 0 the
fit has β = 0 and no risk score. Every w_S = 0 fold has C undefined. The registration (and the
S-9 registration) described this arm wrongly as "β still fitted". The unsupervised comparison
the registration intended is built after the fact from the saved w_S = 0 loadings: a Cox
regression of survival on the training program scores Z = Y E[F], scored on the held-out fold
(`w0cox` in the tables). This arm was added after the registration. It needs no refit.

## Step S3: training CV and programs at the D4 genes

| Fit | Mean fold C (SE) | Pooled out-of-fold C | CV partial log-lik. | Folds converged | Survival-active (all patients) |
|---|---|---|---|---|---|
| w_S = 1 | 0.663 (0.024) | 0.663 | −775 | 3/5 | 3 |
| w_S = 0 then Cox | 0.652 (0.026) | 0.647 | −777 | 3/5 | 4 |

Two w_S = 1 folds stop at the 300-sweep cap (283, 261, 300, 300 and 214 sweeps). The relative
ELBO change was still above 10⁻⁵, and the max |E β| is 0.46. Neither the τ mechanism nor the
β mechanism of the simulation applies: reconstruction is bounded and the coefficients are
small. The fit is moving slowly near its optimum.

Programs of the all-patient w_S = 1 fit (from `programs_D4.csv`; sign convention: positive β
means a higher score raises the hazard; top genes at each pole of the loading):

| Program | E β (SD) | Active | Positive pole | Negative pole |
|---|---|---|---|---|
| 2 | −0.41 (0.08) | yes | ATP2A3, ACACB, METTL7A, ADH1B, ACSS1, C7 | LOXL2, MMP14, P4HA1, LDHA, SLC16A3, PGK1, PLAU |
| 3 | +0.22 (0.07) | yes | RPS7, RPL6, NCAPD2, EEF1D, MCM7 | KIAA1211, SYTL2, CTSE, TSPAN8, TOX3 |
| 5 | +0.23 (0.08) | yes | UHMK1, ARHGAP5, TAOK1, RAD21, DDX21, PRKDC | CHPF, RCN3, MXRA8, LOXL1, SPON2 |
| 4 | +0.15 (0.08) | no (z = 1.9) | CTRB2, CTRB1, CPB1, REG1B, CELA3A (acinar) | BANK1, MS4A1, CD79A, CR2 (B cell) |
| 7 | −0.05 (0.06) | no | ANO1, KRT7, PPL | PRDX4, PPIB, XBP1, TXNDC5 |
| 1 | 0.00 | no | TCF4, FBN1, BNC2 (stromal) | RAB25, CLDN7, PRSS8 (epithelial) |
| 6 | 0.00 | no | CD74, HLA-DRA, CD68, LAPTM5 (myeloid) | IGF1R, PBX1, LTBP4 |

The two-step fit has the same programs and top genes. It also calls program 4 active
(z = 2.3). The w_S = 1 fit's Normal β prior (one prior variance per program, estimated as max(0, x² − s²))
shrinks programs 1 and 6 to exactly zero, while
the unpenalized Cox fit gives them small nonzero coefficients. These biological labels come from
the top genes only. A gene-set analysis has not been run on these programs.

## Feature-set-size curve

Nested outcome-free sets by the per-cohort mean/variance rank. Mean fold C ± SE:

| Genes | w_S = 1 | Two-step | Paired w_S = 1 − two-step (SE) | w_S = 1 folds converged | w_S = 1 match with all-genes programs (|r|) |
|---|---|---|---|---|---|
| 2,218 (D4) | 0.663 ± 0.024 | 0.652 ± 0.026 | +0.011 (0.005) | 3 | 0.78, 0.89 |
| 3,000 | 0.663 ± 0.026 | 0.653 ± 0.027 | +0.010 (0.003) | 3 | 0.72, 0.80 |
| 5,000 | 0.657 ± 0.030 | 0.652 ± 0.031 | +0.005 (0.002) | 4 | 0.69, 0.77 |
| 8,001 | 0.647 ± 0.030 | 0.650 ± 0.033 | −0.003 (0.005) | 5 | 0.74, 0.97 |
| 12,001 | 0.643 ± 0.037 | 0.653 ± 0.037 | −0.010 (0.003) | 5 | 0.96, 0.80 |
| 17,900 (all) | 0.648 ± 0.038 | 0.652 ± 0.037 | −0.004 (0.004) | 5 | 1.00, 1.00 |

`size_elbow.png` plots the first two columns. By the CV partial log-likelihood, which the
registration also lists, the ordering differs: the w_S = 1 value is −775.2 at D4 and −773.2 at
5,000 genes, and the two-step value improves steadily with size (−777.1 to −772.0). The
preferred-size rule uses C and still picks D4, but larger sets are not worse by every metric. The curve has no elbow above the D4 size: the
smallest set does as well as any larger one, so the more parsimonious D4 set is preferred,
which the registered rule confirms. The 5,000-gene fit matches one all-genes program at
r = 0.69, just under the 0.7 threshold. The risk scores of every size correlate with the
all-genes risk at 0.94 or more. The w_S = 1 advantage over the two-step fit shrinks with size
and reverses above 5,000 genes. Fold SEs are 0.024 to 0.038, so all sizes are within one SE of
one another.

## Survival screening before the joint fit

U(N): the N genes with the best outcome-free ranks inside D4. S(N): the D4 genes ranked by a
cohort-stratified univariate Cox test on fold-training patients, top N.

| N | Fit | C, U(N) | C, S(N) | Paired S − U (SE) | S(N) genes also in U(N) |
|---|---|---|---|---|---|
| 250 | w_S = 1 | 0.638 | 0.641 | +0.003 (0.016) | 11% |
| 501 | w_S = 1 | 0.637 | 0.646 | +0.009 (0.009) | 22% |
| 1,000 | w_S = 1 | 0.643 | 0.655 | +0.012 (0.011) | 44% |
| 250 | two-step | 0.626 | 0.647 | +0.021 (0.018) | 11% |
| 501 | two-step | 0.645 | 0.650 | +0.004 (0.017) | 22% |
| 1,000 | two-step | 0.647 | 0.648 | +0.001 (0.006) | 44% |

Screening on survival selects mostly different genes from the outcome-free screen, but no
screened set improves on the w_S = 1 fit with all 2,218 D4 genes (0.663). Gains over U(N) are
within about one paired SE. Smaller sets also converge less reliably: two U(N) w_S = 1 cells
(N = 250 and 1,000) converge in 0 of 5 folds. The October 26 results with the historical double-log, α_F = 0 model
(S − U of 0.038, −0.016 and 0.014) stay on record. They used different folds and preprocessing
and are not directly comparable. Together they point the same way: a survival screen does not
give a consistent gain, and it weakens the claim that the factorization finds the prognostic
programs itself.

## Null calibration check (simulation)

| p | Planted variance | Mean false-active | Converged |
|---|---|---|---|
| 50 | 0.25% / 1% / 5% | 0.1 / 0.1 / 0.2 | 10/10 each |
| 2,000 | 0.25% / 1% / 5% | 0.4 / 0.5 / 0.2 | 10/10 each |

The criterion (at most 1 per cell) holds in every cell. With no survival signal, w_S = 1 does
not create prognostic programs.

## Interpretation and next steps

On the real data at w_S = 1 the joint fit's loadings are effectively those of the unsupervised
factorization, because the real survival effects are small and the survival share of the
loading precision scales with E[β²]. The modest CV gain over the two-step fit comes from fitting
β jointly with a shrinkage prior. The simulation showed that a larger weight breaks convergence
at 2,000 genes, so raising w_S is not the fix. A calibrated tempering of the Cox likelihood
(SafeBayes or Holmes–Walker) is the principled way to set how much survival should move the
loadings, and it is the main open method question.

Choices fixed by training CV: w_S = 1, the D4 gene set (2,218 genes), no survival screen. Next,
under the 2026-11-04 decision: one scoring of this fixed model on the five external cohorts,
registered first and run after Andrew has seen these results. Then the M-modality extension.
