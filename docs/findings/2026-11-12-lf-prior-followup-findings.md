# Follow-up A2 findings: point-exponential L and F with the cohort mean difference removed

Date: 2026-11-12. Author: Andrew Walther.

**Result.** Neither setup makes point-exponential competitive with point-Laplace. Both
recover most of the concordance that task A's point-exponential arm lost. SZ removes the
cohort artifact found in task A. PES leaves it in the programs but absorbs it in the survival
part of the fit. Both still fail the same two registered criteria.

- **Shifted z-score (SZ).** Mean within-cohort CV concordance is 0.635 against 0.658 for
  point-Laplace, a paired difference of −0.023 (SE 0.008), which is 2.9 SEs below. Only 1
  of 5 folds converged in 300 sweeps. No program correlates with the cohort label (largest
  |r| 0.05), and 2 to 4 programs per fold are survival-active. In simulation SZ behaves like
  task A's point-exponential arm. It averages 2.8 to 4.4 false-active programs per fit
  (limit 1), and every one of its 40 fits has at least one.
- **X/sd with cohort strata (PES).** Mean within-cohort CV concordance is 0.628, a
  difference of −0.030 (SE 0.015) from point-Laplace, and no fold converged. The program
  directions are the same as in task A's point-exponential fit (|r| ≥ 0.9999 for every
  program in every fold), so the programs still separate TCGA from CPTAC. What changes is β.
  3 or 4 programs become survival-active, where every coefficient in task A's fit was below
  10⁻⁷.

Against task A's point-exponential arm, within-cohort concordance is higher by 0.062
(SE 0.031, 2.0 SEs) for SZ and by 0.055 (SE 0.035, 1.6 SEs) for PES.

Registration: the 2026-11-12 "A2" entry in the
[attempts log](../plans/2026-10-16-round2-attempts-log.md), a follow-up to
[task A](2026-11-12-lf-prior-comparison-findings.md). Runner:
`results/benchmark_sim/run_lf_prior_followup.R`. Aggregate outputs:
`results/benchmark_sim/outputs/lf_prior_2026-11-09/a2/`.

## Design

Both setups use `fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1)` with point-exponential L and F,
per-program Normal β, K = 7, no pruning, at most 300 sweeps and tolerance 10⁻⁵, the same
call as task A's PE arm apart from the input and the strata.

| Setup | Input | Strata | Cox projection |
|---|---|---|---|
| SZ | per-cohort z-score plus c_j = −min_i z_ij over all 273 patients (nonnegative, minimum 0 per gene) | none | centered with one pooled mean |
| PES | X/sd within cohort (task A's PE input) | `strata_id` = cohort | centered within cohort |
| PL (task A, reused) | per-cohort z-score, point-Laplace L and F | none | as is |
| PE (task A, reused) | X/sd within cohort | none | centered with one pooled mean |

The SZ constant is the same in both cohorts, so each cohort keeps gene mean c_j and SD 1. It
depends on expression only, not on outcomes. It is computed once from all 273 patients,
held-out folds included, as registered.

**Real data.** The 273 TCGA and CPTAC training patients, the 2,218 D4 genes and the five S-11
folds, as in task A, plus a full-data fit for each setup. Held-out risk is
Y[test] E[F] E[β] in the setup's own input. A per-gene constant moves every held-out score by
the same amount, and a per-cohort offset does not change within-cohort pairs, so the primary
metric, within-cohort concordance, is unaffected by either. PES scores keep each cohort's
offset, so its pooled concordance (mean 0.564) is not interpreted.

**Simulation.** SZ only, in task A's four cells (10 seeds each). The synthetic training
z-score is shifted by its own per-gene minimum, and the validation z-score by the training
constants. Scoring is `s9_score()`, as in task A. PES was not refitted in simulation. The
synthetic data have one cohort, and with one stratum the PES fit is task A's PE fit, so A's
PE results stand for PES in criteria 1 to 3. Those criteria therefore say nothing about the
strata, which only matter with two cohorts.

All 52 planned fits returned (40 simulation, 12 real), with no errors and no fit warnings.
With task B and other jobs sharing the machine, real fits took 53 to 214 seconds and
simulation fits 4 to 11 seconds on one core each.

## Simulation (SZ)

| Cell | Arm | Recovered | Mean false-active | Converged | Median best \|r\| | Median validation C | Median calibration slope | Median programs with nonzero loadings |
|---|---|---|---|---|---|---|---|---|
| 0.05% (S-14) | SZ | 6/10 | 2.8 | 10/10 | 0.741 | 0.841 | 1.30 | 7 |
| 0.05% (S-14) | PE | 8/10 | 2.6 | 10/10 | 0.762 | 0.839 | 1.26 | 7 |
| 0.05% (S-14) | PL | 5/10 | 3.0 | 10/10 | 0.711 | 0.832 | 1.58 | 7 |
| 0.25% | SZ | 10/10 | 3.1 | 10/10 | 0.923 | 0.941 | 1.37 | 6 |
| 0.25% | PE | 10/10 | 3.0 | 10/10 | 0.923 | 0.940 | 1.39 | 6 |
| 0.25% | PL | 10/10 | 0.9 | 10/10 | 0.934 | 0.947 | 1.45 | 7 |
| 1% | SZ | 10/10 | 4.1 | 10/10 | 0.973 | 0.970 | 1.70 | 6 |
| 1% | PE | 10/10 | 4.0 | 10/10 | 0.973 | 0.971 | 1.82 | 6 |
| 1% | PL | 10/10 | 0.1 | 10/10 | 0.973 | 0.974 | 2.16 | 7 |
| 5% | SZ | 10/10 | 4.4 | 10/10 | 0.992 | 0.981 | 2.13 | 5.5 |
| 5% | PE | 10/10 | 4.3 | 10/10 | 0.992 | 0.982 | 2.07 | 5.5 |
| 5% | PL | 10/10 | 0.0 | 10/10 | 0.991 | 0.987 | 3.23 | 6 |

In the one-cohort simulation, SZ and task A's PE differ only by a per-gene constant
added to the input, and their results are nearly the same in every cell. Removing the cohort
mean cannot change the false-activity failure here, because the simulation has no cohort
difference to remove. The failure is seen with point-exponential L and F and a centered
projection. No arm separates the prior from the centering.

**Size of the false-active coefficients (described after seeing the results, not a
criterion).** For each fit with a matched program (best |r| ≥ 0.7), the table gives the
largest ratio |E β_j| / |E β_k| over false-active programs j, with k the matched program.

| Cell | Fits with a ratio | Median of the largest ratio | 90th percentile | Largest | Mean of the per-fit mean ratios | Median \|E β\| of the matched program |
|---|---|---|---|---|---|---|
| 0.05% (S-14) | 6 | 0.28 | 0.36 | 0.41 | 0.18 | 2.5 |
| 0.25% | 10 | 0.11 | 0.17 | 0.18 | 0.08 | 5.4 |
| 1% | 10 | 0.20 | 0.30 | 0.32 | 0.09 | 7.0 |
| 5% | 10 | 0.23 | 0.42 | 0.51 | 0.12 | 8.6 |

In the S-14 cell, the 4 fits with no matched program have no reference coefficient and are
left out. The false-active coefficients are mostly small, about a tenth to a quarter of the
matched program's. In a few fits the largest reaches 0.4 to 0.5.
Whether that is acceptable is Andrew's judgment. It does not change criterion 2, which is
applied literally.

## Real training data

| Setup | Mean within-cohort C | Difference from PL (SE) | Difference from PE (SE) | Mean pooled C | Folds converged | Survival-active programs per fold | Max \|cor\| of a program with the cohort label | Max \|cor\| with mean expression | Max within-cohort \|cor\| with mean expression |
|---|---|---|---|---|---|---|---|---|---|
| PL (task A) | 0.658 | | | 0.663 | 3/5 | 2 to 4 | 0.03 to 0.05 | 0.33 to 0.42 | 0.65 to 0.70 |
| PE (task A) | 0.573 | −0.085 (0.031) | | 0.502 | 0/5 | 0 | 0.998 | 0.88 to 0.91 | 0.79 to 0.90 |
| SZ | 0.635 | −0.023 (0.008) | +0.062 (0.031) | 0.638 | 1/5 | 2 to 4 | 0.03 to 0.05 | 0.43 to 0.50 | 0.81 to 0.88 |
| PES | 0.628 | −0.030 (0.015) | +0.055 (0.035) | not interpreted | 0/5 | 3 to 4 | 0.998 (offset, see below) | 0.88 to 0.91 | 0.79 to 0.90 |

The pooled-C difference for SZ is −0.025 (SE 0.004). Fold-wise values are in
`a2_real_fits.csv` and per-program values in `a2_real_programs.csv`. All seven programs keep
nonzero loadings in every setup and fold. The full-data fits have 3 survival-active programs
each. The fold-wise within-cohort C for SZ is 0.656, 0.678, 0.611, 0.670 and 0.560, against
0.666, 0.732, 0.622, 0.689 and 0.581 for PL, so SZ is lower in every fold. Fold 5 has the
lowest TCGA concordance for every arm (PL 0.456, SZ 0.440, PES 0.430). The held-out
within-cohort partial log-likelihood, summed over the five folds, is lower than
point-Laplace's by 2.5 for SZ, 6.7 for PES and 13.9 for task A's PE, the same order as the
concordance.

**The last column was added after the smoke test.** Under PES each program's uncentered
training scores still carry the cohort offset that the strata absorb, so task A's
correlation with the TCGA indicator stays at 0.998 and says nothing about what is left
within cohort. The added column centers both the program score and the patient's mean log
expression within cohort before correlating. It is computed from the saved E[F] for all four
arms (`a2_real_within_cohort_cor.csv`). Within cohort, every arm has a program that tracks
mean expression. That program reaches 0.65 to 0.70 under point-Laplace and 0.79 to 0.90
under the three point-exponential arms.

**SZ programs.** The shifted z-score removes the cohort artifact by construction. Every gene
has mean c_j in both cohorts, so on the full data any program's score has zero correlation
with the TCGA indicator, whatever F is. The fold values of 0.03 to 0.05 come from leaving out
the held-out patients, and point-Laplace on z-scores has the same property.
Point-exponential programs still follow mean expression more closely than point-Laplace
programs do (pooled 0.43 to 0.50 against 0.33 to 0.42, within cohort 0.81 to 0.88 against
0.65 to 0.70). Every program's loadings are nonnegative, so a program's score sums positive
contributions, and a program that loads broadly rises with overall expression.

**PES programs.** The strata change the survival part of the fit, while the reconstruction
term is the same as in task A's PE fit. In every fold, each PES program matches one PE program
in direction, with |r| ≥ 0.9999. Their scale differs. In 6 of the 7 programs the column norm
of E[F] is about 0.07 to 0.17 under PES against 0.007 to 0.017 under PE, so the survival
term does affect F, but not which genes a program loads on. Task A's PE coefficients were all
below 10⁻⁷ in absolute value, which a tenfold scale change cannot explain. With the cohort
offset absorbed by the strata, the within-cohort variation of the same programs reaches
|E β| of 0.2 to 0.7, and 3 or 4 programs are survival-active. The programs still describe the
TCGA versus CPTAC difference first.

**Overlap with the current programs.** Top-100 gene overlap with the current S-11 fit's
programs 2 and 3: the best SZ program shares 16 genes with program 2 and 50 with program 3.
PES shares 8 and 42, the same as task A's PE fit. The best SZ programs share 79 and 76 genes
with current programs 4 and 6.

## Criteria

| Setup | Criterion | Result | Pass |
|---|---|---|---|
| SZ | 1. Recovery in every cell ≥ PL's minus 1 of 10 | 6 vs 5, 10 vs 10, 10 vs 10, 10 vs 10 | yes |
| SZ | 2. Mean false-active programs ≤ 1 in every cell | 2.8, 3.1, 4.1, 4.4 | **no**, every cell |
| SZ | 3. ≥ 9 of 10 converged in every cell | 10/10 in every cell | yes |
| SZ | 4. Within-cohort CV C not lower than PL's by more than one SE, all folds converge | −0.023 (SE 0.008), 1/5 converged | **no** |
| PES | 1 to 3 | task A's PE results (one stratum): 8 vs 5 and 10/10 recovery, 2.6 to 4.3 false-active, 10/10 converged | 1 yes, 2 **no**, 3 yes |
| PES | 4. As above | −0.030 (SE 0.015), 0/5 converged | **no** |

Point-Laplace stays as the L and F prior. Both setups fail criterion 2 in every simulation
cell and fail criterion 4 on both concordance and convergence. If criterion 2 were relaxed
for small false-active coefficients, both setups would still fail criterion 4. Under the
literal rule, SZ concordance is 2.9 SEs below point-Laplace, and the folds do not converge.

## Limits

- The SZ shift constants use all 273 patients, held-out folds included, as registered. They
  depend on expression only, and a per-gene constant does not change held-out concordance.
- The 300-sweep cap is reached in 4 of 5 SZ folds and in every PES fold. Concordance is from
  the state at the cap. The full-data fits also reach the cap.
- PES has no simulation fits of its own. Its criteria 1 to 3 come from task A's PE fits,
  which equal PES with one stratum.
- The within-cohort mean-expression correlation was added after the smoke test, and it is
  descriptive only.
- The paired SE is the SD of the five fold differences divided by √5, as registered for A.
  The folds share training patients, so this SE probably understates the uncertainty.
- Point-Laplace itself converges in 3 of 5 folds, so it would not meet the convergence half of
  criterion 4 either.
- No external cohort was used.

## Text for the shared status files

Proposed for the session that owns DECISIONS.md, PROJECT_STATUS.qmd, the week plan and the
attempts log.

**Attempts log, A2 outcome (append under the A2 entry):**
> Outcome (2026-11-12): neither setup qualifies. The shifted z-score removes the cohort
> artifact (largest program correlation with cohort 0.05), but within-cohort CV C is 0.635
> against 0.658 for point-Laplace (−0.023, SE 0.008), 1/5 folds converge, and the simulation
> gives 2.8 to 4.4 false-active programs per fit. Cohort strata on X/sd keep the directions of task A's
> point-exponential programs (|r| ≥ 0.9999) and make 3 or 4 of them survival-active.
> Within-cohort C is 0.628 (−0.030, SE 0.015), and no fold converges. False-active
> coefficients in the SZ simulation are mostly 0.1 to 0.3 of the matched program's |β|
> (largest 0.51). Findings: `docs/findings/2026-11-12-lf-prior-followup-findings.md`.

**DECISIONS.md (addendum to the task A entry):**
> 2026-11-12 follow-up A2: removing the TCGA/CPTAC mean difference, by a shifted per-cohort
> z-score or by cohort strata on X/sd, does not make point-exponential L and F competitive.
> Both fail the false-activity limit and the within-cohort CV criterion (−0.023 and −0.030 from
> point-Laplace, 1/5 and 0/5 folds converged). Point-Laplace on per-cohort z-scores stays.

**PROJECT_STATUS.qmd / week plan (one line):**
> A2 complete: shifted z-score and cohort strata both fail. Point-Laplace L/F stays
> ([findings](docs/findings/2026-11-12-lf-prior-followup-findings.md)).
