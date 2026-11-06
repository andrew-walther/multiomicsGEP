# Real-data prognostic comparators

The [registration](../plans/2026-11-06-real-data-comparators-plan.md) fixes
the S-11 supervised and two-step references and trains the comparators on
TCGA and CPTAC. The reference factorizations were reused. The two-step Cox
coefficients were reconstructed from the fixed unsupervised program scores.

The clinical comparison requires covariates that can be harmonized across
the training data and every external cohort. None of the requested fields
met that requirement. Age was recorded in both training cohorts and three of
the five external cohorts, and margin in four of seven cohorts, so a model on
a subset of cohorts was possible but outside the registered rule. The field audit preserves the original source column
names. For example, Puleo's FFPE block age is not patient age, and Moffitt's
T, N and M columns are component categories rather than an overall stage.

## Reproduce and trace the results

From the repository root with the local S-11 inputs and the pinned DeSurv
library installed:

```sh
export PDAC_DATA_ROOT=/absolute/path/to/PDAC_data
Rscript results/benchmark_sim/run_prompt6_desurv.R 8
Rscript results/benchmark_sim/run_prompt6_comparators.R
```

The first command uses the authors' reported training settings (rank 3 and
the published α, λ and ν, which the authors selected on these two training
cohorts; the training-only Bayesian optimization in `run_desurv_h2h.R` was not
rerun), their
preprocessing, independent seed fits, and a consensus initialization for
the final fit. Its checkpoints stay under ignored `data/prompt6_desurv/`.
The second command fits the PurIST probability Cox model and cross-validates
elastic-net Cox on the fixed D4 features and S-11 folds. It scores the
external samples, reproduces the S-12 reference concordances and writes the
following aggregate tables to
[`outputs/prompt6_comparators/`](../../results/benchmark_sim/outputs/prompt6_comparators/).

| File | Reported values |
|---|---|
| `external.csv` | Eligible, scored, dropped and event counts, per-cohort C and bootstrap limits, calibration slope and SE |
| `paired.csv` | Supervised-minus-comparator C on common rows, eligible strata, row and donor counts, donor-bootstrap limits and replicate count |
| `pooled_hr.csv` | Cohort-stratified HR per pooled external SD, interval, p-value and denominator |
| `clinical_fields.csv` | Source clinical fields by cohort |
| `training.csv` | Training denominators and fit counts by arm, excluding DeSurv |
| `tuning.csv` | Elastic-net alpha, lambda and CV deviance by candidate |
| `desurv_starts.csv` | One row per author seed, including convergence, errors, warnings and runtime |
| `desurv_final.csv` | Final DeSurv training denominator, convergence and iterations |

External outcomes entered only the performance and calibration analyses.
The pooled C contrast gives each represented cohort equal weight. The paired
bootstrap keeps the two PACA-AU assays together at the donor level. The HR
scaling uses the pooled external score SD as a descriptive unit.

## Results

The paired external analysis used 616 records from 568 donors in each
contrast. The 48 repeated records are the 47 donors measured on both PACA-AU
assays and one donor with two RNA-seq samples. Every interval includes zero.

| Comparator | Supervised minus comparator C | 95% donor bootstrap interval |
|---|---:|---:|
| Two-step | 0.0071 | −0.0029 to 0.0173 |
| PurIST | 0.0211 | −0.0134 to 0.0573 |
| Elastic net | 0.0015 | −0.0146 to 0.0180 |
| DeSurv | 0.0069 | −0.0166 to 0.0300 |

The contrast table and denominators come from `paired.csv`. DeSurv completed
all 100 starts with no errors or convergence failures, and its final fit
converged in 34 iterations on 273 training patients
(`desurv_starts.csv`, `desurv_final.csv`). Elastic net selected the ridge
endpoint (α = 0, λ = 39.4461). Some lower-penalty solutions at nonzero α hit
the iteration limit. The candidate warnings and selected settings are in
`tuning.csv`. All eligible patients were scored by every fitted arm
(`training.csv`, `external.csv`). The clinical arm was excluded before fitting.

The elastic-net tuning used the fixed D4 feature set and already prepared
training matrix, so its CV objective is conditional on that preprocessing.
The external cohorts have a history of project development use, described
in the [external-cohort audit](2026-11-06-external-cohort-uses.md). The new
comparisons are descriptive and the intervals do not establish equivalence.

## Review notes (2026-11-06)

An independent rerun of `run_prompt6_comparators.R` reproduced every committed
CSV exactly. The ridge penalty selected by cross-validation lies inside the
λ path (position 36 of 100), so the elastic-net choice is not a boundary
solution. The glmnet iteration warnings at α ≥ 0.75 (`tuning.csv`) affect only the smallest penalties of those paths, but they may have favoured the ridge choice. The field table originally listed TCGA margin as unavailable,
although TCGA records `Residual.tumor`. The table and runner were corrected.
The exclusion of the clinical arm is unchanged because neither margin nor any
other field is available in all seven cohorts.

The pooled HRs in the registered analysis treat all 616 records as
independent, although 47 donors contribute a record from each PACA-AU assay.
The runner now also reports a robust interval clustered on donor
(`donors`, `donor_lo`, `donor_hi` in `pooled_hr.csv`). The point estimates
and the registered model-based intervals are unchanged, and the clustered
intervals are slightly wider. For the supervised score the HR is 1.37 with a
model-based interval of 1.24 to 1.51 and a donor-clustered interval of 1.22
to 1.53.

The DeSurv arm deviates from plan item 3, which named the training-only
tuning workflow. It used the settings the authors selected on the same 273
training patients with their own folds, not a new search on the S-11 folds.
The per-cohort C intervals resample records. This matters only for PACA-AU
RNA-seq, where 52 records come from 51 donors.
