# Methylation-only programs on matched TCGA patients

Date: 2026-11-12. Author: Andrew Walther.

**K = 5 result.** The full sweep and timing-eligible CV are complete. Mean
methylation concordance is below matched expression in all seven completed
comparisons. No PL program meets the methylation-only rule. Formal PE
classification is unassessable because the expression PE full fit capped.
These are descriptive results from 144 patients, not a model adoption decision.

**K = 7 status: stopped at Andrew's request on November 12 to prioritize K = 5.**
No K = 7 process remains active. Nine methylation full fits completed, six
converged and three capped. The final all-CpG PE fit was interrupted after
46.0 minutes without a final checkpoint. Expression fits, official CV and ICGC
scoring were not run. The approved 600-iteration check is deferred.

**K = 5 follow-up:** all 57 scheduled full/CV fits are saved. Of these,
45 converged and 12 capped. There were no fitting errors or fit warnings.

**Scientific conclusion for the stopped K = 7 run:** this partial run establishes computing cost and
full-fit diagnostics only. It does not answer relative predictive performance
or identify validated methylation-only prognostic programs.

## Question and design

This is the fixed **K = 7** diagnostic study. The working rank has since changed
to **K = 5**. Its matched methylation and expression follow-up is reported below under
“Rerun at K = 5.”

This study asks whether methylation alone identifies survival-active programs
that the expression fit does not, and how concordance and computing cost change
with the number of CpGs. The controlling specification is
[Prompt 1](../prompts/2026-11-12-codex-c-d-e-prompts.md), with the November 12
methylation decision in `DECISIONS.md`. Neither Prompt 2 nor tasks D and E is part
of this analysis.

The training set contains 144 matched TCGA patients with 75 events. Expression
uses the 2,218 fixed D4 genes on those same patients. Both modalities use the same
saved five-fold split, balanced by event with seed 2026.

The prepared methylation universe contains 305,352 probes. The saved mask based
on observed TCGA β means before imputation retains 182,934 probes with means in
[0.1, 0.9], reproducing the registered count. Means computed after imputation
would retain 182,946. Sixteen probes cross a boundary between these conventions.
Eligible probes are ranked by variance of their prepared, imputed β values,
following the existing decreasing-variance ranking convention. The nested sizes
are 10,000, 25,000, 50,000, 100,000 and all 182,934 probes. Screening is
outcome-free and computed once on all TCGA patients before splitting folds, as
registered. The data remain β values throughout.

Each size has two fixed setups. PL uses CpG z-scores and point-Laplace priors on
L and F. PE divides each CpG by its training SD, without subtracting its mean,
and uses point-exponential priors with a centered Cox projection. Scaling
parameters are estimated separately within each training fold and applied
unchanged to its held-out patients. The PE projection center also comes from
that training fold. Full-fit scaling is frozen on TCGA for external scoring.
The fit uses seven columns without pruning, canonicalization, unit survival
weight and Normal β. The iteration cap, tolerance and alpha come from
`config/globals.yml`. No failed fit is tuned.

## Full-fit results

<!-- FULL-FIT-RESULTS -->
The original 300-iteration results are preserved below. Counts from capped
iterates are diagnostics. The stopped PE row has no fitted-state result.

| CpGs | Setup | Status | Iterations | Process minutes | Peak RSS MiB | Nonzero programs | Survival-active |
|---:|:---|:---|---:|---:|---:|---:|---:|
| 10,000 | PL | converged | 123 | 1.18 | 1422 | 7 | 0 |
| 10,000 | PE | converged | 284 | 2.12 | 1156 | 7 | 0 |
| 25,000 | PL | converged | 113 | 4.01 | 1464 | 7 | 0 |
| 25,000 | PE | converged | 288 | 16.14 | 1394 | 7 | 0 |
| 50,000 | PL | capped | 300 | 12.87 | 2272 | 7 | 0 |
| 50,000 | PE | capped | 300 | 16.36 | 2148 | 7 | 0 |
| 100,000 | PL | converged | 150 | 14.98 | 3303 | 7 | 1 |
| 100,000 | PE | capped | 300 | 22.20 | 3248 | 7 | 2 |
| 182,934 | PL | converged | 113 | 18.84 | 3852 | 7 | 0 |
| 182,934 | PE | stopped_by_user | NA | 46.02 | 4282 | NA | NA |

The converged 100k PL fit has one survival-active program. Without expression
fits or fold stability results it cannot be classified as methylation-only.
The two active coefficients in the capped 100k PE iterate are not evidence of
a converged prognostic solution.

The interrupted PE timing and memory describe its partial execution, not a
completed fit. Its CV projection is a lower bound based on the elapsed time
before interruption. That lower bound already exceeds three hours. The 10k
and 25k counts would have qualified for CV under the original timing and
convergence rules. Both 50k fits capped. The 100k combined projection was
3.10 hours, so its CV was already deferred before the stop request.
<!-- /FULL-FIT-RESULTS -->

## Computing and reproducibility

The runner is
[`run_methylation_single_modality.R`](../../results/multimodal_real/run_methylation_single_modality.R).
The following commands describe the implementation. **Do not resume K = 7
without a new instruction.** `run` would restart the interrupted PE fit from
its original initialization, then launch the remaining fits. The methods and
summary sections below describe implemented behavior, not completed stages.
From the repository root:

```sh
Rscript results/multimodal_real/run_methylation_single_modality.R prepare
Rscript results/multimodal_real/run_methylation_single_modality.R smoke
Rscript results/multimodal_real/run_methylation_single_modality.R run
Rscript results/multimodal_real/run_methylation_single_modality.R summary
Rscript tests/test_methylation_single_modality.R
```

`prepare` refuses to replace an existing design. `smoke` runs one held-out fold
on 2,000 CpGs under both setups. `run` checkpoints full fits in increasing size,
projects the combined five-fold cost from measured process elapsed time, and
runs CV only for converged setups at sizes whose projected cost is at most three
hours. Expression receives both full fits and all five CV folds. `summary`
rebuilds descriptions and fixed external scoring from saved fits.

Every fit runs in a fresh R process with one numerical thread. At most three CV
fits run concurrently. Fits at 100,000 and 182,934 CpGs run serially. Peak
resident memory is measured by macOS `/usr/bin/time -l`, so subprocess launches
need access to its system resource query. Fit time includes scaling, fitting and
projection. Process elapsed time also includes model loading and compilation.
The CV projection is five times the sum of the two full-fit process times,
without assuming parallel speedup.

Before fitting, the runner freezes committed model code and configuration in a
private source snapshot, excluding other sessions' uncommitted model edits.
The immutable runner copy is used for child fits. The source revision and
source hashes are saved with the aggregate results. Early smoke and timing
attempts before this freeze are retained privately and excluded from the final
tables. One early completed full fit was repeated after a live runner edit
caused a parsing error at process exit.

Patient matrices, design manifests, posterior fits, patient scores and detailed
logs are under ignored `data/methylation_single_modality_2026-11-12/`. Only
aggregate CSVs are published in
[`methylation_single_modality_2026-11-12`](../../results/multimodal_real/outputs/methylation_single_modality_2026-11-12/).

## Program and prediction summaries

A program has nonzero loadings when the sum of squared posterior mean feature
loadings is positive. It is survival-active when |E[β_k]| is at least 1.96
posterior SD, with positive posterior SD and nonzero loadings. Patient program
scores are the canonical YF projections used by the fitted Cox predictor.

For each converged methylation full fit, the expression comparison uses the
largest absolute correlation with any expression program under the same setup.
Fold labels are aligned to full-fit programs by the one-to-one assignment that
maximizes total absolute score correlation on common training patients. This
handles arbitrary column order and sign. Fold activity is counted only for
converged fits and defined, nonzero matched correlations. The registered study
does not specify a minimum matching correlation. Weak matches therefore remain
a limitation, and all matched correlations are reported in
`fold_program_matching.csv`.

A methylation-only program must be active in the full fit and at least three
folds, with maximum absolute expression correlation below 0.5. This description
is assessed only when all five folds were run, at least four converged and the
expression full fit converged. Deferred CV leaves the classification undefined.

Concordance uses the fixed Cox direction, with higher scores indicating shorter
survival, through `frozen_reverse_cindex()`. Held-out outcomes never determine a
sign flip. Fold-wise paired differences are methylation minus expression under
the same setup. Their SE is the sample SD of the five differences divided by
√5. A missing concordance is retained as undefined, rather than dropped to form
a smaller complete-case average. These comparisons are descriptive at n = 144.

External scoring applies each valid saved fixed full fit to the matched ICGC
donors after the existing donor-level clinical reconciliation. Primary PDAC and
the broader all-matched sensitivity population are labeled separately. Both use
frozen TCGA features, scaling, coefficients and score direction. External
concordance does not choose a setup or CpG count. A capped fit is scored as the
saved unconverged iterate and retains its convergence flag and iteration count.
Its concordance does not make it an accepted model. Program descriptions remain
restricted to converged full fits.

## Implementation walkthrough

The preparation and fitting stages are separate so that every arm receives the
same patient order, feature ranking and fold vector. `c_prepare()` reads the
existing private universe, verifies that its original screened values exactly
match the frozen baseline, then writes a private design containing numeric
patient-by-feature matrices and named survival vectors. It does not read ICGC
outcomes. `c_folds()` returns one integer label per patient. Its event groups
receive nearly equal fold sizes, without searching for a favorable split.

`c_scale()` takes a numeric matrix and a setup label. It returns the transformed
matrix plus numeric center and SD vectors, each of length equal to the number
of features. `c_transform()` applies those vectors to another matrix in the
same feature order. PL subtracts the training mean and divides by the sample
SD. PE subtracts zero and divides by the sample SD. A constant training feature
causes an explicit failed fit rather than an undocumented feature deletion.
Base R is sufficient for these operations.

`c_freeze()` saves the committed `code/` and `config/` trees and a runner copy in
the private analysis directory. Existing snapshots must match their saved hash
manifest. `c_load()` loads model functions and compiles the existing Rcpp sweep
while that frozen tree is the working directory. It then restores the project
working directory. This preserves the model implementation while allowing
other sessions to work in the same checkout.

`c_one()` takes modality, feature count, setup and held-out-fold identifiers.
It fits only the training rows, then predicts all rows with the fitted scaling
and projection center. The saved fit, score matrix and patient indices support
later program matching without refitting. Warnings are echoed and retained in
the private checkpoint. Model errors become explicit failed rows.
`c_launch()` runs this function in a separate R process, measures process wall
time and peak resident memory, and resumes a completed checkpoint without
refitting. A child that terminates before writing a checkpoint is recorded by
its exit status and private log path. Missing memory measurements stay missing.

`c_run()` schedules full fits in increasing size, writes each completed count's
aggregate timing table, and constructs eligible CV jobs. Three orchestration
workers can launch three small fitting children. Large jobs run one at a time.
The orchestrators load no model or numerical libraries before forking. The
existing compiled fitter is therefore loaded only in fresh child processes,
which avoids the project's known macOS fork-after-compilation failure.

`c_activity()` maps a fitted loading matrix and coefficient moments to logical
program indicators. `c_match()` takes two patient-by-program score matrices and
returns an integer permutation of the fold columns. At seven columns, checking
all 7! assignments is small and avoids a new assignment-package dependency.
Constant-score correlations contribute zero to the assignment. `c_paired()`
returns a two-element numeric vector containing the mean paired difference and
its SE, or two undefined values if any fold is missing.

`c_read()` retrieves a private checkpoint or returns NULL for deferred CV.
`c_summary()` joins fold concordances by their fold labels, aligns programs and
writes aggregate descriptions. `c_external()` uses the existing clinical and
molecular data helpers to score the fixed fits, then writes population-level
concordances. `survival` supplies concordance through the shared orientation
helper, `yaml` reads the existing model configuration, and the fitter retains
its established `ebnm` empirical Bayes updates. No new analysis dependency or
model option was introduced.

## Verification and limitations

Both focused synthetic test scripts passed before the stop request:
`tests/test_methylation_single_modality.R` and
`tests/test_methylation_iteration_extension.R`. Coverage includes missing and
duplicate CV results, varied arithmetic/permutation fixtures, isolated cap
configuration, source/config mismatch rejection, divergent ELBO prefixes,
constant-score reporting and checkpoint-preserving resume. No real 600-iteration
fit was run. Independent reviews found no remaining launch blocker after the
reported fixes. The real-data smoke fits converged in 61 PL and 185 PE
iterations. Their single-fold concordances, 0.615 and 0.571, are smoke diagnostics
and are not estimates from the registered five-fold comparison.

No file under `code/` was changed by this session, so the prompt's conditional
full-suite run was not triggered. Other sessions' model edits were excluded
from the private model snapshot.

KNN imputation is inherited from the prepared full-TCGA universe, rather than
refitted within each fold. The registered screening also precedes the folds.
Only scaling and model fitting are estimated within each training complement.
Any future CV results would describe this fixed prepared-data workflow.

Synthetic tests cover scaling arithmetic and application to held-out rows,
reproducible event-balanced folds, one-to-one program matching under column
permutation and sign reversal, paired SE, undefined-fold propagation and exact
nonzero-loading/activity definitions. They also exercise subprocess failure
recording and resumption without a memory log.

This study was specified with Normal β while task B was pending. The
[β-prior comparison](2026-11-12-beta-prior-comparison-findings.md) has since
retained Normal, including its real-data K = 5 rerun. The original results remain
at K = 7. The updated working-rank comparison is reported below. The
registered recovery, false-active-program and convergence gate remains
necessary for scientific adoption. This descriptive comparison does not adopt
a model or alter that gate.

## Would more iterations help?

Andrew asked during the sweep whether the larger fits should receive more than
300 iterations. Saved late-iteration histories give different answers by setup.
The final relative ELBO changes were 1.41 × 10⁻⁵ for 50k PL, 1.54 × 10⁻⁴ for
50k PE and 2.83 × 10⁻⁵ for 100k PE, against the fixed tolerance of 10⁻⁵. None
had an ELBO decrease in its final 50 iterations. The 50k PL change was slowly
declining toward tolerance. The PE changes were roughly flat or increasing
recently, so these histories do not establish how many more iterations they
would need. The 100k PL fit already converged in 150 iterations.

A separately labeled 600-iteration check of 50k PL is a reasonable numerical
follow-up, with the same initialization and all other settings fixed. Compare
its risk scores and program identities with the saved 300-iteration fit, in
addition to the convergence flag. Andrew approved this numerical follow-up,
then prioritized the K = 5
handoff and deferred its execution. The script is prepared and its synthetic
tests pass, but no 600-iteration fit has been run. It will use the frozen K = 7
model in a separate private source/configuration and checkpoint directory. No
additional CV is authorized by this check. The registered 300-iteration results
remain unchanged.
See `iteration_cap_diagnostics.csv` for the saved summaries. If Andrew later
resumes this deferred diagnostic, use:

```sh
Rscript -e 'source("results/multimodal_real/run_methylation_iteration_extension.R")'
```

The extension script copies the frozen input and model code to a separate
private directory and doubles only the iteration cap to 600. It reruns from
the same seed and initialization. Before comparing results, it requires the
first 300 ELBO values to reproduce the baseline trajectory to relative error
at most 10⁻¹⁰. It reports the signed correlation of risk scores, their RMS
change divided by the baseline score SD, and matched program correlations and
activity indicators. Patient-level outputs remain private.

## Handoff to the K = 5 session

Andrew stopped this K = 7 run and prioritized the working-rank K = 5 study.
There are no active processes from this session. The runner and shared files
are available for handoff. Preserve the K = 7 outputs and private checkpoints.

Use `data/methylation_single_modality_2026-11-12/design.rds` for the fixed
patients, D4 genes, eligible probe ranking and saved folds. Model code and
configuration used here are in its private `source/` directory. The original
`k_merged_yfb_desurv` remains seven, while the current working key is five.
K = 5 needs separate checkpoint/output keys and explicit rank labels.

**No count completed official K = 7 CV.** The K = 5 prompt's condition referring
to counts with completed K = 7 CV cannot supply any methylation counts. Its
remaining protocol must account for this user-requested stop rather than
waiting for or treating absent K = 7 fold results as completed. The saved
full-fit timing/convergence table is available for revising that scope. No
K = 5 fit was launched or modified by this session.

## Text for the shared status files

**Attempts log:** the registered K = 7 methylation-only comparison verified
182,934 eligible probes on 144 matched TCGA patients and 75 events. Nine full
fits completed, six converged and three capped. Andrew stopped the final
all-CpG PE fit after 46.0 minutes to free compute for the updated K = 5 study.
No expression fits, official CV or ICGC scoring ran. Preserve completed private
checkpoints and aggregate timing/convergence evidence. The original scientific
comparison is incomplete. The 600-iteration diagnostic is prepared and tested
but deferred at Andrew's request.

**DECISIONS:** retain these results as historical K = 7 numerical diagnostics.
They do not establish relative accuracy or methylation-only programs and do
not change model adoption. Prioritize the separate K = 5 study. Its protocol
must account for the absence of completed K = 7 CV. Normal β was fixed before
task B and has since been retained by that comparison.

Shared status files were not edited by this session.

## Rerun at K = 5

The K = 7 run above was stopped at Andrew's request. This follow-up uses
K = 5, preserving the K = 7 checkpoints and aggregate results. The same
144 matched TCGA patients, D4 genes, eligible probe ranking and saved folds
are reused. Full fits and each training fold retain the original scaling,
prior, initialization, iteration cap and scoring definitions.

### Revised scope and computing screen

No official K = 7 methylation CV completed before the stop request. The
original K = 5 prompt's completed-CV prerequisite therefore yields no eligible
methylation count. Andrew subsequently requested K = 5 full fits through
50k, 100k and all eligible CpGs, followed by the three-hour CV screen. The
resumed scope is therefore all five original CpG counts at K = 5, alongside
the complete matched expression reference. This is a change to the execution
protocol, not a completed-CV claim about K = 7.

The full fits run serially in increasing CpG count. The original computing
rule is then applied to measured K = 5 process times. For each count, five
times the sum of its PL and PE full-fit process times projects the combined
five-fold CV cost. Both setups contribute to that timing screen. CV runs
only for converged full-fit setups at a count whose projection is at most
three hours. A cap is a failure to converge, rather than an invitation to
change the iteration budget. The approved K = 7 iteration extension remains
deferred. No K = 5 extension is included.

### Matched expression reference

All 12 expression fits are saved, with no fitting errors or fit warnings.
PL converged in its full fit and two of five folds. PE's full fit and four
folds reached the fixed 300-iteration cap. Concordance from capped states
remains labeled as such and does not establish convergence.

| K | Setup | Mean CV C | Defined CV folds | Converged CV folds | Mean nonzero programs | Mean survival-active | Full status | Full nonzero | Full survival-active |
|---|---|---|---|---|---|---|---|---|---|
| 5 | PL | 0.591160 | 5/5 | 2/5 | 5 | 1 | converged | 5 | 1 |
| 5 | PE | 0.622361 | 5/5 | 1/5 | 5 | 1 | capped | 5 | 1 |
| 7 | PL | not run | 0/5 | not available | not available | not available | not run | not available | not available |
| 7 | PE | not run | 0/5 | not available | not available | not available | not run | not available | not available |

These estimates describe the matched TCGA subset, not the pooled TCGA+CPTAC
sample used in the β-prior comparison. The higher PE concordance does not
make it a converged solution. Its capped full fit also prevents an assessable
PE methylation-only program classification under the registered rule.

### Methylation full fits at both ranks

All five K = 5 PL full fits converged. PE converged at 100k and all CpGs,
while the three smaller PE full fits capped. None returned a fitting error
or fit warning. Every K = 5 full fit retains five programs with nonzero
loadings. The nine completed K = 7 full fits retained seven. These are
nonzero loading columns, rather than validated biological programs.

| CpGs | Setup | K = 7 status | K = 7 active | K = 5 status | K = 5 iterations | K = 5 active | K = 5 process minutes | K = 5 peak RSS MiB |
|---:|:---|:---|---:|:---|---:|---:|---:|---:|
| 10,000 | PL | converged | 0 | converged | 112 | 1 | 1.11 | 962 |
| 10,000 | PE | converged | 0 | capped | 300 | 0 | 2.37 | 880 |
| 25,000 | PL | converged | 0 | converged | 98 | 0 | 2.32 | 1161 |
| 25,000 | PE | converged | 0 | capped | 300 | 0 | 6.27 | 1183 |
| 50,000 | PL | capped | 0 | converged | 188 | 0 | 6.73 | 2023 |
| 50,000 | PE | capped | 0 | capped | 300 | 0 | 7.44 | 2225 |
| 100,000 | PL | converged | 1 | converged | 122 | 0 | 7.30 | 3038 |
| 100,000 | PE | capped | 2 | converged | 243 | 0 | 11.40 | 3164 |
| 182,934 | PL | converged | 0 | converged | 94 | 0 | 11.01 | 4037 |
| 182,934 | PE | stopped_by_user | NA | converged | 253 | 0 | 24.56 | 4450 |

The K = 7 stopped row has no fitted-state activity count. At K = 5, only
the 10k PL full fit is survival-active, with one program. Full-fit activity
alone does not establish a stable methylation-only program. The CV and
expression-correlation requirements still apply.

| CpGs | K = 5 projected CV hours, both setups | Setups scheduled for CV |
|---:|---:|:---|
| 10,000 | 0.291 | PL |
| 25,000 | 0.716 | PL |
| 50,000 | 1.181 | PL |
| 100,000 | 1.558 | PL, PE |
| 182,934 | 2.964 | PL, PE |

All counts passed the per-count screen. All CpGs passed by only 129.8
seconds in the projected five-fold total. This is an eligibility forecast,
not a guarantee that actual CV will finish within three hours. The scheduled
work comprises 25 PL folds and ten PE folds, with no PE CV at 10k, 25k or
50k because those full fits capped.

### Paired CV results at K = 5

All 35 scheduled methylation CV fits completed, with 34 converged and one
capped at 100k PE. Every scheduled fold has defined concordance. The three
smaller PE setups received no CV because their full fits capped. No fold
was dropped from a reported five-fold mean or paired SE.

| CpGs | Setup | Mean CV C | Difference from expression (paired SE) | Converged CV folds | Methylation-only programs |
|---:|:---|---:|---:|:---|:---|
| 10,000 | PL | 0.545 | −0.046 (0.041) | 5/5 | 0 |
| 25,000 | PL | 0.504 | −0.088 (0.038) | 5/5 | 0 |
| 50,000 | PL | 0.511 | −0.080 (0.044) | 5/5 | 0 |
| 100,000 | PL | 0.502 | −0.089 (0.038) | 5/5 | 0 |
| 182,934 | PL | 0.520 | −0.071 (0.038) | 5/5 | 0 |
| 10,000 | PE | not run | not run | not run | not assessed |
| 25,000 | PE | not run | not run | not run | not assessed |
| 50,000 | PE | not run | not run | not run | not assessed |
| 100,000 | PE | 0.539 | −0.083 (0.053) | 4/5 | not assessed |
| 182,934 | PE | 0.512 | −0.111 (0.052) | 5/5 | not assessed |

The matched expression references are PL C = 0.591 and PE C = 0.622. Only
two PL expression folds and one PE expression fold converged, so the paired
comparisons include predictions from capped expression states.
They describe this run rather than the prediction of fully converged models.
The SE is the sample SD of the five fold differences divided by √5.

No evaluated methylation setup has a higher mean CV C than its expression
reference. The estimates remain noisy at n = 144. This comparison does not
rule out useful methylation information in a different model or establish
that a multimodal model would fail. No rank, prior or CpG count is adopted
from these results.

K = 7 has no official CV or expression fits, so a paired K = 5 versus K = 7
CV comparison is unavailable. The full-fit comparison above preserves the
actual historical results rather than imputing missing K = 7 estimates.

The actual selected CV process times, summed across folds and setups, were
0.061, 0.182, 0.596, 1.348 and 2.255 hours at 10k, 25k, 50k, 100k and all
CpGs. Smaller counts include PL only. The all-CpG two-setup CV completed
below its 2.964-hour projection. These summed process times do not assume
parallel speedup.

### Methylation-only program check

No PL program meets the methylation-only rule at K = 5. The only
survival-active PL full-fit program is program 5 at 10k CpGs, with
E[β] = 0.2639 and posterior SD 0.1313. It is active in only one of the five
converged matched folds, below the required three. Its largest absolute score
correlation with an expression program is 0.5957, above the required limit
of 0.5. It fails both conditions. The smallest matched full-to-fold score
correlation is 0.2122, so program correspondence is weak in at least one fold.
The other PL full fits have no survival-active programs.

The converged PE full fits at 100k and all CpGs also have no survival-active
programs. Formal PE methylation-only classification remains undefined because
the expression PE full fit capped. The three smaller methylation PE full fits
also capped and received no CV. These unassessed cases are kept separate from
the zero qualifying-program count in the assessable PL comparisons.

### Rebuild commands and implementation

From the repository root:

```sh
OMP_NUM_THREADS=1 Rscript results/multimodal_real/run_methylation_single_modality.R rerun --K=5
OMP_NUM_THREADS=1 Rscript results/multimodal_real/run_methylation_single_modality.R sweep --K=5
Rscript tests/test_methylation_k5_rerun.R
```

The first command completes the 12 expression fits. The second resumes those
checkpoints and runs the full methylation sweep with the measured-time CV
screen. Each full count has both setups. Every model process uses one numerical
thread, with at most three concurrent CV fits. Full methylation fits and CV
at 100k/all CpGs run serially. Neither command performs external scoring or
model selection.

- `c_set_rank()` takes an optional integer K. The historical default remains
  `benchmark$k_merged_yfb_desurv`, which is still seven. K = 5 selects separate
  `K05/` directories beneath the private and aggregate study directories. The
  private design is linked to the original and its hash must agree, preserving
  the patient order, screening and folds.
- For this separate rank, `c_freeze()` copies the original hash-verified model
  code and config, together with the runner. The historical config remains
  unchanged. Its source revision and copied source hashes are written beside
  the K = 5 aggregates. The actual child runner is immutable, including its
  rank-aware fitting entry point. Later parent-orchestration edits do not
  change it. Concurrent fitter edits in other sessions cannot enter these
  comparisons.
- `c_launch()` passes `--K=5` to each fresh child. `c_one()` uses that explicit
  rank while retaining the saved model settings. Private checkpoints contain
  posterior fits, scaling vectors and patient scores. Aggregate metric rows
  carry the actual fitted rank. Resume reads only the selected rank's directory.
- `c_rerun()` takes no arguments and completes the matched expression reference.
  It verifies complete, unique coverage of the ten held-out jobs and retains
  failed or undefined results explicitly. `expression_summary.csv` reports
  full-fit and fold program counts, mean CV C and convergence. `sweep` reuses
  `c_run()` for the original five-count design and actual-time CV screen, then
  calls `c_summary()` without external scoring. Failed fits cause a nonzero
  exit after aggregate diagnostics are saved.
- `c_summary()` uses the existing one-to-one assignment of full and fold
  programs by absolute score correlation. A methylation-only program must
  satisfy the full-fit activity criterion, remain active in at least three
  converged matched folds, and have maximum absolute expression-program
  correlation below 0.5. Classification is undefined when the required full
  fits or fold coverage do not support assessment. Constant scores cannot
  substantiate a program match.
- `tests/test_methylation_k5_rerun.R` uses temporary synthetic files to verify
  rank isolation, frozen-source integrity, child rank propagation, all-count
  full-fit coverage, measured-time CV gating, unique fold coverage and
  expression summaries. It starts no model fits. The original helper and
  iteration-extension tests also pass. Base R handles file and matrix
  operations, yaml reads the established config, and parallel orchestrates
  fresh children.

Patient-level K = 5 outputs remain under ignored
`data/methylation_single_modality_2026-11-12/K05/`. Only aggregate CSVs are
written to `results/multimodal_real/outputs/methylation_single_modality_2026-11-12/K05/`.

The full regression suite passed 819 tests with zero failures. The focused
rank-five rerun, original helper and iteration-extension tests also passed.
The rerun fixture uses different PL and PE mean concordances, so a constant
summary cannot satisfy the numerical check. Independent code review found no
blocker in the expanded sweep and confirmed that the frozen child fitting
functions match the reviewed implementation.

The completed run exited with status zero. The final checkpoint audit verified
57 distinct fits, all saved features and fold indices, rank five, risk scores
from the stored program projections, every held-out concordance recomputed
under the frozen Cox direction, the paired means and SEs, and the program rule.
It also verified the actual-time eligibility screen and complete scheduled CV
coverage. It found no fitting errors or fit warnings. The 12 capped fits are
eight expression fits and four methylation fits, including three methylation
full fits and one methylation CV fold.

Aggregate evidence is in `full_fits.csv`, `cv_fits.csv`, `cv_summary.csv`,
`programs.csv`, `fold_program_matching.csv` and `checkpoint_audit.csv` under
`K05/`. `rank_comparison_full.csv` combines shared columns of the preserved
K = 7 and K = 5 full-fit tables. `rerun_input_hashes.csv` records the saved
design and parent runner hashes. No patient-level inputs, fits or scores were
committed. No ICGC scoring was run in this follow-up.

Independent final review found no substantive issue or blocker. It recomputed
all 45 held-out concordances, the seven paired means and SEs, and every
full-to-fold program assignment by independently enumerating the 120 possible
rank-five permutations. It confirmed the program classifications, timing
screen, source and input integrity, K = 7 preservation and aggregate-only
outputs. The review used saved fits and process measurements rather than
refitting models or remeasuring runtime and memory.

## Text for the shared status files (K = 5)

**Attempts log:** the matched-TCGA methylation comparison was repeated at K = 5
on the same 144 patients, 75 events, five saved folds and 182,934 eligible
probes. Andrew expanded the stopped K = 7 handoff to full fits at all five
CpG counts, followed by the original per-count three-hour CV screen using
measured K = 5 timings. All ten methylation full fits completed, with seven
converged and three capped. All counts passed the timing screen. The 35
scheduled methylation CV fits completed, with 34 converged and one capped.
The expression reference contributed 12 fits. All 57 fits were audited, with
no fit errors or warnings. The source snapshot and K = 7 results were preserved.

**Findings:** mean PL methylation CV C was 0.545, 0.504, 0.511, 0.502 and
0.520 across increasing counts, against 0.591 for expression. PE was 0.539
at 100k and 0.512 at all CpGs, against 0.622 for expression. The paired SEs
and convergence limits are in the findings table. No PL program met the
methylation-only rule. The sole survival-active PL full-fit program, 10k
program 5, was active in only one fold and correlated 0.596 with expression.
PE full fits at the two CV-eligible setups had no survival-active programs,
and formal uniqueness assessment is undefined because expression PE capped.

**DECISIONS:** no model, prior or CpG count is adopted from this descriptive
comparison. Normal β remains the specified prior. Expression convergence,
pre-fold screening and imputation, and the small matched sample limit the
interpretation. Reliable supervision still requires the registered recovery,
false-active-program and convergence gate. Shared status files and the
deferred K = 7 iteration-extension study were not changed.
