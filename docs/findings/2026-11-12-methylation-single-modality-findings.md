# Methylation-only programs on matched TCGA patients

Date: 2026-11-12. Author: Andrew Walther.

**K = 5 follow-up.** All eight capped expression fits converged after the approved
1,000-sweep extension. K = 5 expression mean CV C is now 0.582311 for PL and
0.620214 for PE, and methylation remains lower in all seven completed comparisons.
No K = 5 prognostic methylation-only program qualifies. The descriptive K = 5
analysis identifies 19 of 35 factor-at-count instances as methylation-specific
under the fixed correlation, stability and recurrence thresholds. These K = 5
results, purity associations and caveats appear in the final follow-up section.

**K = 5 result at the original 300-sweep cap.** The full sweep and timing-eligible CV are complete. Mean
methylation concordance is below matched expression in all seven completed
comparisons. No PL program meets the methylation-only rule. No PE program can
meet it either, because neither converged PE methylation full fit (100k, all
CpGs) has a survival-active program. Only the expression-correlation step of the
PE rule is unassessable, because the expression PE full fit capped. Of the 57
fits, 12 capped, including 8 of the 12 expression fits, so the paired
differences use capped expression states in most folds.
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

## Full-fit results (K = 7)

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

## Text for the shared status files (K = 7, superseded by the K = 5 text below)

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
the expression PE full fit capped. That limit affects only the
expression-correlation condition. The rule also requires full-fit survival
activity, which neither converged PE full fit has, so no PE program at 100k or
all CpGs could qualify whatever the expression fit. The three smaller
methylation PE full fits capped and received no CV, so they remain
unassessed. These cases are kept separate from the zero qualifying-program
count in the assessable PL comparisons.

In the two folds where the PL expression fit converged (folds 2 and 5),
methylation PL concordance was below expression at every CpG count except
10k in fold 2 (0.708 against 0.685). Fold-level concordance varies widely,
from 0.30 to 0.71 for methylation, with fold 5 low for both modalities. This
check is descriptive and was not registered.

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


## K = 5 expression reruns and factor characterization

### Fixed K = 5 descriptive rules before examining results

At K = 5, a methylation factor is called methylation-specific only when all three
conditions hold: its largest absolute correlation with any expression factor in the
same setup is below 0.5, its one-to-one matched full-to-fold correlation is at least
0.5 in at least three converged CV folds, and its largest absolute correlation with
a full-fit factor at another CpG count in the same setup is at least 0.7. These
thresholds are fixed before this follow-up's results are inspected. At K = 5, this
descriptive definition does not require survival activity and does not replace the
earlier methylation-only prognostic-program rule.

At K = 5, the eight capped expression fits restart from the original initialization
with a 1,000-sweep cap. Every earlier ELBO value must reproduce within relative
error 10⁻¹⁰ before a refit is used. At K = 5, all other frozen model settings and
the saved 144-patient design remain unchanged. The K = 5 characterization uses saved
posterior means, the 2,218-gene D4 expression universe, and no external outcomes.

### K = 5 expression convergence and paired CV

At K = 5, all eight expression restarts converged in 305 to 438 sweeps. The first
300 ELBO values reproduced exactly for every K = 5 restart, with maximum relative
error zero. At K = 5, every final relative ELBO change was below the fixed tolerance
10⁻⁵, and no restart produced a fitting warning or error. Four already-converged
K = 5 fits were reused unchanged, so both K = 5 expression references now have a
converged full fit and five converged CV folds.

| K | Setup | Fold, 0 = full | Sweeps | Final relative change | Risk r with 300 sweeps | Minimum matched program abs(r) | CV C, 300 | CV C, extended |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 5 | PL | 1 | 404 | 9.8416936e-06 | 0.945589 | 0.950295 | 0.694690 | 0.650442 |
| 5 | PL | 3 | 326 | 9.9699771e-06 | 0.998764 | 0.999014 | 0.552124 | 0.552124 |
| 5 | PL | 4 | 305 | 9.9465965e-06 | 0.999999 | 0.999934 | 0.593750 | 0.593750 |
| 5 | PE | 0 | 438 | 9.9789158e-06 | 0.999976 | 0.996190 | NA | NA |
| 5 | PE | 1 | 389 | 9.9097006e-06 | 0.999531 | 0.997946 | 0.663717 | 0.650442 |
| 5 | PE | 3 | 393 | 9.9620203e-06 | 0.999857 | 0.997882 | 0.563707 | 0.571429 |
| 5 | PE | 4 | 328 | 9.9545094e-06 | 0.999923 | 0.999653 | 0.629464 | 0.629464 |
| 5 | PE | 5 | 404 | 9.9751311e-06 | 0.995155 | 0.998895 | 0.611399 | 0.606218 |

At K = 5, the point-Laplace expression mean CV C changed from 0.59116016 to
0.58231061, and the point-exponential expression mean changed from 0.62236103 to
0.62021429. The K = 5 point-Laplace fold 1 change accounts for the point-Laplace mean
change. K = 5 risk-score correlations can remain high while individual comparable
pairs change order, so convergence alone does not imply identical concordance.

| K | Setup | CpGs | Methylation mean C | Paired Δ at 300 sweeps | SE at 300 | Paired Δ after extension | SE after extension |
|---|---|---:|---:|---:|---:|---:|---:|
| 5 | PL | 10,000 | 0.545442 | -0.045718 | 0.040687 | -0.036869 | 0.039548 |
| 5 | PL | 25,000 | 0.503502 | -0.087658 | 0.037697 | -0.078808 | 0.037648 |
| 5 | PL | 50,000 | 0.511267 | -0.079893 | 0.044426 | -0.071043 | 0.044876 |
| 5 | PL | 100,000 | 0.502136 | -0.089025 | 0.038474 | -0.080175 | 0.035598 |
| 5 | PL | 182,934 | 0.519814 | -0.071346 | 0.038021 | -0.062497 | 0.032802 |
| 5 | PE | 100,000 | 0.539498 | -0.082863 | 0.053348 | -0.080716 | 0.052502 |
| 5 | PE | 182,934 | 0.511707 | -0.110654 | 0.052326 | -0.108507 | 0.051225 |

At K = 5, Δ is methylation minus expression on the same five folds and SE is the
sample SD of those five differences divided by √5. All seven completed K = 5
methylation comparisons remain lower than expression. At K = 5, the 10k, 25k and
50k point-exponential full methylation fits remain capped at 300 sweeps and have
no CV. The K = 5 100k point-exponential methylation CV still contains one capped fold.
The K = 5 extension changes expression only, and preserves every 300-sweep table.

At K = 5, the point-exponential expression full fit converged at sweep 438, so the
expression-correlation step of the earlier prognostic-program rule is now assessable
for both converged point-exponential methylation full fits. Neither K = 5 fit has a
survival-active program, and neither can meet that rule regardless of correlation.
Across all 35 K = 5 full-fit methylation factors, zero meet the prognostic-program
rule. See `expression_iteration1000.csv`, `expression_iteration1000_programs.csv`,
`cv_summary_iteration1000.csv` and `programs_iteration1000.csv` under the K05 output
directory for the K = 5 fit-level and per-program values.

### K = 5 factor overlap, annotation and expression biology

At K = 5, 19 of the 35 methylation factors meet the fixed descriptive specificity
rule: nine point-Laplace factors and all ten point-exponential factors. These K = 5
counts include recurring factors at different CpG counts and are not counts of
19 independent biological processes. At K = 5, low pairwise correlations with
individual expression factors do not establish separation from every linear
combination of expression factors. The K = 5 canonical correlations quantify that
distinction. Across the seven K = 5 comparisons, the first canonical correlation
is 0.903 to 0.921 and the second is 0.849 to 0.877, while the fifth is 0.002 to
0.174. At K = 5, the modalities share some score directions and have other directions
with little observed linear overlap. These K = 5 correlations are descriptive
within the same 144 patients, without held-out validation of factor identity.

At K = 5, the annotation package is
`IlluminaHumanMethylation450kanno.ilmn12.hg19` version 0.6.1. All 182,934 eligible
K = 5 CpGs were found in its hg19 annotation. At K = 5, each factor's top 100 CpGs
are ranked by the absolute posterior mean gene-feature loading `EF`. The analogous
K = 5 expression top-100 genes also use `EF`, as Andrew clarified. In the saved
K = 5 model, `EL` is a 144-patient loading matrix and `EF` is a feature-by-five
loading matrix. At K = 5, multi-gene annotation entries remain separate gene-region
records, while annotation fractions and CpG-overlap counts count each CpG once.

At K = 5, the eligible-CpG background is 20.91% island, 25.73% shore, 10.53%
shelf and 42.84% open sea, with 32.89% promoter and 41.30% gene-body annotation.
At K = 5, top-100 island fractions range from 0% to 70% for point-Laplace and
0% to 35% for point-exponential, while promoter fractions range from 11% to 59%
and 13% to 38%, respectively. The K = 5 factors therefore differ substantially
in genomic context rather than following a single common annotation profile.

The K = 5 table below reports every methylation factor. At K = 5, stability is the
number of converged matched folds with |r| ≥ 0.5, recurrence is the largest |r|
with a factor at another CpG count in the same setup, and specific applies the
three fixed conditions. At K = 5, the annotation column gives the most frequent
island category and the promoter fraction among the top 100 CpGs. The K = 5
Hallmark column gives the smallest adjusted p-value within that factor's Hallmark
analysis, with signed NES and q. All K = 5 purity correlations use ABSOLUTE purity
on 143 complete patients because one of the 144 estimates is missing.

| K | Setup | Methylation factor at CpGs | Closest expression factor, abs(r) | Stable folds | Recurrence | Specific | Top-CpG annotation | Leading Hallmark, NES, q | Purity abs(r) |
|---|---|---|---|---:|---:|---|---|---|---:|
| 5 | PL | factor 1 at 10,000 | 1, 0.653 | 5 | 0.993 | no | island 68%, promoter 59% | EPITHELIAL_MESENCHYMAL_TRANSITION, 2.19, 1.9e-06 | 0.749 |
| 5 | PL | factor 2 at 10,000 | 1, 0.548 | 5 | 0.835 | no | open sea 81%, promoter 25% | ALLOGRAFT_REJECTION, 2.55, 3.2e-08 | 0.707 |
| 5 | PL | factor 3 at 10,000 | 1, 0.583 | 5 | 0.967 | no | open sea 58%, promoter 39% | EPITHELIAL_MESENCHYMAL_TRANSITION, 3.11, 2.3e-24 | 0.507 |
| 5 | PL | factor 4 at 10,000 | 2, 0.436 | 5 | 0.837 | yes | open sea 64%, promoter 22% | EPITHELIAL_MESENCHYMAL_TRANSITION, 3.17, 6.7e-25 | 0.107 |
| 5 | PL | factor 5 at 10,000 | 2, 0.596 | 4 | 0.873 | no | open sea 71%, promoter 19% | EPITHELIAL_MESENCHYMAL_TRANSITION, 2.45, 9.3e-29 | 0.089 |
| 5 | PL | factor 1 at 25,000 | 1, 0.698 | 5 | 0.997 | no | island 69%, promoter 47% | EPITHELIAL_MESENCHYMAL_TRANSITION, 2.21, 5.3e-07 | 0.799 |
| 5 | PL | factor 2 at 25,000 | 1, 0.573 | 5 | 0.904 | no | open sea 65%, promoter 27% | EPITHELIAL_MESENCHYMAL_TRANSITION, 2.44, 5.9e-11 | 0.595 |
| 5 | PL | factor 3 at 25,000 | 4, 0.726 | 5 | 0.844 | no | open sea 69%, promoter 39% | EPITHELIAL_MESENCHYMAL_TRANSITION, 3.11, 6.2e-24 | 0.179 |
| 5 | PL | factor 4 at 25,000 | 2, 0.436 | 4 | 0.948 | yes | open sea 67%, promoter 29% | EPITHELIAL_MESENCHYMAL_TRANSITION, -2.23, 6.1e-14 | 0.063 |
| 5 | PL | factor 5 at 25,000 | 3, 0.519 | 4 | 0.910 | no | open sea 67%, promoter 29% | EPITHELIAL_MESENCHYMAL_TRANSITION, -2.97, 9.5e-19 | 0.316 |
| 5 | PL | factor 1 at 50,000 | 1, 0.721 | 5 | 0.997 | no | island 70%, promoter 47% | EPITHELIAL_MESENCHYMAL_TRANSITION, 2.21, 4.7e-07 | 0.829 |
| 5 | PL | factor 2 at 50,000 | 1, 0.414 | 5 | 0.904 | yes | open sea 57%, promoter 12% | ALLOGRAFT_REJECTION, 2.52, 1.4e-09 | 0.480 |
| 5 | PL | factor 3 at 50,000 | 4, 0.710 | 5 | 0.977 | no | open sea 68%, promoter 39% | EPITHELIAL_MESENCHYMAL_TRANSITION, 3.28, 2.5e-26 | 0.236 |
| 5 | PL | factor 4 at 50,000 | 2, 0.401 | 5 | 0.948 | yes | open sea 72%, promoter 31% | EPITHELIAL_MESENCHYMAL_TRANSITION, -2.13, 1.5e-08 | 0.198 |
| 5 | PL | factor 5 at 50,000 | 1, 0.489 | 5 | 0.910 | yes | open sea 67%, promoter 23% | EPITHELIAL_MESENCHYMAL_TRANSITION, 3.41, 4.8e-27 | 0.370 |
| 5 | PL | factor 1 at 100,000 | 1, 0.748 | 5 | 0.994 | no | island 68%, promoter 44% | EPITHELIAL_MESENCHYMAL_TRANSITION, 2.22, 2.8e-07 | 0.861 |
| 5 | PL | factor 2 at 100,000 | 1, 0.125 | 5 | 0.963 | yes | open sea 54%, promoter 11% | ALLOGRAFT_REJECTION, 2.12, 3.8e-05 | 0.198 |
| 5 | PL | factor 3 at 100,000 | 4, 0.684 | 5 | 0.977 | no | open sea 67%, promoter 44% | EPITHELIAL_MESENCHYMAL_TRANSITION, 3.14, 1.6e-23 | 0.319 |
| 5 | PL | factor 4 at 100,000 | 2, 0.546 | 5 | 0.929 | no | open sea 73%, promoter 30% | EPITHELIAL_MESENCHYMAL_TRANSITION, -2.40, 1.3e-11 | 0.268 |
| 5 | PL | factor 5 at 100,000 | 1, 0.359 | 3 | 0.968 | yes | open sea 72%, promoter 22% | EPITHELIAL_MESENCHYMAL_TRANSITION, 2.21, 1.1e-09 | 0.202 |
| 5 | PL | factor 1 at 182,934 | 1, 0.774 | 5 | 0.992 | no | island 45%, promoter 43% | EPITHELIAL_MESENCHYMAL_TRANSITION, 2.27, 3e-08 | 0.872 |
| 5 | PL | factor 2 at 182,934 | 2, 0.096 | 5 | 0.963 | yes | open sea 55%, promoter 16% | TNFA_SIGNALING_VIA_NFKB, 1.91, 9.4e-05 | 0.150 |
| 5 | PL | factor 3 at 182,934 | 4, 0.759 | 5 | 0.965 | no | open sea 69%, promoter 37% | EPITHELIAL_MESENCHYMAL_TRANSITION, 2.95, 1.2e-19 | 0.151 |
| 5 | PL | factor 4 at 182,934 | 1, 0.288 | 5 | 0.968 | yes | open sea 75%, promoter 24% | EPITHELIAL_MESENCHYMAL_TRANSITION, -1.97, 6.5e-06 | 0.160 |
| 5 | PL | factor 5 at 182,934 | 2, 0.581 | 5 | 0.929 | no | open sea 75%, promoter 23% | EPITHELIAL_MESENCHYMAL_TRANSITION, 2.76, 5.7e-32 | 0.048 |
| 5 | PE | factor 1 at 100,000 | 2, 0.267 | 4 | 0.872 | yes | open sea 76%, promoter 26% | ALLOGRAFT_REJECTION, 2.17, 2.5e-05 | 0.248 |
| 5 | PE | factor 2 at 100,000 | 4, 0.406 | 4 | 0.975 | yes | open sea 41%, promoter 31% | EPITHELIAL_MESENCHYMAL_TRANSITION, -2.81, 6e-15 | 0.506 |
| 5 | PE | factor 3 at 100,000 | 4, 0.434 | 4 | 0.941 | yes | open sea 55%, promoter 38% | ALLOGRAFT_REJECTION, 2.69, 1.1e-10 | 0.575 |
| 5 | PE | factor 4 at 100,000 | 3, 0.404 | 4 | 0.941 | yes | open sea 71%, promoter 23% | EPITHELIAL_MESENCHYMAL_TRANSITION, -2.53, 1e-11 | 0.298 |
| 5 | PE | factor 5 at 100,000 | 4, 0.269 | 4 | 0.907 | yes | open sea 46%, promoter 31% | GLYCOLYSIS, 2.38, 1.9e-06 | 0.503 |
| 5 | PE | factor 1 at 182,934 | 4, 0.175 | 5 | 0.830 | yes | open sea 51%, promoter 13% | TNFA_SIGNALING_VIA_NFKB, 1.72, 0.0036 | 0.113 |
| 5 | PE | factor 2 at 182,934 | 4, 0.396 | 5 | 0.975 | yes | island 35%, promoter 30% | EPITHELIAL_MESENCHYMAL_TRANSITION, -2.72, 3.3e-13 | 0.443 |
| 5 | PE | factor 3 at 182,934 | 5, 0.253 | 5 | 0.941 | yes | open sea 34%, promoter 13% | ALLOGRAFT_REJECTION, 2.54, 1.3e-09 | 0.402 |
| 5 | PE | factor 4 at 182,934 | 3, 0.359 | 5 | 0.941 | yes | open sea 59%, promoter 31% | EPITHELIAL_MESENCHYMAL_TRANSITION, -2.61, 1.3e-12 | 0.214 |
| 5 | PE | factor 5 at 182,934 | 4, 0.296 | 5 | 0.870 | yes | open sea 47%, promoter 22% | EPITHELIAL_MESENCHYMAL_TRANSITION, -2.69, 1.4e-13 | 0.396 |

At K = 5, the expression factors receive the same D4 gene-score correlation and
Hallmark analysis, allowing direct descriptive comparison:

| K | Setup | Expression factor | Leading Hallmark | NES | q | Purity abs(r) |
|---|---|---:|---|---:|---:|---:|
| 5 | PL | 1 | EPITHELIAL_MESENCHYMAL_TRANSITION | 2.68 | 1e-14 | 0.734 |
| 5 | PL | 2 | EPITHELIAL_MESENCHYMAL_TRANSITION | 3.08 | 4e-32 | 0.118 |
| 5 | PL | 3 | MYC_TARGETS_V1 | 2.62 | 1.6e-05 | 0.231 |
| 5 | PL | 4 | APICAL_JUNCTION | -2.12 | 0.00032 | 0.267 |
| 5 | PL | 5 | EPITHELIAL_MESENCHYMAL_TRANSITION | -2.24 | 5.7e-26 | 0.278 |
| 5 | PE | 1 | EPITHELIAL_MESENCHYMAL_TRANSITION | 2.25 | 7.2e-17 | 0.181 |
| 5 | PE | 2 | EPITHELIAL_MESENCHYMAL_TRANSITION | 2.16 | 1.1e-17 | 0.505 |
| 5 | PE | 3 | INTERFERON_ALPHA_RESPONSE | -1.87 | 0.015 | 0.224 |
| 5 | PE | 4 | EPITHELIAL_MESENCHYMAL_TRANSITION | -2.48 | 3.2e-10 | 0.499 |
| 5 | PE | 5 | EPITHELIAL_MESENCHYMAL_TRANSITION | 2.18 | 2.1e-13 | 0.379 |

At K = 5, fgsea ranks the Pearson correlation of each factor's saved canonical
YF score with each of the 2,218 D4 genes. At K = 5, the Hallmark and PDAC
collections include retained Moffitt, Bailey and DeSurv sets. Each K = 5 gene set
is intersected with the D4 universe before testing, with the existing helper's
size limits of 10 to 500 genes. At K = 5, adjusted p-values apply within a factor
and collection, not across the complete analysis. The D4 universe and DeSurv
gene sets come from related expression preprocessing, so K = 5 enrichment is
not independent validation of a subtype or a causal mechanism. At K = 5, score
correlations with purity support a mixture-composition interpretation for some
factors, but do not identify a cell type or establish a tumor-intrinsic mechanism.

### K = 5 CpG–gene associations and overlap

At K = 5, the eligible universe contains 9,086 unique promoter CpG–D4-gene pairs.
Of these K = 5 pairs, 6,952 have negative Spearman correlation, a fraction of
0.76513317. At K = 5, 3,001 pairs are below −0.3 and 949 are below −0.5.
The strongest negative K = 5 pairs include cg25373430–GJB5, ρ = −0.864842,
cg02430183–WNK2, ρ = −0.819448, and cg03449125–CAPN5, ρ = −0.813785.
At K = 5, these are cross-patient associations after the saved imputation,
without adjustment for purity, multiple testing or potential confounding.

At K = 5, overlap of a methylation factor's top 100 CpGs with annotations to its
closest expression factor's top 100 genes ranges from zero to four CpGs. The
K = 5 null draws 100 CpGs without replacement from all 182,934 eligible CpGs,
1,000 times with seed 2026. Four K = 5 comparisons have unadjusted upper-tail
Monte Carlo p < 0.05: point-Laplace factor 5 at 10k, factor 4 at 25k, factor 5
at 25k, and factor 5 at all CpGs. The K = 5 observed counts are three, three,
three and four, with p = 0.039960, 0.039960, 0.009990 and 0.005994. None passes
a Bonferroni threshold of 0.05/35. At K = 5, the CpG annotation overlap therefore
does not establish gene-level correspondence beyond this descriptive screen.

At K = 5, `cross_modality_overlap.csv` also reports each factor's promoter-pair
count, negative fraction and counts below −0.3 and −0.5 against the complete
promoter-pair background. At K = 5, all 9,086 pair-level correlations appear in
`promoter_cpg_gene_correlations.csv`. Both are aggregate per-pair or per-factor
outputs, with no patient values.

### K = 5 implementation and rebuild explanation

At K = 5, `results/multimodal_real/run_expression_iteration_extension_k5.R` exists
to distinguish a capped numerical iterate from a converged expression reference.
Its K = 5 `e5_prepare()` block verifies the original source manifest, copies the
frozen source and design, and checks that only `cavi$max_iter` becomes 1000.
At K = 5, `e5_runner()` selects rank five through the existing `c_set_rank(5)`
mechanism and redirects private checkpoints to `iteration1000/`. At K = 5,
`e5_worker()` calls the unchanged fitter, which resets seed 42, and `e5_prefix()`
requires the saved ELBO prefix to reproduce before accepting the new posterior.
At K = 5, `e5_run()` launches batches of at most three child R processes and stops
after a batch if any child fails. At K = 5, `e5_summary()` reuses the four
converged fits, evaluates held-out scores with `frozen_reverse_cindex()`, and uses
`e5_paired()` to retain all five paired differences in the SE calculation.

At K = 5, the extension inputs are private saved fits and the frozen design, and
its outputs are eight new private fit checkpoints plus aggregate diagnostics,
matched program correlations, CV summaries and updated prognostic-program flags.
At K = 5, missing or nonfinite ELBO values, changed source hashes or any config
change beyond the cap are hard failures. At K = 5, this is a restart from the
original seed, not a posterior warm start. The K = 5 R implementation uses
`parallel::mclapply()` only to launch independent child processes, and `yaml` to
read and verify the configuration. No file under `code/` or the live global
configuration was changed by this follow-up.

At K = 5, `results/multimodal_real/characterize_methylation_factors_k5.R` reads
saved full-fit posterior means and the original fold-matching table without
refitting a model. Its K = 5 `classify_specific()` applies the three fixed
thresholds, `expand_annotation()` preserves aligned multi-gene and gene-region
annotations, and `overlap_permutation()` counts unique CpGs under a seeded
1,000-draw null. At K = 5, `safe_cor()` permits an undefined correlation for a
constant factor rather than inventing a value, and `island_class()` maps the
Illumina island labels to island, shore, shelf or open sea. At K = 5,
`characterize_k5()` checks feature and patient alignment, assembles the aggregate
outputs, and aligns local ABSOLUTE purity by unique saved sample ID.

At K = 5, the characterization inputs are seven converged methylation full fits,
two converged expression full fits, the saved 144-patient design and fold matches,
the approved hg19 annotation, the existing gene sets and local purity metadata.
At K = 5, the outputs are factor-level, CpG-annotation-level and CpG–gene-pair-level
CSV rows with no patient values. At K = 5, Pearson correlations describe factor
scores and gene-expression rankings, while Spearman correlations describe
promoter CpG–gene pairs. At K = 5, base R `cancor()` compares the two five-dimensional
score spaces, and the existing `run_fgsea_program()` uses signed rankings with
`scoreType = "std"` and seed 2026. At K = 5, promoter and gene-body categories can
overlap for a CpG with several annotations, and small D4 gene-set intersections
can be excluded by the existing fgsea size limits. The K = 5 annotation package
is installed in ignored `data/`, and package versions and input hashes are saved
alongside the aggregate characterization.

At K = 5, rebuild from the repository root with:

```sh
OMP_NUM_THREADS=1 Rscript results/multimodal_real/run_expression_iteration_extension_k5.R run
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 Rscript results/multimodal_real/characterize_methylation_factors_k5.R
Rscript tests/test_expression_iteration_extension_k5.R
Rscript tests/test_methylation_factor_characterization_k5.R
Rscript tests/run_tests.R
```

### K = 5 recommendation for a future joint fit

At K = 5, the recommended first CpG count for a future authorized joint fit is
10,000 with the point-Laplace setup. At K = 5, its full methylation fit and all
five CV folds converge, it already contains a stable recurrent descriptively
specific factor, and it has the highest methylation-only mean CV C among the
tested counts. At K = 5, larger counts add computational cost without improving
this marginal prediction comparison. This K = 5 recommendation does not establish
the best count for joint prediction, and it does not authorize a joint fit or
adopt a count. At K = 5, the point-exponential full fits at 10k, 25k and 50k remain
nonconverged, so the same recommendation does not extend to that setup.

### K = 5 enrichment warnings and numerical limits

At K = 5, the analysis writes all 2,475 factor-by-gene-set rows, including 18
rows with undefined NES, p-value and adjusted p-value because the fgsea simple
permutation calculation had unbalanced positive and negative gene statistics.
At K = 5, the 19 contextual warning records also include warnings about
overestimated p-values and the default 10⁻⁵⁰ p-value floor. These K = 5 warnings
are retained in `factor_characterization/analysis_warnings.csv`, and the
settings were not tuned to remove them. At K = 5, undefined rows are not counted
as evidence of enrichment or absence of enrichment. At K = 5, factor-specific
promoter-pair subsets can contain zero pairs, in which case their negative
fraction is undefined and their denominator is reported explicitly.

### K = 5 verification and independent review

At K = 5, `Rscript tests/run_tests.R` finishes with 849 passed and zero failed.
An earlier K = 5 full-suite run had 842 passed and five failures in concurrently
modified feature-group simulation tests. After that separate work was committed
in `a0bdb81d`, the full suite passed. The K = 5 verification log retains the
reported bounded-β profile-likelihood fallback diagnostic warning.

At K = 5, the focused extension and characterization tests each pass 15 assertions.
The K = 5 extension tests check matching and changed ELBO prefixes, a short or
nonfinite trajectory, cap-only configuration copying, preserved design contents
and paired means and SEs that change when the inputs change. The K = 5
characterization tests check planted shared and methylation-specific factors,
threshold boundaries, multi-gene annotations and seeded nonconstant permutation
results.

At K = 5, a fresh independent reviewer recomputed all expression concordances and
paired statistics from the private checkpoints, all 175 score correlations, all
35 canonical correlations using a separate QR/SVD calculation, and all 175 fold
assignments and 35 specificity calls. The K = 5 review also verified all top-100
CpG sets and annotations, all 9,086 promoter-pair correlations, all 35 overlap
nulls and all 45 factor-purity correlations. At K = 5, all 2,475 pathway sizes
were checked and 330 representative enrichment rows reproduced, including
undefined results. No substantive K = 5 numerical or code blocker was found.
The K = 5 review identified unescaped absolute-correlation symbols in table
headers, which were corrected.

### Verification of the follow-up and a joint-overlap check

All numbers in this follow-up were recomputed from the committed aggregate CSVs on
2026-11-12: the 12 converged expression fits, the expression means (0.5823 PL,
0.6202 PE), all seven paired differences and SEs, the 35 specificity calls, the
canonical-correlation ranges, the 9,086 promoter pairs and their counts, and the
maximum purity correlations. All reproduce. Two points qualify the specificity result.

First, only the expression-correlation condition decides the calls. Every one of the
35 methylation factors passes the fold-stability and cross-count recurrence conditions,
and the 16 non-specific factors all fail on max |r| ≥ 0.5.

Second, that condition compares a methylation factor with one expression factor at a
time, so a factor can pass while being largely a combination of several expression
factors. To check this, each methylation factor's full-fit scores were regressed on all
five expression factor scores of the same setup (the PL full fit, and the PE full fit
from the 1,000-sweep rerun), on the same 144 patients. The R² values are in
`factor_characterization/joint_expression_r2.csv`.

- Of the 19 factors called specific, 12 have R² < 0.25 and 3 have R² ≥ 0.5
  (PL factor 4 at 25k, 0.517; PL factor 4 at 50k, 0.527; PL factor 5 at 50k, 0.664).
  Those three are mostly explained by expression jointly.
- The clearest methylation-specific factor is PL factor 2 at 100k and at all CpGs
  (R² 0.055 and 0.036, purity |r| 0.20 and 0.15, recurrence |r| 0.963). Its leading
  Hallmark sets are allograft rejection and TNF-α signaling via NF-κB, which suggests
  an immune component, but the enrichment does not establish a cell type.
- PL factor 4 at all CpGs (R² 0.146) and PL factor 5 at 100k (R² 0.193) are also
  largely separate from expression. Most PE factors have R² of 0.04 to 0.24.
- The one specific factor at 10k (PL factor 4) has R² 0.431, so 10k carries no clearly
  separate methylation factor. This weakens the stated reason for recommending 10k for
  a joint fit if the aim of the joint fit is methylation-specific structure. The larger
  counts carry the most separate factors and still finish (full fit 7 to 11 minutes
  under PL).

The regression is in-sample on 144 patients with five predictors, so R² is biased
upward by roughly 5/143 under no association. The check was not part of the
registered follow-up and is descriptive.

Enrichment discriminates poorly between factors. The Hallmark epithelial–mesenchymal
transition set has 159 genes in the 2,218-gene universe and is the leading set for
26 of 35 methylation factors and 7 of 10 expression factors, so it mainly indicates
correlation with stromal genes.

## Text for the documents

All results are at rank K = 5 on the 144 TCGA patients with matched expression,
methylation and survival (75 deaths). The expression comparator uses the 2,218 genes
selected by ranking each gene's mean and variance within each training cohort, keeping
the 3,000 with the smallest rank sum, and intersecting the two cohorts' sets.

**Expression comparator.** The eight expression fits that stopped at 300 sweeps were
restarted from the same initialization with a 1,000-sweep limit. All converged, in 305
to 438 sweeps, and reproduced their first 300 ELBO values exactly. Both expression full
fits and all ten expression folds now converge. Mean CV concordance for expression is
0.582 under point-Laplace (0.591 before) and 0.620 under point-exponential (0.622 before).

**Prediction.** Point-Laplace methylation mean CV concordance was 0.545, 0.504, 0.511,
0.502 and 0.520 at 10k, 25k, 50k, 100k and all 182,934 eligible CpGs. The paired
differences from expression were −0.037 (SE 0.040), −0.079 (0.038), −0.071 (0.045),
−0.080 (0.036) and −0.062 (0.033). Point-exponential methylation reached 0.539 at
100k and 0.512 at all CpGs, with differences −0.081 (0.053) and −0.109 (0.051). The
point-exponential methylation full fits at 10k to 50k did not converge within 300
sweeps and were not cross-validated, and one 100k point-exponential fold did not
converge. Methylation is lower than expression in all seven comparisons.

**Survival-active programs.** No methylation program is survival-active in the full
fit and in at least three folds while correlating below 0.5 with every expression
program.

**Methylation factors compared with expression factors.** The two sets of patient
scores share two strong directions: the first two canonical correlations are 0.90 to
0.92 and 0.85 to 0.88 across the seven comparisons, and the fifth is 0.00 to 0.17.
Under a descriptive rule that does not require survival activity (correlation below
0.5 with every expression factor, matched in at least three folds and recurring at
another CpG count), 19 of the 35 factor fits are methylation-specific, nine under
point-Laplace and ten under point-exponential. These are recurring factors seen at
several counts, not 19 separate programs. Regressing each factor on all five
expression factors gives R² below 0.25 for 12 of the 19 and above 0.5 for three, so
those three are mostly combinations of expression factors. The clearest
methylation-specific factor is point-Laplace factor 2 at 100k and all CpGs (R² 0.06
and 0.04). Its strongest gene-set associations are allograft rejection and TNF-α
signaling, which suggests an immune component. The most purity-related methylation
factor (point-Laplace factor 1) correlates up to 0.87 with ABSOLUTE purity (143 of 144
patients) and overlaps expression factor 1.

**CpG–gene associations.** Of 9,086 pairs of a promoter CpG and one of the 2,218
genes, 76.5% have negative Spearman correlation between methylation and expression,
3,001 are below −0.3 and 949 below −0.5. These are not adjusted for purity. The top
CpGs of a methylation factor rarely sit at the top genes of its closest expression
factor (0 to 4 of 100). Four comparisons have unadjusted permutation p < 0.05, and
none survives a Bonferroni correction over 35.

**Caveats.** The results are descriptive at n = 144. CpG screening and imputation
preceded the fold split. Gene-set analysis is limited to the 2,218 genes, where the
epithelial–mesenchymal transition set (159 genes) leads for most factors in both data
types. 18 of 2,475 enrichment rows are undefined. No ICGC scoring was run, and no
rank, setup, prior or CpG count is adopted. Codex recommended 10k CpGs under
point-Laplace as a first count for a joint fit because it predicts best and is
cheapest, but the clearly methylation-specific factors appear only at 100k and all
CpGs. The choice of count is open.

**Sources.** `results/multimodal_real/outputs/methylation_single_modality_2026-11-12/K05/`:
`expression_iteration1000.csv`, `cv_summary_iteration1000.csv`,
`programs_iteration1000.csv`, and under `factor_characterization/`:
`score_correlations.csv`, `canonical_correlations.csv`, `factor_stability.csv`,
`joint_expression_r2.csv`, `top_cpgs_annotated.csv`, `annotation_enrichment.csv`,
`fgsea_by_factor.csv`, `promoter_cpg_gene_correlations.csv`,
`cross_modality_overlap.csv`, `purity_correlations.csv`, `analysis_warnings.csv`.
Patient-level inputs, scores and fits remain in ignored private storage.
