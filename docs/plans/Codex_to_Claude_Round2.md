# Round 2 progress and review handoff

Authority: `docs/plans/Codex_Handoff_Round2_10_16_26.md`.
Attempt history: `docs/plans/Round2_attempts_log.md`.
Current cumulative branch: `round2/supervised`. No pushes or merges into the base
or main branches are authorized. The benchmark is unchanged and no replacement
configuration has been proposed or approved.

Current checkpoint: Task0 and Task7 are complete. Task1's original800-cell
panel rejects both protected models. The approved 160-cell full-state Single
rerun improves recovery but fails recovery/convergence/rank and is not adopted.
Single K selection needs review before an additional remedy; the approved
K7/no-pruning fallback remains available. Task2's ten-arm training diagnosis is complete: selected/no-start shared-prior
absorption is observed, fixed/shared-Normal do not give a usable fit, and all
five multimodal arms cap. The remedy ladder remains open. Andrew deferred
Task3's real-data head-to-head for review and omitted DeSurv from Task4;
Task3 software preparation is retained. Task4 is complete (40 datasets/200 rows, zero errors; two supervised caps);
Tasks5–6 remain pending. No benchmark replacement is
proposed. All ten benchmark manifest entries and the immutable tag were verified
unchanged at this checkpoint (2026-10-15).
The chronological updates below retain the earlier execution failures; the
complete attempt history is in `Round2_attempts_log.md`.

## Task 0 — integrate Round 1 (complete)

Benchmark manifest before work: ten files passed. Benchmark tag resolves to
`a352f17ca9cff523caa80756aa2af085e62551e1`.

- `2131c30`: merge protected-pruning branch; 553 tests passed, zero failed.
- `eb6b290`: merge simulation comparators; 554 passed, zero failed.
- `f2c2a85`: merge multimodal comparators; 556 passed, zero failed.
- `bf3e050`: merge robustness and initialize Round 2 records; 559 passed, zero failed.
- `b9a3952`: final DeSurv audit merge; 559 passed, zero failed.

Documentation conflicts preserve every dated entry. Round 1's blocked assessments,
conditional comparator denominators, collapsed multimodal fit and D's nonzero
process exit remain explicit. No new scientific fit has started.

## Task 1 — positive-SD protection and rank panel (pre-registered)

The approved SD>0 protection correction is implemented. Dead-program and
credible-program tests pass; default/off fits match the base branch bit for bit.
Full suite before the attempt: **563 passed, 0 failed**.
The corrected 800-cell panel and its adoption rule are committed before running;
see Task 1 attempt 1 in the attempt history. Per-cell errors, warnings and
intermediate pruning caps are retained. No empirical adoption result yet.

## Tasks 2–6

Pending Task 1 assessment. The handoff's corrected positive-SD protection rule
and each model's own convergence rule control Task 1. No criteria or defaults
have been changed. Task 2 selection and diagnosis will use training data only.


### Task 1 execution update

`1d5ae7a` committed the design before launch (563/0 tests). Task1-1 exited before
any cell: `grid <<-` reached a locked graphics binding. It has no scientific
result. Task1-2 fixes local runner scoping/export only, with an isolated entry-point
regression test and the exact same design/criteria in a fresh attempt2 folder.
No tuning or fitter change is involved in that retry.

Task1-2 pre-run verification: **564 passed, 0 failed**. The panel retry is
paused before launch to perform Andrew's added Task 7 at this stopping point.


## Task 7 — write-up process-record housekeeping

Performed at the Task 1 pre-run stopping point, with no scripts running.
No analysis, result, code or wording of the write-up was changed.

Moved with `git mv` into flat `paper/writeup/_process/` (14 files):
planning.qmd, planning.html, implementation_plan.md, drafting_prompt.md,
plan_review_prompt.md, plan_review_findings.md, review_prompt.md, review.md,
review_dispositions.md, final_editorial_check.md, verification.md,
verification_round2.md, verification_round3.md, source_inventory.csv.
The folder has no subfolders. Added `paper/writeup/README.md`; deliverable,
appendices and figures stay in their original locations. No file was deleted.

Updated live pointers in CLAUDE.md's Quick Reference, PROJECT_STATUS.qmd,
README.md, ROADMAP.md and the meeting-note prompts/sequencing table. Process
records retain bare links to each other, while deliverable/asset paths now use
`../` and planning's other relative source links have one extra parent level.
The Round 2 handoff now records Andrew's added Task 7.

Intentional historical references: four source lines in DECISIONS.md retain the
old locations, as explicitly requested. The 2026-10-14 “Writeup format after
author review of the initial draft” entry mentions implementation_plan.md and
review_prompt.md; the 2026-10-14 “Finalize the writeup plan and drafting handoff”
entry mentions implementation_plan.md and drafting_prompt.md. These are dated
audit entries, not live instructions. A new local-date 2026-10-15 relocation entry
lists all files and explains preservation for the Prompt 7 revision. No remaining
live old-path pointer outside those historical entries was found.

Render verification:
- Quarto 1.3.450 failed first under sandbox sysctl restrictions, then under host
  access because its bundled Sass Dart executable has the wrong CPU architecture.
- The existing native RStudio Quarto **1.10.18** succeeded for both sources using
  the requested `QUARTO_PANDOC=/opt/homebrew/bin/pandoc`. No tool was installed.
  Command: `PATH=/Applications/RStudio.app/Contents/Resources/app/quarto/bin:$PATH
  QUARTO_PANDOC=/opt/homebrew/bin/pandoc quarto render <source>`, from
  `paper/writeup/`, for `writeup.qmd` and `_process/planning.qmd` separately.
- Main HTML was byte-identical. The 39-page PDF was byte-identical after
  normalizing only CreationDate, ModDate and document ID; extracted text/layout
  also matched. Restored HTML/PDF with the requested `git checkout -- ...`.
- Planning HTML has **30 changed link destinations** and is included because
  those links changed. All relative Markdown targets checked exist. Bare process
  links still resolve in their shared folder.
- Hash checks confirm all 15 source/appendix/figure files, and both restored
  deliverable outputs, are byte-identical to their pre-housekeeping state.

Full-suite verification: **564 passed, 0 failed**. Task 7 is the standalone
commit containing this record (message starts `Task 7`). Task 1's unchanged
panel retry resumes after this housekeeping commit.


### Task 1 second execution update

Task1-2 was stopped after a reporting-only `!NULL` error affected fits without
pruning histories. All 209 checkpoints are retained (48 metric rows, 161 report
errors, zero malformed files). No scientific adoption conclusion is drawn.
The coordinator and all six workers were stopped before source edits. Task1-3
has the identical scientific design and fixes only this cap-reporting boundary;
a real-fit regression test reproduces the bug. The earlier entry's “no scientific
cell yet” statement describes its pre-run point, not this later partial run.

Task1-3 pre-run verification: **565 passed, 0 failed**, including the real-fit
reporting regression. Its preregistration and two-line reporting fix are committed
before retrying. Task 7 is complete at `639351b`.


### Interim Task 1 single-modality assessment — needs Andrew

All 320 single-modality cells of Task1-3 are complete without errors; the
multimodal portion is still running. Every protected group has 2/20 recovery,
20/20 ELBO convergence, and a passing false-positive margin. Median retained K
is 2 at noise .25 and 1 at noise1 for every K_init {5,7,10,15}. The prescribed
16/20 recovery and true-rank3 criteria fail. See the attempt folder's
`single_summary_interim.csv` for all eight groups and output-rule flags.

The final ladder option is already in use for single modality:
`yf_canonical_elbo()` calls `compute_survival_elbo()`, which subtracts
½ Σw Var_q(η). Repeating a “corrected” arm would duplicate the current objective.
The convergence-only rung does not apply because recovery and rank fail, not
convergence. The false-positive-stop rung is not triggered. No other approved
single-modality K remedy remains; changing the protection threshold, rank rule,
pruning criterion or model requires a decision. No protection or new single K
is adopted. Recommendation: retain K7/no-pruning as the working single model,
keep K selection explicitly unresolved, and have Claude/Andrew agree on any new
structural remedy. Ask Andrew whether Task5 may use that unchanged K7 fallback
while the K-selection issue awaits review. Continue the prescribed multimodal
panel and independent tasks; do not tune from this failed assessment.


### Andrew's clarification while Task1-3 was running

Andrew approved option 1: Task5 may use the unchanged K7/no-pruning single
fallback. He also asked that supervised K selection remain an active objective,
recognizing the strong earlier alpha_F=0 result. The failed protected-pruning
panel does not invalidate that benchmark or the unpruned supervised simulation
result (17/20 and 19/20 recovery). Any viable prior remedy from Task2 can motivate
rechecking supervised K selection under the unchanged criteria; this approval
does not authorize changing a threshold, test bed, default or adoption rule.
The recorded 2026-07-13 real-data K4/K5 warm-start alternatives have K_eff=2 and
mean external C=.6270, versus K7's .6267/K_eff=2. The older K_eff=4 beta-boost
entry is explicitly superseded in DECISIONS.md. Retained rank and survival-active
rank must remain distinct in subsequent reports.


### Task 2 preparation and data authorization (before diagnosis fits)

The frozen 144-patient, 3,000-gene, 10,000-CpG TCGA block has 2,130 missing
entries in the original raw screened CpGs. Re-imputing that subset would change
the data because the original KNN neighborhood used the larger probe universe.
Andrew explicitly approved loading the existing combined cache to extract the
unchanged training block. Saved only $training to the ignored
`data/round2_tcga_training.rds`; patient IDs, feature order and absence of missing
values match the frozen fit metadata. No ICGC outcome or score is used for
selection or diagnosis. The temporary object holding other cohorts was discarded.

Prepared trace tests fail on the current fitter for the expected missing
capabilities. Separate temporary copies pass 3/3: tracing leaves all fitted
quantities identical; the fixed shared Laplace prior is held exactly; and the
shared Normal update agrees with the analytic posterior and warm-up rule.
The original fitter remains untouched while Task1-3 runs. Four diagnosis-runner
helper tests pass, including checkpointing a failed input load. Diagnosis has
not run and still requires its preregistration commit. Longleaf wrappers are
prepared but no jobs have been submitted.

### Task 1 follow-up diagnosis preregistered

At Andrew's request, continue tracing the supervised single-modality pruning
loss while retaining the approved Task5 K7/no-pruning fallback. Task1-4 observes
the first prescribed simulation seed at both noise levels, using exactly the
current protected model. It saves the state/gains at each nullcheck, so loss
before deletion can be separated from loss during a refit. The observer changes
only its own R process's function binding and restores it on success/error;
the original fitter and running Task1-3 panel remain unchanged. No new
configuration or selection rule is proposed. The identity/error-restoration and
no-overwrite tests pass; full suite **567 passed, 0 failed** before committing
this design. Scientific case traces run only after this preregistration commit.

### Task 1 case-trace result and Andrew's additional remedy approval

Task1-4 completed (pre-run commit93954f1,567/0 tests). In the first prescribed
low-noise case the planted program is protected at K = 7 (|cor| .895844), a different
program is deleted, and recovery is lost during the K = 6 refit (best |cor| .390941;
terminal K2:.477163). At noise 1 it remains recovered but rank falls to1
(.958283→.715245 after first refit → .734053 terminal). Every stage converges;
no errors/warnings. The prior is nondegenerate throughout these two cases.
See `attempt4_trace/stages.csv`; no seeds were substituted.

The current single refit carries only L/F means and the shared beta prior,
resetting posterior variances, beta, feature precision and persistent F priors.
Andrew approved one new opt-in full-state refit remedy beyond the handoff's
Task1 ladder. Prepare its tests without changing the running panel; commit its
exact design before rerunning the same adoption criteria. Default fits and
protected-pruning thresholds are unchanged. This is an approved investigation,
not an adopted model or benchmark replacement.


### Completed Task1-3 original panel — neither model adopted

The process exited 0 after 93.23 minutes with 800 unique cells, zero errors and zero
warnings. Every original reference/selected/protected row is retained in
`results/benchmark_sim/outputs/k_selection/round2/attempt3/` (`results.csv`,
`summary.csv`, `adoption.csv`, individual checkpoints). No external cohort was
loaded. Final fit flags and intermediate pruning caps are separate below.

single protected arm:

| K_init | Noise | Recovered/20 | False/fit | Reference+.5 | Final converged/20 | MedianK | Intermediate caps |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 0.25 | 2 | 1.10 | 1.75 | 20 | 2 | 1 |
| 5 | 1 | 2 | 0.25 | 0.90 | 20 | 1 | 0 |
| 7 | 0.25 | 2 | 1.10 | 1.65 | 20 | 2 | 1 |
| 7 | 1 | 2 | 0.25 | 0.85 | 20 | 1 | 1 |
| 10 | 0.25 | 2 | 1.10 | 1.80 | 20 | 2 | 0 |
| 10 | 1 | 2 | 0.30 | 0.90 | 20 | 1 | 1 |
| 15 | 0.25 | 2 | 1.10 | 1.80 | 20 | 2 | 0 |
| 15 | 1 | 2 | 0.25 | 0.90 | 20 | 1 | 0 |

multimodal protected arm:

| K_init | Noise | Recovered/20 | False/fit | Reference+.5 | Final converged/20 | MedianK | Intermediate caps |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 0.25 | 19 | 0.50 | 1.75 | 20 | 3 | 0 |
| 5 | 1 | 20 | 0.15 | 1.05 | 20 | 2 | 0 |
| 7 | 0.25 | 18 | 0.55 | 1.25 | 20 | 3 | 0 |
| 7 | 1 | 20 | 0.15 | 1.15 | 20 | 3 | 0 |
| 10 | 0.25 | 17 | 0.70 | 1.10 | 20 | 3 | 1 |
| 10 | 1 | 19 | 0.30 | 0.95 | 20 | 3 | 0 |
| 15 | 0.25 | 18 | 0.60 | 1.25 | 20 | 3 | 0 |
| 15 | 1 | 20 | 0.50 | 1.25 | 19 | 3 | 0 |

Single fails recovery/rank in all groups. The objective is already variance-
corrected; the approved full-state remedy investigates the refit loss without
altering criteria. Multimodal passes every recovery and false-positive group,
but fails rank at K = 5, noise = 1 and convergence at K = 15, noise = 1. The latter is
seed 20261033 (finalK = 5,500 iterations,loading correlation .934102,output flagFALSE).
Protected MM also has one capped intermediate stage (K = 10/noise .25). The selected
unprotected MM arm has one final cap and one intermediate cap; the 10/1 reference
has no final cap and one intermediate cap. Single protected has four capped
intermediate stages despite every returned final fit converging. These are
reported limitations; no model meets every success criterion.

No new starting K, protection default or benchmark replacement is adopted.
The conditional 1-SE analysis is not eligible from Task1-3. Continue Andrew's
approved Single full-state remedy and subsequently the prescribed Task2
training-only collapse diagnosis. MM's current failures are rank/convergence,
not a failed recovery criterion; do not tune1.96, rank 3 or the output tolerance.

### Task1-5 ready to execute

The approved full-state remedy and unchanged 160-cell design are committed before
running. Fresh full suite: **573 passed, 0 failed**, including exact retained
state transfer, its numerical use in the first update, malformed-state errors,
no-overwrite guards, and bit-for-bit default identity against the base branch.
The new path returns n×K and p×K posterior moments, K beta moments, p precisions,
persistent loading priors and the shared beta prior; all fields use the same
retained indices. Existing coordinate updates and criteria are unchanged.
No posterior-state candidate result is claimed yet; execution follows this commit.


### Task1-5 completed: full-state refits improve recovery but are not adopted

All 160 new fits completed with zero execution errors or warnings. The comparison
retains all 160 original reference and 160 mean-only protected rows. Full-state
recovery is 15–19/20 instead of 2/20, and every false-positive margin passes.
However, K_init=5/noise=.25 misses the recovery threshold (15/20); five final
fits hit the 300-iteration ELBO cap; high-noise median K is 1 at every starting K.
The candidate therefore fails recovery, convergence and rank. No 1-SE run is
eligible, no configuration is adopted, and the benchmark/defaults stay unchanged.

| K_init | Noise | Recovery /20 | False active | Reference + .5 | Final ELBO converged /20 | Median K | Intermediate caps |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 0.25 | 15 | 1.15 | 1.75 | 19 | 3 | 2 |
| 5 | 1 | 18 | 0.20 | 0.90 | 20 | 1 | 0 |
| 7 | 0.25 | 17 | 1.10 | 1.65 | 19 | 3 | 4 |
| 7 | 1 | 19 | 0.15 | 0.85 | 19 | 1 | 5 |
| 10 | 0.25 | 17 | 1.15 | 1.80 | 19 | 3 | 6 |
| 10 | 1 | 19 | 0.20 | 0.90 | 19 | 1 | 8 |
| 15 | 0.25 | 16 | 1.10 | 1.80 | 20 | 3 | 0 |
| 15 | 1 | 19 | 0.30 | 0.90 | 20 | 1 | 0 |

The five final caps are seed20261026 at noise .25 with K_init5/7/10 (final K3),
and seed20261017 at noise1 with K_init7/10 (final K2). Intermediate caps total25,
reported separately from returned-fit flags. Outputs:
`results/benchmark_sim/outputs/k_selection/round2/attempt5_posterior/`;
`adoption.csv` records `adopted=FALSE`. Original reference table MD5:
`d45c207936855fe274ad14209474bf34`.

This completes Andrew's additional full-state experiment. The listed Single
ladder has no remaining applicable remedy: its survival objective is already
variance-corrected, and the failure is not convergence alone. Further changes
require a new decision; the approved K7/no-pruning fallback remains available
for Task5. Continue Task2's separate multimodal diagnosis and prescribed ladder.
No external data were scored or used for this assessment.

Task1-5 result checkpoint: fresh full suite **573 passed, 0 failed**. All ten
benchmark manifest checks pass; annotated tag and peeled commit remain unchanged.


## Task2-1 preparation: exact training-only diagnosis

The prepared opt-in trace is now implemented; no real diagnostic fit has run.
Multimodal traces retain each outer/inner sweep and recursive pruning stage,
including prior before/after, beta moments, raw pseudo-observations/precisions,
canonical scales, live programs, ELBO components and unchanged convergence
quantities. Full and centered molecular predictions and centered eta are extra
observations, not replacement convergence criteria. Single's existing trace
also observes its actual shared-prior refit and raw beta update quantities.

The fixed diagnostic arms are selected, fixed starting shared Laplace, shared
Normal with the same ten-iteration warmup, historical per-program Normal, and
selected without supervised start, for both models at K7/seed42. Inputs are the
unchanged144 TCGA training donors (3000 genes/10000 CpGs), extracted from the
existing imputed cache with Andrew's approval; no external scores select arms.
The runner retains per-arm errors/warnings and appends the three relevant
frozen TCGA benchmark training rows to its comparison. R1 fixed/shared-Normal
controls are opt-in; this diagnostic experiment adopts no remedy.

Teach-back: tracing copies values already computed by the fitter; it does not
feed them back into updates. Raw x/A and canonical scale remain separate because
beta changes units during normalization. Effective slab probability sums only
weights on positive-scale components, so a single weight on scale0 is correctly
recognized as a point mass. Large matrices are not saved in every trace record;
only per-program vectors and per-modality changes are retained. The scalar
coefficient helper accepts a Normal prior for the specified family contrast;
its conjugate posterior is checked analytically. The diagnosis runner has ten
independent cells, captures errors within each arm, and refuses to overwrite
completed checkpoints. A Longleaf array wrapper is provided without submission.

Software verification: the new Single prior test failed before implementation
and passed afterward (5/5 diagnosis tests). The molecular-change test failed
when its fields were absent; the no-intercept test caught NULL intercept handling;
both pass afterward (6/6 multimodal trace tests). One temporary test invocation
used a nonexistent summary helper after all assertions passed; rerunning with
the project's report_results helper exits0. Full-suite verification is pending
before the preregistration commit and any real-data execution.

Task2-1 software verification update (before any real-data fit): the first full
suite returned581 passed/3 failed. All three failures came from the new Single
trace trying to read g_beta in historical per-program fits, where no shared
prior exists. Restricted the additional prior records to shared-prior fits;
existing per-program traces remain unchanged. The affected supervision tests
and new diagnosis tests now pass22/22; the full suite is being rerun. This is a
software correction before the diagnostic attempt, not a scientific retry or
a change to its arms, criteria, inputs or defaults.

Task2-1 final pre-run verification: **584 passed, 0 failed**, including
base-branch/default identity checks, all6 multimodal trace tests and all5
diagnosis tests. Longleaf wrapper syntax passes. No diagnostic real-data fit
has run before this preregistration commit.


### Task2-1 interim: Single TCGA arms finished; multimodal arms pending

All five Single arms finish without errors/warnings. These144 TCGA donors and
3000 expression features differ from the pooled D4 training experiment. This
is diagnosis, not remedy selection; no external cohort has been scored.

| Source | Model / arm | Credible programs | Training C | Converged |
|---|---|---:|---:|---|
| round2 | single / fixed | 0 | 0.638645 | TRUE |
| round2 | single / framework | 3 | 0.651632 | FALSE |
| round2 | single / no_start | 1 | 0.654967 | FALSE |
| round2 | single / normal | 3 | 0.661636 | TRUE |
| round2 | single / selected | 1 | 0.654791 | TRUE |
| benchmark_snapshot | multimodal_framework / pruned | 1 | 0.637768 | FALSE |
| benchmark_snapshot | single_alpha_F0 / pruned_expression_only | NA | 0.593717 | NA |
| benchmark_snapshot | expression_only_framework / pruned_expression_only | 1 | 0.643559 | FALSE |

The selected/no-start/fixed shared-Laplace prior never becomes a point mass
in these Single traces; effective slab probability remains at least.5. Shared
Normal remains positive-scale (effective probability1). The fixed arm has no
credible program despite nonconstant risk. Framework/no-start Single fits hit
the300-iteration ELBO cap; selected/fixed/Normal converge in297/297/298 iterations.
Table: `results/multimodal_real/outputs/round2/diagnosis1/single_interim_comparison.csv`.


### Task3 software preparation while the unchanged Task2 fits run

No head-to-head scientific fit or external scoring has started. Isolated library:
ignored `data/round2_R_library/`. DeSurv1.0.1 is installed from
`rashidlab/DeSurv@cc9d4c2486b3125704ab71eea61f35e040a40080`; paper pipeline pinned
at `0a870a78673296a3802aeddfee0514e4a8b0f9ae`. The first install failed for missing
cvwrapr. Installed cvwrapr1.0 at `9188a9e0b5fb94f227bdd19a8b9e4f46b08a8723` from
its [maintainer's documented repository](https://kjytay.github.io/cvwrapr/).
Compilation then failed to find cmath; a successful compiler probe and retry
used only temporary SDKROOT and CPLUS_INCLUDE_PATH pointing to the existing
CommandLineTools SDK. No package source or system configuration changed.

A tiny nonnegative software-only fit at the reported parameters completed with
finite objective components and nonnegative W/H. The first reporting assertion
incorrectly treated the package's loss list as a scalar; corrected reporting
passes without a fitting change. This is not a simulation assessment or a basis
for choosing parameters. The adapter's three tests failed before implementation
and pass against the authors' preprocessing/projection APIs.

Adapter teach-back: `desurv_author_prepare()` accepts author-loaded genes-by-
patient cohorts, aligns their sample IDs, keeps the same eligible patients,
performs the authors' outer union with zero fill, then calls the package's
per-cohort feature filtering and within-sample rank transform. For validation,
frozen training genes missing from a platform are zero-filled before ranking.
The output is the package's nonnegative matrix and survival metadata. It keeps
the authors' duplicate-gene assignment behavior rather than the project raw
helper's first-duplicate convention. `desurv_author_risk()` takes W, beta and an
author-preprocessed validation matrix, restricts to the union of training-chosen
top genes, and computes X'W beta with fixed Cox orientation. Missing frozen
projection genes fail loudly. No new package dependency beyond the approved
method's required dependency was introduced.

The local PDAC raw loader matches the pinned paper loader byte for byte, and
all required expression/survival/subtype files are present. It already log2(x+1)
transforms the RNA cohorts. The DeSurv path must not apply another log. The
existing D4 wrapper's additional log flag appears to log these cohorts again;
this is a preprocessing concern for review, not an authorized benchmark change.
The tagged pipeline and benchmark remain intact; no new code has made the
benchmark irreproducible. Task3 execution still follows Task2's stopping point.


Task3 runner preparation now passes9 focused tests. No scientific execution yet.
`desurv_reported_controls()` preserves the paper's reported numeric parameters
and full4000/4000 seed limits. `desurv_seed_fit()` calls the unmodified author
optimizer, saves each fit, records convergence/errors/warnings separately, and
refuses overwrite; its software fixture reproduces the direct author's W/H/beta
and objective exactly. `desurv_author_merge()` keeps raw values separate from
screening/ranks so full tuning receives raw inputs and preprocesses within
folds. `desurv_bo_controls()` gives the paper's full50+100-evaluation,30-start,
5-fold specification. `desurv_bo_summary()` sources the unchanged pinned paper
helper from an ignored cache, checks MD5, and follows paper03's one-SE rank
postprocessing literally, retaining the other overall-best parameters. This
surprising rank-only override is the authors' script behavior, not a new choice.

`run_desurv_h2h()` checks the package and raw-loader versions, verifies ordered
training outcomes and sample counts against D4's pipeline, fits all100 seeds
and a consensus final fit, refits the specified Single model with seed42, then
makes the one external comparison with frozen benchmark rows and paired
bootstrap. It has prepare/finalize modes for an independent100-seed Longleaf
array, plus a separate full-BO mode. Missing seed statuses stop finalization;
failed seeds remain explicit and the successful-count consensus frequency is
exactly the paper's. A non-finite final fit stops before scoring. DeSurv W is
genes-by-program and its input X is genes-by-patient; the project Single uses
patients-by-genes. Orientation is fixed by their Cox predictors, and D4 uses
its frozen documented sign. The author package does the matrix optimization;
base R/parallel manages isolated checkpointing, and the project's existing
survival bootstrap supplies the paired intervals. Controlled matrices and fits
are ignored RDS files, never staged.

Longleaf scripts are prepared, not submitted: `round2_desurv_seeds.slurm`
(100 independent full fits,1CPU/4GB/2h each, prepare then finalize) and
`round2_desurv_full_bo.slurm` (one adaptive30CPU/48GB/48h search). Full BO needs
DiceKriging/lhs and the unchanged paper R/bo_helpers.R cache on Longleaf;
its hash is checked before submission work and in R postprocessing. The local
reported-parameter comparison does not require those BO dependencies. The
fixed design is in the attempt log; full-suite verification and a commit still
precede any head-to-head scientific execution.


## Task2-1 completed TCGA diagnosis (2026-10-15)

Commit `ccebdf8` preceded all ten fits; its full suite passed584/0. The process
completed with exit0, ten status checkpoints, no errors or warnings. Single
fits took100–102 seconds; multimodal fits took2384–2556 seconds each. All returned
K7, with no pruning refits. Fits and trace CSVs are in
`results/multimodal_real/outputs/round2/diagnosis1/`; real-data RDS files are
ignored. All new training results are shown beside frozen benchmark rows:

| Source/model | Arm | K | Active | Risk SD | Training C | Converged/iterations |
|---|---|---:|---:|---:|---:|---|
| round2/single | selected | 7 | 1 | 0.350339 | 0.654791 | TRUE/297 |
| round2/single | fixed | 7 | 0 | 0.309121 | 0.638645 | TRUE/297 |
| round2/single | normal | 7 | 3 | 0.396799 | 0.661636 | TRUE/298 |
| round2/single | framework | 7 | 3 | 0.456210 | 0.651632 | FALSE/300 |
| round2/single | no_start | 7 | 1 | 0.346748 | 0.654967 | FALSE/300 |
| round2/multimodal | selected | 7 | 0 | 0.000000 | NA | FALSE/500 |
| round2/multimodal | fixed | 7 | 0 | 0.059893 | 0.635837 | FALSE/500 |
| round2/multimodal | normal | 7 | 0 | 0.000000 | 0.629695 | FALSE/500 |
| round2/multimodal | framework | 7 | 1 | 0.255337 | 0.646367 | FALSE/500 |
| round2/multimodal | no_start | 7 | 0 | 0.000000 | NA | FALSE/500 |
| benchmark_snapshot/multimodal_framework | pruned | 7 | 1 | NA | 0.637768 | FALSE/NA |
| benchmark_snapshot/single_alpha_F0 | pruned_expression_only | NA | NA | NA | 0.593717 | NA/NA |
| benchmark_snapshot/expression_only_framework | pruned_expression_only | 7 | 1 | NA | 0.643559 | FALSE/NA |

The absorbing-prior hypothesis holds in both affected arms. In selected
multimodal, slab probability reaches exactly0 at outer99/inner1; all beta means
and SDs become exactly0 at outer99/inner2 and stay0 thereafter. Without the
supervised start, these events occur at outer54/inner3 and54/inner4. The prior
at extinction retains a nominal Laplace scale1, with zero slab weight; this is
not a numerical scale-zero diagnosis. Fixed Laplace never loses its .5 slab
weight, but has zero credible programs. Shared Normal avoids exact point-mass
absorption but ends with risk SD1.61e-7 and zero credible programs. Historical
per-factor Normal retains one credible program. None of Single's shared-prior
contrasts has a zero slab probability. No configuration is adopted.

Every multimodal arm hits the unchanged500-iteration output cap. Final raw
expression/methylation reconstruction changes are respectively:
selected .000446/.008929; fixed .001263/.015696; Normal .000456/.007998;
framework .001139/.006451; no-start .001139/.006375. Both exceed1e-4 in all
arms. Including intercepts or centering predictions does not remove molecular
movement (see `cap_diagnosis.csv`). Selected and no-start eta changes are0;
fixed/framework eta changes are .000424/.001163, Normal2.23e-10. Every arm's
last100 ELBO differences is positive: no observed ELBO oscillation supports
choosing damping at this point. This does not prove every coordinate is
monotone. No stopping criterion, iteration cap, prior default or benchmark changed.

The newly fitted framework training C .646367 differs from the frozen
framework's .637768. It is a current-code diagnostic refit, not a tagged-code
reproduction; no equality or superiority claim follows. This difference and
the framework's output cap remain explicit. No external score selects a remedy.

Next causally applicable rung: constrain the refitted shared Laplace slab
probability below at1/K (R2); R1 fixed-prior diagnostic already fails the
necessary full-training credible-program and convergence requirements, but
its required simulation/CV evaluation has not yet been performed. Those
unrun parts must not be described as completed. Continue the prescribed
ladder with a committed design and unchanged simulation/CV criteria.

Task2-1 result-checkpoint verification: fresh full suite **584 passed,0 failed**;
all ten benchmark SHA256 entries pass and the annotated/peeled tag IDs remain
unchanged. The result checkpoint adopts no configuration.

## Task3 stopping point — deferred by Andrew (2026-10-15)

Andrew explicitly deferred the separate real-data head-to-head for review.
The pinned installation, author-pipeline adapter, reported-parameter driver
and two Longleaf wrappers are retained as preparation. Nine focused software
tests pass; they are now registered in the master suite. No DeSurv real-data
fit, consensus run, tuning run or external scoring has occurred. The prepared
commands are dormant pending renewed authorization. This task is deferred,
not completed as a scientific comparison; no DeSurv result table exists.

The adapter takes raw author-format expression cohorts (genes x patients),
keeps the same eligible patients, merges by the paper's zero-filled union,
uses the package's own feature filter and within-sample ranks, and projects
with paper05's frozen union-of-top-genes X'W beta rule. The driver preserves
reported parameters at full precision, independent seed checkpoints and failed
seeds, then the authors' consensus procedure. Its full BO path receives raw
training values so fold preprocessing occurs inside CV. Inputs and model RDSs
stay ignored. Neither adapter changes our fitters or benchmark preprocessing.
The source pins, independent patient/outcome alignment checks, fixed orientation
and helper checksum are deliberate safeguards against comparing different
patients or silently modifying the authors' pipeline. Missing package/library
or data files fail loudly. Longleaf wrappers are prepared, never submitted.

Review concerns remain: the existing D4 wrapper appears to apply an additional
log after the author raw loader; the paper's BO summary changes K under its
one-SE rule while retaining the other overall-best parameters. Both are
preserved and reported, with no unilateral correction.

Task3 preparation/defer checkpoint: fresh full suite **593 passed,0 failed**,
including all nine author-adapter/driver tests. Both Longleaf scripts pass
bash syntax checks. No scientific DeSurv run was started.

## Task4 preparation — clean five-comparator rerun

Andrew approved omitting DeSurv. The original five-arm calibrated panel is
unchanged; new per-arm checkpoints prevent an EBMF/Cox failure from discarding
other methods. Rank0 means not recovered, C undefined, retained in n20;
rank1 uses ordinary Cox and marks lasso inapplicable. The runner reports
rank0/rank1/rank>=2 and fitting errors separately. Six focused intent tests
pass after reproducing the original small-rank failures. The attempt history
records the fixed design before execution. No scientific rerun has begun.

The added row/try helpers create explicit failures; the EBMF helper scores
three dependent arms without erasing molecular recovery on a Cox error.
The summary helper reports full sample counts and labelled finite-metric
counts, and the rank helper counts each EBMF dataset once. Inputs are the
unchanged simulation bed and real retained loading matrix; outputs are five
rows per dataset plus fixed-denominator tables. A valid rank0 flash object
cannot pass through ldf(), so its p-by-zero loadings are constructed directly.
Cox warnings and fitting errors remain visible. Lasso's nonzero count is
reported without treating it as Bayesian credibility. No new package is used.

Task4 final pre-run verification: fresh full suite **599 passed,0 failed**,
including all six small-rank/error-isolation intent tests. The five-arm design
is committed before the40-dataset scientific rerun.

Task4's first launch (dd11923) exited1 before any fit because the sandbox
could not open PSOCK port11504. Only design.csv was written. The failure is
preserved in outputs/supervision_comparators/round2/execution_status.md.
An unchanged host-socket retry uses fresh round2_retry1 outputs after the
retry record is committed. This is an environment retry, not tuning.

Task4 socket-retry pre-run verification: fresh **599 passed,0 failed**;
no scientific code changed from dd11923.

## Task2 R1 evaluation preparation

The next committed attempt completes the fixed shared prior's required80
simulation and20 training-CV fits. It changes only beta_prior_fix=TRUE,
uses the original model rules and explicit simulation/framework references,
and adopts nothing unless all criteria pass. The full-TCGA diagnosis already
fails credible-program and convergence requirements, so a favorable partial
table cannot establish adoption. The attempt log records the exact design.

The new runner creates the100-cell grid and fixed event-stratified folds;
its scaling helper estimates means/SDs from training patients only. The cell
function fits one specified model/arm, scores frozen orientation, and saves
metrics or the exact error independently. Summary logic keeps all20 seeds
and all five paired folds; an undefined C cannot pass by deletion. Prepare/
finalize modes support the Longleaf array, with a frozen training-cache hash
and all100 checkpoints required. Inputs are the approved TCGA-only cache or
original simulated beds; outputs are CSV records and ignored fit RDSs.
Eight focused tests pass. No default fitter change or new dependency is used.


## Task4 completed five-comparator panel (2026-10-15)

Pre-run commits dd11923 and d353322; full suite599/0 before each. The first
sandbox launch failed before any dataset because a PSOCK socket could not
open; its design and execution_status.md remain in round2/. The unchanged
host-socket retry completed with exit0:40 datasets,200 method rows, zero
fitting errors. Outputs:
`results/benchmark_sim/outputs/supervision_comparators/round2_retry1/`.
DeSurv was omitted by Andrew's explicit decision; no arm or threshold changed.

| Method | Noise | Recovered/n | Validation C (finite n) | Null active/fit | Median K | ELBO-converged |
|---|---:|---:|---:|---:|---:|---|
| ebmf_alone | 0.25 | 0/20 | NA (0) | NA | 2 | not reported |
| ebmf_cox | 0.25 | 0/20 | 0.634182 (20) | 1.15 | 2 | not reported |
| ebmf_lasso_cox | 0.25 | 0/20 | 0.638527 (19) | NA | 2 | not reported |
| pca_cox | 0.25 | 0/20 | 0.687753 (20) | 4.25 | 7 | not reported |
| supervised | 0.25 | 17/20 | 0.763449 (20) | 1.1 | 7 | 19/20 |
| ebmf_alone | 1 | 0/20 | NA (0) | NA | 1.5 | not reported |
| ebmf_cox | 1 | 0/20 | 0.538504 (18) | 0.45 | 1.5 | not reported |
| ebmf_lasso_cox | 1 | 0/20 | 0.553546 (6) | NA | 1.5 | not reported |
| pca_cox | 1 | 0/20 | 0.596066 (20) | 3.3 | 7 | not reported |
| supervised | 1 | 19/20 | 0.765290 (19) | 0.25 | 7 | 19/20 |

EBMF rank counts: low noise0 zero-factor,0 one-factor,20 at least2; high
noise2 zero-factor,8 one-factor,10 at least2. No EBMF fitting error. Every
method/noise denominator is20, including constant-risk and inapplicable
cases. At high noise lasso is applicable to10 datasets and has finite C in6;
its undefined scores are retained. EBMF-alone has no fitted survival risk,
so its C is undefined by design. The summaries label every finite-C count.
Lasso's nonzero count is not Bayesian credibility; false-active isNA.

The supervised fit caps at300 iterations for seed20261026/noise .25 and
seed20261017/noise1; both recover the planted program. At seed20261018/noise1
it has no active program and undefined constant-risk C, counted as not
recovered. Consequently this panel does not establish universal convergence
or satisfy all project criteria. It is a comparator report, not adoption.
There are30 warning-bearing lasso rows: installed glmnet warns that a future
release will change the Cox tie default from Breslow to Efron. The exact
warning is retained per row; no package/default was changed to silence it.

This completes the bounded Task4 rerun, including rank/denominator reporting.
Task2's remedy ladder remains open and is the priority. No real-data or
external-cohort selection occurred in this simulation task.

Task2-2 final pre-run verification: fresh full suite **608 passed,0 failed**,
including all nine new fixed-prior protocol tests and unchanged-default checks.
Longleaf wrapper syntax passes. No candidate simulation/CV fit ran before
this preregistration commit.

Task4 result-checkpoint verification: fresh full suite **608 passed,0 failed**.
The40-dataset retry exited0; all200 rows and40 cell checkpoints are preserved.
No model/default change follows the two supervised caps.

Task2 R1 interim: all80 simulations complete without errors. Single recovers
17/20 and19/20, with false-positive margins passing, but one low-noise fit
caps (seed20261026,300 iterations). Multimodal recovers19/20 and20/20, with
all false-positive and output-convergence criteria passing. These results do
not override the failed full-TCGA credible-program/convergence requirements.
TCGA CV is still running; no criterion, control or running script is changed.

Training CV's patient-level folds.csv is local/ignored because it contains
patient IDs and event indicators. Only aggregate/per-fit metrics are committed;
reproduction uses the fixed seed42 event-stratified rule and approved local
training cache. Transfer the fold manifest securely with controlled inputs
when reusing references on Longleaf. No running analysis changed.

CV interpretation: this evaluation conditions on the handoff's existing
screened/imputed TCGA cache. It does not redo raw-data screening or methylation
imputation separately inside each fold. Single's centering/scaling is estimated
within the training fold. Both candidate and reference use the same fixed
cache/folds. These are conditional CV comparisons, not fully nested raw-data
pipeline validation; the approved imputation/test bed is not changed.

### Fixed-bed interpretation for review (no analysis change)

The approved generator explicitly uses eta=(Y_expression F_expression+
Y_methylation F_methylation)beta (`code/simulate_multimodal_yfb.R`, cohort
hazard calculation). Its Gaussian molecular residuals therefore enter the
survival predictor too. Raising the molecular noise multiplier need not make
survival recovery harder, because it also changes prognostic variation.
Single's fixed bed takes expression only from that joint generator, while its
stored c_oracle uses the full joint eta; this is not an attainable
expression-only oracle. The original recovery correlation still compares
fitted loadings after expression z-standardization to the original raw F truth.
Those feature units differ. These are limitations/interpretation points of the
approved bed and score, not Round2 changes. Nothing is regenerated, redefined
or tuned in response. Review this framing before making claims in the write-up.

A zero single-program recovery count also does not imply no predictive value:
PCA->Cox has low-noise validation C .687753 in Task4. Its recovery criterion
selects one component by Cox z and applies the specified loading-correlation/
sign rule; a predictive span can exist without one matching axis.


Task2-2 framework CV collapse and Andrew's direction: completed multimodal
framework fold1 has beta means and second moments exactly0, no credible
program, risk SD0, C=NA, and hits500 iterations. Other completed framework
folds2/3 retain one prognostic program. Andrew explicitly requests continued
diagnosis of the constant prediction, not a tied-risk C=.5 substitution.
Keep the committed scorer and all five finite paired-fold requirements.
Compare this failed training fold with the successful full-TCGA framework
diagnosis and the other folds; do not tune on external data. The ongoing
R1 evaluation remains unchanged. No scoring amendment or new remedy is
adopted. This is an additional instability in the reference itself.

### Multimodal framework instability in a training fold

The historical framework's completed fold1 returns zero beta posterior means
and second moments; C remains undefined. Folds2/3 retain a prognostic program.
The saved-state per-factor EB update estimates variance max(0,x²−1/A), with
all seven pseudo-z magnitudes below1 in fold1. Beta KL is0 from outer167
onward. Unlike the shared point-mass prior, these variances are re-estimated
each coordinate and can revive; an unchanged-fold trace will establish the
trajectory. Andrew explicitly requests this diagnosis rather than assigning
C=.5. The reference instability is not resolved, and no configuration is
adopted. The preregistered trace must reproduce the original posterior and
history exactly; only observation changes.

The diagnostic takes the original event-stratified fold and raw imputed
training matrices, enables the existing trace, and records beta moments and
pseudo-observations after each inner sweep. Its summary distinguishes a
temporary all-zero sweep from persistent loss through the last sweep.
It uses no external cohort, changes no prior or convergence rule, and reports
the frozen benchmark rows alongside the repeated fold's metrics.

Task2-3 preregistration verification: new focused checks3/3 passed after
RED0/3; full suite **611 passed,0 failed** before the diagnostic commit.
The trace has not yet been fitted at this checkpoint.
