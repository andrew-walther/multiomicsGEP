# Round 2 progress and review handoff

Authority: `docs/plans/Codex_Handoff_Round2_10_16_26.md`.
Attempt history: `docs/plans/Round2_attempts_log.md`.
Current cumulative branch: `round2/supervised`. No pushes or merges into the base
or main branches are authorized. The benchmark is unchanged and no replacement
configuration has been proposed or approved.

Active feature-study checkpoint: the advisor report (HTML/PDF), full-CpG
variability and fixed-K expression curves through17900 shared genes are
complete. All24 expression fits converge; no new configuration is adopted.
Methylation-only and joint count curves remain next. See the latest result
section below; historical pause/preparation statements are superseded.

Latest direction (2026-10-15, after the pause): Andrew authorizes resumption
and requests explicit prioritization of the10/2 advisor deliverables. Read
`Advisor_Priorities_and_Execution_10_15_26.md` before any further long fit.
It audits what is solved/partial/open and proposes: advisor progress report,
full-CpG variability summaries, increasing nested feature sets at fixed K7,
then matched multimodal comparisons alongside bounded supervision remedies.
Andrew approved the benchmark/fallback sensitivity order in his reply; this
changes original Task5 sequencing, not adoption criteria. No new fit has run.
The limited Round1 feature-count experiment and existing skewness plots are
acknowledged; neither completes the requested full-range feature/SD study.
The earlier supervised Single external non-inferiority failure is explicit.

Historical pause checkpoint (superseded by the resumption above): Task0 integration, Task4's
40-dataset/200-row comparator rerun and Task7 housekeeping are complete.
Task1's original panel and full-state refit remedy are not adopted; further
Single K selection needs review, with approved K7/no-pruning fallback retained.
Task2 diagnosis, fixed-prior evaluation, exact failed-fold trace and the
approved profiled-curvature panel are complete, but Task2 remains unresolved.
Profiling restores nonconstant risk in all five candidate CV folds; all five
still cap, four lack credible programs, and the full multimodal fit also caps.
Single CV passes, but one simulation and full training fail convergence.
Task3's real-data DeSurv comparison is deferred by Andrew; Tasks5–6 are pending.
The next prior-consistency observer is prepared and tested, never executed.
No new scientific fit is to start before Andrew resumes work. Full suite:
631 passed,0 failed. Benchmark tag/snapshot and write-up deliverables are unchanged;
no external Round2 score has run and no replacement is proposed.

Resume priorities: review unresolved K selection and the KL/prior-scale
questions; prospectively log/commit the observational prior audit before
running it; continue the applicable Task2 remedy ladder under the original
criteria and address the undefined original-reference pair with Andrew.
Task5 feature expansion follows the Task1/2 stopping decisions; Task6 is an
optional clean robustness rerun. Task3 requires renewed authorization.
The final sections contain the profiled-panel tables and the relation to the
approved write-up. This is an interim pause checkpoint, not completion of Round2.

The chronological updates below retain the earlier execution failures; the
complete attempt history is in `Round2_attempts_log.md`.

## Commit and verification index at the pause

All work is on `round2/supervised`; none is pushed. The pause-result commit
follows5a56f66 and is identifiable by its `Task 2 Record profiled-curvature`
subject. Historical test counts belong to their checkpoints; the current
suite is631/631.

| Task | Main commits | Verification and status |
|---|---|---|
| 0 | 2131c30, eb6b290, f2c2a85, bf3e050, b9a3952 | five integrations; final559/559; complete |
| 1 | 1d5ae7a, 221b3a6, 860a8ff, 93954f1, 06a18d5, 95b8cc4 | full-state result573/573; no adoption; needs review |
| 2 | ccebdf8, 589656a, d5ccf3f, ef4a929, 51822db, 083f94b, 5a56f66, pause-result commit | diagnosis584/584; fixed panel608 before/611 after; exact trace611/611; profile prereg625/625; pause631/631; unresolved |
| 3 | f37d646 |593/593; software prepared; real-data head-to-head deferred |
| 4 | dd11923, d353322, b941d53 |599/599 before fits,608/608 at results; complete |
| 5 | none | feature expansion pending Task1/2 stopping decisions |
| 6 | none | optional clean robustness rerun pending |
| 7 | 639351b |564/564;14 moves, live references and both renders verified; complete |

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


### Auxiliary mathematical audit during Task2-3 (no real-data fit or remedy)

Canonicalization maps (L,F,beta) to (aL,F/a,a beta), factor by factor.
It transforms L/F priors but has no beta-prior argument; the shared prior
g_beta is not transformed by the caller. Preserving fitted products alone
therefore does not establish preservation of KL(q_beta||g_beta) under a
shared fixed prior. A deterministic two-program software calculation uses
scales2.581989/1.154701 and a shared N(0,1) prior: eta changes0, while the
correct Normal KL changes1.629438/1.209438 to3.627545/1.098930.
The current carried beta KL is not recomputed by canonicalization.

This is a prior/parameterization consistency question for review, distinct
from the historical per-factor Normal CV collapse; that update re-estimates
its own variance at each coordinate. It does not demonstrate that scale
handling causes the observed zero beta. In the saved selected/no-start
traces the per-sweep scales are already close to1 at prior extinction
(selected outer99 range .998530–1.000556; no-start outer54 .998964–1.003347),
although first-outer scales are much larger. No prior, canonicalization,
objective, default or candidate is changed on this evidence. Any proposed
shared-prior configuration needs this formulation checked in addition to
the numerical adoption criteria; product-invariance tests alone are
insufficient mathematical validation. No new real-data fit or score follows.


### Saved framework loading-pattern comparison (descriptive, no new fit)

All comparisons use existing fits and the same frozen feature order. For
benchmark program4 (the frozen fit's sole credible program), the largest
absolute loading correlations by block are:

| Saved fit | Expression r | Methylation r | Combined r | Credible programs | Training C | Fold-validation C |
|---|---:|---:|---:|---:|---:|---:|
| Frozen benchmark full TCGA | 1 | 1 | 1 | 1 | .637768 | NA |
| Current full-TCGA framework diagnosis | .915133 | .898948 | .895973 | 1 | .646367 | NA |
| Current framework CV fold1 | .787112 | .685382 | .752163 | 0 | NA | NA |

Each strongest match above is program4. Relative to the current full-TCGA
program4, fold1's correlations are .938139 (expression), .846790
(methylation), .933739 (combined), yet its beta is exactly0. Thus the
loading pattern remains recognizable; this is not removal of the program
by pruning (fold1 has K7 and no pruning refit). This comparison is a
description of saved fits, not a cross-validation estimate using a
full-training loading matrix, an external selection criterion, or proof of
benchmark reproducibility. No fit, score convention or model is changed.


### Task2-2 final fixed-prior evaluation: not adopted

The100-fit local evaluation exited0:80 simulations and20 TCGA CV fits,
with0 fitting errors and0 nonempty warning records. All100 design keys,
the unchanged training-cache MD5 and the local fold-manifest SHA256 were
verified. Output: results/benchmark_sim/outputs/round2/fixed_prior/.

| Model | Noise | Recovered | Null active / reference | Converged |
|---|---:|---:|---:|---:|
| Single | .25 | 17/20 | 1.00 / 1.15 | 19/20 |
| Single | 1 | 19/20 | .35 / .35 | 20/20 |
| Multimodal | .25 | 19/20 | .35 / .75 | 20/20 |
| Multimodal | 1 | 20/20 | 0 / .65 | 20/20 |

Recovery and false-positive margins pass throughout. Single's low-noise
seed20261026 hits300 iterations, so universal simulation convergence fails.
These are unchanged K7 beds; K selection has not been adopted.

| Source/model | Fixed-prior CV mean C | Framework CV mean C | Delta / paired SE | Converged candidate/reference | Frozen training C |
|---|---:|---:|---:|---:|---:|
| Round2 Single (5 pairs) | .608225 | .623985 | -.015761 / .006968 | 5/5, 5/5 | NA |
| Round2 multimodal (5 specified pairs) | NA | NA | NA / NA | 0/5, 0/5 | NA |
| Frozen multimodal benchmark, full TCGA | NA | NA | NA | NA | .637768 |
| Frozen alpha_F0 expression-only benchmark, full TCGA | NA | NA | NA | NA | .593717 |
| Frozen expression-only framework, full TCGA | NA | NA | NA | NA | .643559 |

Single fails the one-paired-SE CV rule. Multimodal's reference fold1 has
C=NA/risk SD0; no finite five-pair difference exists, so the rule cannot
pass. The candidate's own five C values are all finite, mean .599978,
but none of its five fits has a credible program, and all cap at500.
All five reference fits also cap at500; reference credible-program counts
by fold are0,1,1,0,0. Its finite fold C values2–5 are .541502,.586047,
.485437,.359223. No conditional four-fold mean is substituted.

The frozen training rows above are labelled non-CV references, not
comparisons against a refitted benchmark. Earlier full-TCGA fixed-prior
diagnosis likewise has0 credible programs and a500 cap. Thus R1 fails
necessary real-data criteria despite the multimodal simulation passes.
No fixed prior, K, default or benchmark change is adopted. No external
cohort is scored. The unchanged-fold trace continues, and Andrew approves
the separate prospective profiled-curvature experiment after exact trace
reproduction; its implementation is still private and no scientific
profiled-curvature fit has run.

Task2-2 final checkpoint verification: fresh full suite **611 passed,0 failed**.
All100 fit records are retained; clinical folds.csv and all RDS/data remain
ignored. Only aggregate/per-fit metrics and hash/package provenance are staged.


### Task2-3 result: exact reproduction; loss at outer167/inner1

The unchanged-fold trace exited0 in2142.164 seconds, with no warnings.
All13 checks in identity.csv are TRUE: every saved posterior/prior field
and the complete500-iteration history equal original R1 cell84 exactly.
No numerical fit changed when observation was enabled.

Of2500 inner sweeps, the first all-zero posterior occurs at sweep831,
outer167/inner1. The maximum absolute pseudo-z falls from1.368511,
1.262920 and1.128647 in outer166/inner3–5 to .966195 at the next sweep.
Programs3 and4 are the last two nondegenerate coefficients. All1670
sweeps from831 through2500 have zero means/SDs; there is no later revival.
Earlier revival is real: programs1/3/4/5/6 revive7/8/2/1/2 times,
respectively. The existing full-TCGA framework trace never loses all
coefficients, and its program4 is nondegenerate throughout2500 sweeps.
This supports the per-coordinate variance-boundary diagnosis, not shared
point-mass absorption, a scoring defect or removal by pruning.

The repeated fold retains K7, has no credible program, risk SD0/C=NA,
and reaches the unchanged500 cap. Its comparison_with_benchmark.csv retains
all frozen benchmark training rows. The frozen multimodal benchmark
training C remains .637768 (one credible program); this is not a CV row.

Outputs: results/multimodal_real/outputs/round2/framework_fold_trace/:
identity.csv, fold1_trace.csv, fold1_boundaries.csv,
full_training_boundaries.csv, comparison_with_benchmark.csv,
training_source.txt and warnings.txt. The fit RDS remains ignored.

Andrew's explicit approval of the separate profiled-curvature test was
conditional on this reproduction, and that prerequisite now passes.
No remedy has been fitted at this checkpoint. Source edits still wait for
the fresh full suite/result commit. The prospective test uses the existing
framework priors/initialization and changes only the working-intercept
profile flag; it will not combine the unresolved shared-prior scale question.

Task2-3 result-checkpoint verification: fresh full suite **611 passed,0 failed**.
All13 exact reproduction checks pass; only safe factor-level diagnostic CSVs
and hash/status records are committed, with the model RDS kept ignored.


### What the profiled-curvature implementation changes (scientific results pending)

The helper `code/profile_cox_intercept.R` takes a numeric patient vector or
patient-by-feature matrix and frozen nonnegative Cox weights. It subtracts
the weighted mean, preserving dimensions/names; zero total curvature returns
the original input. Base R `colSums` and `sweep` implement this without a
patient-by-patient covariance matrix or a new package dependency.

Both fitters use this centered design only for survival projection means,
loading uncertainty, beta updates and the survival contribution to F updates.
Genomic residuals remain on the original observations. This is equivalent to
the dense profiled quadratic Hc=diag(w)-ww'/sum(w), including cross-subject
covariance induced by uncertain loadings. The existing compiled scalar sweep
is reused. Canonicalization and raw prediction retain their original inputs.

The option is FALSE by default and disabled objects match the protected
base-branch definitions. Supported experimental priors are per-factor Normal;
shared priors fail explicitly rather than combining the unresolved prior-scale
question. Single additionally requires canonical fitting, no cohort/strata
effects and no pruning. Its original ELBO monitor/rule remains unchanged;
the multimodal output rule also stays fixed. Convergence can still fail.

The evaluation runner records80 new simulation fits,10 new candidate CV fits,
two full-TCGA fits and10 unchanged cached framework CV references. It checks
input/fold hashes, per-fit keys and provenance, preserves undefined C, and
reports each arm's own mean only with all5 finite folds. Error/warning records
are kept independently. Dense precision/linear-term and covariance-objective
checks test the actual model path; separate checks verify raw genomic
residuals, base-branch identity and full-training routing. No scientific
profiled fit or external score has run at this checkpoint.

Task2-4 preregistration verification: focused23/23 (14 new plus9 existing)
and fresh full suite **625 passed,0 failed**. Default/off fits match
protected-base fitter definitions. The SLURM wrapper passes bash syntax
checking. No scientific profiled-curvature fit has run before this commit.


### Additional software invariant found during Task2-4 (no model change)

A deterministic zero-pseudo-observation helper call returns a refitted
point-mass molecular prior while returning posterior moments computed under
the previous nondegenerate prior. For F with Normal prior variance1, the
returned second moment is .1851852 and stored KL .4357921, while the new
prior is delta0. For the actual point-Laplace loading family (pi=.5,rate1),
the returned second moment is .04679413 and KL .2366599, again with a
returned delta0 prior. A shared-score example similarly returns mean
.3246915, second moment .4220161 and KL .1185925 with a new delta0 prior.

These q distributions are not absolutely continuous with respect to the
returned point prior; their KL to that returned prior is infinite. The
finite stored KL instead refers to the prior used before hyperparameter
refitting. This is a concrete prior/posterior consistency gap, predating
Round2, that product-invariance/default-identity tests do not check.
It affects interpretation of the reported ELBO/pruning criterion as the
variational objective of the returned state. It is not established as the
cause of the real-data collapse or Single oscillation. The completed full-TCGA
selected/framework states do not have a point-mass/positive-second-moment
support conflict on inspection; finite prior-lag discrepancies require
further observation. No benchmark result, default, helper or running script
is changed on this finding; record it for methodological review.

Single's completed simulation seed20261026/noise .25 shows oscillation in
both fixed-prior and profiled-framework fits. Over the last50 iterations,
relative ELBO increments are negative16 and17 times, respectively; ranges
are [-.000197505,.000128715] and[-.000180143,.000116265].
Final output changes .143724 and .105693 also fail the output tolerance.
This is not simply a near-convergence cap; the original ELBO rule remains
the criterion. No extra iterations, damping or rule change is applied.


### Pause checkpoint preparation — prior audit prepared, not executed

Andrew requested a pause after the running Task2-4 panel finishes and its
results are committed. No additional scientific fit will start before that
pause. The next diagnostic is prepared only, not an executed attempt and not
a model remedy; it must receive its prospective attempt-log entry and commit
before any execution after resumption.

`round2_prior_consistency_helpers.R` computes posterior slab probabilities
and E[|F|] from each existing point-Laplace quadratic (A,B), then evaluates
E_q[log g] under the prior used for the update and the returned refitted prior.
The difference yields the corresponding F-prior KL correction; a support
conflict is retained as infinity. The observer returns the original update
object unchanged and records aggregate iteration/factor diagnostics. It uses
base R probability functions and the existing truncated-Normal moment helper;
no dependency, fitted prior or posterior is replaced.

`run_round2_prior_consistency_trace.R` is prepared to repeat only the already
completed Single simulation seed20261026/noise .25 with identical controls.
It must reproduce the entire saved fit exactly. It records every warning,
partial loading trace on error, prior-density differences and a diagnostic
F-prior-adjusted monitor. That monitor corrects only the measured F-prior
contribution; it is not a corrected full objective or an adoption criterion.
The analytic/native slab-probability tolerance is a software check only.
Six focused tests cover numerical quadrature, units, incompatible support,
a fixed-q prior-density identity, the no-spike boundary and unchanged update
objects. These tests do not establish a remedy for any scientific failure.


### Relation to the approved interim write-up

The approved `paper/writeup/writeup.pdf` remains the benchmark report, not a
report of these unadopted Round2 experiments. Its single-modality result is
alpha_F=0, starting K7 with four retained biological programs and external
mean C=.627. That evidence remains intact; it does not establish that survival
shapes the loading updates. The new supervised recovery results answer part
of that open question on the fixed simulation panel, but convergence and
rank-selection failures prevent replacing the write-up's evaluated method.
The write-up's four classified retained programs are not evidence that the
new survival-protected deletion/refitting algorithm satisfies its K criterion.

The write-up already reports multimodal convergence limitations, no proven
prediction gain from methylation, low-variance recovery failure, a diagonal
Cox approximation and KL penalties cached against earlier priors. Round2 has
made several of these limitations more precise: collapse can arise from the
shared coefficient prior, the unchanged framework can also lose all risk in
one fixed training-CV fold, the working curvature penalizes a constant-risk
direction, and tiny helper fixtures expose posterior/prior KL inconsistency.
The latter is stronger mathematical evidence about the approximate removal
score, not proof that it caused the real-data failures. The write-up already
qualifies that score; its historical results are not invalidated by finding
an unadopted supervised configuration fails. A future methods revision needs
to address these issues before making stronger Bayesian/convergence claims.

The simulation comparator panel supplies new controlled evidence for
supervised recovery, with every failure and full denominator retained. It
does not yet show external non-inferiority, biological superiority or
multimodal predictive gain. No new Round2 external cohort has been scored.
DeSurv's separate real-data head-to-head is deferred by Andrew; omission from
the signed-Gaussian panel is also his explicit decision. Feature expansion
and final program interpretation remain future work. Process-file moves in
Task7 are housekeeping; the write-up's analysis and wording are unchanged.

Checkpoint software verification (2026-10-15): fresh
`Rscript tests/run_tests.R` completed with **631 passed,0 failed**.
This includes the six prior-consistency fixtures. The prepared observational
simulation rerun has not been executed; no new remedy was tested by those fixtures.


## Task2-4 completed: profiled Cox curvature is not adopted (2026-10-15)

The preregistered panel finished with exit0:102 records,92 new fits and10
unchanged cached framework CV references, zero execution errors and zero
warning-bearing records. All original controls, criteria, fixed simulation
beds, TCGA cache and folds were retained. Output directory:
`results/benchmark_sim/outputs/round2/profiled/`.

| Model | Noise | Recovered /20 | False active | Reference false active | Converged /20 |
|---|---:|---:|---:|---:|---:|
| Single | .25 | 18 | 1.35 | 1.15 | 19 |
| Single | 1 | 20 | .40 | .35 | 20 |
| Multimodal | .25 | 16 | .80 | .75 | 20 |
| Multimodal | 1 | 19 | .40 | .65 | 20 |

Every recovery and false-positive criterion passes. Single seed20261026 at
noise .25 reaches300 iterations; its oscillation is documented above. All40
multimodal simulations converge. These are the approved recovery/false-positive
criteria, not proof that protected K selection is resolved.

| Source/model | Mean candidate CV C | Reference CV C | Paired delta | Paired SE | CV converged | Frozen training C |
|---|---:|---:|---:|---:|---:|---:|
| Round2 Single, profiled | .624776 | .623985 | .000791 | .000791 | 5/5 | — |
| Round2 multimodal, profiled | .563678 | NA | NA | NA | 0/5 | — |
| Frozen multimodal framework | — | — | — | — | — | .637768 |
| Frozen single alpha_F=0, matching TCGA comparison | — | — | — | — | — | .593717 |
| Frozen expression-only framework | — | — | — | — | — | .643559 |

Single passes its five-pair one-SE comparison. Multimodal has five finite
candidate scores, but only four finite original-reference scores: the
predeclared paired comparison remains undefined and fails eligibility. The
candidate's own mean uses all five folds; no undefined reference is replaced
or dropped. All five profiled multimodal folds now have nonconstant risk,
including formerly collapsed fold1; their credible-program counts are
0,1,0,0,0. Every multimodal candidate CV fit reaches500 iterations, as did
all five cached reference fits. All five Single candidate and reference CV
fits converge. Frozen training C values are labelled historical training
quantities, not CV estimates or external scores.

| Source/model | Training C | Credible programs | Nonconstant risk | Converged | Iterations | K_final |
|---|---:|---:|---|---|---:|---:|
| Round2 Single, profiled | .651632 | 3 | yes | no | 300 | 7 |
| Round2 multimodal, profiled | .639172 | 1 | yes | no | 500 | 7 |
| Frozen multimodal framework | .637768 | 1 | — | no | — | 7 |
| Frozen single alpha_F=0, matching TCGA comparison | .593717 | — | — | — | — | — |
| Frozen expression-only framework | .643559 | 1 | — | no | — | 7 |

Both full-training candidates pass the defined no-collapse check and fail
convergence. No full-training factor is deleted. Eight of92 new fits cap:
one Single simulation, the full Single fit, five multimodal CV fits and the
full multimodal fit; another five caps belong to the unchanged cached
multimodal references. No intermediate pruning cap is reported. Detailed
per-fit C, activity counts, risk SD, rank, caps and warnings are in
`results.csv`, `training_cv_with_benchmark.csv` and
`training_full_with_benchmark.csv`; the latter two retain frozen rows.

**Decision: do not adopt.** Profiling removes exact constant-risk collapse
in this candidate's five training-CV folds, but does not deliver stable,
credibly prognostic programs throughout CV or universal convergence. The
original paired-reference issue remains open. No external cohort is scored,
no criterion/default is changed, and the benchmark remains the current best
state. Further Task2 work is paused by Andrew, not resolved or ladder-exhausted.
The prior/posterior KL and common-beta-prior scale questions remain separate
review items, not established explanations for these results.

Final checkpoint verification: full R suite **631 passed,0 failed** (including
six new prior-consistency fixtures); prepared observer syntax parses. The
observer itself has not been run and requires prospective registration and
commit before execution. All ten benchmark SHA256 checks pass after the
panel; tag object0fea14903a75606e495267f2f35712eac10d9fc6 still peels to
a352f17ca9cff523caa80756aa2af085e62551e1. The writeup.qmd/html/pdf files
have no diff against that tag. Patient-level folds.csv and RDS/data are
ignored and excluded from the commit. No push, protected-branch edit or
merge has occurred.


## Approved-priority deliverable1: advisor-facing progress report

Created `docs/reports/Advisor_Progress_since_10_2_26.qmd` and its five-page
PDF from saved evidence only. It covers the failed earlier supervised
external non-inferiority result, recovery/false-positive/convergence tables,
limited feature CV, multimodal training failures and the frozen seven-program
role table. It distinguishes four classified Single programs from deletion
to K4, and says K7 may appropriately represent expression, methylation and
prognostic structure. Roles overlap; fitted shares are not additive PVE or
survival contributions, and frozen labels are not transferred to a new fit.

Andrew approved the revised feature-sensitivity order before implementation.
Report renders with native RStudio Quarto and external Pandoc. Initial
LuaLaTeX font-cache failure is resolved by a temporary writable TEXMF cache;
no package/default was changed. All five pages are visually checked; the
crowded CI table was corrected. Full suite631 passed,0 failed. Benchmark
and original write-up remain unchanged; no fit or external score ran.
Next deliverable: full-TCGA-CpG descriptive inventory and variability plots.


## Task5-QC1: full-probe TCGA methylation variability (prospective)

Andrew approved the revised order before implementation. This descriptive
attempt addresses the10/2 per-probe mean/SD histogram request, before any
feature-count model fits. It changes no model, filter, preprocessing,
imputation, criterion or default and uses no survival outcomes or external
scores. No full-probe mean/SD analysis has run at this preregistration.

Inputs: raw TCGA `tcga_meth_filtered_sex.rds`$meth (305492 probes x150 samples),
existing TCGA-only training cache (144 patients x10000 screened probes), and
the frozen joint fit's training_spec feature/patient manifest as the screened reference.
Patient matching is by cache IDs, preserving its144-patient order. The raw
universe is all305492 available TCGA probes, not a newly approved shared
TCGA/ICGC universe. No ICGC values, outcomes or scores are used. Input MD5
and a non-identifying patient-ID hash are retained; raw/cached inputs and the
snapshot must not change. Cross-reactive/SNP-filter provenance is currently
unverified; the run applies no additional QC.

For each probe on observed raw beta values: count observed/missing, mean,
sample variance/SD using n_observed-1, arcsine-scale SD after asin(2*beta-1),
and proportions <.1, >.9 and outside[.1,.9]. Probes with n<2 have undefined
SD; n=0 also has undefined mean/filter membership. They remain in inventory
and aggregate denominators. No missing value is filled. Constants stay zero.
The unchanged cached screened imputation is summarized separately, without
rerunning imputation.

Descriptive groups: all available; mean filter .1<=mean<=.9; proportion
filter retaining probes with neither extreme fraction>=.95; both; mean-only;
proportion-only; current screened set. The two extreme fractions are tested
separately, not combined. Unknown membership is reported. Screened probe IDs
must match the frozen benchmark exactly as sets. There is no new SD cutoff.

Figures: per-probe SD histograms on beta/arcsine scales,80 bins, separate axes
for each group/scale; mean-beta versus SD binned density. Histograms display
all/mean/proportion/screened groups; finite/missing SD counts remain explicit.
They do not estimate a survival benefit or choose features from outcomes.
Probe inventory RDS stays ignored/local. Commit only aggregate tables,
provenance and plots. Fresh output only:
`results/multimodal_real/outputs/feature_sensitivity/cpg_qc/`.

Software validity: known means/sample variances and transform moments match
independent calculations; block size does not change values; missing/constant
and boundary cases are correct; invalid range/identities fail; filter counts
and histogram bins conserve finite observations; the driver preserves input
hashes and refuses overwrite. All warnings and errors are retained. A failure
is reported, not silently retried or used to change a cutoff. This attempt
has no model-adoption pass rule; it cannot satisfy supervision/pruning goals.

Runner: `Rscript results/multimodal_real/summarize_cpg_variability.R run`.
Serial local run, minutes, several-GB input peak;10k-row summaries bound the
extra working memory. No jobs or parallel workers submitted. Preregister,
run the full suite, record the actual count and commit BEFORE execution.


### CpG software teach-back

`align_cpg_training()` aligns a probe-by-patient raw matrix to the existing
patient-by-probe cache by patient names and retains every raw probe. It fails
on missing/duplicate/incompatible IDs. `cpg_variability_table()` takes raw
beta values and returns one row per probe; two-pass centered sums give sample
variance, while observed counts preserve missing-data denominators. Transforming
individual observations gives arcsine SD. Blocks control memory only.
`summarize_cpg_filters()` counts the predeclared groups, missing SDs, constants
and screened overlap. `cpg_variability_plots()` returns two ggplots; finite
observations enter bins, and aggregate tables retain excluded-undefined counts.
Separate axes reflect different units and group sizes; no cutoff is fitted.
`run_cpg_variability()` reads only TCGA molecular inputs and frozen probe IDs,
verifies the reference set, saves a local inventory plus safe aggregates and
figures, and checks inputs stay unchanged. Existing base-R row arithmetic,
yaml config and ggplot2 avoid new dependencies. Cached post-imputation
reference statistics are labelled separately; full raw data are not imputed.
These functions describe the data and do not fit survival or change the model.

Software preflight clarification before any descriptive run: the frozen
loading matrix has no row names. Its training_spec stores feature and patient
identities; the driver verifies those fields, not row positions. The
end-to-end fixture was tightened to this real archive schema (RED), then the
driver was corrected (GREEN). No result, test bed or statistical rule changes.

Task5-QC1 preregistration verification: focused11/11 pass after real-archive
manifest RED/GREEN check. Fresh full suite **642 passed,0 failed**. Frozen
feature IDs and patient order match the unchanged cache. No full-probe
variability result has been computed before this commit.


## Task5-QC1 result: full-probe methylation variability completed

Preregistered in dc69f6d; serial run completed with exit0, no errors and zero
warning records. All305492 available TCGA probes are described on the unchanged
144 training patients. Frozen screened probe IDs and training patient order
match the cache. No imputation, survival fit, filter/default adoption or
external scoring occurred. Outputs:
`results/multimodal_real/outputs/feature_sensitivity/cpg_qc/`.

| Descriptive group | Probes | Overlap with frozen screened10k |
|---|---:|---:|
| All available TCGA |305492|10000|
| Mean beta in[.1,.9] |182996|9990|
| Neither observed extreme fraction>=.95 |221222|10000|
| Both filters |182996|9990|
| Mean-only |0|0|
| Proportion-only |38226|10|
| Frozen/current screened reference |10000|10000|

Every probe has defined sample SD; no probe is constant. There are29414
missing raw entries, retained as missing in observed-value summaries. Of
these,2130 occur in the screened reference; the unchanged cached imputation
contains none. Median beta SD is .060008 across all probes, .179288 in raw
screened probes and .179186 in the imputed screened cache. Median arcsine SD
is .153757/.417231/.416941, respectively. Median probe mean beta is
.461005/.350031/.350125. These are data descriptions, not survival results.
The current10k screen is strongly enriched for variability; that does not
establish the best feature universe. The mean filter would discard10 of the
currently selected probes and is not adopted. Upstream cross-reactive/SNP
filtering remains unverified; this is an explicit provenance gap.

The requested SD histograms and mean-versus-SD binned plot were visually
checked. No SD threshold was selected. An independent check on21 probes
(including the most-missing probe) agrees with base-R mean/SD to maxima
4.44e-16/4.16e-17, and transformed SD to5.55e-17. Cached values at observed
raw entries agree exactly. Input hashes remain unchanged; zero warnings.
The full probe inventory RDS remains local/ignored; only aggregate tables,
plots and non-identifying hashes/provenance are committed.

The advisor progress PDF is updated to six pages with these results and the
SD figure; changed pages were visually checked. The original write-up is
unchanged. Fresh post-run full suite **642 passed,0 failed**. No new
performance, biological or supervision claim is made. Next: the fixed-K7
expression feature-count study and its training-only factor-agreement curve.


## Task5 expression preflight and agreement software (no feature fit)

The D4 training universe has17900 shared genes on273 patients. Its frozen
baseline gene order is reproduced exactly, and its standardized matrix is
exactly the corresponding subset of the all-gene matrix. Existing per-cohort
combined ranking reaches shared counts2064/5001/10000/17900 at screening
sizes3000/6607/11602/17900. No ranking rule or tie-breaking changes. The nominal
5000 target is crossed at5001 and will be plotted at its actual count.
Safe design tables are in
`results/benchmark_sim/outputs/feature_expansion/design_preflight/`.

The frozen D4 activity matrix has rank4 of7 and exact zero columns1,2,4.
Thus original Task5's inverse(L-transpose-L) continuation is undefined. A
cold-start primary curve at all four counts is proposed, preserving fixed K7,
alpha_F=0 benchmark settings, five training folds and all original criteria;
warm continuation would wait for a rank-deficient policy. This decision is
pending Andrew; no feature-count fit has started.

`code/feature_sensitivity.R` supplies only independent preparation/agreement
utilities. `feature_count_manifest()` consumes cohort-specific log-scale
training matrices and the existing rankings, reports the reached shared counts
and unavailable targets, and returns nested gene sets. It changes neither
preprocessing nor model fitting. `feature_cor()` calculates centered Pearson
correlations, preserving NA for constant columns. `feature_program_agreement()`
uses the existing greedy non-reused activity matching pattern, aligns loading
rows on common gene identities, compares proper normalized projection/risk
scores and reports beta times projection SD. This avoids raw coefficient-unit
comparisons. Matching matrices, unmatched slots and continuous best-match
margins stay visible; no match-quality threshold is selected.

Inputs to agreement must represent the same ordered patients. Named row-order
mismatches fail, and the future runner must retain the ordered identity manifest
for unnamed matrices. Constant factors remain unmatched, rather than forcing
correlation0 or a fabricated match. Signed correlations remain signed; frozen
biological labels are not automatically assigned. Existing base-R operations,
preprocessing and prediction avoid dependencies. Nine focused checks pass
from RED to GREEN, including scale/permutation/common-feature/constant and
patient-order cases. Full-suite count is recorded after completion below.

The matched expression-only multimodal framework itself retains K7 in the
frozen results. The pooled Single classified4 versus joint7 comparison also
changes training population and model settings; it does not isolate a
methylation effect. A retained K7 can be appropriate if its nonredundant
molecular/prognostic roles and stability justify it. Fewer factors alone are
not an adoption criterion for the real-data study.

Feature preparation checkpoint: fresh full suite **651 passed,0 failed**,
including nine new agreement/manifest checks. No feature-performance fit has
run. The original benchmark/write-up are unchanged; the initialization
decision remains pending.


## Task2-5: unchanged Single loading-prior consistency audit (prospective)

This is the bounded observational diagnostic already included in the approved
priority plan, executed while the independent feature-curve initialization
decision is pending. It is not a remedy, new test bed, parameter selection or
claim that debugging replaces the feature deliverables. The report/CpG tasks
are complete; no feature-performance fit is authorized by this entry.

Hypothesis: the F update returns moments under its incoming point-Laplace
prior, then refits and returns a different prior, while retaining KL against
the incoming prior. Quantify that discrepancy during one actual completed
simulation, beyond the existing synthetic counterexample. One observational
change only: a process-local wrapper records old/new priors and moments and
returns the original loading-update object unchanged.

Repeat ONLY Task2-4 cell021: seed20261026, noise .25,400 training/120 simulated
validation subjects,p50,K7; framework per-program Normal beta, no supervised
start, profile_cox_intercept TRUE, w_S1, no pruning, the original300-iteration
ELBO rule/cap. No weight, prior, initialization, bed, convergence criterion or
scorer is changed. No external or patient-level real-data input is read.
The saved original fit is
`results/benchmark_sim/outputs/round2/profiled/cells/cell_021.rds`.
Its original result is recoveryTRUE and convergenceFALSE at300; neither can
be reclassified by this diagnostic. The simulated validation outcomes are
not used to select or repair anything.

Per loading block record incoming/returned pi and rate, reconstructed slab
mass and absolute moment, native returned slab mass and second-moment sum,
old/new expected log prior, cached KL and the F-prior KL correction. Analytic
support conflicts and returned point priors with strictly positive native
second moments are counted separately; a tiny probability discrepancy alone
does not establish the native posterior support. The monitor adjustment is
only the measured F-prior contribution, not a corrected full ELBO.

Validity requires exact identity of the ENTIRE saved fit, plus every named
field and complete history, and max analytic/native slab-probability gap
<=1e-7 (existing software tolerance, not an adoption threshold). Preserve all
warnings, partial rows and errors; any mismatch blocks interpretation. Keep
all finite/nonfinite adjustments and support conflicts; no missing cases,
extra iterations, retry, damping or new prior update is allowed. Report the
original convergence failure unchanged. No adoption is possible in this audit.

Runner: `Rscript results/benchmark_sim/run_round2_prior_consistency_trace.R`.
Fresh output:
`results/benchmark_sim/outputs/round2/prior_consistency_seed26/`.
Local serial, minutes,1GB; no HPC job. Safe CSVs/provenance only may be
committed; fit RDS remains ignored. The already prepared observer is extended
only with native moment telemetry, tested RED/GREEN before this commit.
Run the full suite and record its actual count before committing this entry;
no audit fit runs until the commit exists.

Task2-5 preregistration verification: strengthened native-telemetry fixture
RED5/6 then GREEN6/6; fresh full suite **651 passed,0 failed**. No audit fit
has run before this commit.


## Task2-5 completed; cold-primary feature study authorized

The unchanged simulation audit exited0 in13.544 seconds, with no warnings.
The entire saved fit and all18 field/history checks match exactly. Across2100
loading updates, maximum reconstructed/native probability discrepancy is
6.66e-16. Three early returned point priors have strictly positive native
posterior second moments (iterations2,3,8); KL against those returned point
priors is infinite despite finite cached KL. No such support conflicts occur
in the last50 iterations. The F-prior-only diagnostic adjustment still leaves
17 negative monitor increments in the last50 iterations. It does not remove
late oscillation, correct the full ELBO, establish the cause of collapse, or
change the original convergenceFALSE/300 result. No remedy is adopted.
Aggregate evidence: `results/benchmark_sim/outputs/round2/prior_consistency_seed26/`.
Fresh full suite: **651 passed,0 failed**. No patient data are in this audit.

Andrew authorized judgment on the singular continuation issue and requested
progress. The primary fixed-K7 expression curve will use cold starts at each
of2064/5001/10000/17900 shared genes, preserving benchmark alpha_F0 settings.
Warm continuation is deferred for review; no adoption criterion changes.
Its prospective performance-panel design still precedes any fitting.
Andrew also requested HTML alongside the six-page PDF progress report; the
same source now supports self-contained HTML without changing the PDF text.

Task5-F1 prepares the authorized24-fit cold-primary expression study; see its
prospective entry for exact controls, conditional-CV limitation and agreement
analysis. The source now renders HTML and PDF with Andrew's author line only;
future presentation deliverables/code omit tool attribution at his request.
The protected benchmark write-up remains unchanged.

Task5-F1 verification before fitting: scorer RED9/12 then GREEN12/12; full
suite654 passed,0 failed. Preparation confirms all24 cell definitions, exact
frozen baseline gene order/matrix and the17900-gene maximum. Frozen D4's
historical risk=-eta orientation is retained, as previously documented;
training C0.656307. No external outcome was loaded and no model fitted.

Task5-F1 cell1 has completed (20 iterations, four retained programs, training
C0.656259 versus frozen0.656307). Two reporting-only defects were exposed:
empty-error CSV type inference and an R-version-dependent memory-column
index. Regression tests15/15 pass after RED12/15. The completed fit is
preserved; memory telemetry is unavailable for cell1, with the raw CSV kept.
Remaining23 cells await the verified correction commit.

Reporting correction verified:657 passed,0 failed; completed-cell integration
exits0 and preserves the saved fit without refitting. The summary runner
produces conditional-CV curves and the preregistered full-fit agreement tables
only once all24 cells exist. Remaining23 fits can proceed after this commit.


## Task5-F1 completed: full-range fixed-K expression sensitivity

All24 fits finish with exit0, converge under the unchanged rule and yield
finite C; zero warnings/errors. Five patient folds at each of four counts,
plus full fits, use unchanged benchmark alpha_F0/K7/no-deletion controls.
Screening/scaling precede folds: conditional CV, not nested validation.

| Genes | Training C | Conditional CV C ± SE | CV converged | Retained / active |
|---:|---:|---:|---:|---:|
| 2,064 | 0.656259 | 0.608404 ± 0.046195 | 5/5 | 4 / 2 |
| 5,001 | 0.615915 | 0.577289 ± 0.026584 | 5/5 | 4 / 2 |
| 10,000 | 0.616977 | 0.608789 ± 0.039651 | 5/5 | 4 / 2 |
| 17,900 | 0.619149 | 0.606791 ± 0.018562 | 5/5 | 3 / 2 |
| Frozen D4: 2,064 | .656307 | Not rescored | Historical full fit converged | 4 / 2 |

The full fits take20/25/31/37 iterations; maximum fold211/300. Total measured
fit time553.487 seconds; maximum R allocation estimate819.3MB (not RSS).
First-cell memory is unavailable, with its original invalid telemetry saved.
All-gene input retains240 constant columns under unchanged preprocessing.
No external cohort is scored and no feature count/default is adopted.

Full-training matches to the fresh2064 baseline:4/4/3 defined pairs at
5001/10000/17900; median common-loading correlations .297/.333/.566,
activity .301/.306/.408, projection .485/.467/.534, risk .655/.638/.505.
Adjacent5k-to10k risk correlation .972,10k-toall .667. Preserve all unmatched
and ambiguous comparisons. Greedy matches do not distinguish rotation from
biological change. Loading summaries show dense posterior means in every
nonzero full-fit program; no gene-level threshold or selection claim is made.

Fresh baseline loading/activity maxima differ from archive by about1e-10;
its oriented risk correlation is .9999999587. Small coefficient differences
and training-C delta .000048 remain explicit. The archive's beta second
moments equal squared means, consistent with the documented historical
variance-discarding bug; current uncertainty is retained. Historical D4 uses
fixed risk=-eta. Neither the archive nor its benchmark is overwritten.
The older limited feature study used point-Laplace L/F and beta; its CV means
are not directly comparable to this benchmark point-exponential/Normal curve.

Andrew clarified that avoiding arbitrary screening may favor full-feature
fitting without a C improvement, provided performance/stability are preserved.
Retain this rationale; these five-fold results do not establish noninferiority
or reliable gene-level signal/noise separation. No threshold/gate is changed.

Outputs: `results/benchmark_sim/outputs/feature_expansion/cold_primary/`:
`results_with_benchmark.csv`, `summary_with_benchmark.csv`,
`loading_shrinkage_with_benchmark.csv`, `program_agreement.csv`,
`agreement_curve.csv`, `risk_agreement.csv`, two PNG curves,
`baseline_reproduction.csv`, `independent_validation.txt`, per-cell rows and
program classifications. Patient/fitted RDS stay ignored. Mean/SE arithmetic
independently matches all five fold rows to1e-14. Both plots visually checked;
progress report HTML/PDF updated. Full suite **657 passed,0 failed**.

Next: methylation-only feature curve, then joint one-modality-at-a-time curves.
Read-only probe-ID preflight confirms original shared universe305352 versus
305492 raw TCGA probes (140 extra TCGA-only). No ICGC outcome or score was
used. Reconstruct the original TCGA imputation and verify exact frozen10k
values before fitting; do not silently change eligibility, imputation or mean
filter. Warm continuation remains deferred. Tasks1/2 remain unresolved;
Task3 deferred; Task4/7 complete; Task6 optional/unrun. This is not final goal
completion or a replacement supervised model.

Final expression checkpoint: fresh full suite657/0 after implementation,
24-cell independent checks pass, both report formats render, all eight PDF
pages visually inspected, HTML embeds three figures, and all10 benchmark
SHA256 entries pass. Main write-up QMD/HTML/PDF have no diff versus the tag.
No external scoring, default change or benchmark replacement.

### CpG expansion preparation registered, not yet executed

Task5-P1 registers original whole-universe TCGA reconstruction before any fit.
The planned methylation-only curve is10k/25k/50k/100k/all eligible, with10k
only the reference. Guard tests4/4 pass after RED0/4. The preflight also
retains the tested mean[.1,.9] and95%-extreme-tail descriptive masks for a
later separate filter comparison, as Andrew reiterates. No SD cutoff was
specified or invented. No filter/default is adopted. Exact original10k
values and ranking must pass before further modeling.

Task5-P1 software checkpoint: fresh full suite661 passed,0 failed; four guard
fixtures pass. No reconstruction run before its prospective commit.

Task5-P1 completed after da9bce1:305352 eligible CpGs,144 patients,24.719s,
zero warnings; original29307 missing values in8133 probes imputed, exact
frozen10k value/order/ranking and expression/patient/survival checks pass.
Within this universe the mean mask retains182934 and proportion mask221142.
The unfiltered30-fit methylation-only reference curve is prospectively
registered as Task5-M1; no new model has run before its commit. Filters
will be evaluated separately rather than restricting the reference to10k.

Task5-M1 preregistration checkpoint: full suite664/0, three new scoring tests
pass after RED0/3, all30 design cells prepared and two frozen training rows
reproduced. P1 source/cache fingerprints remain unchanged. No M1 fit yet.

Interim publication checkpoint: advisor report remains eight pages in HTML
and PDF and now records verified full-CpG preparation and the five-count plan.
Fresh full suite664/0 after the preregistration commit. Task5-M1 cell1 (full
10k) is running with its original controls; no result row exists yet. Do not
edit its script or loaded helpers. Larger counts have not run. Andrew's
latest clarification favors reducing arbitrary screening when model-level
performance/structure hold up; this is a motivation, not evidence already
established or a changed adoption gate.

### First methylation-count result and saved-role audit registration

Task5-M1 cell1 finished:10k CpGs, full144-patient training C=.6126711,
risk SD=.1172642, zero credible coefficients, retained classification4,
500-iteration cap,963.827s,1075.7MB R allocation high-water, no warnings/errors.
No convergence or supervised-model success claim; larger counts remain
unrun at this checkpoint. Frozen contextual training C is .637768 joint
and .643559 expression-only, with different molecular inputs.

Task5-Roles1 prospectively registers read-only characterization of frozen
joint versus saved unadopted profiled cell102, both unconverged. No new fit,
survival score, biological label transfer or external outcome. Five new
agreement/role fixtures pass; fresh full suite669 passed,0 failed. Execute
only after this registration commit. Scientific defaults/criteria unchanged.

Task5-Roles1 completed after0c6851e: safe CSVs in
results/multimodal_real/outputs/feature_sensitivity/saved_program_roles/.
Seven activity matches preserve indices (.832–.999), risk correlation.954,
one credible program4 in each fit. Program7 remains methylation-only, with
methylation norm share58.9% benchmark versus68.7% newer; program4 has
expression/methylation shares6.2%/4.5% versus7.6%/2.7%. Shares are normalized
rank-one norms, not additive PVE or survival fractions. Both fits cap500;
this supports distinct molecular roles, not convergence or optimal K.
No external scoring/biological-label transfer; input hashes unchanged.

M1 cell7 (full25k) is now running, unchanged controls and caps. No full
30-cell local launch: remaining cells require their actual timing/resource
assessment. Larger features remain the target; no silent restriction to10k.

Correction to the immediately preceding live status: M1 cell7 completed
with an execution error after27.981s, not convergence. F-update moment guard
rejects a Cox/projection input; current aggregate message cannot identify
which. No score/fit exists; the failed row is preserved with benchmark rows.
No further count cells launched. A bounded unchanged-settings diagnostic
trace versus deferral is pending Andrew's response. No scientific remedy or
filter/criterion/default change is inferred from this numeric failure.

Implementation review: the new agreement helper adapts saved multimodal
fits to the existing activity matcher using prefixed feature IDs and unit
normalizers, preserving raw sum(Y_m F_m) projection. It then replaces pooled
loading correlations by within-modality common-feature correlations, because
blocks have different units. Inputs are fitted posterior objects plus named
raw training matrices with exact patient order; outputs are matched-factor
tables and correlation matrices. Constant factors produce undefined entries,
not fabricated agreement. The role helper computes rank-one squared norms
within each block and coefficient uncertainty/per-projection-SD units; dense
nonzero means are not interpreted as validated selection. Base R matrix
operations suffice; yaml reads the existing configured credibility cutoff.
The saved-role runner only reads two capped fits, checks source fingerprints
and writes aggregate CSVs. No update equation or model default was edited.

Final checkpoint for this work block: fresh full suite669 passed,0 failed;
all10 snapshot SHA256 checks pass; frozen rounded molecular shares reproduce
exactly and four source-by-modality norm shares sum to one. Main write-up
QMD/HTML/PDF have no diff against the immutable tag. Progress HTML/PDF render
successfully; the10-page PDF was visually checked, including the new role and
methylation tables. No scientific script remains running. Only a diagnostic
decision is pending for the count experiment; the overall goal is incomplete.

Andrew subsequently authorized the unchanged-settings25k diagnosis (Task5-M2).
Its observer/guard fixtures pass3/3 after RED0/3, including exact F-update
output identity on/off. The prospective entry must be committed with a fresh
passing full suite before real-data execution. No remedy is authorized by
this diagnostic approval. Original failed rows and frozen benchmark survive.

Task5-M2 preregistration verification: fresh full suite672 passed,0 failed;
three diagnostic fixtures pass. No replay before the prospective commit.
The progress-book chapter now also records the completed expression curve,
full-CpG preparation, saved joint-role audit and both CpG terminal results.

Task5-M2 completed after df455e0: exact original error at outer6/inner1/k1,
176 F calls,32.218s, zero warnings. Inputs/data/source fingerprints unchanged.
EZ/VZ finite and VZ positive; Cox w/h have17 nonfinite values. The saved
finite predictor spans5734 units; global-max exponential shifting underflows
143/144 weights to zero, leaving17 zero risk-set sums (6 event rows).
This identifies Cox risk-set numerical underflow, alongside rapid predictor
expansion; a numeric fix alone would not demonstrate scientific stability.
No remedy or subsequent cell was executed. Safe aggregate tables are under
feature_sensitivity/moment_diagnosis/, with original benchmark-linked rows;
detailed state is ignored/local.

Andrew now explicitly asks to document larger-CpG limitations, include
loading correlations where fits exist, and update the larger write-up where
applicable. This authorizes a bounded manuscript update from saved verified
feature/diagnostic evidence; the benchmark tag/snapshot remain immutable.
Only10k methylation has a saved full fit, so cross-count CpG loading
correlation is currently unavailable, not zero or inferred from a failed run.
The expression common-feature correlations are already complete. Further
scientific numeric remedies require their own prospective design; none are
implied by this editorial authorization.

The larger write-up now incorporates the bounded author-requested evidence:
Evaluation feature-count design; Results expression performance/loading
correlations and partial CpG failure tables; Limitations/Continuing work
reflect the available evidence. No new biological labels or supervised-model
adoption is claimed. Source inventory addsS64–S70. Current manuscript source
and renders intentionally differ from the benchmark after Andrew's request;
the approved interim remains preserved at the unchanged tag. Both new tables
include saved benchmark context. The meeting progress report explicitly
separates unavailable CpG-count correlations from complete expression
correlations and the distinct saved-joint-fit audit.

Final verification for the diagnosis/manuscript work block: fresh full suite
672 passed,0 failed; diagnostic fixture RED0/3 then GREEN3/3. Both report
formats render successfully (manuscript42 pages, meeting progress10 pages).
Visual checks cover all manuscript pages in contact sheets plus readable
new/changed pages, including the corrected unbroken rank table. HTML anchors,
embedded figures, balanced Revisit markers and70 source-inventory IDs pass.
No assistant attribution in either presentation PDF. Snapshot10/10 hashes
pass and benchmark tag object0fea14903a75606e495267f2f35712eac10d9fc6 is unchanged.
Current manuscript intentionally includes the author-requested sensitivity
revision; approved interim remains retrievable from the tag. No numeric
remedy, new external score, default change or model adoption is claimed.

Task2-B1 saved-state likelihood-balance audit is prospectively registered;
five focused fixtures pass and fresh full suite677/0. Four saved capped TCGA
fits only, including two frozen benchmarks, no new fit/score or selected
weight. Execute after its commit. Andrew's latest clarification asks to solve
the25k arithmetic hangup; he does not accept inferring that CpG expansion is
inherently infeasible. Imputation is verified; numeric Cox risk-set underflow
and rapid predictor expansion must be distinguished. Opt-in arithmetic repair
and the same25k case are next, under a separate committed prospective design.

Task2-B1 completes before arithmetic-source edits. Across saved capped fits,
non-point-mass program median survival precision shares are6.30e-11–3.03e-7,
and absolute-numerator shares8.70e-7–1.82e-4. Frozen joint program4 expression/
methylation values are2.03e-7/9.00e-10 (precision),7.27e-5/4.40e-6 (numerator).
These are frozen independent-coordinate likelihood shares excluding priors,
not executed-sweep telemetry, a causal explanation or an adopted weight.

Andrew's clarification explicitly prioritizes resolving the25k expansion
failure. Task5-M3 adds opt-in same-likelihood log-space Cox arithmetic and a
bounded same25k case; no model/default/cap/criterion changes. Seven focused
checks pass after RED0/5; fresh full suite684 passed,0 failed. No captured-
point real-data evaluation or candidate fit has run before this prospective
commit. Remaining data/source identities are checked, and protected-base
original/default full-fit equality is tested. Numeric usability and scientific
convergence are separate outcomes; no blanket claim against large CpG sets.
