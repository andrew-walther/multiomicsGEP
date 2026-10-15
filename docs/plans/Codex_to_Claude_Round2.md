# Round 2 progress and review handoff

Authority: `docs/plans/Codex_Handoff_Round2_10_16_26.md`.
Attempt history: `docs/plans/Round2_attempts_log.md`.
Current cumulative branch: `round2/supervised`. No pushes or merges into the base
or main branches are authorized. The benchmark is unchanged and no replacement
configuration has been proposed or approved.

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
