# Housekeeping and CV fold-option handoff

Date: 2026-11-11. Author: Andrew Walther.

The requested housekeeping, G design and item 0 are complete on main. After
this handoff was first published, Andrew authorized the [combined independent
CV review](2026-11-11-stratified-cv-independent-review.md). It approved `68e3646`
and the new fold option on November 11. All 62 focused tests and additional
independent checks passed. No analysis task A to F was started. The later
approved move of 267 early artifacts is complete, with contents preserved.
Current November results, frozen evidence and checkpoints remain in place.
No scientific default changed.

## Current continuation prompt and results navigation

Completed in `9ed706a`, with independent review and local-link verification.

Use the [November 11 continuation prompt](../prompts/2026-11-11-week-plan-continuation-prompt.md)
for the next session. It identifies completed work and leaves decisions on the
scope and location of A to F to Andrew's discussion with the next owner.
The [results index](../../results/README.md) now separates the pilot, reported
November evidence, frozen benchmark and historical studies. The
[history index](../../results/history/README.md) links existing paths in place.
The first physical batch was subsequently approved: 110 early figures and
157 tables moved to `results/history/early-benchmarks/`. The move is in `872b92b0`. The
[migration findings and manifest](../findings/2026-11-11-early-benchmark-migration-findings.md)
record preservation and path validation. Current reported results, frozen
evidence and checkpoints remain in place. Further physical batches need
separate approval.

## Work and commits

| Item | Delivered evidence | Work commits, followed by checklist commits where separate |
|---|---|---|
| 00 Gemini review | [33-idea table and actionable short list](../findings/2026-11-10-gemini-advice-findings.md), two preserved specifications in this folder. Source folder deleted after Andrew approved closure. | `3f03cf5b`, `ae79b19a`; approved closeout `1e2c73af`, `47690aa2` |
| 0a README cleanup | [README](../../README.md), [retention ledger](../findings/2026-11-10-readme-cleanup-findings.md), quantified simulation statement and simplified documentation links. Andrew accepted it. | `789d38c4`, `904ca349`; feedback `961365cf`, `71750325`; acceptance/deferred document revamp `83f13932` |
| 0b Inventory | [Findings and three indexed tables](../findings/2026-11-11-results-inventory-findings.md): 204 scripts, 251 locations and 82 report sources. Independently reviewed. | `728d0678`, `1a6188de` |
| 0c Branch review and cleanup | [Approved table and preservation record](../findings/2026-11-11-branch-review-findings.md). Removed 33 local and 14 remote branch names. Retained main, publishing branch gh-pages and standalone rebuild/unified-ssbmf. | Review `8746bab`, `a23cdbc`; cleanup `fc8cd2c8`, `0a44c84a` |
| G Design only | [Rebuild design](../plans/2026-11-11-results-rebuild-design.md), with proposed roles for all 204 scripts. Independently reviewed, delivered for Andrew's consideration. | `7db86ef`, `86ee3ead` |
| 0 Code only | [Cohort/event fold generator](../../code/train_test_split.R), four direct CV helper functions and [targeted tests](../../tests/test_cohort_folds.R). | `dbb4113`, `debe6e0`; combined review `681ce91c` |

## Added local pilot and final handoff

Andrew subsequently authorized a small implementation beyond G's design.
`9613ab6` registers its scope and `cc11b34` implements it. The
[pilot guide](../../results/pilot/README.md) gives setup/build/status/catalog
commands. The [findings](../findings/2026-11-11-local-results-pilot-findings.md)
record 14 workflow checks, 801 passing full-suite tests and successful
clean-checkout setup/build/catalog execution. Both aggregate CSVs reproduced
byte for byte under the tested environment. The independent software review
approved the pilot after correcting the CLI/environment issues it exposed.

The six original items and this added pilot are complete. This is ready to hand
off. The next owner has a tested synthetic workflow to extend, not a migrated
scientific results collection. No analysis task A to F was started. The later approved early-artifact
relocation is recorded above. Follow the week plan and obtain Andrew's approval for the next scoped
analysis or workflow extension. The full results migration and private archive
destination remain later decisions.

## Fold option and verification

Use `cohort_id = cohort` with `create_stratified_folds()`, or
`fold_cohort_id = cohort` with `select_K_cv()`, `select_alpha_cv()`,
`cv_survival_loglik()` or `bicv_genomics_loglik()`. The latter uses the option
only for row folds. This split label is separate from molecular offsets,
Cox strata and cohort-specific coefficients. Those still need their own arguments.

The generator shuffles within joint cohort/event cells and assigns them
round-robin, continuing across cohorts within each event status. It checks that
every training complement contains all supplied cohorts. NULL retains the exact
event-only default. A sparse allocation can error even when a different design
would work, and events are not guaranteed within every cohort of each held-out
fold. No alternative-allocation search was added. Inputs, outputs, design choices
and a function-by-function explanation are in the
[code guide](../../code/README.md#cohort-and-event-folds).

`Rscript tests/run_tests.R` completed with exit status 0: **794 passed, 0 failed**.
The 11 new targeted checks also passed standalone and were independently rerun.
They check exact partitions, cohort/event balance, training-cohort coverage,
reproducibility, unchanged historical assignments, invalid/sparse labels and the
actual training row IDs used by each direct helper, including both rank-CV model
paths. The full suite emitted a warning that the bounded β optimizer used its
profile-likelihood fallback, plus zero-factor diagnostic messages. No failure was
suppressed. A fresh review of the new option found no actionable defect. This
bounded check was followed by the separately authorized combined review linked above.

## Gemini ideas retained for future consideration

- Test multimodal scale invariance, modality balance and consistent prior/KL rescaling.
- Derive a coherent point-exponential fixed-q prior update if the planned prior comparison exposes convergence failures.
- Distinguish anchored fixed residual precision from bounded adaptive precision in instability studies.
- Register genomics tempering with Cox weight one only after settling the objective and KL weights.
- Consider regularized nonzero β variance and paired positive/negative Cox starts for documented collapse.
- Retain the planned prior/PIP, cohort-noise and biological checks. None of these suggestions is newly adopted or authorized for analysis here.

## Left for Andrew and the independent reviewer

Andrew can review G's proposed pilot, core release families, layout and retention
rules. Expanding the scientific pipeline, moving further runners/results or selecting an external
checkpoint archive needs a later approved task. The broader cleanup of
PROJECT_STATUS.qmd, DECISIONS.md and ROADMAP.md remains on the ROADMAP.

The combined review of `68e3646` and `dbb4113` is complete. No defects required
code changes. The scale-robust dead-program rule remains a scientific follow-up,
not a blocker for the reviewed CV integration. G now explicitly specifies local
execution, an extensible result catalog and clean-checkout reproduction. Full scientific workflow implementation and physical results migration remain
separate tasks.

The platform-normalization archive is local only because it contains unpublished
patient risk-score history. Its tag and verified, restoration-tested bundle are
recorded in the branch closeout. It has no independent off-machine backup from
this cleanup. The attached old worktree and its untracked/ignored data remain
intact at detached commit `bde88ff8`. No worktree was removed.
