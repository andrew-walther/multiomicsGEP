# Continue the November week plan after housekeeping

Date: 2026-11-11. Author: Andrew Walther.

Paste the following into the next project session.

---

Work in `~/GithubProjects/multiomicsGEP` on `main`. Read AGENTS.md, README.md,
the full `docs/plans/2026-11-09-week-plan.md`, and
`docs/reviews/2026-11-11-housekeeping-and-cv-handoff.md` before proposing work.
Use the week plan's Progress checklist as the authoritative progress record.

The following work is complete. Do not repeat it:

- 00: Gemini advice review, preserved specifications and actionable short list.
  Andrew approved closure and source-folder deletion.
- 0a: README cleanup accepted by Andrew.
- 0b: Indexed script/folder/report inventory, independently reviewed.
- 0c: Approved branch cleanup. Keep main, gh-pages and the standalone
  rebuild/unified-ssbmf reference. Deleted tips remain on main or archive tags.
- G: Full results-rebuild design delivered and reviewed, with local execution,
  an expandable result catalog and model/version provenance as explicit goals.
- 0: Cohort/event folds implemented in `dbb4113`. Combined independent review of
  that change and `68e3646` approved in `681ce91c`, with 62 focused tests and
  additional center, score, missing-stratum and compatibility checks passing.

Andrew also authorized a bounded local software pilot after these items.
`cc11b34` implements it under `results/pilot/`. Read its README and
`docs/findings/2026-11-11-local-results-pilot-findings.md`. The latest full suite
passed 801 tests with zero failures. All 14 workflow checks and the independent
software review passed. Setup, build and catalog commands passed in a clean
checkout, reproducing both aggregate CSVs. The bounded-β optimizer fallback
warning remains documented. A restricted-tool-sandbox OpenMP shared-memory
failure did not reproduce when the unchanged pilot commands ran outside that
sandbox. Do not change model code to work around that execution restriction.

The pilot is synthetic software validation, not a scientific result. It has not
rebuilt the historical scientific analyses or changed a production results
pointer. The results README now separates the pilot, current reported evidence,
frozen evidence and indexed history. `results/history/README.md` indexes existing
locations in place. No historical files or private checkpoints have moved.

Andrew will discuss how and where tasks A to F should be done before starting
them. Begin by confirming the current Git/checklist state and giving a concise
summary of those choices. Present the next proposed item's checklist line,
implementation plan and done condition, then wait for Andrew's go-ahead. Do not
interpret this handoff as authorization to run A to F or a results migration.

For results work, use the tested pilot as a starting point. The full scientific
workflow, real-data result families, richer figures, review/promotion contract,
physical migration and private checkpoint storage are separate scoped tasks.
The proposed first physical batch is described in
`docs/findings/2026-11-11-results-navigation-findings.md`; it is not approved by
this prompt. Preserve current report inputs and the frozen benchmark.

Keep future jobs local where practical. Measure cost before starting a large
grid. Longleaf requires a separate decision and practical setup guidance.
Every result must identify its model/options, data/preprocessing, code version
and evaluation design. Use training data for all selection and external cohorts
only to score fixed models. Save the same paired cohort/event folds where
comparisons require them. Request `fold_cohort_id` and `strata_id` separately.

The scale-robust dead-program rule is still an open scientific follow-up. Review
the week plan's diagnostic requirements before changing it. The Gemini short
list contains candidates for consideration, not adopted methods or permission
to expand an already registered comparison.

Work directly in this checkout. Do not create a worktree or change branches by
default. Do not delete ignored fits or patient inputs. Keep the checklist current
with commit hashes and evidence, index every new document, obtain independent
review, run the relevant tests, and commit/push each completed logical item.

---
