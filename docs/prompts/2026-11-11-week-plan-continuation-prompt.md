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
locations. Andrew subsequently approved the first physical archive batch:
110 early figures and 157 tables now live under
`results/history/early-benchmarks/` in `872b92b0`. The per-file preservation manifest and
verification are in `docs/findings/2026-11-11-early-benchmark-migration-findings.md`.
Current November results, the frozen benchmark and private checkpoints retain
their paths and contents.

Later on November 11, after this prompt was first written, a bios-dissertation
session made these changes. Read them before planning A to F:

- Week plan, task F item 2 (`2a1c3a91`, `c0b45982`, `3b250851`): a fixed
  comparison arm added by Andrew, run whatever prior task A chooses. Input is a
  within-patient rank transform followed by division of each gene by its
  within-cohort SD of the ranks, without centering ("rank transform, per-gene SD
  scaling"), with point-exponential L and F, `center_projection = TRUE` and
  `strata_id = cohort`. It is compared with the current z-score point-Laplace
  model by five-fold training CV only. The plan records why per-gene SD scaling
  was chosen over raw ranks or ranks divided by p (the fitter's fixed floors), and
  the extra five fits if A keeps point-Laplace. The arm follows DeSurv's approach
  without copying it.
- Prelim Chapter 4 (`paper/prelim/project3-ssbmf.qmd`; `581487f8`, `52392829`,
  `0211e846`, `c51faa15`): the Introduction is rewritten as Background & Motivation
  and Aims & Contributions, with no results, dates or Chapter 1 references. Intro
  citations were checked against their sources. `Young2026` now cites Amber Young's
  2026 UNC dissertation (doi:10.17615/jcm5-0a36; DeSurv is its Chapter 2), not
  the paper under review (DECISIONS in bios-dissertation, 2026-11-11). The gene
  filter is described as adapted from DeSurv, because our 2,218 genes are not
  DeSurv's 1,970. For W3, after the body revamp, recheck the Introduction's current
  and planned contributions and its outline against the new sections, and remove the
  history wording listed in `paper/prelim/notes/trimmed-material-2026-11.md`.
  Keep the bib entries in `paper/ssbmf-refs.bib` identical to the bios-dissertation
  master bib, or the sync's bib check fails.

Andrew will discuss how and where tasks A to F should be done before starting
them, including whether to use Codex for any of them. Begin by confirming the current Git/checklist state and giving a concise
summary of those choices. Present the next proposed item's checklist line,
implementation plan and done condition, then wait for Andrew's go-ahead. Do not
interpret this handoff as authorization to run A to F or a results migration.

For results work, use the tested pilot as a starting point. The full scientific
workflow, real-data result families, richer figures, review/promotion contract,
physical migration and private checkpoint storage are separate scoped tasks.
The first physical batch is complete; do not repeat it. Future physical batches
need their own approval and dependency/preservation checks. Preserve current
report inputs and the frozen benchmark.

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
