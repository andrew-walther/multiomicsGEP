# Handoff prompts for 2026-11-06 (usage-limit window)

Three independent tasks. Run them one at a time, in this checkout on branch `fix/yfb-supervision`.
Other sessions may be working in this folder, so commit only the files you changed, never stash,
reset or switch branches, and do not merge into main. Follow AGENTS.md: no tool attribution in
code or deliverables, plain literal prose with no semicolons or dashes as punctuation, never
use external cohorts to choose anything, and never commit patient data. A reviewer will check
each task afterwards, so list every file you changed and every number you report with its
source file.

---

## Prompt 1: 11/6 project writeup, phase 1 (prep and data-driven skeleton)

Prepare the updated full project writeup in `paper/writeup/2026-11-06-writeup/`. A previous
session started it and was stopped mid-turn. Run `git log --oneline -10 -- paper/writeup/2026-11-06-writeup`
and `git status` first, and build on what is there. Do not edit the circulated
`paper/writeup/2026-10-16-writeup/`.

1. Read the 10-16 writeup and its `_process/` notes, PROJECT_STATUS.qmd, the newest DECISIONS.md
   entries, ROADMAP.md and `docs/findings/2026-11-0*`.
2. Write a section-by-section plan in the folder README. For every section, mark what carries
   over, what changes and what is new, and name the evidence file for every number.
3. Set up `writeup.qmd` so every number comes from R chunks that read saved CSVs. No hand-typed
   numbers (GitHub issue #6 item 8). `docs/progress_book/chapters/2026-11-06.qmd` shows the
   pattern.
4. The writeup must address issue #6 (`gh issue view 6`):
   - item 2: use single-log numbers;
   - item 4: state the pre-registered October 15 supervision null result plainly (ΔC = −0.003,
     −0.029 to 0.024; DECISIONS.md around line 1031);
   - item 5: lead with the rank-7 EBMF + Cox comparison, not the 40-factor lasso comparator;
   - item 6: list every use of the external cohorts (K = 7, preprocessing variant, α_F,
     N_frozen, comparator configuration) and disclose the 47 shared PACA-AU donors.
5. Draft only the sections that do not depend on the prelim Chapter 4 figures. Chapter 4 is
   being revised separately. Leave clearly marked placeholders for its figures.
6. Commit and push. Stop there. Do not render the final PDF or run the prose review yet.

---

## Prompt 2: GitHub issue #6 item 5 leftover (risk orientation in other runners)

Issue #6 item 5 flagged `max(c_raw, 1 - c_raw)`, which orients a risk score with the
validation outcome. Commit b85d390 fixed the two named scripts. About ten other runners still
use the pattern (for example `results/benchmark_sim/run_desurv_comparison.R:164`,
`run_joint_bo.R:107`, `run_cohort_lmm_benchmark.R:126`).

1. Find every remaining occurrence: `grep -rn "1 - c" results/ code/ | grep -i max`, plus variants.
2. For each one, decide whether it feeds a number reported in a current document (writeup,
   progress book, prelim, DECISIONS.md, findings).
   - If it does, replace it with the training-fixed orientation `frozen_reverse_cindex()`
     (`code/concordance_ci.R`).
   - If it does not, still replace it, unless the script is in a dated archive or legacy
     folder. Leave those unchanged and list them.
3. Write a short table (script, line, reported anywhere?, action) in
   `docs/findings/2026-11-06-orientation-audit.md` and add its README row.
4. If rescoring is possible from saved outputs without patient data, rescore, report whether
   any number changes, and save the rescore CSV. Do not refit models.
5. Run `Rscript tests/run_tests.R` if you touched anything under `code/`. Commit and push.

---

## Prompt 3: GitHub issue #6 item 6 (external-cohort uses and donor-level resampling)

1. Write `docs/findings/2026-11-06-external-cohort-uses.md`, an exact list of every decision
   that looked at external-cohort results: K = 7, the preprocessing variant, α_F, N_frozen, the
   comparator configuration and anything else. For each one give the source (DECISIONS.md line,
   findings file or commit). Name at least one external cohort that played no part in any
   choice, or say plainly that none qualifies.
2. Add donor-level resampling for the pooled external comparison in
   `results/benchmark_sim/score_yfb_supervised_external.R` (S-12). PACA_AU_array and
   PACA_AU_seq share 47 donors. The current bootstrap resamples records. Resample donors so
   that a shared donor's array and RNA-seq records move together. Find the donor identifiers in
   the loaders (`load_pdac_raw()`). This needs the local patient data (`PDAC_DATA_ROOT`), and
   if it is unavailable, stop and say so. Write only aggregate outputs:
   `external_paired_diff_donor.csv` next to the S-12 outputs. Report the new interval beside
   the record-level one.
3. Do not choose anything from these results. Add a findings README row. Commit and push.
