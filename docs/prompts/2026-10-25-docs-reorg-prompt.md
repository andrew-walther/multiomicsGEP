# Kickoff prompt: reorganize multiomicsGEP `docs/`

Paste the block below into a new Claude Code session started in
`/Users/ajwalther/GithubProjects/multiomicsGEP/`. It applies the convention from the
IncidenceDesign reorganization of 2026-10-25 (`../README.md` shows the result).

---

I want to reorganize `docs/` in this repo the same way I just did for IncidenceDesign in
SpatialCRT. It's hard to see what's in there. Dates are written several different ways,
and most of the folder names don't say what kind of document is inside. Start in plan
mode. Read the documents, classify each one, and show me a mapping table before
moving anything.

**The convention, already settled. Don't re-ask these:**

- Type folders sit directly under `docs/`: `plans/`, `findings/` (results, checkpoints,
  reports from work that ran), `reviews/` (independent reviewer output), `prompts/`
  (paste-in prompts, briefs, handoffs, Codex↔Claude handoffs), `meetings/`. Add a type
  only if a real group of documents fits none of these, and ask first.
- Names are `YYYY-MM-DD-topic-type.ext`, all hyphens, with the date first. The date is the
  day the document was written: take it from the old filename, or failing that from the
  first git commit of the file. Keep the topic words recognizable.
- No archive folder. Each type folder gets a `README.md` table (date, linked doc, one-line
  purpose, status: active / complete / superseded / not started), and `docs/README.md`
  maps the folders and holds an old→new rename table.
- Living documents with no single date (a TODO, a running log such as
  `Round2_attempts_log.md`) keep an undated name. Ask me where each one goes.

**What I already know about this folder:**

- Folders already there: `plans/` (with `plans/archive/`, about 36 files),
  `reviews/` (already date-prefixed), `reports/` (about 56 entries, mostly Quarto),
  `notes/` (derivations plus `notes/gemini/`), `lb_walkthroughs/`, `progress_report/`,
  `progress_book/`, and top-level `PDAC_data_audit.*` and `Makefile`.
- I think many of these are meeting notes: `Meeting_Notes_Followup_*`,
  `meeting_2026-05-29_*`, `rashid_lab_meeting_notes_*`, `rashidlab_meeting_feedback_*`,
  `Advisor_Meeting_Summary_*`. Confirm each one by reading it.
- Dates come in at least five forms: `_10_15_26`, `_4_27_26`, `_07_12_2026`, `_05_08_26`
  and `2026-10-22_` prefixes. Some files have no date at all.
- Many reports are a Quarto set: `.qmd` plus rendered `.html`/`.pdf` and `_files/`,
  `_cache/` folders (`.aux` too). A set moves together under one new stem. Check whether
  the rendered files are tracked or ignored.

**Things to check before planning any moves:**

1. `progress_book/` is a Quarto book published to GitHub Pages (`publish_gh_pages.sh`,
   `.github/`). Its chapter paths feed the site, so propose leaving it in place unless
   you find it's safe to move.
2. `docs/Makefile` lists `lb_walkthroughs/*.qmd` and `PDAC_data_audit.qmd` by path.
   Moving those means editing the Makefile in the same commit.
3. The repo has a `post-commit` hook. Read it and confirm which paths trigger it, so the
   reorganization commits don't fire it by accident.
4. Search the whole repo (code, `AGENTS.md`, `README.md`, `ROADMAP.md`, `DECISIONS.md`,
   `PROJECT_STATUS.qmd`, `paper/`, `.claude/`, `_quarto.yml`) and the `bios-dissertation`
   repo for references to every file you plan to move. Report the counts.
5. Look for files you must not edit: hashed or manifest-checked code, frozen provenance,
   generated output. Leave their stale paths alone and list them in the rename table instead.
6. Check `git status` first. Commit or set aside any uncommitted work in `docs/` before the
   moves (ask me which).

**How to run it, once I approve:**

- Commit 1 is pure renames with `git mv`. Check it with `git show -M --numstat`: every row
  should be 0/0.
- Commit 2 updates references with a script that holds the mapping. Links between docs
  become relative (`../plans/...`). References from elsewhere become `docs/...` paths.
  Dry-run it, show me the change list, then apply it, then confirm a second run finds
  nothing.
- Then write the READMEs. Read each document to get its purpose and status. Don't guess
  from the filename.
- Verify: every relative Markdown link and `docs/...` path resolves, no old filename
  survives outside the excluded files, the progress book still renders if anything it
  reads moved, and the repo's test suite still passes.
- Update `AGENTS.md`'s rule on where documents go, and add a dated entry to the status
  log. Don't push. Ask me first.
