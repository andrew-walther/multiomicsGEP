# presentation/

Talks for multiomicsGEP (dissertation Project 3, SSBMF).

## Contents

| Folder or file | What it is | Status |
|---|---|---|
| `walther_lab_meeting_04_09_2026/` | Rashid lab meeting, 2026-04-09 (V1–V3 drafts and a companion handout in `Notes/`) | Past talk; reference only |
| `walther_lab_meeting_06_18_2026/` | Rashid lab meeting, 2026-06-18 (deck, PDFs, revision notes) | Past talk; reference only |
| `walther_lab_meeting_08_27_2026/` | Rashid lab meeting, 2026-08-27 (deck, figure scripts in `figs/`) | Most recent talk; the best starting point for content, but results have moved on since |
| `walther_lab_meeting_12_03_2026/` | December 3 lab meeting | Placeholder. May reuse the Part III prelim deck. |
| `resources/` | Earlier Quarto template and assets | Superseded by the prelim theme for new decks |
| `prelim-slides/` | Prelim Part III deck: `ssbmf_talk.qmd`, `OUTLINE.md`, `outline.html`, figures and build tools | Draft 1, 2026-11-03 |

The lab-meeting decks use their own CSS and section dividers. They record what was shown
at each meeting, so leave them unchanged.

## Part III of the prelim deck

A 20–25 minute talk on SSBMF for the prelim exam (2026-12-07), the third of three parts:
Part I covers SpatialCRT spillover and Part II covers SpatialCRT incidence. The committee is
four statisticians and one non-statistician. Draft 1 was built on 2026-11-03 from evidence
on file that day. Slide content freezes on 2026-11-30 for practice.

**Rebuild** (run from the multiomicsGEP root):

```bash
Rscript presentation/prelim-slides/tools/make_figures.R
python3 presentation/prelim-slides/tools/make_diagrams.py
cd presentation/prelim-slides && /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render ssbmf_talk.qmd && cd ../..
python3 presentation/prelim-slides/tools/outline_from_deck.py
python3 presentation/prelim-slides/tools/build_outline_html.py
```

`make_figures.R` copies scripted figures from the writeup, Chapter 4 and the progress book
(sources listed in the script) and draws the concordance forest plot. `make_diagrams.py`
draws the schematic diagrams and needs Google Chrome. `OUTLINE.md` keeps its hand-written
preamble; everything below its `SLIDES (generated)` marker comes from the deck.

**Sync to bios-dissertation.** Every commit that touches `prelim-slides/` runs
`prelim-slides/tools/sync_slides_to_prelim.sh` through the post-commit hook
(`paper/prelim/tools/post_commit_hook.sh`). It mirrors the deck into
`bios-dissertation/slides/prelim/part-3/` and commits it there without pushing. Edit the
deck here, not in bios-dissertation.

### Format

Match Parts I and II. The Part II deck in
`SpatialCRT/projects/IncidenceDesign/presentation/` is the working example: see its
`README.md` and `sudden_talk_2026-10-29.qmd`.

- **Rendering.**
  - Quarto revealjs at 1920×1080, rendered to one self-contained HTML file.
  - Use RStudio's Quarto
    (`/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto`).
- **Theme.**
  - Copy `unc-gillings.css`, `prelim-overrides.css` and `footer-part.html` unchanged from
    `bios-dissertation/slides/prelim/test/`.
  - Put deck-specific styles in `theme/deck.css`, with class names that won't clash with
    Parts I and II.
- **Sections and footers.**
  - Use section names in the footer (`data-footer-part`) and no divider slides inside the
    standalone deck.
  - Suggested sections: Background, Methods, Simulation Results, Application (PDAC
    cohorts), Discussion, Backup.
- **Slides.**
  - Each slide title is a claim written as a sentence.
  - Figures carry no in-figure titles.
  - Define every symbol of the model on a notation slide. Part III's notation is its own;
    Parts I and II each define theirs.
- **Speaker notes.** A verbatim read-aloud script with per-slide timing, and short "use
  when asked" notes on backup slides.
- **Review.** Run independent statistician and non-statistician reviews until neither
  reports a blocking issue.

### Content

Build from the project's current state, not from the August deck. Sources:
- `PROJECT_STATUS.qmd`, `DECISIONS.md` and `ROADMAP.md`;
- the [current project state](../docs/progress_book/current-state.qmd) and the October 30 supervision findings;
- the prelim chapter, `bios-dissertation/prelim/project-proposals/project3-ssbmf/`.

Reuse figures from the lab-meeting decks only where the results behind them are still
current.

### Adapting for the prelim

These are the assembly rules recorded in `bios-dissertation/slides/prelim/README.md`
(author, 2026-10-28):
- Drop the standalone title, outline and thank-you slides.
- Open with the `Part III · SSBMF` divider.
- Prefix the footer labels with the part (for example "Part III · Methods") and the menu
  titles with `III:`.
- Background builds on Parts I and II instead of repeating them.
- Name the paper or preprint, if any, on the last Discussion slide.
- Acknowledgements and funding (Dr. Liu, Dr. Rashid, Rashid Lab) go on the single prelim
  "Thank you" slide after Part III.
- Backup slides go in the third backup section, under `Backup · Part III · SSBMF`.
