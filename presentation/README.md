# presentation/

Talks for multiomicsGEP (dissertation Project 3, SSBMF).

## Contents

| Folder or file | What it is | Status |
|---|---|---|
| `walther_lab_meeting_04_09_2026/` | Rashid lab meeting, 2026-04-09 (V1–V3 drafts and a companion handout in `Notes/`) | Past talk; reference only |
| `walther_lab_meeting_06_18_2026/` | Rashid lab meeting, 2026-06-18 (deck, PDFs, revision notes) | Past talk; reference only |
| `walther_lab_meeting_08_27_2026/` | Rashid lab meeting, 2026-08-27 (deck, figure scripts in `figs/`) | Most recent talk; the best starting point for content, but results have moved on since |
| `resources/` | Earlier Quarto template and assets | Superseded by the prelim theme for new decks |
| `OUTLINE.md` | Placeholder outline for the prelim Part III deck | Not started |

The lab-meeting decks use their own CSS and section dividers. They record what was shown
at each meeting, so leave them unchanged.

## Planned: Part III of the prelim deck (not started)

A ~20-minute talk on SSBMF for the prelim exam (2026-12-07). It is the third of three
parts: Part I covers SpatialCRT spillover and Part II covers SpatialCRT incidence. The
committee is four statisticians and one non-statistician. Part III may run shorter than
Parts I and II, since the project is continuing work.

**Build it here as a standalone deck first**, then copy it into
`bios-dissertation/slides/prelim/part-3/` after Andrew approves it.

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
- the plans in `docs/plans/`, including the 2026-10-22 rebuild plan and the 2026-10-26
  parallel-sessions plan;
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
