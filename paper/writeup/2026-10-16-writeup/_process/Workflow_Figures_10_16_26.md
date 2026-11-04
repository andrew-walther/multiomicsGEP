# Real-data workflow figures: implementation and verification

Andrew approved the two-figure plan on October 16, 2026. The figures document
the evaluated historical pipelines, not a proposed new analysis. They are
Figures 1 and 2 in the Data section of `writeup.qmd`, rendered on PDF pages
8 and 9. The complete revised manuscript is 46 pages.

## Delivered changes

- `figures/make_real_data_workflows.R` draws each workflow as vector PDF and
  300-dpi PNG, using base R only.
- The single-modality diagram separates pooled TCGA/CPTAC training from
  the five external series. It shows preprocessing, fixed prediction
  parameters, evaluation, and the PACA-AU donor-overlap caveat.
- The multimodal diagram separates matched TCGA training from ICGC/PACA-AU
  scoring. Its primary50 and broader67 evaluation branches explicitly use
  the same frozen model; the latter includes the same50 plus17 other donors.
- The manuscript explains matching within patient and why the primary
  population remains primary PDAC. Table2 now gives75/64 deaths separately
  for TCGA/CPTAC and139 for pooled training, clarifies ICGC naming/nesting,
  and labels external totals as records rather than independent patients.
- The CPTAC RNA-seq correction was already present in the working tree and
  was preserved. No unrelated existing edits were reverted.
- Root/writeup READMEs document the rebuild command, and CLAUDE.md provides
  the new figure pointer. AGENTS.md remains its symlink.

## Code walkthrough and teach-back

### `stage_box()`

1. **Purpose:** draw one readable stage of either workflow.
2. **Logic:** the `shape` argument selects an input parallelogram, process
   rectangle, model hexagon, evaluation ellipse, or dashed context box.
   `polygon()`/`rect()` draw the boundary; `text()` centers the supplied
   newline-separated label. Ellipse vertices use cosine/sine on a full circle.
3. **Design:** shape and restrained fill identify the stage without relying
   on color alone. Explicit line breaks keep the printed labels predictable
   without adding a graph-layout package. Notes are visually distinct from
   processing steps and have no data-flow arrows.
4. **Inputs/outputs:** four numeric normalized coordinates, a character
   label, an optional numeric text scale and a role-specific shape; drawing side effects and an
   invisible NULL return, with no data changes.
5. **Gotchas:** a longer label can exceed its fixed box. Inspect both output
   formats after text/layout changes; do not just shrink the text indefinitely.
6. **R notes:** `graphics::rect()`, `graphics::polygon()` and `graphics::text()`
   are base plotting functions; no new package is required. `stopifnot()`
   rejects unrecognized shapes instead of silently drawing the wrong role.

### `stage_arrow()`

1. **Purpose:** show processing order and transfer of frozen parameters.
2. **Logic:** `arrows()` connects the supplied endpoints with a short
   arrowhead and a consistent line width.
3. **Design:** explicit coordinates keep the two training/validation lanes
   and the multimodal branching unambiguous.
4. **Inputs/outputs:** four numeric normalized coordinates; a graphical
   arrow and invisible NULL, without modifying fitted objects.
5. **Gotchas:** arrows must be repositioned when their stage boxes move.
   They describe the documented pipeline, not a dependency graph generated
   automatically from model code.
6. **R notes:** base `graphics::arrows()` suffices for this static figure.

### `draw_workflow()`

1. **Purpose:** render one complete historical workflow while keeping model
   estimation separate from scoring and evaluation.
2. **Logic:** validate the selected model name, initialize a unit-square
   plotting canvas, draw the appropriate labels/stages, and connect the
   vertical lanes and parameter-transfer arrow. The multimodal evaluation
   has two branches; neither returns validation outcomes to the fitter.
3. **Design:** the two figures use the same geometry so readers can compare
   their populations, preprocessing and prediction rules. The primary50
   and broader67 appear as nested analysis populations, not separate cohorts.
4. **Inputs/outputs:** character `"single"` or `"multimodal"`; drawing on the
   active graphics device and invisible NULL. No patient matrix is loaded.
5. **Gotchas:** labels intentionally describe frozen historical settings.
   They do not track future defaults automatically. The single benchmark's
   log-transform provenance and development use of external rank evidence
   remain disclosed in its caption. The multimodal fit remains unconverged.
6. **R notes:** `par()`, `plot.new()` and `plot.window()` establish a fixed
   drawing area. `stopifnot()` rejects an unsupported diagram selection.

### Export block

1. **Purpose:** provide matching assets for the existing PDF/HTML manuscript.
2. **Logic:** loop over the two supported workflows, open a PDF device,
   draw/close it, then repeat on a PNG device using the same dimensions.
3. **Design:** vector PDF preserves printed text;300-dpi PNG matches the
   existing embedded HTML-image workflow. Both are6.5 by5.5 inches.
4. **Inputs/outputs:** the fixed drawing specifications; four files under
   `paper/writeup/figures/`, with no model, cache or patient-data input.
5. **Gotchas:** run from the repository root. The first explicit Cairo PNG
   attempt exposed missing local X11 libraries; the final script uses R's
   system-default native PNG device and runs without those warnings.
6. **R notes:** `grDevices::pdf()`, `png()` and `dev.off()` are standard R
   devices. `useDingbats=FALSE` improves PDF portability. Native PNG avoids
   introducing a dependency solely for figure export.

The prose/caption block states the flow in words and defines each cohort;
the cohort-table changes make its denominators agree with the figures.
Neither block alters the historical results, likelihood, selection rules,
or evaluation populations.

## Verification

- `Rscript paper/writeup/figures/make_real_data_workflows.R`: exit0; both
  PDFs and both PNGs created without warnings on the final run. PNG headers
  verify1950x1650 pixels, corresponding to300dpi at the declared dimensions.
- RStudio Quarto HTML and PDF renders: both exit0. The HTML refresh required
  scoped access to Quarto's existing external cache after the sandboxed
  attempt failed to open its Sass database; no dependency was installed.
- HTML inspection: both figure IDs and cross-references exist; all internal
  fragment links resolve. Each embedded workflow PNG matches the current
  generated PNG byte-for-byte.
- PDF inspection: Table2 and both diagrams/captions were visually checked;
  both figures are complete and readable. The initial page-sized gap before
  the multimodal figure was removed by compacting the drawing and placing
  the explanation after that figure. Final figure pages are8 and9, with no
  intervening nearly blank page.
- Fresh `Rscript tests/run_tests.R`: **710 passed,0 failed**, exit0. These
  verify the existing analysis code; graphical correctness is established
  separately by rendering/content checks, not by mirrored layout tests.
- All10 frozen benchmark manifest entries pass SHA256 verification from
  `results/benchmark_snapshot_2026-10-14/`.
- `git diff --check` passes for this work's changed text files. A whole-tree
  check also sees unrelated generated-HTML whitespace in concurrently edited
  `docs/PDAC_data_audit.html`; that file was not changed by this figure work.

The pre-render canonical PDF is preserved, unchanged, at the ignored local
path `data/writeup_review_backups/writeup_before_workflows_2026-10-16.pdf`,
SHA256 `7b6c1d8a34235aafa415cb827bef2f8342f36ebeb20c1ce4f7d8c518c4c4c29c`.
The existing dated circulation PDF was not overwritten.

No new scientific fit, external scoring, statistical reanalysis, default
change, commit or push was performed. Other accuracy-review items remain
in `Accuracy_Review_10_16_26.md`; this bounded diagram change does not resolve
the repeated-donor inference or historical log-transform questions.

## Explanatory callouts and stage shapes

Andrew requested clarifying the formerly unboxed rank-selection and
patient-matching notes and making the data flow visually easier to follow.
The rank callout now states that five-fold training CV selected K7 and
external concordance later informed its assessment; it no longer directs
the reader to an unspecified caption caveat. The matching callout defines
one row as expression, methylation and survival for one person and states
that TCGA trains and ICGC validates, without pooling their patients.

`stage_box()` now accepts a documented `shape` argument. Cohort data use
parallelograms, preprocessing/scoring use rectangles, fitted models use
hexagons, and evaluation uses ellipses. Dashed boxes hold contextual notes
and have no processing arrows. The shape key appears in the manuscript;
restrained fills supplement geometry without requiring color to interpret
the diagram. `draw_workflow()` chooses these roles and positions the
expanded callouts; the export block and statistical pipeline are unchanged.
The updated helper teach-back above covers the polygon/ellipse calculation,
inputs/outputs, validation and layout constraints.

Final revision checks: the figure generator exits0 without warnings;
HTML/PDF render with scoped access to the existing Sass/TeX font caches.
Both figures remain complete and readable on pages8/9 of the46-page PDF.
The updated labels were verified in extracted vector-PDF text, all HTML
fragment links resolve, and both embedded PNGs match the current figure
files byte-for-byte. Fresh full suite:710 passed,0 failed, exit0. All10
benchmark manifest entries pass, and changed manuscript/documentation text
passes `git diff --check`. No scientific fit, result, default or cohort
change, new dependency, commit or push occurred.
