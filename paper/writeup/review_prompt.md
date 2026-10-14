# Codex handoff: independent review of the initial writeup draft

## Assignment

On branch `writeup/draft` (draft commit `0a7cfc9`, monochrome styling `a9f4b47`), review the initial writeup: `paper/writeup/writeup.qmd`, `appendices/single_modality.qmd`, `appendices/multimodal.qmd`, and the rendered `writeup.html` and `writeup.pdf`. This is Task 4 of `paper/writeup/implementation_plan.md`.

Run `git status -sb` and `git branch --show-current` first; the branch must be `writeup/draft`. If it differs, stop and ask. Work in this checkout only (no worktree or second clone).

**Write findings to `paper/writeup/review.md`. Do not edit the draft, appendices, figures, source inventory, code, or saved results.** Rank findings by severity. Give each one a location (section, table/figure number, or appendix equation), the evidence, and a concrete suggested fix: replacement wording, a proposed structure, or the exact correction. Claude Code applies the confirmed findings afterward (Prompt 4), and Andrew approves before merging. Ask before committing `review.md`. Do not push.

Read first:

1. `paper/writeup/implementation_plan.md`: **Author feedback on the initial draft (2026-10-14)** (this section takes precedence over earlier plan requirements where they conflict), Review focus, Global constraints, and the Revisit mapping.
2. `paper/writeup/source_inventory.csv` (62 entries; each records what was verified).
3. The draft sources and both rendered outputs.

## Part A: scientific and numerical correctness

Check:

- **Quantitative claims.** Check every quantitative claim against its recorded source path and selector: values, denominators, cohort sizes and overlap, events, comparator identity (the principal comparator is 40 unsupervised factors with LASSO Cox), risk orientation, and how each uncertainty interval was computed.
- **Equations and algorithms.** Check the algorithm descriptions and appendix equations against the implementation and the evaluated settings: single-modality `code/fit_cox_on_yf.R` and its update files; multimodal `code/fit_multimodal_yfb.R`, `code/multimodal_yfb_updates.R` and `code/multimodal_yfb_helpers.R`; realized settings in the runners named in the inventory.
- **Overstated claims.** Look especially for any suggestion that survival shapes the single-modality loadings when α_F = 0. Also look for biological interpretation that goes beyond what enrichment, subtype concordance, overlap, and unadjusted associations can support.
- **Revisit markers.** Confirm the hidden `<!-- REVISIT Rn START/END -->` markers cover every method-dependent passage, equation, algorithm, exhibit and caption. They must survive the restructuring in Part B. R9 now maps to the Discussion and Continuing work rather than an evidence-gap table.

Drafter findings to verify independently rather than accept:

- The principal 0.026 difference is a C-index difference over 616 concatenated external patients, not a fixed-effect pooled difference.
- All saved per-cohort C-indices equal their frozen-orientation values.
- The established model's C = 0.685 on the 50 matched ICGC donors has no saved source.
- Under the historical fit's convention (risk = −η), the joint coefficients of Programs 3 and 7 agree in sign with their marginal direction labels.
- The stored coefficient variances of that fit are zero.
- Program 3's enrichment includes both the classical and the basal-like DeSurv lists.
- Starting rank 7 has the best ICGC concordance among the multimodal starting ranks.

## Part B: author feedback to review against

Andrew reviewed the rendered draft and asked for the changes below. For each item, identify every affected location in the draft and propose the specific revision. Where an item conflicts with an earlier requirement in the plan or the drafting prompt, the item wins.

**Guiding principle.** The writeup is a standalone account of the current method, results, and biological findings, written so that its text, figures, and tables can move into a journal manuscript with little change. It is not meeting documentation and not an audit of earlier claims.

1. **Title and subtitle.** Title: *Survival-Supervised Bayesian Matrix Factorization for Prognostic Molecular Programs in Pancreatic Cancer* (the prelim chapter's working title with "Gene Expression" replaced by "Molecular", which covers the RNA-seq, proteomics and methylation data). Drop the subtitle.
2. **Authored appearance.** For the finalized document, remove every trace of the drafting process. That means the running header "Initial draft awaiting independent review", the status box, and any mention of drafts, reviewers, Claude, Codex, the source inventory, or "this report could not trace". The document should read as written by Andrew.
3. **No repository references in the rendered text.** Remove:
   - code scripts, file paths, function names, arguments and configuration labels (for example `fit_cox_on_yf()`, `flip = FALSE`, `config == pruned`, `pruned_expression_only`, `prune`, `ebnm_point_laplace`, `plot_real_cindex_by_kinit.R`, output directory names);
   - internal labels such as "D4", "P3", "flash k5" and "expr-only";
   - the visible source IDs (S01–S62).

   Software may be named as in a methods section (for example the R packages ebnm, flashier, fgsea, glmnet, survival). Source IDs are dropped entirely, not moved into hidden comments; the document is a step toward a journal manuscript. `source_inventory.csv` stays in the repository as the verification record but is not referenced by the document. The hidden Revisit markers are unaffected.
4. **No dated or process references.** Remove references to meeting notes, plans, progress chapters, dated decision records, and versions named by date, for example "September model", "September design noise", "the October progress chapter", "the August 21 report", "the September 4 correction", and "the decision record". Describe each item by its content instead: "the base model without an intercept or pruning", "high noise", "an earlier orientation check".
5. **Table 5 noise column** (multimodal simulation). "September" is meaningless to a reader. Label the levels by what they are, for example low and high noise as multiples of the noise standard deviation. Apply the same fix to Figure 2's panel titles and to all prose.
6. **Self-contained exhibits; define abbreviations.** Every table and figure must be readable on its own, as in a journal article.
   - Define each acronym or abbreviation in the caption or a table note. Examples: ELBO, BIC, CV, PL, PVE, HR, SD, CI, C, NES, padj, BH, ORA, fgsea, EBMF, LASSO, KM, PDAC, TCGA, ICGC, CPTAC, PurIST, EMT, ECM, GO-BP, DeSurv.
   - Table 4 (starting-rank sweep) is the clearest case.
   - List the undefined abbreviations exhibit by exhibit.
   - Captions should state the population, the metric and its orientation, and the uncertainty method, without repository selectors.
7. **Algorithms.** Do not use a box that breaks across pages; the break ruins the appearance. Numbered, indented algorithm formatting is sufficient. Recommend a format that renders well in both HTML and PDF.
8. **Figures.** Every figure should be sharp and well sized; Figure 1, the held-out survival curves by risk tertile, is currently fuzzy.
   - For each figure, assess resolution, size, font size relative to body text, and consistency of style.
   - Recommend regeneration where needed, for example vector PDF for the PDF output and high-resolution PNG for HTML. Figures may be redrawn from saved risk scores and aggregate outputs, without refitting any model.
   - Figures keep colour; all text stays black and white.
9. **Drop the evidence-gap table (current Table 11) and the claim-by-claim framing.** The writeup should not hinge on evaluating earlier claims.
   - Findings the current evidence supports belong in the Results and Discussion narrative.
   - Goals still being worked toward belong in Continuing work, written as plain-language aims. Examples: survival-guided recovery of low-variance programs, rank selection that accounts for survival and reconstruction, comparisons with unsupervised factorization and DeSurv, multimodal added value, and characterization of the final programs.
   - Propose where each current gap-table row's substance should go, or whether it should be dropped.
   - This supersedes the six-claim table required by Task 3.
10. **Section 5 becomes a manuscript-style Discussion** with subsections for:
    - what the findings mean and why they matter;
    - limitations;
    - continuing work, with the α_F = 0 statement kept here.

    Propose an outline and identify text to move or cut.
11. **Remove Appendix C.**
    - Convert the biological characterization methods table (current Table 14) into a concise methods subsection in Evaluation: prose or a compact list giving, for each analysis, input, population, statistic, multiplicity and what it supports. The scientific caveats from the drafting prompt stay (shared inputs, post hoc correction, within-training concordance, within-collection multiple-testing control, availability-driven method choice).
    - Place the single-modality simulation figure and table (current Figure 5 and Table 15) in Results Section 3.1 with proper context, or recommend dropping one if it adds little.
12. **No References section** for now. Attribute DeSurv in the text (for example "DeSurv (Young et al., submitted)") without a bibliography.
13. **Length and organization.** The document is long (21-page main text). Identify redundancy (summary box versus Discussion, prose that repeats tables, notation repeated in Table 1 and in the text) and propose cuts or reorganization. Keep important content but organize it as a manuscript would: Introduction-free for now, then Methods, Evaluation, Results, Biology, Discussion.
14. **Appendices A and B: complete, gap-free derivations.** Andrew will not review Appendix A now, so check both appendices step by step:
    - Each derivation starts from the model (likelihood, priors, Cox partial likelihood and its working approximation, variational objective) and reaches the exact pseudo-observation and standard error, $x = B/A$ and $s = A^{-1/2}$, passed to the empirical-Bayes normal-means solver.
    - Every intermediate step is shown, with no jumps a reader would question: expanding the quadratic, taking expectations under the mean-field posterior, collecting terms, completing the square.
    - Flag each place where a step is skipped or an equation is wrong, and supply the missing steps.
    - Also confirm the appendices state where the evaluated algorithm departs from the joint model (plug-in projection, α weight, α_F = 0, diagonal curvature, plug-in survival term in the multimodal removal objective).

## Output and stop

`review.md` should have:

- a short summary;
- Part A findings ranked by severity;
- Part B findings organized by the numbered items above, each with proposed text or structure;
- a list of checks that could not be completed and why.

Mark nothing as fixed. Ask Andrew before committing `review.md`. Do not push or merge.

Next step, for reference only: Claude Code reads `review.md`, verifies each finding against the sources, applies the confirmed ones (including all Part B items), records a disposition for every finding (fixed, not confirmed with reason, or deferred to the method work), re-renders HTML and PDF, and asks Andrew before merging.
