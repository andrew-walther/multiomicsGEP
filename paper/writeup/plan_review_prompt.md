# Claude Code prompt: review the HTML writeup plan before approval

Review the proposed writeup plan and identify gaps or areas to reinforce **before Andrew officially approves it**. This is a review-only task. Do not execute the implementation plan or draft the writeup.

## Read first

- `CLAUDE.md` and the current writeup evidence-boundary and authorship entries in `DECISIONS.md`.
- `paper/writeup/planning.html`, the document Andrew is reviewing, together with its source `paper/writeup/planning.qmd`.
- `paper/writeup/implementation_plan.md` and `paper/writeup/source_inventory.csv`.
- The **Sequencing and handoffs** section of `docs/plans/Meeting_Notes_Followup_10_2_26.md`. It records the newer plan for Claude to draft and Codex to review, which may not yet be reflected throughout the planning materials.

Run `git status -sb` and `git branch --show-current` to establish the current state. Do not switch branches, create branches/worktrees, reset, stash, commit, push, or merge for this review. Preserve any existing changes.

## Review focus

The planned document is an advisor-facing, article-style technical report: a 15–20-page main text covering current single-modality and multimodal methodology, simulation and real-data performance, and biology, plus derivation appendices. It should support eventual manuscript writing and journal discussion. A literature review, bibliography expansion, journal selection, new fitting experiments, and method development are outside this drafting scope.

Andrew's feedback is that **the writeup must explain the methods used for biological characterization, not merely describe the biological function of each factor/program**. Check whether the plan gives sufficient direction to explain:

- Which inputs each characterization uses: posterior gene loadings, latent activities, projected program scores, or gene-signature scores.
- Gene-set enrichment and over-representation: ranking/sign conventions, collections, gene backgrounds, statistical comparisons, multiple testing, thresholds, leading-edge genes, and interpretation.
- Subtype concordance, DeSurv gene-list overlap, external survival associations and controls: how each analysis is performed, what its outputs mean, and how those outputs support program labels.
- Multimodal characterization: program correspondence, modality contributions, signed loadings, and which CpG-level characterization is completed versus still proposed.

Also assess:

1. **Coverage and coherence:** Are models, algorithms, prediction, evaluation, results, biology, and derivations adequately covered and connected? Does the exhibit/page plan support the intended audience?
2. **Theory and implementation:** Does the plan distinguish the evaluated single-modality alpha_F=0 algorithm from the intended joint formulation, and account for the multimodal intercept, signed priors, uncertainty, scaling, and pruning?
3. **Evidence and claims:** Are comparator identities, rank definitions, cohort overlap, uncertainty, convergence, and biological interpretations handled accurately? Identify claims needing more support, especially value over EBMF/NMF alone, two-step EBMF → Cox, and DeSurv; rank selection accounting for both survival and reconstruction; and multimodal added value.
4. **Future updates:** Does the plan make it easy to revisit method-dependent passages after `feature/supervised-F`, while keeping prospective supervision work concentrated in one Continuing work subsection?
5. **Handoff feasibility:** Are drafting versus independent-review responsibilities clear? What minimum reinforcement is needed for the existing October 14 draft / October 16 sharing targets, without adding new analyses or weakening scientific qualifications?

Inspect relevant repository sources only as needed to substantiate a concern. For biological methods, start with `docs/reports/pathway_enrichment_report_08_21_26.qmd`, `code/pathway_enrichment.R`, and `results/multimodal_real/summarize_multimodal_real.R`. Do not assume that an analysis exists merely because a script name suggests it, or transfer settings from nonnegative single-modality loadings to signed multimodal loadings without checking.

## Return findings in chat, then stop

Give a brief overall assessment and at most **five substantive gaps or recommended reinforcements**, prioritized by impact. For each, provide:

- The specific plan section or omission.
- Why it matters and the relevant source evidence, if applicable.
- A concrete, minimal suggested amendment.
- Whether it is a **writing/organization gap** or a **scientific evidence gap for later method work**.

List minor suggestions separately only if useful. Do not manufacture issues or repeat requirements already adequately covered. Missing scientific evidence should usually be described honestly in the writeup, rather than becoming a demand for new experiments before drafting.

Do not edit any files, revise the plan, write the report, execute its tasks, or run new analyses. Andrew will bring your findings back to Codex to discuss and revise the plan before formal approval.
