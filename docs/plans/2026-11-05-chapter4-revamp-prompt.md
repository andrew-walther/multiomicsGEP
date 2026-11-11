# Follow-up prompt: Chapter 4 revamp (supervised YFB model)

> Branch routing updated 2026-11-11: `fix/yfb-supervision` was merged into
> `main` on November 9. Use `main` for new work. Branch instructions below
> preserve the original execution context.

Paste the block below into a new Claude Code session started in
`~/GithubProjects/multiomicsGEP` (branch `fix/yfb-supervision`).

---

> Revamp prelim Chapter 4 (SSBMF) so the supervised single-modality model is its centerpiece. That
> model is `fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1)`, with the YFB parameterization,
> point-Laplace priors on L and F, and a Normal prior on β. Follow
> `docs/plans/2026-11-05-chapter4-revamp-plan.md`, which has the current model and evidence, the
> target outline, the results checklist (tables T1–T5, figures F1–F11, biology B1–B5), the other
> documents the revamp affects, and the order of work. Target: done by about 2026-11-07, so I have
> a week to review before the 2026-11-16 prelim deadline.
>
> Andrew's notes on the current chapter are in plan section 0. Apply all of them:
> - the new title;
> - Introduction first, with no bold Summary paragraph;
> - no bold paragraph leads (usually plain topic sentences; a heading only where navigation needs it);
> - completed and proposed work blended across Intro / Methods / Results (with the PDAC application) / Discussion, in connected prose rather than an outline;
> - every figure legible at printed size, rebuilt as vector PDF where fuzzy, with a final review of all figures in the rendered prelim.
>
> Before planning, read:
> - `AGENTS.md`, `paper/prelim/AGENTS.md` (continuing work as methodology and expected evidence,
>   no dates), `DECISIONS.md` (entries from 2026-11-04 on), and `docs/findings/2026-11-04-*` and
>   `2026-11-05-*`;
> - `bios-dissertation/AGENTS.md` (writing style: no mannered prose; no invented numbers; never
>   delete cut material).
>
> Another session (S-11: configuration review and external scoring) may still be running here.
> Check `git log` and the newest DECISIONS/findings first. Don't rerun, retune or rescore its work,
> and don't write into its output folders; put new outputs in new dated folders. Ask me before
> saving per-patient external risk scores (checklist F5/F6/B5), and whether the prelim abstract
> should change.
>
> Then:
> 1. Go through the checklist with me. Ask which items I've already started, mark them, and confirm
>    the open decisions in plan section 5.
> 2. Produce the missing tables, figures and biology. Tick each checklist row in the plan file as it
>    lands, with its output path. Anything that can't be done in time becomes "planned evaluation"
>    in the chapter.
> 3. Rewrite `paper/prelim/project3-ssbmf.qmd` to the plan's outline: Introduction / Methods /
>    Results (simulation, then PDAC application) / Discussion, matching Chapters 2 and 3.
>    - Put the proposed work inside the sections where it belongs, for example a Methods
>      subsection on the M-modality model with modality-specific priors on F, and others on a
>      data-driven supervision strength, rank selection and the remaining comparators.
>    - Use descriptive headings and no bold run-in paragraph openers.
>    - Multimodal is proposed work only.
>    - State plainly that survival barely moves the real-data loadings at w_S = 1, and that the
>      gain comes from the joint β update.
>    - Use the human-writer skill while drafting.
>    - Move removed text and results to `paper/prelim/notes/trimmed-material-2026-11.md`, dated,
>      with reasons.
> 4. Commit the chapter on `fix/yfb-supervision` (the hook syncs it into bios-dissertation). Update
>    `assembly_manifest.json` for new figures, rebuild the full prelim in bios-dissertation, run the
>    assembly test, and report page counts.
> 5. Rewrite the Chapter 1 passages listed in plan section 4. Exact lines are in
>    `bios-dissertation/prelim/full-draft/notes/human-writer-review/round3/ch1/`.
> 6. Run the Round 3 human-writer review on the new Chapter 4: a writer subagent using
>    `bios-dissertation/prelim/full-draft/notes/human-writer-review/round3/WRITER-BRIEF.md`, then an
>    independent reviewer, then apply and report.
> 7. Update ROADMAP/DECISIONS in both repos.
>
> Present a numbered plan and wait for my approval before writing code. Use subagents where they
> cut cost. Before declaring the task done, check every checklist row and every bullet of this
> prompt.
