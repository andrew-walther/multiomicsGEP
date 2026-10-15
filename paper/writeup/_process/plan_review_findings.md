# Planning-page review findings

Review of `planning.qmd`/`planning.html`, `implementation_plan.md`, and `source_inventory.csv`
against `plan_review_prompt.md`, by Claude Code on 2026-10-14 (branch `writeup/draft`). Review
only: no plan files were edited and no analyses were run.

## Overall assessment

The plan is careful about statistical evidence. It names the principal comparator, the α_F = 0
configuration, CPTAC being proteomics, the 50-donor overlap, incomplete convergence, and the
approximate pruning objective. Those boundaries are well specified.

It is weakest on Andrew's feedback. It lists the biological characterization methods by name but
never says what each one takes as input, how it is computed, or what its output can support. The
multimodal labels come from a different, much thinner procedure than the single-modality labels,
and the plan doesn't say so. The plan also hasn't absorbed the 10/13 sequencing yet: there is no
Revisit list, no evidence-gap task, no `drafting_prompt.md`, and Task 4 still has Claude
dispatching the reviewer.

## Prioritized gaps

### 1. Single-modality characterization: methods are named, not specified

**Where:** §2 Evaluation framework (the "Biological characterization methods belong here"
sentence) and §4 Single-modality.

**Why it matters:** The analyses use four different inputs, and a reader can't judge a label
without knowing which one. This is what the code does (`code/pathway_enrichment.R`,
`results/benchmark_sim/make_km_figures.R`, and §2 of
`docs/reports/pathway_enrichment_report_08_21_26.qmd`):

| Analysis | Input | Procedure as run |
|---|---|---|
| fgsea | Posterior mean loadings E[F], nonnegative, over the 2,064 DeSurv-selected genes | One-sided `scoreType = "pos"`; minSize 10, maxSize 500. Collections: Hallmark, Reactome, KEGG, GO-BP, plus a custom PDAC collection (Moffitt/Bailey and DeSurv's 3 lists). BH within each program × collection, padj < 0.10, seed 1 |
| ORA (sensitivity) | Top 50/100/150 genes | Background is the 2,064-gene universe |
| Subtype concordance | Training latent activities E[L], TCGA n = 144 (not projection scores) | Spearman against PurIST.prob; Kruskal–Wallis against PurIST calls |
| DeSurv overlap | Top 270 genes by E[F] | Jaccard and hypergeometric test on the restricted background |
| External association (a) | Leading-edge gene-signature score: mean of within-cohort z-scored genes | Univariable Cox per cohort (C2/C4) |
| External association (b) | Frozen projection scores of Programs 3 and 7 | Within-cohort tertiles; pooled KM with cohort-stratified log-rank |
| Direction labels | Projection score | Marginal Cox on ZF, which disagrees in sign with the joint β̃ |

**Interpretive limits the plan doesn't carry forward:**

- **Fgsea and the overlap test are not independent.** The 8/21 report calls them "two independent
  methods," but DeSurv's lists sit in the fgsea custom collection and are also the overlap
  targets. The 2,064-gene universe is itself DeSurv's selection.
- **The C4 negative-control rule was chosen after the C2 result.** It should be described as a
  post hoc correction.
- **Subtype concordance is measured within the training cohort.**
- **BH correction is within program × collection only**, so there is no control across programs
  or collections.
- **The report's own caveat goes unmentioned.** It says these methods were "chosen for
  availability" and should be confirmed before the labels go into the manuscript. The S21
  inventory caveat mentions this, but the plan body doesn't.

**Amendment:** Add a compact characterization-methods table to §2 (or an appendix table referenced
from §2). Columns: analysis, input object, population, statistic and settings, multiplicity, what
the output supports. In §4, require each program label to cite which rows support it. State the
non-independence and post hoc points in prose.

**Type:** Writing/organization gap. The availability caveat is honest disclosure, not a reason to
run new analyses.

### 2. Multimodal characterization: labels come from a different procedure, and completed vs. proposed work isn't separated

**Where:** §4 Multimodal, Table 5, inventory rows S31–S32.

**Why it matters:** `results/multimodal_real/characterize_programs.R` runs no enrichment. Nothing
under `results/multimodal_real` calls fgsea or ORA. Here is what it does:

- **Labels** ("basal-like", "classical/gastric", "exocrine", and so on) come from reading the top
  8 positive and top 6 negative expression genes by signed loading.
- **Modality shares** are ‖L_k‖²‖F_mk‖², normalized within each modality.
- **Marginal HRs** use the projection Z_k = Σ_m Y_m F_mk, which includes the methylation term,
  standardized within each cohort. That gives 14 unadjusted tests.
- **Correspondence** with the expression-only fit and the D4 fit uses |Pearson r| (in both
  `characterize_programs.R` and `summarize_multimodal_real.R`). That drops the sign, which matters
  for signed F: an inverted match looks the same as a direct one.
- **CpG level:** no characterization exists.

The plan says "expression-defined labels," but not that the evidence behind them is weaker than
for the single-modality programs. Single-modality enrichment settings (one-sided, nonnegative
ranking) must not carry over. `results/benchmark_sim/run_ebmf_pathway_concordance.R` already shows
the correct pattern for signed loadings: sign alignment, then two-sided `scoreType = "std"`.

**Amendment:** Add a multimodal row group to the methods table from #1. Give Table 5 a "label
basis" column ("top signed genes; correspondence r"). Add a two-column completed/proposed list:

- **Completed:** top-gene reading, within-modality shares, |r| correspondence, unadjusted
  projection HRs.
- **Proposed (M8):** signed two-sided enrichment of expression loadings, CpG annotation and
  genomic-context enrichment, and checks for cell-composition and technical explanations.

Report signed r wherever correspondence is shown, or say that only |r| is available.

**Type:** Writing gap now; scientific evidence gap for M8.

### 3. No evidence-gap table for the manuscript's comparative claims, and two source traps

**Where:** "Gaps that affect interpretation" table, source map, inventory. Prompt 1, item 3 of the
sequencing section in `docs/plans/Meeting_Notes_Followup_10_2_26.md` asks for this table, but no
task produces it.

**Why it matters:** The current gaps table is organized by method issue, not by the claims the
manuscript must make. Evidence for those claims:

- **DeSurv:** There is no head-to-head prediction comparison.
  `results/benchmark_sim/outputs/desurv_comparison/` compares our own models under DeSurv-aligned
  preprocessing (the LB/YFB variants in `desurv_comparison_results.csv`). The directory name
  invites a drafter to call it a DeSurv comparison. The only DeSurv evidence is gene overlap (S35).
- **EBMF alone:** Evidence exists but isn't in the inventory:
  `results/benchmark_sim/outputs/pathway_enrichment/T5_ebmf_factor_correspondence.csv` and
  `T6_ebmf_pathway_concordance.csv` (unsupervised EBMF factors matched to Programs 3 and 7, with
  pathway concordance). Under α_F = 0 this is the honest test of "same programs as EBMF." Watch the
  name collision: `T5_enrichment_kept.csv` (S22) and `T5_ebmf_factor_correspondence.csv` both
  start with T5.
- **NMF alone:** I found no comparator.
- **Rank selection using both survival and reconstruction:** DECISIONS.md (around line 2672)
  records the reconstruction-versus-K explanation as "not yet verified." The K = 7 recommendation
  used external performance.
- **Multimodal added value:** 0.6852 vs. 0.6757 on 50 overlapping donors, with both fits
  unconverged.

**Amendment:** Add the task Prompt 1 specifies, as a table in §5 with columns claim, existing
evidence, source ID, missing. Rows: EBMF/NMF alone, two-step EBMF → Cox, DeSurv, K selection
(survival + reconstruction), multimodal added value, biology. Add T5_ebmf/T6 to the inventory with
an unambiguous selector, and add a caveat to S09/S13 that `desurv_comparison` is not a DeSurv fit.

**Type:** Writing gap (the table); scientific evidence gap for M5 and M7 (the missing rows).

### 4. No "Revisit after feature/supervised-F" mechanism or single Continuing-work subsection

**Where:** `implementation_plan.md` (section absent), §1 Methods briefs, §5.

**Why it matters:** These passages depend on the method work:

- the α / α_F paragraph and the ZF definition (M2)
- the single-modality loading update in the algorithm box and Appendix A (M2–M3)
- the rank-selection paragraph (M5)
- the comparator results (M7)
- all of §4 Biology (M8)

§5 currently mixes interpretation with prospective work. The Methods brief also discusses α_F
inline, so prospective supervision material could end up in three places.

**Amendment:** Add the Revisit section to `implementation_plan.md`, mapping each passage to its
M-step. Require one visible-in-source marker convention, for example `<!-- REVISIT: M2 -->`, that
Codex's step-3 review can grep. Split §5 into "Interpretation" and one "Continuing work"
subsection, which is the only place supervision, C4, and M5 plans are described. Methods states
only the evaluated configuration, plus one sentence pointing to Continuing work.

**Type:** Writing/organization gap.

### 5. The handoff isn't executable as written, and the 10/14 target needs an explicit cut list

**Where:** `implementation_plan.md` Task 4 and its header; `drafting_prompt.md` (doesn't exist).

**Why it matters:**

- **Conflicting review responsibility:** Task 4 has Claude "dispatch one independent scientific
  reviewer." The 10/13 sequencing gives that review to Codex in step 3, and Prompt 2 assumes
  `drafting_prompt.md` exists.
- **Timing:** As of 10/14, Tasks 2–3 haven't started.
- **Unverified sources:** Eleven inventory sources are still "extraction pending." That includes
  every biology source except S23 and S31, plus S14 (rank) and S39 (the established model's
  C = 0.685 on matched donors, whose provenance hasn't been traced).

**Amendment, minimum reinforcement:**

- Reassign Task 4 to Codex, and limit Claude's part to rendering checks plus the disposition of
  findings (step 4).
- In the drafting prompt, name the cuts:
  - reuse existing figures with corrected captions rather than redrawing them
  - reduce supporting simulations (S19, S20) and the threshold grid (S15) to one sentence each
  - Appendix A builds on the existing derivations, with a table separating the mathematical model
    from the evaluated algorithm
  - Appendix B gets the full equation-level reconciliation, since its components are newest
- Add one rule: any claim whose source is still unverified when drafting stops is omitted or
  listed in the evidence-gap table. It never appears in the text as a placeholder number.

None of this weakens a qualification.

**Type:** Writing/organization gap.

## Minor suggestions

- **Program 6 label differs between sources.** The planning page says "immune/quiescent-stroma";
  the 8/21 report says "immune infiltrate with a vascular component." Pick one and cite its source.
- **The saved D4 fit predates the sign-correction fix.** `make_km_figures.R` uses risk = −η as a
  single orientation fixed before validation. State this exactly in the prediction section, not
  just as a general orientation rule.
- **S08 notes that pruning can run after the iteration cap.** Carry this into the multimodal
  convergence text, since the real fits hit the cap.
- **`summarize_multimodal_real.R` also compares against an unsupervised flash reference**
  (`reference = flash`). S32 selects only `d4`. That flash comparison is the multimodal analogue of
  the EBMF-alone evidence in #3.
- **Excluding bibliography development shouldn't drop the DeSurv citation** (Young et al.). Per the
  10/14 authorship decision, the work is cited rather than coauthored, so the citation has to be
  there.
