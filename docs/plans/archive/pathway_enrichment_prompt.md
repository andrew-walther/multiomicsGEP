# Pathway-enrichment planning prompt

Copy the block below into a **new** Claude Code session to build the implementation plan
for the gene-set / pathway enrichment feature. The session is expected to **produce a
plan, not code**, and to save that plan to `docs/plans/pathway_enrichment_plan.md`.

---

```
I want to build a PLAN (not implementation yet) for a gene-set / pathway enrichment
feature on my supervised Bayesian matrix factorization (SBMF) PDAC project. Please
read MEMORY.md, ROADMAP.md ("Pathway enrichment on D4 active factors" item), and
DECISIONS.md (2026-06-16 entry) first, then propose a plan and walk me through it
before writing any code.

GOAL
Assign biology to the two survival-active gene expression programs from the
recommended model by running pathway/gene-set enrichment on their top-weighted genes,
and check concordance with established PDAC molecular axes (Moffitt basal/classical,
Bailey 2016 subtypes) and with the DeSurv programs.

RECOMMENDED MODEL (what to enrich)
- Config D4: YFB projection predictor η=(YF)β, per-platform z-standardization,
  survival-ranked (combined mean+variance rank) gene selection, top-3000 per cohort
  before merging → 2064 genes, no cohort indicator, K=7, K_eff=2.
- Two survival-active programs. IMPORTANT direction convention (see DECISIONS.md
  2026-06-16 and the YFB KM sign/suppression memory): label by the MARGINAL
  (YF)-projection survival direction, NOT the joint-β sign.
    * Program 7 = ADVERSE (MET, ITGA3, BCAR3, glycolytic GAPDH/ENO1/TPI1/PGK1)
    * Program 3 = PROTECTIVE (epithelial: MLPH, SLC45A3, TJP3, CAPN5, ...)
  The joint-β signs (β̂₇=−0.041, β̂₃=+0.011) are opposite the marginal direction
  (suppression among correlated programs) — don't let that flip the labels.

DATA / INPUTS
- Fit object: results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_fits.rds
  (the "D4" element has $EF = p×K non-negative gene weights, $EBeta, $EL).
- Gene names align to the 2064-gene merged training set; recover via
  code/preprocess_desurv.R (needs PDAC_DATA_ROOT). Top genes per active program are
  already exported to
  presentation/walther_lab_meeting_06_18_2026/assets/active_factor_genes.csv.
- The F prior is point-exponential (non-negative), so all gene weights are ≥0 —
  rank genes by weight per program.

QUESTIONS THE PLAN SHOULD ANSWER
1. Enrichment method(s): fgsea (ranked, using the full per-program weight vector) vs
   over-representation (top-N gene list, e.g. clusterProfiler/Enrichr). Recommend one
   as primary and say why. Which gene-set collections (MSigDB Hallmark, KEGG,
   Reactome, GO-BP)? Gene ID mapping (symbols → Entrez) and the background set
   (the 2064 selected genes, not the whole genome).
2. How many genes / what ranking statistic per program; whether to enrich all K=7
   programs or just the 2 active ones (I lean: all 7 as a sanity check, focus on 2).
3. Concordance checks: correlate program loadings L̂_{i,3}, L̂_{i,7} with Moffitt
   basal/classical and Bailey subtype scores in TCGA_PAAD.
4. Cross-reference with DeSurv (Young et al., PNAS 2026 — see the DeSurv reference
   memory): do our adverse/protective programs recover similar biology (e.g. their
   D1 Classical+iCAF coupling)?
5. Outputs & placement: runner script, result tables/figures, and a short dated
   report — follow repo conventions (code/ or results/benchmark_sim/; docs/reports/
   for the report). Reproducibility (seeds, gene-set collection versions).
6. New R package dependencies (fgsea/clusterProfiler/msigdbr/org.Hs.eg.db) — flag for
   discussion before adding, per project conventions.

FIGURES, TABLES & COMPARISONS THE PLAN SHOULD SPECIFY
Figures
- F1. Per-program enrichment dot/bar plot: top enriched gene sets for Program 7
  (adverse) and Program 3 (protective) side by side; x = NES, size = set size,
  color = -log10(padj).
- F2. GSEA running-enrichment-score plots for the 2-3 headline pathways per program.
- F3. Gene-weight × leading-edge heatmap: top genes (rows) × the 7 programs (cols),
  annotated by which active-program pathway they belong to (extends the deck's
  gene-program heatmap with pathway labels).
- F4. Concordance: program loadings L̂_{i,7}, L̂_{i,3} by Moffitt basal/classical and
  by Bailey subtype (boxplots/violins with test p-values), in TCGA_PAAD.
- F5. SBMF-vs-DeSurv overlap: heatmap of leading-edge gene Jaccard (or enriched-
  pathway overlap) between our Programs 3/7 and DeSurv's D1-D3.
Tables
- T1. Top enriched gene sets per active program: collection, set name, size, NES,
  p, padj, leading-edge genes.
- T2. Top-N weighted genes per active program with weights (formalize the CSV).
- T3. Loading-vs-subtype concordance stats: correlations + Kruskal/ANOVA across
  Moffitt/Bailey subtypes, with effect sizes.
- T4. SBMF-vs-DeSurv concordance: pathway/gene overlap counts, Jaccard, and
  hypergeometric p for Program 7↔DeSurv-adverse and Program 3↔DeSurv programs.
Comparisons / controls
- C1. Active vs expression-only: do the 2 survival-active programs enrich for
  distinct, coherent biology while the other 5 programs do not (or enrich for
  generic/stromal sets)? — the enrichment analog of the supervision story.
- C2. Robustness: do the active programs' leading-edge signatures score consistently
  across the 5 external cohorts (sanity check that the biology travels)?
- C3. Direction sanity: confirm the adverse program's enriched biology (e.g.
  glycolysis/EMT/MET signaling) is consistent with worse survival, and the
  protective program's (epithelial/differentiation) with better — i.e. biology
  agrees with the marginal direction convention.
For each figure/table, the plan should state the exact input, the function/package,
and the output path.

CONSTRAINTS / STYLE
- R, matching existing project conventions (script headers, roxygen, section
  dividers); don't introduce new dependencies without flagging them.
- Plan first, get my approval, then implement step-by-step with verification.

DELIVERABLE FOR THIS SESSION
Produce a numbered implementation plan with verification checks per step, the
figure/table/comparison list above scoped concretely (inputs, methods, output paths),
and the open decisions you need me to weigh in on. **Save the finalized plan to
`docs/plans/pathway_enrichment_plan.md`** (Markdown). Do not write feature code yet.
```
