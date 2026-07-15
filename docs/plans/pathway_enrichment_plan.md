# Implementation Plan — Pathway / Gene-Set Enrichment on D4 Active Programs

**Status:** Draft for approval (planning session, 2026-06-16). No feature code written yet.
**Roadmap item:** "Pathway enrichment on D4 active factors" `[Priority: High]` `[Effort: Small]`
**Goal:** Assign biology to the two survival-active gene expression programs of the recommended
model (D4) by running gene-set enrichment on their gene-weight vectors, and test concordance
with established PDAC molecular axes (Moffitt basal/classical, Bailey 2016) and with the DeSurv
programs (Young et al., PNAS 2026).

---

## 0. Grounding facts (verified this session)

| Fact | Value / location | Implication for plan |
|------|------------------|----------------------|
| Fit object | `results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_fits.rds`, element `D4` | `$EF` = 2064×7 non-negative gene weights, `$EBeta` = length-7, `$EL` = n×7 |
| Gene names | `results/benchmark_sim/outputs/desurv_comparison/d4_gene_names.rds` (length 2064, gene symbols) | **Gene names already recovered — `PDAC_DATA_ROOT` is NOT needed for the enrichment step.** It is needed only for the concordance checks (TCGA_PAAD expression + clinical/subtype). |
| Top genes (exported) | `presentation/walther_lab_meeting_06_18_2026/assets/active_factor_genes.csv` | Program 3 (Protective): CAPN5, LPCAT4, MLPH, SLC45A3, TJP3, PLEKHA6, MTMR11, GALNT6. Program 7 (Adverse): ITGA3, MET, BCAR3, FAM3C, GAPDH, ENO1, TPI1, PGK1. |
| Direction convention | **Program 7 = ADVERSE, Program 3 = PROTECTIVE** by the marginal (YF)-projection survival direction (DECISIONS.md 2026-06-16) | Joint-β signs (β̂₇=−0.041, β̂₃=+0.011) are *opposite* and must **not** drive labels (suppression among correlated programs). |
| F prior | point-exponential (non-negative) → all weights ≥ 0 | Rank genes by weight; enrichment is effectively **one-sided** (no "depleted" sets). Drives the method recommendation in §2. |
| Existing doc | `docs/reports/pathway_enrichment_overview.qmd` (2026-05-30, 910 lines) | A **conceptual primer/illustration** ("what it is, why it matters"), NOT the analysis. **It predates the 2026-06-16 correction** — its code comment labels Factor 3 as adverse (β>0). Keep it as a primer; the new report is the analysis. Do not inherit its labeling. |
| Package status | `fgsea`, `clusterProfiler`, `msigdbr`, `org.Hs.eg.db`, `ComplexHeatmap` all **MISSING**; `preprocessCore`, `flashier`, `survival`, `ggplot2`, `pheatmap` installed | New dependencies required — see Open Decision **D1** (must approve before any install). |

---

## 1. Deliverables & placement (repo conventions)

| Artifact | Path | Convention followed |
|----------|------|---------------------|
| Reusable functions (roxygen, headers) | `code/pathway_enrichment.R` | `code/` for reusable, documented functions |
| Orchestration runner | `results/benchmark_sim/run_pathway_enrichment.R` | matches `run_desurv_comparison.R` etc. |
| Outputs (tables, figures, manifest) | `results/benchmark_sim/outputs/pathway_enrichment/` | `outputs/<analysis>/` |
| Dated analysis report | `docs/reports/pathway_enrichment_report_06_DD_26.{qmd,pdf,html}` | `docs/reports/<name>_MM_DD_YY` |
| Custom PDAC gene-set object | `results/benchmark_sim/outputs/pathway_enrichment/pdac_genesets.rds` | reproducibility |
| Gene-set version manifest | `.../pathway_enrichment/genesets_manifest.txt` | seeds + collection versions |

---

## 2. Methodological choices (recommendations + rationale)

**Primary enrichment method — `fgsea` ranked by gene weight (recommended).**
Because the point-exponential F prior makes every weight ≥ 0, the per-program weight vector is a
natural continuous ranking statistic. Ranking all 2064 genes by `EF[, k]` (descending) and running
`fgsea` uses the full weight information and avoids an arbitrary top-N cutoff. The enrichment is
**one-sided by construction** — we report positive NES only; "depletion" is meaningless for a
non-negative weight (state this explicitly in the report). This matches the prompt's lean toward
using the full weight vector.

**Confirmatory method — over-representation analysis (ORA, hypergeometric) on top-N genes.**
Run `clusterProfiler::enricher()` on the top-N weighted genes per program against the **2064 selected
genes as background** (not the genome — using the genome inflates significance because the 2064 are
already survival/variance-selected). N swept over {50, 100, 150} for robustness. ORA is the cross-check
that fgsea's leading edge is not a ranking artifact.

**Gene-ID mapping.** Gene symbols → Entrez via `org.Hs.eg.db` (`clusterProfiler::bitr`), keeping a
mapping table and logging unmapped symbols (fail-loud: report count + names of dropped genes; do not
silently drop). `fgsea` runs on symbols directly when using `msigdbr` symbol sets — prefer symbol-space
to minimize mapping loss; use Entrez only where a collection requires it.

**Gene-set collections (via `msigdbr`, species = "Homo sapiens"):**
- **MSigDB Hallmark (H)** — primary, compact, interpretable (50 sets). EMT, glycolysis, hypoxia expected for Program 7.
- **Reactome (C2:CP:REACTOME)** and **KEGG (C2:CP:KEGG_MEDICUS)** — pathway-level detail.
- **GO Biological Process (C5:GO:BP)** — broad coverage; report top sets only (large, redundant).
- **Custom PDAC sets** (built into `pdac_genesets.rds`): Moffitt basal-like (25) + classical (25)
  signatures (Moffitt et al. 2015); Bailey 2016 subtype signatures; DeSurv programs (Young et al. 2026)
  — pending source availability (Open Decision **D3**).

**Multiple-testing.** BH-adjusted p (`padj`) per collection; report `padj < 0.10` with NES and
leading-edge genes. Min/max set size filters: 10/500.

**Scope of programs.** Enrich **all K=7 programs** as a sanity check (expect the 5 inactive programs
to enrich for generic/stromal/no coherent biology — this is comparison **C1**, the enrichment analog
of the supervision story), then **focus reporting on Programs 3 and 7**.

**Reproducibility.** `set.seed()` before fgsea permutations; record `msigdbr` version, `fgsea` version,
R version, and collection retrieval date in `genesets_manifest.txt`.

---

## 3. Numbered implementation steps (with verification)

> Each step lists: action → **verify:** check. Commit between steps (per project git convention).

**Step 1 — Dependencies & environment.** *(Gated on Open Decision D1.)*
Install `fgsea`, `clusterProfiler`, `msigdbr`, `org.Hs.eg.db` (Bioconductor); optionally `ComplexHeatmap`.
Register them in the report/runner header. → **verify:** `requireNamespace()` returns TRUE for each;
`Rscript tests/run_tests.R` still 246/246 (no namespace clashes introduced).

**Step 2 — `code/pathway_enrichment.R`: load + rank.**
Function `load_d4_weights()` reads the fit + gene names, attaches `rownames(EF)`, returns a list with
`EF` (2064×7, named rows), `EBeta`, and a `program_labels` lookup hard-coding the **corrected** convention
(7→Adverse, 3→Protective; others→inactive). → **verify:** unit test asserts `nrow(EF)==2064`,
`length(gene_names)==2064`, no duplicate gene symbols (or de-dup rule applied + logged), and
`program_labels[["7"]]=="Adverse"`, `program_labels[["3"]]=="Protective"` (guards against the stale label).

**Step 3 — `code/pathway_enrichment.R`: enrichment engines.**
`run_fgsea_program(weights_k, genesets, seed)` and `run_ora_program(top_genes, background, genesets)`.
Both return tidy data frames (collection, set, size, NES/oddsratio, p, padj, leading_edge). Fail loud on
zero overlap or all-NA results. → **verify:** unit test on a tiny synthetic gene set where the answer is
known (e.g. a set = the top-5 genes of a program must enrich with NES>0, p small); a test that passes
when the function returns a constant must FAIL (assert NES differs between a matched vs scrambled set).

**Step 4 — Build custom PDAC gene sets.** *(Partly gated on D3.)*
`build_pdac_genesets()` assembles Moffitt basal/classical, Bailey, and DeSurv signatures into a named list,
saved to `pdac_genesets.rds`. → **verify:** each set non-empty, gene symbols present in MSigDB universe,
manifest records the source citation per set.

**Step 5 — Runner: enrich all 7 programs.**
`run_pathway_enrichment.R` loads weights, retrieves collections, runs fgsea (all 7) + ORA (N∈{50,100,150},
Programs 3 & 7), writes **T1** and **T2** CSVs and the fgsea/ORA result objects (.rds). → **verify:** output
CSVs exist and non-empty; Programs 3 & 7 each return ≥1 set at `padj<0.10`; the 5 inactive programs
summarized for **C1** (count of coherent hallmark hits — expected low / generic).

**Step 6 — Figures F1–F3 (enrichment).**
F1 dot/bar (top sets, Programs 3 vs 7), F2 GSEA running-ES for 2–3 headline sets/program,
F3 gene-weight × leading-edge heatmap (genes × 7 programs, annotated by active-program pathway).
→ **verify:** files written to `outputs/pathway_enrichment/`; F1 x-axis = NES, size = set size,
color = −log10(padj) as specified; visual sanity (Program 7 shows glycolysis/EMT/MET-type sets).

**Step 7 — Concordance with PDAC subtypes (F4, T3).** *(Gated on `PDAC_DATA_ROOT` + Open Decision D2.)*
Load TCGA_PAAD expression + clinical; obtain Moffitt basal/classical score and Bailey subtype labels
(see D2 for label source). Correlate patient loadings `EL[,3]`, `EL[,7]` (continuous molecular scores)
with the Moffitt basal-classical axis (Spearman) and compare across Bailey subtypes (Kruskal–Wallis,
effect size). → **verify:** sample barcodes align (report N matched / N dropped, fail loud if <80% match);
Program 7 loading should track the basal/squamous axis and Program 3 the classical/progenitor axis if
the biology is coherent — state the directional prediction *before* computing.

**Step 8 — Robustness across external cohorts (C2).**
Score each program's leading-edge signature (mean of z-scored leading-edge genes) in each of the 5
held-out cohorts; relate the signature score to survival (Cox HR / C-index per cohort). → **verify:**
adverse-program signature → HR>1 (worse survival) in a majority of cohorts; protective → HR<1; report
the count and fail loud if signatures are unavailable in a cohort (missing genes logged, not silently
zero-filled).

**Step 9 — SBMF vs DeSurv overlap (F5, T4).** *(Gated on D3 — DeSurv gene lists.)*
Compute leading-edge gene Jaccard and enriched-pathway overlap between Programs 3/7 and DeSurv D1–D3;
hypergeometric p for overlap against the 2064 background. → **verify:** Jaccard/overlap matrices written;
hypergeometric p reported with the background size stated.

**Step 10 — Direction sanity (C3).**
Confirm Program 7's enriched biology (glycolysis / EMT / MET signaling) is consistent with *worse*
survival and Program 3's (epithelial / differentiation) with *better* — i.e. biology agrees with the
marginal direction convention. → **verify:** narrative cross-check in the report ties each headline
pathway to its program's survival direction; flag any contradiction explicitly rather than smoothing it.

**Step 11 — Report.**
`docs/reports/pathway_enrichment_report_06_DD_26.qmd` assembling F1–F5, T1–T4, C1–C3, methods (collections,
versions, background, seeds), and a "direction convention" callout. Render PDF+HTML. → **verify:**
`quarto render` succeeds; all figures/tables present; reproducibility manifest embedded.

**Step 12 — Living docs.**
Update `ROADMAP.md` (check off the item), `DECISIONS.md` (enrichment method choice, gene-set versions,
subtype-label source, any new dependency), `PROJECT_STATUS.md`, and CLAUDE.md quick-reference table.
→ **verify:** links resolve; 246/246 tests still pass.

---

## 4. Figures, tables & comparisons — concrete spec

### Figures
| ID | Content | Input | Function/package | Output path |
|----|---------|-------|------------------|-------------|
| F1 | Top enriched sets, Program 7 vs 3, side by side; x=NES, size=set size, color=−log10(padj) | fgsea results (Step 5) | `ggplot2` | `outputs/pathway_enrichment/F1_enrichment_dotplot.{png,pdf}` |
| F2 | GSEA running-ES plots, 2–3 headline sets/program | fgsea + ranked weights | `fgsea::plotEnrichment` | `.../F2_running_es_*.{png,pdf}` |
| F3 | Gene-weight × leading-edge heatmap (top genes × 7 programs), annotated by active-program pathway | `EF` + leading-edge sets | `pheatmap` (or `ComplexHeatmap` if D1 approves) | `.../F3_geneweight_heatmap.{png,pdf}` |
| F4 | `EL[,7]`, `EL[,3]` by Moffitt basal/classical and by Bailey subtype (violin + test p) | TCGA_PAAD loadings + subtypes (Step 7) | `ggplot2` | `.../F4_loading_vs_subtype.{png,pdf}` |
| F5 | SBMF×DeSurv leading-edge Jaccard (or pathway-overlap) heatmap | Step 9 | `pheatmap` | `.../F5_sbmf_vs_desurv.{png,pdf}` |

### Tables
| ID | Content | Input | Output path |
|----|---------|-------|-------------|
| T1 | Top enriched sets/active program: collection, set, size, NES, p, padj, leading-edge | fgsea (Step 5) | `.../T1_enrichment_active.csv` |
| T2 | Top-N weighted genes/active program + weights (formalizes the CSV) | `EF` ranked | `.../T2_top_genes.csv` |
| T3 | Loading-vs-subtype concordance: Spearman ρ + Kruskal/ANOVA across Moffitt/Bailey, effect sizes | Step 7 | `.../T3_concordance_stats.csv` |
| T4 | SBMF-vs-DeSurv overlap: counts, Jaccard, hypergeometric p (Program 7↔DeSurv-adverse, 3↔DeSurv) | Step 9 | `.../T4_sbmf_desurv_overlap.csv` |

### Comparisons / controls
- **C1 — Active vs inactive programs:** do Programs 3 & 7 enrich for distinct coherent biology while
  the other 5 do not? (enrichment analog of the supervision story). Output: summary table of coherent-hit
  counts per program.
- **C2 — External robustness:** leading-edge signatures score consistently (correct survival direction)
  across the 5 external cohorts.
- **C3 — Direction sanity:** enriched biology agrees with the marginal survival direction (adverse=worse,
  protective=better).

---

## 5. Decisions — RESOLVED (2026-06-16)

**D1 — New R dependencies — APPROVED.** Install `fgsea`, `clusterProfiler`, `msigdbr`,
`org.Hs.eg.db` (Bioconductor); `ComplexHeatmap` optional for F3 (`pheatmap` fallback).

**D2 — Subtype labels — RESOLVED (data confirmed local).** Per-cohort subtype files live alongside
the expression data: `$PDAC_DATA_ROOT/original/<cohort>_subtype.csv` and `<cohort>.caf_subtype.rds`
(loaded by `load_data_internal.R`). For TCGA_PAAD:
- **Moffitt basal/classical axis — already computed per-sample.** Use `MS` / `MS_K2` labels and
  **`PurIST`** (label + continuous `PurIST.prob` in `TCGA_PAAD.caf_subtype.rds$Subtype`, keyed by `sampID`).
  `PurIST.prob` gives a continuous basal↔classical score for the Spearman correlation in T3; `MS_K2`
  gives the 2-group split for the violin/Kruskal test in F4. **No signature re-derivation needed.**
- **Bailey 4-subtype — registered schema, labels NOT pre-computed.** `cmbSubtypes.RData$schemaList`
  includes `"Bailey"` (and `"Collisson"`, `"MS_K2"`), but no per-sample Bailey column exists in the
  CSV or `Subtype` df, and `subtypeGeneList` is empty in that file. **Plan:** primary concordance uses
  the Moffitt axis (PurIST/MS_K2, fully available); attempt Bailey via the classifier referenced by the
  schema if it is readily callable, otherwise **defer Bailey to a follow-up** and note it in the report.
  (Decision point retained in Step 7.)

**D3 — DeSurv gene lists — RESOLVED (source confirmed local).**
`UNC Dissertation (Liu)/papers/DeSurv/` contains `si_appendix.pdf`, `paper.pdf`, and
`DeSurv_Paper.html` (links to the R-package + paper git repos). DeSurv program gene lists will be
**extracted from `si_appendix.pdf`** (or pulled from the linked repo if the SI is not machine-readable).
F5/T4 proceed; if the SI lists are not cleanly extractable, fall back to the repo gene lists and log the source.

**D4 — Primary method — CONFIRMED.** fgsea ranked-by-weight primary; ORA top-N confirmatory.

**D5 — Report vs overview — CONFIRMED.** New dated analysis report; keep
`pathway_enrichment_overview.qmd` as the conceptual primer (predates the direction correction).

---

## 6. Constraints / style
- R, matching project conventions: script headers, roxygen on every function, `# Section ----` dividers.
- No new dependencies added until **D1** is approved.
- Fail loud: log dropped/unmapped genes, unmatched barcodes, missing cohort genes — never silently zero-fill.
- Plan first → approval → implement step-by-step with the verification checks above.
