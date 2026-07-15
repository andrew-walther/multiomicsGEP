# C3 — Direction sanity check

**Question:** does the biology enriched in each active program agree with its assigned marginal
survival direction (Program 7 = Adverse, Program 3 = Protective; DECISIONS.md 2026-06-16), rather
than being an artifact of an arbitrary label?

**Method:** for each program, tie every headline pathway (T1, padj < 0.10) to a known, literature-established
direction in PDAC prognosis, and separately check whether the two independent, purely-outcome-based results
from Steps 7-8 (subtype concordance, external-cohort survival) point the same way. Flag agreement explicitly;
flag any disagreement rather than smoothing over it.

## Program 7 (Adverse)

| Headline set | Genes (representative) | Established direction |
|---|---|---|
| DeSurv_D3_BasalLikeTumor (padj=6.1e-4) | ITGA3, GAPDH, CDCP1, KRT7/17/19 | Basal-like/squamous PDAC — established worse-prognosis subtype |
| KEGG GF_RTK_PI3K / RTK_RAS_PI3K / RTK_RAS_ERK signaling (padj~1.1-1.2e-3) | MET, EGFR, AREG, EREG, TGFA, VEGFA | MET/EGFR-driven RTK signaling — canonical aggressive-tumor growth/invasion axis |
| Reactome developmental cell lineages (padj=0.0102) | KRT5/6A/7/15/17/18/19, ITGA6/ITGB1/ITGB4, LAMA3/LAMB3/LAMC2 | Basement-membrane/integrin + basal keratin program — squamous/basal identity |
| Bailey_Squamous (padj=0.0299) | BCAR3, HMGA2, SNAI2, S100A2, KRT5/6A | SNAI2/HMGA2 are canonical EMT drivers; Bailey's own "Squamous" subtype = worse prognosis |

**Verdict:** every headline pathway for Program 7 is glycolysis/EMT/MET/basal-squamous biology —
exactly what the plan predicted before this analysis was run, and exactly the literature's own
basal-like/squamous = worse-prognosis story. **Agrees with "Adverse."**

## Program 3 (Protective)

| Headline set | Genes (representative) | Established direction |
|---|---|---|
| DeSurv_D1_ClassicalTumor (padj=0.0172) | CLDN18, HNF4A, GATA6, TFF2 | GATA6/HNF4A are canonical classical/progenitor PDAC transcription factors |
| Moffitt_Classical (padj=0.0352) | TFF1/2/3, AGR2, AGR3, REG4 | TFF/AGR/REG4 secretory-differentiation markers — Moffitt's own "Classical" subtype = better prognosis |
| KEGG ITGA_B_FAK/RAC/RHOA (padj=0.0228, driven solely by SRC) | SRC | Weak/single-gene-driven; not a coherent program-level signal on its own |
| DeSurv_D3_BasalLikeTumor (padj=0.0172, secondary) | LPCAT4, EPCAM, KRT8, TFF1/3 | See note below |

**Verdict:** Program 3's dominant, most-significant association is Classical/differentiated biology —
**agrees with "Protective."** The secondary DeSurv_D3_BasalLikeTumor hit (also found independently in
Step 9's gene-overlap table, T4: 49/270 genes, p=0.0068, versus 61/270 and p=2.4e-6 for D1) is a real,
reproducible finding, not noise — but it is honestly a nuance, not a contradiction: several of the shared
genes (EPCAM, KRT8, TFF1/TFF3) are general epithelial/tumor-identity markers rather than basal-specific
drivers, plausibly explaining why a "classical, differentiated-epithelial" program picks up some overlap
with any tumor-cell gene list. Flagged here explicitly rather than omitted.

## Cross-check against the two purely outcome-based results (Steps 7-8)

These were derived without reference to any pathway or gene-set annotation, so they are an independent
check on the same direction claim:

- **Step 7 (PurIST subtype concordance, TCGA_PAAD, n=144):** Program 7 loading correlates positively with
  PurIST's basal-likelihood score (Spearman rho=+0.58, p=2.1e-14); Program 3 correlates negatively
  (rho=-0.57, p=4.9e-14). Basal-like PurIST calls are the established worse-prognosis subtype.
- **Step 8 (external cohort survival, 5 held-out cohorts):** Program 7's leading-edge signature has HR>1
  (worse survival) in 5/5 cohorts; Program 3's has HR<1 (better survival) in 5/5 cohorts.

Both independent checks agree with the pathway-based direction above. **No contradictions found across
any of the four methods (Steps 6, 7, 8, 9) run in this plan.**
