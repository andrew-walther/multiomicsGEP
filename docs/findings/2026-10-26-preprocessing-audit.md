# Preprocessing audit: single-modality and multimodal pipelines

Written 2026-10-26. Read-only audit of the code on `main`. No code, configuration or saved result was changed. Numbers in Section 4 come from a rebuild of the training matrix and the saved benchmark fit described there.

## Purpose and scope

The project fits a supervised Bayesian matrix factorization, Y ≈ L Fᵀ, with nonnegative patient activities L and gene loadings F, and a Cox model on the program scores YF. Two real-data analyses exist, and they prepare their data differently:

- **Single modality (expression only).** Training pools TCGA-PAAD (n = 144) and CPTAC-PDAC (n = 129), 273 patients in total. Five external cohorts are scored: Dijk, Moffitt (GEO array), PACA-AU array, PACA-AU RNA-seq and Puleo (array). The saved benchmark fit is the configuration labeled D4 in the code: Cox model on YF, per-cohort z-standardization, DeSurv gene selection, K = 7.
- **Multimodal (expression and methylation).** Training uses TCGA (n = 144). Validation uses ICGC PACA-AU donors with matched expression and methylation: 50 primary PDAC donors in the main set and 67 donors in a broader-histology sensitivity set.

This document records what each pipeline does, step by step and with file and line references. It then quantifies one consequence of the single-modality scaling for a model with nonnegative factors, and ends with one proposed choice per step for Andrew to decide.

Paths are relative to the repository root unless stated otherwise. `PDAC_data/` means the local data directory (`$PDAC_DATA_ROOT`). `data/` is the git-ignored local data folder in the repository.

---

## 1. Pipelines by cohort and data type

### 1.1 Single modality: per-cohort loading

Every cohort is read by `load_pdac_raw()` (`results/benchmark_sim/benchmark_helpers.R:35`). It points a temporary working directory at the data root and sources the per-cohort loader `PDAC_data/load_data_internal.R` (lines 51–52). That file is not in the repository. The copy in `data/PDAC_data/` is byte-identical to the OneDrive copy. `load_pdac_raw()` then keeps the rows flagged `keep == 1` (line 55), transposes to patients × genes (line 60), and keeps the first column of any duplicated gene symbol (lines 73–79).

`load_data_internal()` reads `PDAC_data/original/<cohort>.rds` (expression), `<cohort>.survival_data.rds` and `<cohort>_subtype.csv` (`load_data_internal.R:6–8`). The cohort-specific steps are:

| Cohort | Role | Patients kept (`load_data_internal.R`) | Gene symbols and duplicates | Log transform inside the loader |
|---|---|---|---|---|
| TCGA_PAAD | training | non-missing `Grade` (line 19) | rownames set to `featInfo$SYMBOL` (line 20). Duplicates: first occurrence kept by `load_pdac_raw` | **log2(x+1)** (line 21) |
| CPTAC | training | `histology_diagnosis == "PDAC"` (line 16) | as supplied. First occurrence kept | none (file is already log2, metadata `log.transformed = TRUE`) |
| Dijk | external | all (line 12) | as supplied | **log2(x+1)** (line 13) |
| Moffitt_GEO_array | external | `specimen_type == "Primary"` (line 35) | as supplied | none |
| PACA_AU_array | external | primary tumour and PDAC histology (lines 41–42) | duplicated symbols collapsed by the **median** (lines 44–50) | none |
| PACA_AU_seq | external | primary tumour and PDAC histology (lines 53–54) | duplicated symbols collapsed by the median (lines 57–63). Sample IDs get a `_seq` suffix (lines 67–69) | **log2(x+1)** (line 64) |
| Puleo_array | external | all (line 72) | as supplied | none |

All cohorts then drop patients without positive, non-missing survival time and event (lines 79–80).

Missing expression values are not handled. None are present: `validate_expression_inputs()` stops on any `NA` (`code/preprocess_desurv.R:10–11`), and every cohort passed it in the saved runs.

### 1.2 Single modality: training matrix (TCGA_PAAD + CPTAC, n = 273)

The runner `results/benchmark_sim/run_desurv_comparison.R:180–189` calls `preprocess_merged_cohorts()` (`code/preprocess_desurv.R:276`) with `top_n = 3000` (`config/globals.yml:141`), `selection_per_cohort = TRUE`, `selection_method = "combined_rank"`, `per_platform_standardize = TRUE`, `normalize_method = "none"` and `rank_transform = FALSE`. In order:

1. **Gene universe.** Intersect the gene symbols of the two training cohorts: 20,502 (TCGA, after duplicates are removed) and 28,057 (CPTAC) give 17,900 shared genes (`preprocess_desurv.R:290–292`).
2. **Second log transform.** For cohorts flagged `TRUE` in `cfg$pdac$platform_log_transform` (`config/globals.yml:211–218`), apply `log2_plus1_transform()` (`preprocess_desurv.R:17–21`, called at 304–305). TCGA is flagged, so TCGA receives its second log here. CPTAC is not flagged and is logged zero times by the project. Section 3 covers the double log.
3. **Screening, per cohort.** Within each training cohort, `select_top_variable_genes(..., method = "combined_rank")` ranks genes by mean and by variance separately, sums the two ranks and keeps the 3,000 genes with the smallest sum (`preprocess_desurv.R:57–63`, called at 317–321). The two selected sets are intersected (line 323), which gives **2,064 genes**. Screening sees only the 273 training patients, and each cohort's screen sees only that cohort. The TCGA ranks are computed on the double-logged values.
4. **Scaling.** `per_platform_standardize_cohorts()` (`preprocess_desurv.R:193–200`, called at 340–343) centers and scales each gene within each cohort to mean 0 and SD 1, using that cohort's own patients. The SD is floored at 1e-8.
5. **Merge.** Row-bind the two cohorts (line 345). No further normalization (lines 368–370) and no rank transform (line 389).

The resulting 273 × 2,064 matrix is passed unchanged to `fit_cox_on_yf()` (`run_desurv_comparison.R:262–265`). It has no intercept, and no cohort indicator is used for this configuration.

### 1.3 Single modality: external cohorts

Each external cohort is processed by itself (`run_desurv_comparison.R:297–305`) with `preprocess_desurv_cohort()` (`code/preprocess_desurv.R:100`):

1. Log transform if flagged (lines 110–111). Dijk and PACA_AU_seq are flagged and so receive a second log. The microarray cohorts are not flagged.
2. `select_top_variable_genes()` with `top_n = NULL`, which returns every gene (lines 42–43, called at 113). External cohorts are not screened.
3. Per-gene z-standardization using **the external cohort's own** means and SDs (lines 121–122). Training means and SDs are not used.
4. **Alignment to training** (`run_desurv_comparison.R:313–322`): intersect the 2,064 training genes with the cohort's genes, keep the shared genes in training order, and drop the matching rows of the posterior mean F. Training genes absent from a cohort are dropped, not imputed. Shared genes per cohort were 2,038 (Dijk), 2,000 (Moffitt), 2,050 (PACA-AU array), 2,040 (PACA-AU RNA-seq) and 1,989 (Puleo) (`results/benchmark_snapshot_2026-10-14/single_modality/desurv_comparison_results.csv`).
5. Risk score η = (Y_ext F_sub / ‖F‖) β, from `predict_cox_on_yf()` with the training column norms `EF_norms` (`run_desurv_comparison.R:325–326`).

PACA-AU contributes two external cohorts, the array and RNA-seq series. They share 47 donors, and one RNA-seq donor appears twice, so 616 records come from 568 distinct people. Per-cohort C-indices are unaffected. Any pooled analysis across the two series counts these donors twice unless it keeps one record per donor.

### 1.4 Multimodal: TCGA training (n = 144)

Built by `build_multiomics_cohorts()` (`code/load_multiomics_data.R:274`). Source files sit under `data/multiomicsGEP_code/` (`multiomics_data$root` in `config/globals.yml`):

**Expression**
- Source: `tcga/expr_log_normalized.rds` (`load_multiomics_data.R:155`), 20,501 genes × 150 samples. This is **not** the file the single-modality analysis uses. The generating script `data/multiomicsGEP_code/code_ebmf_tcga.R` computes TMM-normalized counts per million (lines 30–31), rescales them so the median library size matches ICGC (line 39), and takes log2(s · CPM + 1) once (line 43). The `saveRDS` that wrote the file is commented out (line 65). This provenance is inferred from that script and was not reproduced here, because `edgeR` is not installed.
- Transformations in the project: none. The file is used as supplied.
- Gene universe: genes present in both TCGA and ICGC expression (`load_multiomics_data.R:293–294`). This step uses validation gene names, not validation values.
- Screening: zero-variance genes removed (line 231), then the same DeSurv combined mean + variance rank, top 3,000 (lines 232–233, `top_n_genes_desurv`). TCGA patients only (line 298).
- Scaling and centering: **none** by default. The data stay on the nonnegative log2 scale, which `.validate_multimodal_yfb_blocks()` requires for point-exponential loadings (`code/preprocess_multimodal_yfb.R:35`).

**Methylation**
- Source: `tcga/tcga_meth_filtered_sex.rds$meth`, β-values in [0, 1] (`load_multiomics_data.R:156`), 305,492 probes × 150 samples. The file name indicates that sex-chromosome probes were removed upstream. Whether cross-reactive, SNP-affected or detection-p-value probes were also filtered is **not documented**.
- Probe universe: probes present in both the TCGA and ICGC sex-filtered files, matched by name (lines 284–286). This leaves 305,352 probes, 140 fewer than TCGA alone.
- Scale: raw β-values. No M-value or arcsine transform is applied by default.
- Missing values: KNN imputation within the cohort, `impute::impute.knn` with default arguments (k = 10), applied by `impute_methylation_knn()` (lines 130–142, call at 137). Only probes with at least one missing value are passed to `impute.knn` (line 131), so the nearest neighbours are found among those probes only. For TCGA, 29,307 values in 8,133 probes (0.07%) were imputed over the 144 patients, before screening (`results/multimodal_real/outputs/feature_sensitivity/universe_preflight/validation.csv`).
- Screening: the 10,000 probes with the largest β variance among probes with positive variance (lines 234–235, `n_cpgs: 10000` at `config/globals.yml:47`). TCGA patients only.
- Scaling and centering: none by default.

**Patients and survival**
- Survival comes from `PDAC_data/original/TCGA_PAAD.survival_data.rds` (line 157), the same file the single-modality loader uses.
- `filter_usable_survival()` drops 6 of 150 patients for missing or non-positive time (line 168).
- The 144 patients are exactly the 144 TCGA patients in the single-modality training set (identical barcodes), although the two analyses reach them by different filters.

### 1.5 Multimodal: ICGC validation (n = 50 primary, 67 all)

- Source: `icgc/icgc_matched_data.rds` (`load_multiomics_data.R:193`), with `$expr` unlogged, `$meth` β-values and `$info.expr` clinical data.
- Expression: **log2(x+1) once** (line 208). The source scale is presumably TMM-normalized, since the TCGA file was normalized to match it, but this is not documented.
- Methylation: same probe universe, then KNN imputation within ICGC on all 67 usable donors (line 209). The 17 non-primary donors therefore contribute to imputing values for the 50 primary donors. The number of imputed ICGC values was not recorded for this audit.
- Clinical: rows collapsed to one per donor (line 198), with a stop if survival conflicts. 12 of 79 donors are dropped for unusable survival (line 203). Primary PDAC means every sample of the donor is a primary tumour with PDAC histology (lines 102–103, 301).
- Screening: none. The TCGA-selected genes and probes are taken by name (lines 246–252).
- Scaling: none by default. ICGC keeps its own log2 scale, which differs from TCGA's (median 4.4 vs 4.9 on the 3,000 screened genes).
- Alignment: `align_multimodal_yfb_prediction()` (`code/preprocess_multimodal_yfb.R:97–118`) reorders columns to the training features and **stops** if any training feature is missing (lines 110–114). The single-modality pipeline drops missing genes instead.
- The ICGC donors are PACA-AU donors, so they overlap the PACA-AU external cohorts of the single-modality analysis. They do not overlap its training data.

### 1.6 Multimodal: optional input variants

`results/multimodal_real/run_multimodal_real_fit.R` reads the cache `data/processed_multiomics_tcga_icgc.rds` (lines 74–75), which holds the output of `build_multiomics_cohorts()`. Through `MMYFB_INPUT` (lines 76–93) it can, as alternatives to the default raw input:

- center every cohort on the **TCGA training means**. This differs from the single-modality pipeline, which centers each external cohort on its own means.
- apply asin(2β − 1) to methylation before centering.

Separately, `fit_multimodal_yfb()` has an optional per-feature intercept μ_m, re-estimated each sweep (`code/fit_multimodal_yfb.R:393–398`). It is off unless the "pruned" framework is chosen. The expression-only arm of the multimodal runner uses the raw log2 TCGA input, not the z-standardized input of the single-modality benchmark (runner header, lines 27–31).

---

## 2. Side-by-side comparison

| Step | Single modality (TCGA + CPTAC → 5 external) | Multimodal (TCGA → ICGC) | Same? |
|---|---|---|---|
| TCGA expression source | `PDAC_data/original/TCGA_PAAD.rds`, RSEM-type values up to 6.5 × 10⁶, 20,531 rows | `multiomicsGEP_code/tcga/expr_log_normalized.rds`, log2(scaled TMM-CPM + 1), 20,501 genes | **No** |
| TCGA patients | non-missing Grade + usable survival: 144 | matched methylation + usable survival: 144 | Same 144 patients, different rules |
| Log transforms on TCGA | two: `load_data_internal.R:21` and `preprocess_desurv.R:304–305` | zero in the project (one upstream) | **No** |
| Log transforms on other RNA-seq | Dijk and PACA-AU RNA-seq: two. CPTAC: zero (supplied logged) | ICGC: one (`load_multiomics_data.R:208`) | **No** |
| Duplicate gene symbols | first occurrence kept (TCGA, CPTAC, Dijk, Moffitt, Puleo). Median collapse in the loader (PACA-AU array and RNA-seq) | upstream files have unique symbols. No step in the project | Not comparable |
| Gene universe | TCGA ∩ CPTAC (17,900) | TCGA ∩ ICGC (validation gene names) | **No** |
| Gene screen rule | DeSurv combined mean + variance rank, top 3,000 per training cohort, then intersection → 2,064 | the same rule, top 3,000 on TCGA alone → 3,000 (zero-variance genes removed first) | Rule same, application differs |
| Patients seen by the screen | training only (each training cohort separately) | training only | Yes |
| Expression scaling | per-gene z-score within each cohort, training and external alike | none (raw log2). Optional: center on TCGA means | **No** |
| Negative values in the expression input | yes, 47% of training entries | no (nonnegativity enforced) | **No** |
| External/validation scaling | each external cohort z-scored on its own means and SDs | none. Optional centering uses training means | **No** |
| Alignment of validation genes | intersection. Missing genes dropped and the matching F rows removed | exact match required. Stops if a feature is missing | **No** |
| Missing expression values | none present. Code stops on any `NA` | none handled | Yes |
| Methylation source | n/a | `tcga_meth_filtered_sex.rds`, `icgc_meth_filtered_sex.rds` | n/a |
| Methylation scale | n/a | raw β (optional asin(2β − 1)). No M-values | n/a |
| Methylation imputation | n/a | KNN (k = 10) within cohort, before screening. Neighbours limited to probes with missing values | n/a |
| Probe filters | n/a | sex chromosomes (upstream). TCGA ∩ ICGC probe names. Positive variance. Cross-reactive and SNP filtering undocumented | n/a |
| CpG screen | n/a | top 10,000 by β variance, TCGA only | n/a |
| Intercept in the genomics model | none (z-scoring stands in for one) | none by default. Optional per-feature μ_m | **No** |

Of the 2,064 single-modality genes, 1,781 are among the 3,000 multimodal genes.

---

## 3. The double log

**Where it happens.**

| Cohort | First log2(x+1) | Second log2(x+1) |
|---|---|---|
| TCGA_PAAD (training) | `PDAC_data/load_data_internal.R:21` | `code/preprocess_desurv.R:304–305` (`preprocess_merged_cohorts`), flag at `config/globals.yml:212` |
| Dijk (external) | `load_data_internal.R:13` | `code/preprocess_desurv.R:110–111` (`preprocess_desurv_cohort`), flag at `globals.yml:214` |
| PACA_AU_seq (external) | `load_data_internal.R:64` | `code/preprocess_desurv.R:110–111`, flag at `globals.yml:217` |

CPTAC and the three microarray cohorts are logged once, upstream of the project. The flags in `globals.yml` were written on the assumption that `load_pdac_raw()` returns raw RNA-seq values. It does not, because `load_data_internal()` has already logged those three cohorts. Every script that loads data through `load_pdac_raw()` and uses `PLATFORM_LOG_TRANSFORM` inherits the double log. A rebuild confirmed that the loader output for TCGA equals one log2(x+1) of the source file.

**What it does to the scale (TCGA, 144 training patients).**

| Quantity | Min | 25% | Median | 75% | 99% | Max |
|---|---|---|---|---|---|---|
| Source file, all genes | 0 | 9.7 | 224 | 853 | 12,020 | 6.48 × 10⁶ |
| One log, all genes | 0 | 3.42 | 7.81 | 9.74 | 13.55 | 22.63 |
| Two logs, all genes | 0 | 2.14 | 3.14 | 3.43 | 3.86 | 4.56 |
| One log, the 2,064 selected genes | 0 | 9.01 | 10.40 | 11.67 | 15.36 | 22.63 |
| Two logs, the 2,064 selected genes | 0 | 3.32 | 3.51 | 3.66 | 4.03 | 4.56 |
| CPTAC (one log, as supplied), 2,064 genes | 0 | 9.78 | 11.13 | 12.32 | 15.99 | 22.53 |

The second log compresses the upper range strongly. Values from 9 to 23 on the single-log scale map to 3.3–4.6. On the selected genes the median per-gene SD falls from 0.97 to 0.12. The per-cohort z-score removes the change in location and spread, but the shape change is nonlinear and stays: high-expression patients are pulled in, and low-expression patients are spread out relative to them. The TCGA half of the training matrix is therefore on a different transformed scale from the CPTAC half, even after z-scoring.

**Effect on screening.** The ranks of per-gene SD under one and two logs agree closely on the selected genes (Spearman 0.97), but the screen still changes. With TCGA logged once, the same per-cohort rule selects 2,218 shared genes instead of 2,064, and 1,985 of the 2,064 frozen genes are among them. The effect on the fitted benchmark is measured in [the single-log check](2026-10-26-single-log-check.md): external mean C 0.625 versus 0.627, with risk scores correlated at 0.99 or more.

---

## 4. Negative entries and reconstruction error in the frozen single-modality fit

**Method.** The training matrix was rebuilt with the same `preprocess_merged_cohorts()` call as the `prepare` block of `results/benchmark_sim/run_feature_expansion.R`. The script asserts that the gene names equal the frozen 2,064 (`d4_gene_names.rds`) and that n = 273. The frozen fit `results/benchmark_snapshot_2026-10-14/single_modality/desurv_comparison_fits.rds`, element `$D4`, stores the posterior means `EL` (273 × 7) and `EF` (2,064 × 7). `fit_cox_on_yf()` measures reconstruction as Y − EL EFᵀ on the input Y (`code/fit_cox_on_yf.R:1020`). This configuration has no cohort columns, so Ŷ = EL EFᵀ. The rebuilt RMSE, 0.932962, equals the fit's last recorded `history$rmse`, which confirms that the matrix and the reconstruction are on the scale the fit used. Four of the seven factors have nonzero loadings. Only training data were used.

**Results.**

| Quantity | Value |
|---|---|
| Share of training entries with Y < 0 | **47.2%** (TCGA 45.5%, CPTAC 49.1%) |
| Minimum of EL and EF | 0 (both nonnegative). Minimum of Ŷ: 1.2 × 10⁻⁷ |
| Mean Ŷ on entries with Y < 0 / Y ≥ 0 | 0.070 / 0.266 |
| Total sum of squares of Y | 559,344 |
| Residual sum of squares (SSE) | 490,456 (proportion of variance explained 0.123) |
| **Share of SSE on entries with Y < 0** | **68.7%** |
| **Mean residual Y − Ŷ, entries with Y < 0** | **−0.886** |
| **Mean residual Y − Ŷ, entries with Y ≥ 0** | **+0.465** |
| Lowest possible SSE on the Y < 0 entries for any Ŷ ≥ 0, which is Σ Y² over those entries | 310,032, or 63.2% of the fit's SSE and 55.4% of the total sum of squares of Y |
| SSE on Y < 0 entries above that floor | 5.5% of SSE |

**What this implies.** A product of nonnegative L and F cannot take a value below zero. After per-gene centering, about half of every gene's entries are negative by construction, so the reconstruction can only approach zero there. On those entries the fit does nearly the best a nonnegative model can do, with Ŷ close to 0, and its error is almost exactly the floor Σ Y². That floor alone is 55% of the total variation in Y and 63% of the fit's residual error. The remaining error is spread over the positive entries, where the mean residual is +0.47.

Several consequences follow.

- The genomics part of the model describes only above-average expression. A program defined by genes below their mean in some patients can be represented only by an activity near zero, not by a negative contribution.
- The proportion of variance explained (12%) is capped well below what a signed model of the same rank could reach. Most of the residual reflects the sign constraint, not noise.
- The estimated noise precisions are set by this structured error. Tau is about 1.0 per gene, and the mean squared residual is 0.87, so most of what the model treats as noise comes from the sign floor. A smaller precision gives the reconstruction term less weight relative to the survival term, and the same error enters any ELBO-based choice of K.
- Risk prediction uses YF, which includes the negative entries, so the survival scores see the full centered data. The fitted F, however, was learned under the distorted reconstruction.

The input scale and the loading prior disagree. Either the data should stay nonnegative, for example single-logged expression with a per-feature intercept modeled explicitly, or the loadings should be signed.

---

## 5. Proposed alignments

The recommendations below were proposed on October 26. The approved October 30
dispositions follow this historical list; they supersede its pending-decision labels. The saved benchmark numbers must not be quoted for inputs that differ from those that produced them.

1. **Log transform: exactly once for every RNA-seq cohort.** Remove the `platform_log_transform` step for cohorts that `load_data_internal()` already logs (TCGA, Dijk, PACA-AU RNA-seq), or move all logging into the project code and leave the loader to read files only. Add a test that a raw input comes out as log2(x+1) applied once. *For Andrew's decision.*
2. **One TCGA expression source in both analyses.** Either use `TCGA_PAAD.rds` (RSEM, logged once) for the multimodal analysis as well, or use the TMM-CPM file for the single-modality analysis. The first keeps the five-cohort benchmark intact. The second matches the normalization intended for ICGC. *For Andrew's decision.*
3. **The same expression scaling in both analyses.** Options: (a) single-logged, nonnegative expression with a per-feature intercept μ estimated in the model and no z-scoring, which suits nonnegative loadings, (b) per-cohort centering and scaling, which requires signed loadings, or (c) keep z-scoring with nonnegative loadings and accept the sign floor quantified in Section 4. Section 4 argues against (c). *For Andrew's decision.*
4. **How validation cohorts are scaled.** Choose between each cohort's own means and SDs (current single-modality practice, which removes cohort-level shifts but uses the validation cohort's distribution) and the training means and SDs (current multimodal centering option, which is a strict frozen transform but carries platform offsets into the validation data). Apply the choice to both analyses. *For Andrew's decision.*
5. **Gene universe and screen.** Keep the DeSurv combined rank. Decide whether the gene universe may use validation gene names, as the multimodal pipeline does, or only training cohorts, as the single-modality pipeline does. Decide whether pooled training should screen per cohort and intersect, or screen once on the pooled data. *For Andrew's decision.*
6. **Missing validation features.** Use one rule in both analyses: drop missing features and the matching F rows (single modality), or require all training features (multimodal). The second is stricter and reproducible. The first is needed for external cohorts with smaller gene sets. *For Andrew's decision.*
7. **Duplicate gene symbols.** Use one rule for all cohorts, either the first occurrence or the median across duplicates. *For Andrew's decision.*
8. **Methylation scale.** Keep β-values (nonnegative, bounded), or use M-values or asin(2β − 1) with a signed methylation loading prior. *For Andrew's decision.*
9. **Methylation imputation.** Keep within-cohort KNN, but pass the full probe matrix, or a fixed neighbour pool, to `impute.knn` so that neighbours are not limited to probes with missing values. Impute the validation cohort on the reported validation set (50 primary donors) or state that all 67 donors are used. *For Andrew's decision.*
10. **Probe filters.** Establish whether cross-reactive and SNP-affected probes were removed upstream. If not, decide whether to add those filters before screening. *For Andrew's decision.*
11. **PACA-AU donor accounting.** Keep per-cohort C-indices. In any pooled external analysis, use one record per donor (568 people from 616 records). *For Andrew's decision.*

---

## Not determined in this audit

- The upstream processing of the microarray cohorts (Moffitt, PACA-AU array, Puleo) and of CPTAC, beyond the fact that they are supplied on a log scale.
- The exact provenance of `tcga/expr_log_normalized.rds`. It was inferred from `code_ebmf_tcga.R` and not reproduced, because `edgeR` is not installed.
- Upstream methylation probe QC beyond sex-chromosome removal.
- The number of imputed ICGC methylation values.
- The effect of any proposed change other than the single log (measured separately in [the single-log check](2026-10-26-single-log-check.md)) on fitted programs or external C-index.

## October 30 dispositions for Section 5
The approved realignment plan adopts single logging, the existing TCGA RSEM
source (`original/TCGA_PAAD.rds`), and training-only gene names selected by
per-cohort combined mean/variance ranks with intersection. For the new fitter
only, divide nonnegative single-logged expression by each cohort's feature SD
without centering; center inside the survival projection. Validation uses its
own within-cohort SD, a distribution-adaptive transform that requires a cohort
rather than a single new patient. No external outcomes enter scaling or choice.

Keep the D4 missing-feature rule (intersect genes and subset matching loading
rows). Keep the wrapper's first-occurrence duplicate-symbol rule. The internal
PACA-AU loaders already median-collapse symbols before that wrapper; adopting
one rule across raw sources will require removing that earlier aggregation
when those pipelines are next changed. No such loader change is authorized here.
Keep methylation β-values and within-cohort KNN imputation; the primary set has
50 donors, while the existing imputation pool contains all 67 usable donors.
Cross-reactive and SNP probe filtering is required when multimodal work resumes.
Pooled external analyses must use one record per donor (568 unique individuals),
with donor-level or cohort-stratified resampling. These records do not authorize
a new external analysis. Survival pre-screening remains an ablation only.
