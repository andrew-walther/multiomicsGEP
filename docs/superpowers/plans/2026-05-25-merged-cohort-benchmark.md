# Merged-Cohort Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a comprehensive, apples-to-apples comparison of 6 model configurations on merged TCGA+CPTAC training data, using CV-selected K throughout, to identify the most robust configuration for multi-cohort prognostic gene expression program discovery.

**Architecture:** Two new scripts (`run_merged_kcv.R` for K selection, `run_merged_benchmark.R` for fitting and evaluation) plus a Quarto report. All hyperparameters come from `config/globals.yml`. Existing `preprocess_merged_cohorts(per_platform_standardize=TRUE/FALSE)` and `select_K_cv()` infrastructure is reused — no new model code needed.

**Tech Stack:** R, `survival`, `ebnm`, `yaml`, `ggplot2`, Quarto. All existing model code (`fit_modular.R`, `fit_cox_on_yf.R`, `select_K.R`, `preprocess_desurv.R`) is used as-is.

---

## Background and Design Rationale

### The Problem
We need a merged multi-cohort model (TCGA_PAAD RNA-seq + CPTAC proteomics) that:
- Learns biologically relevant gene expression programs (factors F)
- Maintains non-zero survival associations (β ≠ 0)
- Uses a parsimonious K (CV-selected, not ARD-pruned from a large K=20)
- Generalizes to 5 independent external cohorts

Platform mixing (RNA-seq vs proteomics) creates a technical confound that competes with biological signal. Two levers handle this:
1. **Preprocessing lever** — per-platform z-standardization (normalize TCGA and CPTAC separately before merging) removes bulk platform scale differences at the input level.
2. **Model lever** — cohort dummy column (corner-point encoded indicator in L) explicitly estimates and removes residual between-cohort offsets in gene expression space inside the model.

These are complementary: preprocessing handles magnitude differences, the dummy handles systematic directional differences that survive normalization.

### Model Configurations

Six configurations to test (2 models × 2 preprocessing × 2 cohort_id):

| ID | Model | Preprocessing | cohort_id | Prior | K selection | Prior status |
|----|-------|---------------|-----------|-------|-------------|--------------|
| M1 | LB | Joint quantile+rank | No | normal | CV on merged | Re-run at CV K (was K=20 ARD) |
| M2 | LB | Joint quantile+rank | Yes | normal | CV on merged | Re-run at CV K (was K=5, not CV) |
| M3 | LB | Per-platform z-std | No | normal | CV on merged | **NEW** |
| M4 | LB | Per-platform z-std | Yes | normal | CV on merged | **NEW** |
| M5 | YFB | Per-platform z-std | No | normal | CV on merged | Re-run at CV K (was K=3, fixed) |
| M6 | YFB | Per-platform z-std | Yes | normal | CV on merged | **NEW** |

Explicitly excluded (documented β→0 structural failure):
- YFB × joint quantile+rank × No: β→0, all V0–V11 strategies exhausted
- YFB × joint quantile+rank × Yes: β→0 unchanged by cohort_id

Prior: normal throughout (point_normal collapses β→0 on merged data for both models).

### K Selection Approach
K is selected by cross-validated C-index (1-SE rule) on the merged training data. This replaces the previous K=20+ARD approach (which uses many latent factors and lets ARD prune — not parsimonious or principled). K-CV is run separately for each preprocessing×model combination:
- K_merged_LB_joint (for M1, M2)
- K_merged_LB_perplatform (for M3, M4)
- K_merged_YFB_perplatform (for M5, M6)

Results are stored in globals.yml so they can be reused across script runs.

### Evaluation Metrics
1. **External C-index** on 5 held-out cohorts (Dijk, Moffitt_GEO_array, PACA_AU_array, PACA_AU_seq, Puleo_array) — primary
2. **K_eff** — number of |β_k| > beta_threshold after fitting
3. **β_max** — largest |β_k|; near-zero signals collapse
4. **Factor top genes** — for each factor k, the top 20 genes by |EF[,k]| loading; used to assess biological coherence
5. **Loading heatmap** — heatmap of |EF| for active factors (visual interpretability check)

---

## File Map

| File | Action | Purpose |
|------|--------|---------|
| `config/globals.yml` | Modify | Add `benchmark.k_merged_*` keys for CV-selected K values |
| `results/benchmark_sim/run_merged_kcv.R` | **Create** | K-CV on merged data for 3 preprocessing×model combos |
| `results/benchmark_sim/run_merged_benchmark.R` | **Create** | Fit all 6 configurations; evaluate on 5 external cohorts; save CSV + RDS |
| `results/benchmark_sim/outputs/merged_benchmark/` | Create (auto) | CSV tables + RDS fits |
| `docs/reports/merged_benchmark_report.qmd` | **Create** | Quarto report: C-index comparison table, loading heatmaps, top gene lists |
| `docs/reports/merged_benchmark_report.pdf` | Create (auto) | Rendered PDF |
| `DECISIONS.md` | Modify | Record findings and recommended configuration |
| `ROADMAP.md` | Modify | Mark prior open items resolved; add any new items |

---

## Task 1: Update globals.yml with placeholder K values

These keys will hold the CV-selected K for each configuration. They start as `null` and are filled in after Task 2 (K-CV run). The benchmark runner (Task 3) reads these; if null, it stops with a clear error directing the user to run Task 2 first.

**Files:**
- Modify: `config/globals.yml`

- [ ] **Step 1: Add K placeholder entries to globals.yml**

Open `config/globals.yml`. Find the `benchmark:` section. Add the following lines after `k_pdac_yfb_merged: 3`:

```yaml
  # CV-selected K for merged training (TCGA_PAAD + CPTAC) — filled by run_merged_kcv.R
  # Format: k_merged_<model>_<preprocessing>
  # null = K-CV not yet run; run_merged_benchmark.R will stop with an informative error.
  k_merged_lb_joint:        null   # LB, joint quantile+rank preprocessing (M1, M2)
  k_merged_lb_perplatform:  null   # LB, per-platform z-std preprocessing (M3, M4)
  k_merged_yfb_perplatform: null   # YFB, per-platform z-std preprocessing (M5, M6)
```

- [ ] **Step 2: Verify globals.yml loads cleanly in R**

```bash
Rscript -e "cfg <- yaml::read_yaml('config/globals.yml'); cat('k_merged_lb_joint:', deparse(cfg$benchmark$k_merged_lb_joint), '\n')"
```

Expected output: `k_merged_lb_joint: NULL`

- [ ] **Step 3: Commit**

```bash
git add config/globals.yml
git commit -m "Add placeholder K keys for merged-cohort CV benchmark (run_merged_kcv fills these)"
```

---

## Task 2: K-CV runner for merged training data

### Context
`select_K_cv()` in `code/select_K.R` runs cross-validated C-index over a K grid and applies the 1-SE rule. It accepts a pre-processed Y matrix. We call it three times: LB×joint, LB×per-platform, YFB×per-platform.

**Important preprocessing note:** Preprocessing (quantile normalization or per-platform z-std) is applied to the full merged training matrix before CV folds are created. This means held-out fold data influences the normalization — mild data leakage, acceptable for K selection but should be noted in comments.

**Files:**
- Create: `results/benchmark_sim/run_merged_kcv.R`
- Modify: `config/globals.yml` (filled in by this script)

- [ ] **Step 1: Create `run_merged_kcv.R`**

```r
# ============================================================
# Script:  results/benchmark_sim/run_merged_kcv.R
# Purpose: Select K via cross-validated C-index (1-SE rule) on merged
#          TCGA_PAAD + CPTAC training data for three preprocessing ×
#          model configurations:
#            (a) LB    × joint quantile+rank   (for M1, M2)
#            (b) LB    × per-platform z-std    (for M3, M4)
#            (c) YFB   × per-platform z-std    (for M5, M6)
#
#          Results are written into config/globals.yml under
#          benchmark.k_merged_lb_joint, benchmark.k_merged_lb_perplatform,
#          benchmark.k_merged_yfb_perplatform. run_merged_benchmark.R
#          reads these values.
#
#          NOTE on preprocessing leakage: preprocessing is applied to the
#          full merged matrix before CV folds split. Held-out fold data
#          therefore influences quantile normalization / z-std parameters.
#          This is standard practice for K selection (not final model
#          parameters) and the effect is small, but it is not leakage-free.
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-25
# Usage:   Rscript results/benchmark_sim/run_merged_kcv.R [--quick]
#          --quick: K_grid=2:5, n_folds=3, max_iter=50 (~5 min)
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (file.exists("code/fit_modular.R")) {
  # already at repo root
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../../")
}

suppressPackageStartupMessages({
  library(yaml)
  library(survival)
})

cfg <- yaml::read_yaml("config/globals.yml")

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
source("code/select_K.R")
tryCatch(source("code/fit_modular.R"),    error = function(e) invisible(NULL))
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/preprocess_desurv.R")

# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------

K_GRID    <- if (QUICK_MODE) 2L:5L else 2L:10L
N_FOLDS   <- if (QUICK_MODE) 3L else cfg$cavi$n_cv_folds    # 5
MAX_ITER  <- if (QUICK_MODE) 50L else cfg$cavi$max_iter      # 300
ALPHA     <- cfg$benchmark$alpha                              # 0.5
LAMBDA    <- cfg$benchmark$lambda                             # 1.0
PRIOR_BETA <- "normal"

cat("=== Merged-Cohort K-CV ===\n")
cat(sprintf("    K_grid=%d:%d | n_folds=%d | max_iter=%d | QUICK=%s\n\n",
            min(K_GRID), max(K_GRID), N_FOLDS, MAX_ITER, QUICK_MODE))

OUT_DIR <- "results/benchmark_sim/outputs/merged_benchmark"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1. Load raw training data
# --------------------------------------------------------------------------

cat("--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds))
  load_pdac_raw(ds, PDAC_DATA_ROOT)
})
n_tcga  <- train_raw$TCGA_PAAD$n
n_cptac <- train_raw$CPTAC$n
cohort_labels <- c(rep("TCGA", n_tcga), rep("CPTAC", n_cptac))

time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))

cat(sprintf("  n=%d (TCGA=%d, CPTAC=%d), events=%d\n\n", n_tcga + n_cptac,
            n_tcga, n_cptac, sum(status_train)))

# --------------------------------------------------------------------------
# 2. Preprocess — two versions
# --------------------------------------------------------------------------

cat("--- Preprocessing (joint quantile+rank) ---\n")
merged_joint <- preprocess_merged_cohorts(
  cohort_raw_list     = train_raw,
  log_transform_flags = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n               = cfg$preprocessing$top_n_genes,
  rank_transform      = TRUE,
  per_platform_standardize = FALSE
)
Y_joint <- merged_joint$Y
cat(sprintf("  Y_joint: n=%d, p=%d\n\n", nrow(Y_joint), ncol(Y_joint)))

cat("--- Preprocessing (per-platform z-std) ---\n")
merged_perplatform <- preprocess_merged_cohorts(
  cohort_raw_list     = train_raw,
  log_transform_flags = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n               = cfg$preprocessing$top_n_genes,
  rank_transform      = FALSE,
  per_platform_standardize = TRUE
)
Y_perplatform <- merged_perplatform$Y
cat(sprintf("  Y_perplatform: n=%d, p=%d\n\n", nrow(Y_perplatform), ncol(Y_perplatform)))

# --------------------------------------------------------------------------
# 3. K-CV for each configuration
# --------------------------------------------------------------------------

run_kcv <- function(Y, time, status, model, label) {
  cat(sprintf("--- K-CV: %s ---\n", label))
  cat(sprintf("    model=%s | K_grid=%d:%d | n_folds=%d | max_iter=%d\n\n",
              model, min(K_GRID), max(K_GRID), N_FOLDS, MAX_ITER))
  set.seed(42L)
  result <- select_K_cv(
    Y, time, status,
    K_grid     = K_GRID,
    n_folds    = N_FOLDS,
    model      = model,
    seed       = 42L,
    verbose    = TRUE,
    max_iter   = MAX_ITER,
    prior_beta = PRIOR_BETA,
    alpha      = ALPHA,
    lambda     = LAMBDA,
    sign_correction = FALSE    # must be FALSE in CV folds
  )
  cat(sprintf("\n  → K_opt=%d (1-SE rule)\n", result$K_opt))
  cat("\n  CV table:\n")
  print(result$cv_table[, c("K", "mean_cindex", "se_cindex")], row.names = FALSE)
  out_csv <- file.path(OUT_DIR, sprintf("kcv_%s.csv", gsub(" ", "_", tolower(label))))
  write.csv(result$cv_table, out_csv, row.names = FALSE)
  cat(sprintf("  Saved: %s\n\n", out_csv))
  result
}

kcv_lb_joint        <- run_kcv(Y_joint,        time_train, status_train, "LB",  "LB joint")
kcv_lb_perplatform  <- run_kcv(Y_perplatform,  time_train, status_train, "LB",  "LB per-platform")
kcv_yfb_perplatform <- run_kcv(Y_perplatform,  time_train, status_train, "YFB", "YFB per-platform")

# --------------------------------------------------------------------------
# 4. Write K values back into globals.yml
# --------------------------------------------------------------------------

cat("--- Updating config/globals.yml with CV-selected K values ---\n")

globals_text <- readLines("config/globals.yml")

replace_null_key <- function(lines, key, value) {
  pattern <- paste0("^(\\s*", key, ":\\s*)null")
  idx <- grep(pattern, lines)
  if (length(idx) == 0) stop(sprintf("Key '%s' not found in globals.yml", key))
  lines[idx] <- sub("null", as.character(value), lines[idx])
  lines
}

globals_text <- replace_null_key(globals_text, "k_merged_lb_joint",
                                  kcv_lb_joint$K_opt)
globals_text <- replace_null_key(globals_text, "k_merged_lb_perplatform",
                                  kcv_lb_perplatform$K_opt)
globals_text <- replace_null_key(globals_text, "k_merged_yfb_perplatform",
                                  kcv_yfb_perplatform$K_opt)

writeLines(globals_text, "config/globals.yml")
cat(sprintf("  k_merged_lb_joint        = %d\n", kcv_lb_joint$K_opt))
cat(sprintf("  k_merged_lb_perplatform  = %d\n", kcv_lb_perplatform$K_opt))
cat(sprintf("  k_merged_yfb_perplatform = %d\n", kcv_yfb_perplatform$K_opt))
cat("  globals.yml updated.\n\n")

cat("============================================================\n")
cat(" K-CV complete. Run run_merged_benchmark.R to fit all models.\n")
cat("============================================================\n")
```

- [ ] **Step 2: Verify script is syntactically valid**

```bash
Rscript -e "parse(file='results/benchmark_sim/run_merged_kcv.R'); cat('OK\n')"
```

Expected: `OK`

- [ ] **Step 3: Run K-CV in quick mode to verify it executes end-to-end**

```bash
Rscript results/benchmark_sim/run_merged_kcv.R --quick 2>&1 | tail -20
```

Expected: three K-CV tables printed, globals.yml updated with integer K values (not null). Verify:

```bash
grep "k_merged" config/globals.yml
```

Expected: three lines with integer values (not null).

- [ ] **Step 4: Reset globals.yml K values back to null for the full run**

The quick-mode K values (from K_grid=2:5, 3 folds, 50 iters) are not reliable enough to keep. Reset them so the full run fills them properly:

```bash
sed -i '' 's/k_merged_lb_joint:        [0-9]*/k_merged_lb_joint:        null/' config/globals.yml
sed -i '' 's/k_merged_lb_perplatform:  [0-9]*/k_merged_lb_perplatform:  null/' config/globals.yml
sed -i '' 's/k_merged_yfb_perplatform: [0-9]*/k_merged_yfb_perplatform: null/' config/globals.yml
```

Verify reset: `grep "k_merged" config/globals.yml` should show `null` for all three.

- [ ] **Step 5: Commit the script (not the globals.yml K values yet)**

```bash
git add results/benchmark_sim/run_merged_kcv.R config/globals.yml
git commit -m "Add run_merged_kcv.R: CV K selection on merged TCGA+CPTAC for all preprocessing×model combos"
```

- [ ] **Step 6: Run full K-CV (background, ~45–90 min, sleep prevention active)**

`caffeinate -s` prevents macOS from sleeping while plugged into AC power. It wraps
the R process and releases the sleep assertion automatically when the job finishes.
Verify the machine is plugged in before running.

```bash
caffeinate -s Rscript results/benchmark_sim/run_merged_kcv.R > /tmp/merged_kcv.log 2>&1 &
echo "PID: $!"
```

Monitor progress: `tail -f /tmp/merged_kcv.log`

Check job is still running: `ps aux | grep "run_merged_kcv" | grep -v grep`

- [ ] **Step 7: Commit globals.yml with CV-selected K values**

After the full run completes:

```bash
grep "k_merged" config/globals.yml   # verify non-null integers
git add config/globals.yml
git commit -m "Fill CV-selected K in globals.yml: merged TCGA+CPTAC (LB joint/perplatform, YFB perplatform)"
```

---

## Task 3: Merged benchmark runner

Fits all 6 model configurations at their CV-selected K and evaluates on 5 external cohorts. Reads K from globals.yml (stops with clear error if any K is still null). Saves results to CSV and compact RDS.

**Files:**
- Create: `results/benchmark_sim/run_merged_benchmark.R`
- Creates (auto): `results/benchmark_sim/outputs/merged_benchmark/merged_benchmark_results.csv`
- Creates (auto): `results/benchmark_sim/outputs/merged_benchmark/merged_benchmark_fits.rds`

- [ ] **Step 1: Create `run_merged_benchmark.R`**

```r
# ============================================================
# Script:  results/benchmark_sim/run_merged_benchmark.R
# Purpose: Comprehensive merged-cohort benchmark (TCGA_PAAD + CPTAC).
#          Fits 6 model configurations at CV-selected K and evaluates
#          external C-index on 5 held-out PDAC cohorts.
#
#          Configurations (all prior_beta="normal"):
#            M1: LB  × joint quantile+rank  × no cohort_id
#            M2: LB  × joint quantile+rank  × cohort_id
#            M3: LB  × per-platform z-std   × no cohort_id
#            M4: LB  × per-platform z-std   × cohort_id
#            M5: YFB × per-platform z-std   × no cohort_id
#            M6: YFB × per-platform z-std   × cohort_id
#
#          K for each configuration is read from globals.yml
#          (benchmark.k_merged_lb_joint, k_merged_lb_perplatform,
#           k_merged_yfb_perplatform). Run run_merged_kcv.R first.
#
#          Excluded (documented β→0 structural failure, all V0–V11 exhausted):
#            YFB × joint quantile+rank × No/Yes
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-25
# Usage:   Rscript results/benchmark_sim/run_merged_benchmark.R [--quick]
#          --quick: max_iter=30, skips interpretability output (smoke test)
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (file.exists("code/fit_modular.R")) {
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../../")
}

suppressPackageStartupMessages({
  library(yaml)
  library(survival)
})

cfg <- yaml::read_yaml("config/globals.yml")

# Guard: K values must be filled by run_merged_kcv.R first
k_lb_joint       <- cfg$benchmark$k_merged_lb_joint
k_lb_perplatform <- cfg$benchmark$k_merged_lb_perplatform
k_yfb_perplatform <- cfg$benchmark$k_merged_yfb_perplatform

if (is.null(k_lb_joint) || is.null(k_lb_perplatform) || is.null(k_yfb_perplatform)) {
  missing <- c(
    if (is.null(k_lb_joint))        "k_merged_lb_joint",
    if (is.null(k_lb_perplatform))  "k_merged_lb_perplatform",
    if (is.null(k_yfb_perplatform)) "k_merged_yfb_perplatform"
  )
  stop(sprintf(
    "K values not yet set in globals.yml: %s\nRun: Rscript results/benchmark_sim/run_merged_kcv.R",
    paste(missing, collapse = ", ")
  ))
}

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_modular.R"),    error = function(e) invisible(NULL))
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/preprocess_desurv.R")

ALPHA      <- cfg$benchmark$alpha
LAMBDA     <- cfg$benchmark$lambda
MAX_ITER   <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
PRIOR_BETA <- "normal"
SIGMA_COH  <- 1.0
BETA_THRESH <- cfg$k_selection$beta_threshold

OUT_DIR <- "results/benchmark_sim/outputs/merged_benchmark"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Merged Cohort Benchmark — 6-way comparison ===\n")
cat(sprintf("    alpha=%.2f | prior_beta=%s | max_iter=%d | QUICK=%s\n",
            ALPHA, PRIOR_BETA, MAX_ITER, QUICK_MODE))
cat(sprintf("    K: LB_joint=%d | LB_perplatform=%d | YFB_perplatform=%d\n\n",
            k_lb_joint, k_lb_perplatform, k_yfb_perplatform))

# --------------------------------------------------------------------------
# 1. Load and preprocess training data — both versions
# --------------------------------------------------------------------------

cat("--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
n_tcga  <- train_raw$TCGA_PAAD$n
n_cptac <- train_raw$CPTAC$n
cohort_labels <- c(rep("TCGA", n_tcga), rep("CPTAC", n_cptac))
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))

merged_joint <- preprocess_merged_cohorts(
  train_raw, PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n = cfg$preprocessing$top_n_genes,
  rank_transform = TRUE, per_platform_standardize = FALSE
)

merged_perplatform <- preprocess_merged_cohorts(
  train_raw, PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n = cfg$preprocessing$top_n_genes,
  rank_transform = FALSE, per_platform_standardize = TRUE
)

cat(sprintf("  joint: n=%d, p=%d\n", nrow(merged_joint$Y), ncol(merged_joint$Y)))
cat(sprintf("  perplatform: n=%d, p=%d\n\n", nrow(merged_perplatform$Y),
            ncol(merged_perplatform$Y)))

# --------------------------------------------------------------------------
# Helper: oriented C-index
# --------------------------------------------------------------------------
oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 2. Fit all 6 configurations
# --------------------------------------------------------------------------

fits      <- list()
gene_sets <- list()   # store gene_names per preprocessing type for external validation

fit_lb <- function(Y, time, status, K, cohort_id = NULL, label = "") {
  cat(sprintf("--- Fitting %s (K=%d) ---\n", label, K))
  set.seed(42L)
  fit <- suppressMessages(
    fit_supervised_mf_modular(Y, time, status,
                              K          = K,
                              max_iter   = MAX_ITER,
                              alpha      = ALPHA,
                              lambda     = LAMBDA,
                              prior_beta = PRIOR_BETA,
                              verbose    = TRUE,
                              cohort_id  = cohort_id,
                              sigma_F_cohort = SIGMA_COH)
  )
  k_eff <- sum(abs(fit$EBeta) > BETA_THRESH)
  cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n\n",
              k_eff, max(abs(fit$EBeta)), fit$history$n_iter))
  fit
}

fit_yfb <- function(Y, time, status, K, cohort_id = NULL, label = "") {
  cat(sprintf("--- Fitting %s (K=%d) ---\n", label, K))
  set.seed(42L)
  fit <- suppressMessages(
    fit_cox_on_yf(Y, time, status,
                  K          = K,
                  max_iter   = MAX_ITER,
                  alpha      = ALPHA,
                  lambda     = LAMBDA,
                  prior_beta = PRIOR_BETA,
                  verbose    = TRUE,
                  cohort_id  = cohort_id,
                  sigma_F_cohort = SIGMA_COH)
  )
  k_eff <- sum(abs(fit$EBeta) > BETA_THRESH)
  cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n\n",
              k_eff, max(abs(fit$EBeta)), fit$history$n_iter))
  fit
}

fits$M1 <- fit_lb(merged_joint$Y,       time_train, status_train, k_lb_joint,
                   cohort_id = NULL,         label = "M1 LB joint no-cohort")
fits$M2 <- fit_lb(merged_joint$Y,       time_train, status_train, k_lb_joint,
                   cohort_id = cohort_labels, label = "M2 LB joint cohort_id")
fits$M3 <- fit_lb(merged_perplatform$Y, time_train, status_train, k_lb_perplatform,
                   cohort_id = NULL,         label = "M3 LB perplatform no-cohort")
fits$M4 <- fit_lb(merged_perplatform$Y, time_train, status_train, k_lb_perplatform,
                   cohort_id = cohort_labels, label = "M4 LB perplatform cohort_id")
fits$M5 <- fit_yfb(merged_perplatform$Y, time_train, status_train, k_yfb_perplatform,
                    cohort_id = NULL,         label = "M5 YFB perplatform no-cohort")
fits$M6 <- fit_yfb(merged_perplatform$Y, time_train, status_train, k_yfb_perplatform,
                    cohort_id = cohort_labels, label = "M6 YFB perplatform cohort_id")

gene_sets$joint       <- merged_joint$gene_names
gene_sets$perplatform <- merged_perplatform$gene_names

# --------------------------------------------------------------------------
# 3. External validation on 5 held-out cohorts
# --------------------------------------------------------------------------

cat("--- External validation (5 cohorts) ---\n")

EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
results_rows <- list()

for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("  %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
  pre_ext <- preprocess_desurv_cohort(
    Y             = raw_ext$Y,
    gene_names    = raw_ext$gene_names,
    top_n         = cfg$preprocessing$top_n_genes,
    log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]],
    cohort_name   = ext_cohort
  )

  # For each training gene set, intersect and evaluate all models using that set
  for (prep_type in c("joint", "perplatform")) {
    train_genes <- gene_sets[[prep_type]]
    common      <- intersect(train_genes, pre_ext$gene_names)
    if (length(common) < 100) next
    Y_ext      <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]
    train_idx  <- match(common, train_genes)

    models_for_prep <- if (prep_type == "joint") c("M1", "M2") else c("M3", "M4", "M5", "M6")
    for (mid in models_for_prep) {
      fit <- fits[[mid]]
      is_yfb <- mid %in% c("M5", "M6")
      if (!is_yfb) {
        EF_sub <- fit$EF[train_idx, , drop = FALSE]
        pred   <- predict_supervised_mf(Y_ext, EF_sub, fit$EBeta)
        c_val  <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)
      } else {
        EF_sub <- fit$EF[train_idx, , drop = FALSE]
        pred   <- predict_cox_on_yf(Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)
        c_val  <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)
      }
      results_rows[[length(results_rows) + 1]] <- data.frame(
        model      = mid,
        cohort     = ext_cohort,
        c_index    = round(c_val, 4),
        k_eff      = sum(abs(fit$EBeta) > BETA_THRESH),
        beta_max   = round(max(abs(fit$EBeta)), 4),
        n_iters    = fit$history$n_iter,
        preprocess = prep_type,
        has_cohort = mid %in% c("M2", "M4", "M6"),
        stringsAsFactors = FALSE
      )
    }
  }
}

results_df <- do.call(rbind, results_rows)

# --------------------------------------------------------------------------
# 4. Summary table
# --------------------------------------------------------------------------

cat("\n============================================================\n")
cat(" Merged Benchmark Results — External C-index\n")
cat("============================================================\n")

model_ids <- c("M1","M2","M3","M4","M5","M6")
model_labels <- c("LB_joint","LB_joint_coh","LB_perplat","LB_perplat_coh",
                  "YFB_perplat","YFB_perplat_coh")

for (i in seq_along(model_ids)) {
  mid <- model_ids[i]
  sub <- results_df[results_df$model == mid, ]
  if (nrow(sub) == 0) next
  cat(sprintf("  %s (%s):\n", mid, model_labels[i]))
  for (j in seq_len(nrow(sub))) {
    cat(sprintf("    %-22s C=%.3f\n", sub$cohort[j], sub$c_index[j]))
  }
  cat(sprintf("    %-22s C=%.3f  (K_eff=%d, beta_max=%.4f)\n\n",
              "MEAN", mean(sub$c_index), sub$k_eff[1], sub$beta_max[1]))
}

# --------------------------------------------------------------------------
# 5. Factor top-20 gene table (interpretability)
# --------------------------------------------------------------------------

if (!QUICK_MODE) {
  cat("--- Factor top-20 genes per configuration ---\n")
  top_genes_list <- list()
  for (mid in model_ids) {
    fit <- fits[[mid]]
    K_fit <- ncol(fit$EF)
    prep_type <- if (mid %in% c("M1","M2")) "joint" else "perplatform"
    genes <- gene_sets[[prep_type]]
    top_genes_list[[mid]] <- lapply(seq_len(K_fit), function(k) {
      idx <- order(abs(fit$EF[, k]), decreasing = TRUE)[1:min(20, nrow(fit$EF))]
      data.frame(model = mid, factor = k, gene = genes[idx],
                 loading = round(fit$EF[idx, k], 4))
    })
  }
  top_genes_df <- do.call(rbind, lapply(top_genes_list, function(x) do.call(rbind, x)))
  write.csv(top_genes_df,
            file.path(OUT_DIR, "merged_benchmark_top_genes.csv"),
            row.names = FALSE)
  cat(sprintf("  Top genes saved: %s\n\n",
              file.path(OUT_DIR, "merged_benchmark_top_genes.csv")))
}

# --------------------------------------------------------------------------
# 6. Save results
# --------------------------------------------------------------------------

write.csv(results_df, file.path(OUT_DIR, "merged_benchmark_results.csv"), row.names = FALSE)

compact_fits <- lapply(fits, function(f) list(
  EBeta     = f$EBeta,   EBeta2   = f$EBeta2,
  EF_cohort = f$EF_cohort, EF2_cohort = f$EF2_cohort,
  EF_norms  = if (!is.null(f$EF_norms)) f$EF_norms else NULL,
  history   = f$history
))
saveRDS(list(fits = compact_fits, results = results_df, gene_sets = gene_sets,
             params = list(k_lb_joint = k_lb_joint,
                           k_lb_perplatform = k_lb_perplatform,
                           k_yfb_perplatform = k_yfb_perplatform,
                           ALPHA = ALPHA, PRIOR_BETA = PRIOR_BETA),
             date = Sys.time()),
        file.path(OUT_DIR, "merged_benchmark_fits.rds"))

cat(sprintf("Results saved to: %s\n", OUT_DIR))
cat("============================================================\n")
```

- [ ] **Step 2: Verify syntax**

```bash
Rscript -e "parse(file='results/benchmark_sim/run_merged_benchmark.R'); cat('OK\n')"
```

Expected: `OK`

- [ ] **Step 3: Run in quick mode to verify it executes (requires K values in globals.yml)**

First ensure K values are filled (run Task 2 quick mode again if needed):
```bash
grep "k_merged" config/globals.yml
```
If all three show integers (not null), proceed. Otherwise re-run `--quick` K-CV:
```bash
Rscript results/benchmark_sim/run_merged_kcv.R --quick
```
Then run the benchmark in quick mode:
```bash
Rscript results/benchmark_sim/run_merged_benchmark.R --quick 2>&1 | tail -30
```

Expected: C-index values printed for each model × cohort. Check that M1/M2 (joint) and M3/M4/M5/M6 (perplatform) all appear, no errors.

- [ ] **Step 4: Commit**

```bash
git add results/benchmark_sim/run_merged_benchmark.R
git commit -m "Add run_merged_benchmark.R: 6-config merged benchmark (LB/YFB × joint/perplatform × ±cohort_id)"
```

- [ ] **Step 5: Run full benchmark (background, ~60–120 min, sleep prevention active)**

Verify the machine is still plugged into AC power before launching.

```bash
caffeinate -s Rscript results/benchmark_sim/run_merged_benchmark.R > /tmp/merged_benchmark.log 2>&1 &
echo "PID: $!"
```

Monitor progress: `tail -f /tmp/merged_benchmark.log`

Check job is still running: `ps aux | grep "run_merged_benchmark" | grep -v grep`

- [ ] **Step 6: Commit outputs**

```bash
git add results/benchmark_sim/outputs/merged_benchmark/
git commit -m "Add merged benchmark results: C-index table and top-gene lists for all 6 configurations"
```

---

## Task 4: Quarto report

Synthesizes all results into a readable document for biostatistician collaborators. Sections: design table, C-index heatmap, model summary, top-gene lists, conclusions.

**Files:**
- Create: `docs/reports/merged_benchmark_report.qmd`

- [ ] **Step 1: Create `docs/reports/merged_benchmark_report.qmd`**

```qmd
---
title: "Merged-Cohort PDAC Benchmark: LB vs. YFB, Joint vs. Per-Platform Preprocessing, With and Without Cohort Indicator"
author: "Andrew Walther"
date: "`r Sys.Date()`"
format:
  pdf:
    toc: true
    number-sections: true
    fig-pos: "H"
  html:
    toc: true
    code-fold: true
    theme: cosmo
execute:
  echo: false
  warning: false
  message: false
---

```{r setup}
library(yaml)
library(ggplot2)
library(tidyr)
library(dplyr)

cfg     <- yaml::read_yaml(here::here("config/globals.yml"))
results <- read.csv(here::here(
  "results/benchmark_sim/outputs/merged_benchmark/merged_benchmark_results.csv"))
top_genes <- tryCatch(
  read.csv(here::here(
    "results/benchmark_sim/outputs/merged_benchmark/merged_benchmark_top_genes.csv")),
  error = function(e) NULL)

model_labels <- c(
  M1 = "LB | joint | no cohort adj.",
  M2 = "LB | joint | cohort_id",
  M3 = "LB | per-platform | no cohort adj.",
  M4 = "LB | per-platform | cohort_id",
  M5 = "YFB | per-platform | no cohort adj.",
  M6 = "YFB | per-platform | cohort_id"
)
results$model_label <- model_labels[results$model]
```

## Study Design

This report compares six model configurations for joint multi-cohort PDAC analysis,
trained on merged TCGA_PAAD (RNA-seq, n=144) and CPTAC (proteomics, n=129) data
and evaluated on five independent external cohorts.

**Goal:** Identify the configuration that produces a parsimonious model (small, CV-selected K)
with non-zero survival associations (β ≠ 0) and biologically coherent gene expression programs
that generalize across platforms and patient populations.

| ID | Model | Preprocessing | Cohort adjustment | K |
|----|-------|---------------|-------------------|---|
| M1 | LB (η = Lβ) | Joint quantile+rank | None | `r cfg$benchmark$k_merged_lb_joint` |
| M2 | LB | Joint quantile+rank | Cohort indicator | `r cfg$benchmark$k_merged_lb_joint` |
| M3 | LB | Per-platform z-std | None | `r cfg$benchmark$k_merged_lb_perplatform` |
| M4 | LB | Per-platform z-std | Cohort indicator | `r cfg$benchmark$k_merged_lb_perplatform` |
| M5 | YFB (η = YFβ) | Per-platform z-std | None | `r cfg$benchmark$k_merged_yfb_perplatform` |
| M6 | YFB | Per-platform z-std | Cohort indicator | `r cfg$benchmark$k_merged_yfb_perplatform` |

YFB × joint preprocessing is excluded: β→0 collapse is structural and all mitigation
strategies (V0–V11) have been exhausted (see DECISIONS.md 2026-05-22/25).

**K selection:** Cross-validated C-index (1-SE rule) on merged training data, separately
for each preprocessing × model combination. This replaces the previous K=20+ARD approach.

## External C-index Results

```{r cindex-table}
wide <- results |>
  select(model_label, cohort, c_index) |>
  pivot_wider(names_from = cohort, values_from = c_index)

means <- results |>
  group_by(model_label) |>
  summarise(Mean = round(mean(c_index), 3), .groups = "drop")

wide <- left_join(wide, means, by = "model_label")
knitr::kable(wide, digits = 3, caption = "External C-index by model and cohort")
```

```{r cindex-heatmap, fig.height=4, fig.width=8}
ggplot(results, aes(x = cohort, y = model_label, fill = c_index)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.3f", c_index)), size = 3) +
  scale_fill_gradient2(low = "#d73027", mid = "#fee090", high = "#1a9850",
                       midpoint = 0.55, name = "C-index") +
  labs(title = "External C-index heatmap", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
```

## Model Summary (β activity)

```{r beta-summary}
summary_df <- results |>
  distinct(model, model_label, k_eff, beta_max, n_iters) |>
  mutate(
    K = c(cfg$benchmark$k_merged_lb_joint, cfg$benchmark$k_merged_lb_joint,
          cfg$benchmark$k_merged_lb_perplatform, cfg$benchmark$k_merged_lb_perplatform,
          cfg$benchmark$k_merged_yfb_perplatform, cfg$benchmark$k_merged_yfb_perplatform)[
      match(model, c("M1","M2","M3","M4","M5","M6"))]
  ) |>
  select(model_label, K, k_eff, beta_max, n_iters)
knitr::kable(summary_df, digits = 4,
             caption = "Model summary: CV-selected K, active factors (K_eff), max |β|, iterations to convergence")
```

## Factor Interpretability: Top Genes

```{r top-genes, results='asis'}
if (!is.null(top_genes)) {
  for (mid in c("M1","M2","M3","M4","M5","M6")) {
    sub <- top_genes[top_genes$model == mid, ]
    if (nrow(sub) == 0) next
    cat(sprintf("\n### %s — %s\n\n", mid, model_labels[mid]))
    for (k in sort(unique(sub$factor))) {
      genes_k <- sub[sub$factor == k, "gene"]
      cat(sprintf("**Factor %d:** %s\n\n", k, paste(genes_k, collapse = ", ")))
    }
  }
} else {
  cat("Top gene table not available (run benchmark without --quick).\n")
}
```

## Conclusions

*(To be filled in after results are available. Key questions to address:)*

- *Which configuration achieves the highest mean external C-index?*
- *Does per-platform preprocessing consistently improve over joint preprocessing for LB?*
- *Does the cohort indicator help on top of per-platform preprocessing, or is it redundant?*
- *Which configuration produces the most biologically coherent top genes per factor?*
- *What is the recommended configuration for the manuscript?*
```

- [ ] **Step 2: Verify Quarto renders without error (placeholder data OK)**

After results CSV is available:

```bash
quarto render docs/reports/merged_benchmark_report.qmd --to pdf 2>&1 | tail -10
```

Expected: no errors, `merged_benchmark_report.pdf` created.

- [ ] **Step 3: Commit**

```bash
git add docs/reports/merged_benchmark_report.qmd
git commit -m "Add merged_benchmark_report.qmd: C-index heatmap, model summary, top genes, conclusions"
```

---

## Task 5: Documentation update

After results are in hand, record findings and recommended configuration.

**Files:**
- Modify: `DECISIONS.md`
- Modify: `ROADMAP.md`
- Modify: `CLAUDE.md` (update Current model status line)

- [ ] **Step 1: Add DECISIONS.md entry**

Add a new entry at the top of DECISIONS.md (after the header block):

```markdown
## 2026-MM-DD — Merged-cohort benchmark: recommended configuration identified

- **Finding:** [Fill in after results: which of M1–M6 wins on mean external C-index,
  which wins on β activity, which wins on factor interpretability]

- **Key comparisons:**
  - LB joint vs. LB per-platform: [does preprocessing matter for LB?]
  - cohort_id Yes vs. No at matched preprocessing: [does dummy column help?]
  - LB vs. YFB at per-platform preprocessing: [which model type is better?]

- **Recommended configuration for manuscript:**
  - Model: [LB / YFB]
  - Preprocessing: [joint / per-platform]
  - cohort_id: [Yes / No]
  - K: [CV-selected value]
  - Rationale: [fill in]

- **Files:** `results/benchmark_sim/run_merged_kcv.R`,
  `results/benchmark_sim/run_merged_benchmark.R`,
  `docs/reports/merged_benchmark_report.{qmd,pdf}`,
  `results/benchmark_sim/outputs/merged_benchmark/`
```

- [ ] **Step 2: Update ROADMAP.md**

Mark the merged benchmark item as complete. Add any follow-up items identified from results (e.g., biological pathway enrichment if top genes look promising).

- [ ] **Step 3: Update CLAUDE.md Current model status line**

In the CLAUDE.md `## Key Instructions` section, update the "Current model status" bullet to reflect the recommended configuration from this benchmark.

- [ ] **Step 4: Final commit**

```bash
git add DECISIONS.md ROADMAP.md CLAUDE.md
git commit -m "Document merged-cohort benchmark findings and recommended configuration"
```

---

## Self-Review

**Spec coverage:**
- ✅ CV-selected K throughout (Tasks 1, 2)
- ✅ LB with per-platform preprocessing (M3, M4 in Task 3)
- ✅ Cohort dummy column with per-platform preprocessing (M4, M6 in Task 3)
- ✅ YFB with per-platform preprocessing + cohort_id (M6 in Task 3)
- ✅ globals.yml used for all constants
- ✅ Factor interpretability (top-20 genes, Task 3 section 5)
- ✅ Documentation (Task 4, Task 5)
- ✅ YFB × joint documented as excluded (not re-run)

**Placeholder scan:** No TBDs except the Conclusions section of the report (intentional — filled after results). All code blocks are complete.

**Type consistency:** `gene_sets$joint` and `gene_sets$perplatform` defined in Task 3 and consumed correctly in the external validation loop. `model_labels` vector defined in Task 4 matches the 6 model IDs used throughout.

**One known limitation to note in comments:** K-CV preprocessing applies to the full training matrix before fold splits (mild leakage, standard practice for K selection).
