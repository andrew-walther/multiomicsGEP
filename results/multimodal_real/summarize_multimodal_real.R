# ============================================================
# Script: summarize_multimodal_real.R
# Purpose: Collect every cached TCGA -> ICGC multimodal fit into a K_init
#          sweep table and figures, and compare the joint expression
#          programs with the single-modality YFB programs and with Yusha's
#          unsupervised TCGA flash fit.
# Author: Andrew Walther
# Created: 2026-10-01
# Dependencies: survival, ggplot2, yaml; code/ multimodal and single-modality
#               YFB code; results/multimodal_real/outputs/fits/*.rds
# Run: Rscript results/multimodal_real/summarize_multimodal_real.R
#      (safe to rerun while fits are still running; it uses whatever is cached)
# ============================================================

suppressPackageStartupMessages({ library(survival); library(ggplot2) })
for (f in c("code/multimodal_yfb_helpers.R", "code/preprocess_multimodal_yfb.R",
            "code/multimodal_yfb_updates.R", "code/fit_multimodal_yfb.R",
            "code/predict_multimodal_yfb.R", "code/simulate_multimodal_yfb.R",
            "code/update_beta.R", "code/update_L.R", "code/update_F.R",
            "code/update_tau.R", "code/compute_elbo.R")) source(f)
suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL)))
source("code/predict_cox_on_yf.R"); source("code/run_multimodal_yfb_simulation.R")
source("code/concordance_ci.R"); source("code/pathway_enrichment.R")

out_dir <- "results/multimodal_real/outputs"
d <- readRDS("data/processed_multiomics_tcga_icgc.rds")
files <- list.files(file.path(out_dir, "fits"), pattern = "^fits_K.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("No cached fits in ", file.path(out_dir, "fits"))

#' Parse K_init and configuration from a cache file name
#' fits_K07.rds -> original; fits_K07_point_laplace-point_laplace_pruned.rds -> pruned
parse_name <- function(path) {
  b <- sub("\\.rds$", "", basename(path))
  list(K_init = as.integer(sub("^fits_K(\\d+).*", "\\1", b)),
       config = if (grepl("_pruned$", b)) "pruned" else if (b == sub("_.*", "", b) ||
                  grepl("^fits_K\\d+$", b)) "original" else sub("^fits_K\\d+_", "", b))
}

#' Inputs for a configuration: centered/asin variants transform the data as
#' the runner did; raw and pruned use the cached data as-is (the pruned fit's
#' intercept handles centering internally)
cohorts_for <- function(config) {
  x <- d[c("training", "validation_primary", "validation_all")]
  if (grepl("centered", config)) {
    if (grepl("asin", config)) for (nm in names(x)) x[[nm]]$Y$methylation <- asin(2 * x[[nm]]$Y$methylation - 1)
    mu <- lapply(x$training$Y, colMeans)
    for (nm in names(x)) x[[nm]]$Y <- Map(function(y, m) sweep(y, 2, m), x[[nm]]$Y, mu)
  }
  names(x) <- c("tcga_train", "icgc_primary", "icgc_all")
  x
}

score <- function(risk_fn, cohorts) {
  do.call(rbind, lapply(names(cohorts), function(nm) {
    x <- cohorts[[nm]]; r <- risk_fn(x$Y)
    # Frozen orientation: every score here comes from a training-data Cox fit,
    # so it is already a risk score (flip = FALSE). flip = NULL would re-pick
    # the direction from the validation outcomes (circular).
    ci <- bootstrap_concordance_ci(r, x$time, x$event, B = 1000, seed = 1, flip = FALSE)
    data.frame(cohort = nm, c_index = multimodal_yfb_survival_metrics(r, x$time, x$event)$c_index,
               ci_lower = ci$lower, ci_upper = ci$upper)
  }))
}

# Sweep table ----
rows <- list(); joint_fits <- list()
for (path in files) {
  info <- parse_name(path); fits <- readRDS(path); cohorts <- cohorts_for(info$config)
  j <- fits$joint; dj <- j$diagnostics
  joint_fits[[basename(path)]] <- c(info, list(fit = j))
  cls <- dj$factor_class
  rows[[length(rows) + 1]] <- cbind(
    data.frame(config = info$config, K_init = info$K_init, method = "joint_multimodal_yfb",
               K_final = ncol(j$EL), n_variance_only = sum(cls == "reconstruction_only"),
               n_survival_only = sum(cls == "survival_only"), n_both = sum(cls == "both"),
               n_survival_active = sum(dj$survival_active),
               converged = dj$converged, elbo = as.numeric(dj$elbo %||% NA)),
    score(function(Y) predict_multimodal_yfb(j, Y)$risk_scores, cohorts))
  base <- list(
    expression_only_yfb = function(Y) predict_multimodal_yfb_single_modality(fits$expr_only, Y)$risk_scores,
    two_step_ebmf_cox = function(Y) predict_multimodal_yfb_ebmf_cox(fits$two_step, Y))
  for (m in names(base)) {
    rows[[length(rows) + 1]] <- cbind(
      data.frame(config = info$config, K_init = info$K_init, method = m, K_final = NA,
                 n_variance_only = NA, n_survival_only = NA, n_both = NA,
                 n_survival_active = NA, converged = NA, elbo = NA),
      score(base[[m]], cohorts))
  }
}
sweep_tab <- do.call(rbind, rows)
sweep_tab <- sweep_tab[order(sweep_tab$config, sweep_tab$method, sweep_tab$K_init, sweep_tab$cohort), ]
write.csv(sweep_tab, file.path(out_dir, "real_sweep_summary.csv"), row.names = FALSE)
print(sweep_tab[sweep_tab$method == "joint_multimodal_yfb" & sweep_tab$cohort == "icgc_primary",
                c("config", "K_init", "K_final", "n_variance_only", "n_survival_only", "n_both",
                  "converged", "c_index", "ci_lower", "ci_upper")], digits = 3, row.names = FALSE)

# Figures ----
theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                                legend.position = "bottom"))
joint <- sweep_tab[sweep_tab$method == "joint_multimodal_yfb" & sweep_tab$cohort == "tcga_train", ]
counts <- rbind(
  data.frame(joint[c("config", "K_init")], class = "Variance only", n = joint$n_variance_only),
  data.frame(joint[c("config", "K_init")], class = "Survival only", n = joint$n_survival_only),
  data.frame(joint[c("config", "K_init")], class = "Variance and survival", n = joint$n_both))
counts$class <- factor(counts$class, levels = c("Variance only", "Variance and survival", "Survival only"))
p1 <- ggplot(counts, aes(factor(K_init), n, fill = class)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.5) +
  facet_wrap(~config, nrow = 1) +
  scale_fill_manual(values = c("Variance only" = "#2a78d6", "Variance and survival" = "#eb6834",
                               "Survival only" = "#1baf7a"), name = NULL) +
  labs(x = expression(K[init]), y = "Retained factors",
       title = "TCGA multimodal fits: retained factors by role")
ggsave(file.path(out_dir, "real_sweep_factor_roles.png"), p1, width = 9, height = 4, dpi = 150)

cvals <- sweep_tab[sweep_tab$cohort != "tcga_train", ]
cvals$cohort <- factor(cvals$cohort, levels = c("icgc_primary", "icgc_all"),
                       labels = c("ICGC primary PDAC (n = 50)", "ICGC all donors (n = 67)"))
p2 <- ggplot(cvals, aes(K_init, c_index, colour = method, shape = method)) +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey60") +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), position = position_dodge(width = 0.6)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.6)) +
  facet_grid(cohort ~ config) +
  scale_colour_manual(values = c(joint_multimodal_yfb = "#2a78d6", expression_only_yfb = "#eb6834",
                                 two_step_ebmf_cox = "#1baf7a"), name = NULL) +
  scale_shape_manual(values = c(joint_multimodal_yfb = 16, expression_only_yfb = 17,
                                two_step_ebmf_cox = 15), name = NULL) +
  scale_x_continuous(breaks = sort(unique(cvals$K_init))) +
  labs(x = expression(K[init]), y = "External C-index (95% bootstrap CI)",
       title = "Trained on TCGA, validated on ICGC")
ggsave(file.path(out_dir, "real_sweep_cindex.png"), p2, width = 9, height = 6, dpi = 150)

# Factor comparison ----
# Joint expression loadings vs single-modality YFB programs (D4 fit) and
# Yusha's unsupervised TCGA flash fit, on shared genes (Pearson |cor|)
d4 <- load_d4_weights()
flash <- readRDS("data/multiomicsGEP_code/tcga/tcga_flash_K14.rds")$F_pm
ref <- list(d4 = d4$EF[, c(3, 7, 5, 6), drop = FALSE], flash = flash)
colnames(ref$d4) <- c("P3 protective", "P7 adverse", "P5 genomics", "P6 genomics")
colnames(ref$flash) <- paste0("flash k", seq_len(ncol(flash)))
comp <- list()
for (nm in names(joint_fits)) {
  jf <- joint_fits[[nm]]; if (jf$config != "pruned") next
  E <- jf$fit$EF$expression
  rownames(E) <- jf$fit$training_spec$feature_names$expression
  cls <- jf$fit$diagnostics$factor_class
  colnames(E) <- paste0("k", seq_len(ncol(E)), " (", sub("reconstruction_only", "variance",
                        sub("survival_only", "survival", cls)), ", beta=",
                        sprintf("%.2f", jf$fit$EBeta), ")")
  for (r in names(ref)) {
    g <- intersect(rownames(E), rownames(ref[[r]]))
    if (length(g) < 50) next
    cm <- abs(stats::cor(E[g, , drop = FALSE], ref[[r]][g, , drop = FALSE]))
    comp[[paste(nm, r)]] <- data.frame(fit = nm, K_init = jf$K_init, reference = r,
      joint_factor = rep(rownames(cm), ncol(cm)), ref_factor = rep(colnames(cm), each = nrow(cm)),
      abs_cor = as.vector(cm), n_genes = length(g))
  }
}
if (length(comp)) {
  comp <- do.call(rbind, comp)
  write.csv(comp, file.path(out_dir, "real_factor_comparison.csv"), row.names = FALSE)
  best <- do.call(rbind, lapply(split(comp, list(comp$fit, comp$joint_factor, comp$reference), drop = TRUE),
                                function(x) x[which.max(x$abs_cor), ]))
  print(best[order(best$K_init, best$reference, best$joint_factor),
             c("K_init", "reference", "joint_factor", "ref_factor", "abs_cor", "n_genes")],
        digits = 2, row.names = FALSE)
}
cat("Wrote real_sweep_summary.csv, real_sweep_*.png, real_factor_comparison.csv to", out_dir, "\n")
