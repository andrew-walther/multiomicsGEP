# ============================================================
# Script: make_real_data_workflows.R
# Purpose: Draw the evaluated real-data training and validation workflows.
# Created: 2026-10-16
# Dependencies: base R (graphics, grDevices)
# Run: Rscript paper/writeup/figures/make_real_data_workflows.R
# ============================================================
# Cohort counts are the audited historical analysis populations, not values
# inferred from a consortium's total enrollment. Sources: load_multiomics_data.R,
# benchmark_helpers.R, and _process/Accuracy_Review_10_16_26.md. This script
# reads no patient data, changes no configuration, and fits/scores no model.
# Displayed feature counts describe the frozen benchmark, not future defaults.

# Drawing helpers ----

#' Draw a workflow stage with a shape indicating its role
#' @param x1,y1,x2,y2 Numeric box edges in normalized figure coordinates.
#' @param label Character stage description, with newline-separated lines.
#' @param cex Numeric text scale relative to the 12-point device font.
#' @param shape Character role: process, data, model, result, or note.
#' @return Invisibly NULL; draws the role-specific shape and centered text.
#' @examples
#' # After plot.new(): stage_box(.05, .8, .45, .95, "Training\n144 patients")
stage_box <- function(x1, y1, x2, y2, label, cex = 0.79,
                      shape = "process") {
  stopifnot(shape %in% c("process", "data", "model", "result", "note"))
  # Role-specific geometry: input parallelograms, process rectangles,
  # model hexagons, evaluation ellipses; dashed notes are not processing steps.
  if (shape == "data") {
    polygon(c(x1, x2 - .025, x2, x1 + .025), c(y1, y1, y2, y2),
            col = "#EEF4F8", border = "grey40", lwd = .8)
  } else if (shape == "model") {
    polygon(c(x1, x1 + .03, x2 - .03, x2, x2 - .03, x1 + .03),
            c((y1 + y2) / 2, y1, y1, (y1 + y2) / 2, y2, y2),
            col = "#F0EDF5", border = "grey40", lwd = .8)
  } else if (shape == "result") {
    angle <- seq(0, 2 * pi, length.out = 101)
    polygon((x1 + x2) / 2 + (x2 - x1) / 2 * cos(angle),
            (y1 + y2) / 2 + (y2 - y1) / 2 * sin(angle),
            col = "#EDF5EF", border = "grey40", lwd = .8)
  } else {
    rect(x1, y1, x2, y2, col = if (shape == "note") "white" else "grey97",
         border = "grey40", lwd = .8, lty = if (shape == "note") 2 else 1)
  }
  text((x1 + x2) / 2, (y1 + y2) / 2, label, cex = cex)
  invisible(NULL)
}

#' Connect two workflow stages with a directional arrow
#' @param x1,y1,x2,y2 Numeric endpoints in normalized figure coordinates.
#' @return Invisibly NULL; draws an arrow without changing the plot limits.
#' @examples
#' # After plot.new(): stage_arrow(.25, .8, .25, .75)
stage_arrow <- function(x1, y1, x2, y2) {
  arrows(x1, y1, x2, y2, length = 0.065, lwd = 0.9, col = "grey25")
  invisible(NULL)
}

#' Draw one of the two fixed historical real-data workflows
#' @param model Character, either "single" or "multimodal".
#' @return Invisibly NULL; draws the selected workflow on the active device.
#' @examples
#' # pdf(tempfile(fileext = ".pdf")); draw_workflow("single"); dev.off()
draw_workflow <- function(model) {
  stopifnot(model %in% c("single", "multimodal"))
  par(mar = c(0, 0, 0, 0), family = "sans", xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  text(.25, .975, "TRAINING", font = 2, cex = .9)
  text(.75, .975, "EXTERNAL VALIDATION", font = 2, cex = .9)

  if (model == "single") {
    stage_box(.02, .755, .47, .94,
              "TCGA: 144 patients, 75 deaths\nCPTAC: 129 patients, 64 deaths\nBoth RNA-seq\nExpression + survival", shape = "data")
    stage_box(.53, .755, .98, .94,
              "Dijk: 90 records (RNA-seq)\nMoffitt: 123 (microarray)\nPACA-AU: 63 (microarray)\nPACA-AU: 52 (RNA-seq)\nPuleo: 288 (microarray)", cex = .76, shape = "data")
    stage_box(.02, .50, .47, .70,
              "Historical log preprocessing\nTop 3,000 genes per cohort\n2,064 genes in the intersection\nStandardize within each cohort\nPool: 273 patients x 2,064 genes")
    stage_box(.53, .50, .98, .70,
              "Historical log preprocessing\nUse measured training genes\nStandardize within each cohort\n1,989-2,050 genes per series\nNo external outcome screening")
    stage_box(.02, .255, .47, .44,
              "Fit single-block YFB, starting K = 7\nLoading survival weight = 0\nFreeze loadings, stored norms,\ncoefficients and risk direction", cex = .76, shape = "model")
    stage_box(.53, .255, .98, .44,
              "Project each external expression\nmatrix onto frozen loadings\nApply stored norms and coefficients\nObtain risk scores without refitting", cex = .76)
    stage_box(.02, .085, .47, .24,
              "RANK SELECTION\n5-fold training CV selected K = 7.\nExternal concordance later informed\nits assessment; this was not a\nfully untouched external test.", cex = .70, shape = "note")
    stage_box(.53, .10, .98, .205,
              "Compare risk with observed survival\nConcordance in each series\nand descriptive survival curves", cex = .73, shape = "result")
    stage_box(.02, .003, .98, .073,
              "616 records; 568 unique people if other cohorts are distinct\nPACA-AU series share 47 donors; RNA-seq has 52 records from 51 donors",
              cex = .70, shape = "note")
  } else {
    stage_box(.02, .755, .47, .94,
              "TCGA: 150 matched patients\nExclude 6 without usable survival\n144 patients, 75 deaths\nExpression + methylation + survival", shape = "data")
    stage_box(.53, .755, .98, .94,
              "ICGC = PACA-AU in this study\n79 matched donors\nExclude 12 without usable survival\n67 donors, 40 deaths", shape = "data")
    stage_box(.02, .50, .47, .70,
              "Shared assay feature universe\nImpute methylation within TCGA\nScreen on TCGA only:\n3,000 genes + 10,000 CpGs\nLog expression; no z-standardization")
    stage_box(.53, .50, .98, .70,
              "Same training-selected features\nImpute methylation within ICGC\nLog expression; methylation values\nNo z-standardization\nNo outcome-based feature selection")
    stage_box(.02, .255, .47, .44,
              "Fit joint expression-methylation YFB\nShared activities; separate loadings\nFreeze modality loadings,\ncoefficients and risk direction", cex = .74, shape = "model")
    stage_box(.53, .255, .98, .44,
              "Project both assays for each donor\nSum their program projections\nApply frozen survival coefficients\nObtain risk without refitting", cex = .76)
    stage_box(.02, .085, .47, .24,
              "PATIENT MATCHING\nEach row: expression, methylation\nand survival for one person.\nTCGA trains; ICGC validates.\nPatients are not pooled across cohorts.", cex = .70, shape = "note")
    stage_box(.515, .095, .74, .225,
              "PRIMARY PDAC\n50 donors; 30 deaths\nConcordance", cex = .70, shape = "result")
    stage_box(.755, .095, .995, .225,
              "BROADER SET\n67 donors; 40 deaths\nIncludes primary 50", cex = .70, shape = "result")
    stage_box(.02, .003, .98, .073,
              "Broader set = primary 50 + 17 donors with other tumor histologies\nThe same fitted model is evaluated in both sets; these are not independent cohorts",
              cex = .69, shape = "note")
  }
  # Separate lanes prevent the evaluation outcomes from appearing as fit inputs.
  for (x in c(.245, .755)) {
    stage_arrow(x, .755, x, .70)
    stage_arrow(x, .50, x, .44)
  }
  stage_arrow(.47, .35, .53, .35)
  if (model == "single") {
    stage_arrow(.755, .255, .755, .205)
  } else {
    stage_arrow(.70, .255, .6275, .225)
    stage_arrow(.81, .255, .875, .225)
  }
  invisible(NULL)
}

# Export identical content for print and HTML ----
# Vector PDF keeps text sharp in the manuscript; PNG avoids browser-specific
# SVG/font behavior. Explicit devices have identical dimensions and font size.
out_dir <- "paper/writeup/figures"
for (model in c("single", "multimodal")) {
  stem <- file.path(out_dir, paste0("real_data_workflow_", model))
  pdf(paste0(stem, ".pdf"), width = 6.5, height = 5.5, pointsize = 12,
      family = "Helvetica", useDingbats = FALSE)
  draw_workflow(model)
  dev.off()
  png(paste0(stem, ".png"), width = 6.5, height = 5.5, units = "in",
      res = 300, pointsize = 12)
  draw_workflow(model)
  dev.off()
  cat("Wrote", paste0(stem, ".{pdf,png}"), "\n")
}
