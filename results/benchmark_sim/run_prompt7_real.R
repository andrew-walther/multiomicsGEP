# ============================================================
# Script: run_prompt7_real.R
# Purpose: Evaluate frozen risk groups, program HRs and ranked EBMF comparators.
# Created: 2026-11-06
# Dependencies: survival, ggplot2, flashier and Prompt 6 helpers
# ============================================================

source("results/benchmark_sim/run_prompt6_comparators.R")
source("code/fit_ssbmf_yfb.R")
source("code/simulate_multimodal_yfb.R")
OUT7 <- exhibit_out("results/benchmark_sim/outputs/prompt7_desurv_style", "prompt7_desurv_style")  # EXHIBIT_K

#' Derive one maximally separated cutoff from standardized OOF training scores.
#' @param risk,time,event OOF score and aligned training survival vectors.
#' @return One row with cutoff, training denominator, and the log-rank statistic.
#' @examples
#' # prompt7_cutpoint(oof, time, event)
prompt7_cutpoint <- function(risk, time, event) {
  stopifnot(all(is.finite(risk)))
  candidates <- sort(unique(risk))
  candidates <- (head(candidates, -1L) + tail(candidates, -1L)) / 2
  low <- vapply(candidates, function(x) sum(risk <= x), integer(1))
  keep <- pmin(low, length(risk) - low) >= .2 * length(risk)
  x <- candidates[keep]
  stat <- vapply(x, function(cut) survdiff(Surv(time, event) ~ I(risk > cut))$chisq,
                 numeric(1))
  best <- order(-stat, abs(x), x)[1]
  stopifnot(length(best) == 1L, !is.na(best))
  data.frame(cutpoint = x[best], n = length(risk), low = sum(risk <= x[best]),
             high = sum(risk > x[best]), logrank_chisq = stat[best])
}

#' Fit a descriptive external Cox HR and preserve estimability and sample counts.
#' @param dat Data frame with x, time, event, cohort and optionally subtype.
#' @param pooled Logical. Use cohort-specific baseline hazards when pooled.
#' @param adjusted Logical. Include subtype if requested.
#' @return One row with HR, limits, p-value, and fit status.
#' @examples
#' # prompt7_hr(dat, FALSE)
prompt7_hr <- function(dat, pooled, adjusted = FALSE) {
  # Review addition: pooled fits add cluster(donor) for a robust interval,
  # because 47 PACA-AU donors contribute records from both assays. The
  # registered model-based interval (lo, hi) is unchanged.
  formula <- reformulate(c("x", if (adjusted) "subtype",
                           if (pooled) c("strata(cohort)", "cluster(donor)")), "Surv(time, event)")
  output <- data.frame(n = nrow(dat), events = sum(dat$event), hr = NA_real_,
                       lo = NA_real_, hi = NA_real_, p = NA_real_,
                       donor_lo = NA_real_, donor_hi = NA_real_, status = "unestimable")
  if (length(unique(dat$x)) < 2L || sum(dat$event) < 2L ||
      (adjusted && length(unique(dat$subtype)) < 2L)) return(output)
  warnings <- character()
  fit <- withCallingHandlers(coxph(formula, data = dat, na.action = na.fail),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w)))
  if (is.finite(coef(fit)["x"]) && !length(warnings)) {
    b <- coef(fit)[["x"]]
    se <- sqrt((if (pooled) fit$naive.var else vcov(fit))[1, 1])
    output[1, c("hr", "lo", "hi", "p")] <- c(exp(b), exp(b + c(-1, 1) * qnorm(.975) * se),
                                            2 * pnorm(-abs(b / se)))
    if (pooled) {
      rse <- sqrt(vcov(fit)[1, 1])
      output[1, c("donor_lo", "donor_hi")] <- exp(b + c(-1, 1) * qnorm(.975) * rse)
    }
    output$status <- "ok"
  }
  if (length(warnings)) output$status <- paste(warnings, collapse = " | ")
  output
}

#' Compute coefficient ratios and signed risk variance allocations.
#' @param Y Training matrix. @param F Loading matrix. @param beta Coefficient vector.
#' @return One row per loading with coefficient ratio and signed variance share.
#' @examples
#' # prompt7_contribution(Y, F, beta)
prompt7_contribution <- function(Y, F, beta) {
  component <- sweep(Y %*% F, 2, beta, "*")
  eta <- rowSums(component)
  sb <- apply(component, 2, sd)
  data.frame(program = seq_along(beta), ratio = if (max(sb) > 0) sb / max(sb) else NA_real_,
    variance_share = if (var(eta) > 0) drop(cov(component, eta)) / var(eta) else NA_real_)
}

#' Identify live posterior-active coefficients.
#' @param fit Fixed YFB checkpoint state with EF and coefficient moments.
#' @return Logical vector identifying active programs.
#' @examples
#' # prompt7_active(d$w1)
prompt7_active <- function(fit) {
  se <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
  colSums(fit$EF^2) > 1e-8 & se > 0 & abs(fit$EBeta) >= 1.96 * se
}

#' Match live loading directions one to one by absolute correlation.
#' @param ref,estimate Gene by program loading matrices.
#' @return Matches in original column indices, or NULL if no live columns.
#' @examples
#' # prompt7_matches(d$w1$EF, fold$w1$EF)
prompt7_matches <- function(ref, estimate) {
  i <- which(apply(ref, 2, sd) > 0)
  j <- which(apply(estimate, 2, sd) > 0)
  if (!length(i) || !length(j)) return(NULL)
  m <- multimodal_yfb_match_factors(ref[, i, drop = FALSE], estimate[, j, drop = FALSE])
  m$true_factor <- i[m$true_factor]
  m$estimated_factor <- j[m$estimated_factor]
  m
}

if (sys.nframe() == 0L) {
  d <- prompt6_reference()
  cohorts <- cfg$pdac$external_cohorts
  folds <- sort(unique(d$fold))
  stopifnot(identical(as.integer(folds), 1:5))
  arms <- c("w1", "twostep")
  oof <- matrix(NA_real_, nrow(d$Y), length(arms), dimnames = list(NULL, arms))
  matches <- list()
  for (f in folds) {
    df <- prompt6_reference(f)
    stopifnot(identical(df$genes, d$genes), identical(df$fold, d$fold))
    r <- prompt6_risk(d$Y, df)
    hold <- d$fold == f
    for (arm in arms) {
      stopifnot(all(is.na(oof[hold, arm])), sd(r[[arm]][!hold]) > 0)
      oof[hold, arm] <- (r[[arm]][hold] - mean(r[[arm]][!hold])) / sd(r[[arm]][!hold])
    }
    m <- prompt7_matches(d$w1$EF, df$w1$EF)
    matches[[f]] <- transform(m, fold = f, active_match = abs_correlation >= .7 &
                                prompt7_active(df$w1)[estimated_factor])
  }
  cutoff <- do.call(rbind, lapply(arms, function(arm)
    cbind(model = arm, prompt7_cutpoint(oof[, arm], d$time, d$status))))
  matches <- do.call(rbind, matches)
  active <- prompt7_active(d$w1)
  program <- prompt7_contribution(d$Y, sweep(d$w1$EF, 2, d$w1$EF_norms, "/"), d$w1$EBeta)
  stopifnot(abs(sum(program$variance_share) - 1) < 1e-10)
  program$active <- active
  program$active_folds <- vapply(program$program, function(k)
    sum(matches$true_factor == k & matches$active_match), integer(1))
  program$main <- with(program, active & ratio >= .25 & variance_share >= .05 & active_folds >= 4)
  program$pve <- colSums(d$w1$EL^2) * colSums(d$w1$EF^2) / sum(d$Y^2)
  # E[L²] is unavailable in the checkpoint. E[L]² gives a lower bound on
  # genomics precision and hence an upper bound on the median survival share.
  risk <- prompt6_risk(d$Y, d)$w1
  w <- calc_cox_taylor_yf(risk, d$time, d$status)$w
  survival_precision <- outer(colSums(d$Y^2 * w), d$w1$EBeta2 / d$w1$EF_norms^2)
  genomics_lower <- outer(d$w1$Tau, colSums(d$w1$EL^2))
  program$precision_share_upper <- apply(survival_precision /
    (survival_precision + genomics_lower), 2, median)
  program$old_survival_only <- active & program$pve < .01
  program$precision_survival_only <- ifelse(!active | program$precision_share_upper <= .5,
                                            FALSE, NA)
  result <- list(cutpoints = cutoff, program_rule = program, fold_matches = matches)

  training_risk <- prompt6_risk(d$Y, d)
  external <- setNames(lapply(cohorts, prompt6_external, genes = d$genes), cohorts)
  group_hr <- program_hr <- adjusted_hr <- enrichment <- km <- at_risk <- list()
  for (arm in arms) {
    F <- if (arm == "w1") d$w1$EF else d$w0$EF
    use <- if (arm == "w1") active else {
      live <- which(colSums(F^2) > 1e-8)
      seq_len(ncol(F)) %in% live[abs(coef(d$cox0) / sqrt(diag(vcov(d$cox0)))) >= 1.96]
    }
    program_sd <- apply(d$Y %*% F, 2, sd)
    train_mu <- mean(training_risk[[arm]])
    train_sd <- sd(training_risk[[arm]])
    data <- do.call(rbind, lapply(cohorts, function(co) {
      e <- external[[co]]
      risk <- (prompt6_risk(e$Y, d)[[arm]] - train_mu) / train_sd
      data.frame(time = e$time, event = e$status, cohort = co, risk = risk,
        high = as.integer(risk > cutoff$cutpoint[cutoff$model == arm]),
        e$meta[, c("PurIST", "DeCAF", "donor")])
    }))
    P <- do.call(rbind, lapply(external, function(e) e$Y %*% F))
    for (co in c(cohorts, "pooled")) {
      pooled <- co == "pooled"
      i <- if (pooled) seq_len(nrow(data)) else which(data$cohort == co)
      a <- data[i, ]
      group_hr[[length(group_hr) + 1L]] <- cbind(model = arm, cohort = co,
        low = sum(a$high == 0), high = sum(a$high == 1),
        prompt7_hr(transform(a, x = high), pooled))
      for (k in which(use)) {
        program_hr[[length(program_hr) + 1L]] <- cbind(model = arm, cohort = co, program = k,
          prompt7_hr(transform(a, x = P[i, k] / program_sd[k]), pooled))
      }
      for (nm in c("PurIST", "DeCAF")) {
        target <- if (nm == "PurIST") "Basal-like" else "permCAF"
        labels <- if (nm == "PurIST") c("Classical", "Basal-like") else c("restCAF", "permCAF")
        keep <- a[[nm]] %in% labels
        b <- transform(a[keep, ], subtype = factor(a[[nm]][keep], levels = labels), x = risk)
        for (adjusted in c(FALSE, TRUE)) {
          adjusted_hr[[length(adjusted_hr) + 1L]] <- cbind(model = arm, cohort = co,
            subtype = nm, adjusted = adjusted, dropped = sum(!keep),
            prompt7_hr(b, pooled, adjusted))
        }
        if (!pooled) {
          ct <- table(factor(b$high, levels = 0:1), factor(b$subtype == target, c(FALSE, TRUE)))
          n <- rowSums(ct)
          enrichment[[length(enrichment) + 1L]] <- data.frame(model = arm, cohort = co,
            subtype = nm, target = target, dropped = sum(!keep), low_n = n[1], high_n = n[2],
            low_positive = ct[1, 2], high_positive = ct[2, 2], low_prop = ct[1, 2] / n[1],
            high_prop = ct[2, 2] / n[2], p = if (all(n > 0)) fisher.test(ct)$p.value else NA_real_)
        }
      }
      if (!pooled) for (g in sort(unique(a$high))) {
        b <- a[a$high == g, ]
        sf <- survfit(Surv(time, event) ~ 1, data = b)
        # Evaluate the KM estimate on a monthly grid before publishing.
        # This preserves the survival summary without exposing exact event times.
        months <- seq.int(0, min(60, floor(max(b$time))))
        km[[length(km) + 1L]] <- data.frame(model = arm, cohort = co, high = g,
          month = months, survival = c(1, sf$surv)[findInterval(months, sf$time) + 1L])
        ticks <- seq(0, 60, 12)
        at_risk[[length(at_risk) + 1L]] <- data.frame(model = arm, cohort = co, high = g,
          month = ticks, n = vapply(ticks, function(t) sum(b$time >= t), integer(1)))
      }
    }
  }
  result <- c(result, lapply(list(group_hr = group_hr, program_hr = program_hr,
    adjusted_hr = adjusted_hr, enrichment = enrichment, at_risk = at_risk, km_monthly = km), function(x) do.call(rbind, x)))
  dir.create(OUT7, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(result)) write.csv(result[[nm]], file.path(OUT7, paste0(nm, ".csv")), row.names = FALSE)
  # Publish monthly KM estimates and prespecified at-risk counts.
  library(ggplot2)
  for (arm in arms) {
    curves <- result$km_monthly[result$km_monthly$model == arm, ]
    counts <- result$at_risk[result$at_risk$model == arm, ]
    p <- ggplot(curves, aes(month, survival, colour = factor(high))) + geom_step() +
      geom_text(data = counts, aes(month, -.08 - .08 * high, label = n),
                size = 3, show.legend = FALSE) +
      facet_wrap(~ cohort, ncol = 2) + scale_color_discrete("Risk group", labels = c("Low", "High")) +
      scale_x_continuous("Months (counts at risk below each panel)", breaks = seq(0, 60, 12)) +
      scale_y_continuous("Survival probability", breaks = seq(0, 1, .25)) + theme_bw()
    ggsave(file.path(OUT7, paste0("km_", arm, ".pdf")), p, width = 9, height = 8)
  }
  forest <- result$program_hr[result$program_hr$status == "ok", ]
  forest$label <- paste(forest$model, "program", forest$program)
  p <- ggplot(forest, aes(hr, label)) + geom_vline(xintercept = 1, linetype = "dashed") +
    geom_pointrange(aes(xmin = lo, xmax = hi), orientation = "y") + scale_x_log10() +
    facet_wrap(~ cohort, ncol = 2) + labs(x = "Hazard ratio per training SD", y = NULL) + theme_bw()
  ggsave(file.path(OUT7, "program_forest.pdf"), p, width = 10, height = 8)
}
