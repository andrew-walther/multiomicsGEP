# ============================================================
# Script: test_methylation_single_modality.R
# Purpose: Check methylation study preprocessing and summary helpers.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: base R
# ============================================================

local({
source("results/multimodal_real/run_methylation_single_modality.R", local = TRUE)

# Training-derived scaling must also be used for held-out observations. ----
train <- cbind(a = c(1, 2, 3), b = c(2, 4, 6))
held <- matrix(c(10, 20), nrow = 1, dimnames = list(NULL, colnames(train)))
pl <- c_scale(train, "PL")
pe <- c_scale(train, "PE")
stopifnot(isTRUE(all.equal(unname(pl$center), c(2, 4))),
          isTRUE(all.equal(unname(pl$sd), c(1, 2))),
          isTRUE(all.equal(unname(pl$Y), unname(scale(train)), check.attributes = FALSE)),
          isTRUE(all.equal(unname(c_transform(held, pl)), matrix(c(8, 8), 1))),
          isTRUE(all.equal(unname(pe$center), c(0, 0))),
          isTRUE(all.equal(unname(pe$Y), unname(sweep(train, 2, c(1, 2), "/")))),
          isTRUE(all.equal(unname(c_transform(held, pe)), matrix(c(10, 10), 1))),
          all(pe$Y >= 0))

train2 <- cbind(a = c(4, 8, 12), b = c(5, 8, 11))
pl2 <- c_scale(train2, "PL")
pe2 <- c_scale(train2, "PE")
stopifnot(isTRUE(all.equal(unname(pl2$center), c(8, 8))),
          isTRUE(all.equal(unname(pl2$sd), c(4, 3))),
          isTRUE(all.equal(unname(pl2$Y), matrix(rep(c(-1, 0, 1), 2), 3))),
          isTRUE(all.equal(unname(pe2$Y), unname(sweep(train2, 2, c(4, 3), "/")))))

# A fixed seed provides the same event-balanced folds for every study arm. ----
event <- rep(c(0L, 1L), each = 10)
f <- c_folds(event, seed = 2026L)
stopifnot(identical(f, c_folds(event, seed = 2026L)), length(f) == length(event),
          setequal(f, 1:5), all(table(event, f) == 2))

event2 <- c(rep(0L, 13), rep(1L, 17))
f2 <- c_folds(event2, seed = 37L)
counts2 <- table(event2, factor(f2, levels = 1:5))
stopifnot(length(f2) == length(event2), all(f2 %in% 1:5),
          identical(f2, c_folds(event2, seed = 37L)),
          all(apply(counts2, 1, function(x) max(x) - min(x) <= 1)))

# Match permuted and sign-reversed programs without reusing a fold program. ----
a <- c(-2, -1, 0, 1, 2, 3)
b <- c(1, -1, 0, 0, -1, 1)
full <- cbind(a, b, rep(1, length(a)))
fold <- cbind(-b, rep(4, length(a)), -a)
matched <- c_match(full, fold)
stopifnot(identical(as.integer(matched), c(3L, 1L, 2L)),
          length(unique(matched)) == 3L)

stopifnot(identical(c_match(cbind(a, b), cbind(a, b)), c(1L, 2L)),
          identical(c_match(cbind(a, b), cbind(-b, a)), c(2L, 1L)))

# Pair concordances by fold and propagate undefined folds. ----
meth <- c(0.6, 0.7, 0.5, 0.8, 0.6)
expression <- c(0.5, 0.65, 0.6, 0.7, 0.55)
d <- meth - expression
paired <- c_paired(meth, expression)
stopifnot(isTRUE(all.equal(unname(paired), c(mean(d), sd(d) / sqrt(5)))),
          all(is.na(c_paired(replace(meth, 3, NA_real_), expression))))

second <- c_paired(c(.2, .4, .6, .8, 1), rep(.2, 5))
stopifnot(isTRUE(all.equal(unname(second), c(.4, sqrt(.1 / 5)))))

# Alive means exactly nonzero EF, and activity requires positive posterior SD. ----
fit <- list(EF = matrix(c(1e-12, 0, 0, 0, 1, 0, 1, 0), nrow = 2),
            EBeta = c(2, 3, 2, 4), EBeta2 = c(5, 10, 4, 20))
activity <- c_activity(fit)
stopifnot(identical(as.logical(activity$alive), c(TRUE, FALSE, TRUE, TRUE)),
          identical(as.logical(activity$active), c(TRUE, FALSE, FALSE, TRUE)))

fit2 <- fit
fit2$EBeta <- c(0, 3, 3, 0)
fit2$EBeta2 <- c(1, 10, 10, 1)
activity2 <- c_activity(fit2)
stopifnot(identical(as.logical(activity2$active), c(FALSE, FALSE, TRUE, FALSE)))

cat("Methylation single-modality helper tests passed.\n")

# A child startup failure must yield a failed row, and resume can lack RSS. ----
saved_local <- C_LOCAL
C_LOCAL <- tempfile("methylation-test-")
dir.create(C_LOCAL)
# Mock just the subprocess boundary, leaving the failure-recording code real.
#' Simulate an unsuccessful subprocess for failure-recording checks
#' @param ... Ignored command arguments.
#' @return Integer exit code 9.
#' @examples
#' # system2('unused')
system2 <- function(...) 9L
failed <- c_launch("methylation", 10L, "PL", 0L)
stopifnot(identical(failed$status, "failed"), !failed$converged,
          grepl("Child exit 9", failed$error), !failed$memory_measured,
          is.na(failed$peak_rss_MiB),
          file.exists(file.path(C_LOCAL, "methylation_10_PL_f0.rds")))
rm(system2)
resumed <- c_launch("methylation", 10L, "PL", 0L)
stopifnot(identical(resumed$status, "failed"), !resumed$memory_measured)
unlink(C_LOCAL, recursive = TRUE)
C_LOCAL <- saved_local
cat("Methylation subprocess failure/resume checks passed.\n")


# CV coverage detects missing and repeated worker outputs without starting fits. ----
C_LOCAL <- tempfile("methylation-coverage-local-")
C_OUT <- tempfile("methylation-coverage-out-")
dir.create(C_LOCAL)
dir.create(C_OUT)
saveRDS(list(counts = 10L), file.path(C_LOCAL, "design.rds"))
coverage_case <- "success"

#' Return synthetic fit metrics for a requested study job
#' @param modality,count,arm,fold Study job identifiers.
#' @return One metric row, or NULL for the designated missing worker.
#' @examples
#' # c_launch("methylation", 10, "PL", 1)
c_launch <- function(modality, count, arm, fold) {
  if (fold > 0L && modality == "expression" && arm == "PL" && fold == 1L &&
      coverage_case == "missing") return(NULL)
  if (fold > 0L && modality == "expression" && arm == "PL" && fold == 1L &&
      coverage_case == "duplicate") fold <- 2L
  data.frame(modality = modality, count = count, arm = arm, fold = fold,
             converged = TRUE, seconds = 1, process_seconds = 1)
}

#' Execute a mock parallel map sequentially in the current process
#' @param X,FUN Inputs and worker function.
#' @param ... Ignored parallel options.
#' @return List of worker results, retaining NULL entries.
#' @examples
#' # sequential_map(1:2, identity)
sequential_map <- function(X, FUN, ...) lapply(X, FUN)
parallel_ns <- asNamespace("parallel")
original_map <- get("mclapply", envir = parallel_ns)
unlockBinding("mclapply", parallel_ns)
assign("mclapply", sequential_map, envir = parallel_ns)
lockBinding("mclapply", parallel_ns)
tryCatch({
  c_run()
  cv_rows <- read.csv(file.path(C_OUT, "cv_fits.csv"))
  keys <- with(cv_rows, paste(modality, count, arm, fold))
  expected <- expand.grid(modality = c("expression", "methylation"),
                          arm = c("PL", "PE"), fold = 1:5,
                          stringsAsFactors = FALSE)
  expected$count <- ifelse(expected$modality == "expression", 2218L, 10L)
  expected_keys <- with(expected, paste(modality, count, arm, fold))
  stopifnot(nrow(cv_rows) == 20L, !anyDuplicated(keys), setequal(keys, expected_keys))
  coverage_case <- "missing"
  missing_error <- tryCatch({ c_run(); "" }, error = conditionMessage)
  stopifnot(grepl("Missing/failed CV worker result", missing_error, fixed = TRUE))
  coverage_case <- "duplicate"
  duplicate_error <- tryCatch({ c_run(); "" }, error = conditionMessage)
  stopifnot(grepl("CV result coverage differs from the registered task list",
                  duplicate_error, fixed = TRUE))
}, finally = {
  unlockBinding("mclapply", parallel_ns)
  assign("mclapply", original_map, envir = parallel_ns)
  lockBinding("mclapply", parallel_ns)
  unlink(c(C_LOCAL, C_OUT), recursive = TRUE)
})
cat("Methylation CV coverage checks passed.\n")

})
