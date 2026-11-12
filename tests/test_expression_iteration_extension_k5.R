# ============================================================
# Script: test_expression_iteration_extension_k5.R
# Purpose: Verify K = 5 trajectory reproduction, cap isolation and paired CV arithmetic.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: yaml; expression extension helpers
# ============================================================
local({
source("results/multimodal_real/run_expression_iteration_extension_k5.R", local = TRUE)
# A changed prefix, short history or nonfinite history must fail.
x <- -c(100, 90, 82)
stopifnot(e5_prefix(x, c(x, -80)) == 0,
          e5_prefix(x, c(x + 1e-12, -80)) < 1e-10)
for (bad in list(c(x + .01, -80), x[1:2], c(x[1:2], NA)))
  stopifnot(inherits(try(e5_prefix(x, bad), silent = TRUE), "try-error"))
cfg <- list(cavi = list(max_iter = 300L, tol = 1e-5), benchmark = list(alpha = .7))
new <- e5_config(cfg)
stopifnot(new$cavi$max_iter == 1000L, cfg$cavi$max_iter == 300L)
new$cavi$max_iter <- cfg$cavi$max_iter
stopifnot(identical(cfg, new))
m <- c(.5, .7, .4, .8, .65); e <- c(.6, .5, .6, .7, .55)
a <- e5_paired(m, e); b <- e5_paired(m, e + .1)
stopifnot(abs(a[1] - mean(m-e)) < 1e-14,
          abs(a[2] - sd(m-e)/sqrt(5)) < 1e-14,
          abs(b[1] - (a[1]-.1)) < 1e-14,
          is.na(e5_paired(replace(m, 2, NA), e)[1]))
root <- tempfile("e5-"); dir.create(root)
on.exit(unlink(root, recursive = TRUE))
E5_BASE <- file.path(root, "private"); E5_OUT <- file.path(root, "aggregate")
dir.create(file.path(E5_BASE, "source/config"), recursive = TRUE); dir.create(E5_OUT)
saveRDS(list(fold = 1:5), file.path(E5_BASE, "design.rds"))
yaml::write_yaml(cfg, file.path(E5_BASE, "source/config/globals.yml"))
writeLines("# frozen", file.path(E5_BASE, "source/runner.R"))
paths <- c("config/globals.yml", "runner.R")
write.csv(data.frame(path = paths, md5 = unname(tools::md5sum(file.path(E5_BASE, "source", paths)))),
          file.path(E5_OUT, "source_hashes.csv"), row.names = FALSE)
e5_prepare()
extended_cfg <- file.path(E5_BASE, "iteration1000/source/config/globals.yml")
stopifnot(yaml::read_yaml(extended_cfg)$cavi$max_iter == 1000L,
          identical(readRDS(file.path(E5_BASE, "iteration1000/design.rds")), list(fold = 1:5)))
changed <- yaml::read_yaml(extended_cfg); changed$cavi$tol <- 1e-4
yaml::write_yaml(changed, extended_cfg)
stopifnot(inherits(try(e5_prepare(), silent = TRUE), "try-error"))
cat("K = 5 expression extension: 15 assertions passed.\n")
})
