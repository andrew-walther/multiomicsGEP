# ============================================================
# Script: test_multimodal_yfb_helpers.R
# Purpose: Test mathematical helpers for the isolated multimodal YFB model.
# Author: Andrew Walther
# Created: 2026-09-17
# Dependencies: base R
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
if (!exists("multimodal_yfb_cox_working")) {
  source("code/multimodal_yfb_helpers.R")
}

cat("\n========================================\n")
cat("  Tests: multimodal_yfb_helpers.R\n")
cat("========================================\n\n")

# T1: Cox working quantities ----

run_test("MMYFB-Helpers-T1: Breslow working quantities handle tied events", {
  working <- multimodal_yfb_cox_working(
    time = c(1, 2, 2, 4), event = c(1, 1, 1, 0),
    eta = c(0.2, -0.3, 0.4, -0.1)
  )
  assert_near(working$u,
              c(0.719790, 0.357806, -0.293219, -0.784377), tol = 1e-5)
  assert_near(working$w,
              c(0.280210, 0.642194, 1.293219, 0.784377), tol = 1e-5)
  assert_near(working$h, working$u + working$w * c(0.2, -0.3, 0.4, -0.1),
              tol = 1e-12)
})

# T2: Projection and predictor moments ----

run_test("MMYFB-Helpers-T2: projection moments pool modalities and loading uncertainty", {
  Y <- list(
    expression = matrix(c(1, 2, 3, 4), nrow = 2,
                        dimnames = list(c("s1", "s2"), c("g1", "g2"))),
    methylation = matrix(c(2, 1), nrow = 2,
                          dimnames = list(c("s1", "s2"), "cg1"))
  )
  EF <- list(expression = matrix(c(0.5, 1.0), ncol = 1),
             methylation = matrix(2, ncol = 1))
  EF2 <- list(expression = matrix(c(0.29, 1.04), ncol = 1),
              methylation = matrix(4.25, ncol = 1))
  moments <- multimodal_yfb_projection_moments(Y, EF, EF2)
  assert_near(moments$EZ[, 1], c(7.5, 7), tol = 1e-12)
  assert_near(moments$VZ[, 1], c(1.4, 1.05), tol = 1e-12)
  assert_near(moments$EZ2[, 1], c(57.65, 50.05), tol = 1e-12)
})

run_test("MMYFB-Helpers-T3: predictor variance uses beta second moments", {
  predictor <- multimodal_yfb_predictor_moments(
    EZ = matrix(c(2, 3), ncol = 1), VZ = matrix(c(0.5, 1), ncol = 1),
    EBeta = 2, EBeta2 = 4.25
  )
  assert_near(predictor$mean, c(4, 6), tol = 1e-12)
  assert_near(predictor$variance, c(3.125, 6.5), tol = 1e-12)
})

# T3: Residuals and moment contracts ----

run_test("MMYFB-Helpers-T4: factor-excluded residual removes only other factors", {
  Y <- list(expression = matrix(c(5, 7), ncol = 1),
            methylation = matrix(c(6, 8), ncol = 1))
  EL <- cbind(c(1, 2), c(3, 4))
  EF <- list(expression = rbind(c(2, 1)), methylation = rbind(c(1, 2)))
  residual <- multimodal_yfb_residual_minus_k(Y, EL, EF, k = 1)
  assert_near(residual$expression[, 1], c(2, 3), tol = 1e-12)
  assert_near(residual$methylation[, 1], c(0, 0), tol = 1e-12)
})

run_test("MMYFB-Helpers-T5: invalid posterior second moments fail loudly", {
  error <- tryCatch({
    multimodal_yfb_validate_moments(mean = c(1, 2), second = c(0.9, 4),
                                    name = "EL")
    NULL
  }, error = conditionMessage)
  assert_true(grepl("second", error, ignore.case = TRUE))
})

if (sys.nframe() == 0L) report_results("test_multimodal_yfb_helpers.R")
