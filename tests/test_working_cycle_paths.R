# ============================================================
# Script: test_working_cycle_paths.R
# Purpose: Keep the invalid and corrected registered diagnostic outputs separate.
# Created: 2026-10-16
# Dependencies: existing test helpers and audit driver
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
source("results/benchmark_sim/run_working_cycle_audit.R")
run_test("CyclePaths-T1: corrected attempt uses fresh fixed output and telemetry paths", {
  original<-cycle_audit_paths(character())
  corrected<-cycle_audit_paths("corrected")
  assert_equal(original$out,"results/benchmark_sim/outputs/round2/working_cycle_seed26")
  assert_equal(corrected$out,"results/benchmark_sim/outputs/round2/working_cycle_seed26_corrected")
  assert_equal(corrected$detailed,"data/round2_working_cycle_seed26_corrected.rds")
  assert_false(identical(original,corrected),msg="invalid attempt cannot be overwritten")
})
run_test("CyclePaths-T2: arbitrary output modes are rejected", {
  for(mode in list("retry2",c("corrected","extra"))) {
    error<-tryCatch(cycle_audit_paths(mode),error=function(e)e)
    assert_true(inherits(error,"error") && grepl("Only the original or registered corrected attempt",
      conditionMessage(error),fixed=TRUE),msg="reject the mode itself, not an unrelated error")
  }
})
