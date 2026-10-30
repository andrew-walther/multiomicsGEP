# The generated Cox direction and reconstruction loading must share units. ----
source("results/benchmark_sim/run_yfb_supervision_decision.R")
run_test("YFBS simulation: prepared reconstruction and Cox use one loading", local({
  d<-yfb_decision_data(1,50,.01,truth_mode="shared")
  expected<-d$raw_F3/d$training_sd
  expected<-expected/sqrt(sum(expected^2))
  assert_near(d$truth,expected,1e-12)
  eta<-drop(scale(d$Y,scale=FALSE)%*%expected)*3.3
  assert_near(d$training_eta,eta,1e-12)
  assert_true(sd(eta)>.1)
  component<-outer(d$activity3,d$raw_F3/d$training_sd)
  assert_near(d$training_component3,component,1e-12)
  assert_near(d$actual_share,.01,1e-12)
}))
