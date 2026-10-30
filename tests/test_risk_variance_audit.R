# Regression checks for the documented Cox hazard convention and β moments. ----
run_test("C6: positive Cox coefficient denotes adverse risk", {
  x <- seq(-2, 2, length.out = 40)
  time <- exp(-x)
  model <- survival::coxph(survival::Surv(time, rep(1, 40)) ~ x,
                          control = survival::coxph.control(iter.max = 5))
  assert_true(unname(stats::coef(model)) > 0)
  assert_near(frozen_reverse_cindex(x, time, rep(1, 40)), 1)
})
run_test("C6: fixed nondegenerate β prior retains uncertainty at zero signal", {
  ans <- update_beta_k(rep(1, 20), rep(0, 20), rep(1, 20), rep(1.1, 20),
                       g_fixed = ashr::normalmix(1, 0, 1))
  assert_near(ans$second - ans$mean^2, 1 / (1 + 0.5 * 22))
  assert_true(ans$second > ans$mean^2)
  assert_near(ans$second - (-ans$mean)^2, ans$sd^2)
})
