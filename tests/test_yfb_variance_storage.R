# Separate variances must remain representable when second moments round. ----
source("code/fit_yfb_supervised.R")
run_test("YFBS: strong point-normal messages retain variance and exact KL", local({
  g<-ashr::normalmix(c(.8,.2),c(0,0),c(0,1e5))
  q<-yfb_ebnm(5e4,1e8,"point_normal",g)
  v<-1/(1e8+1e-10);m<-v*1e8*5e4
  exact<--log(.2)+.5*(v/1e10+m^2/1e10-1+log(1e10/v))
  assert_length(q$variance,1);assert_true(q$variance>0)
  assert_near(q$variance,v,1e-18);assert_near(q$mean,m,1e-10)
  assert_near(q$kl,exact,1e-8)
}))
run_test("YFBS: strong exponential messages retain separate variance", local({
  q<-yfb_stable_exp(1e12,1e20,.4,1)
  assert_length(q$variance,1);assert_true(q$variance>0)
  assert_near(q$variance,1e-12,1e-20)
  exact<-1e8-log(.4)-.5*log(2*pi*exp(1)/1e12)
  assert_near(q$kl,exact,1e-6)
}))
run_test("YFBS: predictor variance does not subtract large second moments", {
  out<-yfb_predictor_variance(matrix(c(1,2),2),matrix(c(.1,.2),2),1e12,1e-8)
  assert_near(out,c(1e23,2e23),1e9)
  assert_true(all(out>=0))
  tiny<-yfb_predictor_variance(matrix(c(1,2),2),matrix(0,2),1e12,1e-8)
  assert_near(tiny,c(1e-8,4e-8),1e-18)
})
