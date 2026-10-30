# Behavioral checks for the opt-in normalized-projection engine. ----
source("code/fit_yfb_supervised.R")
run_test("YFBS: exact Normal KL from ebnm output", {
  g <- ashr::normalmix(1,0,2)
  q <- yfb_ebnm(c(-1,2), c(2,3), "normal", g)
  v <- q$second-q$mean^2
  exact <- sum(.5*((v+q$mean^2)/4-1+log(4/v)))
  assert_near(q$kl,exact,1e-8)
})
run_test("YFBS: normalized projection is scale invariant", {
  set.seed(19); Y<-matrix(rexp(120),20); F<-matrix(rexp(18),6)
  z1<-yfb_projection(scale(Y,scale=FALSE),F,F^2+.01)
  cc<-c(.2,3,5); G<-sweep(F,2,cc,"*")
  z2<-yfb_projection(scale(Y,scale=FALSE),G,sweep(F^2+.01,2,cc^2,"*"))
  assert_near(z1$Z,z2$Z,1e-7);assert_near(z1$V,z2$V,1e-7)
  b<-c(.1,-.4,.3);assert_near(drop(z1$Z%*%b),drop(z2$Z%*%b),1e-7)
})
run_test("YFBS: bounded point-normal remains nondegenerate under null signal", {
  g<-yfb_beta_prior(rep(0,7),rep(1,7),7)
  assert_true(g$pi[2]>=1/7 && g$pi[2]<=.5 && g$sd[2]>=.05)
  q<-yfb_ebnm(rep(0,7),rep(1,7),"point_normal",g)
  assert_positive(q$second-q$mean^2)
})
set.seed(17)
ys <- matrix(rexp(40*12),40,12)
ts <- rexp(40); ds <- rep(c(0,1),20)
run_test("YFBS: alpha zero equals explicitly removed survival message", {
  a<-fit_yfb_supervised(list(expression=ys),ts,ds,K=3,alpha=0,max_iter=8)
  b<-fit_yfb_supervised(list(expression=ys),ts,ds,K=3,alpha=.5,gamma=0,max_iter=8)
  assert_near(a$EL,b$EL,1e-12);assert_near(a$EF[[1]],b$EF[[1]],1e-12)
  assert_near(a$history$objective,b$history$objective,1e-12)
})
run_test("YFBS: objective decreases are explicitly reported", {
  a<-fit_yfb_supervised(list(expression=ys),ts,ds,K=3,gamma=10,max_iter=8)
  decrease<-c(FALSE,diff(a$history$objective)< -1e-8*pmax(1,abs(head(a$history$objective,-1))))
  assert_equal(a$history$objective_decreased,decrease)
  assert_true(!a$converged || !tail(decrease,1))
  assert_finite(a$eta)
})
run_test("YFBS: cap and degenerate input cannot be silent convergence", {
  a<-fit_yfb_supervised(list(expression=ys),ts,ds,K=3,max_iter=1)
  assert_false(a$converged);assert_true("max_iter" %in% a$flags)
  err<-try(fit_yfb_supervised(list(expression=ys*0),ts,ds,K=3),silent=TRUE)
  assert_true(inherits(err,"try-error"))
  err<-try(yfb_cox(ts,ds,rep(NaN,40)),silent=TRUE)
  assert_true(inherits(err,"try-error"))
})
run_test("YFBS: collapsed point-exponential columns retain exact zero KL", {
  set.seed(4)
  x<-matrix(rexp(400*50),400,50)
  a<-fit_yfb_supervised(list(expression=x),rexp(400),rep(c(0,1),200),K=7,
                         gamma=1,max_iter=3)
  assert_finite(a$EF[[1]]);assert_true(all(a$EF2[[1]]>=a$EF[[1]]^2-1e-10))
})
