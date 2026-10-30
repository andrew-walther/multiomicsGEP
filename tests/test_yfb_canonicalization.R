# Canonicalization must preserve the same predictor and posterior prior units. ----
source("code/fit_yfb_supervised.R")
run_test("YFBS: tiny live activities are normalized rather than skipped", local({
  yc<-matrix(c(1,-1,0,2,-2,0),3,2)
  L<-matrix(rep(.5e-12,3),3);L2<-matrix(rep(1e-24,3),3)
  F<-matrix(c(3e12,4e12),2);F2<-F^2+1e24
  g<-ebnm::ebnm_point_exponential(c(1,3),1)$fitted_g
  g$pi<-c(.5,.5);g$scale<-c(0,1)
  st<-list(L=L,L2=L2,F=F,F2=F2,gL=list(g),gF=list(g),
    b=.3,b2=.1,gB=ashr::normalmix(c(.5,.5),c(0,0),c(0,1)))
  before<-yfb_projection(yc,F,F2)
  q<-yfb_canonicalize(st)
  after<-yfb_projection(yc,q$F,q$F2)
  assert_near(colMeans(q$L2),1,1e-12)
  assert_near(q$L%*%t(q$F),L%*%t(F),1e-12)
  assert_near(before$Z,after$Z,1e-12);assert_near(before$V,after$V,1e-12)
  assert_equal(q$b,st$b);assert_equal(q$gB,st$gB)
}))
run_test("YFBS: zero-activity loading has an invariant finite norm and KL", local({
  g<-ebnm::ebnm_point_exponential(c(1,3),1)$fitted_g
  g$pi<-c(.6,.4);g$scale<-c(0,2)
  x<-c(1,3);A<-c(2,4);post<-yfb_ebnm(x,A,"point_exponential",g)
  yc<-matrix(c(1,-1,0,2,-2,0),3,2)
  st<-list(L=matrix(0,3,1),L2=matrix(0,3,1),F=matrix(post$mean,2),
    F2=matrix(post$second,2),gL=list(g),gF=list(g),b=.3,b2=.1,gB=g)
  before<-yfb_projection(yc,st$F,st$F2);q<-yfb_canonicalize(st)
  after<-yfb_projection(yc,q$F,q$F2)
  r<-sqrt(sum(st$F^2))
  check<-yfb_ebnm(x/r,A*r^2,"point_exponential",q$gF[[1]])
  assert_near(sum(q$F^2),1,1e-12)
  assert_near(post$kl,check$kl,1e-10)
  assert_near(q$F,check$mean,1e-10)
  assert_near(before$Z,after$Z,1e-12);assert_near(before$V,after$V,1e-12)
  assert_equal(q$b,st$b);assert_equal(q$gB,st$gB)
}))
run_test("YFBS: live projection normalization has no positive scale floor", local({
  yc<-matrix(c(1,-1,0,2,-2,0),3,2)
  F<-matrix(c(3,4),2);F2<-F^2+1
  a<-yfb_projection(yc,F,F2)
  b<-yfb_projection(yc,F*1e-12,F2*1e-24)
  assert_near(a$Z,b$Z,1e-12);assert_near(a$V,b$V,1e-12)
}))
run_test("YFBS: reconstruction-only fitting precedes every coefficient update", local({
  old<-yfb_ebnm;families<-character()
  assign("yfb_ebnm",function(x,A,family,g=NULL) {
    families<<-c(families,family);old(x,A,family,g)
  },envir=.GlobalEnv)
  set.seed(68);X<-matrix(rexp(40*12),40,12)
  ans<-tryCatch(fit_yfb_supervised(list(expression=X),rexp(40),rep(c(0,1),20),
    K=3,gamma=0,max_iter=1),finally=assign("yfb_ebnm",old,envir=.GlobalEnv))
  assert_equal(families[1],"point_exponential")
}))
