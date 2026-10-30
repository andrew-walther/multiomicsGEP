// ============================================================
// Script: yfb_exponential_stable.cpp
// Purpose: Stable point-exponential posterior and opt-in coupled loading sweep.
// Created: 2026-10-30
// Dependencies: Rcpp; used only by the experimental normalized-projection fitter
// ============================================================
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
struct ExpPost {double mean, second, variance, kl;};

// Cancel N(x;0,1/A) before computing spike/slab odds or KL. The slab
// integral is I=integral λ exp((B-λ)t - At²/2)dt, t>=0. This avoids
// subtracting large Gaussian terms when A is tiny and x=B/A is huge.
static ExpPost exp_post(double A,double B,double pi,double rate) {
  if(pi==0) return {0,0,0,0};
  if(!(A>0) || !std::isfinite(A+B+pi+rate) || pi<0 || pi>1 || rate<=0)
    stop("Invalid stable exponential message");
  double se=1/std::sqrt(A), h=(B-rate)*se, sm, ss, logI;
  if(h < -20) {
    double a=-h, v=1/(a*a);
    double delta=(1/a)*(1+v*(-2+v*(10+v*(-74+v*(706+v*(-8162+v*110410))))));
    sm=se*delta;
    ss=se*se*v*(2+v*(-10+v*(74+v*(-706+v*(8162-v*110410)))));
    double series=1+v*(-1+v*(3+v*(-15+v*(105+v*(-945+v*(10395-v*135135))))));
    logI=std::log(rate)+std::log(se)-std::log(a)+std::log(series);
  } else {
    double mills=std::exp(R::dnorm(h,0,1,1)-R::pnorm(h,0,1,1,1));
    sm=se*(h+mills);ss=se*se*(1+h*h+h*mills);
    logI=std::log(rate)+std::log(se)+.5*std::log(2*M_PI)+.5*h*h+
      R::pnorm(h,0,1,1,1);
  }
  double l0=pi<1 ? std::log1p(-pi) : R_NegInf;
  double l1=std::log(pi)+logI, high=std::max(l0,l1);
  double logN=high+std::log(std::exp(l0-high)+std::exp(l1-high));
  double prob=std::exp(l1-logN), mean=prob*sm;
  double slabvar;
  if(h < -20) slabvar=ss-sm*sm;
  else {
    double mills=std::exp(R::dnorm(h,0,1,1)-R::pnorm(h,0,1,1,1));
    slabvar=se*se*(1-h*mills-mills*mills);
  }
  double variance=prob*slabvar+prob*(1-prob)*sm*sm;
  double second=mean*mean+variance;
  double kl=B*mean-.5*A*second-logN;
  // For a well-separated positive Gaussian slab, use its entropy form;
  // the quadratic terms otherwise cancel to machine zero at huge A.
  if(h>20 && prob==1) kl=-std::log(pi)-std::log(rate)+rate*mean-
    .5*std::log(2*M_PI*std::exp(1.0)/A);
  if(!std::isfinite(mean+second+kl) || second<mean*mean-1e-8*std::max(1.0,second))
    stop("Invalid stable exponential posterior");
  return {mean,second,variance,kl};
}

// [[Rcpp::export]]
List yfb_exp_vector_cpp(NumericVector A,NumericVector B,double pi,double rate) {
  int p=A.size();NumericVector mean(p),second(p),variance(p),kl(p);
  for(int j=0;j<p;++j) {auto q=exp_post(A[j],B[j],pi,rate);
    mean[j]=q.mean;second[j]=q.second;variance[j]=q.variance;kl[j]=q.kl;}
  return List::create(_["mean"]=mean,_["second"]=second,_["variance"]=variance,_["kl"]=kl);
}

// [[Rcpp::export]]
List yfb_exp_sweep_cpp(NumericMatrix Y,NumericVector tau,double sumL2,
 NumericVector wy2,NumericVector LR,NumericVector yh,NumericVector w,
 NumericVector Z0,NumericVector V0,NumericVector F0,NumericVector F20,NumericVector VF0,
 double beta,double beta2,double pi,double rate) {
 int n=Y.nrow(),p=Y.ncol();
 NumericVector Z=clone(Z0),V=clone(V0),F=clone(F0),F2=clone(F20),VF=clone(VF0);
 NumericVector x(p),s2(p),A(p),B(p),kl(p);
 for(int j=0;j<p;++j) {
   double coupling=0;
   for(int i=0;i<n;++i)coupling+=w[i]*Y(i,j)*Z[i];
   coupling-=F[j]*wy2[j];
   A[j]=tau[j]*sumL2+beta2*wy2[j];
   B[j]=tau[j]*LR[j]+beta*yh[j]-beta2*coupling;
   auto q=exp_post(A[j],B[j],pi,rate);
   double dm=q.mean-F[j],dv=q.variance-VF[j];
   for(int i=0;i<n;++i) {Z[i]+=Y(i,j)*dm;V[i]+=Y(i,j)*Y(i,j)*dv;}
   F[j]=q.mean;F2[j]=q.second;VF[j]=q.variance;kl[j]=q.kl;
   x[j]=A[j]>0 ? B[j]/A[j] : NA_REAL;s2[j]=A[j]>0 ? 1/A[j] : R_PosInf;
 }
 return List::create(_["EF"]=F,_["EF2"]=F2,_["VF"]=VF,_["EZ"]=Z,_["VZ"]=V,
   _["x"]=x,_["s2"]=s2,_["A"]=A,_["B"]=B,_["kl"]=kl);
}
