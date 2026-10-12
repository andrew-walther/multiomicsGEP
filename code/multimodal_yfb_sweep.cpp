// ============================================================
// Script: multimodal_yfb_sweep.cpp
// Purpose: C++ version of the per-feature Gauss-Seidel loading sweep in
//          multimodal_yfb_update_F_mk() (the runtime bottleneck). It mirrors
//          the R code line by line; tests check agreement to 1e-10.
// Author: Andrew Walther
// Created: 2026-10-01
// Dependencies: Rcpp
// ============================================================
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Moments of N(location, s^2) truncated to [0, inf); inverse-Mills expansion
// for h = location / s < -8 (see multimodal_yfb_positive_truncnorm_moments()).
static void truncnorm_pos(double location, double s, double &mean, double &second) {
  double h = location / s;
  if (h < -8) {
    double a = -h, ia = 1.0 / a;
    double delta = ia - 2 * std::pow(ia, 3) + 10 * std::pow(ia, 5) - 74 * std::pow(ia, 7);
    mean = s * delta;
    second = s * s * (1 - a * delta);
  } else {
    double mills = std::exp(R::dnorm(h, 0.0, 1.0, 1) - R::pnorm(h, 0.0, 1.0, 1, 1));
    mean = location + s * mills;
    second = location * location + s * s + location * s * mills;
  }
}

static double logsumexp2(double a, double b) {
  double m = std::max(a, b);
  if (!std::isfinite(m)) return m;
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

// family: 0 point_exponential, 1 point_laplace, 2 normal
struct Post { double mean, second, slab_prob, x, s2, log_ml; };

static Post posterior(double A, double B, int family, bool point_mass,
                      double pi, double rate, double variance) {
  Post p;
  if (point_mass) {
    p.mean = 0; p.second = 0; p.slab_prob = 0; p.log_ml = NA_REAL;
    if (family == 0) { p.x = NA_REAL; p.s2 = NA_REAL; }
    else { p.x = A > 0 ? B / A : NA_REAL; p.s2 = A > 0 ? 1 / A : R_PosInf; }
    return p;
  }
  if (A == 0) {
    if (std::fabs(B) <= std::sqrt(DBL_MIN)) B = 0;
    if (B != 0) stop("A = 0 with B != 0 is not a valid quadratic update.");
    p.x = NA_REAL; p.s2 = R_PosInf; p.log_ml = NA_REAL;
    if (family == 0) { p.mean = pi / rate; p.second = 2 * pi / (rate * rate); p.slab_prob = pi; }
    else if (family == 2) { p.mean = 0; p.second = variance; p.slab_prob = 1; }
    else { p.mean = 0; p.second = 2 * pi / (rate * rate); p.slab_prob = pi; }
    return p;
  }
  double x = B / A, s2 = 1 / A, s = std::sqrt(s2);
  p.x = x; p.s2 = s2;
  if (family == 2) {
    double v = 1 / (A + 1 / variance), m = v * B;
    p.mean = m; p.second = v + m * m; p.slab_prob = 1;
    p.log_ml = R::dnorm(x, 0.0, std::sqrt(s2 + variance), 1);
    return p;
  }
  double log_spike = pi < 1 ? std::log1p(-pi) + R::dnorm(x, 0.0, s, 1) : R_NegInf;
  if (family == 0) {
    double loc = x - rate * s2;
    double log_slab = std::log(pi) + std::log(rate) - rate * x + 0.5 * rate * rate * s2 +
      R::pnorm(loc / s, 0.0, 1.0, 1, 1);
    double log_total = logsumexp2(log_spike, log_slab);
    double slab_prob = std::exp(log_slab - log_total);
    double sm, ss; truncnorm_pos(loc, s, sm, ss);
    p.mean = slab_prob * sm;
    p.second = std::max(slab_prob * ss, p.mean * p.mean);
    p.slab_prob = slab_prob; p.log_ml = log_total;
    return p;
  }
  // point-Laplace
  double lam = rate;
  double c = std::log(pi) + std::log(lam / 2) + 0.5 * lam * lam * s2;
  double log_pos = c - lam * x + R::pnorm((x - lam * s2) / s, 0.0, 1.0, 1, 1);
  double log_neg = c + lam * x + R::pnorm((-x - lam * s2) / s, 0.0, 1.0, 1, 1);
  double mx = std::max(log_spike, std::max(log_pos, log_neg));
  double log_ml = mx + std::log(std::exp(log_spike - mx) + std::exp(log_pos - mx) + std::exp(log_neg - mx));
  double w_spike = std::exp(log_spike - log_ml), w_pos = std::exp(log_pos - log_ml),
         w_neg = std::exp(log_neg - log_ml);
  double pm, ps, nm, ns;
  truncnorm_pos(x - lam * s2, s, pm, ps);
  truncnorm_pos(-x - lam * s2, s, nm, ns);
  p.mean = w_pos * pm - w_neg * nm;
  p.second = std::max(w_pos * ps + w_neg * ns, p.mean * p.mean);
  p.slab_prob = 1 - w_spike; p.log_ml = log_ml;
  return p;
}

// [[Rcpp::export]]
List mmyfb_F_sweep_cpp(NumericMatrix Y, NumericVector Tau, double sum_EL2,
                       NumericVector sum_wy2, NumericVector sum_LR, NumericVector sum_yh,
                       NumericVector w, NumericVector EZ, NumericVector VZ,
                       NumericVector EF, NumericVector EF2, double EBeta, double EBeta2,
                       int family, bool point_mass, double pi, double rate, double variance) {
  int n = Y.nrow(), p = Y.ncol();
  NumericVector EF_new = clone(EF), EF2_new = clone(EF2), EZ_new = clone(EZ), VZ_new = clone(VZ);
  for (int i = 0; i < n; i++) if (VZ_new[i] < 0) VZ_new[i] = 0;
  NumericVector A_out(p), B_out(p), x_out(p), s2_out(p), mean_out(p), second_out(p),
                slab_out(p), logml_out(p);
  for (int j = 0; j < p; j++) {
    double old_mean = EF_new[j];
    double old_var = std::max(EF2_new[j] - old_mean * old_mean, 0.0);
    double wyEZ = 0;
    for (int i = 0; i < n; i++) wyEZ += w[i] * Y(i, j) * EZ_new[i];
    double sum_wyEZ_without_j = wyEZ - old_mean * sum_wy2[j];
    double A = Tau[j] * sum_EL2 + EBeta2 * sum_wy2[j];
    double B = Tau[j] * sum_LR[j] + EBeta * sum_yh[j] - EBeta2 * sum_wyEZ_without_j;
    Post q = posterior(A, B, family, point_mass, pi, rate, variance);
    double new_var = std::max(q.second - q.mean * q.mean, 0.0);
    EF_new[j] = q.mean;
    EF2_new[j] = q.mean * q.mean + new_var;
    for (int i = 0; i < n; i++) {
      double y = Y(i, j);
      EZ_new[i] += y * (q.mean - old_mean);
      VZ_new[i] += y * y * (new_var - old_var);
    }
    A_out[j] = A; B_out[j] = B; x_out[j] = q.x; s2_out[j] = q.s2; mean_out[j] = q.mean;
    second_out[j] = q.second; slab_out[j] = q.slab_prob; logml_out[j] = q.log_ml;
  }
  for (int i = 0; i < n; i++) if (VZ_new[i] < 0) VZ_new[i] = 0;
  return List::create(_["EF"] = EF_new, _["EF2"] = EF2_new, _["EZ"] = EZ_new, _["VZ"] = VZ_new,
                      _["A"] = A_out, _["B"] = B_out, _["x"] = x_out, _["s2"] = s2_out,
                      _["mean"] = mean_out, _["second"] = second_out,
                      _["slab_prob"] = slab_out, _["log_ml"] = logml_out);
}

// Vectorized posterior for independent coordinates (used by the score update,
// where every coordinate shares the same A): same formulas as posterior().
// [[Rcpp::export]]
List mmyfb_posterior_vec_cpp(NumericVector A, NumericVector B, int family, bool point_mass,
                             double pi, double rate, double variance) {
  int n = B.size();
  NumericVector mean(n), second(n), slab(n), x(n), s2(n), log_ml(n);
  for (int i = 0; i < n; i++) {
    Post q = posterior(A[A.size() == 1 ? 0 : i], B[i], family, point_mass, pi, rate, variance);
    mean[i] = q.mean; second[i] = q.second; slab[i] = q.slab_prob;
    x[i] = q.x; s2[i] = q.s2; log_ml[i] = q.log_ml;
  }
  return List::create(_["mean"] = mean, _["second"] = second, _["slab_prob"] = slab,
                      _["x"] = x, _["s2"] = s2, _["log_ml"] = log_ml);
}
