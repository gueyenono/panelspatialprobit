// Exact univariate truncated normal sampling using R's RNG.
//
// Algorithms: naive normal rejection when the support holds enough mass,
// Robert (1995) exponential rejection for one-sided tails, and Robert's
// uniform rejection for short two-sided intervals.
#ifndef PSP_TRUNCNORM_H
#define PSP_TRUNCNORM_H

#include <Rcpp.h>
#include <cmath>
#include <limits>

namespace psp {

// Standard normal truncated to [a, inf), a finite.
inline double rtnorm_std_lower(double a) {
  if (a <= 0.0) {
    // acceptance probability >= 0.5
    for (;;) {
      double x = norm_rand();
      if (x >= a) return x;
    }
  }
  // Robert (1995): translated exponential proposal with optimal rate
  const double lambda = 0.5 * (a + std::sqrt(a * a + 4.0));
  for (;;) {
    double x = a + exp_rand() / lambda;
    double d = x - lambda;
    if (unif_rand() <= std::exp(-0.5 * d * d)) return x;
  }
}

// Standard normal truncated to [a, b], both finite, a < b.
inline double rtnorm_std_two_sided(double a, double b) {
  if (a < 0.0 && b > 0.0) {
    if (b - a >= 2.5066282746310002) {  // sqrt(2 pi): normal rejection is efficient
      for (;;) {
        double x = norm_rand();
        if (x >= a && x <= b) return x;
      }
    }
    for (;;) {  // uniform rejection, mode inside the interval
      double x = a + (b - a) * unif_rand();
      if (unif_rand() <= std::exp(-0.5 * x * x)) return x;
    }
  }
  // interval on one side of zero: work with a >= 0 and mirror if needed
  bool flip = false;
  if (b <= 0.0) {
    double tmp = a;
    a = -b;
    b = -tmp;
    flip = true;
  }
  // Robert (1995) rule for choosing uniform versus exponential proposal
  const double sq = std::sqrt(a * a + 4.0);
  const double bound = 2.0 * std::sqrt(M_E) / (a + sq) *
                       std::exp(0.25 * (a * a - a * sq));
  double x;
  if (b - a <= bound) {
    for (;;) {
      x = a + (b - a) * unif_rand();
      if (unif_rand() <= std::exp(0.5 * (a * a - x * x))) break;
    }
  } else {
    for (;;) {
      x = rtnorm_std_lower(a);
      if (x <= b) break;
    }
  }
  return flip ? -x : x;
}

// N(mean, sd^2) truncated to [lo, hi]; lo/hi may be infinite.
inline double rtnorm(double mean, double sd, double lo, double hi) {
  const double a = (lo - mean) / sd;
  const double b = (hi - mean) / sd;
  const bool fin_a = std::isfinite(a), fin_b = std::isfinite(b);
  double x;
  if (!fin_a && !fin_b) {
    x = norm_rand();
  } else if (fin_a && !fin_b) {
    x = rtnorm_std_lower(a);
  } else if (!fin_a && fin_b) {
    x = -rtnorm_std_lower(-b);
  } else {
    if (b <= a) return mean + sd * a;  // degenerate interval
    x = rtnorm_std_two_sided(a, b);
  }
  return mean + sd * x;
}

}  // namespace psp

#endif
