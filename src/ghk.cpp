// GHK simulator (Geweke-Hajivassiliou-Keane) for rectangle probabilities of a
// Gaussian vector given by its sparse precision (plan/02, approximate
// likelihood). For x ~ N(0, Q^{-1}) with Q = L L' (L lower triangular, the
// coordinates already in the factor's order), L'x ~ N(0, I), so the
// coordinates can be drawn one at a time from the last to the first.
#include <RcppEigen.h>
#include <cmath>
#include "truncnorm.h"

// [[Rcpp::depends(RcppEigen)]]

namespace {

// log(Phi(b) - Phi(a)) for a < b, accurate in both tails
double log_interval(double a, double b) {
  if (a > 0.0) {  // upper tail: Q(a) - Q(b)
    const double la = R::pnorm(a, 0.0, 1.0, 0, 1), lb = R::pnorm(b, 0.0, 1.0, 0, 1);
    return la + std::log1p(-std::exp(lb - la));
  }
  if (b < 0.0) {
    const double la = R::pnorm(a, 0.0, 1.0, 1, 1), lb = R::pnorm(b, 0.0, 1.0, 1, 1);
    return lb + std::log1p(-std::exp(la - lb));
  }
  return std::log(R::pnorm(b, 0.0, 1.0, 1, 0) - R::pnorm(a, 0.0, 1.0, 1, 0));
}

}  // namespace

// Returns the log importance weights of 'reps' GHK replications (rows), split
// by the independent blocks of coordinates given in 'block' (0-based; columns).
// The probability estimate of a block is the mean of the exponentials of its
// column, and independent blocks multiply.
// [[Rcpp::export]]
Rcpp::NumericMatrix psp_ghk_cpp(SEXP L_sexp, const Eigen::VectorXd& lo,
                                const Eigen::VectorXd& hi, const Eigen::VectorXi& block,
                                int nblocks, int reps) {
  typedef Eigen::Map<Eigen::SparseMatrix<double>> MapSp;
  const MapSp L(Rcpp::as<MapSp>(L_sexp));
  const int n = L.cols();
  Rcpp::NumericMatrix logw(reps, nblocks);
  Eigen::VectorXd x(n);
  for (int r = 0; r < reps; ++r) {
    for (int i = n - 1; i >= 0; --i) {
      double d = 0.0, s = 0.0;
      for (MapSp::InnerIterator it(L, i); it; ++it) {
        if (it.row() == i) d = it.value();
        else s += it.value() * x[it.row()];
      }
      // (L'x)_i = d x_i + s = e_i ~ N(0, 1), with lo_i < x_i < hi_i
      const double a = std::isinf(lo[i]) ? lo[i] : d * lo[i] + s;
      const double b = std::isinf(hi[i]) ? hi[i] : d * hi[i] + s;
      if (!(std::isinf(a) && std::isinf(b))) logw(r, block[i]) += log_interval(a, b);
      x[i] = (psp::rtnorm(0.0, 1.0, a, b) - s) / d;
    }
  }
  return logw;
}
