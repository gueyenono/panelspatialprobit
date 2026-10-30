// MCMC sampler for panel spatial (ordered) probit models.
//
// Latent model, observations stacked with time slow and unit fast:
//   R_t S_t z_t = R_t (X_t beta + alpha) + eps_t,  eps_t ~ N(0, I_N),
//   S_t = I - rho W_t (spatial lag), R_t = I - lambda M_t (spatial error),
// i.e. the kernel (A z - b)'(A z - b) with block-diagonal A = R S (plan/02
// s2.2). Weights are time-invariant (one matrix) or given per period.
// Spatial random effects: alpha = phi W alpha + a, a ~ N(0, sigma_alpha^2 I).
// Latent AR(1) dynamics: gamma z_{t-1} enters period t >= 2 as the design
// column 'lag_col' (refreshed from z), so gamma is drawn with beta. The full
// A is then block lower-bidiagonal with blocks A_t and -gamma R_t below the
// diagonal; log|A| is unchanged.
// Student-t errors (fixed nu) are the scale mixture eps_it ~ N(0, v_it),
// v_it ~ IG(nu/2, nu/2): the kernel becomes sum_i w_i r_i^2 with w = 1/v.
// One call runs one chain. Block order per iteration (see plan/02 s2.1):
//   z | . -> rho | z, alpha, lambda (beta integrated out) -> beta | .
//   -> lambda | . -> alpha | . -> phi | alpha, sigma_alpha -> sigma_alpha^2 | alpha
//   -> cut-points | . -> scale move
#include <RcppEigen.h>
#include <cmath>
#include <limits>
#include <vector>
#include "truncnorm.h"

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<SpMat> MapSp;

namespace {

const double INF = std::numeric_limits<double>::infinity();

struct Model {
  int N, T, n, k, J;
  bool ordered, has_rho, has_lambda, has_re, spatial_re;
};

// Per-period sparse operators: one matrix (time-invariant) or one per period
typedef std::vector<SpMat> Ops;

inline const SpMat& at(const Ops& ops, int t) { return ops[ops.size() == 1 ? 0 : t]; }

// Sparse operator applied period by period: out_t = Op_t v_t
void apply_op(const Ops& ops, const VectorXd& v, VectorXd& out, int N, int T) {
  for (int t = 0; t < T; ++t) {
    out.segment(t * N, N).noalias() = at(ops, t) * v.segment(t * N, N);
  }
}

// e = S z - X beta - Delta alpha and r = R e = A z - b
void compute_resid(const Model& m, const Ops& W, const Ops& M, const MatrixXd& X,
                   const VectorXd& z, const VectorXd& beta, const VectorXd& alpha,
                   double rho, double lambda, VectorXd& work, VectorXd& e, VectorXd& r) {
  e.noalias() = z - X * beta;
  if (m.has_rho) {
    apply_op(W, z, work, m.N, m.T);
    e -= rho * work;
  }
  if (m.has_re) {
    for (int t = 0; t < m.T; ++t) e.segment(t * m.N, m.N) -= alpha;
  }
  r = e;
  if (m.has_lambda) {
    apply_op(M, e, work, m.N, m.T);
    r -= lambda * work;
  }
}

// A_t = (I - lambda M_t)(I - rho W_t), with the column sums of squares of
// each A_t (columns of 'colsq'). One A when the weights do not vary.
void build_A(const Model& m, const Ops& W, const Ops& M, double rho, double lambda, Ops& A,
             MatrixXd& colsq) {
  const bool varies = (m.has_rho && W.size() > 1) || (m.has_lambda && M.size() > 1);
  const int nA = varies ? m.T : 1;
  A.resize(nA);
  colsq.setZero(m.N, nA);
  SpMat I(m.N, m.N);
  I.setIdentity();
  for (int a = 0; a < nA; ++a) {
    SpMat S = I;
    if (m.has_rho) S = I - rho * at(W, a);
    if (m.has_lambda) {
      const SpMat R = I - lambda * at(M, a);
      A[a] = R * S;
    } else {
      A[a] = S;
    }
    A[a].makeCompressed();
    for (int c = 0; c < m.N; ++c)
      for (SpMat::InnerIterator it(A[a], c); it; ++it) colsq(c, a) += it.value() * it.value();
  }
}

// Truncation bounds of z_idx. Missing outcomes (code 0 for ordered, -1 for
// binary) are unrestricted, which imputes them (plan/02 s2.3).
void bounds(const Model& m, const Eigen::VectorXi& y, const std::vector<double>& kappa,
            int idx, double& lo, double& hi) {
  if (m.ordered) {
    if (y[idx] == 0) {
      lo = -INF;
      hi = INF;
    } else {
      lo = kappa[y[idx] - 1];
      hi = kappa[y[idx]];
    }
  } else if (y[idx] == 1) {
    lo = 0.0;
    hi = INF;
  } else if (y[idx] == 0) {
    lo = -INF;
    hi = 0.0;
  } else {
    lo = -INF;
    hi = INF;
  }
}

// One Gibbs sweep over all latent values using residual updates (plan/02 s2.3):
// for coordinate k with column a_k of A, z_k | . ~ TN(z_k - a_k'Dr / h_k, 1 / h_k)
// with h_k = a_k'D a_k, then r += a_k (z_k' - z_k). D = diag(wt) for t errors
// ('wt' null: D = I and h_k comes from 'colsqs').
// With latent AR (gamma != 0) z_{k,t} also enters period t + 1 through column
// k of -gamma R_{t+1} ('Rs'; identity when empty).
void latent_sweep(const Model& m, const Ops& As, const MatrixXd& colsqs,
                  const Eigen::VectorXi& y, const std::vector<double>& kappa,
                  const VectorXd* wt, double gamma, const Ops& Rs, VectorXd& z,
                  VectorXd& r) {
  for (int t = 0; t < m.T; ++t) {
    const int off = t * m.N;
    const SpMat& A = at(As, t);
    const auto colsq = colsqs.col(As.size() == 1 ? 0 : t);
    for (int k = 0; k < m.N; ++k) {
      const int idx = off + k;
      double ar = 0.0, h = 0.0;
      if (wt) {
        for (SpMat::InnerIterator it(A, k); it; ++it) {
          const double wa = (*wt)[off + it.row()] * it.value();
          ar += wa * r[off + it.row()];
          h += wa * it.value();
        }
      } else {
        for (SpMat::InnerIterator it(A, k); it; ++it) ar += it.value() * r[off + it.row()];
        h = colsq[k];
      }
      const bool lagged = gamma != 0.0 && t < m.T - 1;
      const int off1 = off + m.N;
      if (lagged) {
        if (Rs.empty()) {
          const double wv = wt ? (*wt)[off1 + k] : 1.0;
          ar -= gamma * wv * r[off1 + k];
          h += gamma * gamma * wv;
        } else {
          for (SpMat::InnerIterator it(at(Rs, t + 1), k); it; ++it) {
            const double wv = wt ? (*wt)[off1 + it.row()] : 1.0;
            ar -= gamma * wv * it.value() * r[off1 + it.row()];
            h += gamma * gamma * wv * it.value() * it.value();
          }
        }
      }
      const double mean = z[idx] - ar / h;
      double lo, hi;
      bounds(m, y, kappa, idx, lo, hi);
      const double znew = psp::rtnorm(mean, 1.0 / std::sqrt(h), lo, hi);
      const double delta = znew - z[idx];
      z[idx] = znew;
      for (SpMat::InnerIterator it(A, k); it; ++it) r[off + it.row()] += it.value() * delta;
      if (lagged) {
        if (Rs.empty()) {
          r[off1 + k] -= gamma * delta;
        } else {
          for (SpMat::InnerIterator it(at(Rs, t + 1), k); it; ++it)
            r[off1 + it.row()] -= gamma * it.value() * delta;
        }
      }
    }
  }
}

double lndet_interp(const VectorXd& grid, const VectorXd& lndet, double x) {
  const int G = grid.size();
  const double step = (grid[G - 1] - grid[0]) / (G - 1);
  double pos = (x - grid[0]) / step;
  int i = static_cast<int>(std::floor(pos));
  if (i < 0) i = 0;
  if (i > G - 2) i = G - 2;
  const double w = pos - i;
  return (1.0 - w) * lndet[i] + w * lndet[i + 1];
}

// Draw from a density known on a grid: pick a cell by trapezoid mass, then
// sample exactly from the linear density inside it.
double draw_from_grid(const std::vector<double>& x, const std::vector<double>& lp) {
  const int G = x.size();
  double mx = -INF;
  for (int i = 0; i < G; ++i) if (lp[i] > mx) mx = lp[i];
  std::vector<double> f(G), cum(G - 1);
  for (int i = 0; i < G; ++i) f[i] = std::exp(lp[i] - mx);
  double tot = 0.0;
  for (int i = 0; i < G - 1; ++i) {
    tot += 0.5 * (f[i] + f[i + 1]) * (x[i + 1] - x[i]);
    cum[i] = tot;
  }
  const double u = unif_rand() * tot;
  int c = 0;
  while (c < G - 2 && cum[c] < u) ++c;
  const double prev = (c == 0) ? 0.0 : cum[c - 1];
  const double h = x[c + 1] - x[c];
  const double f0 = f[c], f1 = f[c + 1];
  const double mass = 0.5 * (f0 + f1) * h;
  const double v = (mass > 0.0) ? (u - prev) / mass : unif_rand();
  const double slope = (f1 - f0) / h;
  double s;
  if (std::fabs(slope) * h < 1e-10 * (f0 + f1)) {
    s = v * h;
  } else {
    s = (-f0 + std::sqrt(f0 * f0 + 2.0 * slope * v * mass)) / slope;
  }
  if (s < 0.0) s = 0.0;
  if (s > h) s = h;
  return x[c] + s;
}

// Collapsed griddy Gibbs draw of rho (plan/02 s2.4), two-stage grid: a coarse
// pass (every 'stride'-th grid point) locates the region within 30 log units
// of the mode, which is then refined on 400 points.
double draw_rho(const VectorXd& grid, const VectorXd& lndet, int T, double q00,
                double q01, double q11, double a, double b, double rmin, double rmax) {
  const int G = grid.size();
  const int stride = std::max(1, G / 500);
  std::vector<int> idx;
  for (int i = 0; i < G; i += stride) idx.push_back(i);
  if (idx.back() != G - 1) idx.push_back(G - 1);
  const int C = idx.size();
  std::vector<double> lp(C);
  double mx = -INF;
  for (int c = 0; c < C; ++c) {
    const double x = grid[idx[c]];
    lp[c] = T * lndet[idx[c]] - 0.5 * (q00 - 2.0 * x * q01 + x * x * q11) +
            (a - 1.0) * std::log(x - rmin) + (b - 1.0) * std::log(rmax - x);
    if (lp[c] > mx) mx = lp[c];
  }
  int first = C - 1, last = 0;
  for (int c = 0; c < C; ++c) {
    if (lp[c] > mx - 30.0) {
      if (c < first) first = c;
      if (c > last) last = c;
    }
  }
  first = idx[std::max(0, first - 2)];
  last = idx[std::min(C - 1, last + 2)];
  const int R = 400;
  std::vector<double> xs(R), lps(R);
  const double lo = grid[first], hi = grid[last];
  for (int i = 0; i < R; ++i) {
    const double x = lo + (hi - lo) * i / (R - 1.0);
    xs[i] = x;
    lps[i] = T * lndet_interp(grid, lndet, x) -
             0.5 * (q00 - 2.0 * x * q01 + x * x * q11) +
             (a - 1.0) * std::log(x - rmin) + (b - 1.0) * std::log(rmax - x);
  }
  return draw_from_grid(xs, lps);
}

double rinvgamma(double shape, double scale) { return scale / R::rgamma(shape, 1.0); }

double log_kappa_prior(double x, bool flat, double mean, double sd) {
  if (flat) return 0.0;
  const double d = (x - mean) / sd;
  return -0.5 * d * d;
}

// Albert-Chib Gibbs step for the free cut-points given z.
void cut_ac_step(const Model& m, const Eigen::VectorXi& y, const VectorXd& z,
                 std::vector<double>& kappa, bool flat, double mean, double sd) {
  for (int j = 2; j <= m.J - 1; ++j) {
    double lo = kappa[j - 1], hi = kappa[j + 1];
    for (int i = 0; i < m.n; ++i) {
      if (y[i] == j && z[i] > lo) lo = z[i];
      if (y[i] == j + 1 && z[i] < hi) hi = z[i];
    }
    if (hi > lo) {
      kappa[j] = flat ? lo + (hi - lo) * unif_rand() : psp::rtnorm(mean, sd, lo, hi);
    }
  }
}

// Joint (kappa_j, z) transformation move (plan/02 s2.7): kappa_j moves on the
// logit scale of (kappa_{j-1}, kappa_{j+1}) (log scale for the top cut-point)
// and the latent values of categories j and j + 1 are rescaled with it, so
// every z stays inside its category. Updates kappa, z and the residual
// r = A z - b in place on acceptance; returns acceptance.
bool cut_transform_move(const Model& m, const Ops& A, const Eigen::VectorXi& y,
                        const std::vector<int>& ncat, int j, double step, bool flat,
                        double kmean, double ksd, const VectorXd* wt, double gamma,
                        const Ops& Rs, std::vector<double>& kappa, VectorXd& z, VectorXd& r,
                        VectorXd& Ad) {
  const double L = kappa[j - 1], U = kappa[j + 1], kj = kappa[j];
  const bool top = (j == m.J - 1);
  const double eps = step * norm_rand();
  double knew, s1, s2 = 1.0, logjac;
  if (!top) {
    const double pp = (kj - L) / (U - L);
    const double lg = std::log(pp / (1.0 - pp)) + eps;
    const double pn = 1.0 / (1.0 + std::exp(-lg));
    if (pn <= 0.0 || pn >= 1.0) return false;
    knew = L + pn * (U - L);
    s1 = (knew - L) / (kj - L);
    s2 = (U - knew) / (U - kj);
    logjac = ncat[j] * std::log(s1) + ncat[j + 1] * std::log(s2) +
             std::log(pn * (1.0 - pn)) - std::log(pp * (1.0 - pp));
  } else {
    s1 = std::exp(eps);
    knew = L + (kj - L) * s1;
    logjac = (ncat[j] + 1.0) * eps;
  }
  VectorXd delta(m.n);
  for (int i = 0; i < m.n; ++i) {
    if (y[i] == j) {
      delta[i] = (L + (z[i] - L) * s1) - z[i];
    } else if (y[i] == j + 1) {
      delta[i] = top ? knew - kj : (U - (U - z[i]) * s2) - z[i];
    } else {
      delta[i] = 0.0;
    }
  }
  apply_op(A, delta, Ad, m.N, m.T);
  if (gamma != 0.0) {  // latent AR: -gamma R_t delta_{t-1} in period t
    for (int t = 1; t < m.T; ++t) {
      const auto dprev = delta.segment((t - 1) * m.N, m.N);
      if (Rs.empty()) Ad.segment(t * m.N, m.N) -= gamma * dprev;
      else Ad.segment(t * m.N, m.N) -= gamma * (at(Rs, t) * dprev);
    }
  }
  const double dq = wt ? (wt->cwiseProduct(2.0 * r + Ad)).dot(Ad)
                      : 2.0 * r.dot(Ad) + Ad.squaredNorm();
  const double logr = -0.5 * dq + logjac + log_kappa_prior(knew, flat, kmean, ksd) -
                      log_kappa_prior(kj, flat, kmean, ksd);
  if (std::log(unif_rand()) < logr) {
    z += delta;
    r += Ad;
    kappa[j] = knew;
    return true;
  }
  return false;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List psp_sampler_cpp(const Eigen::VectorXi& y, const Eigen::Map<Eigen::MatrixXd>& X,
                           Rcpp::List W_list, Rcpp::List M_list, int N, int T, int J,
                           bool ordered,
                           bool has_rho, bool has_lambda, bool has_re, bool spatial_re,
                           const Eigen::VectorXd& beta_mean, const Eigen::MatrixXd& beta_prec,
                           double rho_a, double rho_b, double rmin, double rmax,
                           const Eigen::VectorXd& rho_grid, const Eigen::VectorXd& lndet_grid,
                           double lam_a, double lam_b, double lmin, double lmax,
                           const Eigen::VectorXd& lam_grid, const Eigen::VectorXd& lndetM_grid,
                           double phi_a, double phi_b, double nu, int lag_col,
                           int sigma_prior, double sigma_p1, double sigma_p2, bool kappa_flat,
                           double kappa_mean, double kappa_sd, int iter, int warmup,
                           int thin, int z_sweeps, int cut_sampler, bool scale_move,
                           const Eigen::VectorXd& beta_start, double rho_start,
                           double lambda_start, double sigma2_start,
                           const Eigen::VectorXd& kappa_start, int alpha_every, bool progress,
                           int chain) {
  Ops W, M;
  for (int i = 0; i < W_list.size(); ++i) W.push_back(SpMat(Rcpp::as<MapSp>(W_list[i])));
  for (int i = 0; i < M_list.size(); ++i) M.push_back(SpMat(Rcpp::as<MapSp>(M_list[i])));
  Model m;
  m.N = N;
  m.T = T;
  m.n = N * T;
  m.k = X.cols();
  m.J = J;
  m.ordered = ordered;
  m.has_rho = has_rho;
  m.has_lambda = has_lambda;
  m.has_re = has_re;
  m.spatial_re = spatial_re;
  const int n = m.n, k = m.k;

  // state
  VectorXd beta = beta_start;
  double rho = has_rho ? rho_start : 0.0;
  double lambda = has_lambda ? lambda_start : 0.0;
  double sigma2 = sigma2_start;
  double phi = 0.0;
  // Student-t errors: precision weights w = 1 / v (all 1 for normal errors)
  const bool has_t = nu > 0.0;
  VectorXd wt = VectorXd::Ones(n), v_sum = VectorXd::Zero(n);
  const VectorXd* wtp = has_t ? &wt : nullptr;
  double xi = 1.0;
  VectorXd alpha = VectorXd::Zero(N);
  std::vector<double> kappa(J + 1);
  kappa[0] = -INF;
  kappa[J] = INF;
  if (ordered) {
    kappa[1] = 0.0;
    for (int j = 2; j < J; ++j) kappa[j] = kappa_start[j - 2];
  }
  std::vector<int> ncat(J + 1, 0);
  if (ordered) for (int i = 0; i < n; ++i) if (y[i] > 0) ncat[y[i]]++;

  // A = R S and its column norms; rebuilt whenever rho or lambda changes
  Ops A;
  MatrixXd colsq;
  build_A(m, W, M, rho, lambda, A, colsq);

  // latent AR: working copy of the design with the lag column, R_t = I -
  // lambda M_t for the coupling between periods (empty: identity)
  const bool latent_ar = lag_col >= 0;
  MatrixXd Xm = X;
  Ops Rs;
  auto build_R = [&]() {
    Rs.clear();
    if (!(latent_ar && has_lambda)) return;
    SpMat I(N, N);
    I.setIdentity();
    for (std::size_t i = 0; i < M.size(); ++i) Rs.push_back(SpMat(I - lambda * M[i]));
  };
  build_R();
  auto gamma_now = [&]() { return latent_ar ? beta[lag_col] : 0.0; };

  // Transformed design Xt = R X and the posterior covariance of beta. Without
  // an error term R = I and both are constant.
  MatrixXd Xt = Xm;
  MatrixXd B, LB;
  const VectorXd Pc = beta_prec * beta_mean;
  auto refresh_design = [&]() {
    if (has_lambda) {
      MatrixXd MX(n, k);
      for (int t = 0; t < T; ++t) MX.middleRows(t * N, N) = at(M, t) * Xm.middleRows(t * N, N);
      Xt = Xm - lambda * MX;
    } else if (latent_ar) {
      Xt = Xm;
    }
    B = has_t ? MatrixXd((Xt.transpose() * wt.asDiagonal() * Xt + beta_prec).inverse())
              : MatrixXd((Xt.transpose() * Xt + beta_prec).inverse());
    LB = B.llt().matrixL();
  };
  refresh_design();

  // sparse Cholesky for the random effects under an error term or with
  // spatial random effects (plan/02 s2.6): with B = I - phi W,
  // P = sum_t R_t'R_t + B'B / sigma2
  //   = T I - lambda sum_t (M_t + M_t') + lambda^2 sum_t M_t'M_t
  //     + (I - phi (W + W') + phi^2 W'W) / sigma2.
  // The pattern is fixed, so it is analysed once and refactorised per draw.
  const bool sparse_re = has_re && (has_lambda || spatial_re);
  Eigen::SimplicialLLT<SpMat> re_chol;
  SpMat I_N(N, N);
  I_N.setIdentity();
  SpMat MpMt(N, N), MtM(N, N), WpWt(N, N), WtW(N, N);
  if (has_lambda && has_re) {
    for (std::size_t i = 0; i < M.size(); ++i) {
      const SpMat Mt = M[i].transpose();
      MpMt += M[i] + Mt;
      MtM += Mt * M[i];
    }
    if (M.size() == 1) {
      MpMt *= double(T);
      MtM *= double(T);
    }
  }
  if (spatial_re) {
    const SpMat Wt = W[0].transpose();
    WpWt = W[0] + Wt;
    WtW = Wt * W[0];
  }
  // structural pattern of P (values irrelevant here)
  if (sparse_re) re_chol.analyzePattern(SpMat(I_N + MpMt + MtM + WpWt + WtW));

  // initial latent values consistent with the observed outcomes
  VectorXd z(n);
  for (int i = 0; i < n; ++i) {
    if (!ordered) {
      z[i] = (y[i] == 1) ? 0.5 : ((y[i] == 0) ? -0.5 : 0.0);
    } else {
      const int c = y[i];
      if (c == 0) z[i] = 0.0;
      else if (c == 1) z[i] = -0.5;
      else if (c == J) z[i] = kappa[J - 1] + 0.5;
      else z[i] = 0.5 * (kappa[c - 1] + kappa[c]);
    }
  }

  // latent AR: the lag column holds z_{t-1} (zero in the first period)
  auto refresh_lag = [&]() {
    if (!latent_ar) return;
    Xm.col(lag_col).head(N).setZero();
    Xm.col(lag_col).tail(n - N) = z.head(n - N);
  };
  refresh_lag();

  VectorXd Wz = VectorXd::Zero(n), r(n), a0(n), a1(n), e(n), work(n);

  const int n_kappa = ordered ? J - 2 : 0;
  const int p = k + (has_rho ? 1 : 0) + (has_lambda ? 1 : 0) + (has_re ? 1 : 0) +
                (spatial_re ? 1 : 0) + n_kappa;
  const int nsave = iter / thin;
  Rcpp::NumericMatrix draws(nsave, p);
  VectorXd z_sum = VectorXd::Zero(n), a_sum = VectorXd::Zero(N), a_sq = VectorXd::Zero(N);
  // alpha draws kept every 'alpha_every' saved draws (for unit-specific impacts)
  const int n_alpha = (has_re && alpha_every > 0) ? (nsave + alpha_every - 1) / alpha_every : 0;
  Rcpp::NumericMatrix alpha_draws(n_alpha, has_re ? N : 0);

  // transformation-move tuning (one step size per free cut-point)
  std::vector<double> step(J + 1, 0.5);
  std::vector<int> acc(J + 1, 0), att(J + 1, 0), acc_w(J + 1, 0), att_w(J + 1, 0);

  int save = 0;
  const int total = warmup + iter;
  for (int it = 0; it < total; ++it) {
    if (it % 100 == 0) Rcpp::checkUserInterrupt();
    if (progress && total >= 10 && (it + 1) % (total / 10) == 0) {
      Rprintf("Chain %d: %3d%% (%s)\n", chain, 10 * ((it + 1) / (total / 10)),
              it < warmup ? "warm-up" : "sampling");
    }

    // 1. latent block
    refresh_lag();
    compute_resid(m, W, M, Xm, z, beta, alpha, rho, lambda, work, e, r);
    for (int s = 0; s < z_sweeps; ++s) {
      latent_sweep(m, A, colsq, y, kappa, wtp, gamma_now(), Rs, z, r);
    }
    if (latent_ar) {
      refresh_lag();
      refresh_design();
    }

    // 2. rho with beta integrated out (given lambda), then 3. beta.
    //    a0 = R (z - Delta alpha), a1 = R W z; the regression is
    //    a0 - rho a1 = Xt beta + eps.
    a0 = z;
    if (has_re) for (int t = 0; t < T; ++t) a0.segment(t * N, N) -= alpha;
    if (has_lambda) {
      apply_op(M, a0, work, N, T);
      a0 -= lambda * work;
    }
    if (has_rho) {
      apply_op(W, z, a1, N, T);
      if (has_lambda) {
        apply_op(M, a1, work, N, T);
        a1 -= lambda * work;
      }
      const VectorXd c0 = a0 - Xt * beta_mean;
      const VectorXd wc0 = has_t ? VectorXd(wt.cwiseProduct(c0)) : c0;
      const VectorXd wa1 = has_t ? VectorXd(wt.cwiseProduct(a1)) : a1;
      const VectorXd u0 = Xt.transpose() * wc0;
      const VectorXd u1 = Xt.transpose() * wa1;
      const double q00 = c0.dot(wc0) - u0.dot(B * u0);
      const double q01 = wa1.dot(c0) - u1.dot(B * u0);
      const double q11 = a1.dot(wa1) - u1.dot(B * u1);
      rho = draw_rho(rho_grid, lndet_grid, T, q00, q01, q11, rho_a, rho_b, rmin, rmax);
      a0 -= rho * a1;
    }
    {
      const VectorXd mean =
          B * (Xt.transpose() * (has_t ? VectorXd(wt.cwiseProduct(a0)) : a0) + Pc);
      VectorXd zn(k);
      for (int j = 0; j < k; ++j) zn[j] = norm_rand();
      beta = mean + LB * zn;
    }

    // 4. lambda | z, rho, beta, alpha (plan/02 s2.5): with e = S z - X beta -
    //    Delta alpha, the kernel is |e - lambda M e|^2.
    if (has_lambda) {
      compute_resid(m, W, M, Xm, z, beta, alpha, rho, lambda, work, e, r);
      apply_op(M, e, work, N, T);
      const VectorXd we = has_t ? VectorXd(wt.cwiseProduct(e)) : e;
      const VectorXd wMe = has_t ? VectorXd(wt.cwiseProduct(work)) : work;
      const double q00 = e.dot(we), q01 = we.dot(work), q11 = work.dot(wMe);
      lambda = draw_rho(lam_grid, lndetM_grid, T, q00, q01, q11, lam_a, lam_b, lmin, lmax);
      refresh_design();
      build_R();
    }
    build_A(m, W, M, rho, lambda, A, colsq);

    // 5. random effects and their variance
    if (has_re) {
      // u_t = S_t z_t - X_t beta, summed over periods
      compute_resid(m, W, M, Xm, z, beta, alpha, rho, lambda, work, e, r);
      VectorXd usum = VectorXd::Zero(N);
      for (int t = 0; t < T; ++t) usum += e.segment(t * N, N) + alpha;
      if (!sparse_re) {
        // diagonal precision sum_t w_it + 1 / sigma2
        VectorXd prec = VectorXd::Constant(N, 1.0 / sigma2), wu = VectorXd::Zero(N);
        for (int t = 0; t < T; ++t) {
          const auto wt_t = wt.segment(t * N, N);
          prec += wt_t;
          wu += wt_t.cwiseProduct(e.segment(t * N, N) + alpha);
        }
        for (int i = 0; i < N; ++i) {
          alpha[i] = wu[i] / prec[i] + norm_rand() / std::sqrt(prec[i]);
        }
      } else {
        // precision P = sum_t R_t'D_t R_t + B'B / sigma2 (D_t = I for normal
        // errors), mean P^{-1} sum_t R_t'D_t R_t u_t
        SpMat P;
        VectorXd rhs = VectorXd::Zero(N);
        if (!has_t) {
          P = (T + 1.0 / sigma2) * I_N - lambda * MpMt + (lambda * lambda) * MtM -
              (phi / sigma2) * WpWt + (phi * phi / sigma2) * WtW;
        } else {
          // same structural pattern as analysed (the zero-weighted terms keep it)
          P = (1.0 / sigma2) * I_N + 0.0 * (MpMt + MtM) - (phi / sigma2) * WpWt +
              (phi * phi / sigma2) * WtW;
          for (int t = 0; t < T; ++t) {
            const SpMat Dt = SpMat(wt.segment(t * N, N).asDiagonal());
            const VectorXd u = e.segment(t * N, N) + alpha;
            if (has_lambda) {
              const SpMat Rt = I_N - lambda * at(M, t);
              const SpMat RtD = SpMat(Rt.transpose()) * Dt;
              P += RtD * Rt;
              rhs += RtD * (Rt * u);
            } else {
              P += Dt;
              rhs += wt.segment(t * N, N).cwiseProduct(u);
            }
          }
        }
        re_chol.factorize(P);
        if (re_chol.info() != Eigen::Success) Rcpp::stop("random-effects Cholesky failed");
        if (has_t) {
          // rhs computed above
        } else if (!has_lambda) {
          rhs = usum;
        } else if (M.size() == 1) {
          const VectorXd w = usum - lambda * (M[0] * usum);
          rhs = w - lambda * (M[0].transpose() * w);
        } else {
          for (int t = 0; t < T; ++t) {
            const VectorXd u = e.segment(t * N, N) + alpha;
            const VectorXd w = u - lambda * (M[t] * u);
            rhs += w - lambda * (M[t].transpose() * w);
          }
        }
        const VectorXd mean = re_chol.solve(rhs);
        VectorXd zn(N);
        for (int i = 0; i < N; ++i) zn[i] = norm_rand();
        const VectorXd u = re_chol.matrixU().solve(zn);
        alpha = mean + re_chol.permutationPinv() * u;
      }
      // phi | alpha, sigma2: kernel |alpha - phi W alpha|^2 / sigma2 with log|B|
      if (spatial_re) {
        const VectorXd Wa = W[0] * alpha;
        phi = draw_rho(rho_grid, lndet_grid, 1, alpha.squaredNorm() / sigma2,
                       alpha.dot(Wa) / sigma2, Wa.squaredNorm() / sigma2, phi_a, phi_b, rmin,
                       rmax);
      }
      const double ss = spatial_re ? (alpha - phi * (W[0] * alpha)).squaredNorm()
                                   : alpha.squaredNorm();
      if (sigma_prior == 0) {  // half-Cauchy(0, A) via inverse-gamma mixture
        sigma2 = rinvgamma(0.5 * (N + 1.0), 0.5 * ss + 1.0 / xi);
        xi = rinvgamma(1.0, 1.0 / sigma2 + 1.0 / (sigma_p1 * sigma_p1));
      } else {  // inverse gamma(shape, scale) on sigma^2
        sigma2 = rinvgamma(sigma_p1 + 0.5 * N, sigma_p2 + 0.5 * ss);
      }
    }

    // 5b. mixing variances of t errors: v_i | . ~ IG((nu + 1)/2, (nu + r_i^2)/2)
    if (has_t) {
      compute_resid(m, W, M, Xm, z, beta, alpha, rho, lambda, work, e, r);
      for (int i = 0; i < n; ++i) wt[i] = R::rgamma(0.5 * (nu + 1.0), 2.0 / (nu + r[i] * r[i]));
      refresh_design();
    }

    // 6. cut-points
    if (ordered) {
      compute_resid(m, W, M, Xm, z, beta, alpha, rho, lambda, work, e, r);
      cut_ac_step(m, y, z, kappa, kappa_flat, kappa_mean, kappa_sd);
      if (cut_sampler == 0) {
        for (int j = 2; j <= J - 1; ++j) {
          const bool ok = cut_transform_move(m, A, y, ncat, j, step[j], kappa_flat,
                                             kappa_mean, kappa_sd, wtp, gamma_now(), Rs, kappa,
                                             z, r, work);
          att[j]++;
          att_w[j]++;
          if (ok) {
            acc[j]++;
            acc_w[j]++;
          }
          // adapt during warm-up only
          if (it < warmup && att_w[j] == 50) {
            const double rate = acc_w[j] / 50.0;
            step[j] *= (rate > 0.35) ? 1.25 : 0.8;
            acc_w[j] = 0;
            att_w[j] = 0;
          }
        }
      }
    }

    // 7. scale move (generalised Gibbs, Liu & Sabatti 2000): x -> g x for
    //    x = (z, beta, alpha, kappa). The outcome is invariant; the target
    //    for g is p(g) ~ g^(d-1) exp(-C g^2 / 2 + L g) (Haar measure dg/g).
    if (scale_move) {
      refresh_lag();  // z may have moved with the cut-points
      compute_resid(m, W, M, Xm, z, beta, alpha, rho, lambda, work, e, r);
      // gamma (latent AR) is not rescaled: its column z_{t-1} already is
      VectorXd bs = beta;
      if (latent_ar) bs[lag_col] = 0.0;
      double C = (has_t ? r.dot(wt.cwiseProduct(r)) : r.squaredNorm()) + bs.dot(beta_prec * bs);
      double L = bs.dot(Pc);
      double d = n + k - (latent_ar ? 1 : 0);
      if (has_re) {
        C += (spatial_re ? (alpha - phi * (W[0] * alpha)).squaredNorm() : alpha.squaredNorm()) /
             sigma2;
        d += N;
      }
      if (ordered) {
        for (int j = 2; j <= J - 1; ++j) {
          if (!kappa_flat) {
            C += kappa[j] * kappa[j] / (kappa_sd * kappa_sd);
            L += kappa[j] * kappa_mean / (kappa_sd * kappa_sd);
          }
          d += 1.0;
        }
      }
      double g = 1.0;
      if (L == 0.0) {
        g = std::sqrt(R::rgamma(0.5 * d, 2.0 / C));  // g^2 ~ Gamma(d/2, rate C/2)
      } else {
        // random walk on log g, symmetric; acceptance includes the Jacobian g^d
        const double prop = std::exp(norm_rand() / std::sqrt(2.0 * d));
        const double logr = -0.5 * (prop * prop - 1.0) * C + (prop - 1.0) * L +
                            d * std::log(prop);
        if (std::log(unif_rand()) < logr) g = prop;
      }
      z *= g;
      const double gam = gamma_now();
      beta *= g;
      if (latent_ar) beta[lag_col] = gam;
      alpha *= g;
      for (int j = 2; j <= J - 1 && ordered; ++j) kappa[j] *= g;
    }

    // store
    if (it >= warmup && (it - warmup) % thin == 0 && save < nsave) {
      int c = 0;
      for (int j = 0; j < k; ++j) draws(save, c++) = beta[j];
      if (has_rho) draws(save, c++) = rho;
      if (has_lambda) draws(save, c++) = lambda;
      if (has_re) draws(save, c++) = std::sqrt(sigma2);
      if (spatial_re) draws(save, c++) = phi;
      for (int j = 2; j <= J - 1 && ordered; ++j) draws(save, c++) = kappa[j];
      z_sum += z;
      if (has_t) v_sum += wt.cwiseInverse();
      a_sum += alpha;
      a_sq += alpha.cwiseProduct(alpha);
      if (n_alpha > 0 && save % alpha_every == 0) {
        for (int i = 0; i < N; ++i) alpha_draws(save / alpha_every, i) = alpha[i];
      }
      ++save;
    }
    if (it == warmup - 1) {
      std::fill(acc.begin(), acc.end(), 0);
      std::fill(att.begin(), att.end(), 0);
    }
  }

  Rcpp::NumericVector acc_rate(n_kappa), steps(n_kappa);
  for (int j = 2; j <= J - 1 && ordered; ++j) {
    acc_rate[j - 2] = att[j] > 0 ? static_cast<double>(acc[j]) / att[j] : NA_REAL;
    steps[j - 2] = step[j];
  }
  const VectorXd a_mean = a_sum / save;
  const VectorXd a_var = (a_sq / save - a_mean.cwiseProduct(a_mean)).cwiseMax(0.0);
  return Rcpp::List::create(
      Rcpp::Named("draws") = draws, Rcpp::Named("z_mean") = z_sum / save,
      Rcpp::Named("alpha_mean") = a_mean, Rcpp::Named("alpha_sd") = a_var.cwiseSqrt(),
      Rcpp::Named("cut_accept") = acc_rate, Rcpp::Named("cut_step") = steps,
      Rcpp::Named("alpha_draws") = alpha_draws, Rcpp::Named("v_mean") = v_sum / save);
}

// Expose the truncated normal sampler for testing.
// [[Rcpp::export]]
Rcpp::NumericVector psp_rtnorm_cpp(int n, double mean, double sd, double lo, double hi) {
  Rcpp::NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = psp::rtnorm(mean, sd, lo, hi);
  return out;
}

// Toy harness for the cut-point updates (no spatial term, unit variance):
// latent z_i ~ N(mu_i, 1) with categories y and free cut-points kappa_2.. (as
// in the sampler, prior N(kmean, ksd) restricted to the ordering). Applies
// 'reps' rounds of the Albert-Chib step (ac = TRUE) and/or the transformation
// move (move = TRUE) and returns the final kappa (free ones) and z.
// [[Rcpp::export]]
Rcpp::List psp_cut_kernel_cpp(const Eigen::VectorXi& y, const Eigen::VectorXd& mu,
                              Eigen::VectorXd z, const Eigen::VectorXd& kappa_free, int J,
                              double kmean, double ksd, double step, int reps, bool ac,
                              bool move) {
  Model m;
  m.N = y.size();
  m.T = 1;
  m.n = m.N;
  m.k = 0;
  m.J = J;
  m.ordered = true;
  m.has_rho = false;
  m.has_lambda = false;
  m.has_re = false;
  m.spatial_re = false;
  std::vector<double> kappa(J + 1);
  kappa[0] = -INF;
  kappa[1] = 0.0;
  kappa[J] = INF;
  for (int j = 2; j < J; ++j) kappa[j] = kappa_free[j - 2];
  std::vector<int> ncat(J + 1, 0);
  for (int i = 0; i < m.n; ++i) ncat[y[i]]++;
  Ops A(1, SpMat(m.N, m.N));
  A[0].setIdentity();
  VectorXd r = z - mu, work(m.n);
  for (int rep = 0; rep < reps; ++rep) {
    if (ac) cut_ac_step(m, y, z, kappa, false, kmean, ksd);
    if (move) {
      for (int j = 2; j < J; ++j) {
        cut_transform_move(m, A, y, ncat, j, step, false, kmean, ksd, nullptr, 0.0, Ops(),
                           kappa, z, r, work);
      }
    }
  }
  Rcpp::NumericVector kf(J - 2);
  for (int j = 2; j < J; ++j) kf[j - 2] = kappa[j];
  return Rcpp::List::create(Rcpp::Named("kappa") = kf, Rcpp::Named("z") = z);
}
