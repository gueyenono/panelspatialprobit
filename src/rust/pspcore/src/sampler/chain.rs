//! One MCMC chain: a port of `psp_sampler_cpp` in the C++ engine
//! (src/sampler.cpp at tag `cpp-final`), block for block and in the same
//! order, consuming random numbers in the same order (plan/06 §6.7).
//!
//! Every model is the kernel |A z - b|^2 with block-diagonal A = R S,
//! S_t = I - rho W_t, R_t = I - lambda M_t (plan/02 §2.2).

use super::griddy::draw_grid_param;
use super::input::{ChainOutput, CutSampler, GridPrior, KappaPrior, SamplerInput, SigmaPrior};
use super::re::ReBlock;
use crate::dense::{cholesky, dot, l_solve, lt_solve, sq_matvec, x_mul, xt_mul};
use crate::dist::rtnorm;
use crate::par::{self, TaskCtx};
use crate::rng::{RMrg32k3a, Rng};
use crate::sparse::{Csc, FilterPattern};
use crate::PspError;

const INF: f64 = f64::INFINITY;

fn rinvgamma<R: Rng>(rng: &mut R, shape: f64, scale: f64) -> f64 {
    scale / rng.gamma(shape, 1.0)
}

fn log_kappa_prior(x: f64, prior: KappaPrior) -> f64 {
    match prior {
        KappaPrior::Flat => 0.0,
        KappaPrior::Normal { mean, sd } => {
            let d = (x - mean) / sd;
            -0.5 * d * d
        }
    }
}

/// The filters A_t, one per distinct pair (W_t, M_t), and the factors R_t =
/// I - lambda M_t of the latent AR coupling, one per distinct M_t.
struct Filters {
    pats: Vec<FilterPattern>,
    mats: Vec<Csc>,
    colsq: Vec<Vec<f64>>,
    pair: Vec<usize>,
    r_pats: Vec<FilterPattern>,
    r_mats: Vec<Csc>,
}

impl Filters {
    fn new(inp: &SamplerInput, latent_ar: bool) -> Result<Self, PspError> {
        let n = inp.n_units;
        let mut keys: Vec<(usize, usize)> = Vec::new();
        let mut pair = Vec::with_capacity(inp.periods);
        for t in 0..inp.periods {
            let key = (
                if inp.has_rho { inp.w.period[t] } else { 0 },
                if inp.has_lambda { inp.m.period[t] } else { 0 },
            );
            let p = match keys.iter().position(|&k| k == key) {
                Some(p) => p,
                None => {
                    keys.push(key);
                    keys.len() - 1
                }
            };
            pair.push(p);
        }
        let mut pats = Vec::with_capacity(keys.len());
        for &(wu, mu) in &keys {
            pats.push(FilterPattern::new(
                n,
                inp.has_rho.then(|| &inp.w.mats[wu]),
                inp.has_lambda.then(|| &inp.m.mats[mu]),
            )?);
        }
        let mats = pats.iter().map(|p| p.new_matrix()).collect();
        let colsq = vec![vec![0.0; n]; pats.len()];
        let (mut r_pats, mut r_mats) = (Vec::new(), Vec::new());
        if latent_ar && inp.has_lambda {
            for m in &inp.m.mats {
                let p = FilterPattern::new(n, None, Some(m))?;
                r_mats.push(p.new_matrix());
                r_pats.push(p);
            }
        }
        Ok(Self {
            pats,
            mats,
            colsq,
            pair,
            r_pats,
            r_mats,
        })
    }

    fn build_a(&mut self, rho: f64, lambda: f64) {
        for ((p, a), cs) in self.pats.iter().zip(&mut self.mats).zip(&mut self.colsq) {
            p.update(rho, lambda, a, cs);
        }
    }

    fn build_r(&mut self, lambda: f64) {
        let mut scratch = Vec::new();
        for (p, r) in self.r_pats.iter().zip(&mut self.r_mats) {
            scratch.resize(p.n, 0.0);
            p.update(0.0, lambda, r, &mut scratch);
        }
    }

    fn a(&self, t: usize) -> &Csc {
        &self.mats[self.pair[t]]
    }

    /// R_t of the latent AR coupling (`None`: the identity).
    fn r<'a>(&'a self, inp: &SamplerInput, t: usize) -> Option<&'a Csc> {
        (!self.r_mats.is_empty()).then(|| &self.r_mats[inp.m.period[t]])
    }
}

/// Truncation bounds of the latent values given the outcome codes: binary
/// 1 / 0 (-1 missing), ordered 1..J (0 missing) with cut-points `kappa`.
#[derive(Clone, Copy)]
struct Bounds<'a> {
    ordered: bool,
    kappa: &'a [f64],
}

impl Bounds<'_> {
    #[inline]
    fn of(&self, y: i32) -> (f64, f64) {
        if self.ordered {
            if y == 0 {
                (-INF, INF)
            } else {
                (self.kappa[y as usize - 1], self.kappa[y as usize])
            }
        } else if y == 1 {
            (0.0, INF)
        } else if y == 0 {
            (-INF, 0.0)
        } else {
            (-INF, INF)
        }
    }
}

/// What one period's sweep reads: A_t with its squared column norms, R_{t+1}
/// of the latent AR coupling (`None`: the identity), the outcomes and the t
/// weights of periods t and t + 1.
struct Period<'a> {
    a: &'a Csc,
    colsq: &'a [f64],
    r_next: Option<&'a Csc>,
    y: &'a [i32],
    wt: &'a [f64],
    wt1: &'a [f64],
}

impl Period<'_> {
    /// One sweep over z_t, updating r_t and, with the latent AR coupling
    /// (`LAG`), r_{t+1}. `T_ERR`: weights D_t from t errors.
    #[allow(clippy::too_many_arguments)]
    fn sweep<const T_ERR: bool, const LAG: bool>(
        &self,
        rng: &mut RMrg32k3a,
        bounds: Bounds,
        gamma: f64,
        z: &mut [f64],
        r: &mut [f64],
        r1: &mut [f64],
    ) {
        let a = self.a;
        for k in 0..z.len() {
            let (cols, vals) = {
                let (p0, p1) = (a.colptr[k], a.colptr[k + 1]);
                (&a.rowind[p0..p1], &a.val[p0..p1])
            };
            let (mut ar, mut h) = (0.0, 0.0);
            if T_ERR {
                for (&row, &v) in cols.iter().zip(vals) {
                    let wa = self.wt[row] * v;
                    ar += wa * r[row];
                    h += wa * v;
                }
            } else {
                for (&row, &v) in cols.iter().zip(vals) {
                    ar += v * r[row];
                }
                h = self.colsq[k];
            }
            if LAG {
                match self.r_next {
                    None => {
                        let wv = if T_ERR { self.wt1[k] } else { 1.0 };
                        ar -= gamma * wv * r1[k];
                        h += gamma * gamma * wv;
                    }
                    Some(rm) => {
                        for (row, v) in rm.col(k) {
                            let wv = if T_ERR { self.wt1[row] } else { 1.0 };
                            ar -= gamma * wv * v * r1[row];
                            h += gamma * gamma * wv * v * v;
                        }
                    }
                }
            }
            let mean = z[k] - ar / h;
            let (lo, hi) = bounds.of(self.y[k]);
            let znew = rtnorm(rng, mean, 1.0 / h.sqrt(), lo, hi);
            let delta = znew - z[k];
            z[k] = znew;
            for (&row, &v) in cols.iter().zip(vals) {
                r[row] += v * delta;
            }
            if LAG {
                match self.r_next {
                    None => r1[k] -= gamma * delta,
                    Some(rm) => {
                        for (row, v) in rm.col(k) {
                            r1[row] -= gamma * v * delta;
                        }
                    }
                }
            }
        }
    }
}

struct Chain<'a> {
    inp: &'a SamplerInput,
    n: usize,
    nn: usize,
    k: usize,
    latent_ar: bool,
    has_t: bool,
    // state
    beta: Vec<f64>,
    rho: f64,
    lambda: f64,
    sigma2: f64,
    phi: f64,
    xi: f64,
    alpha: Vec<f64>,
    kappa: Vec<f64>,
    // observations by outcome category (ordered models)
    cat_idx: Vec<Vec<usize>>,
    z: Vec<f64>,
    wt: Vec<f64>,
    // design: Xm (X with the latent AR lag column), Xt = R Xm, the
    // precision P = Xt' D Xt + P0 of beta and its Cholesky factor
    xm: Vec<f64>,
    xt: Vec<f64>,
    prec: Vec<f64>,
    lp: Vec<f64>,
    pc: Vec<f64>,
    mx: Vec<f64>,
    gram: Option<[Vec<f64>; 3]>,
    filters: Filters,
    // work: e and r as left by `compute_resid`; r is kept current by the
    // blocks that change z or alpha in place, and `resid_ok` says whether it
    // is (e is current only right after `compute_resid`)
    e: Vec<f64>,
    r: Vec<f64>,
    resid_ok: bool,
    work: Vec<f64>,
}

impl<'a> Chain<'a> {
    /// |alpha|^2, or |alpha - phi W alpha|^2 for spatial random effects.
    fn re_ss(&self) -> f64 {
        if self.inp.spatial_re {
            let mut wa = vec![0.0; self.nn];
            self.inp.w.mul_one(0, &self.alpha, &mut wa);
            self.alpha
                .iter()
                .zip(&wa)
                .map(|(a, w)| (a - self.phi * w).powi(2))
                .sum()
        } else {
            dot(&self.alpha, &self.alpha)
        }
    }

    fn gamma_now(&self) -> f64 {
        self.inp.lag_col.map_or(0.0, |c| self.beta[c])
    }

    /// The latent AR lag column holds z_{t-1} (zero in the first period).
    fn refresh_lag(&mut self) {
        if let Some(c) = self.inp.lag_col {
            let (n, nn) = (self.n, self.nn);
            let col = &mut self.xm[c * n..(c + 1) * n];
            col[..nn].fill(0.0);
            col[nn..].copy_from_slice(&self.z[..n - nn]);
        }
    }

    /// Column j of M Xm, period by period.
    fn compute_mx_col(&mut self, j: usize) {
        let (n, nn) = (self.n, self.nn);
        for t in 0..self.inp.periods {
            let r = j * n + t * nn..j * n + (t + 1) * nn;
            self.inp
                .m
                .at(t)
                .mul_vec(&self.xm[r.clone()], &mut self.mx[r]);
        }
    }

    /// M Xm (computed once: only the latent AR lag column changes) and, with
    /// normal errors and a fixed design, the Gram matrices X'X, X'MX and
    /// (MX)'MX, so that Xt'Xt costs O(k^2) after each lambda draw.
    fn init_design(&mut self) {
        let (n, k) = (self.n, self.k);
        if !self.inp.has_lambda {
            return;
        }
        self.mx = vec![0.0; n * k];
        for j in 0..k {
            self.compute_mx_col(j);
        }
        if !self.has_t && !self.latent_ar {
            let col = |a: &[f64], j: usize| a[j * n..(j + 1) * n].to_vec();
            let mut g = [vec![0.0; k * k], vec![0.0; k * k], vec![0.0; k * k]];
            for j in 0..k {
                for i in 0..k {
                    let (xi, mxi) = (col(&self.xm, i), col(&self.mx, i));
                    let (xj, mxj) = (col(&self.xm, j), col(&self.mx, j));
                    g[0][j * k + i] = dot(&xi, &xj);
                    g[1][j * k + i] = dot(&xi, &mxj);
                    g[2][j * k + i] = dot(&mxi, &mxj);
                }
            }
            self.gram = Some(g);
        }
    }

    /// Xt = R Xm and the precision P = Xt' D Xt + P0 with its Cholesky
    /// factor.
    fn refresh_design(&mut self) -> Result<(), PspError> {
        let (n, k, inp) = (self.n, self.k, self.inp);
        if inp.has_lambda {
            if let Some(c) = inp.lag_col {
                self.compute_mx_col(c);
            }
            for x in 0..n * k {
                self.xt[x] = self.xm[x] - self.lambda * self.mx[x];
            }
        } else if self.latent_ar {
            self.xt.copy_from_slice(&self.xm);
        }
        let mut xtx = inp.beta_prec.clone();
        if let Some(g) = &self.gram {
            // (X - l MX)'(X - l MX) = X'X - l (X'MX + (X'MX)') + l^2 (MX)'MX
            let (l, l2) = (self.lambda, self.lambda * self.lambda);
            for j in 0..k {
                for i in 0..k {
                    xtx[j * k + i] += g[0][j * k + i] - l * (g[1][j * k + i] + g[1][i * k + j])
                        + l2 * g[2][j * k + i];
                }
            }
        } else {
            let mut wcol = vec![0.0; n];
            for j in 0..k {
                let cj = &self.xt[j * n..(j + 1) * n];
                if self.has_t {
                    for i in 0..n {
                        wcol[i] = self.wt[i] * cj[i];
                    }
                } else {
                    wcol.copy_from_slice(cj);
                }
                for i in 0..=j {
                    let s = dot(&self.xt[i * n..(i + 1) * n], &wcol);
                    xtx[j * k + i] += s;
                    if i != j {
                        xtx[i * k + j] += s;
                    }
                }
            }
        }
        self.lp = cholesky(&xtx, k)?;
        self.prec = xtx;
        Ok(())
    }

    /// e = S z - Xm beta - Delta alpha and r = R e = A z - b; with an error
    /// term, leaves M e in `work` (the lambda block uses it).
    fn compute_resid(&mut self) {
        let (n, nn, inp) = (self.n, self.nn, self.inp);
        x_mul(&self.xm, n, &self.beta, &mut self.work);
        for i in 0..n {
            self.e[i] = self.z[i] - self.work[i];
        }
        if inp.has_rho {
            inp.w.apply(&self.z, &mut self.work);
            for i in 0..n {
                self.e[i] -= self.rho * self.work[i];
            }
        }
        if inp.has_re {
            for t in 0..inp.periods {
                for i in 0..nn {
                    self.e[t * nn + i] -= self.alpha[i];
                }
            }
        }
        self.r.copy_from_slice(&self.e);
        if inp.has_lambda {
            inp.m.apply(&self.e, &mut self.work);
            for i in 0..n {
                self.r[i] -= self.lambda * self.work[i];
            }
        }
        self.resid_ok = true;
    }

    /// r = A z - b, recomputed only if a block has changed parameters since.
    fn ensure_resid(&mut self) {
        if !self.resid_ok {
            self.compute_resid();
        } else if cfg!(debug_assertions) {
            // the r kept in place must match a fresh computation
            let kept = self.r.clone();
            self.compute_resid();
            for (a, b) in kept.iter().zip(&self.r) {
                assert!((a - b).abs() <= 1e-9 * (1.0 + b.abs()), "stale residual");
            }
            self.r = kept;
        }
    }

    /// r -= R_t d for every period: alpha has moved by d.
    fn shift_resid_alpha(&mut self, d: &[f64]) {
        let (nn, inp) = (self.nn, self.inp);
        let mut rd = vec![0.0; nn];
        let mut last = None;
        for t in 0..inp.periods {
            if inp.has_lambda {
                let u = inp.m.period[t];
                if last != Some(u) {
                    inp.m.mul_one(u, d, &mut rd);
                    for i in 0..nn {
                        rd[i] = d[i] - self.lambda * rd[i];
                    }
                    last = Some(u);
                }
            } else if last.is_none() {
                rd.copy_from_slice(d);
                last = Some(0);
            }
            for i in 0..nn {
                self.r[t * nn + i] -= rd[i];
            }
        }
    }

    /// The latent block: `z_sweeps` Gibbs sweeps over z with residual
    /// updates (plan/02 §2.3), z_k | . ~ TN(z_k - a_k'D r / h_k, 1 / h_k),
    /// h_k = a_k'D a_k, then r += a_k (z_k' - z_k).
    ///
    /// Each period has its own random-number substream (`prngs`), so the
    /// draws do not depend on `inp.threads` (plan/06 R6). Without latent AR
    /// coupling the periods are conditionally independent: each runs all its
    /// sweeps as one task. With it, z_{k,t} also enters period t + 1 through
    /// column k of -gamma R_{t+1}, so each sweep visits the even periods and
    /// then the odd ones (chromatic Gibbs; a period writes only r_t and
    /// r_{t+1}).
    fn latent_block(&mut self, prngs: &mut [RMrg32k3a]) {
        let (nn, inp) = (self.nn, self.inp);
        let tt = inp.periods;
        let gamma = self.gamma_now();
        let coupled = gamma != 0.0 && tt > 1;
        let bounds = Bounds {
            ordered: inp.ordered,
            kappa: &self.kappa,
        };
        let filters = &self.filters;
        let (wt, has_t) = (&self.wt[..], self.has_t);
        let task = |t: usize, z: &mut [f64], r: &mut [f64], r1: &mut [f64], rng: &mut RMrg32k3a| {
            let p = Period {
                a: filters.a(t),
                colsq: &filters.colsq[filters.pair[t]],
                r_next: if r1.is_empty() {
                    None
                } else {
                    filters.r(inp, t + 1)
                },
                y: &inp.y[t * nn..(t + 1) * nn],
                wt: &wt[t * nn..(t + 1) * nn],
                wt1: if r1.is_empty() {
                    &[]
                } else {
                    &wt[(t + 1) * nn..(t + 2) * nn]
                },
            };
            let sweeps = if coupled { 1 } else { inp.z_sweeps };
            for _ in 0..sweeps {
                match (has_t, !r1.is_empty()) {
                    (false, false) => p.sweep::<false, false>(rng, bounds, gamma, z, r, r1),
                    (true, false) => p.sweep::<true, false>(rng, bounds, gamma, z, r, r1),
                    (false, true) => p.sweep::<false, true>(rng, bounds, gamma, z, r, r1),
                    (true, true) => p.sweep::<true, true>(rng, bounds, gamma, z, r, r1),
                }
            }
        };
        if !coupled {
            let tasks: Vec<_> = self
                .z
                .chunks_mut(nn)
                .zip(self.r.chunks_mut(nn))
                .zip(prngs.iter_mut())
                .enumerate()
                .collect();
            par::for_each(tasks, inp.threads, |(t, ((z, r), rng))| {
                task(t, z, r, &mut [], rng)
            });
            return;
        }
        for _ in 0..inp.z_sweeps {
            for parity in 0..2 {
                // period t = parity + 2 i writes r_t and r_{t+1}
                let tasks: Vec<_> = self.z[parity * nn..]
                    .chunks_mut(2 * nn)
                    .zip(self.r[parity * nn..].chunks_mut(2 * nn))
                    .zip(prngs[parity..].iter_mut().step_by(2))
                    .enumerate()
                    .map(|(i, ((z, r), rng))| {
                        let (r, r1) = r.split_at_mut(nn);
                        (parity + 2 * i, &mut z[..nn], r, r1, rng)
                    })
                    .collect();
                par::for_each(tasks, inp.threads, |(t, z, r, r1, rng)| {
                    task(t, z, r, r1, rng)
                });
            }
        }
    }

    /// Albert-Chib Gibbs step for the free cut-points given z.
    fn cut_ac_step<R: Rng>(&mut self, rng: &mut R) {
        let inp = self.inp;
        for j in 2..inp.n_cat {
            let (mut lo, mut hi) = (self.kappa[j - 1], self.kappa[j + 1]);
            for &i in &self.cat_idx[j] {
                lo = lo.max(self.z[i]);
            }
            for &i in &self.cat_idx[j + 1] {
                hi = hi.min(self.z[i]);
            }
            if hi > lo {
                self.kappa[j] = match inp.kappa_prior {
                    KappaPrior::Flat => lo + (hi - lo) * rng.unif(),
                    KappaPrior::Normal { mean, sd } => rtnorm(rng, mean, sd, lo, hi),
                };
            }
        }
    }

    /// Joint (kappa_j, z) transformation move (plan/02 §2.7). Updates kappa,
    /// z and r in place on acceptance; returns acceptance.
    fn cut_transform_move<R: Rng>(
        &mut self,
        rng: &mut R,
        ncat: &[usize],
        j: usize,
        step: f64,
        delta: &mut [f64],
        ad: &mut [f64],
    ) -> bool {
        let (n, nn, inp) = (self.n, self.nn, self.inp);
        let (l, u, kj) = (self.kappa[j - 1], self.kappa[j + 1], self.kappa[j]);
        let top = j == inp.n_cat - 1;
        let eps = step * rng.norm();
        let (knew, s1, s2, logjac);
        if !top {
            let pp = (kj - l) / (u - l);
            let lg = (pp / (1.0 - pp)).ln() + eps;
            let pn = 1.0 / (1.0 + (-lg).exp());
            if pn <= 0.0 || pn >= 1.0 {
                return false;
            }
            knew = l + pn * (u - l);
            s1 = (knew - l) / (kj - l);
            s2 = (u - knew) / (u - kj);
            logjac =
                ncat[j] as f64 * s1.ln() + ncat[j + 1] as f64 * s2.ln() + (pn * (1.0 - pn)).ln()
                    - (pp * (1.0 - pp)).ln();
        } else {
            s1 = eps.exp();
            s2 = 1.0;
            knew = l + (kj - l) * s1;
            logjac = (ncat[j] as f64 + 1.0) * eps;
        }
        // only categories j and j + 1 move: delta and A delta come from
        // their columns alone (ad is dense, zero elsewhere)
        ad.fill(0.0);
        let gamma = self.gamma_now();
        let tt = inp.periods;
        for (c, &i) in self.cat_idx[j]
            .iter()
            .map(|i| (j, i))
            .chain(self.cat_idx[j + 1].iter().map(|i| (j + 1, i)))
        {
            let zi = self.z[i];
            let d = if c == j {
                (l + (zi - l) * s1) - zi
            } else if top {
                knew - kj
            } else {
                (u - (u - zi) * s2) - zi
            };
            delta[i] = d;
            let (t, k) = (i / nn, i % nn);
            let off = t * nn;
            for (row, v) in self.filters.a(t).col(k) {
                ad[off + row] += v * d;
            }
            // latent AR: -gamma R_{t+1} delta_t in period t + 1
            if gamma != 0.0 && t + 1 < tt {
                let off1 = off + nn;
                match self.filters.r(inp, t + 1) {
                    None => ad[off1 + k] -= gamma * d,
                    Some(rm) => {
                        for (row, v) in rm.col(k) {
                            ad[off1 + row] -= gamma * v * d;
                        }
                    }
                }
            }
        }
        let dq = if self.has_t {
            (0..n)
                .map(|i| self.wt[i] * (2.0 * self.r[i] + ad[i]) * ad[i])
                .sum::<f64>()
        } else {
            2.0 * dot(&self.r, ad) + dot(ad, ad)
        };
        let logr = -0.5 * dq + logjac + log_kappa_prior(knew, inp.kappa_prior)
            - log_kappa_prior(kj, inp.kappa_prior);
        if rng.unif().ln() < logr {
            for &i in self.cat_idx[j].iter().chain(&self.cat_idx[j + 1]) {
                self.z[i] += delta[i];
            }
            for (r, a) in self.r.iter_mut().zip(ad.iter()) {
                *r += a;
            }
            self.kappa[j] = knew;
            true
        } else {
            false
        }
    }
}

/// The observations in each category 1..=`n_cat` (index 0 is empty).
fn category_index(y: &[i32], n_cat: usize) -> Vec<Vec<usize>> {
    let mut idx = vec![Vec::new(); n_cat + 1];
    for (i, &c) in y.iter().enumerate() {
        if c >= 1 && (c as usize) <= n_cat {
            idx[c as usize].push(i);
        }
    }
    idx
}

/// Runs one chain. `ctx` gives the cancel flag (checked every 100 iterations)
/// and the progress counter (iterations done).
pub fn run_chain(
    inp: &SamplerInput,
    rng: &mut RMrg32k3a,
    ctx: &TaskCtx,
) -> Result<ChainOutput, PspError> {
    inp.validate()?;
    // the latent sweep of period t draws from substream t + 1 of the chain's
    // stream; everything else from the stream itself
    let mut prngs = Vec::with_capacity(inp.periods);
    let mut sub = rng.next_substream();
    for _ in 0..inp.periods {
        let next = sub.next_substream();
        prngs.push(sub);
        sub = next;
    }
    let mut re_block = if inp.sparse_re() {
        Some(ReBlock::new(inp)?)
    } else {
        None
    };
    // phi uses the rho grid (log|I - phi W|) with its own Beta prior
    let phi_grid = GridPrior {
        a: inp.phi_a,
        b: inp.phi_b,
        ..inp.rho.clone()
    };
    let (n, nn, k, tt) = (inp.n(), inp.n_units, inp.k, inp.periods);
    let j_cat = inp.n_cat;
    let latent_ar = inp.lag_col.is_some();
    let has_t = inp.nu > 0.0;

    let mut kappa = vec![0.0; j_cat + 1];
    kappa[0] = -INF;
    kappa[j_cat] = INF;
    if inp.ordered {
        kappa[1] = 0.0;
        kappa[2..j_cat].copy_from_slice(&inp.kappa_start);
    }
    let mut ncat = vec![0usize; j_cat + 1];
    if inp.ordered {
        for &y in &inp.y {
            if y > 0 {
                ncat[y as usize] += 1;
            }
        }
    }

    let mut pc = vec![0.0; k];
    sq_matvec(&inp.beta_prec, &inp.beta_mean, &mut pc);
    let mut ch = Chain {
        inp,
        n,
        nn,
        k,
        latent_ar,
        has_t,
        beta: inp.beta_start.clone(),
        rho: if inp.has_rho { inp.rho_start } else { 0.0 },
        lambda: if inp.has_lambda {
            inp.lambda_start
        } else {
            0.0
        },
        sigma2: inp.sigma2_start,
        phi: 0.0,
        xi: 1.0,
        alpha: vec![0.0; nn],
        kappa,
        cat_idx: category_index(&inp.y, if inp.ordered { j_cat } else { 0 }),
        z: vec![0.0; n],
        wt: vec![1.0; n],
        xm: inp.x.clone(),
        xt: inp.x.clone(),
        prec: Vec::new(),
        lp: Vec::new(),
        pc,
        mx: Vec::new(),
        gram: None,
        filters: Filters::new(inp, latent_ar)?,
        e: vec![0.0; n],
        r: vec![0.0; n],
        resid_ok: false,
        work: vec![0.0; n],
    };
    ch.filters.build_a(ch.rho, ch.lambda);
    ch.filters.build_r(ch.lambda);
    ch.init_design();
    ch.refresh_design()?;

    // initial latent values consistent with the observed outcomes
    for i in 0..n {
        let c = inp.y[i];
        ch.z[i] = if !inp.ordered {
            if c == 1 {
                0.5
            } else if c == 0 {
                -0.5
            } else {
                0.0
            }
        } else if c == 0 {
            0.0
        } else if c == 1 {
            -0.5
        } else if c as usize == j_cat {
            ch.kappa[j_cat - 1] + 0.5
        } else {
            0.5 * (ch.kappa[c as usize - 1] + ch.kappa[c as usize])
        };
    }
    ch.refresh_lag();

    let mut a0 = vec![0.0; n];
    let mut a1 = vec![0.0; n];
    let mut tmp = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut kv = vec![0.0; k];
    let mut kv2 = vec![0.0; k];
    let mut kv3 = vec![0.0; k];
    let mut wv = vec![0.0; k];
    let mut a_old = vec![0.0; nn];

    let n_kappa = if inp.ordered { j_cat - 2 } else { 0 };
    let p = inp.n_par();
    let nsave = inp.iter / inp.thin;
    let mut draws = vec![0.0; nsave * p];
    let mut z_sum = vec![0.0; n];
    let mut v_sum = vec![0.0; n];
    let mut a_sum = vec![0.0; nn];
    let mut a_sq = vec![0.0; nn];
    let n_alpha = if inp.has_re && inp.alpha_every > 0 {
        nsave.div_ceil(inp.alpha_every)
    } else {
        0
    };
    let mut alpha_draws = vec![0.0; n_alpha * nn];

    let mut step = vec![0.5; j_cat + 1];
    let (mut acc, mut att) = (vec![0usize; j_cat + 1], vec![0usize; j_cat + 1]);
    let (mut acc_w, mut att_w) = (vec![0usize; j_cat + 1], vec![0usize; j_cat + 1]);

    let mut save = 0;
    let total = inp.warmup + inp.iter;
    for it in 0..total {
        if it % 100 == 0 {
            ctx.check()?;
        }
        ctx.set_progress(it + 1);

        // 1. latent block
        ch.refresh_lag();
        ch.ensure_resid();
        ch.latent_block(&mut prngs);
        if latent_ar {
            ch.refresh_lag();
            ch.refresh_design()?;
        }

        // 2. rho with beta integrated out (given lambda), then 3. beta.
        //    a0 = R (z - Delta alpha), a1 = R W z; the regression is
        //    a0 - rho a1 = Xt beta + eps. Everything needed is in the
        //    k-vectors g_i = Xt' D a_i and the products a_i' D a_j: with
        //    c0 = a0 - Xt beta_mean and G = Xt' D Xt = P - P0,
        //    Xt' D c0 = g0 - G beta_mean.
        a0.copy_from_slice(&ch.z);
        if inp.has_re {
            for t in 0..tt {
                for i in 0..nn {
                    a0[t * nn + i] -= ch.alpha[i];
                }
            }
        }
        if inp.has_lambda {
            inp.m.apply(&a0, &mut ch.work);
            for i in 0..n {
                a0[i] -= ch.lambda * ch.work[i];
            }
        }
        // D a0 (tmp), g0 (kv)
        let da0: &[f64] = if has_t {
            for i in 0..n {
                tmp[i] = ch.wt[i] * a0[i];
            }
            &tmp
        } else {
            &a0
        };
        xt_mul(&ch.xt, n, da0, &mut kv);
        if inp.has_rho {
            inp.w.apply(&ch.z, &mut a1);
            if inp.has_lambda {
                inp.m.apply(&a1, &mut ch.work);
                for i in 0..n {
                    a1[i] -= ch.lambda * ch.work[i];
                }
            }
            // D a1 (tmp2), g1 (kv2)
            let da1: &[f64] = if has_t {
                for i in 0..n {
                    tmp2[i] = ch.wt[i] * a1[i];
                }
                &tmp2
            } else {
                &a1
            };
            xt_mul(&ch.xt, n, da1, &mut kv2);
            let (s00, s01, s11) = (dot(&a0, da0), dot(&a0, da1), dot(&a1, da1));
            // G beta_mean = P beta_mean - P0 beta_mean (kv3)
            sq_matvec(&ch.prec, &inp.beta_mean, &mut kv3);
            for j in 0..k {
                kv3[j] -= ch.pc[j];
            }
            let bm = &inp.beta_mean;
            let c00 = s00 - 2.0 * dot(bm, &kv) + dot(bm, &kv3);
            let c01 = s01 - dot(bm, &kv2);
            // w0 = Lp^{-1} Xt' D c0 (kv3), w1 = Lp^{-1} g1 (wv)
            for j in 0..k {
                kv3[j] = kv[j] - kv3[j];
            }
            l_solve(&ch.lp, &mut kv3);
            wv.copy_from_slice(&kv2);
            l_solve(&ch.lp, &mut wv);
            let q00 = c00 - dot(&kv3, &kv3);
            let q01 = c01 - dot(&kv3, &wv);
            let q11 = s11 - dot(&wv, &wv);
            ch.rho = draw_grid_param(rng, &inp.rho, tt as f64, q00, q01, q11);
            // Xt' D (a0 - rho a1)
            for j in 0..k {
                kv[j] -= ch.rho * kv2[j];
            }
        }
        {
            // beta ~ N(P^{-1} (Xt' D y + P0 beta_mean), P^{-1})
            for j in 0..k {
                kv[j] += ch.pc[j];
            }
            l_solve(&ch.lp, &mut kv);
            for (j, v) in kv2.iter_mut().enumerate() {
                *v = kv[j] + rng.norm();
            }
            lt_solve(&ch.lp, &mut kv2);
            ch.beta.copy_from_slice(&kv2);
        }
        ch.resid_ok = false;

        // 4. lambda | z, rho, beta, alpha (plan/02 §2.5): with e = S z - X beta
        //    - Delta alpha, the kernel is |e - lambda M e|^2.
        if inp.has_lambda {
            // compute_resid leaves M e in work
            ch.compute_resid();
            for i in 0..n {
                tmp[i] = if has_t { ch.wt[i] * ch.e[i] } else { ch.e[i] };
                tmp2[i] = if has_t {
                    ch.wt[i] * ch.work[i]
                } else {
                    ch.work[i]
                };
            }
            let q00 = dot(&ch.e, &tmp);
            let q01 = dot(&tmp, &ch.work);
            let q11 = dot(&ch.work, &tmp2);
            ch.lambda = draw_grid_param(rng, &inp.lambda, tt as f64, q00, q01, q11);
            ch.resid_ok = false;
            ch.refresh_design()?;
            ch.filters.build_r(ch.lambda);
        }
        ch.filters.build_a(ch.rho, ch.lambda);

        // 5. random effects and their variance
        if inp.has_re {
            ch.compute_resid();
            a_old.copy_from_slice(&ch.alpha);
            if let Some(re) = re_block.as_mut() {
                // sparse precision: u_t = S_t z_t - X_t beta = e_t + alpha
                for t in 0..tt {
                    for i in 0..nn {
                        tmp[t * nn + i] = ch.e[t * nn + i] + ch.alpha[i];
                    }
                }
                let wt = has_t.then_some(&ch.wt[..]);
                re.draw(
                    inp,
                    &tmp,
                    wt,
                    ch.lambda,
                    ch.phi,
                    ch.sigma2,
                    &mut || rng.norm(),
                    &mut ch.alpha,
                )?;
            } else {
                // diagonal precision sum_t w_it + 1 / sigma2
                let mut prec = vec![1.0 / ch.sigma2; nn];
                let mut wu = vec![0.0; nn];
                for t in 0..tt {
                    for i in 0..nn {
                        let w = ch.wt[t * nn + i];
                        prec[i] += w;
                        wu[i] += w * (ch.e[t * nn + i] + ch.alpha[i]);
                    }
                }
                for i in 0..nn {
                    ch.alpha[i] = wu[i] / prec[i] + rng.norm() / prec[i].sqrt();
                }
            }
            // keep r current: alpha moved by alpha - a_old
            for i in 0..nn {
                a_old[i] = ch.alpha[i] - a_old[i];
            }
            ch.shift_resid_alpha(&a_old);
            // phi | alpha, sigma2: kernel |alpha - phi W alpha|^2 / sigma2 with log|B|
            if inp.spatial_re {
                let mut wa = vec![0.0; nn];
                inp.w.mul_one(0, &ch.alpha, &mut wa);
                let s = ch.sigma2;
                ch.phi = draw_grid_param(
                    rng,
                    &phi_grid,
                    1.0,
                    dot(&ch.alpha, &ch.alpha) / s,
                    dot(&ch.alpha, &wa) / s,
                    dot(&wa, &wa) / s,
                );
            }
            let ss = ch.re_ss();
            match inp.sigma_prior {
                SigmaPrior::HalfCauchy { scale } => {
                    // half-Cauchy(0, A) via the inverse-gamma mixture
                    ch.sigma2 = rinvgamma(rng, 0.5 * (nn as f64 + 1.0), 0.5 * ss + 1.0 / ch.xi);
                    ch.xi = rinvgamma(rng, 1.0, 1.0 / ch.sigma2 + 1.0 / (scale * scale));
                }
                SigmaPrior::InvGamma { shape, scale } => {
                    ch.sigma2 = rinvgamma(rng, shape + 0.5 * nn as f64, scale + 0.5 * ss);
                }
            }
        }

        // 5b. mixing variances of t errors: v_i | . ~ IG((nu + 1)/2, (nu + r_i^2)/2)
        if has_t {
            ch.ensure_resid();
            for i in 0..n {
                ch.wt[i] = rng.gamma(0.5 * (inp.nu + 1.0), 2.0 / (inp.nu + ch.r[i] * ch.r[i]));
            }
            ch.refresh_design()?;
        }

        // 6. cut-points
        if inp.ordered {
            ch.ensure_resid();
            ch.cut_ac_step(rng);
            if inp.cut_sampler == CutSampler::Transform {
                for j in 2..j_cat {
                    let ok = ch.cut_transform_move(rng, &ncat, j, step[j], &mut tmp, &mut tmp2);
                    att[j] += 1;
                    att_w[j] += 1;
                    if ok {
                        acc[j] += 1;
                        acc_w[j] += 1;
                    }
                    // adapt during warm-up only
                    if it < inp.warmup && att_w[j] == 50 {
                        let rate = acc_w[j] as f64 / 50.0;
                        step[j] *= if rate > 0.35 { 1.25 } else { 0.8 };
                        acc_w[j] = 0;
                        att_w[j] = 0;
                    }
                }
            }
        }

        // 7. scale move (generalised Gibbs, Liu & Sabatti 2000): x -> g x for
        //    x = (z, beta, alpha, kappa).
        if inp.scale_move {
            ch.refresh_lag();
            ch.ensure_resid();
            // gamma (latent AR) is not rescaled: its column z_{t-1} already is
            kv.copy_from_slice(&ch.beta);
            if let Some(c) = inp.lag_col {
                kv[c] = 0.0;
            }
            sq_matvec(&inp.beta_prec, &kv, &mut kv2);
            let rr = if has_t {
                (0..n).map(|i| ch.r[i] * (ch.wt[i] * ch.r[i])).sum::<f64>()
            } else {
                dot(&ch.r, &ch.r)
            };
            let mut c = rr + dot(&kv, &kv2);
            let mut l = dot(&kv, &ch.pc);
            let mut d = (n + k) as f64 - if latent_ar { 1.0 } else { 0.0 };
            if inp.has_re {
                c += ch.re_ss() / ch.sigma2;
                d += nn as f64;
            }
            if inp.ordered {
                for j in 2..j_cat {
                    if let KappaPrior::Normal { mean, sd } = inp.kappa_prior {
                        c += ch.kappa[j] * ch.kappa[j] / (sd * sd);
                        l += ch.kappa[j] * mean / (sd * sd);
                    }
                    d += 1.0;
                }
            }
            let mut g = 1.0;
            if l == 0.0 {
                g = rng.gamma(0.5 * d, 2.0 / c).sqrt();
            } else {
                // random walk on log g; the acceptance includes the Jacobian g^d
                let prop = (rng.norm() / (2.0 * d).sqrt()).exp();
                let logr = -0.5 * (prop * prop - 1.0) * c + (prop - 1.0) * l + d * prop.ln();
                if rng.unif().ln() < logr {
                    g = prop;
                }
            }
            // r = A z - b is linear in (z, beta, alpha) jointly, and gamma's
            // column z_{t-1} is rescaled with z, so r scales with g too
            for (z, r) in ch.z.iter_mut().zip(&mut ch.r) {
                *z *= g;
                *r *= g;
            }
            let gam = ch.gamma_now();
            for v in ch.beta.iter_mut() {
                *v *= g;
            }
            if let Some(c) = inp.lag_col {
                ch.beta[c] = gam;
            }
            for v in ch.alpha.iter_mut() {
                *v *= g;
            }
            if inp.ordered {
                for j in 2..j_cat {
                    ch.kappa[j] *= g;
                }
            }
        }

        // store
        if it >= inp.warmup && (it - inp.warmup) % inp.thin == 0 && save < nsave {
            let mut c = 0;
            let mut put = |v: f64| {
                draws[c * nsave + save] = v;
                c += 1;
            };
            for &b in &ch.beta {
                put(b);
            }
            if inp.has_rho {
                put(ch.rho);
            }
            if inp.has_lambda {
                put(ch.lambda);
            }
            if inp.has_re {
                put(ch.sigma2.sqrt());
            }
            if inp.spatial_re {
                put(ch.phi);
            }
            if inp.ordered {
                for j in 2..j_cat {
                    put(ch.kappa[j]);
                }
            }
            for i in 0..n {
                z_sum[i] += ch.z[i];
                if has_t {
                    v_sum[i] += 1.0 / ch.wt[i];
                }
            }
            for i in 0..nn {
                a_sum[i] += ch.alpha[i];
                a_sq[i] += ch.alpha[i] * ch.alpha[i];
            }
            if n_alpha > 0 && save % inp.alpha_every == 0 {
                let row = save / inp.alpha_every;
                for i in 0..nn {
                    alpha_draws[i * n_alpha + row] = ch.alpha[i];
                }
            }
            save += 1;
        }
        if inp.warmup > 0 && it == inp.warmup - 1 {
            acc.fill(0);
            att.fill(0);
        }
    }

    let s = save as f64;
    let alpha_mean: Vec<f64> = a_sum.iter().map(|a| a / s).collect();
    let alpha_sd = a_sq
        .iter()
        .zip(&alpha_mean)
        .map(|(q, m)| (q / s - m * m).max(0.0).sqrt())
        .collect();
    Ok(ChainOutput {
        draws,
        nsave,
        n_par: p,
        z_mean: z_sum.iter().map(|z| z / s).collect(),
        alpha_mean,
        alpha_sd,
        cut_accept: (2..2 + n_kappa)
            .map(|j| {
                if att[j] > 0 {
                    acc[j] as f64 / att[j] as f64
                } else {
                    f64::NAN
                }
            })
            .collect(),
        cut_step: (2..2 + n_kappa).map(|j| step[j]).collect(),
        alpha_draws,
        n_alpha,
        v_mean: v_sum.iter().map(|v| v / s).collect(),
    })
}

/// The cut-point updates on a toy ordered model, for the stationarity test
/// (`test-cutmove.R`; the C++ engine's `psp_cut_kernel_cpp`): latent
/// z_i ~ N(mu_i, 1), categories `y` (1..J), free cut-points `kappa_free`
/// with prior N(kmean, ksd) restricted to the ordering. Applies `reps` rounds
/// of the Albert-Chib step (`ac`) and/or the transformation move (`mv`), and
/// returns the final free cut-points and z.
#[allow(clippy::too_many_arguments)]
pub fn cut_kernel<R: Rng>(
    rng: &mut R,
    y: &[i32],
    mu: &[f64],
    z: &[f64],
    kappa_free: &[f64],
    n_cat: usize,
    kmean: f64,
    ksd: f64,
    step: f64,
    reps: usize,
    ac: bool,
    mv: bool,
) -> Result<(Vec<f64>, Vec<f64>), PspError> {
    use crate::sampler::input::GridPrior;
    use crate::sparse::Weights;
    let n = y.len();
    if mu.len() != n || z.len() != n || n_cat < 3 || kappa_free.len() != n_cat - 2 {
        return Err(PspError::InvalidInput(
            "cut kernel: inputs do not match".to_string(),
        ));
    }
    if y.iter().any(|&c| c < 1 || c as usize > n_cat) {
        return Err(PspError::InvalidInput(
            "cut kernel: categories must be in 1..J".to_string(),
        ));
    }
    let eye = Weights::new(vec![Csc::identity(n)], vec![0])?;
    let grid = GridPrior {
        a: 1.0,
        b: 1.0,
        min: -1.0,
        max: 1.0,
        grid: vec![-0.5, 0.5],
        lndet: vec![0.0, 0.0],
    };
    let inp = SamplerInput {
        y: y.to_vec(),
        x: Vec::new(),
        k: 0,
        n_units: n,
        periods: 1,
        n_cat,
        ordered: true,
        has_rho: false,
        has_lambda: false,
        has_re: false,
        spatial_re: false,
        w: eye.clone(),
        m: eye,
        beta_mean: Vec::new(),
        beta_prec: Vec::new(),
        rho: grid.clone(),
        lambda: grid,
        phi_a: 1.0,
        phi_b: 1.0,
        nu: 0.0,
        lag_col: None,
        sigma_prior: SigmaPrior::InvGamma {
            shape: 1.0,
            scale: 1.0,
        },
        kappa_prior: KappaPrior::Normal {
            mean: kmean,
            sd: ksd,
        },
        iter: 1,
        warmup: 0,
        thin: 1,
        z_sweeps: 1,
        threads: 1,
        cut_sampler: CutSampler::Transform,
        scale_move: false,
        beta_start: Vec::new(),
        rho_start: 0.0,
        lambda_start: 0.0,
        sigma2_start: 1.0,
        kappa_start: kappa_free.to_vec(),
        alpha_every: 0,
    };
    let mut kappa = vec![0.0; n_cat + 1];
    kappa[0] = -INF;
    kappa[n_cat] = INF;
    kappa[2..n_cat].copy_from_slice(kappa_free);
    let mut ncat = vec![0usize; n_cat + 1];
    for &c in y {
        ncat[c as usize] += 1;
    }
    let mut ch = Chain {
        inp: &inp,
        n,
        nn: n,
        k: 0,
        latent_ar: false,
        has_t: false,
        beta: Vec::new(),
        rho: 0.0,
        lambda: 0.0,
        sigma2: 1.0,
        phi: 0.0,
        xi: 1.0,
        alpha: Vec::new(),
        kappa,
        cat_idx: category_index(y, n_cat),
        z: z.to_vec(),
        wt: vec![1.0; n],
        xm: Vec::new(),
        xt: Vec::new(),
        prec: Vec::new(),
        lp: Vec::new(),
        pc: Vec::new(),
        mx: Vec::new(),
        gram: None,
        filters: Filters::new(&inp, false)?,
        e: vec![0.0; n],
        r: z.iter().zip(mu).map(|(z, m)| z - m).collect(),
        resid_ok: true,
        work: vec![0.0; n],
    };
    ch.filters.build_a(0.0, 0.0);
    let (mut delta, mut ad) = (vec![0.0; n], vec![0.0; n]);
    for _ in 0..reps {
        if ac {
            ch.cut_ac_step(rng);
        }
        if mv {
            for j in 2..n_cat {
                ch.cut_transform_move(rng, &ncat, j, step, &mut delta, &mut ad);
            }
        }
    }
    Ok((ch.kappa[2..n_cat].to_vec(), ch.z))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::par;
    use crate::rng::RMrg32k3a;
    use crate::sparse::Weights;
    use std::time::Duration;

    fn rook(nr: usize, nc: usize, skip: usize) -> Csc {
        // rook neighbours, minus every `skip`-th link (a second W for
        // time-varying weights), row-standardised
        let mut trip = Vec::new();
        for r in 0..nr {
            for c in 0..nc {
                let i = r * nc + c;
                let mut nb = Vec::new();
                for (ok, j) in [
                    (r > 0, i.wrapping_sub(nc)),
                    (r + 1 < nr, i + nc),
                    (c > 0, i.wrapping_sub(1)),
                    (c + 1 < nc, i + 1),
                ] {
                    if ok && (skip == 0 || (i + j) % skip != 0 || nb.is_empty()) {
                        nb.push(j);
                    }
                }
                for &j in &nb {
                    trip.push((i, j, 1.0 / nb.len() as f64));
                }
            }
        }
        Csc::from_triplets(nr * nc, nr * nc, &trip).unwrap()
    }

    #[derive(Clone, Copy, Default)]
    struct Variant {
        rho: bool,
        lambda: bool,
        re: bool,
        spatial_re: bool,
        ordered: bool,
        t_errors: bool,
        latent_ar: bool,
        varying_w: bool,
        scale_move: bool,
    }

    fn input(v: Variant) -> SamplerInput {
        let (nn, tt, k) = (20, 3, 3);
        let n = nn * tt;
        let mut rng = RMrg32k3a::from_seed([7; 6]).unwrap();
        let mut x = vec![1.0; n * k];
        for e in x[n..].iter_mut() {
            *e = rng.norm();
        }
        let y = (0..n)
            .map(|i| {
                let z = x[n + i] + 0.5 * rng.norm();
                if v.ordered {
                    1 + (z > -0.5) as i32 + (z > 0.5) as i32
                } else {
                    (z > 0.0) as i32
                }
            })
            .collect();
        let w = if v.varying_w {
            Weights::new(vec![rook(4, 5, 0), rook(4, 5, 3)], vec![0, 1, 0]).unwrap()
        } else {
            Weights::new(vec![rook(4, 5, 0)], vec![0; tt]).unwrap()
        };
        let grid: Vec<f64> = (0..=18).map(|i| -0.9 + 0.1 * i as f64).collect();
        let gp = GridPrior {
            a: 1.0,
            b: 1.0,
            min: -1.0,
            max: 1.0,
            lndet: grid.iter().map(|g| -0.3 * g * g).collect(),
            grid,
        };
        SamplerInput {
            y,
            x,
            k,
            n_units: nn,
            periods: tt,
            n_cat: if v.ordered { 3 } else { 2 },
            ordered: v.ordered,
            has_rho: v.rho,
            has_lambda: v.lambda,
            has_re: v.re,
            spatial_re: v.spatial_re,
            m: w.clone(),
            w,
            beta_mean: vec![0.1; k],
            beta_prec: (0..k * k)
                .map(|i| if i % (k + 1) == 0 { 0.1 } else { 0.0 })
                .collect(),
            rho: gp.clone(),
            lambda: gp,
            phi_a: 1.0,
            phi_b: 1.0,
            nu: if v.t_errors { 5.0 } else { 0.0 },
            lag_col: v.latent_ar.then_some(k - 1),
            sigma_prior: SigmaPrior::HalfCauchy { scale: 1.0 },
            kappa_prior: if v.ordered {
                KappaPrior::Normal { mean: 0.0, sd: 3.0 }
            } else {
                KappaPrior::Flat
            },
            iter: 60,
            warmup: 60,
            thin: 1,
            z_sweeps: 2,
            threads: 1,
            cut_sampler: CutSampler::Transform,
            scale_move: v.scale_move,
            beta_start: vec![0.0; k],
            rho_start: 0.0,
            lambda_start: 0.0,
            sigma2_start: 1.0,
            kappa_start: if v.ordered { vec![1.0] } else { Vec::new() },
            alpha_every: 0,
        }
    }

    /// Every block that keeps r = A z - b current in place agrees with a
    /// fresh computation (`ensure_resid` checks this in debug builds).
    #[test]
    fn residual_kept_current() {
        let base = Variant {
            rho: true,
            re: true,
            scale_move: true,
            ..Default::default()
        };
        let variants = [
            base,
            Variant {
                scale_move: false,
                ..base
            },
            Variant { re: false, ..base },
            Variant {
                ordered: true,
                ..base
            },
            Variant {
                t_errors: true,
                ordered: true,
                ..base
            },
            Variant {
                rho: false,
                lambda: true,
                ..base
            },
            Variant {
                lambda: true,
                t_errors: true,
                ..base
            },
            Variant {
                lambda: true,
                varying_w: true,
                ..base
            },
            Variant {
                latent_ar: true,
                ordered: true,
                ..base
            },
            Variant {
                latent_ar: true,
                lambda: true,
                t_errors: true,
                ..base
            },
            Variant {
                spatial_re: true,
                ..base
            },
            Variant {
                varying_w: true,
                ordered: true,
                scale_move: false,
                ..base
            },
        ];
        for (i, &v) in variants.iter().enumerate() {
            let inp = input(v);
            let out = par::run(
                1,
                1,
                Duration::from_millis(50),
                |_, ctx| {
                    let mut rng = RMrg32k3a::from_seed([1, 2, 3, 4, 5, 6])?;
                    run_chain(&inp, &mut rng, ctx)
                },
                |_| Ok::<(), ()>(()),
            )
            .unwrap_or_else(|_| panic!("variant {i} failed"));
            assert!(out[0].draws.iter().all(|d| d.is_finite()), "variant {i}");
        }
    }
}
