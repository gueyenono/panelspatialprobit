//! One MCMC chain: a port of `psp_sampler_cpp` in the C++ engine
//! (src/sampler.cpp at tag `cpp-final`), block for block and in the same
//! order, consuming random numbers in the same order (plan/06 §6.7).
//!
//! Every model is the kernel |A z - b|^2 with block-diagonal A = R S,
//! S_t = I - rho W_t, R_t = I - lambda M_t (plan/02 §2.2).

use super::griddy::draw_grid_param;
use super::input::{ChainOutput, CutSampler, GridPrior, KappaPrior, SamplerInput, SigmaPrior};
use super::re::ReBlock;
use crate::dense::{cholesky, dot, spd_inverse, sq_matvec, x_mul, xt_mul};
use crate::dist::rtnorm;
use crate::par::TaskCtx;
use crate::rng::Rng;
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
    z: Vec<f64>,
    wt: Vec<f64>,
    // design: Xm (X with the latent AR lag column), Xt = R Xm, B, chol(B)
    xm: Vec<f64>,
    xt: Vec<f64>,
    b: Vec<f64>,
    lb: Vec<f64>,
    pc: Vec<f64>,
    mx: Vec<f64>,
    gram: Option<[Vec<f64>; 3]>,
    filters: Filters,
    // work
    e: Vec<f64>,
    r: Vec<f64>,
    work: Vec<f64>,
}

impl<'a> Chain<'a> {
    /// |alpha|^2, or |alpha - phi W alpha|^2 for spatial random effects.
    fn re_ss(&self) -> f64 {
        if self.inp.spatial_re {
            let mut wa = vec![0.0; self.nn];
            self.inp.w.mats[0].mul_vec(&self.alpha, &mut wa);
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

    /// Xt = R Xm and B = (Xt' D Xt + P0)^{-1} with its Cholesky factor.
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
        self.b = spd_inverse(&xtx, k)?;
        self.lb = cholesky(&self.b, k)?;
        Ok(())
    }

    /// e = S z - Xm beta - Delta alpha and r = R e = A z - b.
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
    }

    fn bounds(&self, idx: usize) -> (f64, f64) {
        let y = self.inp.y[idx];
        if self.inp.ordered {
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

    /// One Gibbs sweep over the latent values with residual updates
    /// (plan/02 §2.3): z_k | . ~ TN(z_k - a_k'D r / h_k, 1 / h_k), h_k =
    /// a_k'D a_k, then r += a_k (z_k' - z_k). With latent AR, z_{k,t} also
    /// enters period t + 1 through column k of -gamma R_{t+1}.
    fn latent_sweep<R: Rng>(&mut self, rng: &mut R) {
        let (nn, inp) = (self.nn, self.inp);
        let gamma = self.gamma_now();
        for t in 0..inp.periods {
            let off = t * nn;
            let a = self.filters.a(t);
            let colsq = &self.filters.colsq[self.filters.pair[t]];
            let lagged = gamma != 0.0 && t < inp.periods - 1;
            let off1 = off + nn;
            let r_next = if lagged {
                self.filters.r(inp, t + 1)
            } else {
                None
            };
            for k in 0..nn {
                let idx = off + k;
                let (mut ar, mut h) = (0.0, 0.0);
                if self.has_t {
                    for (row, v) in a.col(k) {
                        let wa = self.wt[off + row] * v;
                        ar += wa * self.r[off + row];
                        h += wa * v;
                    }
                } else {
                    for (row, v) in a.col(k) {
                        ar += v * self.r[off + row];
                    }
                    h = colsq[k];
                }
                if lagged {
                    match r_next {
                        None => {
                            let wv = if self.has_t { self.wt[off1 + k] } else { 1.0 };
                            ar -= gamma * wv * self.r[off1 + k];
                            h += gamma * gamma * wv;
                        }
                        Some(rm) => {
                            for (row, v) in rm.col(k) {
                                let wv = if self.has_t { self.wt[off1 + row] } else { 1.0 };
                                ar -= gamma * wv * v * self.r[off1 + row];
                                h += gamma * gamma * wv * v * v;
                            }
                        }
                    }
                }
                let mean = self.z[idx] - ar / h;
                let (lo, hi) = self.bounds(idx);
                let znew = rtnorm(rng, mean, 1.0 / h.sqrt(), lo, hi);
                let delta = znew - self.z[idx];
                self.z[idx] = znew;
                for (row, v) in a.col(k) {
                    self.r[off + row] += v * delta;
                }
                if lagged {
                    match r_next {
                        None => self.r[off1 + k] -= gamma * delta,
                        Some(rm) => {
                            for (row, v) in rm.col(k) {
                                self.r[off1 + row] -= gamma * v * delta;
                            }
                        }
                    }
                }
            }
        }
    }

    /// Albert-Chib Gibbs step for the free cut-points given z.
    fn cut_ac_step<R: Rng>(&mut self, rng: &mut R) {
        let inp = self.inp;
        for j in 2..inp.n_cat {
            let (mut lo, mut hi) = (self.kappa[j - 1], self.kappa[j + 1]);
            for i in 0..self.n {
                let y = inp.y[i];
                if y == j as i32 && self.z[i] > lo {
                    lo = self.z[i];
                }
                if y == j as i32 + 1 && self.z[i] < hi {
                    hi = self.z[i];
                }
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
        for i in 0..n {
            let y = inp.y[i];
            delta[i] = if y == j as i32 {
                (l + (self.z[i] - l) * s1) - self.z[i]
            } else if y == j as i32 + 1 {
                if top {
                    knew - kj
                } else {
                    (u - (u - self.z[i]) * s2) - self.z[i]
                }
            } else {
                0.0
            };
        }
        for t in 0..inp.periods {
            let r = t * nn..(t + 1) * nn;
            self.filters.a(t).mul_vec(&delta[r.clone()], &mut ad[r]);
        }
        let gamma = self.gamma_now();
        if gamma != 0.0 {
            // latent AR: -gamma R_t delta_{t-1} in period t
            let mut tmp = vec![0.0; nn];
            for t in 1..inp.periods {
                let prev = &delta[(t - 1) * nn..t * nn];
                match self.filters.r(inp, t) {
                    None => {
                        for i in 0..nn {
                            ad[t * nn + i] -= gamma * prev[i];
                        }
                    }
                    Some(rm) => {
                        rm.mul_vec(prev, &mut tmp);
                        for i in 0..nn {
                            ad[t * nn + i] -= gamma * tmp[i];
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
            for i in 0..n {
                self.z[i] += delta[i];
                self.r[i] += ad[i];
            }
            self.kappa[j] = knew;
            true
        } else {
            false
        }
    }
}

/// Runs one chain. `ctx` gives the cancel flag (checked every 100 iterations)
/// and the progress counter (iterations done).
pub fn run_chain<R: Rng>(
    inp: &SamplerInput,
    rng: &mut R,
    ctx: &TaskCtx,
) -> Result<ChainOutput, PspError> {
    inp.validate()?;
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
        z: vec![0.0; n],
        wt: vec![1.0; n],
        xm: inp.x.clone(),
        xt: inp.x.clone(),
        b: Vec::new(),
        lb: Vec::new(),
        pc,
        mx: Vec::new(),
        gram: None,
        filters: Filters::new(inp, latent_ar)?,
        e: vec![0.0; n],
        r: vec![0.0; n],
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
        ch.compute_resid();
        for _ in 0..inp.z_sweeps {
            ch.latent_sweep(rng);
        }
        if latent_ar {
            ch.refresh_lag();
            ch.refresh_design()?;
        }

        // 2. rho with beta integrated out (given lambda), then 3. beta.
        //    a0 = R (z - Delta alpha), a1 = R W z; the regression is
        //    a0 - rho a1 = Xt beta + eps.
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
        if inp.has_rho {
            inp.w.apply(&ch.z, &mut a1);
            if inp.has_lambda {
                inp.m.apply(&a1, &mut ch.work);
                for i in 0..n {
                    a1[i] -= ch.lambda * ch.work[i];
                }
            }
            // c0 = a0 - Xt beta_mean (tmp), D c0 (tmp2), D a1 (work)
            x_mul(&ch.xt, n, &inp.beta_mean, &mut tmp);
            for i in 0..n {
                tmp[i] = a0[i] - tmp[i];
                tmp2[i] = if has_t { ch.wt[i] * tmp[i] } else { tmp[i] };
                ch.work[i] = if has_t { ch.wt[i] * a1[i] } else { a1[i] };
            }
            xt_mul(&ch.xt, n, &tmp2, &mut kv); // u0
            xt_mul(&ch.xt, n, &ch.work, &mut kv2); // u1
            sq_matvec(&ch.b, &kv, &mut kv3); // B u0
            let q00 = dot(&tmp, &tmp2) - dot(&kv, &kv3);
            let q01 = dot(&ch.work, &tmp) - dot(&kv2, &kv3);
            sq_matvec(&ch.b, &kv2, &mut kv3); // B u1
            let q11 = dot(&a1, &ch.work) - dot(&kv2, &kv3);
            ch.rho = draw_grid_param(rng, &inp.rho, tt as f64, q00, q01, q11);
            for i in 0..n {
                a0[i] -= ch.rho * a1[i];
            }
        }
        {
            for i in 0..n {
                tmp[i] = if has_t { ch.wt[i] * a0[i] } else { a0[i] };
            }
            xt_mul(&ch.xt, n, &tmp, &mut kv);
            for j in 0..k {
                kv[j] += ch.pc[j];
            }
            sq_matvec(&ch.b, &kv, &mut kv2); // mean
            for v in kv.iter_mut() {
                *v = rng.norm();
            }
            sq_matvec(&ch.lb, &kv, &mut kv3);
            for j in 0..k {
                ch.beta[j] = kv2[j] + kv3[j];
            }
        }

        // 4. lambda | z, rho, beta, alpha (plan/02 §2.5): with e = S z - X beta
        //    - Delta alpha, the kernel is |e - lambda M e|^2.
        if inp.has_lambda {
            ch.compute_resid();
            inp.m.apply(&ch.e, &mut ch.work);
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
            ch.refresh_design()?;
            ch.filters.build_r(ch.lambda);
        }
        ch.filters.build_a(ch.rho, ch.lambda);

        // 5. random effects and their variance
        if inp.has_re {
            ch.compute_resid();
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
            // phi | alpha, sigma2: kernel |alpha - phi W alpha|^2 / sigma2 with log|B|
            if inp.spatial_re {
                let mut wa = vec![0.0; nn];
                inp.w.mats[0].mul_vec(&ch.alpha, &mut wa);
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
            ch.compute_resid();
            for i in 0..n {
                ch.wt[i] = rng.gamma(0.5 * (inp.nu + 1.0), 2.0 / (inp.nu + ch.r[i] * ch.r[i]));
            }
            ch.refresh_design()?;
        }

        // 6. cut-points
        if inp.ordered {
            ch.compute_resid();
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
            ch.compute_resid();
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
            for v in ch.z.iter_mut() {
                *v *= g;
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
        z: z.to_vec(),
        wt: vec![1.0; n],
        xm: Vec::new(),
        xt: Vec::new(),
        b: Vec::new(),
        lb: Vec::new(),
        pc: Vec::new(),
        mx: Vec::new(),
        gram: None,
        filters: Filters::new(&inp, false)?,
        e: vec![0.0; n],
        r: z.iter().zip(mu).map(|(z, m)| z - m).collect(),
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
