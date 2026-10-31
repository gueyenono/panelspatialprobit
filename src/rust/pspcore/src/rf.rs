//! Reduced form of the latent outcome for one posterior draw (plan/06 §6.8),
//! the port of `.psp_reduced_form`, `.psp_rf_block` and `.psp_rf_latent_ar`
//! (kept in `R/reference_postproc.R` as the test oracle).
//!
//! With `S = I - rho W`, `R = I - lambda M` and `B = I - phi W`, the latent
//! outcome of a period has mean `m = S^{-1} X beta` and covariance
//! `sigma^2 S^{-1} V S^{-T} + K K'`, `K = S^{-1} R^{-1}`, `V = (B'B)^{-1}`
//! (`V = I` without spatial random effects). [`Rf`] holds `m`, the sd `s`
//! and the diagonals the impacts need. Matrices are column-major; vectors of
//! length n = N T are stacked time-slow, unit-fast.

use crate::dense::x_mul;
use crate::sparse::{Csc, FilterPattern, SparseLu, Weights};
use crate::PspError;

/// The data a reduced form needs.
pub struct RfModel {
    pub n_units: usize,
    pub periods: usize,
    /// n x k design matrix.
    pub x: Vec<f64>,
    pub k: usize,
    pub w: Weights,
    /// Error weights (W when the model has no lambda).
    pub m: Weights,
}

/// The parameters of one draw. `sig2` is `sigma_alpha^2` (0 when pooled).
#[derive(Debug, Clone, Copy)]
pub struct RfPar<'a> {
    pub beta: &'a [f64],
    pub sig2: f64,
    pub rho: f64,
    pub lambda: f64,
    pub phi: f64,
    pub gamma: f64,
}

/// An N x `cols` column-major matrix (Rademacher probes, or the N x Q mixing
/// variances of t errors).
#[derive(Debug, Clone, Copy)]
pub struct Cols<'a> {
    pub val: &'a [f64],
    pub cols: usize,
}

/// Reduced-form quantities, each of length n (`s`: n x `q`, one column per
/// mixture draw for t errors, else one). `d1 = diag(S^{-1})`,
/// `d2 = diag(S^{-1} W)`, `t1 = S^{-1} 1`, `t2 = S^{-1} W 1`, filled only
/// when the impacts are needed.
#[derive(Debug, Clone, Default)]
pub struct Rf {
    pub m: Vec<f64>,
    pub s: Vec<f64>,
    pub q: usize,
    pub d1: Vec<f64>,
    pub d2: Vec<f64>,
    pub t1: Vec<f64>,
    pub t2: Vec<f64>,
}

// Per-unit results of one block (N rows).
struct Block {
    m: Vec<f64>,
    s: Vec<f64>,
    q: usize,
    d: Option<[Vec<f64>; 4]>,
}

// A filter `I - r A` on its pattern, with its LU (of `I - cur A`).
struct Filter {
    pat: FilterPattern,
    mat: Csc,
    lu: SparseLu,
    colsq: Vec<f64>,
    cur: Option<f64>,
}

impl Filter {
    fn new(a: &Csc) -> Result<Self, PspError> {
        let pat = FilterPattern::new(a.n_rows, Some(a), None)?;
        let mat = pat.new_matrix();
        let lu = SparseLu::new(&mat)?;
        Ok(Self {
            colsq: vec![0.0; a.n_rows],
            pat,
            mat,
            lu,
            cur: None,
        })
    }

    // Factorises I - r A, unless it already holds that factorisation.
    fn factorize(&mut self, r: f64) -> Result<(), PspError> {
        if self.cur == Some(r) {
            return Ok(());
        }
        self.cur = None;
        self.pat.update(r, 0.0, &mut self.mat, &mut self.colsq);
        self.lu.factorize(&self.mat.val)?;
        self.cur = Some(r);
        Ok(())
    }

    fn solve(&mut self, b: &mut [f64], cols: usize, transpose: bool) -> Result<(), PspError> {
        self.lu.solve(b, cols, transpose)
    }
}

// The filters of one distinct matrix: S = I - rho W_u and B = I - phi W_u,
// or R = I - lambda M_u.
struct UnitW {
    w: Csc,
    w1: Vec<f64>,
    s: Filter,
    b: Filter,
}

/// The patterns and LUs of every distinct W and M, built once per call and
/// reused for every draw.
pub struct RfContext<'a> {
    model: &'a RfModel,
    ws: Vec<UnitW>,
    ms: Vec<Filter>,
}

fn identity(n: usize) -> Vec<f64> {
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        a[i * n + i] = 1.0;
    }
    a
}

// Row sums of squares of an n x cols matrix.
fn row_sq(a: &[f64], n: usize) -> Vec<f64> {
    let mut out = vec![0.0; n];
    for col in a.chunks_exact(n) {
        for (o, v) in out.iter_mut().zip(col) {
            *o += v * v;
        }
    }
    out
}

// Row means of probes * a (both n x p): the stochastic diagonal estimate.
fn probe_diag(probes: &[f64], a: &[f64], n: usize, p: usize) -> Vec<f64> {
    let mut out = vec![0.0; n];
    for (pc, ac) in probes.chunks_exact(n).zip(a.chunks_exact(n)) {
        for i in 0..n {
            out[i] += pc[i] * ac[i];
        }
    }
    for o in &mut out {
        *o /= p as f64;
    }
    out
}

/// Which filter of period t: `S_t = I - r W_t`, `R_t = I - r M_t` or
/// `B_t = I - r W_t` (spatial random effects).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterKind {
    S,
    R,
    B,
}

impl<'a> RfContext<'a> {
    /// Overwrites the N x `cols` matrix `x` with `F^{-1} x` (or `F^{-T} x`)
    /// for the filter `F = I - coef A` of period `t`. Each filter is
    /// factorised once per new coefficient.
    pub fn solve_filter(
        &mut self,
        kind: FilterKind,
        t: usize,
        coef: f64,
        x: &mut [f64],
        cols: usize,
        transpose: bool,
    ) -> Result<(), PspError> {
        let f = match kind {
            FilterKind::S => &mut self.ws[self.model.w.period[t]].s,
            FilterKind::B => &mut self.ws[self.model.w.period[t]].b,
            FilterKind::R => &mut self.ms[self.model.m.period[t]],
        };
        f.factorize(coef)?;
        f.solve(x, cols, transpose)
    }

    pub fn model(&self) -> &RfModel {
        self.model
    }

    /// The reduced-form sd of period `t` alone, without random effects
    /// (`.psp_rf_block(W_t, M_t, ., 0, rho, lambda, probes, FALSE, mix =
    /// mix)$s`): N x q. It does not depend on the linear index; the means are
    /// `S_t^{-1} xb` ([`Self::solve_filter`]).
    pub fn period_sd(
        &mut self,
        t: usize,
        rho: f64,
        lambda: f64,
        probes: Option<Cols>,
        mix: Option<Cols>,
    ) -> Result<(Vec<f64>, usize), PspError> {
        let (uw, um) = (self.model.w.period[t], self.model.m.period[t]);
        let b = self.block(uw, um, &[], 0, 0.0, rho, lambda, 0.0, probes, false, mix)?;
        Ok((b.s, b.q))
    }

    /// The indices of the distinct (W, M) pair of period `t`.
    pub fn period_pair(&self, t: usize) -> (usize, usize) {
        (self.model.w.period[t], self.model.m.period[t])
    }

    pub fn new(model: &'a RfModel) -> Result<Self, PspError> {
        let n = model.n_units;
        if model.w.n() != n || model.m.n() != n {
            return Err(PspError::InvalidInput(
                "reduced form: W and M must be N x N".to_string(),
            ));
        }
        if model.w.periods() != model.periods
            || model.m.periods() != model.periods
            || model.x.len() != n * model.periods * model.k
        {
            return Err(PspError::InvalidInput(
                "reduced form: X, W and M disagree on N or T".to_string(),
            ));
        }
        let ws = model
            .w
            .mats
            .iter()
            .map(|w| {
                let mut w1 = vec![0.0; n];
                w.mul_vec(&vec![1.0; n], &mut w1);
                Ok(UnitW {
                    w: w.clone(),
                    w1,
                    s: Filter::new(w)?,
                    b: Filter::new(w)?,
                })
            })
            .collect::<Result<_, PspError>>()?;
        let ms = model
            .m
            .mats
            .iter()
            .map(Filter::new)
            .collect::<Result<_, _>>()?;
        Ok(Self { model, ws, ms })
    }

    /// The reduced form for one draw. With `alpha` (the N random effects of
    /// the draw) the mean includes them and `s` leaves `sigma_alpha` out.
    pub fn reduced_form(
        &mut self,
        par: &RfPar,
        alpha: Option<&[f64]>,
        need_impacts: bool,
        probes: Option<Cols>,
        mix: Option<Cols>,
    ) -> Result<Rf, PspError> {
        let (nu, tt) = (self.model.n_units, self.model.periods);
        let n = nu * tt;
        if par.beta.len() != self.model.k {
            return Err(PspError::InvalidInput(
                "reduced form: beta has the wrong length".to_string(),
            ));
        }
        let mut xb = vec![0.0; n];
        x_mul(&self.model.x, n, par.beta, &mut xb);
        let mut sig2 = par.sig2;
        if let Some(a) = alpha {
            for t in 0..tt {
                for i in 0..nu {
                    xb[t * nu + i] += a[i];
                }
            }
            sig2 = 0.0;
        }
        if par.gamma != 0.0 {
            return self.latent_ar(&xb, sig2, par, probes, need_impacts, mix);
        }
        let varies = (par.rho != 0.0 && self.model.w.varies())
            || (par.lambda != 0.0 && self.model.m.varies());
        let (wp, mp) = (&self.model.w.period, &self.model.m.period);
        if !varies {
            let (uw, um) = (wp[0], mp[0]);
            let b = self.block(
                uw,
                um,
                &xb,
                tt,
                sig2,
                par.rho,
                par.lambda,
                par.phi,
                probes,
                need_impacts,
                mix,
            )?;
            return Ok(expand(b, nu, tt));
        }
        let mut blocks = Vec::with_capacity(tt);
        for t in 0..tt {
            let (uw, um) = (wp[t], mp[t]);
            blocks.push(self.block(
                uw,
                um,
                &xb[t * nu..(t + 1) * nu],
                1,
                sig2,
                par.rho,
                par.lambda,
                par.phi,
                probes,
                need_impacts,
                mix,
            )?);
        }
        Ok(stack(blocks, nu))
    }

    // `.psp_rf_block`: periods that share W_uw and M_um; `xb` is N x c.
    #[allow(clippy::too_many_arguments)]
    fn block(
        &mut self,
        uw: usize,
        um: usize,
        xb: &[f64],
        c: usize,
        sig2: f64,
        rho: f64,
        lambda: f64,
        phi: f64,
        probes: Option<Cols>,
        need_impacts: bool,
        mix: Option<Cols>,
    ) -> Result<Block, PspError> {
        let n = self.model.n_units;
        let phi = if sig2 == 0.0 { 0.0 } else { phi };
        let UnitW { w, w1, s, b } = &mut self.ws[uw];
        let r = &mut self.ms[um];
        if rho == 0.0 && lambda == 0.0 && phi == 0.0 {
            let (sv, q) = match mix {
                None => (vec![(1.0 + sig2).sqrt(); n], 1),
                Some(mx) => (mx.val.iter().map(|v| (sig2 + v).sqrt()).collect(), mx.cols),
            };
            let d = need_impacts.then(|| [vec![1.0; n], vec![0.0; n], vec![1.0; n], w1.clone()]);
            return Ok(Block {
                m: xb.to_vec(),
                s: sv,
                q,
                d,
            });
        }
        s.factorize(rho)?;
        if lambda != 0.0 {
            r.factorize(lambda)?;
        }
        if phi != 0.0 {
            b.factorize(phi)?;
        }
        let mut m = xb.to_vec();
        s.solve(&mut m, c, false)?;

        let Some(pr) = probes else {
            // exact: dense S^{-1}, K = S^{-1} R^{-1}, G = S^{-1} B^{-1}
            let mut sinv = identity(n);
            s.solve(&mut sinv, n, false)?;
            let kk = if lambda != 0.0 {
                let mut k = identity(n);
                r.solve(&mut k, n, false)?;
                s.solve(&mut k, n, false)?;
                Some(k)
            } else {
                None
            };
            let g_sq = if phi != 0.0 {
                let mut g = identity(n);
                b.solve(&mut g, n, false)?;
                s.solve(&mut g, n, false)?;
                row_sq(&g, n)
            } else {
                row_sq(&sinv, n)
            };
            let k = kk.as_deref().unwrap_or(&sinv);
            let (sv, q) = match mix {
                None => {
                    let k_sq = row_sq(k, n);
                    (
                        (0..n).map(|i| (sig2 * g_sq[i] + k_sq[i]).sqrt()).collect(),
                        1,
                    )
                }
                Some(mx) => {
                    // sig2 rowSums(G^2) + K^2 %*% mix
                    let mut out = vec![0.0; n * mx.cols];
                    for qq in 0..mx.cols {
                        let o = &mut out[qq * n..(qq + 1) * n];
                        let mv = &mx.val[qq * n..(qq + 1) * n];
                        for (j, kc) in k.chunks_exact(n).enumerate() {
                            for i in 0..n {
                                o[i] += kc[i] * kc[i] * mv[j];
                            }
                        }
                        for i in 0..n {
                            o[i] = (sig2 * g_sq[i] + o[i]).sqrt();
                        }
                    }
                    (out, mx.cols)
                }
            };
            let d = if need_impacts {
                let d1: Vec<f64> = (0..n).map(|i| sinv[i * n + i]).collect();
                let d2: Vec<f64> = (0..n)
                    .map(|i| w.col(i).map(|(j, v)| sinv[j * n + i] * v).sum())
                    .collect();
                let mut t1 = vec![0.0; n];
                let mut t2 = vec![0.0; n];
                for (j, col) in sinv.chunks_exact(n).enumerate() {
                    for i in 0..n {
                        t1[i] += col[i];
                        t2[i] += col[i] * w1[j];
                    }
                }
                Some([d1, d2, t1, t2])
            } else {
                None
            };
            return Ok(Block { m, s: sv, q, d });
        };

        // stochastic diagonals: diag(G) ~ mean over probes v of v * (G v).
        // One multi-RHS solve with S of the columns needed: (1, W 1, v, W v)
        // for the impacts, then UG, UK (when it differs from UG) and RV.
        let p = pr.cols;
        let probes = pr.val;
        let mut u = probes.to_vec();
        s.solve(&mut u, p, true)?; // S^{-T} v
        let uk = if lambda != 0.0 {
            let mut t = u.clone();
            r.solve(&mut t, p, true)?;
            r.solve(&mut t, p, false)?;
            Some(t)
        } else {
            None // UK = U
        };
        let ug = if phi != 0.0 {
            let mut t = u.clone();
            b.solve(&mut t, p, true)?;
            b.solve(&mut t, p, false)?;
            t
        } else {
            u.clone()
        };
        // with neither R nor B, UK = UG = U: solve it once
        let uk = match uk {
            Some(t) => Some(t),
            None if phi != 0.0 => Some(u),
            None => None,
        };
        let rv = match mix {
            None => None,
            Some(_) => {
                let mut t = probes.to_vec();
                if lambda != 0.0 {
                    r.solve(&mut t, p, false)?;
                }
                Some(t)
            }
        };
        let mut sol = Vec::with_capacity(n * (2 + 5 * p));
        if need_impacts {
            let mut wv = vec![0.0; n * p];
            for (src, dst) in probes.chunks_exact(n).zip(wv.chunks_exact_mut(n)) {
                w.mul_vec(src, dst);
            }
            sol.extend(std::iter::repeat_n(1.0, n));
            sol.extend_from_slice(w1);
            sol.extend_from_slice(probes);
            sol.extend_from_slice(&wv);
        }
        let at_ug = sol.len() / n;
        sol.extend_from_slice(&ug);
        let at_uk = if let Some(uk) = &uk {
            sol.extend_from_slice(uk);
            at_ug + p
        } else {
            at_ug
        };
        let at_rv = sol.len() / n;
        if let Some(rv) = &rv {
            sol.extend_from_slice(rv);
        }
        let ncols = sol.len() / n;
        s.solve(&mut sol, ncols, false)?;
        let col = |start: usize, len: usize| &sol[start * n..(start + len) * n];
        let d_ss = probe_diag(probes, col(at_ug, p), n, p);
        let d_kk = if at_uk == at_ug {
            d_ss.clone()
        } else {
            probe_diag(probes, col(at_uk, p), n, p)
        };
        let floor = |v: f64| v.max(1e-12).sqrt();
        let (sv, q) = match mix {
            None => ((0..n).map(|i| floor(sig2 * d_ss[i] + d_kk[i])).collect(), 1),
            Some(mx) => {
                // t errors: own error term exact, the others at their mean variance
                let d_k = probe_diag(probes, col(at_rv, p), n, p);
                let ev = mx.val.iter().sum::<f64>() / mx.val.len() as f64;
                let mut out = vec![0.0; n * mx.cols];
                for qq in 0..mx.cols {
                    for i in 0..n {
                        let v = mx.val[qq * n + i];
                        out[qq * n + i] = floor(
                            sig2 * d_ss[i]
                                + d_k[i] * d_k[i] * v
                                + (d_kk[i] - d_k[i] * d_k[i]).max(0.0) * ev,
                        );
                    }
                }
                (out, mx.cols)
            }
        };
        let d = need_impacts.then(|| {
            [
                probe_diag(probes, col(2, p), n, p),
                probe_diag(probes, col(2 + p, p), n, p),
                col(0, 1).to_vec(),
                col(1, 1).to_vec(),
            ]
        });
        Ok(Block { m, s: sv, q, d })
    }

    // `.psp_rf_latent_ar`: z_t = S^{-1}(X_t beta + alpha + gamma z_{t-1} +
    // R^{-1} eps_t), time-invariant weights; probes have P (T + 1) columns.
    fn latent_ar(
        &mut self,
        xb: &[f64],
        sig2: f64,
        par: &RfPar,
        probes: Option<Cols>,
        need_impacts: bool,
        mix: Option<Cols>,
    ) -> Result<Rf, PspError> {
        let (nu, tt) = (self.model.n_units, self.model.periods);
        let (rho, lambda, gamma) = (par.rho, par.lambda, par.gamma);
        let (uw, um) = (self.model.w.period[0], self.model.m.period[0]);
        let p = probes.map_or(0, |pr| pr.cols / (tt + 1));
        let pr0 = probes.map(|pr| Cols {
            val: &pr.val[..nu * p],
            cols: p,
        });
        let base = self.block(
            uw,
            um,
            &xb[..nu],
            1,
            0.0,
            rho,
            lambda,
            0.0,
            pr0,
            need_impacts,
            mix,
        )?;
        let has_b = par.phi != 0.0 && sig2 > 0.0;
        let UnitW { s, b, .. } = &mut self.ws[uw];
        let r = &mut self.ms[um];
        s.factorize(rho)?;
        if lambda != 0.0 {
            r.factorize(lambda)?;
        }
        if has_b {
            b.factorize(par.phi)?;
        }
        let mut m = vec![0.0; nu * tt];
        for t in 0..tt {
            let mut v: Vec<f64> = xb[t * nu..(t + 1) * nu].to_vec();
            if t > 0 {
                for i in 0..nu {
                    v[i] += gamma * m[(t - 1) * nu + i];
                }
            }
            s.solve(&mut v, 1, false)?;
            m[t * nu..(t + 1) * nu].copy_from_slice(&v);
        }
        let mut lag_var = vec![0.0; nu * tt];
        let mut re_var = vec![0.0; nu * tt];
        if p == 0 {
            let mut ph = identity(nu);
            if lambda != 0.0 {
                r.solve(&mut ph, nu, false)?;
            }
            s.solve(&mut ph, nu, false)?;
            let mut acc = vec![0.0; nu];
            for t in 0..tt {
                if t > 0 {
                    s.solve(&mut ph, nu, false)?;
                    let g2 = gamma.powi(2 * t as i32);
                    for (a, v) in acc.iter_mut().zip(row_sq(&ph, nu)) {
                        *a += g2 * v;
                    }
                }
                lag_var[t * nu..(t + 1) * nu].copy_from_slice(&acc);
            }
            if sig2 > 0.0 {
                let mut vh = identity(nu);
                if has_b {
                    b.solve(&mut vh, nu, false)?;
                }
                let mut h = vec![0.0; nu * nu];
                for t in 0..tt {
                    for (hv, v) in h.iter_mut().zip(&vh) {
                        *hv = v + gamma * *hv;
                    }
                    s.solve(&mut h, nu, false)?;
                    re_var[t * nu..(t + 1) * nu].copy_from_slice(&row_sq(&h, nu));
                }
            }
        } else {
            let pv = probes.map_or(&[][..], |pr| pr.val);
            let blk = |g: usize| &pv[g * nu * p..(g + 1) * nu * p];
            let mut x = vec![0.0; nu * p];
            let mut yv = vec![0.0; nu * p];
            let mut c0 = blk(1).to_vec();
            if has_b {
                b.solve(&mut c0, p, false)?;
            }
            let mean_sq = |a: &[f64]| -> Vec<f64> {
                row_sq(a, nu).into_iter().map(|v| v / p as f64).collect()
            };
            for t in 0..tt {
                if t > 0 {
                    let mut kw = blk(t + 1).to_vec();
                    if lambda != 0.0 {
                        r.solve(&mut kw, p, false)?;
                    }
                    s.solve(&mut kw, p, false)?;
                    for (xv, k) in x.iter_mut().zip(&kw) {
                        *xv += k;
                    }
                    s.solve(&mut x, p, false)?;
                    for xv in &mut x {
                        *xv *= gamma;
                    }
                }
                lag_var[t * nu..(t + 1) * nu].copy_from_slice(&mean_sq(&x));
                if sig2 > 0.0 {
                    for (y, c) in yv.iter_mut().zip(&c0) {
                        *y = c + gamma * *y;
                    }
                    s.solve(&mut yv, p, false)?;
                    re_var[t * nu..(t + 1) * nu].copy_from_slice(&mean_sq(&yv));
                }
            }
        }
        let ev = mix.map_or(1.0, |mx| mx.val.iter().sum::<f64>() / mx.val.len() as f64);
        let n = nu * tt;
        let q = base.q;
        let mut sv = vec![0.0; n * q];
        for qq in 0..q {
            for t in 0..tt {
                for i in 0..nu {
                    let k = t * nu + i;
                    let s0 = base.s[qq * nu + i];
                    sv[qq * n + k] = (s0 * s0 + ev * lag_var[k] + sig2 * re_var[k]).sqrt();
                }
            }
        }
        let mut out = Rf {
            m,
            s: sv,
            q,
            ..Default::default()
        };
        if let Some([d1, d2, t1, t2]) = base.d {
            out.d1 = d1.repeat(tt);
            out.d2 = d2.repeat(tt);
            out.t1 = t1.repeat(tt);
            out.t2 = t2.repeat(tt);
        }
        Ok(out)
    }
}

// One block for all T periods: per-unit quantities repeated over periods.
fn expand(b: Block, nu: usize, tt: usize) -> Rf {
    let n = nu * tt;
    let mut s = vec![0.0; n * b.q];
    for qq in 0..b.q {
        for t in 0..tt {
            s[qq * n + t * nu..qq * n + (t + 1) * nu].copy_from_slice(&b.s[qq * nu..(qq + 1) * nu]);
        }
    }
    let mut out = Rf {
        m: b.m,
        s,
        q: b.q,
        ..Default::default()
    };
    if let Some([d1, d2, t1, t2]) = b.d {
        out.d1 = d1.repeat(tt);
        out.d2 = d2.repeat(tt);
        out.t1 = t1.repeat(tt);
        out.t2 = t2.repeat(tt);
    }
    out
}

// One block per period, stacked.
fn stack(blocks: Vec<Block>, nu: usize) -> Rf {
    let tt = blocks.len();
    let n = nu * tt;
    let q = blocks[0].q;
    let mut out = Rf {
        m: Vec::with_capacity(n),
        s: vec![0.0; n * q],
        q,
        ..Default::default()
    };
    for (t, b) in blocks.into_iter().enumerate() {
        out.m.extend_from_slice(&b.m);
        for qq in 0..q {
            out.s[qq * n + t * nu..qq * n + (t + 1) * nu]
                .copy_from_slice(&b.s[qq * nu..(qq + 1) * nu]);
        }
        if let Some([d1, d2, t1, t2]) = b.d {
            out.d1.extend(d1);
            out.d2.extend(d2);
            out.t1.extend(t1);
            out.t2.extend(t2);
        }
    }
    out
}
