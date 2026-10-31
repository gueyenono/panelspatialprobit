//! Simulation from the fitted model (plan/06 §6.8): replicated outcomes for
//! posterior predictive checks (`psp_ppc()`) and outcome histories for
//! response paths (`psp_response_path()`). Each draw has its own random
//! number stream.

use crate::dense::x_mul;
use crate::rf::{Cols, FilterKind, RfContext, RfPar};
use crate::rng::Rng;
use crate::PspError;

/// Standard latent errors: normal, or Student t with `df` degrees of freedom.
fn rerr<R: Rng>(rng: &mut R, df: Option<f64>) -> f64 {
    match df {
        None => rng.norm(),
        Some(df) => rng.t(df),
    }
}

/// The outcome codes of latent values: binary 0/1 (`z > 0`), or ordered
/// categories 1..J from the full cut-point vector `(-Inf, 0, ..., Inf)`.
fn code(z: f64, kappa: Option<&[f64]>, left_open: bool) -> i32 {
    match kappa {
        None => (z > 0.0) as i32,
        Some(kap) => {
            let inner = &kap[1..kap.len() - 1];
            let below = if left_open {
                inner.iter().filter(|&&b| b < z).count()
            } else {
                inner.iter().filter(|&&b| b <= z).count()
            };
            1 + below as i32
        }
    }
}

/// One replicate of the outcomes (length n) for one draw: z = mu + e, with
/// `e_t = R_t^{-1} eps_t`, plus `gamma z_{t-1}` (latent AR), then
/// `z_t = S_t^{-1} z_t`.
pub fn ppc_draw<R: Rng>(
    ctx: &mut RfContext,
    par: &RfPar,
    alpha: Option<&[f64]>,
    kappa: Option<&[f64]>,
    df: Option<f64>,
    rng: &mut R,
) -> Result<Vec<i32>, PspError> {
    let (nu, tt) = (ctx.model().n_units, ctx.model().periods);
    let n = nu * tt;
    let mut mu = vec![0.0; n];
    x_mul(&ctx.model().x, n, par.beta, &mut mu);
    if let Some(a) = alpha {
        for t in 0..tt {
            for i in 0..nu {
                mu[t * nu + i] += a[i];
            }
        }
    }
    let mut z: Vec<f64> = mu.iter().map(|m| m + rerr(rng, df)).collect();
    let mut tmp = vec![0.0; nu];
    for t in 0..tt {
        let r = t * nu..(t + 1) * nu;
        if par.lambda != 0.0 {
            for i in 0..nu {
                tmp[i] = z[t * nu + i] - mu[t * nu + i];
            }
            ctx.solve_filter(FilterKind::R, t, par.lambda, &mut tmp, 1, false)?;
            for i in 0..nu {
                z[t * nu + i] = mu[t * nu + i] + tmp[i];
            }
        }
        if t > 0 && par.gamma != 0.0 {
            for i in 0..nu {
                z[t * nu + i] += par.gamma * z[(t - 1) * nu + i];
            }
        }
        if par.rho != 0.0 {
            ctx.solve_filter(FilterKind::S, t, par.rho, &mut z[r], 1, false)?;
        }
    }
    Ok(z.iter().map(|&v| code(v, kappa, false)).collect())
}

/// What a response path needs besides the model and the draw.
pub struct ResponseSpec<'a> {
    /// Design columns of the lagged outcomes (rebuilt from the simulated
    /// histories).
    pub dynamic: &'a [usize],
    /// Per outcome level (ordered: 2..J; binary: 1), the columns of
    /// `lag_*` and `W_lag_*`, if in the design.
    pub lag: &'a [Option<usize>],
    pub wlag: &'a [Option<usize>],
    /// The shocked regressor, and its spatial lag `W_x` if in the design.
    pub var: usize,
    pub wvar: Option<usize>,
    /// The shock per unit.
    pub shock: &'a [f64],
    /// First period (0-based) and whether the shock is permanent.
    pub start: usize,
    pub permanent: bool,
    /// Unit groups (0-based) whose mean responses are reported.
    pub groups: &'a [Vec<usize>],
    pub nsim: usize,
    pub n_cat: usize,
    pub probes: Option<Cols<'a>>,
    pub mix: Option<Cols<'a>>,
    pub df: Option<f64>,
}

// N x ncat outcome probabilities, averaged over the nsim columns of m
// (N x nsim) and the q columns of s (N x q).
fn sim_probs(m: &[f64], s: &[f64], q: usize, nu: usize, kappa: Option<&[f64]>) -> Vec<f64> {
    use crate::dist::pnorm;
    let cdf = |x: f64| pnorm(x, 0.0, 1.0, true, false);
    let nsim = m.len() / nu;
    let ncat = kappa.map_or(1, |k| k.len() - 1);
    let mut out = vec![0.0; nu * ncat];
    for j in 0..ncat {
        for qq in 0..q {
            let sc = &s[qq * nu..(qq + 1) * nu];
            for i in 0..nu {
                let mut acc = 0.0;
                for c in 0..nsim {
                    let mi = m[c * nu + i];
                    acc += match kappa {
                        None => cdf(mi / sc[i]),
                        Some(k) => cdf((k[j + 1] - mi) / sc[i]) - cdf((k[j] - mi) / sc[i]),
                    };
                }
                out[j * nu + i] += acc / nsim as f64;
            }
        }
        for i in 0..nu {
            out[j * nu + i] /= q as f64;
        }
    }
    out
}

/// Mean differences in outcome probabilities between the shocked and the
/// baseline histories, for one draw: layout `h + H (c + C g)` for horizon h,
/// category c and unit group g. Both histories use the same random numbers
/// (common random numbers); each period's probabilities are exact given the
/// simulated history of the previous periods.
pub fn response_draw<R: Rng>(
    ctx: &mut RfContext,
    par: &RfPar,
    alpha: Option<&[f64]>,
    kappa: Option<&[f64]>,
    spec: &ResponseSpec,
    rng: &mut R,
) -> Result<Vec<f64>, PspError> {
    let (nu, tt, k) = (ctx.model().n_units, ctx.model().periods, ctx.model().k);
    let n = nu * tt;
    let nsim = spec.nsim;
    let horizons = tt - spec.start;
    let ngroups = spec.groups.len();
    let mut out = vec![0.0; horizons * spec.n_cat * ngroups];
    let beta = par.beta;
    let mut is_dyn = vec![false; k];
    for &j in spec.dynamic {
        is_dyn[j] = true;
    }
    let levels: Vec<i32> = if kappa.is_some() {
        (2..=spec.n_cat as i32).collect()
    } else {
        vec![1]
    };
    // histories of the previous period (N x nsim codes), baseline and shocked
    let mut hist: [Vec<i32>; 2] = [Vec::new(), Vec::new()];
    let mut ind = vec![0.0; nu];
    let mut wind = vec![0.0; nu];
    let mut we = vec![0.0; nu];
    // reduced-form sd per distinct (W, M) pair: the same in every period
    // that uses it, and for both histories
    let mut sds: Vec<((usize, usize), Vec<f64>, usize)> = Vec::new();
    for h in 0..horizons {
        let p = spec.start + h;
        let x = &ctx.model().x;
        let xrow = |j: usize, i: usize| x[j * n + p * nu + i];
        let mut fixed = vec![0.0; nu];
        let mut lagged0 = vec![0.0; nu];
        for j in 0..k {
            for i in 0..nu {
                if is_dyn[j] {
                    lagged0[i] += xrow(j, i) * beta[j];
                } else {
                    fixed[i] += xrow(j, i) * beta[j];
                }
            }
        }
        if let Some(a) = alpha {
            for i in 0..nu {
                fixed[i] += a[i];
            }
        }
        let mut shift = vec![0.0; nu];
        if h == 0 || spec.permanent {
            for i in 0..nu {
                shift[i] = spec.shock[i] * beta[spec.var];
            }
            if let Some(wv) = spec.wvar {
                ctx.model().w.at(p).mul_vec(spec.shock, &mut we);
                for i in 0..nu {
                    shift[i] += we[i] * beta[wv];
                }
            }
        }
        let mut eta: [Vec<f64>; 2] = [vec![0.0; nu * nsim], vec![0.0; nu * nsim]];
        for (sc, e) in eta.iter_mut().enumerate() {
            for c in 0..nsim {
                let col = &mut e[c * nu..(c + 1) * nu];
                for i in 0..nu {
                    col[i] = fixed[i] + sc as f64 * shift[i];
                }
                if h == 0 {
                    for i in 0..nu {
                        col[i] += lagged0[i];
                    }
                } else {
                    let y = &hist[sc][c * nu..(c + 1) * nu];
                    for (l, &lev) in levels.iter().enumerate() {
                        for i in 0..nu {
                            ind[i] = (y[i] == lev) as i32 as f64;
                        }
                        if let Some(j) = spec.lag[l] {
                            for i in 0..nu {
                                col[i] += beta[j] * ind[i];
                            }
                        }
                        if let Some(j) = spec.wlag[l] {
                            ctx.model().w.at(p).mul_vec(&ind, &mut wind);
                            for i in 0..nu {
                                col[i] += beta[j] * wind[i];
                            }
                        }
                    }
                }
            }
        }
        let pair = ctx.period_pair(p);
        let k_sd = match sds.iter().position(|(pr, _, _)| *pr == pair) {
            Some(k) => k,
            None => {
                let (s, q) = ctx.period_sd(p, par.rho, par.lambda, spec.probes, spec.mix)?;
                sds.push((pair, s, q));
                sds.len() - 1
            }
        };
        let mut probs: [Vec<f64>; 2] = [Vec::new(), Vec::new()];
        for sc in 0..2 {
            let mut m = eta[sc].clone();
            if par.rho != 0.0 || par.lambda != 0.0 {
                ctx.solve_filter(FilterKind::S, p, par.rho, &mut m, nsim, false)?;
            }
            let (_, s, q) = &sds[k_sd];
            probs[sc] = sim_probs(&m, s, *q, nu, kappa);
        }
        for c in 0..spec.n_cat {
            for (g, units) in spec.groups.iter().enumerate() {
                let mut acc = 0.0;
                for &i in units {
                    acc += probs[1][c * nu + i] - probs[0][c * nu + i];
                }
                out[h + horizons * (c + spec.n_cat * g)] = acc / units.len() as f64;
            }
        }
        if h + 1 < horizons {
            // common latent errors: z = S^{-1}(eta + R^{-1} eps)
            let mut eps: Vec<f64> = (0..nu * nsim).map(|_| rerr(rng, spec.df)).collect();
            if par.lambda != 0.0 {
                ctx.solve_filter(FilterKind::R, p, par.lambda, &mut eps, nsim, false)?;
            }
            let full = kappa.map(|k| k.to_vec());
            for sc in 0..2 {
                let mut z: Vec<f64> = eta[sc].iter().zip(&eps).map(|(a, b)| a + b).collect();
                if par.rho != 0.0 {
                    ctx.solve_filter(FilterKind::S, p, par.rho, &mut z, nsim, false)?;
                }
                hist[sc] = z.iter().map(|&v| code(v, full.as_deref(), true)).collect();
            }
        }
    }
    Ok(out)
}
