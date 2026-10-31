//! Per-draw post-processing on top of the reduced form (plan/06 §6.8):
//! category probabilities (`predict()`), their derivatives and the impacts
//! (`psp_impacts()`). Ports of `.psp_density_weights`, `.psp_category_probs`
//! and the draw loop of `psp_impacts()`.

use crate::dist::{dnorm, pnorm};
use crate::rf::Rf;

/// The full cut-point vector `(-Inf, 0, free..., Inf)` of an ordered model.
pub fn cut_points(free: &[f64]) -> Vec<f64> {
    let mut k = Vec::with_capacity(free.len() + 3);
    k.push(f64::NEG_INFINITY);
    k.push(0.0);
    k.extend_from_slice(free);
    k.push(f64::INFINITY);
    k
}

fn phi(x: f64) -> f64 {
    dnorm(x, 0.0, 1.0, false)
}

fn cdf(x: f64) -> f64 {
    pnorm(x, 0.0, 1.0, true, false)
}

// Mean over the q mixture columns of f(m_i, s_iq), for every i.
fn mix_mean(m: &[f64], s: &[f64], q: usize, f: impl Fn(f64, f64) -> f64) -> Vec<f64> {
    let n = m.len();
    let mut out = vec![0.0; n];
    for col in s.chunks_exact(n).take(q) {
        for i in 0..n {
            out[i] += f(m[i], col[i]);
        }
    }
    if q > 1 {
        for o in &mut out {
            *o /= q as f64;
        }
    }
    out
}

/// n x ncat matrix of `d P(category) / d(index)` given the reduced form
/// (`kappa`: the full cut-point vector of an ordered model; `None` for a
/// binary one).
pub fn density_weights(m: &[f64], s: &[f64], q: usize, kappa: Option<&[f64]>) -> Vec<f64> {
    match kappa {
        None => mix_mean(m, s, q, |m, s| phi(m / s) / s),
        Some(kap) => (0..kap.len() - 1)
            .flat_map(|j| {
                mix_mean(m, s, q, |m, s| {
                    (phi((kap[j] - m) / s) - phi((kap[j + 1] - m) / s)) / s
                })
            })
            .collect(),
    }
}

/// n x ncat matrix of outcome probabilities given the reduced form.
pub fn category_probs(m: &[f64], s: &[f64], q: usize, kappa: Option<&[f64]>) -> Vec<f64> {
    match kappa {
        None => mix_mean(m, s, q, |m, s| cdf(m / s)),
        Some(kap) => (0..kap.len() - 1)
            .flat_map(|j| {
                mix_mean(m, s, q, |m, s| {
                    cdf((kap[j + 1] - m) / s) - cdf((kap[j] - m) / s)
                })
            })
            .collect(),
    }
}

/// A regressor with impacts: its column in the design (`beta`) and the
/// column of its spatial lag `W x` in Durbin models (`theta`).
#[derive(Debug, Clone, Copy)]
pub struct ImpactVar {
    pub beta: usize,
    pub theta: Option<usize>,
}

/// Direct, indirect and total impacts of one draw, given the reduced form
/// (with its impact diagonals) and the density weights `g` (n x ncat).
/// Layout: `v + V (c + C (e + 3 p))` for variable v, category c, effect e
/// (direct, indirect, total) and period p (0 = all periods, then each period
/// when `by_period`).
pub fn impacts(
    rf: &Rf,
    g: &[f64],
    n_cat: usize,
    beta: &[f64],
    vars: &[ImpactVar],
    n_units: usize,
    by_period: bool,
) -> Vec<f64> {
    let n = rf.m.len();
    let periods = n / n_units;
    let nv = vars.len();
    let np = 1 + if by_period { periods } else { 0 };
    let mut out = vec![0.0; nv * n_cat * 3 * np];
    let mut put = |v: usize, c: usize, p: usize, dir: f64, tot: f64| {
        let at = |e: usize| v + nv * (c + n_cat * (e + 3 * p));
        out[at(0)] = dir;
        out[at(1)] = tot - dir;
        out[at(2)] = tot;
    };
    let mut dir_i = vec![0.0; n];
    let mut tot_i = vec![0.0; n];
    for (v, var) in vars.iter().enumerate() {
        let b = beta[var.beta];
        let th = var.theta.map_or(0.0, |j| beta[j]);
        for i in 0..n {
            dir_i[i] = b * rf.d1[i] + th * rf.d2[i];
            tot_i[i] = b * rf.t1[i] + th * rf.t2[i];
        }
        for c in 0..n_cat {
            let gc = &g[c * n..(c + 1) * n];
            let (mut dir, mut tot) = (0.0, 0.0);
            for i in 0..n {
                dir += gc[i] * dir_i[i];
                tot += gc[i] * tot_i[i];
            }
            put(v, c, 0, dir / n as f64, tot / n as f64);
            if by_period {
                for t in 0..periods {
                    let (mut dir, mut tot) = (0.0, 0.0);
                    for i in t * n_units..(t + 1) * n_units {
                        dir += gc[i] * dir_i[i];
                        tot += gc[i] * tot_i[i];
                    }
                    put(v, c, t + 1, dir / n_units as f64, tot / n_units as f64);
                }
            }
        }
    }
    out
}
