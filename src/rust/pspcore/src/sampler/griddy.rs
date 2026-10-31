//! Griddy Gibbs draws of rho, lambda and phi (plan/02 §2.4), ported from the
//! C++ engine's `lndet_interp`, `draw_from_grid` and `draw_rho`.

use super::input::GridPrior;
use crate::rng::Rng;

/// Linear interpolation of the log-determinant grid at `x`.
pub fn lndet_interp(grid: &[f64], lndet: &[f64], x: f64) -> f64 {
    let g = grid.len();
    let step = (grid[g - 1] - grid[0]) / (g - 1) as f64;
    let pos = (x - grid[0]) / step;
    let i = (pos.floor() as i64).clamp(0, g as i64 - 2) as usize;
    let w = pos - i as f64;
    (1.0 - w) * lndet[i] + w * lndet[i + 1]
}

/// A draw from a density known on a grid: a cell by trapezoid mass, then an
/// exact draw from the linear density inside it.
pub fn draw_from_grid<R: Rng>(rng: &mut R, x: &[f64], lp: &[f64]) -> f64 {
    let g = x.len();
    let mx = lp.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let f: Vec<f64> = lp.iter().map(|&l| (l - mx).exp()).collect();
    let mut cum = vec![0.0; g - 1];
    let mut tot = 0.0;
    for i in 0..g - 1 {
        tot += 0.5 * (f[i] + f[i + 1]) * (x[i + 1] - x[i]);
        cum[i] = tot;
    }
    let u = rng.unif() * tot;
    let mut c = 0;
    while c < g - 2 && cum[c] < u {
        c += 1;
    }
    let prev = if c == 0 { 0.0 } else { cum[c - 1] };
    let h = x[c + 1] - x[c];
    let (f0, f1) = (f[c], f[c + 1]);
    let mass = 0.5 * (f0 + f1) * h;
    let v = if mass > 0.0 {
        (u - prev) / mass
    } else {
        rng.unif()
    };
    let slope = (f1 - f0) / h;
    let s = if slope.abs() * h < 1e-10 * (f0 + f1) {
        v * h
    } else {
        (-f0 + (f0 * f0 + 2.0 * slope * v * mass).sqrt()) / slope
    };
    x[c] + s.clamp(0.0, h)
}

/// Collapsed griddy Gibbs draw of a spatial parameter whose kernel is
/// `T log|I - x W| - (q00 - 2 x q01 + x^2 q11) / 2` plus the Beta prior. A
/// coarse pass (every `stride`-th grid point) finds the region within 30 log
/// units of the mode, which is refined on 400 points.
pub fn draw_grid_param<R: Rng>(
    rng: &mut R,
    g: &GridPrior,
    t: f64,
    q00: f64,
    q01: f64,
    q11: f64,
) -> f64 {
    let lp = |x: f64, ld: f64| {
        t * ld - 0.5 * (q00 - 2.0 * x * q01 + x * x * q11)
            + (g.a - 1.0) * (x - g.min).ln()
            + (g.b - 1.0) * (g.max - x).ln()
    };
    let n = g.grid.len();
    let stride = (n / 500).max(1);
    let mut idx: Vec<usize> = (0..n).step_by(stride).collect();
    if *idx.last().unwrap_or(&0) != n - 1 {
        idx.push(n - 1);
    }
    let c = idx.len();
    let coarse: Vec<f64> = idx.iter().map(|&i| lp(g.grid[i], g.lndet[i])).collect();
    let mx = coarse.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let (mut first, mut last) = (c - 1, 0);
    for (j, &l) in coarse.iter().enumerate() {
        if l > mx - 30.0 {
            first = first.min(j);
            last = last.max(j);
        }
    }
    let first = idx[first.saturating_sub(2)];
    let last = idx[(last + 2).min(c - 1)];
    const R: usize = 400;
    let (lo, hi) = (g.grid[first], g.grid[last]);
    let xs: Vec<f64> = (0..R)
        .map(|i| lo + (hi - lo) * i as f64 / (R as f64 - 1.0))
        .collect();
    let lps: Vec<f64> = xs
        .iter()
        .map(|&x| lp(x, lndet_interp(&g.grid, &g.lndet, x)))
        .collect();
    draw_from_grid(rng, &xs, &lps)
}
