//! Log-determinant grids `log |I - r W|` for the griddy-Gibbs draws of rho,
//! lambda and phi (plan/06 §6.5 items 4-5). The grid itself, the choice of
//! method and the spline through the sparse points stay in R.

use crate::sparse::{Csc, FilterPattern, SparseLu};
use crate::PspError;

/// `log |I - r W|` at each of `points`, from one sparse LU per point on a
/// pattern ordered once.
pub fn lndet_sparse(w: &Csc, points: &[f64]) -> Result<Vec<f64>, PspError> {
    let n = w.n_rows;
    let pat = FilterPattern::new(n, Some(w), None)?;
    let mut s = pat.new_matrix();
    let mut colsq = vec![0.0; n];
    let mut lu = SparseLu::new_simplicial(&s)?;
    points
        .iter()
        .map(|&r| {
            pat.update(r, 0.0, &mut s, &mut colsq);
            lu.factorize(&s.val)?;
            lu.log_abs_det()
        })
        .collect()
}

/// `log |I - r W| = sum_i log |1 - r e_i|` over the eigenvalues `re + i im`
/// of W, at each of `grid`.
pub fn lndet_eigen(re: &[f64], im: &[f64], grid: &[f64]) -> Vec<f64> {
    grid.iter()
        .map(|&r| {
            let s: f64 = re
                .iter()
                .zip(im)
                .map(|(&a, &b)| ((1.0 - r * a).powi(2) + (r * b).powi(2)).ln())
                .sum();
            0.5 * s
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sparse_and_eigen_agree_on_a_ring() {
        // Row-standardised ring of n units: eigenvalues cos(2 pi k / n).
        let n = 9;
        let (mut colptr, mut rowind, mut val) = (vec![0], Vec::new(), Vec::new());
        for j in 0..n {
            let mut rows = [(j + n - 1) % n, (j + 1) % n];
            rows.sort();
            for i in rows {
                rowind.push(i);
                val.push(0.5);
            }
            colptr.push(rowind.len());
        }
        let w = Csc::from_parts(n, n, colptr, rowind, val).unwrap();
        let re: Vec<f64> = (0..n)
            .map(|k| (2.0 * std::f64::consts::PI * k as f64 / n as f64).cos())
            .collect();
        let im = vec![0.0; n];
        let pts = [-0.95, -0.3, 0.0, 0.5, 0.99];
        let a = lndet_sparse(&w, &pts).unwrap();
        let b = lndet_eigen(&re, &im, &pts);
        for (x, y) in a.iter().zip(&b) {
            assert!((x - y).abs() < 1e-12, "{x} vs {y}");
        }
    }
}
