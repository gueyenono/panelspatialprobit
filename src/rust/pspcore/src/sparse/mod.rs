//! Sparse matrices and spatial operators (plan/06 §6.4, §6.7).

mod chol;
mod csc;
mod pattern;
mod weights;

pub use chol::SparseChol;
pub use csc::Csc;
pub use pattern::FilterPattern;
pub use weights::Weights;

#[cfg(test)]
mod tests {
    use super::*;

    // Row-standardised rook weights on an r x c grid.
    fn rook(r: usize, c: usize) -> Csc {
        let n = r * c;
        let mut colptr = vec![0];
        let mut rowind = Vec::new();
        let mut val = Vec::new();
        let deg = |i: usize| {
            let (a, b) = (i / c, i % c);
            (a > 0) as usize + (a + 1 < r) as usize + (b > 0) as usize + (b + 1 < c) as usize
        };
        for j in 0..n {
            let (a, b) = (j / c, j % c);
            let mut nb = Vec::new();
            if a > 0 {
                nb.push(j - c);
            }
            if b > 0 {
                nb.push(j - 1);
            }
            if b + 1 < c {
                nb.push(j + 1);
            }
            if a + 1 < r {
                nb.push(j + c);
            }
            for i in nb {
                rowind.push(i);
                val.push(1.0 / deg(i) as f64);
            }
            colptr.push(rowind.len());
        }
        Csc::from_parts(n, n, colptr, rowind, val).unwrap()
    }

    fn dense_mul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
        let mut c = vec![0.0; n * n];
        for j in 0..n {
            for k in 0..n {
                for i in 0..n {
                    c[j * n + i] += a[k * n + i] * b[j * n + k];
                }
            }
        }
        c
    }

    #[test]
    fn validation_rejects_bad_parts() {
        assert!(Csc::from_parts(2, 2, vec![0, 1, 2], vec![1, 0], vec![1.0, 1.0]).is_ok());
        assert!(Csc::from_parts(2, 2, vec![0, 2, 2], vec![1, 0], vec![1.0, 1.0]).is_err());
        assert!(Csc::from_parts(2, 2, vec![0, 1, 2], vec![2, 0], vec![1.0, 1.0]).is_err());
        assert!(Csc::from_parts(2, 2, vec![0, 1], vec![0], vec![1.0]).is_err());
        assert!(Csc::from_parts(2, 2, vec![0, 1, 2], vec![1, 0], vec![1.0, f64::NAN]).is_err());
    }

    #[test]
    fn products_match_dense_algebra() {
        let w = rook(3, 4);
        let m = rook(4, 3);
        let n = w.n_rows;
        let (wd, md) = (w.to_dense(), m.to_dense());
        let mw = m.matmul(&w).unwrap().to_dense();
        let want = dense_mul(&md, &wd, n);
        for (a, b) in mw.iter().zip(&want) {
            assert!((a - b).abs() < 1e-15);
        }
        let x: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
        let (mut y, mut yt) = (vec![0.0; n], vec![0.0; n]);
        w.mul_vec(&x, &mut y);
        w.tmul_vec(&x, &mut yt);
        for i in 0..n {
            let yi: f64 = (0..n).map(|k| wd[k * n + i] * x[k]).sum();
            let yti: f64 = (0..n).map(|k| wd[i * n + k] * x[k]).sum();
            assert!((y[i] - yi).abs() < 1e-15 && (yt[i] - yti).abs() < 1e-15);
        }
    }

    #[test]
    fn filter_updates_match_dense_algebra() {
        let w = rook(3, 4);
        let m = rook(4, 3);
        let n = w.n_rows;
        let eye = Csc::identity(n).to_dense();
        let (wd, md) = (w.to_dense(), m.to_dense());
        for (with_m, rho, lambda) in [(false, 0.4, 0.0), (true, 0.4, -0.3), (true, -0.9, 0.7)] {
            let pat =
                FilterPattern::new(n, Some(&w), if with_m { Some(&m) } else { None }).unwrap();
            let mut a = pat.new_matrix();
            let mut colsq = vec![0.0; n];
            pat.update(rho, lambda, &mut a, &mut colsq);
            let s: Vec<f64> = eye.iter().zip(&wd).map(|(i, w)| i - rho * w).collect();
            let r: Vec<f64> = eye.iter().zip(&md).map(|(i, m)| i - lambda * m).collect();
            let want = if with_m { dense_mul(&r, &s, n) } else { s };
            let got = a.to_dense();
            for (g, w) in got.iter().zip(&want) {
                assert!((g - w).abs() < 1e-14);
            }
            for j in 0..n {
                let cs: f64 = (0..n).map(|i| want[j * n + i].powi(2)).sum();
                assert!((colsq[j] - cs).abs() < 1e-13);
            }
        }
    }

    #[test]
    fn weights_apply_per_period() {
        let w1 = rook(3, 3);
        let w2 = Csc::identity(9);
        let wts = Weights::new(vec![w1.clone(), w2], vec![0, 1, 0]).unwrap();
        let v: Vec<f64> = (0..27).map(|i| i as f64).collect();
        let mut out = vec![0.0; 27];
        wts.apply(&v, &mut out);
        let mut y = vec![0.0; 9];
        w1.mul_vec(&v[18..27], &mut y);
        assert_eq!(&out[18..27], &y[..]);
        assert_eq!(&out[9..18], &v[9..18]);
        assert!(Weights::new(vec![w1], vec![1]).is_err());
    }
}
