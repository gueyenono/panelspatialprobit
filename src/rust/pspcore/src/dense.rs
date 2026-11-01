//! Small dense helpers for the k x k regression block and the n x k design
//! (plan/06 §6.4). Matrices are column-major `Vec<f64>`.

use crate::PspError;

/// Lower Cholesky factor of a symmetric positive-definite k x k matrix.
pub fn cholesky(a: &[f64], k: usize) -> Result<Vec<f64>, PspError> {
    let mut l = vec![0.0; k * k];
    for j in 0..k {
        let mut d = a[j * k + j];
        for p in 0..j {
            d -= l[p * k + j] * l[p * k + j];
        }
        // also rejects NaN
        if d.is_nan() || d <= 0.0 {
            return Err(PspError::InvalidInput(
                "matrix is not positive definite".to_string(),
            ));
        }
        let d = d.sqrt();
        l[j * k + j] = d;
        for i in j + 1..k {
            let mut s = a[j * k + i];
            for p in 0..j {
                s -= l[p * k + i] * l[p * k + j];
            }
            l[j * k + i] = s / d;
        }
    }
    Ok(l)
}

/// Inverse of a symmetric positive-definite k x k matrix, via its Cholesky
/// factor (symmetric to rounding).
pub fn spd_inverse(a: &[f64], k: usize) -> Result<Vec<f64>, PspError> {
    let l = cholesky(a, k)?;
    // Columns of L^{-1}, then inv = L^{-T} L^{-1}.
    let mut linv = vec![0.0; k * k];
    for c in 0..k {
        for i in c..k {
            let mut s = if i == c { 1.0 } else { 0.0 };
            for p in c..i {
                s -= l[p * k + i] * linv[c * k + p];
            }
            linv[c * k + i] = s / l[i * k + i];
        }
    }
    let mut inv = vec![0.0; k * k];
    for j in 0..k {
        for i in j..k {
            let mut s = 0.0;
            for p in i..k {
                s += linv[i * k + p] * linv[j * k + p];
            }
            inv[j * k + i] = s;
            inv[i * k + j] = s;
        }
    }
    Ok(inv)
}

/// Solves `L x = b` in place, `L` lower triangular k x k.
pub fn l_solve(l: &[f64], x: &mut [f64]) {
    let k = x.len();
    for j in 0..k {
        x[j] /= l[j * k + j];
        let xj = x[j];
        for i in j + 1..k {
            x[i] -= l[j * k + i] * xj;
        }
    }
}

/// Solves `L' x = b` in place, `L` lower triangular k x k.
pub fn lt_solve(l: &[f64], x: &mut [f64]) {
    let k = x.len();
    for j in (0..k).rev() {
        let mut s = x[j];
        for i in j + 1..k {
            s -= l[j * k + i] * x[i];
        }
        x[j] = s / l[j * k + j];
    }
}

/// `y = A x` for a k x k matrix.
pub fn sq_matvec(a: &[f64], x: &[f64], y: &mut [f64]) {
    let k = x.len();
    y.fill(0.0);
    for (j, &xj) in x.iter().enumerate() {
        for i in 0..k {
            y[i] += a[j * k + i] * xj;
        }
    }
}

/// `y = X b` for an n x k design.
pub fn x_mul(x: &[f64], n: usize, b: &[f64], y: &mut [f64]) {
    y.fill(0.0);
    for (j, &bj) in b.iter().enumerate() {
        let col = &x[j * n..(j + 1) * n];
        for i in 0..n {
            y[i] += col[i] * bj;
        }
    }
}

/// `u = X' v` for an n x k design.
pub fn xt_mul(x: &[f64], n: usize, v: &[f64], u: &mut [f64]) {
    for (j, uj) in u.iter_mut().enumerate() {
        *uj = dot(&x[j * n..(j + 1) * n], v);
    }
}

/// `a . b`, with four independent partial sums (which the compiler
/// vectorises; Eigen's dot product also sums in lanes).
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len().min(b.len());
    let (a, b) = (&a[..n], &b[..n]);
    let mut acc = [0.0; 4];
    let (ca, cb) = (a.chunks_exact(4), b.chunks_exact(4));
    let (ra, rb) = (ca.remainder(), cb.remainder());
    for (x, y) in ca.zip(cb) {
        for l in 0..4 {
            acc[l] += x[l] * y[l];
        }
    }
    let mut s = (acc[0] + acc[2]) + (acc[1] + acc[3]);
    for (x, y) in ra.iter().zip(rb) {
        s += x * y;
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse_and_cholesky() {
        let a = [4.0, 2.0, 0.6, 2.0, 5.0, 1.0, 0.6, 1.0, 3.0];
        let inv = spd_inverse(&a, 3).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let s: f64 = (0..3).map(|p| a[p * 3 + i] * inv[j * 3 + p]).sum();
                assert!((s - if i == j { 1.0 } else { 0.0 }).abs() < 1e-14);
            }
        }
        let l = cholesky(&a, 3).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let s: f64 = (0..3).map(|p| l[p * 3 + i] * l[p * 3 + j]).sum();
                assert!((s - a[j * 3 + i]).abs() < 1e-14);
            }
        }
        assert!(cholesky(&[1.0, 2.0, 2.0, 1.0], 2).is_err());
        // L L' x = b through the two triangular solves
        let b = [1.0, -2.0, 0.5];
        let mut x = b;
        l_solve(&l, &mut x);
        lt_solve(&l, &mut x);
        let mut ax = [0.0; 3];
        sq_matvec(&a, &x, &mut ax);
        for i in 0..3 {
            assert!((ax[i] - b[i]).abs() < 1e-14);
        }
    }
}
