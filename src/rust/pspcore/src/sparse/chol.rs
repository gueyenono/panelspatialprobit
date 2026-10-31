//! Sparse Cholesky factorisation with a fixed pattern (plan/06 §6.7): the
//! symbolic analysis (AMD ordering) is done once by faer; each new set of
//! values is factorised numerically into the same buffers. The factor is
//! simplicial, so L is an ordinary CSC matrix and the triangular solves are
//! done here.

use super::Csc;
use crate::PspError;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::perm::PermRef;
use faer::sparse::linalg::cholesky::{
    factorize_symbolic_cholesky, CholeskySymbolicParams, SymbolicCholesky, SymbolicCholeskyRaw,
    SymmetricOrdering,
};
use faer::sparse::linalg::SupernodalThreshold;
use faer::sparse::{SparseColMatRef, SymbolicSparseColMat};
use faer::{Par, Side};

/// A symmetric positive-definite matrix given by its lower triangle on a
/// fixed pattern, and its Cholesky factor `Q A Q' = L L'`.
pub struct SparseChol {
    n: usize,
    lower: SymbolicSparseColMat<usize>,
    symbolic: SymbolicCholesky<usize>,
    l_colptr: Vec<usize>,
    l_rowind: Vec<usize>,
    l_val: Vec<f64>,
    perm_fwd: Vec<usize>,
    perm_inv: Vec<usize>,
    buf: MemBuffer,
    work: Vec<f64>,
}

impl SparseChol {
    /// Analyses the pattern of the lower triangle `lower` (row >= column
    /// entries only; the values are not used), with an AMD ordering.
    pub fn new(lower: &Csc) -> Result<Self, PspError> {
        Self::build(lower, SymmetricOrdering::Amd)
    }

    /// As [`Self::new`], keeping the given order (`Q = I`): for callers that
    /// order the coordinates themselves (GHK). The identity is passed as a
    /// custom ordering: faer 0.24.4's `SymmetricOrdering::Identity` skips the
    /// conversion of a lower-triangle input and misreads it.
    pub fn new_unordered(lower: &Csc) -> Result<Self, PspError> {
        let id: Vec<usize> = (0..lower.n_rows).collect();
        let perm = PermRef::new_checked(&id, &id, lower.n_rows);
        Self::build(lower, SymmetricOrdering::Custom(perm))
    }

    /// The factor L (simplicial CSC; the diagonal is the first entry of each
    /// column): column pointers, row indices and values.
    pub fn factor(&self) -> (&[usize], &[usize], &[f64]) {
        (&self.l_colptr, &self.l_rowind, &self.l_val)
    }

    fn build(lower: &Csc, ordering: SymmetricOrdering<usize>) -> Result<Self, PspError> {
        let n = lower.n_rows;
        let fail =
            |e: &dyn std::fmt::Debug| PspError::InvalidInput(format!("sparse Cholesky: {e:?}"));
        let pat = SymbolicSparseColMat::new_checked(
            n,
            n,
            lower.colptr.clone(),
            None,
            lower.rowind.clone(),
        );
        let symbolic = factorize_symbolic_cholesky(
            pat.as_ref(),
            Side::Lower,
            ordering,
            CholeskySymbolicParams {
                supernodal_flop_ratio_threshold: SupernodalThreshold::FORCE_SIMPLICIAL,
                ..Default::default()
            },
        )
        .map_err(|e| fail(&e))?;
        let (l_colptr, l_rowind) = match symbolic.raw() {
            SymbolicCholeskyRaw::Simplicial(s) => {
                let f = s.factor();
                (f.col_ptr().to_vec(), f.row_idx().to_vec())
            }
            SymbolicCholeskyRaw::Supernodal(_) => {
                return Err(PspError::InvalidInput(
                    "sparse Cholesky: expected a simplicial factor".to_string(),
                ))
            }
        };
        let (perm_fwd, perm_inv) = match symbolic.perm() {
            Some(p) => {
                let (f, i) = p.arrays();
                (f.to_vec(), i.to_vec())
            }
            None => ((0..n).collect(), (0..n).collect()),
        };
        let l_val = vec![0.0; symbolic.len_val()];
        let buf = MemBuffer::try_new(
            symbolic.factorize_numeric_llt_scratch::<f64>(Par::Seq, Default::default()),
        )
        .map_err(|e| fail(&e))?;
        Ok(Self {
            n,
            lower: pat,
            symbolic,
            l_colptr,
            l_rowind,
            l_val,
            perm_fwd,
            perm_inv,
            buf,
            work: vec![0.0; n],
        })
    }

    /// Factorises the matrix whose lower-triangle values, on the analysed
    /// pattern, are `values`.
    pub fn factorize(&mut self, values: &[f64]) -> Result<(), PspError> {
        let a = SparseColMatRef::new(self.lower.as_ref(), values);
        self.symbolic
            .factorize_numeric_llt::<f64>(
                &mut self.l_val,
                a,
                Side::Lower,
                Default::default(),
                Par::Seq,
                MemStack::new(&mut self.buf),
                Default::default(),
            )
            .map_err(|_| PspError::InvalidInput("random-effects Cholesky failed".to_string()))?;
        Ok(())
    }

    /// `x = A^{-1} b`.
    pub fn solve(&mut self, b: &[f64], x: &mut [f64]) {
        let n = self.n;
        // y = Q b, then L y' = y, then L' y'' = y', then x = Q' y''.
        for i in 0..n {
            self.work[i] = b[self.perm_fwd[i]];
        }
        for j in 0..n {
            let a = self.l_colptr[j];
            let d = self.work[j] / self.l_val[a];
            self.work[j] = d;
            for k in a + 1..self.l_colptr[j + 1] {
                self.work[self.l_rowind[k]] -= self.l_val[k] * d;
            }
        }
        self.lt_solve_in_place();
        for i in 0..n {
            x[i] = self.work[self.perm_inv[i]];
        }
    }

    /// `x = Q' L'^{-1} z`: for z standard normal, x ~ N(0, A^{-1}). This is the
    /// C++ engine's `permutationPinv() * matrixU().solve(z)`.
    pub fn lt_solve_perm(&mut self, z: &[f64], x: &mut [f64]) {
        self.work.copy_from_slice(z);
        self.lt_solve_in_place();
        for i in 0..self.n {
            x[i] = self.work[self.perm_inv[i]];
        }
    }

    fn lt_solve_in_place(&mut self) {
        for j in (0..self.n).rev() {
            let a = self.l_colptr[j];
            let mut s = self.work[j];
            for k in a + 1..self.l_colptr[j + 1] {
                s -= self.l_val[k] * self.work[self.l_rowind[k]];
            }
            self.work[j] = s / self.l_val[a];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unordered_factor_keeps_the_order() {
        // tridiagonal(-1, 3, -1): L is lower bidiagonal in the given order
        let n = 5;
        let (mut colptr, mut rowind, mut val) = (vec![0], Vec::new(), Vec::new());
        for j in 0..n {
            rowind.push(j);
            val.push(3.0);
            if j + 1 < n {
                rowind.push(j + 1);
                val.push(-1.0);
            }
            colptr.push(rowind.len());
        }
        let lower = Csc::from_parts(n, n, colptr, rowind, val.clone()).unwrap();
        let mut ch = SparseChol::new_unordered(&lower).unwrap();
        ch.factorize(&val).unwrap();
        let (lp, li, lv) = ch.factor();
        let mut d = 3.0f64.sqrt();
        for j in 0..n {
            assert_eq!(li[lp[j]], j);
            assert!((lv[lp[j]] - d).abs() < 1e-14);
            if j + 1 < n {
                assert_eq!(&li[lp[j] + 1..lp[j + 1]], &[j + 1]);
                let off = -1.0 / d;
                assert!((lv[lp[j] + 1] - off).abs() < 1e-14);
                d = (3.0 - off * off).sqrt();
            }
        }
    }

    #[test]
    fn solves_and_draw_covariance() {
        // A = tridiagonal(-1, 4, -1) plus a corner coupling, n = 7.
        let n = 7;
        let mut dense = vec![0.0; n * n];
        for i in 0..n {
            dense[i * n + i] = 4.0 + i as f64 * 0.1;
            if i + 1 < n {
                dense[i * n + i + 1] = -1.0;
                dense[(i + 1) * n + i] = -1.0;
            }
        }
        dense[6] = 0.5; // (6, 0)
        dense[6 * n] = 0.5;
        let (mut colptr, mut rowind, mut val) = (vec![0], Vec::new(), Vec::new());
        for j in 0..n {
            for i in j..n {
                if dense[j * n + i] != 0.0 {
                    rowind.push(i);
                    val.push(dense[j * n + i]);
                }
            }
            colptr.push(rowind.len());
        }
        let lower = Csc::from_parts(n, n, colptr, rowind, val.clone()).unwrap();
        let mut ch = SparseChol::new(&lower).unwrap();
        ch.factorize(&val).unwrap();
        let b: Vec<f64> = (0..n).map(|i| (i as f64).cos()).collect();
        let mut x = vec![0.0; n];
        ch.solve(&b, &mut x);
        for i in 0..n {
            let ax: f64 = (0..n).map(|j| dense[j * n + i] * x[j]).sum();
            assert!((ax - b[i]).abs() < 1e-13);
        }
        // Columns of Q' L'^{-T}: X X' must be A^{-1}, i.e. A X X' = I.
        let mut xm = vec![0.0; n * n];
        for c in 0..n {
            let mut e = vec![0.0; n];
            e[c] = 1.0;
            ch.lt_solve_perm(&e, &mut xm[c * n..(c + 1) * n]);
        }
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for p in 0..n {
                    let xxt: f64 = (0..n).map(|c| xm[c * n + p] * xm[c * n + j]).sum();
                    s += dense[p * n + i] * xxt;
                }
                assert!((s - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
    }
}
