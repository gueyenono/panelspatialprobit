//! Sparse Cholesky factorisation with a fixed pattern (plan/06 §6.7): the
//! symbolic analysis is done once by faer; each new set of values is
//! factorised numerically into the same buffers.
//!
//! [`SparseChol::new`] (AMD ordering) makes a supernodal factor, which is
//! faster at every size the random-effects block meets (1.8x at N = 1,000,
//! 3.9x at N = 10,000 on a rook grid; plan/06 R6).
//! [`SparseChol::new_unordered`] keeps the factor simplicial, so that L is an
//! ordinary CSC matrix ([`SparseChol::factor`], for GHK).

use super::Csc;
use crate::PspError;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::perm::PermRef;
use faer::sparse::linalg::cholesky::supernodal::SupernodalLltRef;
use faer::sparse::linalg::cholesky::{
    factorize_symbolic_cholesky, CholeskySymbolicParams, SymbolicCholesky, SymbolicCholeskyRaw,
    SymmetricOrdering,
};
use faer::sparse::linalg::SupernodalThreshold;
use faer::sparse::{SparseColMatRef, SymbolicSparseColMat};
use faer::{Conj, MatMut, Par, Side};

/// A symmetric positive-definite matrix given by its lower triangle on a
/// fixed pattern, and its Cholesky factor `Q A Q' = L L'`.
pub struct SparseChol {
    n: usize,
    lower: SymbolicSparseColMat<usize>,
    symbolic: SymbolicCholesky<usize>,
    // L as CSC, for a simplicial factor (empty for a supernodal one)
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
    /// entries only; the values are not used), with an AMD ordering and a
    /// supernodal factor.
    pub fn new(lower: &Csc) -> Result<Self, PspError> {
        Self::build(
            lower,
            SymmetricOrdering::Amd,
            SupernodalThreshold::FORCE_SUPERNODAL,
        )
    }

    /// As [`Self::new`], keeping the given order (`Q = I`) and a simplicial
    /// factor: for callers that order the coordinates themselves and read L
    /// (GHK). The identity is passed as a custom ordering: faer 0.24.4's
    /// `SymmetricOrdering::Identity` skips the conversion of a lower-triangle
    /// input and misreads it.
    pub fn new_unordered(lower: &Csc) -> Result<Self, PspError> {
        let id: Vec<usize> = (0..lower.n_rows).collect();
        let perm = PermRef::new_checked(&id, &id, lower.n_rows);
        Self::build(
            lower,
            SymmetricOrdering::Custom(perm),
            SupernodalThreshold::FORCE_SIMPLICIAL,
        )
    }

    /// The factor L (simplicial CSC; the diagonal is the first entry of each
    /// column): column pointers, row indices and values. Needs a factor from
    /// [`Self::new_unordered`].
    pub fn factor(&self) -> (&[usize], &[usize], &[f64]) {
        assert!(
            !self.l_colptr.is_empty(),
            "sparse Cholesky: factor() needs a simplicial factor"
        );
        (&self.l_colptr, &self.l_rowind, &self.l_val)
    }

    fn build(
        lower: &Csc,
        ordering: SymmetricOrdering<usize>,
        threshold: SupernodalThreshold,
    ) -> Result<Self, PspError> {
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
                supernodal_flop_ratio_threshold: threshold,
                ..Default::default()
            },
        )
        .map_err(|e| fail(&e))?;
        let (l_colptr, l_rowind) = match symbolic.raw() {
            SymbolicCholeskyRaw::Simplicial(s) => {
                let f = s.factor();
                (f.col_ptr().to_vec(), f.row_idx().to_vec())
            }
            SymbolicCholeskyRaw::Supernodal(_) => (Vec::new(), Vec::new()),
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
            symbolic
                .factorize_numeric_llt_scratch::<f64>(Par::Seq, Default::default())
                .or(symbolic.solve_in_place_scratch::<f64>(1, Par::Seq)),
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
        // y = Q b, then L y' = y, then L' y'' = y', then x = Q' y''.
        for i in 0..self.n {
            self.work[i] = b[self.perm_fwd[i]];
        }
        self.tri_solve_in_place(false);
        self.tri_solve_in_place(true);
        for i in 0..self.n {
            x[i] = self.work[self.perm_inv[i]];
        }
    }

    /// `x = Q' L'^{-1} z`: for z standard normal, x ~ N(0, A^{-1}). This is the
    /// C++ engine's `permutationPinv() * matrixU().solve(z)`.
    pub fn lt_solve_perm(&mut self, z: &[f64], x: &mut [f64]) {
        self.work.copy_from_slice(z);
        self.tri_solve_in_place(true);
        for i in 0..self.n {
            x[i] = self.work[self.perm_inv[i]];
        }
    }

    /// `work = L^{-1} work`, or `L'^{-1} work` if `transpose`.
    fn tri_solve_in_place(&mut self, transpose: bool) {
        let n = self.n;
        if let SymbolicCholeskyRaw::Supernodal(s) = self.symbolic.raw() {
            let l = SupernodalLltRef::<usize, f64>::new(s, &self.l_val);
            let rhs = MatMut::from_column_major_slice_mut(&mut self.work, n, 1);
            let stack = MemStack::new(&mut self.buf);
            if transpose {
                l.l_transpose_solve_with_conj(Conj::No, rhs, Par::Seq, stack);
            } else {
                l.l_solve_with_conj(Conj::No, rhs, Par::Seq, stack);
            }
            return;
        }
        let (cp, ri, lv, w) = (&self.l_colptr, &self.l_rowind, &self.l_val, &mut self.work);
        if transpose {
            for j in (0..n).rev() {
                let a = cp[j];
                let mut s = w[j];
                for k in a + 1..cp[j + 1] {
                    s -= lv[k] * w[ri[k]];
                }
                w[j] = s / lv[a];
            }
        } else {
            for j in 0..n {
                let a = cp[j];
                let d = w[j] / lv[a];
                w[j] = d;
                for k in a + 1..cp[j + 1] {
                    w[ri[k]] -= lv[k] * d;
                }
            }
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
        // supernodal (AMD) and simplicial (given order) factors
        for mut ch in [
            SparseChol::new(&lower).unwrap(),
            SparseChol::new_unordered(&lower).unwrap(),
        ] {
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
}
