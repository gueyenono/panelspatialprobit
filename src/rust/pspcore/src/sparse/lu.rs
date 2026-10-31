//! Sparse LU factorisation with a fixed pattern (plan/06 §6.8): the column
//! ordering (COLAMD) is computed once by faer; each new set of values is
//! factorised numerically, with partial pivoting, into the same structure.
//! One factorisation serves both `A x = b` and `A' x = b`.
//!
//! [`SparseLu::new`] lets faer choose a simplicial or a supernodal factor
//! (supernodal is much faster for many right-hand sides on large grids);
//! [`SparseLu::new_simplicial`] keeps the factor simplicial, so that `U` is
//! an ordinary CSC matrix whose diagonal gives the log-determinant.

use super::Csc;
use crate::PspError;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::perm::PermRef;
use faer::sparse::linalg::lu::simplicial::{
    factorize_simplicial_numeric_lu, factorize_simplicial_numeric_lu_scratch,
    solve_in_place_scratch, SimplicialLu,
};
use faer::sparse::linalg::lu::{
    factorize_symbolic_lu, LuRef, LuSymbolicParams, NumericLu, SymbolicLu,
};
use faer::sparse::linalg::SupernodalThreshold;
use faer::sparse::{SparseColMatRef, SymbolicSparseColMat};
use faer::{Conj, MatMut, Par};

// The numeric factor: simplicial (our own, with U accessible) or faer's
// choice.
enum Factor {
    Simplicial {
        col_fwd: Vec<usize>,
        col_inv: Vec<usize>,
        row_fwd: Vec<usize>,
        row_inv: Vec<usize>,
        lu: SimplicialLu<usize, f64>,
    },
    Auto {
        symbolic: SymbolicLu<usize>,
        numeric: NumericLu<usize, f64>,
    },
}

/// The LU factors `P A Q = L U` of a square matrix on a fixed pattern.
pub struct SparseLu {
    n: usize,
    pattern: SymbolicSparseColMat<usize>,
    factor: Factor,
    fact_buf: MemBuffer,
    solve_buf: MemBuffer,
    solve_cols: usize,
    factorized: bool,
}

fn fail(e: &dyn std::fmt::Debug) -> PspError {
    PspError::InvalidInput(format!("sparse LU: {e:?}"))
}

impl SparseLu {
    /// Orders the columns of the pattern of `a` (sorted row indices within
    /// each column; the values are not used); faer chooses the factor type.
    pub fn new(a: &Csc) -> Result<Self, PspError> {
        if !a.is_square() {
            return Err(PspError::InvalidInput(
                "sparse LU: the matrix must be square".to_string(),
            ));
        }
        let n = a.n_rows;
        let pattern =
            SymbolicSparseColMat::new_checked(n, n, a.colptr.clone(), None, a.rowind.clone());
        let symbolic =
            factorize_symbolic_lu(pattern.as_ref(), Default::default()).map_err(|e| fail(&e))?;
        let fact_buf = MemBuffer::try_new(
            symbolic.factorize_numeric_lu_scratch::<f64>(Par::Seq, Default::default()),
        )
        .map_err(|e| fail(&e))?;
        let solve_buf = MemBuffer::try_new(symbolic.solve_in_place_scratch::<f64>(1, Par::Seq))
            .map_err(|e| fail(&e))?;
        Ok(Self {
            n,
            pattern,
            factor: Factor::Auto {
                symbolic,
                numeric: NumericLu::new(),
            },
            fact_buf,
            solve_buf,
            solve_cols: 1,
            factorized: false,
        })
    }

    /// As [`Self::new`], with a simplicial factor (needed by
    /// [`Self::log_abs_det`]).
    pub fn new_simplicial(a: &Csc) -> Result<Self, PspError> {
        if !a.is_square() {
            return Err(PspError::InvalidInput(
                "sparse LU: the matrix must be square".to_string(),
            ));
        }
        let n = a.n_rows;
        let pattern =
            SymbolicSparseColMat::new_checked(n, n, a.colptr.clone(), None, a.rowind.clone());
        let symbolic = factorize_symbolic_lu(
            pattern.as_ref(),
            LuSymbolicParams {
                supernodal_flop_ratio_threshold: SupernodalThreshold::FORCE_SIMPLICIAL,
                ..Default::default()
            },
        )
        .map_err(|e| fail(&e))?;
        let (f, i) = symbolic.col_perm().arrays();
        let (col_fwd, col_inv) = (f.to_vec(), i.to_vec());
        let fact_buf =
            MemBuffer::try_new(factorize_simplicial_numeric_lu_scratch::<usize, f64>(n, n))
                .map_err(|e| fail(&e))?;
        let solve_buf = MemBuffer::try_new(solve_in_place_scratch::<usize, f64>(n, 1, Par::Seq))
            .map_err(|e| fail(&e))?;
        Ok(Self {
            n,
            pattern,
            factor: Factor::Simplicial {
                col_fwd,
                col_inv,
                row_fwd: vec![0; n],
                row_inv: vec![0; n],
                lu: SimplicialLu::new(),
            },
            fact_buf,
            solve_buf,
            solve_cols: 1,
            factorized: false,
        })
    }

    pub fn n(&self) -> usize {
        self.n
    }

    /// Factorises the matrix whose values, on the analysed pattern, are
    /// `values`. Fails if the matrix is singular.
    pub fn factorize(&mut self, values: &[f64]) -> Result<(), PspError> {
        self.factorized = false;
        let a = SparseColMatRef::new(self.pattern.as_ref(), values);
        let singular = |_| PspError::InvalidInput("sparse LU: the matrix is singular".to_string());
        let stack = MemStack::new(&mut self.fact_buf);
        match &mut self.factor {
            Factor::Simplicial {
                col_fwd,
                col_inv,
                row_fwd,
                row_inv,
                lu,
            } => {
                let col_perm = PermRef::new_checked(col_fwd, col_inv, self.n);
                factorize_simplicial_numeric_lu(row_fwd, row_inv, lu, a, col_perm, stack)
                    .map_err(singular)?;
            }
            Factor::Auto { symbolic, numeric } => {
                symbolic
                    .factorize_numeric_lu(numeric, a, Par::Seq, stack, Default::default())
                    .map_err(singular)?;
            }
        }
        self.factorized = true;
        Ok(())
    }

    /// `log |det A|`, from the diagonal of `U` (L has a unit diagonal).
    /// Needs a factor from [`Self::new_simplicial`].
    pub fn log_abs_det(&self) -> Result<f64, PspError> {
        debug_assert!(self.factorized);
        let Factor::Simplicial { lu, .. } = &self.factor else {
            return Err(PspError::InvalidInput(
                "sparse LU: log-determinant needs a simplicial factor".to_string(),
            ));
        };
        let u = lu.u_factor_unsorted();
        let mut s = 0.0;
        for j in 0..self.n {
            for (i, v) in u.row_idx_of_col(j).zip(u.val_of_col(j)) {
                if i == j {
                    s += v.abs().ln();
                }
            }
        }
        Ok(s)
    }

    fn grow_solve_buf(&mut self, cols: usize) -> Result<(), PspError> {
        if cols > self.solve_cols {
            let req = match &self.factor {
                Factor::Simplicial { .. } => {
                    solve_in_place_scratch::<usize, f64>(self.n, cols, Par::Seq)
                }
                Factor::Auto { symbolic, .. } => {
                    symbolic.solve_in_place_scratch::<f64>(cols, Par::Seq)
                }
            };
            self.solve_buf = MemBuffer::try_new(req).map_err(|e| fail(&e))?;
            self.solve_cols = cols;
        }
        Ok(())
    }

    /// Overwrites the n x `cols` column-major matrix `b` with `A^{-1} b`
    /// (`transpose = false`) or `A^{-T} b`.
    pub fn solve(&mut self, b: &mut [f64], cols: usize, transpose: bool) -> Result<(), PspError> {
        debug_assert!(self.factorized);
        debug_assert_eq!(b.len(), self.n * cols);
        if cols == 0 {
            return Ok(());
        }
        self.grow_solve_buf(cols)?;
        let rhs = MatMut::from_column_major_slice_mut(b, self.n, cols);
        let stack = MemStack::new(&mut self.solve_buf);
        match &self.factor {
            Factor::Simplicial {
                col_fwd,
                col_inv,
                row_fwd,
                row_inv,
                lu,
            } => {
                let row_perm = PermRef::new_checked(row_fwd, row_inv, self.n);
                let col_perm = PermRef::new_checked(col_fwd, col_inv, self.n);
                if transpose {
                    lu.solve_transpose_in_place_with_conj(
                        row_perm,
                        col_perm,
                        Conj::No,
                        rhs,
                        Par::Seq,
                        stack,
                    );
                } else {
                    lu.solve_in_place_with_conj(row_perm, col_perm, Conj::No, rhs, Par::Seq, stack);
                }
            }
            Factor::Auto { symbolic, numeric } => {
                // `numeric` is the output of `symbolic.factorize_numeric_lu`
                // on this pattern (`factorized`), as new_unchecked requires
                let lu = LuRef::new_unchecked(symbolic, numeric);
                if transpose {
                    lu.solve_transpose_in_place_with_conj(Conj::No, rhs, Par::Seq, stack);
                } else {
                    lu.solve_in_place_with_conj(Conj::No, rhs, Par::Seq, stack);
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solves_and_log_det_match_dense_algebra() {
        // A nonsymmetric matrix that needs pivoting: zero on the diagonal.
        let n = 6;
        let mut dense = vec![0.0; n * n]; // column-major
        for i in 0..n {
            dense[i * n + i] = if i == 2 { 0.0 } else { 2.0 + i as f64 };
            dense[((i + 1) % n) * n + i] = -1.0 + 0.1 * i as f64;
            dense[i * n + (i + 3) % n] = 0.7;
        }
        let (mut colptr, mut rowind, mut val) = (vec![0], Vec::new(), Vec::new());
        for j in 0..n {
            for i in 0..n {
                if dense[j * n + i] != 0.0 {
                    rowind.push(i);
                    val.push(dense[j * n + i]);
                }
            }
            colptr.push(rowind.len());
        }
        let a = Csc::from_parts(n, n, colptr, rowind, val.clone()).unwrap();
        let mut lu = SparseLu::new_simplicial(&a).unwrap();
        lu.factorize(&val).unwrap();
        let mut auto = SparseLu::new(&a).unwrap();
        auto.factorize(&val).unwrap();
        // log|det| by dense Gaussian elimination with partial pivoting
        let mut m = dense.clone();
        let mut want = 0.0;
        for k in 0..n {
            let p = (k..n)
                .max_by(|&x, &y| m[k * n + x].abs().total_cmp(&m[k * n + y].abs()))
                .unwrap();
            for j in 0..n {
                m.swap(j * n + k, j * n + p);
            }
            want += m[k * n + k].abs().ln();
            for i in k + 1..n {
                let f = m[k * n + i] / m[k * n + k];
                for j in k..n {
                    m[j * n + i] -= f * m[j * n + k];
                }
            }
        }
        assert!((lu.log_abs_det().unwrap() - want).abs() < 1e-13);
        assert!(auto.log_abs_det().is_err());
        let b: Vec<f64> = (0..2 * n).map(|i| (i as f64).cos()).collect();
        for (k, transpose) in [(0, false), (0, true), (1, false), (1, true)] {
            let f = if k == 0 { &mut lu } else { &mut auto };
            let mut x = b.clone();
            f.solve(&mut x, 2, transpose).unwrap();
            for c in 0..2 {
                for i in 0..n {
                    let ax: f64 = (0..n)
                        .map(|j| {
                            let aij = if transpose {
                                dense[i * n + j]
                            } else {
                                dense[j * n + i]
                            };
                            aij * x[c * n + j]
                        })
                        .sum();
                    assert!((ax - b[c * n + i]).abs() < 1e-13);
                }
            }
        }
    }
}
