use super::Csc;
use crate::PspError;

/// The spatial filter `A = (I - lambda M)(I - rho W)`, that is
/// `I - rho W - lambda M + rho lambda M W`, on a fixed sparsity pattern
/// (plan/06 §6.7, plan/02 §2.2). W and M are each optional (a model without
/// rho or without lambda).
///
/// The pattern is the union of the patterns of I, W, M and M W, computed once.
/// Each of the four matrices is stored as a value vector on that pattern, so
/// an update for new (rho, lambda) is a single O(nnz) pass with no allocation,
/// which also gives the column sums of squares that the latent sweep needs.
#[derive(Debug, Clone)]
pub struct FilterPattern {
    pub n: usize,
    pub colptr: Vec<usize>,
    pub rowind: Vec<usize>,
    e_i: Vec<f64>,
    e_w: Option<Vec<f64>>,
    e_m: Option<Vec<f64>>,
    e_mw: Option<Vec<f64>>,
}

impl FilterPattern {
    /// The pattern of `(I - lambda M)(I - rho W)` for N x N matrices, without
    /// the W term (`w = None`) or the M term (`m = None`).
    pub fn new(n: usize, w: Option<&Csc>, m: Option<&Csc>) -> Result<Self, PspError> {
        let ok = |a: Option<&Csc>| a.is_none_or(|a| a.is_square() && a.n_rows == n);
        if !ok(w) || !ok(m) {
            return Err(PspError::InvalidInput(
                "spatial filter: W and M must be square and of the same size".to_string(),
            ));
        }
        let eye = Csc::identity(n);
        let mw = match (m, w) {
            (Some(m), Some(w)) => Some(m.matmul(w)?),
            _ => None,
        };
        let mut parts: Vec<&Csc> = vec![&eye];
        parts.extend(w);
        parts.extend(m);
        parts.extend(mw.as_ref());

        // Union of the patterns, column by column.
        let mut colptr = Vec::with_capacity(n + 1);
        let mut rowind = Vec::new();
        colptr.push(0);
        let mut rows = Vec::new();
        for j in 0..n {
            rows.clear();
            for p in &parts {
                rows.extend(p.col(j).map(|(i, _)| i));
            }
            rows.sort_unstable();
            rows.dedup();
            rowind.extend_from_slice(&rows);
            colptr.push(rowind.len());
        }

        // Each matrix as values on the union pattern.
        let scatter = |a: &Csc| -> Vec<f64> {
            let mut out = vec![0.0; rowind.len()];
            for j in 0..n {
                let mut k = colptr[j];
                for (i, v) in a.col(j) {
                    while rowind[k] != i {
                        k += 1;
                    }
                    out[k] = v;
                }
            }
            out
        };
        let e_i = scatter(&eye);
        let e_w = w.map(scatter);
        let e_m = m.map(scatter);
        let e_mw = mw.as_ref().map(scatter);
        Ok(Self {
            n,
            colptr,
            rowind,
            e_i,
            e_w,
            e_m,
            e_mw,
        })
    }

    pub fn nnz(&self) -> usize {
        self.rowind.len()
    }

    /// A matrix with this pattern (values to be set by [`Self::update`]).
    pub fn new_matrix(&self) -> Csc {
        Csc {
            n_rows: self.n,
            n_cols: self.n,
            colptr: self.colptr.clone(),
            rowind: self.rowind.clone(),
            val: vec![0.0; self.nnz()],
        }
    }

    /// Sets `a` (from [`Self::new_matrix`]) to the filter at (rho, lambda) and
    /// `colsq` to its column sums of squares. With one of W and M, the values
    /// are `1 - rho w` (or `1 - lambda m`), the same arithmetic as Eigen's
    /// `I - rho W` in the C++ engine.
    pub fn update(&self, rho: f64, lambda: f64, a: &mut Csc, colsq: &mut [f64]) {
        debug_assert_eq!(a.val.len(), self.nnz());
        debug_assert_eq!(colsq.len(), self.n);
        let rl = rho * lambda;
        for j in 0..self.n {
            let mut s = 0.0;
            for k in self.colptr[j]..self.colptr[j + 1] {
                let mut v = self.e_i[k];
                if let Some(e_w) = &self.e_w {
                    v -= rho * e_w[k];
                }
                if let Some(e_m) = &self.e_m {
                    v -= lambda * e_m[k];
                }
                if let Some(e_mw) = &self.e_mw {
                    v += rl * e_mw[k];
                }
                a.val[k] = v;
                s += v * v;
            }
            colsq[j] = s;
        }
    }
}
