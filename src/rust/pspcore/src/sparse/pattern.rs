use super::Csc;
use crate::PspError;

/// The spatial filter `A = (I - lambda M)(I - rho W) = I - rho W - lambda M
/// + rho lambda M W` on a fixed sparsity pattern (plan/06 §6.7, plan/02 §2.2).
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
    e_w: Vec<f64>,
    e_m: Option<(Vec<f64>, Vec<f64>)>, // (M, M W)
}

impl FilterPattern {
    /// The pattern of `I - rho W` (`m = None`) or of `(I - lambda M)(I - rho W)`.
    pub fn new(w: &Csc, m: Option<&Csc>) -> Result<Self, PspError> {
        let n = w.n_rows;
        if !w.is_square() || m.is_some_and(|m| !m.is_square() || m.n_rows != n) {
            return Err(PspError::InvalidInput(
                "spatial filter: W and M must be square and of the same size".to_string(),
            ));
        }
        let eye = Csc::identity(n);
        let mw = match m {
            Some(m) => Some(m.matmul(w)?),
            None => None,
        };
        let mut parts: Vec<&Csc> = vec![&eye, w];
        if let (Some(m), Some(mw)) = (m, mw.as_ref()) {
            parts.push(m);
            parts.push(mw);
        }

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
        let e_w = scatter(w);
        let e_m = match (m, mw.as_ref()) {
            (Some(m), Some(mw)) => Some((scatter(m), scatter(mw))),
            _ => None,
        };
        Ok(Self {
            n,
            colptr,
            rowind,
            e_i,
            e_w,
            e_m,
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
    /// `colsq` to its column sums of squares. Without M, the values are
    /// computed as `1 - rho w`, the same arithmetic as Eigen's `I - rho W`.
    pub fn update(&self, rho: f64, lambda: f64, a: &mut Csc, colsq: &mut [f64]) {
        debug_assert_eq!(a.val.len(), self.nnz());
        debug_assert_eq!(colsq.len(), self.n);
        let rl = rho * lambda;
        for j in 0..self.n {
            let mut s = 0.0;
            for k in self.colptr[j]..self.colptr[j + 1] {
                let mut v = self.e_i[k] - rho * self.e_w[k];
                if let Some((e_m, e_mw)) = &self.e_m {
                    v = v - lambda * e_m[k] + rl * e_mw[k];
                }
                a.val[k] = v;
                s += v * v;
            }
            colsq[j] = s;
        }
    }
}
