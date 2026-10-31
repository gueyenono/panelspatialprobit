use crate::PspError;

/// A sparse matrix in compressed sparse column form, as R's `dgCMatrix`
/// (`p`, `i`, `x`): row indices are sorted and unique within each column.
#[derive(Debug, Clone, PartialEq)]
pub struct Csc {
    pub n_rows: usize,
    pub n_cols: usize,
    pub colptr: Vec<usize>,
    pub rowind: Vec<usize>,
    pub val: Vec<f64>,
}

impl Csc {
    /// Builds a matrix from its parts, checking the CSC invariants.
    pub fn from_parts(
        n_rows: usize,
        n_cols: usize,
        colptr: Vec<usize>,
        rowind: Vec<usize>,
        val: Vec<f64>,
    ) -> Result<Self, PspError> {
        let bad = |msg: &str| Err(PspError::InvalidInput(format!("sparse matrix: {msg}")));
        if colptr.len() != n_cols + 1 || colptr[0] != 0 {
            return bad("column pointers have the wrong length or do not start at 0");
        }
        let nnz = colptr[n_cols];
        if rowind.len() != nnz || val.len() != nnz {
            return bad("row indices or values do not match the column pointers");
        }
        for j in 0..n_cols {
            let (a, b) = (colptr[j], colptr[j + 1]);
            if b < a {
                return bad("column pointers decrease");
            }
            for k in a..b {
                if rowind[k] >= n_rows {
                    return bad("row index out of range");
                }
                if k > a && rowind[k] <= rowind[k - 1] {
                    return bad("row indices are not sorted and unique within a column");
                }
            }
        }
        if val.iter().any(|v| !v.is_finite()) {
            return bad("non-finite value");
        }
        Ok(Self {
            n_rows,
            n_cols,
            colptr,
            rowind,
            val,
        })
    }

    /// The n x n identity.
    /// A matrix from (row, column, value) triplets; duplicates are summed.
    pub fn from_triplets(
        n_rows: usize,
        n_cols: usize,
        trip: &[(usize, usize, f64)],
    ) -> Result<Self, PspError> {
        let mut t: Vec<(usize, usize, f64)> = trip.to_vec();
        t.sort_unstable_by_key(|&(i, j, _)| (j, i));
        let mut colptr = vec![0usize; n_cols + 1];
        let mut rowind: Vec<usize> = Vec::with_capacity(t.len());
        let mut val: Vec<f64> = Vec::with_capacity(t.len());
        let mut last: Option<(usize, usize)> = None;
        for (i, j, v) in t {
            if i >= n_rows || j >= n_cols {
                return Err(PspError::InvalidInput(
                    "sparse matrix: triplet out of range".to_string(),
                ));
            }
            if last == Some((i, j)) {
                if let Some(x) = val.last_mut() {
                    *x += v;
                }
            } else {
                rowind.push(i);
                val.push(v);
                colptr[j + 1] += 1;
                last = Some((i, j));
            }
        }
        for j in 0..n_cols {
            colptr[j + 1] += colptr[j];
        }
        Self::from_parts(n_rows, n_cols, colptr, rowind, val)
    }

    pub fn identity(n: usize) -> Self {
        Self {
            n_rows: n,
            n_cols: n,
            colptr: (0..=n).collect(),
            rowind: (0..n).collect(),
            val: vec![1.0; n],
        }
    }

    pub fn nnz(&self) -> usize {
        self.val.len()
    }

    pub fn is_square(&self) -> bool {
        self.n_rows == self.n_cols
    }

    /// The (row, value) entries of column `j`.
    pub fn col(&self, j: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let (a, b) = (self.colptr[j], self.colptr[j + 1]);
        self.rowind[a..b]
            .iter()
            .copied()
            .zip(self.val[a..b].iter().copied())
    }

    /// `y = A x`. Accumulates column by column in storage order, as Eigen's
    /// sparse-times-dense product does, so results agree bit for bit.
    pub fn mul_vec(&self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.n_cols);
        debug_assert_eq!(y.len(), self.n_rows);
        y.fill(0.0);
        for (j, &xj) in x.iter().enumerate() {
            for k in self.colptr[j]..self.colptr[j + 1] {
                y[self.rowind[k]] += self.val[k] * xj;
            }
        }
    }

    /// `y = A' x`.
    pub fn tmul_vec(&self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.n_rows);
        debug_assert_eq!(y.len(), self.n_cols);
        for (j, yj) in y.iter_mut().enumerate() {
            let mut s = 0.0;
            for k in self.colptr[j]..self.colptr[j + 1] {
                s += self.val[k] * x[self.rowind[k]];
            }
            *yj = s;
        }
    }

    /// The product `self * other` (Gustavson's algorithm), with sorted row
    /// indices. Structural entries are kept even when they cancel to zero.
    pub fn matmul(&self, other: &Csc) -> Result<Csc, PspError> {
        if self.n_cols != other.n_rows {
            return Err(PspError::InvalidInput(
                "sparse product: non-conformable matrices".to_string(),
            ));
        }
        let n = self.n_rows;
        let mut colptr = Vec::with_capacity(other.n_cols + 1);
        let mut rowind = Vec::new();
        let mut val = Vec::new();
        let mut acc = vec![0.0; n];
        let mut mark = vec![usize::MAX; n];
        let mut rows = Vec::new();
        colptr.push(0);
        for j in 0..other.n_cols {
            rows.clear();
            for (k, b) in other.col(j) {
                for (i, a) in self.col(k) {
                    if mark[i] != j {
                        mark[i] = j;
                        acc[i] = 0.0;
                        rows.push(i);
                    }
                    acc[i] += a * b;
                }
            }
            rows.sort_unstable();
            for &i in &rows {
                rowind.push(i);
                val.push(acc[i]);
            }
            colptr.push(rowind.len());
        }
        Ok(Csc {
            n_rows: n,
            n_cols: other.n_cols,
            colptr,
            rowind,
            val,
        })
    }

    /// The transpose, with sorted row indices.
    pub fn transpose(&self) -> Csc {
        let mut count = vec![0usize; self.n_rows + 1];
        for &i in &self.rowind {
            count[i + 1] += 1;
        }
        for i in 0..self.n_rows {
            count[i + 1] += count[i];
        }
        let colptr = count.clone();
        let mut next = count;
        let mut rowind = vec![0; self.nnz()];
        let mut val = vec![0.0; self.nnz()];
        for j in 0..self.n_cols {
            for (i, v) in self.col(j) {
                let p = next[i];
                rowind[p] = j;
                val[p] = v;
                next[i] += 1;
            }
        }
        Csc {
            n_rows: self.n_cols,
            n_cols: self.n_rows,
            colptr,
            rowind,
            val,
        }
    }

    /// The dense column-major copy, for tests.
    pub fn to_dense(&self) -> Vec<f64> {
        let mut d = vec![0.0; self.n_rows * self.n_cols];
        for j in 0..self.n_cols {
            for (i, v) in self.col(j) {
                d[j * self.n_rows + i] = v;
            }
        }
        d
    }
}
