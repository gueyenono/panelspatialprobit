use super::Csc;
use crate::PspError;

/// Spatial weights for a panel: the distinct N x N matrices and, for each
/// period t, the index of the matrix it uses (`period[t]`). Time-invariant
/// weights are one matrix used by every period.
#[derive(Debug, Clone, PartialEq)]
pub struct Weights {
    pub mats: Vec<Csc>,
    pub period: Vec<usize>,
    // the transposes, so that products are row-wise gathers (about 1.5x
    // faster than scattering column by column)
    tmats: Vec<Csc>,
}

impl Weights {
    pub fn new(mats: Vec<Csc>, period: Vec<usize>) -> Result<Self, PspError> {
        let bad = |msg: &str| Err(PspError::InvalidInput(format!("weights: {msg}")));
        let Some(first) = mats.first() else {
            return bad("no matrices");
        };
        let n = first.n_rows;
        if mats.iter().any(|m| !m.is_square() || m.n_rows != n) {
            return bad("matrices must be square and of the same size");
        }
        if period.is_empty() || period.iter().any(|&u| u >= mats.len()) {
            return bad("period map is empty or points past the matrices");
        }
        let tmats = mats.iter().map(Csc::transpose).collect();
        Ok(Self {
            mats,
            period,
            tmats,
        })
    }

    /// The number of units N.
    pub fn n(&self) -> usize {
        self.mats[0].n_rows
    }

    /// The number of periods T.
    pub fn periods(&self) -> usize {
        self.period.len()
    }

    /// The matrix of period `t`.
    pub fn at(&self, t: usize) -> &Csc {
        &self.mats[self.period[t]]
    }

    pub fn varies(&self) -> bool {
        self.mats.len() > 1
    }

    /// `out = W v` for the `u`-th distinct matrix (length N).
    pub fn mul_one(&self, u: usize, v: &[f64], out: &mut [f64]) {
        self.tmats[u].tmul_vec(v, out);
    }

    /// `out_t = W_t v_t` for every period, with `v` stacked time-slow,
    /// unit-fast (length N T).
    pub fn apply(&self, v: &[f64], out: &mut [f64]) {
        let n = self.n();
        for t in 0..self.periods() {
            let r = t * n..(t + 1) * n;
            self.tmats[self.period[t]].tmul_vec(&v[r.clone()], &mut out[r]);
        }
    }
}
