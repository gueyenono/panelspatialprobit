//! The random-effects block with a sparse precision (plan/02 §2.6): an error
//! term (sem, sdem, sarar) or spatial random effects. With B = I - phi W,
//!   P = sum_t R_t' D_t R_t + B'B / sigma^2,
//! D_t = I for normal errors. P has a fixed pattern, analysed once; each draw
//! sets its values and refactorises (plan/06 §6.7).

use super::input::SamplerInput;
use crate::sparse::{Csc, FilterPattern, SparseChol};
use crate::PspError;

/// Union of the lower-triangle patterns of `parts` (all N x N).
fn lower_union(n: usize, parts: &[&Csc]) -> Csc {
    let mut colptr = vec![0];
    let mut rowind = Vec::new();
    let mut rows = Vec::new();
    for j in 0..n {
        rows.clear();
        rows.push(j);
        for p in parts {
            rows.extend(p.col(j).map(|(i, _)| i).filter(|&i| i >= j));
        }
        rows.sort_unstable();
        rows.dedup();
        rowind.extend_from_slice(&rows);
        colptr.push(rowind.len());
    }
    let nnz = rowind.len();
    Csc {
        n_rows: n,
        n_cols: n,
        colptr,
        rowind,
        val: vec![0.0; nnz],
    }
}

/// Position of entry (i, j), i >= j, in a lower pattern.
fn pos(lower: &Csc, i: usize, j: usize) -> usize {
    let (a, b) = (lower.colptr[j], lower.colptr[j + 1]);
    a + lower.rowind[a..b]
        .binary_search(&i)
        .unwrap_or_else(|_| unreachable!("entry outside the analysed pattern"))
}

/// The lower-triangle values of a symmetric matrix on the pattern.
fn scatter(lower: &Csc, a: &Csc, scale: f64, out: &mut [f64]) {
    for j in 0..a.n_cols {
        for (i, v) in a.col(j).filter(|&(i, _)| i >= j) {
            out[pos(lower, i, j)] += scale * v;
        }
    }
}

fn add(a: &Csc, b: &Csc) -> Csc {
    // a + b on the union pattern
    let mut colptr = vec![0];
    let (mut rowind, mut val) = (Vec::new(), Vec::new());
    for j in 0..a.n_cols {
        let mut e: Vec<(usize, f64)> = a.col(j).chain(b.col(j)).collect();
        e.sort_by_key(|x| x.0);
        for (i, v) in e {
            if rowind.len() > colptr[j] && rowind.last() == Some(&i) {
                if let Some(last) = val.last_mut() {
                    *last += v;
                }
            } else {
                rowind.push(i);
                val.push(v);
            }
        }
        colptr.push(rowind.len());
    }
    Csc {
        n_rows: a.n_rows,
        n_cols: a.n_cols,
        colptr,
        rowind,
        val,
    }
}

/// For t errors: the contributions d_i R_ia R_ib to P[a, b] (a >= b), one
/// list per distinct M, with the value indices of R_ia and R_ib in R's CSC
/// values.
type Plan = Vec<(usize, usize, usize, usize)>; // (row i, value of R_ia, value of R_ib, position)

pub(super) struct ReBlock {
    e_i: Vec<f64>,
    e_mpmt: Vec<f64>,
    e_mtm: Vec<f64>,
    e_wpwt: Vec<f64>,
    e_wtw: Vec<f64>,
    diag: Vec<usize>,
    r_pats: Vec<FilterPattern>,
    r_mats: Vec<Csc>,
    plans: Vec<Plan>,
    chol: SparseChol,
    vals: Vec<f64>,
    rhs: Vec<f64>,
    mean: Vec<f64>,
    u: Vec<f64>,
    v: Vec<f64>,
    v2: Vec<f64>,
}

impl ReBlock {
    pub(super) fn new(inp: &SamplerInput) -> Result<Self, PspError> {
        let n = inp.n_units;
        let tt = inp.periods as f64;
        let eye = Csc::identity(n);
        let mut parts: Vec<Csc> = Vec::new();
        let (mut mpmt, mut mtm) = (Vec::new(), Vec::new());
        if inp.has_lambda {
            for m in &inp.m.mats {
                let mt = m.transpose();
                mpmt.push(add(m, &mt));
                mtm.push(mt.matmul(m)?);
            }
            parts.extend(mpmt.iter().cloned());
            parts.extend(mtm.iter().cloned());
        }
        let (mut wpwt, mut wtw) = (None, None);
        if inp.spatial_re {
            let w = &inp.w.mats[0];
            let wt = w.transpose();
            wpwt = Some(add(w, &wt));
            wtw = Some(wt.matmul(w)?);
            parts.extend(wpwt.iter().cloned());
            parts.extend(wtw.iter().cloned());
        }
        let refs: Vec<&Csc> = parts.iter().collect();
        let lower = lower_union(n, &refs);
        let nz = lower.nnz();
        let mut e_i = vec![0.0; nz];
        scatter(&lower, &eye, 1.0, &mut e_i);
        let (mut e_mpmt, mut e_mtm) = (vec![0.0; nz], vec![0.0; nz]);
        if inp.has_lambda {
            // summed over periods; one matrix: times T, as in the C++ engine
            if inp.m.mats.len() == 1 {
                scatter(&lower, &mpmt[0], 1.0, &mut e_mpmt);
                scatter(&lower, &mtm[0], 1.0, &mut e_mtm);
                for v in e_mpmt.iter_mut().chain(e_mtm.iter_mut()) {
                    *v *= tt;
                }
            } else {
                for t in 0..inp.periods {
                    let u = inp.m.period[t];
                    scatter(&lower, &mpmt[u], 1.0, &mut e_mpmt);
                    scatter(&lower, &mtm[u], 1.0, &mut e_mtm);
                }
            }
        }
        let (mut e_wpwt, mut e_wtw) = (vec![0.0; nz], vec![0.0; nz]);
        if let (Some(a), Some(b)) = (&wpwt, &wtw) {
            scatter(&lower, a, 1.0, &mut e_wpwt);
            scatter(&lower, b, 1.0, &mut e_wtw);
        }
        let diag = (0..n).map(|j| pos(&lower, j, j)).collect();

        let (mut r_pats, mut r_mats, mut plans) = (Vec::new(), Vec::new(), Vec::new());
        if inp.has_lambda {
            for m in &inp.m.mats {
                let p = FilterPattern::new(n, None, Some(m))?;
                let r = p.new_matrix();
                // rows of R with the value index of each entry
                let mut rows: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
                for j in 0..n {
                    for k in r.colptr[j]..r.colptr[j + 1] {
                        rows[r.rowind[k]].push((j, k));
                    }
                }
                let mut plan = Vec::new();
                for (i, row) in rows.iter().enumerate() {
                    for &(a, ka) in row {
                        for &(b, kb) in row {
                            if a >= b {
                                plan.push((i, ka, kb, pos(&lower, a, b)));
                            }
                        }
                    }
                }
                r_pats.push(p);
                r_mats.push(r);
                plans.push(plan);
            }
        }
        let chol = SparseChol::new(&lower)?;
        Ok(Self {
            vals: vec![0.0; nz],
            e_i,
            e_mpmt,
            e_mtm,
            e_wpwt,
            e_wtw,
            diag,
            r_pats,
            r_mats,
            plans,
            chol,
            rhs: vec![0.0; n],
            mean: vec![0.0; n],
            u: vec![0.0; n],
            v: vec![0.0; n],
            v2: vec![0.0; n],
        })
    }

    /// Draws the random effects given the period residuals u_t = e_t + alpha
    /// (`u`, stacked), the t weights (`wt`, None for normal errors), lambda,
    /// phi and sigma^2. `normals` supplies N standard normal draws.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn draw(
        &mut self,
        inp: &SamplerInput,
        u: &[f64],
        wt: Option<&[f64]>,
        lambda: f64,
        phi: f64,
        sigma2: f64,
        normals: &mut dyn FnMut() -> f64,
        alpha: &mut [f64],
    ) -> Result<(), PspError> {
        let n = inp.n_units;
        let tt = inp.periods;
        for (rp, r) in self.r_pats.iter().zip(&mut self.r_mats) {
            rp.update(0.0, lambda, r, &mut self.v2);
        }
        let s = 1.0 / sigma2;
        match wt {
            None => {
                let c = tt as f64 + s;
                let l2 = lambda * lambda;
                for k in 0..self.vals.len() {
                    let mut v = c * self.e_i[k];
                    if inp.has_lambda {
                        v = v - lambda * self.e_mpmt[k] + l2 * self.e_mtm[k];
                    }
                    if inp.spatial_re {
                        v = v - (phi * s) * self.e_wpwt[k] + (phi * phi * s) * self.e_wtw[k];
                    }
                    self.vals[k] = v;
                }
            }
            Some(wt) => {
                for k in 0..self.vals.len() {
                    let mut v = s * self.e_i[k];
                    if inp.spatial_re {
                        v = v - (phi * s) * self.e_wpwt[k] + (phi * phi * s) * self.e_wtw[k];
                    }
                    self.vals[k] = v;
                }
                for t in 0..tt {
                    let d = &wt[t * n..(t + 1) * n];
                    if inp.has_lambda {
                        let u_ = inp.m.period[t];
                        let rv = &self.r_mats[u_].val;
                        for &(i, ka, kb, p) in &self.plans[u_] {
                            self.vals[p] += d[i] * rv[ka] * rv[kb];
                        }
                    } else {
                        for i in 0..n {
                            self.vals[self.diag[i]] += d[i];
                        }
                    }
                }
            }
        }
        self.chol.factorize(&self.vals)?;

        // rhs = sum_t R_t' D_t R_t u_t; with normal errors and one M, R'R is
        // applied once to sum_t u_t
        self.rhs.fill(0.0);
        let once = wt.is_none() && inp.has_lambda && inp.m.mats.len() == 1;
        if once {
            self.v2.fill(0.0);
            for t in 0..tt {
                for i in 0..n {
                    self.v2[i] += u[t * n + i];
                }
            }
            let r = &self.r_mats[0];
            r.mul_vec(&self.v2, &mut self.v);
            r.tmul_vec(&self.v, &mut self.rhs);
        }
        for t in 0..if once { 0 } else { tt } {
            let ut = &u[t * n..(t + 1) * n];
            if inp.has_lambda {
                let r = &self.r_mats[inp.m.period[t]];
                r.mul_vec(ut, &mut self.v);
                if let Some(wt) = wt {
                    for i in 0..n {
                        self.v[i] *= wt[t * n + i];
                    }
                }
                r.tmul_vec(&self.v, &mut self.v2);
                for i in 0..n {
                    self.rhs[i] += self.v2[i];
                }
            } else {
                for i in 0..n {
                    self.rhs[i] += wt.map_or(1.0, |w| w[t * n + i]) * ut[i];
                }
            }
        }
        self.chol.solve(&self.rhs, &mut self.mean);
        for v in self.v.iter_mut() {
            *v = normals();
        }
        self.chol.lt_solve_perm(&self.v, &mut self.u);
        for i in 0..n {
            alpha[i] = self.mean[i] + self.u[i];
        }
        Ok(())
    }
}
