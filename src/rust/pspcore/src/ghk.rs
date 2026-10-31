//! GHK (Geweke-Hajivassiliou-Keane) log-likelihood for one draw (plan/06
//! §6.8), the port of `.psp_ghk_loglik` and `psp_ghk_cpp`.
//!
//! The latent vector is `(z_1, ..., z_T)` and, for the marginal likelihood
//! of a random-effects model, alpha. With `A_t = R_t S_t` the centred kernel
//! is `sum_t |A_t z_t - R_t alpha|^2 + |B alpha|^2 / sigma^2`, whose
//! precision Q is sparse; given alpha only the z blocks remain. For
//! `x ~ N(0, Q^{-1})` with `Q = L L'`, `L'x ~ N(0, I)`, so the coordinates
//! are drawn one at a time from the last to the first, each truncated to its
//! interval.

use crate::dist::{pnorm, rtnorm};
use crate::rf::{FilterKind, RfContext, RfPar};
use crate::rng::Rng;
use crate::sparse::{Csc, FilterPattern, SparseChol};
use crate::PspError;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::linalg::amd;
use faer::sparse::SymbolicSparseColMat;

/// The GHK estimate for one draw: `ll` and its delta-method Monte Carlo se.
#[derive(Debug, Clone, Copy)]
pub struct Ghk {
    pub ll: f64,
    pub mcse: f64,
}

// log(Phi(b) - Phi(a)) for a < b, accurate in both tails. A one-sided
// interval needs one tail only (the other term is exactly zero).
fn log_interval(a: f64, b: f64) -> f64 {
    if b == f64::INFINITY {
        return pnorm(a, 0.0, 1.0, false, true);
    }
    if a == f64::NEG_INFINITY {
        return pnorm(b, 0.0, 1.0, true, true);
    }
    if a > 0.0 {
        let (la, lb) = (
            pnorm(a, 0.0, 1.0, false, true),
            pnorm(b, 0.0, 1.0, false, true),
        );
        return la + (-(lb - la).exp()).ln_1p();
    }
    if b < 0.0 {
        let (la, lb) = (
            pnorm(a, 0.0, 1.0, true, true),
            pnorm(b, 0.0, 1.0, true, true),
        );
        return lb + (-(la - lb).exp()).ln_1p();
    }
    (pnorm(b, 0.0, 1.0, true, false) - pnorm(a, 0.0, 1.0, true, false)).ln()
}

// The filter of period t on its own pattern: W only if rho != 0, M only if
// lambda != 0 (`R_t S_t`, or `R_t` with rho = 0 for `r_only`).
fn filter(ctx: &RfContext, t: usize, rho: f64, lambda: f64) -> Result<Csc, PspError> {
    let md = ctx.model();
    let n = md.n_units;
    let w = (rho != 0.0).then(|| md.w.at(t));
    let m = (lambda != 0.0).then(|| md.m.at(t));
    let pat = FilterPattern::new(n, w, m)?;
    let mut a = pat.new_matrix();
    let mut colsq = vec![0.0; n];
    pat.update(rho, lambda, &mut a, &mut colsq);
    Ok(a)
}

fn push(trip: &mut Vec<(usize, usize, f64)>, a: &Csc, r0: usize, c0: usize, scale: f64) {
    for j in 0..a.n_cols {
        for (i, v) in a.col(j) {
            trip.push((r0 + i, c0 + j, scale * v));
        }
    }
}

// Fill-reducing order of the pattern of the symmetric matrix `q` restricted
// to `idx` (and, with `extra`, of q_rr + |C||C|' via the free columns):
// positions -> entries of idx.
fn amd_order(q: &Csc, idx: &[usize], extra: Option<&[usize]>) -> Result<Vec<usize>, PspError> {
    let n = idx.len();
    if n <= 1 {
        return Ok(idx.to_vec());
    }
    let mut pos = vec![usize::MAX; q.n_rows];
    for (k, &i) in idx.iter().enumerate() {
        pos[i] = k;
    }
    let mut cols: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (k, &j) in idx.iter().enumerate() {
        cols[k].push(k);
        for (i, _) in q.col(j) {
            if pos[i] != usize::MAX {
                cols[k].push(pos[i]);
            }
        }
    }
    if let Some(free) = extra {
        // |C||C|': rest entries coupled through a common free coordinate
        for &f in free {
            let rows: Vec<usize> = q
                .col(f)
                .filter_map(|(i, _)| (pos[i] != usize::MAX).then_some(pos[i]))
                .collect();
            for &a in &rows {
                cols[a].extend_from_slice(&rows);
            }
        }
    }
    let mut colptr = vec![0];
    let mut rowind = Vec::new();
    for c in &mut cols {
        c.sort_unstable();
        c.dedup();
        rowind.extend_from_slice(c);
        colptr.push(rowind.len());
    }
    let nnz = rowind.len();
    let pat = SymbolicSparseColMat::new_checked(n, n, colptr, None, rowind);
    let (mut fwd, mut inv) = (vec![0usize; n], vec![0usize; n]);
    let mut buf = MemBuffer::try_new(amd::order_scratch::<usize>(n, nnz))
        .map_err(|e| PspError::InvalidInput(format!("AMD: {e:?}")))?;
    amd::order(
        &mut fwd,
        &mut inv,
        pat.as_ref(),
        amd::Control::default(),
        MemStack::new(&mut buf),
    )
    .map_err(|e| PspError::InvalidInput(format!("AMD: {e:?}")))?;
    Ok(fwd.into_iter().map(|k| idx[k]).collect())
}

/// What the likelihood needs besides the model and the draw.
pub struct GhkSpec<'a> {
    /// Outcome codes (binary 0/1, ordered 1..J; `None` = missing).
    pub y: &'a [Option<i32>],
    pub ordered: bool,
    pub reps: usize,
}

/// The GHK log-likelihood of one draw. With `alpha` it is conditional on
/// the random effects; otherwise, for a random-effects model (`par.sig2 >
/// 0`), the random effects are integrated out (marginal). `kappa` is the
/// full cut-point vector of an ordered model. Every draw should get a clone
/// of the same generator (common random numbers).
pub fn ghk_loglik<R: Rng>(
    ctx: &mut RfContext,
    par: &RfPar,
    alpha: Option<&[f64]>,
    kappa: Option<&[f64]>,
    spec: &GhkSpec,
    rng: &mut R,
) -> Result<Ghk, PspError> {
    let (nu, tt) = (ctx.model().n_units, ctx.model().periods);
    let n = nu * tt;
    let marginal = par.sig2 > 0.0 && alpha.is_none();
    let (rho, lambda, gamma) = (par.rho, par.lambda, par.gamma);

    // latent mean: m_t = S_t^{-1}(xb_t + gamma m_{t-1})
    let mut xb = vec![0.0; n];
    crate::dense::x_mul(&ctx.model().x, n, par.beta, &mut xb);
    if let Some(a) = alpha {
        for t in 0..tt {
            for i in 0..nu {
                xb[t * nu + i] += a[i];
            }
        }
    }
    let mut m = xb;
    for t in 0..tt {
        if t > 0 && gamma != 0.0 {
            for i in 0..nu {
                m[t * nu + i] += gamma * m[(t - 1) * nu + i];
            }
        }
        if rho != 0.0 {
            ctx.solve_filter(
                FilterKind::S,
                t,
                rho,
                &mut m[t * nu..(t + 1) * nu],
                1,
                false,
            )?;
        }
    }

    // A: blocks R_t S_t, and -gamma R_t below the diagonal (latent AR)
    let mut ta = Vec::new();
    let mut tr = Vec::new(); // stacked R_t (n x N)
    for t in 0..tt {
        push(&mut ta, &filter(ctx, t, rho, lambda)?, t * nu, t * nu, 1.0);
        let r = filter(ctx, t, 0.0, lambda)?;
        if gamma != 0.0 && t > 0 {
            push(&mut ta, &r, t * nu, (t - 1) * nu, -gamma);
        }
        if marginal {
            push(&mut tr, &r, t * nu, 0, 1.0);
        }
    }
    let a = Csc::from_triplets(n, n, &ta)?;
    let at = a.transpose();
    let qzz = at.matmul(&a)?;
    let dim = if marginal { n + nu } else { n };
    let q = if marginal {
        let rall = Csc::from_triplets(n, nu, &tr)?;
        let c = at.matmul(&rall)?; // -C
        let rt = rall.transpose();
        let mut qaa = rt.matmul(&rall)?;
        let b = filter(ctx, 0, par.phi, 0.0)?;
        let btb = b.transpose().matmul(&b)?;
        let mut tq = Vec::new();
        push(&mut tq, &qzz, 0, 0, 1.0);
        push(&mut tq, &c, 0, n, -1.0);
        push(&mut tq, &c.transpose(), n, 0, -1.0);
        push(&mut tq, &qaa, n, n, 1.0);
        push(&mut tq, &btb, n, n, 1.0 / par.sig2);
        qaa = Csc::from_triplets(dim, dim, &tq)?;
        qaa
    } else {
        qzz
    };

    // bounds; missing outcomes and alpha are unrestricted
    let kap: Vec<f64> = match kappa {
        Some(k) => k.to_vec(),
        None => vec![f64::NEG_INFINITY, 0.0, f64::INFINITY],
    };
    let mut lo = vec![f64::NEG_INFINITY; dim];
    let mut hi = vec![f64::INFINITY; dim];
    for (i, y) in spec.y.iter().enumerate() {
        if let Some(y) = *y {
            let code = if spec.ordered { y } else { y + 1 };
            let c = (code.max(1) as usize).min(kap.len() - 1);
            lo[i] = kap[c - 1];
            hi[i] = kap[c];
        }
    }
    let mut mean = m;
    mean.resize(dim, 0.0);

    // order: unrestricted coordinates first (drawn last, integrated out
    // exactly), each part by a fill-reducing order
    let free: Vec<usize> = (0..dim)
        .filter(|&i| lo[i].is_infinite() && hi[i].is_infinite())
        .collect();
    let rest: Vec<usize> = (0..dim)
        .filter(|&i| !(lo[i].is_infinite() && hi[i].is_infinite()))
        .collect();
    let p: Vec<usize> = if free.is_empty() {
        amd_order(&q, &rest, None)?
    } else {
        let mut p = amd_order(&q, &free, None)?;
        p.extend(amd_order(&q, &rest, Some(&free))?);
        p
    };
    let mut inv = vec![0; dim];
    for (k, &i) in p.iter().enumerate() {
        inv[i] = k;
    }
    // lower triangle of Q[p, p]
    let mut tl = Vec::with_capacity(q.nnz());
    for j in 0..dim {
        for (i, v) in q.col(j) {
            let (pi, pj) = (inv[i], inv[j]);
            if pi >= pj {
                tl.push((pi, pj, v));
            }
        }
    }
    let lower = Csc::from_triplets(dim, dim, &tl)?;
    let mut ch = SparseChol::new_unordered(&lower)?;
    ch.factorize(&lower.val)?;
    let (lp, li, lv) = ch.factor();

    // independent periods without random effects to integrate out
    let blocks_of = |i: usize| {
        if marginal || gamma != 0.0 || i >= n {
            0
        } else {
            i / nu
        }
    };
    let nb = if marginal || gamma != 0.0 { 1 } else { tt };
    let blk: Vec<usize> = p.iter().map(|&i| blocks_of(i)).collect();
    let plv: Vec<f64> = p.iter().map(|&i| lo[i] - mean[i]).collect();
    let phv: Vec<f64> = p.iter().map(|&i| hi[i] - mean[i]).collect();
    let reps = spec.reps;
    let mut logw = vec![0.0; reps * nb];
    let mut x = vec![0.0; dim];
    for r in 0..reps {
        for i in (0..dim).rev() {
            let (mut d, mut s) = (0.0, 0.0);
            for k in lp[i]..lp[i + 1] {
                if li[k] == i {
                    d = lv[k];
                } else {
                    s += lv[k] * x[li[k]];
                }
            }
            let a = if plv[i].is_infinite() {
                plv[i]
            } else {
                d * plv[i] + s
            };
            let b = if phv[i].is_infinite() {
                phv[i]
            } else {
                d * phv[i] + s
            };
            if !(a.is_infinite() && b.is_infinite()) {
                logw[blk[i] * reps + r] += log_interval(a, b);
            }
            x[i] = (rtnorm(rng, 0.0, 1.0, a, b) - s) / d;
        }
    }
    let (mut ll, mut rel2) = (0.0, 0.0);
    for b in 0..nb {
        let lw = &logw[b * reps..(b + 1) * reps];
        let mx = lw.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let w: Vec<f64> = lw.iter().map(|v| (v - mx).exp()).collect();
        let mean_w = w.iter().sum::<f64>() / reps as f64;
        ll += mx + mean_w.ln();
        if reps > 1 {
            let var = w.iter().map(|v| (v - mean_w).powi(2)).sum::<f64>() / (reps - 1) as f64;
            rel2 += var / (mean_w * mean_w) / reps as f64;
        }
    }
    Ok(Ghk {
        ll,
        mcse: rel2.sqrt(),
    })
}
