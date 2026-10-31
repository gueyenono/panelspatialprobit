//! Post-processing from R (plan/06 §6.5, §6.8): log-determinant grids,
//! reduced form, impacts and predictions. Each call works through the draws
//! it is given on one worker thread (the main thread checks for interrupts);
//! R spreads blocks of draws over the user's `future::plan()` (R-D11).

use crate::convert::{self, flag, index, reals};
use crate::run;
use pspcore::ghk::{ghk_loglik, GhkSpec};
use pspcore::lndet::{lndet_eigen, lndet_sparse};
use pspcore::par::CHUNK;
use pspcore::post::{category_probs, cut_points, density_weights, impacts, ImpactVar};
use pspcore::rf::{Cols, Rf, RfContext, RfModel, RfPar};
use pspcore::rng::RMrg32k3a;
use pspcore::simulate::{ppc_draw, response_draw, ResponseSpec};
use pspcore::PspError;
use savvy::{
    savvy, IntegerSexp, ListSexp, NotAvailableValue, OwnedIntegerSexp, OwnedListSexp,
    OwnedRealSexp, RealSexp, Sexp,
};
use std::sync::Mutex;

/// A fit's data for post-processing, from `.psp_rs_post_model()`.
struct PostModel {
    rf: RfModel,
    n_cat: usize,
    ordered: bool,
    probes: Vec<f64>,
    n_probes: usize,
    mix: Vec<f64>,
    n_mix: usize,
    df: Option<f64>,
}

impl PostModel {
    fn from_r(model: &ListSexp) -> savvy::Result<Self> {
        let m = "model";
        let n_units = index(model, m, "N")?;
        let periods = index(model, m, "T")?;
        let rf = RfModel {
            n_units,
            periods,
            x: reals(model, m, "X")?,
            k: index(model, m, "k")?,
            w: convert::weights(&convert::list_field(model, m, "W")?, "model$W")?,
            m: convert::weights(&convert::list_field(model, m, "M")?, "model$M")?,
        };
        let probes = reals(model, m, "probes")?;
        let mix = reals(model, m, "mix")?;
        let (n_probes, n_mix) = (index(model, m, "n_probes")?, index(model, m, "n_mix")?);
        if probes.len() != n_units * n_probes || mix.len() != n_units * n_mix {
            return Err(savvy::Error::new(
                "model: `probes` and `mix` must have N rows",
            ));
        }
        Ok(Self {
            rf,
            n_cat: index(model, m, "J")?,
            ordered: flag(model, m, "ordered")?,
            probes,
            n_probes,
            mix,
            n_mix,
            df: Some(convert::real1(model, m, "df")?).filter(|v| !v.is_nan()),
        })
    }

    fn probes(&self) -> Option<Cols<'_>> {
        (self.n_probes > 0).then_some(Cols {
            val: &self.probes,
            cols: self.n_probes,
        })
    }

    fn mix(&self) -> Option<Cols<'_>> {
        (self.n_mix > 0).then_some(Cols {
            val: &self.mix,
            cols: self.n_mix,
        })
    }
}

/// Posterior draws from `.psp_rs_draws()`: one column per draw.
struct Draws {
    n: usize,
    k: usize,
    beta: Vec<f64>,
    sig2: Vec<f64>,
    rho: Vec<f64>,
    lambda: Vec<f64>,
    phi: Vec<f64>,
    gamma: Vec<f64>,
    kappa: Vec<f64>,
    alpha: Vec<f64>,
}

impl Draws {
    fn from_r(draws: &ListSexp, pm: &PostModel) -> savvy::Result<Self> {
        let d = "draws";
        let beta = reals(draws, d, "beta")?;
        let k = pm.rf.k;
        let n = beta.len() / k.max(1);
        let out = Self {
            n,
            k,
            beta,
            sig2: reals(draws, d, "sig2")?,
            rho: reals(draws, d, "rho")?,
            lambda: reals(draws, d, "lambda")?,
            phi: reals(draws, d, "phi")?,
            gamma: reals(draws, d, "gamma")?,
            kappa: reals(draws, d, "kappa")?,
            alpha: reals(draws, d, "alpha")?,
        };
        let nk = if pm.ordered { pm.n_cat - 2 } else { 0 };
        let ok = out.beta.len() == n * k
            && [&out.sig2, &out.rho, &out.lambda, &out.phi, &out.gamma]
                .iter()
                .all(|v| v.len() == n)
            && out.kappa.len() == n * nk
            && (out.alpha.is_empty() || out.alpha.len() == n * pm.rf.n_units);
        if !ok {
            return Err(savvy::Error::new(
                "draws: fields disagree on the number of draws",
            ));
        }
        Ok(out)
    }

    fn par(&self, d: usize) -> RfPar<'_> {
        RfPar {
            beta: &self.beta[d * self.k..(d + 1) * self.k],
            sig2: self.sig2[d],
            rho: self.rho[d],
            lambda: self.lambda[d],
            phi: self.phi[d],
            gamma: self.gamma[d],
        }
    }

    fn alpha(&self, d: usize) -> Option<&[f64]> {
        if self.alpha.is_empty() {
            return None;
        }
        let nu = self.alpha.len() / self.n;
        Some(&self.alpha[d * nu..(d + 1) * nu])
    }

    fn kappa(&self, d: usize, pm: &PostModel) -> Option<Vec<f64>> {
        pm.ordered.then(|| {
            let nk = pm.n_cat - 2;
            cut_points(&self.kappa[d * nk..(d + 1) * nk])
        })
    }
}

fn reduced_form(
    ctx: &Mutex<RfContext>,
    pm: &PostModel,
    dr: &Draws,
    d: usize,
    need_impacts: bool,
) -> Result<Rf, PspError> {
    let mut ctx = ctx
        .lock()
        .map_err(|_| PspError::InvalidInput("reduced form: worker failed".to_string()))?;
    ctx.reduced_form(&dr.par(d), dr.alpha(d), need_impacts, pm.probes(), pm.mix())
}

/// Sparse log-determinants `log |I - r W|` at `points`.
///
/// @noRd
#[savvy]
fn psp_rs_lndet_sparse(w: ListSexp, points: RealSexp) -> savvy::Result<Sexp> {
    let w = convert::csc(&w, "W")?;
    let v = lndet_sparse(&w, points.as_slice()).map_err(savvy::Error::new)?;
    Ok(OwnedRealSexp::try_from_slice(&v)?.into())
}

/// Log-determinants on `grid` from the eigenvalues `re + i im` of W.
///
/// @noRd
#[savvy]
fn psp_rs_lndet_eigen(re: RealSexp, im: RealSexp, grid: RealSexp) -> savvy::Result<Sexp> {
    if re.len() != im.len() {
        return Err(savvy::Error::new("`re` and `im` must have the same length"));
    }
    let v = lndet_eigen(re.as_slice(), im.as_slice(), grid.as_slice());
    Ok(OwnedRealSexp::try_from_slice(&v)?.into())
}

/// The reduced form of the first draw in `draws`, for tests:
/// `list(m, s, q, d1, d2, t1, t2)` (`s` is n x q, column-major).
///
/// @noRd
#[savvy]
fn psp_rs_reduced_form(
    model: ListSexp,
    draws: ListSexp,
    need_impacts: bool,
) -> savvy::Result<Sexp> {
    let pm = PostModel::from_r(&model)?;
    let dr = Draws::from_r(&draws, &pm)?;
    if dr.n == 0 {
        return Err(savvy::Error::new("draws: no draws"));
    }
    let ctx = Mutex::new(RfContext::new(&pm.rf).map_err(savvy::Error::new)?);
    let rf = reduced_form(&ctx, &pm, &dr, 0, need_impacts).map_err(savvy::Error::new)?;
    let mut out = OwnedListSexp::new(7, true)?;
    out.set_name_and_value(0, "m", OwnedRealSexp::try_from_slice(&rf.m)?)?;
    out.set_name_and_value(1, "s", OwnedRealSexp::try_from_slice(&rf.s)?)?;
    out.set_name_and_value(2, "q", OwnedRealSexp::try_from_scalar(rf.q as f64)?)?;
    for (i, (nm, v)) in [
        ("d1", &rf.d1),
        ("d2", &rf.d2),
        ("t1", &rf.t1),
        ("t2", &rf.t2),
    ]
    .into_iter()
    .enumerate()
    {
        out.set_name_and_value(3 + i, nm, OwnedRealSexp::try_from_slice(v)?)?;
    }
    Ok(out.into())
}

/// Impacts for each draw: a vector of length draws x L, draw fastest, with
/// L = vars x categories x 3 x periods laid out as in
/// `pspcore::post::impacts`.
///
/// @param vars_beta 0-based design columns of the regressors.
/// @param vars_theta 0-based design columns of their spatial lags (-1: none).
/// @param lesage Use the LeSage-Pace scale (no reduced-form sd).
/// @noRd
#[savvy]
fn psp_rs_impacts(
    model: ListSexp,
    draws: ListSexp,
    vars_beta: IntegerSexp,
    vars_theta: IntegerSexp,
    by_period: bool,
    lesage: bool,
) -> savvy::Result<Sexp> {
    let pm = PostModel::from_r(&model)?;
    let dr = Draws::from_r(&draws, &pm)?;
    let (vb, vt) = (vars_beta.as_slice(), vars_theta.as_slice());
    if vb.len() != vt.len() {
        return Err(savvy::Error::new(
            "`vars_beta` and `vars_theta` must have the same length",
        ));
    }
    let vars = vb
        .iter()
        .zip(vt)
        .map(|(&b, &t)| {
            let b = usize::try_from(b).ok().filter(|&b| b < pm.rf.k);
            let t = (t >= 0).then(|| usize::try_from(t).ok().filter(|&t| t < pm.rf.k));
            match (b, t) {
                (Some(beta), None) => Ok(ImpactVar { beta, theta: None }),
                (Some(beta), Some(Some(t))) => Ok(ImpactVar {
                    beta,
                    theta: Some(t),
                }),
                _ => Err(savvy::Error::new("impact variables out of range")),
            }
        })
        .collect::<savvy::Result<Vec<_>>>()?;
    let n_cat = if pm.ordered { pm.n_cat } else { 1 };
    let ctx = Mutex::new(RfContext::new(&pm.rf).map_err(savvy::Error::new)?);
    let task = |d: usize, _: &pspcore::par::TaskCtx| {
        let rf = reduced_form(&ctx, &pm, &dr, d, true)?;
        let kap = dr.kappa(d, &pm);
        let g = if lesage {
            density_weights(&rf.m, &vec![1.0; rf.m.len()], 1, kap.as_deref())
        } else {
            density_weights(&rf.m, &rf.s, rf.q, kap.as_deref())
        };
        Ok(impacts(
            &rf,
            &g,
            n_cat,
            dr.par(d).beta,
            &vars,
            pm.rf.n_units,
            by_period,
        ))
    };
    let res = run::run_tasks(dr.n, 1, task, None)?;
    let len = res.first().map_or(0, |r| r.len());
    let mut out = OwnedRealSexp::new(dr.n * len)?;
    for (d, r) in res.iter().enumerate() {
        for (l, v) in r.iter().enumerate() {
            out[d + dr.n * l] = *v;
        }
    }
    Ok(out.into())
}

/// Sums of the link (reduced-form mean) and of the category probabilities
/// over each consecutive block of `pspcore::par::CHUNK` draws:
/// `list(link = n x blocks, prob = n x categories x blocks)`. R adds the
/// blocks in order, so the totals do not depend on how the draws were split
/// between calls, as long as every call starts at a block boundary.
///
/// @noRd
#[savvy]
fn psp_rs_predict(model: ListSexp, draws: ListSexp) -> savvy::Result<Sexp> {
    let pm = PostModel::from_r(&model)?;
    let dr = Draws::from_r(&draws, &pm)?;
    let ctx = Mutex::new(RfContext::new(&pm.rf).map_err(savvy::Error::new)?);
    let task = |d: usize, _: &pspcore::par::TaskCtx| {
        let rf = reduced_form(&ctx, &pm, &dr, d, false)?;
        let kap = dr.kappa(d, &pm);
        let p = category_probs(&rf.m, &rf.s, rf.q, kap.as_deref());
        Ok((rf.m, p))
    };
    let res = run::run_tasks(dr.n, 1, task, None)?;
    let n = pm.rf.n_units * pm.rf.periods;
    let n_cat = if pm.ordered { pm.n_cat } else { 1 };
    let blocks = dr.n.div_ceil(CHUNK);
    let mut link = OwnedRealSexp::new(n * blocks)?;
    let mut prob = OwnedRealSexp::new(n * n_cat * blocks)?;
    for (d, (m, p)) in res.iter().enumerate() {
        let b = d / CHUNK;
        for (i, v) in m.iter().enumerate() {
            link[b * n + i] += v;
        }
        for (i, v) in p.iter().enumerate() {
            prob[b * n * n_cat + i] += v;
        }
    }
    let mut out = OwnedListSexp::new(2, true)?;
    out.set_name_and_value(0, "link", link)?;
    out.set_name_and_value(1, "prob", prob)?;
    Ok(out.into())
}

// One generator per draw from a 6 x draws matrix of L'Ecuyer-CMRG seeds.
fn draw_rngs(seeds: &IntegerSexp, n: usize) -> savvy::Result<Vec<RMrg32k3a>> {
    let s = seeds.as_slice();
    if s.len() != 6 * n {
        return Err(savvy::Error::new(
            "`seeds` must have 6 rows and one column per draw",
        ));
    }
    s.chunks(6)
        .map(|c| {
            let c: [i32; 6] = c.try_into().map_err(|_| savvy::Error::new("bad seed"))?;
            RMrg32k3a::from_seed(c).map_err(savvy::Error::new)
        })
        .collect()
}

fn with_ctx<'a, T>(
    ctx: &Mutex<RfContext<'a>>,
    f: impl FnOnce(&mut RfContext<'a>) -> Result<T, PspError>,
) -> Result<T, PspError> {
    let mut c = ctx
        .lock()
        .map_err(|_| PspError::InvalidInput("post-processing: worker failed".to_string()))?;
    f(&mut c)
}

/// Replicated outcomes for posterior predictive checks: an integer n x draws
/// matrix (column-major), each draw from its own stream.
///
/// @noRd
#[savvy]
fn psp_rs_ppc(model: ListSexp, draws: ListSexp, seeds: IntegerSexp) -> savvy::Result<Sexp> {
    let pm = PostModel::from_r(&model)?;
    let dr = Draws::from_r(&draws, &pm)?;
    let rngs = draw_rngs(&seeds, dr.n)?;
    let ctx = Mutex::new(RfContext::new(&pm.rf).map_err(savvy::Error::new)?);
    let task = |d: usize, _: &pspcore::par::TaskCtx| {
        let mut rng = rngs[d].clone();
        let kap = dr.kappa(d, &pm);
        with_ctx(&ctx, |c| {
            ppc_draw(c, &dr.par(d), dr.alpha(d), kap.as_deref(), pm.df, &mut rng)
        })
    };
    let res = run::run_tasks(dr.n, 1, task, None)?;
    let n = pm.rf.n_units * pm.rf.periods;
    let mut out = OwnedIntegerSexp::new(n * dr.n)?;
    for (d, y) in res.iter().enumerate() {
        for (i, v) in y.iter().enumerate() {
            out[d * n + i] = *v;
        }
    }
    Ok(out.into())
}

fn opt_cols(v: Vec<i32>) -> Vec<Option<usize>> {
    v.into_iter().map(|j| usize::try_from(j).ok()).collect()
}

/// Response paths: for each draw, the mean differences in outcome
/// probabilities (horizon x category x unit group, as in
/// `pspcore::simulate::response_draw`), as a vector of length draws x L,
/// draw fastest.
///
/// @param spec `list(dynamic, lag, wlag, var, wvar, shock, start, permanent,
///   groups, nsim)`: 0-based design columns (-1: absent), the shock per unit,
///   the 0-based first period and the 0-based unit groups.
/// @noRd
#[savvy]
fn psp_rs_response(
    model: ListSexp,
    draws: ListSexp,
    seeds: IntegerSexp,
    spec: ListSexp,
) -> savvy::Result<Sexp> {
    let pm = PostModel::from_r(&model)?;
    let dr = Draws::from_r(&draws, &pm)?;
    let rngs = draw_rngs(&seeds, dr.n)?;
    let s = "spec";
    let dynamic: Vec<usize> = opt_cols(convert::ints(&spec, s, "dynamic")?)
        .into_iter()
        .flatten()
        .collect();
    let lag = opt_cols(convert::ints(&spec, s, "lag")?);
    let wlag = opt_cols(convert::ints(&spec, s, "wlag")?);
    let groups_r = convert::list_field(&spec, s, "groups")?;
    let mut groups = Vec::new();
    for g in groups_r.values_iter() {
        let g = IntegerSexp::try_from(g)
            .map_err(|_| savvy::Error::new("spec: `groups` must hold integer vectors"))?;
        groups.push(
            opt_cols(g.as_slice().to_vec())
                .into_iter()
                .flatten()
                .collect(),
        );
    }
    let shock = reals(&spec, s, "shock")?;
    let n_cat = if pm.ordered { pm.n_cat } else { 1 };
    let rs = ResponseSpec {
        dynamic: &dynamic,
        lag: &lag,
        wlag: &wlag,
        var: index(&spec, s, "var")?,
        wvar: usize::try_from(convert::int1(&spec, s, "wvar")?).ok(),
        shock: &shock,
        start: index(&spec, s, "start")?,
        permanent: flag(&spec, s, "permanent")?,
        groups: &groups,
        nsim: index(&spec, s, "nsim")?,
        n_cat,
        probes: pm.probes(),
        mix: pm.mix(),
        df: pm.df,
    };
    let lev = if pm.ordered { pm.n_cat - 1 } else { 1 };
    if shock.len() != pm.rf.n_units
        || rs.start >= pm.rf.periods
        || lag.len() != lev
        || wlag.len() != lev
        || rs.var >= pm.rf.k
        || groups
            .iter()
            .any(|g: &Vec<usize>| g.is_empty() || g.iter().any(|&i| i >= pm.rf.n_units))
    {
        return Err(savvy::Error::new(
            "spec: inconsistent response-path settings",
        ));
    }
    let ctx = Mutex::new(RfContext::new(&pm.rf).map_err(savvy::Error::new)?);
    let task = |d: usize, _: &pspcore::par::TaskCtx| {
        let mut rng = rngs[d].clone();
        let kap = dr.kappa(d, &pm);
        with_ctx(&ctx, |c| {
            response_draw(c, &dr.par(d), dr.alpha(d), kap.as_deref(), &rs, &mut rng)
        })
    };
    let res = run::run_tasks(dr.n, 1, task, None)?;
    let len = res.first().map_or(0, |r| r.len());
    let mut out = OwnedRealSexp::new(dr.n * len)?;
    for (d, r) in res.iter().enumerate() {
        for (l, v) in r.iter().enumerate() {
            out[d + dr.n * l] = *v;
        }
    }
    Ok(out.into())
}

/// GHK log-likelihood for each draw: `list(ll, mcse)`. Conditional on the
/// random effects when `draws$alpha` is given, marginal otherwise (pooled
/// models have none). Every draw uses the same stream `seed` (common random
/// numbers).
///
/// @param y Outcome codes, NA for missing.
/// @noRd
#[savvy]
fn psp_rs_ghk(
    model: ListSexp,
    draws: ListSexp,
    y: IntegerSexp,
    reps: i32,
    seed: IntegerSexp,
) -> savvy::Result<Sexp> {
    let pm = PostModel::from_r(&model)?;
    let dr = Draws::from_r(&draws, &pm)?;
    let base = draw_rngs(&seed, 1)?.remove(0);
    let y: Vec<Option<i32>> = y
        .as_slice()
        .iter()
        .map(|&v| (!v.is_na()).then_some(v))
        .collect();
    if y.len() != pm.rf.n_units * pm.rf.periods {
        return Err(savvy::Error::new("`y` must have length N T"));
    }
    let reps = usize::try_from(reps)
        .ok()
        .filter(|&r| r > 0)
        .ok_or_else(|| savvy::Error::new("`reps` must be positive"))?;
    let spec = GhkSpec {
        y: &y,
        ordered: pm.ordered,
        reps,
    };
    let ctx = Mutex::new(RfContext::new(&pm.rf).map_err(savvy::Error::new)?);
    let task = |d: usize, _: &pspcore::par::TaskCtx| {
        let mut rng = base.clone();
        let kap = dr.kappa(d, &pm);
        with_ctx(&ctx, |c| {
            ghk_loglik(c, &dr.par(d), dr.alpha(d), kap.as_deref(), &spec, &mut rng)
        })
    };
    let res = run::run_tasks(dr.n, 1, task, None)?;
    let ll: Vec<f64> = res.iter().map(|g| g.ll).collect();
    let mcse: Vec<f64> = res.iter().map(|g| g.mcse).collect();
    let mut out = OwnedListSexp::new(2, true)?;
    out.set_name_and_value(0, "ll", OwnedRealSexp::try_from_slice(&ll)?)?;
    out.set_name_and_value(1, "mcse", OwnedRealSexp::try_from_slice(&mcse)?)?;
    Ok(out.into())
}

/// Category probabilities for each draw: n x categories x draws
/// (column-major).
///
/// @noRd
#[savvy]
fn psp_rs_probs(model: ListSexp, draws: ListSexp) -> savvy::Result<Sexp> {
    let pm = PostModel::from_r(&model)?;
    let dr = Draws::from_r(&draws, &pm)?;
    let ctx = Mutex::new(RfContext::new(&pm.rf).map_err(savvy::Error::new)?);
    let task = |d: usize, _: &pspcore::par::TaskCtx| {
        let rf = reduced_form(&ctx, &pm, &dr, d, false)?;
        let kap = dr.kappa(d, &pm);
        Ok(category_probs(&rf.m, &rf.s, rf.q, kap.as_deref()))
    };
    let res = run::run_tasks(dr.n, 1, task, None)?;
    let len = res.first().map_or(0, |r| r.len());
    let mut out = OwnedRealSexp::new(dr.n * len)?;
    for (d, r) in res.iter().enumerate() {
        out.as_mut_slice()[d * len..(d + 1) * len].copy_from_slice(r);
    }
    Ok(out.into())
}
