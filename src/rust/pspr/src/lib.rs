//! R binding of pspcore (plan/06 §6.4): `#[savvy]` functions only.
//!
//! After changing a `#[savvy]` function, run `savvy-cli update . src/rust/pspr`
//! from the package root to regenerate `R/000-wrappers.R`, `src/init.c` and
//! `src/rust/api.h`, then update `src/registration.c` to match.

// Output must go through R's own printing, and errors must reach R as
// messages rather than panics.
#![deny(
    clippy::print_stdout,
    clippy::print_stderr,
    clippy::unwrap_used,
    clippy::expect_used
)]

mod convert;
mod interrupt;
mod post;
mod run;

use pspcore::dist::{dnorm, pnorm, qnorm_std, rtnorm};
use pspcore::rng::{RMrg32k3a, Rng};
use pspcore::sampler::{cut_kernel, run_chain, ChainOutput};
use pspcore::sparse::FilterPattern;
use savvy::{
    savvy, IntegerSexp, ListSexp, NotAvailableValue, OwnedIntegerSexp, OwnedListSexp,
    OwnedRealSexp, RealSexp, Sexp,
};

fn rng_from_seed(seed: IntegerSexp) -> savvy::Result<RMrg32k3a> {
    let seed: [i32; 6] = seed
        .as_slice()
        .try_into()
        .map_err(|_| savvy::Error::new("`seed` must have 6 elements"))?;
    RMrg32k3a::from_seed(seed).map_err(savvy::Error::new)
}

fn count(n: i32) -> savvy::Result<usize> {
    usize::try_from(n).map_err(|_| savvy::Error::new("`n` must be non-negative"))
}

fn real_map(x: RealSexp, f: impl Fn(f64) -> f64) -> savvy::Result<Sexp> {
    let mut out = OwnedRealSexp::new(x.len())?;
    for (i, &v) in x.as_slice().iter().enumerate() {
        out[i] = f(v);
    }
    Ok(out.into())
}

/// Draws from the Rust ports of R's generators, for tests (plan/06 §6.6).
///
/// @param kind `"unif"`, `"norm"`, `"exp"`, `"gamma"` (`par = c(shape, scale)`)
///   or `"t"` (`par = df`).
/// @param n Number of draws.
/// @param par Distribution parameters (see `kind`).
/// @param seed The six L'Ecuyer-CMRG seed integers (`.Random.seed[-1]`).
/// @noRd
#[savvy]
fn psp_rs_rng(kind: &str, n: i32, par: RealSexp, seed: IntegerSexp) -> savvy::Result<Sexp> {
    let n = count(n)?;
    let par = par.as_slice();
    let need = match kind {
        "unif" | "norm" | "exp" => 0,
        "t" => 1,
        "gamma" => 2,
        _ => return Err(savvy::Error::new(format!("unknown `kind`: \"{kind}\""))),
    };
    if par.len() != need {
        return Err(savvy::Error::new(format!(
            "`par` must have {need} element(s) for kind \"{kind}\""
        )));
    }
    let mut rng = rng_from_seed(seed)?;
    let mut out = OwnedRealSexp::new(n)?;
    for i in 0..n {
        out[i] = match kind {
            "unif" => rng.unif(),
            "norm" => rng.norm(),
            "exp" => rng.exp(),
            "t" => rng.t(par[0]),
            _ => rng.gamma(par[0], par[1]),
        };
    }
    Ok(out.into())
}

/// The seed of the next stream or substream (`parallel::nextRNGStream()`,
/// `nextRNGSubStream()`), for tests.
///
/// @param seed The six L'Ecuyer-CMRG seed integers.
/// @param substream `TRUE` for the next substream.
/// @noRd
#[savvy]
fn psp_rs_next_stream(seed: IntegerSexp, substream: bool) -> savvy::Result<Sexp> {
    let rng = rng_from_seed(seed)?;
    let next = if substream {
        rng.next_substream()
    } else {
        rng.next_stream()
    };
    Ok(OwnedIntegerSexp::try_from_slice(next.seed())?.into())
}

/// Truncated normal draws, N(mean, sd^2) on `[lo, hi]`.
///
/// @param n Number of draws.
/// @param mean,sd,lo,hi Parameters and bounds (bounds may be infinite).
/// @param seed The six L'Ecuyer-CMRG seed integers.
/// @noRd
#[savvy]
fn psp_rs_rtnorm(
    n: i32,
    mean: f64,
    sd: f64,
    lo: f64,
    hi: f64,
    seed: IntegerSexp,
) -> savvy::Result<Sexp> {
    let n = count(n)?;
    let mut rng = rng_from_seed(seed)?;
    let mut out = OwnedRealSexp::new(n)?;
    for i in 0..n {
        out[i] = rtnorm(&mut rng, mean, sd, lo, hi);
    }
    Ok(out.into())
}

/// The ported normal distribution functions, for tests: `pnorm(x)`,
/// `dnorm(x)` and `qnorm(x)` for the standard normal.
///
/// @param fun `"pnorm"`, `"dnorm"` or `"qnorm"`.
/// @param x Values (probabilities for `"qnorm"`).
/// @param lower_tail,log Tail and log flags (`pnorm`; `log` also for `dnorm`).
/// @noRd
#[savvy]
fn psp_rs_norm_fun(fun: &str, x: RealSexp, lower_tail: bool, log: bool) -> savvy::Result<Sexp> {
    match fun {
        "pnorm" => real_map(x, |v| pnorm(v, 0.0, 1.0, lower_tail, log)),
        "dnorm" => real_map(x, |v| dnorm(v, 0.0, 1.0, log)),
        "qnorm" => real_map(x, qnorm_std),
        _ => Err(savvy::Error::new(format!("unknown `fun`: \"{fun}\""))),
    }
}

/// The spatial filter of period `t` (1-based), `A_t = (I - lambda M_t)(I - rho
/// W_t)` (or `I - rho W_t` without M), on its merged pattern, with the column
/// sums of squares, for tests (plan/06 step 12).
///
/// @param w,m Weights from `.psp_rs_weights()` (`m` is ignored unless `has_m`).
/// @param has_m Whether the filter has the error term `I - lambda M`.
/// @param rho,lambda Parameters.
/// @param t Period.
/// @return `list(p, i, x, colsq)`.
/// @noRd
#[savvy]
fn psp_rs_filter(
    w: ListSexp,
    m: ListSexp,
    has_m: bool,
    rho: f64,
    lambda: f64,
    t: i32,
) -> savvy::Result<Sexp> {
    let w = convert::weights(&w, "W")?;
    let m = if has_m {
        Some(convert::weights(&m, "M")?)
    } else {
        None
    };
    let t = usize::try_from(t - 1)
        .ok()
        .filter(|&t| t < w.periods())
        .ok_or_else(|| savvy::Error::new("`t` is out of range"))?;
    let pat = FilterPattern::new(w.n(), Some(w.at(t)), m.as_ref().map(|m| m.at(t)))
        .map_err(savvy::Error::new)?;
    let mut a = pat.new_matrix();
    let mut colsq = vec![0.0; pat.n];
    pat.update(rho, lambda, &mut a, &mut colsq);
    let to_i32 = |v: &[usize]| -> savvy::Result<OwnedIntegerSexp> {
        let v: Vec<i32> = v
            .iter()
            .map(|&x| i32::try_from(x))
            .collect::<Result<_, _>>()?;
        OwnedIntegerSexp::try_from_slice(v)
    };
    let mut out = OwnedListSexp::new(4, true)?;
    out.set_name_and_value(0, "p", to_i32(&a.colptr)?)?;
    out.set_name_and_value(1, "i", to_i32(&a.rowind)?)?;
    out.set_name_and_value(2, "x", OwnedRealSexp::try_from_slice(&a.val)?)?;
    out.set_name_and_value(3, "colsq", OwnedRealSexp::try_from_slice(&colsq)?)?;
    Ok(out.into())
}

/// `W_t v_t` for every period, `v` stacked time-slow, unit-fast, for tests.
///
/// @param w Weights from `.psp_rs_weights()`.
/// @param v Vector of length N T.
/// @noRd
#[savvy]
fn psp_rs_weights_apply(w: ListSexp, v: RealSexp) -> savvy::Result<Sexp> {
    let w = convert::weights(&w, "W")?;
    let v = v.as_slice();
    if v.len() != w.n() * w.periods() {
        return Err(savvy::Error::new("`v` must have length N T"));
    }
    let mut out = vec![0.0; v.len()];
    w.apply(v, &mut out);
    Ok(OwnedRealSexp::try_from_slice(&out)?.into())
}

/// A dummy multi-thread job for testing threads, progress and interrupts
/// (plan/06 step 13): task `i` sums `iters` blocks of 1,000 normal draws from
/// the `i`-th stream after `seed`.
///
/// @param tasks Number of tasks.
/// @param threads Worker threads.
/// @param iters Blocks per task.
/// @param seed The six L'Ecuyer-CMRG seed integers.
/// @param progress Print progress.
/// @return The task sums.
/// @noRd
#[savvy]
fn psp_rs_dummy_job(
    tasks: i32,
    threads: i32,
    iters: i32,
    seed: IntegerSexp,
    progress: bool,
) -> savvy::Result<Sexp> {
    let (tasks, threads, iters) = (count(tasks)?, count(threads)?, count(iters)?);
    let base = rng_from_seed(seed)?;
    let task = |i: usize, ctx: &pspcore::par::TaskCtx| {
        let mut rng = base.clone();
        for _ in 0..=i {
            rng = rng.next_stream();
        }
        let mut s = 0.0;
        for k in 0..iters {
            ctx.check()?;
            ctx.set_progress(k + 1);
            for _ in 0..1000 {
                s += rng.norm();
            }
        }
        Ok(s)
    };
    let printer = progress.then(|| run::Printer::new("Task", tasks, iters));
    let out = run::run_tasks(tasks, threads, task, printer)?;
    Ok(OwnedRealSexp::try_from_slice(&out)?.into())
}

fn real_matrix(x: &[f64], nrow: usize, ncol: usize) -> savvy::Result<OwnedRealSexp> {
    let mut out = OwnedRealSexp::try_from_slice(x)?;
    out.set_dim(&[nrow, ncol])?;
    Ok(out)
}

fn chain_list(o: ChainOutput, has_re: bool) -> savvy::Result<OwnedListSexp> {
    let mut out = OwnedListSexp::new(8, true)?;
    let n_units = o.alpha_mean.len();
    out.set_name_and_value(0, "draws", real_matrix(&o.draws, o.nsave, o.n_par)?)?;
    out.set_name_and_value(1, "z_mean", OwnedRealSexp::try_from_slice(&o.z_mean)?)?;
    out.set_name_and_value(
        2,
        "alpha_mean",
        OwnedRealSexp::try_from_slice(&o.alpha_mean)?,
    )?;
    out.set_name_and_value(3, "alpha_sd", OwnedRealSexp::try_from_slice(&o.alpha_sd)?)?;
    let acc: Vec<f64> = o
        .cut_accept
        .iter()
        .map(|&a| if a.is_nan() { f64::na() } else { a })
        .collect();
    out.set_name_and_value(4, "cut_accept", OwnedRealSexp::try_from_slice(&acc)?)?;
    out.set_name_and_value(5, "cut_step", OwnedRealSexp::try_from_slice(&o.cut_step)?)?;
    let cols = if has_re { n_units } else { 0 };
    let ad = if has_re { &o.alpha_draws[..] } else { &[][..] };
    out.set_name_and_value(6, "alpha_draws", real_matrix(ad, o.n_alpha, cols)?)?;
    out.set_name_and_value(7, "v_mean", OwnedRealSexp::try_from_slice(&o.v_mean)?)?;
    Ok(out)
}

/// Runs the MCMC chains: chain c starts from column c of `seeds`, and up to `threads` chains run at once.
/// The draws do not depend on `threads`.
///
/// @param model Data, weights and flags from `.psp_rs_model()`.
/// @param spec Priors, grids, MCMC settings and starting values.
/// @param seeds 6 x chains integer matrix of L'Ecuyer-CMRG seeds.
/// @param threads Chains run at once.
/// @param progress Print progress.
/// @return A list of chains, each with `draws`, `z_mean`, `alpha_mean`,
/// `alpha_sd`, `cut_accept`, `cut_step`, `alpha_draws` and `v_mean`.
/// @noRd
#[savvy]
fn psp_rs_sample(
    model: ListSexp,
    spec: ListSexp,
    seeds: IntegerSexp,
    threads: i32,
    progress: bool,
) -> savvy::Result<Sexp> {
    let inp = convert::sampler_input(&model, &spec)?;
    let seeds = seeds.as_slice();
    if seeds.is_empty() || seeds.len() % 6 != 0 {
        return Err(savvy::Error::new("`seeds` must have 6 rows"));
    }
    let chains = seeds.len() / 6;
    let rngs = seeds
        .chunks(6)
        .map(|s| {
            let s: [i32; 6] = s.try_into().map_err(|_| savvy::Error::new("bad seed"))?;
            RMrg32k3a::from_seed(s).map_err(savvy::Error::new)
        })
        .collect::<savvy::Result<Vec<_>>>()?;
    let task = |c: usize, ctx: &pspcore::par::TaskCtx| {
        let mut rng = rngs[c].clone();
        run_chain(&inp, &mut rng, ctx)
    };
    let printer = progress.then(|| run::Printer::mcmc(chains, inp.warmup, inp.iter));
    let outs = run::run_tasks(chains, count(threads)?, task, printer)?;
    let mut res = OwnedListSexp::new(chains, false)?;
    for (c, o) in outs.into_iter().enumerate() {
        res.set_value(c, chain_list(o, inp.has_re)?)?;
    }
    Ok(res.into())
}

/// The cut-point updates on a toy ordered model, for `test-cutmove.R`.
///
/// @noRd
#[savvy]
#[allow(clippy::too_many_arguments)]
fn psp_rs_cut_kernel(
    y: IntegerSexp,
    mu: RealSexp,
    z: RealSexp,
    kappa_free: RealSexp,
    n_cat: i32,
    kmean: f64,
    ksd: f64,
    step: f64,
    reps: i32,
    ac: bool,
    mv: bool,
    seed: IntegerSexp,
) -> savvy::Result<Sexp> {
    let mut rng = rng_from_seed(seed)?;
    let (kappa, z) = cut_kernel(
        &mut rng,
        y.as_slice(),
        mu.as_slice(),
        z.as_slice(),
        kappa_free.as_slice(),
        count(n_cat)?,
        kmean,
        ksd,
        step,
        count(reps)?,
        ac,
        mv,
    )
    .map_err(savvy::Error::new)?;
    let mut out = OwnedListSexp::new(2, true)?;
    out.set_name_and_value(0, "kappa", OwnedRealSexp::try_from_slice(&kappa)?)?;
    out.set_name_and_value(1, "z", OwnedRealSexp::try_from_slice(&z)?)?;
    Ok(out.into())
}
