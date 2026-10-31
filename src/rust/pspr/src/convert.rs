//! R lists into pspcore structures, with errors that name the field
//! (plan/06 §6.5, R-D10). The R side builds the lists in `R/rust.R`.

use pspcore::sampler::{CutSampler, GridPrior, KappaPrior, SamplerInput, SigmaPrior};
use pspcore::sparse::{Csc, Weights};
use savvy::{IntegerSexp, ListSexp, LogicalSexp, RealSexp, Sexp};

fn field(list: &ListSexp, what: &str, name: &str) -> savvy::Result<Sexp> {
    list.get(name)
        .ok_or_else(|| savvy::Error::new(format!("{what}: missing field `{name}`")))
}

pub fn ints(list: &ListSexp, what: &str, name: &str) -> savvy::Result<Vec<i32>> {
    let x = IntegerSexp::try_from(field(list, what, name)?)
        .map_err(|_| savvy::Error::new(format!("{what}: `{name}` must be an integer vector")))?;
    Ok(x.as_slice().to_vec())
}

pub fn reals(list: &ListSexp, what: &str, name: &str) -> savvy::Result<Vec<f64>> {
    let x = RealSexp::try_from(field(list, what, name)?)
        .map_err(|_| savvy::Error::new(format!("{what}: `{name}` must be a double vector")))?;
    Ok(x.as_slice().to_vec())
}

pub fn index(list: &ListSexp, what: &str, name: &str) -> savvy::Result<usize> {
    let v = ints(list, what, name)?;
    match v.as_slice() {
        [x] if *x >= 0 => Ok(*x as usize),
        _ => Err(savvy::Error::new(format!(
            "{what}: `{name}` must be one non-negative integer"
        ))),
    }
}

fn indices(v: Vec<i32>, what: &str, name: &str) -> savvy::Result<Vec<usize>> {
    v.into_iter()
        .map(usize::try_from)
        .collect::<Result<_, _>>()
        .map_err(|_| savvy::Error::new(format!("{what}: `{name}` has negative entries")))
}

/// A matrix from `.psp_csc()`: `list(n_row, n_col, p, i, x)`.
pub fn csc(list: &ListSexp, what: &str) -> savvy::Result<Csc> {
    let n_row = index(list, what, "n_row")?;
    let n_col = index(list, what, "n_col")?;
    let p = indices(ints(list, what, "p")?, what, "p")?;
    let i = indices(ints(list, what, "i")?, what, "i")?;
    let x = reals(list, what, "x")?;
    Csc::from_parts(n_row, n_col, p, i, x).map_err(|e| savvy::Error::new(format!("{what}: {e}")))
}

/// Weights from `.psp_rs_weights()`: `list(mats, period)`.
pub fn weights(list: &ListSexp, what: &str) -> savvy::Result<Weights> {
    let mats = ListSexp::try_from(field(list, what, "mats")?)
        .map_err(|_| savvy::Error::new(format!("{what}: `mats` must be a list")))?;
    let mut out = Vec::with_capacity(mats.len());
    for (u, m) in mats.values_iter().enumerate() {
        let m = ListSexp::try_from(m).map_err(|_| {
            savvy::Error::new(format!("{what}: `mats[[{}]]` must be a list", u + 1))
        })?;
        out.push(csc(&m, &format!("{what}$mats[[{}]]", u + 1))?);
    }
    let period = indices(ints(list, what, "period")?, what, "period")?;
    Weights::new(out, period).map_err(|e| savvy::Error::new(format!("{what}: {e}")))
}

pub fn real1(list: &ListSexp, what: &str, name: &str) -> savvy::Result<f64> {
    match reals(list, what, name)?.as_slice() {
        [x] => Ok(*x),
        _ => Err(savvy::Error::new(format!(
            "{what}: `{name}` must be one number"
        ))),
    }
}

pub fn int1(list: &ListSexp, what: &str, name: &str) -> savvy::Result<i32> {
    match ints(list, what, name)?.as_slice() {
        [x] => Ok(*x),
        _ => Err(savvy::Error::new(format!(
            "{what}: `{name}` must be one integer"
        ))),
    }
}

pub fn flag(list: &ListSexp, what: &str, name: &str) -> savvy::Result<bool> {
    let x = LogicalSexp::try_from(field(list, what, name)?)
        .map_err(|_| savvy::Error::new(format!("{what}: `{name}` must be TRUE or FALSE")))?;
    match x.as_slice_raw() {
        [v] if *v != i32::MIN => Ok(*v != 0),
        _ => Err(savvy::Error::new(format!(
            "{what}: `{name}` must be TRUE or FALSE"
        ))),
    }
}

fn list_field(list: &ListSexp, what: &str, name: &str) -> savvy::Result<ListSexp> {
    ListSexp::try_from(field(list, what, name)?)
        .map_err(|_| savvy::Error::new(format!("{what}: `{name}` must be a list")))
}

fn grid(spec: &ListSexp, name: &str) -> savvy::Result<GridPrior> {
    let what = &format!("spec${name}");
    let g = list_field(spec, "spec", name)?;
    Ok(GridPrior {
        a: real1(&g, what, "a")?,
        b: real1(&g, what, "b")?,
        min: real1(&g, what, "min")?,
        max: real1(&g, what, "max")?,
        grid: reals(&g, what, "grid")?,
        lndet: reals(&g, what, "lndet")?,
    })
}

/// The sampler input from `.psp_rs_model()` (`model`: data, weights and
/// flags) and the sampler settings (`spec`: priors, grids, MCMC settings and
/// starting values), as built in `R/rust.R`.
pub fn sampler_input(model: &ListSexp, spec: &ListSexp) -> savvy::Result<SamplerInput> {
    let (m, s) = ("model", "spec");
    let n_units = index(model, m, "N")?;
    let periods = index(model, m, "T")?;
    let x = reals(model, m, "X")?;
    let k = index(model, m, "k")?;
    let lag = int1(model, m, "lag_col")?;
    let sigma = list_field(spec, s, "sigma")?;
    let sigma_prior = match int1(&sigma, "spec$sigma", "family")? {
        0 => SigmaPrior::HalfCauchy {
            scale: real1(&sigma, "spec$sigma", "p1")?,
        },
        _ => SigmaPrior::InvGamma {
            shape: real1(&sigma, "spec$sigma", "p1")?,
            scale: real1(&sigma, "spec$sigma", "p2")?,
        },
    };
    let kappa_prior = if flag(spec, s, "kappa_flat")? {
        KappaPrior::Flat
    } else {
        KappaPrior::Normal {
            mean: real1(spec, s, "kappa_mean")?,
            sd: real1(spec, s, "kappa_sd")?,
        }
    };
    let inp = SamplerInput {
        y: ints(model, m, "y")?,
        x,
        k,
        n_units,
        periods,
        n_cat: index(model, m, "J")?,
        ordered: flag(model, m, "ordered")?,
        has_rho: flag(model, m, "has_rho")?,
        has_lambda: flag(model, m, "has_lambda")?,
        has_re: flag(model, m, "has_re")?,
        spatial_re: flag(model, m, "spatial_re")?,
        w: weights(&list_field(model, m, "W")?, "model$W")?,
        m: weights(&list_field(model, m, "M")?, "model$M")?,
        beta_mean: reals(spec, s, "beta_mean")?,
        beta_prec: reals(spec, s, "beta_prec")?,
        rho: grid(spec, "rho")?,
        lambda: grid(spec, "lambda")?,
        phi_a: real1(spec, s, "phi_a")?,
        phi_b: real1(spec, s, "phi_b")?,
        nu: real1(spec, s, "nu")?,
        lag_col: usize::try_from(lag).ok(),
        sigma_prior,
        kappa_prior,
        iter: index(spec, s, "iter")?,
        warmup: index(spec, s, "warmup")?,
        thin: index(spec, s, "thin")?,
        z_sweeps: index(spec, s, "z_sweeps")?,
        cut_sampler: if int1(spec, s, "cut_sampler")? == 0 {
            CutSampler::Transform
        } else {
            CutSampler::AlbertChib
        },
        scale_move: flag(spec, s, "scale_move")?,
        beta_start: reals(spec, s, "beta_start")?,
        rho_start: real1(spec, s, "rho_start")?,
        lambda_start: real1(spec, s, "lambda_start")?,
        sigma2_start: real1(spec, s, "sigma2_start")?,
        kappa_start: reals(spec, s, "kappa_start")?,
        alpha_every: index(spec, s, "alpha_every")?,
    };
    inp.validate()
        .map_err(|e| savvy::Error::new(e.to_string()))?;
    Ok(inp)
}
