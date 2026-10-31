//! R lists into pspcore structures, with errors that name the field
//! (plan/06 §6.5, R-D10). The R side builds the lists in `R/rust.R`.

use pspcore::sparse::{Csc, Weights};
use savvy::{IntegerSexp, ListSexp, RealSexp, Sexp};

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
