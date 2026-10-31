//! R binding of pspcore (plan/06 §6.4): `#[savvy]` functions only.
//!
//! After changing a `#[savvy]` function, run `savvy-cli update .` from the
//! package root to regenerate `R/000-wrappers.R`, `src/init.c` and
//! `src/rust/api.h`, then add the new routine to `src/registration.c`.

// Output must go through R's own printing, and errors must reach R as
// messages rather than panics.
#![deny(
    clippy::print_stdout,
    clippy::print_stderr,
    clippy::unwrap_used,
    clippy::expect_used
)]

use pspcore::rng::RMrg32k3a;
use savvy::{savvy, IntegerSexp, OwnedRealSexp, Sexp};

fn rng_from_seed(seed: IntegerSexp) -> savvy::Result<RMrg32k3a> {
    let seed: [i32; 6] = seed
        .as_slice()
        .try_into()
        .map_err(|_| savvy::Error::new("`seed` must have 6 elements"))?;
    RMrg32k3a::from_seed(seed).map_err(savvy::Error::new)
}

/// Draws from the Rust ports of R's generators, for tests (plan/06 §6.6).
///
/// @param kind Distribution: `"unif"`.
/// @param n Number of draws.
/// @param seed The six L'Ecuyer-CMRG seed integers (`.Random.seed[-1]`).
/// @noRd
#[savvy]
fn psp_rs_rng(kind: &str, n: i32, seed: IntegerSexp) -> savvy::Result<Sexp> {
    let n = usize::try_from(n).map_err(|_| savvy::Error::new("`n` must be non-negative"))?;
    let mut rng = rng_from_seed(seed)?;
    let mut out = OwnedRealSexp::new(n)?;
    match kind {
        "unif" => {
            for i in 0..n {
                out[i] = rng.unif();
            }
        }
        _ => return Err(savvy::Error::new(format!("unknown `kind`: \"{kind}\""))),
    }
    Ok(out.into())
}
