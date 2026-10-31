//! Normal distribution functions ported from R's nmath, so that values agree
//! with R's own, and the truncated normal sampler (plan/06 §6.5, §6.8).

mod dnorm;
mod pnorm;
mod qnorm;
pub mod truncnorm;

pub use dnorm::dnorm;
pub use pnorm::{pnorm, pnorm_both, Tail};
pub use qnorm::qnorm_std;
pub use truncnorm::rtnorm;
