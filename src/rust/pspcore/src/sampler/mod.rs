//! The MCMC sampler (plan/06 §6.7): one chain per task, chains in parallel
//! through [`crate::par`].

mod chain;
mod griddy;
mod input;
mod re;

pub use chain::{cut_kernel, run_chain};
pub use griddy::{draw_from_grid, draw_grid_param, lndet_interp};
pub use input::{ChainOutput, CutSampler, GridPrior, KappaPrior, SamplerInput, SigmaPrior};
