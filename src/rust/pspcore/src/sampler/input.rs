use crate::sparse::Weights;
use crate::PspError;

/// A parameter drawn on a log-determinant grid (rho, lambda; phi uses the
/// rho grid): Beta(a, b) prior on (min, max), and log|I - x W| averaged over
/// periods at the grid points.
#[derive(Debug, Clone)]
pub struct GridPrior {
    pub a: f64,
    pub b: f64,
    pub min: f64,
    pub max: f64,
    pub grid: Vec<f64>,
    pub lndet: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SigmaPrior {
    /// Half-Cauchy(0, scale) on sigma_alpha.
    HalfCauchy { scale: f64 },
    /// Inverse gamma(shape, scale) on sigma_alpha^2.
    InvGamma { shape: f64, scale: f64 },
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KappaPrior {
    Flat,
    Normal { mean: f64, sd: f64 },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CutSampler {
    /// Albert-Chib step plus the joint transformation move (plan/02 §2.7).
    Transform,
    /// Albert-Chib step only.
    AlbertChib,
}

/// Everything one chain needs: the model (`psp_sampler_cpp`'s arguments in
/// the C++ engine). Observations are stacked time-slow, unit-fast.
#[derive(Debug, Clone)]
pub struct SamplerInput {
    /// Outcome codes: binary 1/0 (-1 missing); ordered 1..J (0 missing).
    pub y: Vec<i32>,
    /// n x k design, column-major.
    pub x: Vec<f64>,
    pub k: usize,
    pub n_units: usize,
    pub periods: usize,
    pub n_cat: usize,
    pub ordered: bool,
    pub has_rho: bool,
    pub has_lambda: bool,
    pub has_re: bool,
    pub spatial_re: bool,
    pub w: Weights,
    pub m: Weights,
    pub beta_mean: Vec<f64>,
    /// k x k prior precision, column-major.
    pub beta_prec: Vec<f64>,
    pub rho: GridPrior,
    pub lambda: GridPrior,
    pub phi_a: f64,
    pub phi_b: f64,
    /// Student-t degrees of freedom (0: normal errors).
    pub nu: f64,
    /// Column of the latent AR coefficient (`lag_z`), if any.
    pub lag_col: Option<usize>,
    pub sigma_prior: SigmaPrior,
    pub kappa_prior: KappaPrior,
    pub iter: usize,
    pub warmup: usize,
    pub thin: usize,
    pub z_sweeps: usize,
    pub cut_sampler: CutSampler,
    pub scale_move: bool,
    pub beta_start: Vec<f64>,
    pub rho_start: f64,
    pub lambda_start: f64,
    pub sigma2_start: f64,
    /// Starting values of the free cut-points kappa_2, ..., kappa_{J-1}.
    pub kappa_start: Vec<f64>,
    /// Keep the random effects of every `alpha_every`-th saved draw (0: none).
    pub alpha_every: usize,
}

impl SamplerInput {
    pub fn n(&self) -> usize {
        self.n_units * self.periods
    }

    /// Random effects that need the sparse precision block: an error term or
    /// spatial random effects (plan/02 §2.6).
    pub fn sparse_re(&self) -> bool {
        self.has_re && (self.has_lambda || self.spatial_re)
    }

    /// The number of saved parameters: beta, rho, lambda, sigma_alpha, phi
    /// and the free cut-points, in that order.
    pub fn n_par(&self) -> usize {
        self.k
            + self.has_rho as usize
            + self.has_lambda as usize
            + self.has_re as usize
            + self.spatial_re as usize
            + if self.ordered { self.n_cat - 2 } else { 0 }
    }

    pub fn validate(&self) -> Result<(), PspError> {
        let bad = |m: &str| Err(PspError::InvalidInput(format!("sampler input: {m}")));
        let n = self.n();
        if self.y.len() != n || self.x.len() != n * self.k {
            return bad("y and X do not match N T");
        }
        if self.w.n() != self.n_units || self.m.n() != self.n_units {
            return bad("weights do not match N");
        }
        if self.w.periods() != self.periods || self.m.periods() != self.periods {
            return bad("weights do not match T");
        }
        if self.beta_mean.len() != self.k
            || self.beta_prec.len() != self.k * self.k
            || self.beta_start.len() != self.k
        {
            return bad("beta prior or start does not match k");
        }
        if self.ordered && (self.n_cat < 2 || self.kappa_start.len() != self.n_cat - 2) {
            return bad("cut-point starts do not match J");
        }
        if self.lag_col.is_some_and(|c| c >= self.k) {
            return bad("lag column out of range");
        }
        if self.thin == 0 || self.iter / self.thin == 0 {
            return bad("no draws to save");
        }
        for g in [&self.rho, &self.lambda] {
            if g.grid.len() < 2 || g.grid.len() != g.lndet.len() {
                return bad("log-determinant grid too short or mismatched");
            }
        }
        if self.spatial_re && self.w.varies() {
            return bad("spatial random effects need time-invariant W");
        }
        Ok(())
    }
}

/// One chain's results, laid out as the C++ engine's list.
#[derive(Debug, Clone)]
pub struct ChainOutput {
    /// Saved draws, nsave x n_par, column-major.
    pub draws: Vec<f64>,
    pub nsave: usize,
    pub n_par: usize,
    pub z_mean: Vec<f64>,
    pub alpha_mean: Vec<f64>,
    pub alpha_sd: Vec<f64>,
    /// Transformation-move acceptance rates (NaN if never attempted).
    pub cut_accept: Vec<f64>,
    pub cut_step: Vec<f64>,
    /// Kept random-effect draws, n_alpha x N, column-major.
    pub alpha_draws: Vec<f64>,
    pub n_alpha: usize,
    pub v_mean: Vec<f64>,
}
