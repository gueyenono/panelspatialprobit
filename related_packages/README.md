# Related packages

Each subfolder holds the CRAN source tarball and its unpacked source (`<pkg>/<pkg>/R`, `src`, `man`, ...). These are third-party code for reference only, not part of this package. All are GPL-licensed, so code copied from them brings GPL terms with it.

| Package | Version | Language / method | Models | Panel? | Ordered? |
|---|---|---|---|---|---|
| `spatialprobit/` | 1.0.4 | R, Bayesian MCMC (Gibbs + truncated MVN, sparse `Matrix`) | `sarprobit`, `semprobit`, `sarorderedprobit`, `sartobit`; marginal effects | No | Yes (`sarorderedprobit`) |
| `ProbitSpatial/` | 1.1 | R + Rcpp/RcppEigen, approximate ML (Mendell–Elston) | SAR, SEM, SARAR probit (`ProbitSpatialFit`, conditional or full likelihood); Katrina data | No | No |
| `pmlsp/` | 1.0.1 | R, partial ML (Billé & Leorato) | SAR and SARAR probit | No | No |
| `McSpatial/` | 2.0 (archived on CRAN) | R, GMM (Pinkse–Slade / Klier–McMillen) and ML | `gmmprobit`, `spprobit`, `spprobitml` | No | No |
| `spmixW/` | 0.2.2 | R, Bayesian MCMC | SAR/SDM/SEM/SDEM **panel** models with fixed effects and convex combinations of W | Yes | No (continuous outcome only) |

Notes:
- **`sarorderedprobit` is a function in `spatialprobit`**, not a separate package (`spatialprobit/spatialprobit/R/sarorderedprobit.R`).
- `spatialprobit` is the natural base to extend: the Gibbs sampler is in `R/SpatialProbit-MCMC.R`, `R/semprobit.R` and `R/sarorderedprobit.R`. The panel version mostly needs W expanded to I_T ⊗ W_N (or time-specific W_t), individual effects added (with their own Gibbs block), and possibly time dynamics.
- `spmixW` is a useful reference for the panel structure (demeaning, fixed effects, Bayesian rho sampling via log-determinants) to combine with the probit data-augmentation step.

## Non-R implementations (not downloaded)
- **Stata**: `spatbinary` (Spinelli 2022) fits cross-sectional SAR logit and probit. Official `sp*` commands support binary outcomes only through linear probability models. See `../literature/Spinelli_2022_spatbinary_Stata.pdf`.
- **MATLAB**: LeSage's Spatial Econometrics Toolbox (spatial-econometrics.com) is where the `spatialprobit` code was ported from (SAR probit, SAR ordered probit). Whether it has a panel probit routine was not checked.
- **Python**: PySAL `spreg` has a probit with spatial diagnostics but no spatial-lag probit estimator, and no panel version was found.

## Gap
No R package (or other package found in this survey) estimates a **panel** spatial probit or spatial ordered probit. The methods exist in the literature (Kakamu & Wago 2005; Baltagi, Egger & Kesina 2016/2018; Wang & Kockelman 2009 for the ordered/dynamic case; Arbia, Billé & Leorato 2023 for fixed-effects ML), but none has a released implementation.
