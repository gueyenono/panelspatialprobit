# Literature: panel spatial (ordered) probit

## Downloaded (open access / working-paper versions)

| File | Citation | Why it matters |
|---|---|---|
| `Kakamu_Wago_2005_bayesian_spatial_panel_probit_MODSIM.pdf` | Kakamu, K. & Wago, H. (2005). Bayesian spatial panel probit model with an application to business cycle in Japan. MODSIM 2005. | Early Bayesian (MCMC) panel spatial probit; 47 prefectures. Closest early template for the target model. |
| `Baltagi_Egger_Kesina_2016_bayesian_spatial_bivariate_panel_probit_WP187.pdf` | Baltagi, B. H., Egger, P. H. & Kesina, M. (2016). Bayesian spatial bivariate panel probit estimation. *Advances in Econometrics* 37, 119–144. doi:10.1108/S0731-905320160000037011 (Syracuse CPR WP 187 version) | Bayesian panel spatial probit with random effects, extended to a bivariate system. |
| `Wang_Kockelman_2008_dynamic_spatial_ordered_probit_TRB.pdf` | Wang, X. & Kockelman, K. M. (2009). Bayesian inference for ordered response data with a dynamic spatial-ordered probit model. *Journal of Regional Science* 49(5), 877–913. doi:10.1111/j.1467-9787.2009.00622.x (TRB 2008 version) | **Key paper for the ordered case**: spatial ordered probit on panel data (spatial errors, AR(1) over time, random effects), Gibbs sampling with data augmentation. |
| `GaoPing_thesis_bayesian_spatial_econometric_models_multivariate_ordered_probit.pdf` | Gao, P. Bayesian analysis on spatial econometric models. PhD thesis, Hokkaido University. | Multivariate and dynamic spatial ordered probit (MSOP/DMSOP) by MCMC, applied to China Family Panel Studies. |
| `Elhorst_Heijnen_Samarina_Jacobs_2013_transitions_spatial_probit_SOM.pdf` | Elhorst, J. P., Heijnen, P., Samarina, A. & Jacobs, J. P. A. M. (2017). Transitions at different moments in time: A spatial probit approach. *Journal of Applied Econometrics* 32(2), 422–439. doi:10.1002/jae.2505 (2013 SOM report version) | Spatial probit on a country × year panel (adoption of inflation targeting, 58 countries, 1985–2008). Replication data are on the ZBW Journal Data Archive. |
| `Wilhelm_deMatos_2013_RJournal_spatialprobit.pdf` | Wilhelm, S. & de Matos, M. G. (2013). Estimating spatial probit models in R. *The R Journal* 5(1), 130–143. doi:10.32614/RJ-2013-013 | Explains the `spatialprobit` package: sparse matrices and a Gibbs sampler with truncated multivariate normal draws. This is the code base to extend. |
| `Spinelli_2022_spatbinary_Stata.pdf` | Spinelli, D. (2022). Fitting spatial autoregressive logit and probit models using Stata: The spatbinary command. *Stata Journal* 22(2), 293–318. doi:10.1177/1536867X221106373 | Documents the Stata implementation (cross-sectional only). |

## To obtain (paywalled; citations checked against Crossref)

Cross-sectional foundations:
- LeSage, J. P. & Pace, R. K. (2009). *Introduction to Spatial Econometrics*. CRC Press. doi:10.1201/9781420064254 (Ch. 10 covers spatial probit and ordered probit; `spatialprobit` is ported from its MATLAB code).
- LeSage, J. P. (2000). Bayesian estimation of limited dependent variable spatial autoregressive models. *Geographical Analysis* 32(1). (original Gibbs sampler; DOI not checked)
- LeSage, J. P., Pace, R. K., Lam, N., Campanella, R. & Liu, X. (2011). New Orleans business recovery in the aftermath of Hurricane Katrina. *JRSS-A*. doi:10.1111/j.1467-985x.2011.00712.x (spatial ordered probit application; the source of the Katrina dataset)
- Wang, H., Iglesias, E. M. & Wooldridge, J. M. (2013). Partial maximum likelihood estimation of spatial probit models. *Journal of Econometrics* 172(1), 77–89. doi:10.1016/j.jeconom.2012.08.005
- Martinetti, D. & Geniaux, G. (2017). Approximate likelihood estimation of spatial probit models. *Regional Science and Urban Economics*. doi:10.1016/j.regsciurbeco.2017.02.002 (the method behind `ProbitSpatial`)
- Pace, R. K. & LeSage, J. P. (2016). Fast simulated maximum likelihood estimation of the spatial probit model capable of handling large samples. *Advances in Econometrics* 37, 3–34. doi:10.1108/s0731-905320160000037008
- Billé, A. G. & Leorato, S. (2019). Partial ML estimation for spatial autoregressive nonlinear probit models with autoregressive disturbances. *Econometric Reviews*. doi:10.1080/07474938.2019.1682314 (the method behind `pmlsp`)

Panel / dynamic:
- Baltagi, B. H., Egger, P. H. & Kesina, M. (2018). Generalized spatial autocorrelation in a panel-probit model with an application to exporting in China. *Empirical Economics* 55(1), 193–211. doi:10.1007/s00181-017-1409-0 (SAR probit with individual random effects, MCMC)
- Arbia, G., Billé, A. G. & Leorato, S. (2023). Concentrated partial ML estimation for dynamic spatial panel data probit models with fixed effects. SSRN 4514415. doi:10.2139/ssrn.4514415 (a frequentist alternative with fixed effects)
- Wang, X. & Kockelman, K. M. (2009). Application of the dynamic spatial ordered probit model: Patterns of land development change in Austin, Texas. *Papers in Regional Science*. doi:10.1111/j.1435-5957.2009.00249.x
- Kakamu, K. & Wago, H. (2008). Small-sample properties of panel spatial autoregressive models: Bayesian vs ML. *Spatial Economic Analysis*. doi:10.1080/17421770802353725 (continuous outcome, but provides the Bayesian panel SAR building block)
- Zhou (2018). Essays on estimation for nonlinear spatial models. PhD dissertation, LSU. https://repository.lsu.edu/gradschool_dissertations/4538 (spatial panel probit: GHK simulated ML, and random effects with a Gibbs sampler)
