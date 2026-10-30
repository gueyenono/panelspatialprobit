# 5. Roadmap

The work is split into phases, each with concrete tasks and an **exit criterion**. Sizes are relative effort (S < M < L), not calendar time. A phase starts only when the previous phase's exit criterion is met, except that Phase 4 and Phase 5 can be swapped.

## Phase 0: Project setup (S)

- [~] `git init`; create the GitHub repository (decision D5: local working copy outside Google Drive). *Status: repo initialised; it needs `git config --global --add safe.directory <path>` (NTFS mount, root-owned files). No remote yet.*
- [~] `usethis::create_package()` in place. *Status: written by hand; the author in `Authors@R` is still a placeholder.* Set DESCRIPTION:
  - Title: "Bayesian Panel Spatial Probit and Ordered Probit Models";
  - `Authors@R` with ORCID;
  - License GPL (≥ 3) (D2).
- [~] `use_roxygen_md()`, `use_testthat(3)`, `use_rcpp_eigen()` (or `use_rcpp()` + LinkingTo), `use_news_md()`, `use_readme_rmd()`, `use_pkgdown()`. *Status: roxygen markdown, testthat 3, RcppEigen, NEWS and vignettes are set up by hand; README is plain Markdown; no pkgdown yet.*
- [x] Add a "hello sparse" C++ function: maps a `dgCMatrix` and returns `W %*% x`, to prove the toolchain works. *Status: superseded by the sampler itself, which maps `dgCMatrix` weights.*
- [x] `.Rbuildignore` / `.gitignore` as in `03` §3.5.
- [~] GitHub Actions: `R-CMD-check` (3 OS × 3 R versions), `test-coverage`, `pkgdown`. *Status: all three workflows are written under `.github/workflows/` (untested until there is a GitHub remote).*
- [x] Re-run `/init` to generate CLAUDE.md. It should cover build/test commands, the stacking convention, the $Az - b$ abstraction and the SBC gate.
- [ ] Obtain the paywalled papers listed in "Reading to obtain" below.

**Exit:** `R CMD check --as-cran` passes with 0 errors, warnings and notes on CI.

## Phase 1: Core engine, binary SAR probit with pooled and random effects (L)

- [x] 1.1 Data pipeline (`model_frame.R`, `weights.R`). Balanced panels, time-invariant $W$, strict ID matching, internal (time, unit) sort with the inverse permutation kept.
- [x] 1.2 Log-determinant grids: eigenvalue and sparse-LU methods, cubic-spline interpolation, tests against `determinant()` and `spatialreg`.
- [x] 1.3 **Pure-R reference sampler** (dense, small $N$). This is the oracle for everything that follows.
- [x] 1.4 C++ engine:
  - `truncnorm`;
  - `latent_sweep` (residual update);
  - `griddy` (two-stage, collapsed $\rho$);
  - $\beta$ block;
  - $\alpha$ block;
  - half-Cauchy $\sigma_\alpha$.
- [x] 1.5 `psp_fit()`: chains, seeds (L'Ecuyer), progress bar, warm-up and thinning, the fit object. `psprobit()` front-end.
- [x] 1.6 Methods: `print`, `summary`, `coef`, `vcov`, `confint`, `nobs`, `as_draws`, basic `plot`.
- [x] 1.7 `sim_psprobit()` for SAR or none × pooled or random, binary.
- [x] 1.8 Unit tests (`04` §4.1). Geweke test. **SBC** for {none, sar} × {pooled, random}. *Status: unit tests and SBC pass (`validation/sbc/results/`). Instead of a full Geweke test there is a stationarity test of the cut-point kernels (`test-cutmove.R`); together with SBC across all variants and the pure-R oracle, it covers what the Geweke test would catch.*
- [x] 1.9 SBC experiment on `spatialprobit::sarprobit` (`04` §4.2), to confirm or retract `02` §2.10. *Status: done; the result is recorded in `02` §2.10.*
- [x] 1.10 First benchmark: $N = 1{,}000$, $T = 10$. Record ESS per second for $\rho$, $\sigma_\alpha$ and the intercept.

**Exit:**
- SBC passes for the four variants;
- the C++ engine and the R oracle agree distributionally on fixed datasets;
- the Katrina $T = 1$ cross-check is documented;
- the benchmark is recorded.

## Phase 2: Interpretation and richer specifications (M)

- [x] 2.1 `impacts()` / `psp_impacts()`:
  - population-averaged and unit-specific;
  - averaged over periods and per period;
  - exact or stochastic diagonals;
  - finite-difference tests;
  - `scale = "lesage_pace"` compatibility option.
- [x] 2.2 CRE (Mundlak and Chamberlain) and time effects. Rank checks.
- [x] 2.3 Durbin terms: `sdm`, `slx`, and impacts with $\theta$.
- [x] 2.4 `fitted()` and `predict(type = "link" | "prob" | "class")` in-sample. Hit rate and AUC (`psp_classification()`).
- [x] 2.5 `tidy()` and `glance()`.
- [x] 2.6 SBC for `sdm`, `slx` and `cre` (Mundlak and Chamberlain).
- [x] 2.7 Vignette "Getting started with panel spatial probit".

**Exit:** impacts tests pass, Durbin and CRE SBC pass, and the vignette builds. Tag **v0.1.0** as a GitHub release (not CRAN).

## Phase 3: Ordered probit (M)

- [x] 3.1 `psoprobit()` front-end: response coding, cut-point parameterisation, the `polr`-style display in `summary`.
- [x] 3.2 Albert–Chib cut-point Gibbs (baseline).
- [x] 3.3 Transformation move on $(\kappa_j, z)$ (`02` §2.7), with the detailed-balance unit test. *Status: the unit test checks stationarity (invariance of the joint distribution of $(\kappa, z)$ on a toy model, `test-cutmove.R`) rather than enumerating detailed balance.*
- [x] 3.4 SBC for both samplers ($J = 3, 5$; pooled and random). Compare ESS per second and **choose the default**. *Status: both pass SBC (J = 3, 4). At N = 400, T = 5 the transformation move gives 4x (J = 3) and 18-78x (J = 5, interior cut-points) the ESS per second of Albert-Chib alone, at about 15% more time per iteration, so `transform` stays the default (`validation/bench/cutpoints_2026-10-27.txt`).*
- [x] 3.5 Category-specific impacts and `predict(type = "prob")` returning an $n \times J$ matrix.
- [~] 3.6 `sim_psprobit(..., cutpoints =)`; Monte Carlo subset; cross-check against `sarorderedprobit` and `MCMCoprobit`. *Status: `cutpoints =` exists; agreement with `MCMCpack::MCMCoprobit` within 0.03 posterior sd (`validation/crosscheck/`). Against `spatialprobit::sarorderedprobit` (1.0.4; 20 x 20 grid, T = 1, J = 4, `validation/crosscheck/spatialprobit_ordered.R`, 2026-10-29): all posterior means within 0.07 sd and sds within 2%; rho differs by 0.004 (beyond MC error, consistent with its concentrated rho kernel); its Albert-Chib cut-points reach bulk ESS about 400 from 4 x 50,000 draws in 4.2 min, against 6,000-7,000 here in 2 s. Monte Carlo (`validation/montecarlo/`, 200 replications, N = 100, T = 10, J = 4): biases of at most 0.017, 90% coverage 0.88-0.92 for coefficients, rho, sigma_alpha, cut-points and the impact of x1.*
- [x] 3.7 Vignette "Ordered outcomes".

**Exit:** SBC passes for both samplers and the default is justified by ESS-per-second evidence. **This covers the package's stated purpose. First CRAN submission candidate (v0.5.0)**, after a docs pass (see Phase 7.4 checklist).

## Phase 4: Error-structure variants (M)

- [x] 4.1 SEM: $\lambda$ griddy step, and the $\alpha$ block via sparse Cholesky.
- [x] 4.2 SARAR ($A = RS$) and SDEM. Identification warning when $M = W$.
- [x] 4.3 Student-$t$ errors (scale mixture). *Status: `errors = psp_t(df)` with fixed df (robit; df is barely identified with discrete outcomes). Precision weights enter every quadratic form of the sampler. Impacts and predictions average the conditional-normal formulas over `mixture_draws` simulated sets of mixing variances, fixed per fit (above `impacts_exact_max` the neighbours' error variances enter at their mean). Checked against the reference sampler (SAR binary; SARAR ordered), the exact t probabilities without spatial terms, and finite differences. SBC: `sar_t_re_bin` and `sem_t_re_ord3` pass. `sarar_t_re_ord3` showed lambda p = 0.024 / 0.030 in two runs, but with chains twice as long lambda passes (p = 0.775) and rho is at p = 0.027, like normal-error SARAR (p = 0.019): slow rho/lambda mixing of SARAR with M = W, not a t-specific error.*
- [x] 4.4 (Optional) `effects = "spatial_random"`. *Status: $\alpha = \phi W\alpha + a$; the $\alpha$ block uses the sparse Cholesky of $\sum_t R_t'R_t + B'B/\sigma^2$ and $\phi$ a griddy step on the $W$ grid. Needs time-invariant $W$. Checked against the reference sampler (SAR and SEM) and finite-difference impacts; SBC variants `sar_sre_bin`, `sem_sre_ord3`.*
- [x] 4.5 SBC for each, and impacts adapted (SEM has direct effects only). *Status: SBC passes for sem (pooled, RE, ordered), sdem and sarar; impacts checked against finite differences.*

**Exit:** SBC passes for every new variant.

## Phase 5: Dynamics (M–L)

- [x] 5.1 `psp_dynamic()`: $y_{t-1}$ and $Wy_{t-1}$ columns, first-period handling, Wooldridge initial conditions (default).
- [x] 5.2 Short-run impacts. Simulated dynamic response paths via posterior prediction. *Status: `psp_response_path()`: temporary or permanent shocks to some or all units; outcome histories simulated with common random numbers, probabilities computed exactly given each history; horizon 0 equals the unit-specific impact (tested).*
- [x] 5.3 (Stretch) Latent AR(1) $\gamma z_{t-1}$: block-bidiagonal $A$, a griddy or Gaussian draw for $\gamma$, treatment of $z_0$. Replicate the Wang & Kockelman simulation. *Status: `psp_dynamic(latent = TRUE)`. $\gamma$ is the coefficient of a design column `lag_z` refreshed from $z$, so it is drawn with $\beta$ (Gaussian) and rho stays collapsed over it; the latent sweep and the cut move include the coupling to period $t+1$; the scale move leaves $\gamma$ unscaled. The first period has no lag (it is the initial condition and shares $\alpha$). Time-invariant weights only; not combined with lagged outcomes. Reduced form: exact recursions for small $N$, simulated probe recursions otherwise; checked against a brute-force stacked solve, the dense reference sampler (incl. SARAR + ordered + t errors), `mvtnorm` for the likelihood, and finite differences. SBC: `sem_lat_re_ord3` passes; `sar_lat_re_bin` fails at the default chain length (lag_z p = 0.001, concentrated in the 45% of replications with max R-hat >= 1.02, U-shaped ranks) and passes with chains three times as long (min p = 0.090): gamma mixes slowly against rho given the latent values, so latent dynamics need longer chains (the convergence warning flags it). The Wang & Kockelman simulation was replicated (`validation/wang_kockelman/`, 2026-10-29): their 12 cells (lambda in {0.1, 0.5, 0.9} x rho in {0.1, 0.6, 0.7, 0.9}, 50 data sets each), adapted to one effect per unit (N = 300, queen grid), unit error variances and an intercept. 600 fits, none failed; overall 90% coverage 0.89, largest average bias 0.12 (sigma^2), against drifts in the paper of up to +0.83 in sigma^2. With rho = 0.9 about two thirds of fits keep R-hat >= 1.05 even with 15,000 draws per chain (the intercept / thresholds), so strongly autocorrelated effects need longer chains. Summarised in `vignette("validation")`.*
- [x] 5.4 SBC and a Monte Carlo study (state dependence versus spatial dependence versus heterogeneity). *Status: SBC passes (`sar_dyn_bin`). Monte Carlo (`validation/montecarlo/results/mc_2026-10-27.txt`, 200 replications, N = 100, T = 10, gamma = 0.8, sigma_alpha = 0.8): the Wooldridge RE model recovers gamma (bias 0.01, 90% coverage 0.93); a pooled dynamic model confounds heterogeneity with state dependence (gamma bias +0.31, coverage 0.11; slopes biased towards zero); a static RE model leaves the other parameters roughly unbiased. Static binary SAR at N = 100 and 400 also gives 90% coverage of 0.90-0.92.*

**Exit:** SBC passes, and the dynamic vignette section is written.

## Phase 6: Scale, robustness and model assessment (L)

- [x] 6.1 Time-varying $W_t$ (list input): per-period log-determinant grids with caching, plus Chebyshev and MC log-determinants for large $N$. *Status: `W` and `M` accept a list with one matrix per period (named by period). The C++ engine keeps one $A_t$ per period only when the weights vary; the log-determinant grid is the average over periods, each distinct matrix computed once. Durbin terms, `W_lag_y`, impacts, predictions and PPC are per period. Checked against the reference sampler and finite-difference impacts; SBC variants `sar_tv_re_bin`, `sarar_tv_re_ord3`. Chebyshev / MC log-determinants are not added: the sparse-LU spline grid already handles large $N$.*
- [x] 6.2 Unbalanced panels with missing $y$: latent imputation with unrestricted bounds. A clear error remains when $X$ is missing. *Rows must still be present (with `NA` outcomes); dynamic models need complete outcomes.*
- [x] 6.3 Parallel chains through `parallel`. Verify reproducibility across `cores`. *Status: done early, in Phase 1.*
- [~] 6.4 Mixing upgrades justified by ESS per second (`02` §2.6): a joint $(\beta, \alpha)$ block, centred/non-centred switching, PX-DA. *Status: a generalised-Gibbs scale move (the PX-DA idea) was added in Phase 1; it roughly doubled ESS for the slopes and sigma_alpha. Three further moves were implemented, verified against the reference sampler, measured and reverted for lack of any ESS gain: a latent-AR time-profile move, a scale move including sigma_alpha, and centred/non-centred interweaving for sigma_alpha. The remaining slowness (rare outcomes, near-separation, strong latent dynamics) is the data-augmentation slowness of Johndrow et al. (2019); see `02` §2.6, "Mixing experiments". It is documented in `?psprobit` (section Convergence) and in the vignette.*
- [x] 6.5 Posterior predictive checks (`ppc.R`). Approximate `logLik` (Mendell–Elston or GHK, labelled experimental) for model comparison. *Status: `psp_ppc()`; `psp_loglik()`, `logLik()` and `psp_dic()` use GHK on the sparse precision of the latent vector (with $\alpha$ appended and ordered so that it is integrated out exactly). Exact in the non-spatial case, checked against quadrature (RE) and `mvtnorm` (SAR). Independent periods are simulated separately, so pooled and conditional likelihoods stay precise for large panels; the marginal RE likelihood is one $NT$-dimensional problem and is noisy beyond a few hundred observations, which a warning reports.*
- [x] 6.6 `predict(newdata = )` for the same units and periods with new $X$, i.e. counterfactual or policy scenarios.
- [x] 6.7 Large benchmark: $N = 10{,}000$, $T = 20$. *Status: 39 ms per iteration, so 10,000 iterations take about 6.5 minutes, against the 30-minute target (`validation/bench/large_2026-10-27.txt`). Rerun 2026-10-30 on the current version: 32.9 ms per iteration (about 5.5 minutes), impacts with 50 draws 19 s (`validation/bench/large_2026-10-30.txt`).*

**Exit:** performance targets met (`02` §2.13) or re-baselined with a written rationale. PPC vignette section written.

## Phase 7: Documentation, paper and release (M)

- [x] 7.1 Vignettes: *Status: all five are written: getting started, ordered outcomes, interpreting impacts (`impacts`), dynamics (`dynamics`) and validation (`validation`). The validation vignette reads two small CSV summaries in `vignettes/`, regenerated from `validation/` results by `validation/vignette_data.R`; rerun it after new SBC or Monte Carlo runs.*
  - getting started;
  - ordered outcomes;
  - interpreting impacts;
  - dynamics;
  - "How we validated the samplers" (an SBC summary with figures generated from `validation/`).
- [ ] 7.2 Empirical example data in `data/`, built from `data-raw/` scripts with licences noted (`04` §4.6). *Open: needs the choice of a real dataset whose licence allows redistribution.*
- [~] 7.3 pkgdown site and `inst/CITATION`. *Status: `_pkgdown.yml` (reference index grouped by task) and `inst/CITATION` are written; the site URL and the citation author wait for the repository and the author details.*
- [~] 7.4 CRAN readiness:
  - examples run in under 5 s each (`\donttest{}` for longer ones);
  - total check time under 10 min;
  - no non-CRAN Suggests;
  - `urlchecker` and `spelling` pass;
  - `rhub` / win-builder checks;
  - reverse dependencies: none yet.

  *Status: the slowest example takes 1.5 s; `R CMD check --as-cran` takes a few minutes and gives only the new-submission NOTE; every Suggests is on CRAN; `urlchecker` and `spelling` pass (`Language: en-GB`, `inst/WORDLIST`, `tests/spelling.R`). Not done: the rhub / win-builder runs, which email the maintainer.*
- [ ] 7.5 Paper for the R Journal or JSS: models, the sampler design (the $Az - b$ engine and the cut-point move), SBC evidence, Monte Carlo results, an empirical application. *Status: draft abstract, outline and Sections 1–6, 8 and 9 (all but the application) in `paper/paper.md` (2026-10-29); venue and application (Section 7) still to do.*

**Exit:** accepted on CRAN; paper submitted.

## Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| Slow mixing of $\rho$, $\alpha$, $\sigma_\alpha$ or $\kappa$ at large $n$ | High | Collapsed $\rho$ draw; transformation cut-point move; ESS-per-second benchmarks drive the Phase 6 upgrades; multiple chains with $\hat R$ and ESS warnings. |
| Weak identification of $\rho$ with binary data | Medium | Prior-sensitivity guidance; report the effective prior-to-posterior shift; Monte Carlo documents the sample sizes needed. |
| Subtle sampler bugs | Medium | Pure-R oracle; Geweke test; **SBC gate per variant**. |
| Impacts too costly for large $N$ | Medium | Post-hoc computation on thinned draws; stochastic diagonals; per-period caching of LU factors. |
| The transformation cut-point move is harder than expected | Medium | Albert–Chib baseline ships first; the move stays experimental until SBC passes. |
| Data and $W$ mis-alignment by users | High | Strict ID matching; informative errors; a vignette section on building $W$. |
| CRAN check time or size | Low | Small tests; validation kept outside the package; heavy examples in `\donttest{}`. |
| Scope creep (frequentist estimators, multivariate models, FE) | Medium | §1.6 out-of-scope list; v2 backlog in GitHub issues. |

## Reading to obtain before or while working on each phase

| Reference | Needed for |
|---|---|
| LeSage & Pace (2009), *Introduction to Spatial Econometrics*, ch. 5 and 10 | Phase 1 (griddy $\rho$, log-determinants), Phase 3 (§10.2, ordered) |
| Baltagi, Egger & Kesina (2018), *Empirical Economics* 55(1) | Phase 1 and 4: their exact "generalized" spatial RE specification; Monte Carlo design |
| Lacombe & LeSage (2018), *Annals of Regional Science* (SSRN 2013, doi:10.2139/ssrn.2314127) | Phase 2: interpretation and impacts in SAR probit |
| LeSage, Pace, Lam, Campanella & Liu (2011), *JRSS-A* | Phase 2 and 3: impacts and the ordered application |
| Martinetti & Geniaux (2017), *RSUE* / Pace & LeSage (2016) | Phase 6: approximate likelihood |
| Wang & Kockelman (2009), *J. Regional Science* (published version) | Phase 5: latent dynamics; simulation design |
| Wooldridge (2005), *J. Applied Econometrics* 20(1) | Phase 5: initial conditions |
| van Dyk & Park (2008), *JASA*, doi:10.1198/016214508000000409 | Phase 1: collapsed-sampler ordering |
| Liu & Sabatti (2000), *Biometrika* 87(2), doi:10.1093/biomet/87.2.353; Cowles (1996), *Statistics and Computing* | Phase 3: cut-point moves |
| Modrák et al. (2025), *Bayesian Analysis*, doi:10.1214/23-BA1404; Talts et al. (2018) | Phase 1+: SBC |
| Arbia, Billé & Leorato (2023), SSRN 4514415 | Positioning versus the frequentist fixed-effects alternative (paper introduction) |
