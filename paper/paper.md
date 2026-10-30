---
title: "panelspatialprobit: Bayesian Panel Spatial Probit and Ordered Probit Models in R"
author: "TODO (author details, see CLAUDE.md item 2)"
status: "Draft: abstract and full text of Sections 1-6, 8 and 9 (all but the application, Section 7)"
target: "R Journal (rjtools Rmd) or JSS (jss.cls LaTeX); venue not yet chosen"
---

<!--
Drafting notes
- Written in Markdown with LaTeX math so it converts to either venue's format.
- Citations are author-year in the text; the reference list at the end is
  partial and every entry needs checking against the published version
  (several papers are in the "Reading to obtain" list in plan/05_roadmap.md).
- Numbers come from validation/ and vignettes/validation-*.csv as of
  2026-10-29. Regenerate them (validation/vignette_data.R) before submission.
- [TODO] marks things to settle or check.
-->

# Abstract

Binary and ordered outcomes observed for the same spatial units over several
periods are common in economics, geography and transport research, but no
released software estimates spatial probit models with a panel dimension.
The R package **panelspatialprobit** fits Bayesian panel spatial probit and
ordered probit models. They can combine seven spatial specifications;
random, correlated random or spatially autocorrelated unit effects; state
dependence through lagged outcomes with Wooldridge initial conditions or
through latent AR(1) dynamics; Student-$t$ errors; time-varying weights; and
missing outcomes. Every model reduces to one sparse latent Gaussian kernel,
which a single C++ Gibbs sampler updates without solving linear systems.
Ordered models add a joint transformation move on the cut-points and the
latent values. The package reports impacts on probabilities, split into
direct and spillover components, and response paths for dynamic models.
Every supported variant passes simulation-based calibration, and a Monte
Carlo study and the replication of a published simulation confirm the
coverage of its intervals. The same calibration test shows that an existing
cross-sectional sampler understates spatial dependence. [TODO: one sentence
on the empirical application, once Section 7 exists.]

# Outline

1. **Introduction** (drafted below)

2. **Models** (drafted below)
   - 2.1 Latent model and observation rules
   - 2.2 Spatial specifications
   - 2.3 Unobserved heterogeneity
   - 2.4 Dynamics
   - 2.5 Error distribution, weights and missing outcomes
   - 2.6 Identification and priors

3. **Posterior sampling** (drafted below)
   - 3.1 One kernel for all models
   - 3.2 Latent block
   - 3.3 Spatial parameters and coefficients
   - 3.4 Random effects
   - 3.5 Cut-points
   - 3.6 A scale move
   - 3.7 Latent dynamics and $t$ errors
   - 3.8 Computational cost

4. **Interpretation** (drafted below)
   - 4.1 Probabilities
   - 4.2 Impacts
   - 4.3 Dynamic responses
   - 4.4 Model checking and comparison

5. **Validation** (drafted below)
   - 5.1 Components and a reference sampler
   - 5.2 Simulation-based calibration
   - 5.3 The same test applied to `spatialprobit`
   - 5.4 Monte Carlo study
   - 5.5 Replication of the Wang and Kockelman simulation
   - 5.6 Cross-checks against other packages
   - 5.7 Performance

6. **Using the package** (drafted below)
   - 6.1 Data and weights; 6.2 Fitting; 6.3 Output and convergence;
     6.4 Impacts; 6.5 Prediction and counterfactuals; 6.6 Dynamic models;
     6.7 Ordered outcomes; 6.8 Model checking and comparison;
     6.9 Practical advice

7. **Empirical application** [needs the example data, CLAUDE.md item 4]
   - Binary: inflation-targeting adoption (Elhorst et al. 2017), 58
     countries, 1985–2008. Dynamic SAR/SDM with Wooldridge initial
     conditions; state dependence versus spatial diffusion; short-run
     impacts and response paths.
   - Ordered (optional): US Drought Monitor categories.
   - Comparison with the misspecified alternatives (pooled, non-spatial,
     static) via PPC persistence and DIC.

8. **Limitations and extensions** (drafted below)
   - 8.1 Mixing
   - 8.2 Approximations
   - 8.3 Restrictions of the implemented models
   - 8.4 Out of scope
   - 8.5 Limits of the validation
   - 8.6 Extensions

9. **Conclusion** (drafted below)

Computational details (R version, package versions, hardware, seeds), code
availability and replication materials.

---

# 1. Introduction

Many outcomes that economists, geographers and transport researchers study
are binary or ordered, observed for the same units over several periods, and
spatially dependent. Examples include:

- whether a country has adopted inflation targeting (Elhorst et al. 2017);
- whether a firm exports (Baltagi, Egger & Kesina 2018);
- whether a business has reopened after a disaster (LeSage et al. 2011);
- the category of land development of a parcel (Wang & Kockelman 2009).

Two kinds of dependence are present at once. Across space, a unit's outcome
depends on what happens to its neighbours: policies diffuse, firms learn from
nearby exporters, and land use spreads along a frontier. Over time, outcomes
persist, for two reasons with different implications. Some units are
persistently more likely to have the outcome (unobserved heterogeneity),
and the outcome itself may raise the probability of having it again (state
dependence).

Models for these data combine three established ingredients:

- the spatial probit and ordered probit models (LeSage 2000; LeSage & Pace
  2009, chap. 10), estimated by Gibbs sampling with data augmentation
  (Albert & Chib 1993);
- random or correlated random effects for unobserved heterogeneity (Mundlak
  1978; Chamberlain 1982), with Wooldridge's (2005) solution to the initial
  conditions problem in dynamic models;
- the spatial lag, spatial error and Durbin specifications of spatial
  econometrics.

Bayesian panel versions have been developed for spatial probit models with
random effects (Kakamu & Wago 2005; Baltagi, Egger & Kesina 2016, 2018)
and for a dynamic spatial ordered probit (Wang & Kockelman 2009). None of
them, to our knowledge, comes with a released software implementation.

The available software covers cross-sections:

- **`spatialprobit`** (Wilhelm & de Matos 2013) estimates SAR and SEM
  probit and SAR ordered probit models by MCMC, ported from LeSage's MATLAB
  code. It has no panel structure, random effects or multiple chains, and it
  offers no convergence diagnostics.
- **`ProbitSpatial`** (Martinetti & Geniaux 2017) estimates SAR, SEM and
  SARAR probit models by approximate maximum likelihood.
- **`pmlsp`** implements partial maximum likelihood (Billé & Leorato 2020).
- **`spldv`** (Piras & Sarrias 2023) implements GMM estimators.
- **`McSpatial`**, now archived on CRAN, provides GMM and ML estimators
  (Klier & McMillen 2008).
- **Stata:** `spatbinary` (Spinelli 2022) fits cross-sectional SAR logit
  and probit models.
- **Continuous outcomes:** `spmixW` estimates spatial panel models, but
  only for continuous outcomes.

A frequentist estimator for dynamic spatial panel probit models with fixed
effects exists (Arbia, Billé & Leorato 2023), but it needs long panels.
Applied researchers with a short panel of binary or ordered spatial data
therefore either pool the periods and ignore heterogeneity, or drop the
spatial dimension. Both choices bias the estimates: Section 5.4 shows that a
pooled dynamic model overstates state dependence by 40%.

**`panelspatialprobit`** fills this gap. It estimates Bayesian panel spatial
probit and ordered probit models with:

- seven spatial specifications (SAR, SDM, SEM, SARAR, SDEM, SLX and a
  non-spatial baseline);
- pooled, random, correlated random (Mundlak or Chamberlain) and spatially
  autocorrelated unit effects;
- state dependence through lagged outcomes, with Wooldridge initial
  conditions, or through latent AR(1) dynamics;
- Student-$t$ errors, time-varying weights and missing outcomes.

The paper makes four contributions.

1. **A single sampler for the whole model family.** Given its spatial
   parameters, every model reduces to a latent Gaussian kernel
   $\exp(-\tfrac12\|Az - b\|^2)$ with a sparse operator $A$. One C++ engine
   serves all variants:
   - an exact latent update that costs $O(\operatorname{nnz}(A))$ per sweep
     and solves no linear systems;
   - a draw of $\rho$ from its exact conditional with $\beta$ integrated
     out, ordered so that it respects the partially collapsed Gibbs rule
     (van Dyk & Park 2008);
   - for ordered outcomes, a joint transformation move on the cut-points
     and the latent values (Liu & Sabatti 2000).

   The sampler handles $N = 10{,}000$ units over $T = 20$ periods in
   minutes.
2. **Interpretable effects.** Impacts on probabilities account for the
   reduced-form variance that the spatial filter and the random effects
   induce, which existing software ignores. They are reported as
   population-averaged or unit-specific effects, with direct and spillover
   components. For dynamic models, simulated response paths trace effects
   over time.
3. **Validation as a release criterion.** No model variant is documented
   as supported until it passes simulation-based calibration (Talts et al.
   2018; Modrák et al. 2025). The evidence is complemented by an
   independent reference sampler, a Monte Carlo study, a replication of
   the Wang & Kockelman (2009) simulation and cross-checks against other
   packages.
4. **Evidence on existing software.** The same calibration test shows that
   `spatialprobit`'s SAR probit sampler fails on about half of the
   simulated data sets and understates $\rho$ on the rest (Section 5.3).
   Its ordered sampler, where it works, agrees with ours, but per second
   it gives about 1/1,800 of the effective sample size for the cut-points
   (Section 5.6).

The rest of the paper is organised as follows:

- **Section 2** defines the models;
- **Section 3** describes the sampler;
- **Section 4** covers probabilities, impacts, dynamic responses and model
  checking;
- **Section 5** presents the validation;
- **Section 6** shows how to use the package;
- **Section 7** applies it to [TODO: the empirical application];
- **Section 8** discusses limitations and extensions;
- **Section 9** concludes.


# 2. Models

## 2.1 Latent model and observation rules

We observe units $i = 1, \dots, N$ (regions, firms, countries) in periods
$t = 1, \dots, T$, so there are $n = NT$ observations. Stacking the latent
values period by period, $z = (z_1', \dots, z_T')'$ with $z_t \in
\mathbb{R}^N$, the general model is

$$
z_t = \rho W_t z_t + X_t\beta + W_t X^{(d)}_t\theta + \alpha + u_t,
\qquad u_t = \lambda M_t u_t + \varepsilon_t,
\qquad \varepsilon_{it} \sim N(0, v_{it}),
\tag{1}
$$

where $W_t$ and $M_t$ are $N \times N$ spatial weight matrices with zero
diagonals (usually row-standardised, and usually equal), $X_t$ holds the
regressors, $X^{(d)}_t$ the subset with spatially lagged (Durbin) terms,
and $\alpha \in \mathbb{R}^N$ the unit effects. Dynamic terms are added in
Section 2.4. In the default probit case $v_{it} = 1$.

The outcome is observed through the latent value. For a **binary** outcome,
$y_{it} = 1\{z_{it} > 0\}$. For an **ordered** outcome with $J \geq 3$
categories,

$$
y_{it} = j \iff \kappa_{j-1} < z_{it} \leq \kappa_j,
\qquad \kappa_0 = -\infty,\ \kappa_1 = 0,\ \kappa_J = \infty,
$$

with free cut-points $\kappa_2 < \dots < \kappa_{J-1}$. Fixing $\kappa_1 = 0$
identifies the location, because $X$ contains an intercept. The summary also
reports the equivalent parameterisation without an intercept used by
`MASS::polr` and `ordinal::clm`.

Write $S_t(\rho) = I_N - \rho W_t$ and $R_t(\lambda) = I_N - \lambda M_t$.
Multiplying (1) through gives

$$
R_t\big(S_t z_t - X_t\beta - W_tX^{(d)}_t\theta - \alpha\big) = \varepsilon_t,
\tag{2}
$$

which is the form the sampler works with (Section 3.1).

## 2.2 Spatial specifications

The argument `model` selects which of $\rho$, $\lambda$ and $\theta$ are
free. Table 1 lists the seven specifications and the kind of spillover each
implies for the impacts of Section 4.2.

Table 1: Spatial specifications.

| `model` | $\rho$ | $\lambda$ | $\theta$ | Spillovers |
|---|---|---|---|---|
| `"sar"` | free | 0 | 0 | global, through $S^{-1}$ |
| `"sdm"` | free | 0 | free | global, through $S^{-1}(\beta_r I + \theta_r W)$ |
| `"sem"` | 0 | free | 0 | none |
| `"sarar"` | free | free | 0 | global |
| `"sdem"` | 0 | free | free | local, through $\theta_r W$ |
| `"slx"` | 0 | 0 | free | local |
| `"none"` | 0 | 0 | 0 | none |

The spatial lag $\rho$ makes a unit's propensity depend on its neighbours'
contemporaneous propensities, so a change anywhere propagates through the
whole network. The spatial error $\lambda$ correlates the unobserved shocks
without transmitting changes in the regressors. The Durbin terms $\theta$
let neighbours' regressors enter directly. In SARAR, $\rho$ and $\lambda$
are identified separately only through the regressors when $M = W$; the
package warns when $M$ is identical to $W$ and there are no Durbin terms
(LeSage & Pace 2009).

## 2.3 Unobserved heterogeneity

The argument `effects` specifies $\alpha$.

- **Pooled:** $\alpha = 0$. A pooled SAR probit is the model of
  `spatialprobit` applied to the block-diagonal weights $I_T \otimes W$.
- **Random effects** (default): $\alpha_i \sim N(0, \sigma_\alpha^2)$
  independently (Kakamu & Wago 2005; Baltagi, Egger & Kesina 2018).
- **Correlated random effects:** the random effect may depend on the
  regressors, $\alpha_i = \bar x_i'\pi + a_i$ (Mundlak 1978), or on all
  their period values (Chamberlain 1982). Both are implemented by appending
  columns to $X$, so $\pi$ is part of $\beta$ and the sampler is unchanged.
- **Spatial random effects:** $\alpha = \phi W\alpha + a$ with
  $a \sim N(0, \sigma_\alpha^2 I_N)$, the spatially autocorrelated regional
  effects of Wang & Kockelman (2009).

Time effects enter as period dummies in $X$; with $N$ observations per
period they raise no incidental-parameter problem.

We deliberately do not offer fixed effects. With $T$ fixed, unit intercepts
in a probit model are inconsistent (the incidental-parameter problem of
Neyman & Scott 1948), and putting flat priors on them does not remove the
bias. Correlated random effects are the standard alternative in nonlinear
panel models (Wooldridge 2010, chap. 15).

## 2.4 Dynamics

The argument `dynamic = psp_dynamic()` adds persistence in one of two forms.

**State dependence.** The previous outcome enters the latent model,

$$
z_t = \rho W_t z_t + \gamma\, y_{t-1} + \eta\, W_t y_{t-1} + X_t\beta + \dots,
$$

with $\gamma$ the coefficient of the own lagged outcome and $\eta$ that of
the neighbours' lagged outcomes (a lagged diffusion channel; optional). For
ordered outcomes the lag enters as one indicator per category above the
first. The lagged terms are data columns, so the sampler is unchanged, and
the first period serves only as the initial condition.

Because $\alpha_i$ also shaped the first outcome $y_{i1}$, the lagged
outcome is correlated with $\alpha_i$, and treating $y_{i1}$ as exogenous
attributes part of the heterogeneity to state dependence. Following
Wooldridge (2005), the default models the random effect given the initial
condition,

$$
\alpha_i = \xi_0\, y_{i1} + \bar x_i'\xi + a_i,
$$

again by appending columns (`init_y` and the unit means of the regressors).
The Monte Carlo results of Section 5.4 show the size of the bias this
removes.

**Latent dynamics.** With `psp_dynamic(latent = TRUE)` the persistence is in
the latent propensity instead (Wang & Kockelman 2009):

$$
z_t = \rho W z_t + \gamma\, z_{t-1} + X_t\beta + \alpha + u_t,
\qquad t \geq 2,
$$

and the first period has no lagged term. The outcome does not feed back;
what persists is the unobserved propensity, including the part of the
previous shock that the outcome did not reveal. The weights must be time
invariant, and latent dynamics cannot be combined with lagged outcomes.

## 2.5 Error distribution, weights and missing outcomes

**Student-$t$ errors.** With `errors = psp_t(df)` the errors follow the scale
mixture $\varepsilon_{it} \mid v_{it} \sim N(0, v_{it})$, $v_{it} \sim
\text{IG}(\nu/2, \nu/2)$, which gives a robit-type link that is robust to
outlying latent values (Geweke 1993; Liu 2004). The degrees of freedom
$\nu$ are fixed by the user, which keeps the scale identified.

**Time-varying weights.** $W_t$ and $M_t$ may differ across periods, for
example when borders or trade links change. A list of identical matrices is
treated as one matrix.

**Missing outcomes.** Missing values of $y_{it}$ are allowed in static
models. Their latent values are left unrestricted and sampled with the
others, which imputes them and keeps the spatial structure intact. The
regressors must be complete, and dynamic models need complete outcomes,
because a missing outcome would make the next period's lagged regressor
missing.

## 2.6 Identification and priors

The error variance is fixed at one (the probit normalisation), and the
location is fixed by the intercept and $\kappa_1 = 0$. For row-standardised
weights the admissible range of $\rho$ is $(1/\omega_{\min}, 1)$, with
$\omega_{\min}$ the smallest eigenvalue of $W$; the default range is
$(-1, 1)$, as in `spatialprobit`, and can be changed, and the same applies
to $\lambda$ and $\phi$.

All default priors are proper, which simulation-based calibration requires
(Table 2). They can be changed with `psp_priors()`.

Table 2: Default priors.

| Parameter | Prior |
|---|---|
| $\beta$, $\theta$, $\gamma$, $\eta$ | $N(0, 10^2 I)$ |
| $\rho$, $\lambda$, $\phi$ | Beta$(1.01, 1.01)$ rescaled to the admissible range |
| $\sigma_\alpha$ | half-Cauchy$(0, 1)$ (Gelman 2006); inverse gamma on $\sigma_\alpha^2$ as an option |
| $\kappa_2 < \dots < \kappa_{J-1}$ | independent $N(0, 10^2)$ restricted to the ordering; flat as an option |

The Beta$(1.01, 1.01)$ prior is the convention of LeSage & Pace (2009):
nearly flat, but with no mass at the boundary. The inverse gamma option on
$\sigma_\alpha^2$ and the flat cut-point option allow comparison with
Kakamu & Wago (2005), Baltagi, Egger & Kesina (2018) and Albert & Chib
(1993).

Binary outcomes carry little information about $\rho$. The panel dimension
helps, but with very unbalanced outcomes the posterior of $\rho$ can be
sensitive to its prior, a point we return to in Section 6.4.

# 3. Posterior sampling

The sampler is a Gibbs sampler with data augmentation (Albert & Chib 1993;
LeSage 2000), extended to panels as in Kakamu & Wago (2005). Its design
departs from existing implementations in three ways:

- every model is handled by one sparse kernel;
- the spatial parameter $\rho$ is drawn exactly from its conditional with
  $\beta$ integrated out;
- the cut-points move jointly with the latent values.

## 3.1 One kernel for all models

Given the parameters, (2) says that $Az - b = \varepsilon$, with

$$
A = \operatorname{blockdiag}_t(R_t S_t),
\qquad b = R\,(X\beta + WX^{(d)}\theta + \Delta\alpha),
$$

where $R = \operatorname{blockdiag}_t(R_t)$ and $\Delta = \iota_T \otimes
I_N$ maps the unit effects to observations. The latent values therefore
have the density

$$
p(z \mid \cdot) \propto |A| \exp\!\Big(-\tfrac12 (Az - b)'V^{-1}(Az - b)\Big),
\qquad V = \operatorname{diag}(v_{it}),
\tag{3}
$$

restricted to the box that the outcomes define, with
$\log|A| = \sum_t \big(\log|S_t| + \log|R_t|\big)$.

Every specification in Table 1, every effects option and both dynamic forms
are special cases of (3): they differ only in how $A$ and $b$ are built.
The sampler therefore needs four operations:

1. a sparse $A$, whose values change with $\rho$ and $\lambda$ but whose
   sparsity pattern is fixed;
2. the residual $r = Az - b$;
3. access to the columns of $A$;
4. precomputed log-determinant grids.

Adding a model variant means saying how to build $A$ and $b$, not writing a
new sampler. Neither the $n \times n$ precision $A'V^{-1}A$ nor the
Kronecker product $I_T \otimes W$ is ever formed.

One iteration updates the blocks in the following order:

1. the latent values $z$;
2. $\rho$, with $\beta$ integrated out;
3. $\beta$;
4. $\lambda$;
5. the effects $\alpha$, then $\phi$;
6. $\sigma_\alpha^2$;
7. the cut-points (ordered outcomes);
8. a scale move.

With $t$ errors, the mixing variances $v$ are drawn after $\alpha$. The
order matters: step 2 marginalises $\beta$, so $\beta$ must be redrawn
before any other block conditions on it. This is the partially collapsed
Gibbs rule of van Dyk & Park (2008); breaking it changes the stationary
distribution.

## 3.2 Latent block

The latent values are updated one coordinate at a time from their exact
full conditionals. For coordinate $k$, with $a_k$ the $k$-th column of $A$,
the conditional distribution is a truncated normal with

$$
h_k = a_k'V^{-1}a_k, \qquad
m_k = z_k - \frac{a_k'V^{-1}r}{h_k}, \qquad
z_k \mid \cdot \sim TN\big(m_k, h_k^{-1}; \ell_k, u_k\big),
$$

where $(\ell_k, u_k)$ are the bounds implied by $y_k$ (and the cut-points).
After the draw, the residual is updated in place, $r \leftarrow r +
a_k(z_k^{\text{new}} - z_k)$.

Three properties follow:

- **Cost.** A full sweep costs $O(\operatorname{nnz}(A))$: no linear system
  is solved and no precision matrix is formed. Existing implementations
  compute a sparse QR factorisation and a solve in every iteration
  (Wilhelm & de Matos 2013).
- **Exactness.** Each update starts from the current state and draws from
  the exact conditional, so one sweep is a valid Gibbs step.
  `spatialprobit` with more than one inner iteration restarts its
  truncated-normal chain in each iteration, which makes its latent update
  only approximately a draw from the conditional.
- **Missing outcomes** get infinite bounds and need no special case.

The univariate truncated normal draws use normal rejection when the
interval holds enough mass, and otherwise the exponential (one-sided tails)
or uniform (short intervals) rejection samplers of Robert (1995). They are
stable for means many standard deviations
beyond a bound.

## 3.3 Spatial parameters and coefficients

**$\rho$.** With the prior $\beta \sim N(c, \Sigma_\beta)$ integrated out,
and $\tilde X = RX$ (including the Durbin columns), define the stacked
vectors $a_0 = R(z - \Delta\alpha) - \tilde Xc$ and $a_1 = RWz$, and

$$
q_{ab} = a_a'V^{-1}a_b - (\tilde X'V^{-1}a_a)'\,B\,(\tilde X'V^{-1}a_b),
\qquad B = (\tilde X'V^{-1}\tilde X + \Sigma_\beta^{-1})^{-1}.
$$

Then

$$
\log p(\rho \mid z, \alpha, \lambda, y) =
\sum_t \log|I - \rho W_t| - \tfrac12\big(q_{00} - 2\rho\,q_{01} + \rho^2 q_{11}\big)
+ \log \pi(\rho) + \text{const}.
\tag{4}
$$

The quadratic form needs three scalars, computed in $O(nk)$. We draw from
(4) on a two-stage grid: a coarse grid locates the mode, and a fine grid
over the mode ± 8 standard deviations gives the draw, by the inverse CDF
with linear interpolation within cells, so the draw is continuous. The
log-determinant is interpolated by cubic splines from a grid computed once
before sampling (Section 3.8).

Expression (4) is exact because the error variance is one. The
implementation in `spatialprobit` instead uses the concentrated form for a
linear SAR model, $-\tfrac{n-k}{2}\log Q(\rho)$, which integrates out an
error variance that a probit model does not have. Section 5.3 shows that
`spatialprobit`'s posterior of $\rho$ is miscalibrated; this kernel is one of
three candidate causes, and which of them is responsible has not been
isolated.

**$\beta$.** Given $\rho$, the conditional of $\beta$ is Gaussian,
$N\big(B(\tilde X'V^{-1}\tilde y + \Sigma_\beta^{-1}c), B\big)$ with
$\tilde y = R(Sz - \Delta\alpha)$.

**$\lambda$.** Given $\beta$, let $e = Sz - X\beta - \Delta\alpha$. The
kernel is again quadratic in $\lambda$,
$e'V^{-1}e - 2\lambda\,e'V^{-1}Me + \lambda^2(Me)'V^{-1}(Me)$, plus
$\sum_t\log|I - \lambda M_t|$, and is drawn on the same kind of grid. It is
not collapsed over $\beta$, because the design $RX$ depends on $\lambda$.

## 3.4 Random effects

**Non-spatial errors** (SAR, SDM, SLX, none). Without the error operator
$R$, the effects are conditionally independent across units. With
$e_{it} = (S_t z_t - X_t\beta)_i$,

$$
\alpha_i \mid \cdot \sim N\!\left(
\frac{\sum_t e_{it}/v_{it}}{\sum_t 1/v_{it} + \sigma_\alpha^{-2}},\
\frac{1}{\sum_t 1/v_{it} + \sigma_\alpha^{-2}}\right),
$$

at cost $O(n)$.

**Spatial errors or spatial effects.** Otherwise $\alpha$ is drawn jointly
from $N(P^{-1}h, P^{-1})$ with the sparse precision

$$
P = \sum_t R_t'V_t^{-1}R_t + \sigma_\alpha^{-2} B_\phi'B_\phi,
\qquad B_\phi = I - \phi W
$$

($B_\phi = I$ without spatial effects). We use a sparse Cholesky
factorisation whose symbolic analysis is done once, because the pattern of
$P$ is fixed. Under spatial effects, $\phi$ is drawn on the grid of $W$ from
its conditional given $\alpha$ and $\sigma_\alpha$.

**$\sigma_\alpha^2$.** The half-Cauchy prior is the mixture
$\sigma_\alpha^2 \mid \xi \sim \text{IG}(1/2, 1/\xi)$,
$\xi \sim \text{IG}(1/2, 1/A^2)$ (Wand et al. 2011), which gives two
conjugate inverse gamma draws:
$\sigma_\alpha^2 \mid \alpha, \xi \sim \text{IG}\big(\tfrac{N+1}{2},
\tfrac12\|B_\phi\alpha\|^2 + 1/\xi\big)$ and
$\xi \mid \sigma_\alpha^2 \sim \text{IG}(1, 1/\sigma_\alpha^2 + 1/A^2)$.

## 3.5 Cut-points

In the Albert & Chib (1993) Gibbs step, $\kappa_j$ is drawn uniformly (times
its prior) between the largest latent value in category $j$ and the
smallest in category $j + 1$. That interval shrinks like $O(1/n)$, so the
cut-points barely move in large panels (Cowles 1996). The collapsed
Metropolis step of Cowles (1996) needs $P(y \mid \kappa, \cdot)$ with $z$
integrated out, which in a spatial model is an $n$-dimensional normal
rectangle probability and is not available exactly.

We therefore keep the Albert–Chib step and add, for each free cut-point, a
joint transformation move on $(\kappa_j, z)$ in the generalised Gibbs
framework of Liu & Sabatti (2000). The move proposes a new value $\kappa_j'$
and maps the latent values piecewise-linearly so that every $z$ stays in
its category.

- **Interior cut-points** ($j < J - 1$). Random walk on the logit of the
  relative position $p = (\kappa_j - \kappa_{j-1})/(\kappa_{j+1} -
  \kappa_{j-1})$. Category $j$ is rescaled from $(\kappa_{j-1}, \kappa_j]$
  onto $(\kappa_{j-1}, \kappa_j']$ by the factor $s_1$, and category
  $j + 1$ from $(\kappa_j, \kappa_{j+1}]$ onto $(\kappa_j', \kappa_{j+1}]$ by
  the factor $s_2$.
- **The top cut-point.** Random walk on $\log(\kappa_{J-1} -
  \kappa_{J-2})$. Category $J - 1$ is rescaled and category $J$ is
  translated.

Let $\delta = z' - z$ and $n_j$ be the number of observations in category
$j$. For an interior cut-point the move is accepted with probability
$\min(1, r)$, where

$$
\log r = -\tfrac12\Big(2\,r'V^{-1}A\delta + (A\delta)'V^{-1}(A\delta)\Big)
+ n_j \log s_1 + n_{j+1}\log s_2
+ \log\frac{p'(1 - p')}{p(1 - p)}
+ \log\frac{\pi(\kappa_j')}{\pi(\kappa_j)}.
$$

The first term is the change in the kernel (3), because the
log-determinant does not depend on $\kappa$. The next two are the Jacobian
of the latent map, and the fourth that of the logit reparameterisation.
For the top cut-point the Jacobian is $(n_{J-1} + 1)\log s_1$. One move
costs $O(\operatorname{nnz}(A))$. The step sizes adapt during warm-up only,
towards an acceptance rate of about 0.35, and are then frozen.

Both samplers pass SBC. At $N = 400$, $T = 5$ the combined move gives 4
times the ESS per second of Albert–Chib alone for $J = 3$, and 18 to 78
times for the interior cut-points with $J = 5$, at about 15% more time per
iteration. It is the default; `control = psp_control(cutpoint_sampler =
"albert_chib")` selects the Albert–Chib step alone.

## 3.6 A scale move

Data augmentation for probit models mixes slowly along the overall scale
of $(z, \beta, \alpha)$ (Liu & Wu 1999; Imai & van Dyk 2005). After each
sweep we apply the generalised Gibbs move $x \mapsto gx$ to
$x = (z, \beta, \alpha, \kappa)$. This leaves the outcomes unchanged, so
only the prior and the kernel (3) determine the target of $g$. With the
Haar measure $dg/g$, the target is

$$
p(g) \propto g^{d-1}\exp\!\big(-\tfrac12 C g^2 + L g\big),
$$

where $d$ is the number of rescaled coordinates and $C$ and $L$ collect the
quadratic and linear terms of the kernel and the priors. When the prior
means are zero, $L = 0$ and $g^2 \sim \text{Gamma}(d/2, \text{rate } C/2)$
exactly. Otherwise a random-walk Metropolis step on $\log g$ is used. The
move roughly doubled the ESS of the slopes and of $\sigma_\alpha$.

## 3.7 Latent dynamics and $t$ errors

**Latent AR(1).** The coefficient $\gamma$ is treated as the coefficient of
a design column `lag_z` that holds $z_{t-1}$ and is refreshed from the
current latent values before every residual computation. Hence $\gamma$ is
drawn with $\beta$, and $\rho$ remains collapsed over it. The full operator
$A$ gains blocks $-\gamma R_t$ below the diagonal. It is block
lower-triangular, so $\log|A|$ is unchanged. The latent sweep and the
cut-point move include the coupling of $z_t$ to period $t + 1$, and the
scale move leaves $\gamma$ unscaled, because its column is already
rescaled.

**Student-$t$ errors.** The precision weights $w_{it} = 1/v_{it}$ enter
every quadratic form. The mixing variances are drawn as
$v_{it} \mid \cdot \sim \text{IG}\big(\tfrac{\nu+1}{2},
\tfrac{\nu + r_{it}^2}{2}\big)$, after which the weighted cross-products
are refreshed.

## 3.8 Computational cost

**Log-determinants.** $\log|I - \rho W|$ is computed once on a grid, before
sampling. For $N$ up to about 3,000 it uses the eigenvalues of $W$
(exact, using the real part of the complex logarithm for asymmetric $W$).
For larger $N$ it uses sparse LU factorisations at grid points (Pace &
Barry 1997), interpolated by cubic splines. With time-varying weights the
grids are summed over periods, and identical periods are computed once.

**Per-iteration cost.** The cost is $O(T\,\operatorname{nnz}(W) + nk)$
without spatial errors or spatial effects, plus the numeric factorisation
of the sparse $N \times N$ matrix $P$ otherwise. No $n \times n$ object is
formed. On one core, an iteration takes 1.3 ms for a SAR random-effects
probit with $N = 1{,}000$, $T = 10$ and $k = 5$, and 33 ms with
$N = 10{,}000$ and $T = 20$ ($n = 200{,}000$), on the computer described in
the computational details (Section 9). The factorisation roughly doubles
the cost for a SEM model (Table 5). [TODO: rerun on the release version.]

**Chains and reproducibility.** Each chain is one C++ call with its own
L'Ecuyer-CMRG stream, so results do not depend on the number of cores, and
the user's random number state is restored afterwards. By default four
chains are run. The fit reports the rank-normalised split-$\hat R$ and the
bulk and tail ESS (Vehtari et al. 2021), and it warns when convergence is
doubtful.

# 4. Interpretation

## 4.1 Probabilities

Let $m_t = S_t^{-1}(X_t\beta + W_tX^{(d)}_t\theta)$ be the systematic part
of the reduced form, excluding $\alpha$. The reduced-form variance of the
latent value is

$$
s_{it}^2 = \sigma_\alpha^2\,[S_t^{-1}B_\phi^{-1}B_\phi^{-\top}S_t^{-\top}]_{ii}
+ [K_t V_t K_t']_{ii},
\qquad K_t = S_t^{-1}R_t^{-1}.
$$

For a SAR model with probit errors this reduces to
$(1 + \sigma_\alpha^2)[(S'S)^{-1}]_{ii}$.

- The **population-averaged** probability integrates over the effects:
  $P(y_{it} = 1) = \Phi(m_{it}/s_{it})$, and, for ordered outcomes,
  $P(y_{it} = j) = \Phi\big((\kappa_j - m_{it})/s_{it}\big) -
  \Phi\big((\kappa_{j-1} - m_{it})/s_{it}\big)$.
- The **unit-specific** probability conditions on the sampled effects. It
  adds $S_t^{-1}\alpha$ to the mean and drops the $\sigma_\alpha^2$ term
  from the variance.

With $t$ errors, the probabilities average the normal formula over
simulated mixing variances. With latent dynamics, the mean follows the
recursion $m_t = S^{-1}(X_t\beta + \gamma m_{t-1})$, and the variance
accumulates the lagged shocks, so it grows with $t$.

## 4.2 Impacts

The effect of regressor $r$ of unit $j$ on the probability of unit $i$ in
the same period is

$$
\frac{\partial P(y_{it} = 1)}{\partial x_{jt,r}}
= \frac{\phi(m_{it}/s_{it})}{s_{it}}
\big[S_t^{-1}(\beta_r I + \theta_r W_t)\big]_{ij}.
\tag{5}
$$

Following LeSage & Pace (2009), we summarise the $N \times N$ matrix in (5)
by three scalars, averaged over units and periods:

- the **average direct impact**, the mean of the diagonal;
- the **average total impact**, the mean row sum;
- the **average indirect impact** (the spillover), their difference.

For ordered outcomes the factor $\phi(\cdot)/s$ becomes
$\big[\phi\big((\kappa_{j-1} - m)/s\big) - \phi\big((\kappa_j - m)/s\big)\big]/s$
for category $j$, and the category impacts sum to zero. Each impact is
computed for every retained posterior draw, so the reports give posterior
means, standard deviations and credible intervals. Per-period impacts are
available as well.

Two choices differ from common practice:

- **The reduced-form scale.** `spatialprobit` evaluates $\phi(m_{it})$
  without dividing by $s_{it}$ (Wilhelm & de Matos 2013), which ignores the
  heteroskedasticity that $S^{-1}$ induces and the random effects. We use
  (5) by default and offer the other convention as
  `scale = "lesage_pace"`, for comparison with published results. With
  random effects the two can differ substantially.
- **Population-averaged versus unit-specific.** By default the effects are
  integrated out (`type = "population"`). With `type = "unit"` the impacts
  condition on the sampled effects, whose draws the sampler stores for a
  subset of iterations. These are the effects for the units in the data
  rather than for a new draw from the population, and they are larger.

Computing (5) needs the diagonals of $S^{-1}$, $S^{-1}W$ and of the
reduced-form covariance for every draw. For $N$ up to 500 (by default) they
are computed exactly from a dense inverse. Above that we use a stochastic
diagonal estimator with Rademacher probes (Bekas, Kokiopoulou & Saad 2007),
each needing one sparse solve. The total impact needs only one solve,
$N^{-1} g'S^{-1}(\beta_r\iota + \theta_r W\iota)$ with
$g_i = \phi(m_i/s_i)/s_i$.

For a discrete change, for example a dummy regressor switched on,
`predict(newdata = )` computes the probabilities under any counterfactual
values of the regressors for the same units and periods.

## 4.3 Dynamic responses

With lagged outcomes, (5) is the short-run effect: a change in period $t$
also changes $y_t$, and through it the probabilities in later periods. In a
nonlinear spatial model there is no closed-form long-run multiplier on
probabilities. `psp_response_path()` therefore simulates outcome histories
forward from a chosen period, with and without the change, using common
random numbers. Given each history, the probabilities are computed exactly
from the reduced form, so only the history is simulated. The difference is
averaged over histories and posterior draws. The change can be temporary
(one period) or permanent, and it can apply to all units or to some of
them. In the latter case the effects on the changed units and on the
others are reported separately.

With latent dynamics the responses follow exactly from the reduced-form
recursion of Section 4.1, so they are computed without simulation.

## 4.4 Model checking and comparison

**Posterior predictive checks.** `psp_ppc()` simulates replicated outcomes
from the posterior and compares them with the data, period by period. The
statistics are prevalence or category shares, persistence (the share of
units whose outcome repeats from the previous period) and Moran's I of the
outcome, a measure of spatial clustering. It reports posterior predictive
p-values.

**Likelihood.** The observed-data likelihood is an $n$-dimensional normal
rectangle probability. We approximate it with the GHK simulator (Geweke
1991; Hajivassiliou & McFadden 1998; Keane 1994) applied to the sparse
precision of $(z, \alpha)$. The coordinates of $\alpha$ and of missing
outcomes are unrestricted, so we order them first and GHK integrates them
out exactly. Periods that are independent given the parameters (pooled
models) are simulated as separate blocks. The approximation is exact in
the non-spatial case, which we check against quadrature. `psp_loglik()`,
`logLik()` and `psp_dic()` (Spiegelhalter et al. 2002) use it and are
labelled approximate. For random-effects models the marginal likelihood is
a single $n$-dimensional problem, and the estimate becomes noisy beyond a
few hundred observations, which a warning reports.

# 5. Validation

An MCMC sampler can run without errors, produce plausible numbers and still
target the wrong distribution: a conditional density off by one term, a
block order that breaks a collapsed Gibbs step, or a move that does not
preserve the target. Such errors produce no error message and can be found
only by testing the output against known answers. We use four layers of
evidence, from the most to the least specific:

1. tests of the individual components;
2. an independent reference sampler;
3. simulation-based calibration (SBC) of every supported variant;
4. the frequentist behaviour of the estimates: a Monte Carlo study, the
   replication of a published simulation, and cross-checks against other
   packages.

A model variant is documented as supported only after it passes SBC. The
scripts and full results are in the `validation/` directory of the source
repository, and `vignette("validation")` reproduces the figures and tables
of this section.

## 5.1 Components and a reference sampler

The unit tests, which run on CRAN, check each component against an
independent computation:

- the truncated normal sampler against analytic moments, including bounds
  many standard deviations from the mean;
- the collapsed kernel (4) against a brute-force dense marginalisation of
  $\beta$;
- the random-effects conditionals against dense Gaussian conditioning;
- the log-determinant grids against `determinant()` and `spatialreg`;
- the impacts (5) against finite differences of the exact probabilities;
- the stochastic diagonal estimator against the exact diagonals;
- the GHK likelihood against quadrature (random effects) and `mvtnorm`
  (spatial models).

The cut-point moves have their own stationarity test: started from draws of
a known joint distribution of $(\kappa, z)$ on a toy model, one step of the
kernel must leave that distribution unchanged.

The test suite also contains a second sampler written in plain R that
shares no sampling code with the C++ engine and deliberately uses different
algorithms. It uses dense precision matrices instead of residual updates,
random-walk Metropolis steps with exact log-determinants for $\rho$,
$\lambda$, $\phi$ and the latent AR coefficient instead of grid draws, and
no scale or cut-point moves. Both samplers are run on the same small data
sets, covering SAR, SEM and SARAR models, spatial random effects,
Student-$t$ errors, time-varying weights and latent dynamics, with binary
and ordered outcomes. Every posterior mean must agree within four Monte
Carlo standard errors. For a mistake in one sampler to go undetected, the
other would have to make the same mistake by a different route. All the
sampler moves of Section 3, and the three moves that were tried and removed
(Section 8), passed this comparison.

## 5.2 Simulation-based calibration

**Method.** SBC (Talts et al. 2018; Modrák et al. 2025) checks a sampler
against its own model. Parameters are drawn from the prior and a data set is
simulated from them. The model is then fitted, and the rank of each true
value among the posterior draws is recorded. If the sampler targets the
correct posterior, the ranks are uniformly distributed whatever the data.
An overconfident sampler gives U-shaped rank histograms, and a biased one
gives skewed histograms. Because the test averages over the prior, the
priors must be proper and the simulated data must be fitted with exactly
the model that generated them.

**Design.** We used a 4 × 4 rook grid ($N = 16$) and $T = 4$ periods.
The priors were proper and moderately informative, so that the simulated
data sets were realistic:

- $\beta \sim N(0, 1)$;
- $(\rho + 1)/2 \sim \text{Beta}(2, 2)$, and the same for $\lambda$ and
  $\phi$;
- $\sigma_\alpha^2 \sim \text{IG}(3, 1)$;
- cut-points $\sim N(1, 1)$, restricted to the ordering.

Each variant was run for 500 replications. Each fit used two chains and
kept 400 thinned draws. Uniformity of each parameter's ranks was tested
with a $\chi^2$ test on 10 bins.

**Variants.** The 30 variants cover:

- every spatial specification with random effects, binary outcomes;
- pooled SAR and SEM;
- Mundlak and Chamberlain correlated random effects;
- spatial random effects;
- ordered outcomes with $J = 3$ and 4, under both cut-point samplers;
- lagged outcomes with Wooldridge initial conditions, and latent AR(1)
  dynamics;
- Student-$t$ errors;
- time-varying weights;
- 20% missing outcomes, binary and ordered;
- combinations such as SARAR with ordered outcomes, $t$ errors and
  time-varying weights.

**Results.** In all, 32 runs of the 30 variants gave 199 rank tests.
Under correct calibration their p-values are themselves uniform, so about
10 should fall below 0.05 and 2 below 0.01 by chance. We observed 19 and 3.
[FIGURE: histogram and uniform QQ plot of the 199 p-values, from
`vignette("validation")`.]

The three p-values below 0.01 all come from known slow mixing rather than
a wrong target:

- **Latent AR(1) with SAR** ($p = 0.001$ for $\gamma$ and $0.008$ for
  $x_2$, at the default chain length). The ranks of $\gamma$ are U-shaped,
  and the departure is concentrated in the replications whose chains had
  not converged (largest $\hat R \geq 1.02$, about 45% of them). The
  coefficient is drawn given the latent values and mixes slowly against
  $\rho$. With chains three times as long, the same variant passes
  (smallest $p = 0.090$).
- **A cut-point of SARAR with time-varying weights** ($p = 0.010$). A
  replication with different seeds passes ($p = 0.874$ for the same
  parameter).

SARAR models with $M = W$ also give borderline p-values for $\rho$ and
$\lambda$ (around 0.02–0.03) at the default chain length. The two
parameters are nearly interchangeable when the weights coincide, so they
mix slowly, and doubling the chain length removes the departure for
$\lambda$. The documentation therefore asks for longer chains for latent
dynamics and for SARAR models with $M = W$, and the convergence warning
flags such fits.

[FIGURE: rank histograms of $\gamma$ at the default and at three times
the chain length, and of six central parameters.]

[TABLE: per variant, the number of parameters, the smallest p-value and the
median of the largest $\hat R$; from `vignette("validation")`. With several
parameters per variant, the smallest p-value falls below 0.05 for roughly a
fifth to a third of the variants under correct calibration.]

## 5.3 The same test applied to `spatialprobit`

To see whether the design choices of Section 3 matter in practice, we
applied the same test to `spatialprobit::sarprobit` (version 1.0.4). We used
an 8 × 8 rook grid, $T = 1$ and a pooled SAR probit, with priors
$\beta \sim N(0, I)$ and $(\rho + 1)/2 \sim \text{Beta}(2, 2)$ set through
its `prior` argument, for 500 replications.

The results:

- **Failures.** The sampler returned undefined (NaN) draws in 264 of the
  500 replications.
- **Miscalibration.** Among the other 236, the ranks of $\rho$ are clearly
  non-uniform ($\chi^2$ $p < 0.0001$). They pile up in the top bin (76
  against an expected 24; mean normalised rank 0.61), so the posterior
  understates $\rho$. The intercept is borderline ($p = 0.037$).
- **Control.** `panelspatialprobit` on the identical design had no failures,
  and all its rank tests passed (all $p > 0.08$; $p = 0.55$ for $\rho$).

[FIGURE: rank histograms of $\rho$ for both packages.]

Reading the source of `spatialprobit` suggests three candidate causes, each
corresponding to a design choice in Section 3:

1. its draw of $\rho$ uses the concentrated likelihood of the linear SAR
   model, which integrates out an error variance that the probit fixes at
   one (Section 3.3);
2. its block order draws $\beta$, then $\rho$ with $\beta$ integrated out,
   and then the latent values using the earlier $\beta$, which breaks the
   partially collapsed Gibbs rule (Section 3.1);
3. its truncated normal step restarts its inner chain in every iteration,
   so the latent update is not an exact conditional draw (Section 3.2).

We have not isolated which of these causes the bias. Keeping only the
successful replications could itself distort the ranks, so the size of the
miscalibration is uncertain, but its presence is not.

## 5.4 Monte Carlo study

SBC averages over the prior. A Monte Carlo study instead fixes the true
parameters and asks how the estimates behave in repeated samples: the bias
and root mean squared error (RMSE) of the posterior mean, and the coverage
of 90% credible intervals. Coverage should be close to 90% when the priors
are weak relative to the data.

Each design was replicated 200 times on rook grids, with $T = 10$ periods
and 2 chains of 2,000 iterations after 1,000 of warm-up. The static designs
are:

- a binary SAR probit with random effects at $N = 100$ and $N = 400$;
- an ordered SAR probit with $J = 4$ at $N = 100$.

In all three, $\beta = (0, 1, -1)$, $\rho = 0.5$ and $\sigma_\alpha = 0.7$.
For the impacts, the target is the total impact of $x_1$ computed from the
true parameters on the same design.

The dynamic designs share one data-generating process: a SAR probit with
state dependence ($\gamma = 0.8$), random effects ($\sigma_\alpha = 0.8$)
and $\rho = 0.4$, on a 10 × 10 grid with 10 periods after the initial one.
It is fitted three ways:

- **Wooldridge:** with random effects and Wooldridge initial conditions,
  the correct specification;
- **Pooled:** pooled, without initial conditions, which confounds state
  dependence with heterogeneity;
- **Static:** a static random-effects model, which omits the state
  dependence.

Table 3: Monte Carlo results (200 replications per design).

| Design | Parameter | True | Bias | RMSE | 90% coverage |
|---|---|---|---|---|---|
| Binary SAR, $N = 100$ | $x_1$ | 1 | 0.012 | 0.079 | 0.91 |
| | $\rho$ | 0.5 | −0.005 | 0.046 | 0.90 |
| | $\sigma_\alpha$ | 0.7 | 0.015 | 0.095 | 0.91 |
| | total impact of $x_1$ | 0.377 | −0.003 | 0.028 | 0.91 |
| Binary SAR, $N = 400$ | $x_1$ | 1 | 0.004 | 0.037 | 0.91 |
| | $\rho$ | 0.5 | −0.004 | 0.022 | 0.91 |
| | $\sigma_\alpha$ | 0.7 | 0.003 | 0.048 | 0.91 |
| | total impact of $x_1$ | 0.380 | −0.002 | 0.015 | 0.90 |
| Ordered SAR, $J = 4$, $N = 100$ | $x_1$ | 1 | 0.007 | 0.057 | 0.92 |
| | $\rho$ | 0.5 | −0.003 | 0.034 | 0.89 |
| | $\kappa_2$ | 0.8 | 0.011 | 0.060 | 0.91 |
| | $\kappa_3$ | 1.6 | 0.017 | 0.091 | 0.89 |
| | total impact of $x_1$ | 0.283 | −0.001 | 0.025 | 0.88 |
| Dynamic, Wooldridge | $\gamma$ | 0.8 | 0.010 | 0.119 | 0.93 |
| | $x_1$ | 1 | 0.040 | 0.089 | 0.89 |
| | $\rho$ | 0.4 | −0.012 | 0.053 | 0.89 |
| Dynamic, pooled | $\gamma$ | 0.8 | 0.313 | 0.333 | 0.11 |
| | $x_1$ | 1 | −0.176 | 0.187 | 0.14 |
| | $\rho$ | 0.4 | −0.025 | 0.058 | 0.85 |
| Dynamic, static | $x_1$ | 1 | −0.016 | 0.075 | 0.92 |
| | $\rho$ | 0.4 | 0.003 | 0.055 | 0.88 |

Table 3 omits $x_2$ and, outside the first design, $\sigma_\alpha$;
`vignette("validation")` lists them. They behave like $x_1$: in the static
designs every bias is at most 0.012 in absolute value and every coverage is
between 0.885 and 0.92. The one exception is $x_2$ in the Wooldridge
dynamic design, with a coverage of 0.85 (bias −0.036).

Table 3 shows three results:

- **Static models.** Biases are at most 0.02 in absolute value, and
  coverage is between 0.88 and 0.92 for every parameter and for the
  impacts. Quadrupling $N$ halves the RMSE and the interval widths, as
  expected.
- **Correct dynamic model.** The Wooldridge random-effects model recovers
  the state dependence: the bias of $\gamma$ is 0.01, with 93% coverage.
- **Misspecified dynamic models.** The pooled dynamic model attributes the
  heterogeneity to state dependence. It overstates $\gamma$ by 0.31 and
  biases the slopes towards zero by about 18%, so its intervals cover the
  truth in only 11–16% of samples. The static random-effects model omits
  the lag but still estimates the slopes and $\rho$ well, so here the
  heterogeneity, not the omitted dynamics, is what distorts the slopes.

## 5.5 Replication of the Wang and Kockelman simulation

Wang & Kockelman (2009) proposed the dynamic spatial ordered probit with
latent AR(1) dynamics and spatially autocorrelated effects, and evaluated
their sampler on simulated data. Their design has:

- three categories, with thresholds 0 and 2.1;
- four $U(0, 1)$ regressors with $\beta = (-1.7, 2, 1, 0.5)$, and 8 periods;
- latent persistence $\lambda \in \{0.1, 0.5, 0.9\}$ crossed with spatial
  dependence of the effects $\rho \in \{0.1, 0.6, 0.7, 0.9\}$;
- 50 data sets per cell.

In this subsection we follow their notation: their $\lambda$ is our
$\gamma$ (`lag_z`), their $\rho$ is our $\phi$, and their $\sigma^2$ is our
$\sigma_\alpha^2$. We fitted the design with

```r
psoprobit(y ~ x1 + x2 + x3 + x4, data, index, W, model = "none",
          effects = "spatial_random", dynamic = psp_dynamic(latent = TRUE))
```

Three adaptations were needed, because the package has one effect per unit
and unit error variance:

- **Units.** There are 300 units on a 15 × 20 queen-contiguity grid, each
  with its own effect. Their model allows this case; their own study used
  30 regions of 10 individuals sharing an effect.
- **Error variances.** All are equal to one, a special case of their
  model; theirs vary by region.
- **Thresholds.** Our fit has an intercept and a free second cut-point.
  Their thresholds follow as $\gamma_1 = -\text{intercept}$ and
  $\gamma_2 = \kappa_2 - \text{intercept}$, so the estimates are directly
  comparable.

Each fit used 2 chains of 15,000 iterations after 5,000 of warm-up, thinned
by 5. They used one chain of 2,000 iterations.

Table 4: Replication of Wang & Kockelman (2009), Table 1. Averages of
posterior means over 50 data sets per cell; true $\sigma^2 = 1$. "Slope
bias" is the largest absolute average bias of the four slopes. Coverage is
that of our 90% intervals, averaged over the nine parameters.

| $\lambda$ | $\rho$ | $\lambda$: ours | $\lambda$: W&K | $\rho$: ours | $\rho$: W&K | $\sigma^2$: ours | $\sigma^2$: W&K | Slope bias: ours | Slope bias: W&K | Coverage | Share $\hat R \geq 1.05$ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.1 | 0.095 | 0.097 | 0.087 | 0.048 | 1.04 | 1.05 | 0.031 | 0.035 | 0.89 | 0.00 |
| 0.1 | 0.6 | 0.101 | 0.101 | 0.595 | 0.452 | 1.00 | 1.09 | 0.024 | 0.028 | 0.92 | 0.00 |
| 0.1 | 0.7 | 0.098 | 0.099 | 0.695 | 0.572 | 1.03 | 1.12 | 0.022 | 0.043 | 0.91 | 0.00 |
| 0.1 | 0.9 | 0.104 | 0.097 | 0.886 | 0.845 | 1.04 | 1.83 | 0.024 | 0.188 | 0.88 | 0.66 |
| 0.5 | 0.1 | 0.499 | 0.492 | 0.088 | 0.039 | 1.01 | 1.16 | 0.055 | 0.046 | 0.89 | 0.00 |
| 0.5 | 0.6 | 0.495 | 0.507 | 0.570 | 0.494 | 1.08 | 1.22 | 0.017 | 0.076 | 0.90 | 0.00 |
| 0.5 | 0.7 | 0.490 | 0.514 | 0.672 | 0.623 | 1.11 | 1.30 | 0.009 | 0.101 | 0.88 | 0.02 |
| 0.5 | 0.9 | 0.498 | 0.511 | 0.875 | 0.863 | 1.12 | 1.77 | 0.081 | 0.278 | 0.87 | 0.64 |
| 0.9 | 0.1 | 0.902 | 0.919 | 0.093 | −0.001 | 1.05 | 1.29 | 0.037 | 0.191 | 0.91 | 0.00 |
| 0.9 | 0.6 | 0.898 | 0.921 | 0.568 | 0.498 | 1.12 | 1.23 | 0.027 | 0.140 | 0.90 | 0.02 |
| 0.9 | 0.7 | 0.905 | 0.909 | 0.681 | 0.616 | 1.04 | 1.31 | 0.043 | 0.181 | 0.88 | 0.04 |
| 0.9 | 0.9 | 0.914 | 0.863 | 0.881 | 0.855 | 1.10 | 1.50 | 0.062 | 0.154 | 0.89 | 0.70 |

None of the 600 fits failed. Our estimates are close to the truth in every
cell:

- the overall coverage of 90% intervals is 0.89;
- the largest average bias of any parameter is 0.12, for $\sigma^2$;
- $\lambda$ is within 0.015 of its true value in every cell;
- $\rho$ is slightly understated at 0.6, by about 0.03.

The estimates reported by Wang & Kockelman drift with $\lambda$ and
$\rho$:

- $\sigma^2$ reaches 1.83 against a true value of 1;
- $\rho$ is understated when it is small (0.048, 0.039 and −0.001 against
  0.1);
- the slope biases reach 0.28.

Their design (shared regional effects and unequal error variances) and
their chain length both differ from ours. The comparison therefore shows
that our sampler recovers the parameters of their model, not what caused
the drift in theirs. [TODO: check the design and Table 1 against the
published version; we used the 2008 working paper.]

When $\rho = 0.9$, about two thirds of the fits still have $\hat R \geq
1.05$ at our chain length. The slopes and $\lambda$ are not affected; the
intercept, and with it the thresholds, are. With such strong spatial
dependence, the effects are close to a level common to all units, which
trades off against the intercept. The RMSE of the thresholds is about 0.62
in these cells, against about 0.2 elsewhere. In a test fit of the hardest
cell at their chain length of 2,000 iterations, the effective sample sizes
of $\lambda$ and of the intercept were only 10–13. Strongly autocorrelated
effects therefore need long chains, a point the documentation makes.

## 5.6 Cross-checks against other packages

Where another implementation of a special case exists, the posteriors should
agree up to Monte Carlo error:

- **`MCMCpack::MCMCoprobit`** (Martin, Quinn & Park 2011). We fitted a
  non-spatial pooled ordered probit with $N = 300$, $T = 4$ and $J = 4$,
  using the same priors in both packages. Every posterior mean agrees
  within 0.03 posterior standard deviations, and the posterior standard
  deviations agree within 2%.
- **`spatialprobit::sarprobit`** on the Katrina data that ship with it: a
  cross-section of 673 businesses after Hurricane Katrina (LeSage et al.
  2011), fitted as a pooled SAR probit with near-diffuse priors. The
  posterior means of all ten parameters agree within 0.65 posterior
  standard deviations. `spatialprobit` had to be started from good values, because
  its defaults fail on these data. Given Section 5.3, exact agreement is not
  expected. On one core, `panelspatialprobit` took 2.0 s for two chains of
  4,000 iterations, against 14 s for one chain of 7,000 in `spatialprobit`.
- **`spatialprobit::sarorderedprobit`**, on a simulated pooled SAR ordered
  probit: a 20 × 20 rook grid, $T = 1$, $J = 4$, with near-diffuse priors on
  both sides. Every posterior mean agrees within 0.07 posterior standard
  deviations, and the posterior standard deviations agree within 2%. The
  largest difference is for $\rho$: 0.004, or 0.07 standard deviations,
  which is small but beyond Monte Carlo error. This is consistent with the
  different $\rho$ kernel (Section 3.3).

  The difference in mixing is large. `spatialprobit` updates its
  cut-points with the Albert–Chib step alone. With 4 chains of 50,000 draws
  it reached a bulk ESS of only about 400 for the cut-points, in 4.2
  minutes. `panelspatialprobit` reached about 6,000–7,000 with 4 chains of
  8,000 draws, in 2 s. Per second, that is about 1,800 times as much ESS
  for the cut-points.

## 5.7 Performance

Table 5 reports the time per iteration and the effective sample size (ESS)
per second for three random-effects models with $k = 5$ regressors: a SAR
probit, a SAR ordered probit with $J = 4$ and a SEM probit
(`validation/bench/bench_models.R`).

Table 5: Speed and mixing, on one core of an Intel Core i7-4790K
(Section 9). [TODO: rerun on the release version.]

| Model | $N$ | $T$ | $n$ | ms per iteration | ESS per second (range over parameters) |
|---|---|---|---|---|---|
| SAR probit | 1,000 | 10 | 10,000 | 1.3 | 28–68 |
| SAR ordered probit | 1,000 | 10 | 10,000 | 1.7 | 24–84 |
| SEM probit | 1,000 | 10 | 10,000 | 2.7 | 9–37 |
| SAR probit | 10,000 | 20 | 200,000 | 33 | 0.7–2.8 |
| SAR ordered probit | 10,000 | 20 | 200,000 | 41 | 0.5–4.0 |
| SEM probit | 10,000 | 20 | 200,000 | 74 | 0.7–1.1 |

Each model was fitted with one chain: 5,000 iterations after 5,000 of
warm-up at $N = 1{,}000$, and 1,000 after 1,000 at $N = 10{,}000$. The ESS
per second divides the bulk ESS by the whole run time, warm-up and the
log-determinant grid included. At $N = 10{,}000$ the chains are short, so
these figures are only a rough guide. The slowest parameter is the
intercept or $\sigma_\alpha$ in the SAR models, and $\lambda$ in the SEM
model.

The ordered model costs about a quarter more per iteration than the binary
one, mainly for the cut-point updates. The SEM model costs two to two and
a half times as much, mainly because it factorises the sparse $N \times N$ precision of the random
effects in every iteration, since that precision depends on $\lambda$.
At $N = 10{,}000$ and $T = 20$, 10,000 iterations take about 5.5 minutes
for the SAR probit, 7 minutes for the ordered model and 12.5 minutes for
the SEM model, including the sparse log-determinant grid. Impacts are computed after
sampling. With 200 draws at $N = 1{,}000$ they take 4.3 s, using the
stochastic diagonals. With 50 draws at $N = 10{,}000$ they take 19 s.

For ordered models, the cut-point move of Section 3.5 is what makes large
panels feasible. At $N = 400$, $T = 5$ it raises the ESS per second of the
cut-points 4-fold for $J = 3$, and 19- to 79-fold for the interior
cut-points with $J = 5$, compared with the Albert–Chib step alone.

# 6. Using the package

<!-- The output below comes from running the code with panelspatialprobit
0.0.0.9000 (2026-10-29) on simulated data. [TODO: rerun with the release
version; consider switching the example to the application data of
Section 7.] -->

This section walks through an analysis with the package, using data
simulated from a known model so that estimates can be compared with the
truth. Section 7 applies the same steps to real data.

## 6.1 Data and weights

The data are a balanced panel in long format: one row per unit and period,
with columns identifying both. The spatial weights carry the unit
identifiers as row names, and every unit in the data must match a row of
`W`; the order of the rows in the data does not matter. `psp_weights()`
accepts dense or sparse matrices and `spdep` neighbour (`nb`) or `listw`
objects, and `style = "W"` row-standardises. With polygons, a typical call
is `psp_weights(spdep::poly2nb(polygons), ids = polygons$id, style = "W")`.
A list with one matrix per period gives time-varying weights.

We use a 10 × 10 grid with rook contiguity, and simulate six periods from a
SAR probit with random effects:
$\beta = (-0.3, 1, -0.8)$, $\rho = 0.5$ and $\sigma_\alpha = 0.7$.

```r
library("panelspatialprobit")
W <- psp_weights(rook_weights(10, 10), style = "W")
d <- sim_psprobit(W, T = 6, beta = c(-0.3, 1, -0.8), rho = 0.5,
                  sigma_alpha = 0.7, seed = 42)
head(d, 3)
```

```
  id time y         x1         x2
1  1    1 1  1.3709584 -0.2484829
2  2    1 0 -0.5646982  0.4223204
3  3    1 0  0.3631284  0.9876533
```

## 6.2 Fitting a model

`psprobit()` fits binary models and `psoprobit()` ordered ones. Both take a
formula, the data, the names of the unit and time columns, and the weights.
The other arguments select the model of Section 2:

| Argument | Values | Section |
|---|---|---|
| `model` | `"sar"` (default), `"sdm"`, `"slx"`, `"sem"`, `"sdem"`, `"sarar"`, `"none"` | 2.2 |
| `effects` | `"random"` (default), `"pooled"`, `"spatial_random"` | 2.3 |
| `cre` | a one-sided formula (Mundlak), or a list with `type = "chamberlain"` | 2.3 |
| `durbin` | a formula restricting the Durbin terms to some regressors | 2.2 |
| `time_effects` | `TRUE` adds period dummies | 2.3 |
| `M` | the error weights, by default `W` | 2.2 |
| `dynamic` | `psp_dynamic(lag_y, lag_Wy, initial, latent)` | 2.4 |
| `errors` | `"normal"` or `psp_t(df)` | 2.5 |
| `priors` | `psp_priors(...)` with `prior_normal()`, `prior_beta()`, ... | 2.6 |
| `control` | `psp_control(iter, warmup, thin, chains, cores, seed, ...)` | 3 |

```r
fit <- psprobit(y ~ x1 + x2, data = d, index = c("id", "time"), W = W,
                model = "sar", effects = "random",
                control = psp_control(iter = 3000, warmup = 1000,
                                      chains = 4, cores = 4, seed = 1))
```

Chains run in parallel on `cores` cores. Each has its own random number
stream derived from `seed`, so the results do not depend on the number of
cores.

## 6.3 Output and convergence

```r
summary(fit)
```

```
Bayesian panel spatial probit (model = sar, effects = random)
Units: 100, periods: 6, observations: 600
Chains: 4, posterior draws: 12000

Posterior summary:
              mean    sd   q2.5  q97.5  rhat ess_bulk ess_tail
(Intercept) -0.285 0.096 -0.473 -0.101 1.002     1462     2683
x1           0.843 0.094  0.663  1.031 1.004     1452     2994
x2          -0.762 0.082 -0.925 -0.608 1.001     1828     3401
rho          0.468 0.063  0.336  0.583 1.002     1595     3091
sigma_alpha  0.728 0.124  0.500  0.991 1.007      599     1328

In-sample classification: hit rate 0.763, AUC 0.837
```

[The output is abridged: the call, the outcome distribution and the
impacts, shown in Section 6.4, are omitted.]

Every true value lies inside its 95% interval. The summary reports the
rank-normalised $\hat R$ and the bulk and tail ESS for each parameter
(Vehtari et al. 2021). The fit warns when $\hat R > 1.01$ or the bulk ESS
falls below 100 per chain. `plot(fit)` draws trace and density plots, and
`posterior::as_draws(fit)` returns the draws for packages such as
`bayesplot`. `coef()`, `vcov()`, `confint()`, `nobs()`, and the `broom`
methods `tidy()` and `glance()`, work as for other model objects.

As Section 8.1 explains, binary data carry little information per
observation, so the chains are autocorrelated. With rare outcomes, strong
dynamics or SARAR with $M = W$, the remedy for a convergence warning is
longer chains: increase `iter` and `warmup`, with `thin` to save memory.

## 6.4 Impacts

The impacts of Section 4.2 are computed after fitting and stored in
`fit$impacts`:

```
Impacts on P(y = 1) (posterior mean, sd, 95% interval; 500 draws)

direct:
      mean     sd    q2.5   q97.5
x1  0.1814 0.0152  0.1498  0.2083
x2 -0.1632 0.0142 -0.1913 -0.1364

indirect:
      mean     sd    q2.5   q97.5
x1  0.1425 0.0340  0.0844  0.2146
x2 -0.1276 0.0282 -0.1910 -0.0765

total:
      mean     sd    q2.5   q97.5
x1  0.3239 0.0410  0.2528  0.4047
x2 -0.2909 0.0331 -0.3619 -0.2320
```

A one-unit increase in $x_1$ for every unit raises the average probability
by 0.32. Of this, 0.18 is the unit's own effect (including feedback through
its neighbours) and 0.14 the spillover from its neighbours' changes. The
spillover is far less precisely estimated, because it depends on $\rho$.

`psp_impacts()` recomputes the impacts with other options:

- `type = "unit"` conditions on the sampled random effects;
- `by_period = TRUE` adds per-period impacts;
- `scale = "lesage_pace"` uses the convention of `spatialprobit`;
- `draws` sets the number of posterior draws used.

For ordered outcomes, impacts are reported for every category and sum to
zero across categories. The draws are in `$draws`, an array indexed by draw,
variable, category and effect, for derived quantities such as the share of
the total impact that is spillover.

## 6.5 Prediction and counterfactuals

`predict()` returns posterior mean probabilities, the latent index or the
most probable outcome, in the row order of the data. With `newdata`, it
evaluates a counterfactual: the same units and periods with changed
regressors. A discrete change, rather than the derivative behind the
impacts, is then a difference of two predictions:

```r
nd <- d
nd$x1 <- nd$x1 + 1
mean(predict(fit, newdata = nd) - predict(fit))
```

```
[1] 0.3260489
```

The average change in probability is close to the total impact of 0.32,
as expected, because the probit curve is nearly linear over a change of
this size at these probabilities. `psp_classification()` reports hit rates,
the AUC and a confusion matrix.

## 6.6 Dynamic models

`dynamic = psp_dynamic()` adds the lagged outcome, with Wooldridge initial
conditions by default (Section 2.4). We simulate eight periods with state
dependence $\gamma = 0.8$:

```r
dd <- sim_psprobit(W, T = 8, beta = c(-0.4, 1, -0.8), rho = 0.4,
                   sigma_alpha = 0.6, gamma = 0.8, seed = 3)
fd <- psprobit(y ~ x1 + x2, data = dd, index = c("id", "time"), W = W,
               dynamic = psp_dynamic(), control = psp_control(...))
round(coef(fd), 2)
```

```
(Intercept)          x1          x2       lag_y      init_y     mean_x1
      -0.47        1.09       -0.79        0.73        0.31       -0.07
    mean_x2         rho sigma_alpha
      -0.35        0.34        0.57
```

`lag_y` estimates the state dependence. `init_y` and the `mean_*` terms
describe how the random effect depends on the initial outcome and the
regressors; they are not causal effects and receive no impacts. In this
fit the sampler warned about the ESS of $\sigma_\alpha$, which the
Wooldridge terms leave little to explain. In an application the chains
should be run longer until the warning disappears.

`psp_response_path()` traces how a change works through the outcome
history (Section 4.3). Here $x_1$ is raised permanently by one unit for all
units from the first modelled period:

```r
psp_response_path(fd, "x1", type = "permanent", draws = 100, nsim = 10)
```

```
Response of P(y = 1) to a permanent shock of 1 in x1 (100 draws x 10 histories)

Average change over all units:
 horizon period   mean     sd   q2.5  q97.5
       0      2 0.2811 0.0262 0.2325 0.3322
       1      3 0.3261 0.0298 0.2799 0.3817
       2      4 0.3142 0.0279 0.2587 0.3704
       ...
       6      8 0.3840 0.0401 0.3132 0.4682
```

The short-run effect of 0.28 grows as the higher outcomes feed through
`lag_y`. It also varies from period to period, because the slope of the
probit depends on where each period's probabilities lie. Options:

- `type = "temporary"` gives a one-period change;
- `units =` changes only some units, and reports the effects on them and the
  spillovers to the others separately;
- `start =` sets the period of the change.

`psp_dynamic(latent = TRUE)` fits the latent AR(1) model instead.

## 6.7 Ordered outcomes

`psoprobit()` takes an ordered factor or integer codes and the same
arguments. The output adds the free cut-points `cut[2|3]`, ... (the first
is fixed at 0), and impacts for each category. For a four-category outcome
simulated with cut-points 0.8 and 1.6:

```r
fo <- psoprobit(y ~ x1 + x2, data = do, index = c("id", "time"), W = W,
                control = psp_control(...))
round(coef(fo), 2)
```

```
(Intercept)          x1          x2         rho sigma_alpha    cut[2|3]
       0.35        0.89       -0.81        0.40        0.62        0.79
   cut[3|4]
       1.66
```

The total impacts of $x_1$ are −0.31, −0.01, 0.05 and 0.27 on the four
category probabilities: raising $x_1$ moves probability from the lowest
category to the highest, and the four impacts sum to zero.

## 6.8 Model checking and comparison

`psp_ppc()` compares the data with replicated data sets (Section 4.4) on
three statistics, period by period: the outcome shares, Moran's I and
persistence. An excerpt:

```
persistence:
            observed rep_mean rep_q2.5 rep_q97.5   ppp
persist t=2     0.68    0.573     0.48      0.67 0.025
persist t=3     0.49    0.551     0.45      0.65 0.915
...
```

Posterior predictive p-values near 0 or 1 signal misfit. With 17
statistics, an occasional extreme value is expected by chance (here 2 of
17 are at or below 0.025, for a model that is correct by construction).
What matters is a systematic pattern across periods. `vignette("dynamics")`
shows such a pattern: a pooled static model underpredicts persistence in
most periods.

`psp_dic()` and `psp_loglik()` give the approximate DIC and log-likelihood
for comparing specifications. For random-effects models they rest on a
single $n$-dimensional GHK problem. On this example ($n = 600$) the
estimate was already noisy enough for the package to warn, so comparisons
of random-effects models need more GHK replications (`reps`) and should be
read with care (Section 8.2).

## 6.9 Practical advice

- **Start simple.** Fit a SAR model with random effects before adding
  Durbin terms, spatial errors or dynamics, and compare the impacts rather
  than the coefficients across specifications.
- **Check convergence first.** Treat any convergence warning as a reason to
  run longer chains, not as noise. Expect to need longer chains with:
  - rare outcomes (roughly under 10% in a category);
  - near-separation;
  - latent dynamics;
  - SARAR with $M = W$;
  - strongly autocorrelated spatial random effects.
- **Check prior sensitivity for $\rho$.** With binary outcomes and
  moderate $N$, refit with a different `prior_beta()` for $\rho$ and
  compare.
- **Handle the initial condition properly.** With lagged outcomes, keep the
  Wooldridge initial condition unless there is a reason not to. Pooled
  dynamic models confound heterogeneity with state dependence
  (Section 5.4).
- **Use impacts for interpretation, not coefficients.** For ordered
  outcomes, read the category impacts together, because they move
  probability between categories.

<!-- Section 7 (Empirical application) is not drafted yet; see the
outline. -->

# 8. Limitations and extensions

## 8.1 Mixing

The sampler mixes well when the outcomes are balanced, and slowly in four
situations:

- rare outcomes (roughly under 10% in one category);
- near-separation (outcomes that barely vary within units or periods);
- strong latent dynamics;
- SARAR models with $M = W$.

In these cases the intercept, $\rho$ and $\sigma_\alpha$ are strongly
correlated a posteriori, and their chains are highly autocorrelated. In a
test case with 4.6% ones ($N = 100$, $T = 5$), the posterior correlations
were 0.84 between the intercept and $\rho$, and −0.69 between the intercept
and $\sigma_\alpha$. The autocorrelation times were 160 for the intercept
and 174 for $\sigma_\alpha$, against 39 for a slope. Even the
best-identified combination, the mean latent level
$\text{intercept}/(1 - \rho)$, had an ESS of only 135 out of 16,000 draws.
So the slow direction is not a single pair of parameters that a
reparameterisation could separate. It is the known slowness of data
augmentation with imbalanced categorical outcomes (Johndrow et al. 2019).
Given the latent values of the many zeros, the intercept and the spatial
level are tightly pinned down, while given the outcomes alone they are not.

We tried three further moves, besides the scale move of Section 3.6. Each
passed the comparison with the reference sampler (Section 5.1) and none
improved the ESS, so all three were removed:

- **Scale move including $\sigma_\alpha$:** $(z, \beta, \alpha, \kappa,
  \sigma_\alpha)$ rescaled jointly. The unit error variance pins the overall
  latent scale, so the scale factor stays within about $1 \pm 1/\sqrt{2n}$
  and the move does almost nothing.
- **Interweaving for $\sigma_\alpha$** (Yu & Meng 2011): $\sigma_\alpha$
  redrawn given the standardised effects $\alpha/\sigma_\alpha$.
- **Time-profile move for latent dynamics:** the generalised Gibbs map
  $z_t \mapsto g^{t-1}z_t$, $\gamma \mapsto g\gamma$.

Running four latent sweeps per iteration did not help either. That result
also rules out a joint update of each unit's latent path $z_{i,1..T}$,
however it is implemented. The best such a move can achieve is an exact draw
of $z$ given the parameters, which repeated sweeps already approximate. The
slow direction couples $z$ with the parameters, and no move on $z$ alone can
remove it.

The remedies in the literature work on the marginal posterior of the
parameters, either through calibrated data augmentation (Duan, Johndrow &
Dunson 2018) or through Metropolis or Hamiltonian moves. Both need the
likelihood $p(y \mid \theta)$, which here is an $n$-dimensional rectangle
probability. Only a noisy GHK estimate of it is available (Section 4.4),
too noisy to drive an exact sampler without pseudo-marginal machinery. For
now the package does three things:

- it documents the four situations;
- it warns when the split-$\hat R$ exceeds 1.01 or the bulk ESS falls
  below its threshold;
- it supports long thinned chains without excessive memory use.

The SBC (Section 5.2) and the Wang & Kockelman replication (Section 5.5)
show that, run long enough, the sampler is calibrated in these cases too.

## 8.2 Approximations

The posterior draws are exact MCMC: no step approximates its target beyond
floating-point and grid resolution. A few quantities are approximated:

- **Grid draws.** $\rho$, $\lambda$ and $\phi$ are drawn on two-stage grids
  with continuous interpolation. For $N$ above about 3,000 the
  log-determinant comes from a spline through sparse LU evaluations, with
  an interpolation error that is negligible next to the posterior spread.
- **Impacts and probabilities for large $N$.** Above `impacts_exact_max`
  (500 units by default), the diagonals in (5) are estimated with
  Rademacher probes. At $N = 1{,}000$ the error was about 5% of the
  posterior standard deviation of the impacts with the default 50 probes.
  With $t$ errors the probabilities average over simulated mixing variances,
  and with latent dynamics for large $N$ the variance recursions use
  probes. In each case the number of simulations can be raised.
- **Response paths** with lagged outcomes simulate outcome histories, so
  they carry simulation error that falls with the number of histories.
- **The likelihood** is a GHK estimate. It is exact in the non-spatial case
  and precise when periods can be simulated separately. For random-effects
  models it is a single $n$-dimensional problem and becomes noisy beyond a
  few hundred observations. It is not available with $t$ errors. DIC and
  likelihood comparisons are therefore labelled approximate.

## 8.3 Restrictions of the implemented models

- **Balanced panels.** Every unit in $W$ must appear in every period with
  complete regressors. Missing outcomes are allowed in static models, which
  covers panels where only the outcome is unobserved. Units that enter or
  leave the sample, or have missing regressors, are not supported.
- **Dynamics.**
  - Missing outcomes are not allowed in dynamic models.
  - Latent dynamics need time-invariant weights and cannot be combined with
    lagged outcomes.
  - The first period of the latent AR(1) model has no lagged term. It is
    not drawn from the stationary distribution of the process, which
    matters mainly for short panels with strong persistence.
  - The Wooldridge initial condition assumes that the random effect depends
    linearly on the initial outcome and on the unit means of the
    regressors.
  - With lagged outcomes, `predict(newdata = )` is not available, because
    new regressors would change the outcome history; response paths take
    its place.
- **Spatial random effects** need time-invariant weights.
- **Student-$t$ errors** have fixed degrees of freedom. Estimating $\nu$
  would compete with the probit scale normalisation.
- **Ordered models** require an intercept, because the first cut-point is
  fixed at zero.
- **SARAR with $M = W$** identifies $\rho$ and $\lambda$ only through the
  regressors. The package warns, but it cannot repair weak identification.
- **Binary outcomes carry little information about $\rho$.** With
  unbalanced outcomes the posterior of $\rho$ can depend on its prior, and a
  prior-sensitivity check is advisable.

## 8.4 Out of scope

Some models were left out deliberately:

- **Fixed effects**, for the reasons in Section 2.3. A large-$T$
  frequentist alternative exists (Arbia, Billé & Leorato 2023).
- **Bivariate, multivariate and multinomial outcomes**, such as the system
  of Baltagi, Egger & Kesina (2016).
- **Spatial Tobit models.**
- **Convex combinations of several weight matrices** (Debarsy & LeSage
  2018) and estimated or endogenous weights.
- **Frequentist estimators**: approximate maximum likelihood
  (`ProbitSpatial`), partial maximum likelihood and GMM.

## 8.5 Limits of the validation

- The SBC uses a small grid ($N = 16$, $T = 4$) so that thousands of fits
  are affordable. It checks the correctness of the sampler, which does not
  depend on $N$, but not its mixing at realistic sizes. That is covered by
  the Monte Carlo study and the benchmarks.
- The Monte Carlo study covers SAR designs only, with three static designs
  and one dynamic one. It does not cover SEM, SDM, SARAR, $t$ errors or
  time-varying weights; their correctness rests on SBC.
- The Wang & Kockelman replication had to adapt their design (Section 5.5).
  It was set up from the 2008 working paper, not the published article.
- The cause of `spatialprobit`'s miscalibration (Section 5.3) has not been
  isolated.

## 8.6 Extensions

- **Faster mixing** for rare outcomes and strong dynamics, through
  calibrated data augmentation or moves on the marginal posterior, possibly
  with a pseudo-marginal GHK likelihood.
- **Unbalanced panels**, with units entering and leaving. This needs weights
  restricted to the units present in each period and a matching
  log-determinant per period.
- **Bivariate outcomes**, which fit the $Az - b$ engine with a $2 \times 2$
  error covariance per observation.
- **A stationary initial condition** for the latent AR(1) model.
- **Parallel sweeps within a chain**, using a colouring of the graph of
  $A$, for very large $N$. This needs care to keep runs reproducible.

# 9. Conclusion

Short panels of binary and ordered outcomes with spatial dependence are
common, but applied researchers have had to choose between spatial
cross-section software and non-spatial panel models. `panelspatialprobit`
removes that choice. It estimates spatial probit and ordered probit models
with unit effects, dynamics and the usual spatial specifications in one
interface, and reports their results as impacts on probabilities rather
than coefficients.

Three points emerge from building and testing the package.

- **One kernel is enough.** Writing every model as the sparse kernel
  $\|Az - b\|^2$ keeps the sampler small: the spatial specification, the
  error process, the random effects and the latent dynamics change only
  the operator $A$ and the offset $b$. The latent update costs
  $O(\operatorname{nnz}(A))$, so a panel of $N = 10{,}000$ units and
  $T = 20$ periods takes minutes (Section 5.7). New variants enter by
  defining $A$, not by writing new conditionals.
- **Small design choices decide correctness.** The order of the blocks
  around the collapsed draw of $\rho$, the kernel used for $\rho$ and the
  exactness of the truncated normal draws are easy to get wrong, and none
  of them produces an error message. The same calibration test that our
  sampler passes shows that `spatialprobit`'s SAR probit sampler fails on
  about half of the simulated data sets and understates $\rho$ on the rest
  (Section 5.3). Making simulation-based calibration a condition for
  documenting a variant, rather than a one-off check, is what makes the
  breadth of the package trustworthy: 30 variants, 199 rank tests, with no
  more departures than chance and slow mixing explain (Section 5.2).
- **Ignoring the panel structure is costly.** In the Monte Carlo study, a
  pooled dynamic model overstated state dependence by about 40% and
  covered the true value in only 11% of samples (Section 5.4). Random
  effects with Wooldridge initial conditions removed the bias.

[TODO: two or three sentences on what the empirical application of
Section 7 shows, e.g. how state dependence and spatial diffusion compare
once heterogeneity is modelled.]

The main open problem is mixing when outcomes are rare or nearly separated,
when latent dynamics are strong, and in SARAR models with $M = W$
(Section 8.1). The slow direction couples the latent values with the
parameters, so no further move on the latent values alone can remove it. A
real gain will need moves on the marginal posterior of the parameters,
perhaps driven by a pseudo-marginal version of the GHK likelihood already
in the package. Until then, the package warns when chains have not
converged, and the calibration results show that long enough chains give
correct answers. Unbalanced panels and bivariate outcomes fit the same
engine and are the natural next extensions (Section 8.6).

## Computational details

The results in this paper were obtained with R 4.6.0 (R Core Team 2026)
and the development version 0.0.0.9000 of `panelspatialprobit` [TODO: rerun
with the release version and update], together with `posterior` 1.7.0
(Bürkner et al. 2026), `Matrix` 1.7-5 (Bates, Maechler & Jagan 2026),
`Rcpp` 1.1.1-1.1 (Eddelbuettel & François 2011), `RcppEigen` 0.3.4.0.2
(Bates & Eddelbuettel 2013), `spatialprobit` 1.0.4 and `MCMCpack` 1.7.1. The
computations ran between 27 and 29 October 2026 on a desktop computer with an
Intel Core i7-4790K processor (4 cores, 8 threads, 4.0 GHz) and 16 GB of
memory, under Pop!_OS 24.04 LTS (Linux). Timings are for a single core
unless stated otherwise; the simulation studies ran their replications in
parallel on up to 7 cores. Every fit sets its seed through
`psp_control(seed = )`, and results do not depend on the number of cores.
The package is available from CRAN [TODO: once accepted] and from
[TODO: repository URL]. The scripts that reproduce the validation results
are in the `validation/` directory of the source repository, and
`vignette("validation")` reproduces the figures and tables of Section 5.
[TODO: replication materials as required by the venue.]

# References (partial; check every entry)

Albert, J. H. and Chib, S. (1993). Bayesian analysis of binary and
polychotomous response data. *Journal of the American Statistical
Association* 88, 669–679.

Baltagi, B. H., Egger, P. H. and Kesina, M. (2016). Bayesian spatial
bivariate panel probit estimation. *Advances in Econometrics* 37, 119–144.

Baltagi, B. H., Egger, P. H. and Kesina, M. (2018). Generalized spatial
autocorrelation in a panel-probit model with an application to exporting in
China. *Empirical Economics* 55, 193–211. [TODO: check the random-effects
specification and priors against this paper, which we have not read.]

Bekas, C., Kokiopoulou, E. and Saad, Y. (2007). An estimator for the
diagonal of a matrix. *Applied Numerical Mathematics* 57, 1214–1229.

Chamberlain, G. (1982). Multivariate regression models for panel data.
*Journal of Econometrics* 18, 5–46.

Cowles, M. K. (1996). Accelerating Monte Carlo Markov chain convergence for
cumulative-link generalized linear models. *Statistics and Computing* 6,
101–111.

Elhorst, J. P., Heijnen, P., Samarina, A. and Jacobs, J. P. A. M. (2017).
Transitions at different moments in time: a spatial probit approach.
*Journal of Applied Econometrics* 32, 422–439.

Gelman, A. (2006). Prior distributions for variance parameters in
hierarchical models. *Bayesian Analysis* 1, 515–534.

Imai, K. and van Dyk, D. A. (2005). A Bayesian analysis of the multinomial
probit model using marginal data augmentation. *Journal of Econometrics*
124, 311–334.

Johndrow, J. E., Smith, A., Pillai, N. and Dunson, D. B. (2019). MCMC for
imbalanced categorical data. *Journal of the American Statistical
Association* 114, 1394–1403.

Kakamu, K. and Wago, H. (2005). Bayesian spatial panel probit model with an
application to business cycle in Japan. In A. Zerger and R. M. Argent
(eds), *MODSIM 2005 International Congress on Modelling and Simulation*,
856–863. Modelling and Simulation Society of Australia and New Zealand.

LeSage, J. P. (2000). Bayesian estimation of limited dependent variable
spatial autoregressive models. *Geographical Analysis* 32, 19–35.

LeSage, J. and Pace, R. K. (2009). *Introduction to Spatial Econometrics*.
CRC Press.

Liu, J. S. and Sabatti, C. (2000). Generalised Gibbs sampler and multigrid
Monte Carlo for Bayesian computation. *Biometrika* 87, 353–369.

Liu, J. S. and Wu, Y. N. (1999). Parameter expansion for data augmentation.
*Journal of the American Statistical Association* 94, 1264–1274.

Martinetti, D. and Geniaux, G. (2017). Approximate likelihood estimation of
spatial probit models. *Regional Science and Urban Economics* 64, 30–45.

Modrák, M., Moon, A. H., Kim, S., Bürkner, P., Huurre, N., Faltejsková, K.,
Gelman, A. and Vehtari, A. (2025). Simulation-based calibration checking
for Bayesian computation: the choice of test quantities shapes sensitivity.
*Bayesian Analysis*. doi:10.1214/23-BA1404.

Mundlak, Y. (1978). On the pooling of time series and cross section data.
*Econometrica* 46, 69–85.

Neyman, J. and Scott, E. L. (1948). Consistent estimates based on partially
consistent observations. *Econometrica* 16, 1–32.

Pace, R. K. and Barry, R. (1997). Quick computation of spatial
autoregressive estimators. *Geographical Analysis* 29, 232–247.

Robert, C. P. (1995). Simulation of truncated normal variables. *Statistics
and Computing* 5, 121–125.

Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A.
(2002). Bayesian measures of model complexity and fit. *Journal of the
Royal Statistical Society B* 64, 583–639.

Talts, S., Betancourt, M., Simpson, D., Vehtari, A. and Gelman, A. (2018).
Validating Bayesian inference algorithms with simulation-based
calibration. arXiv:1804.06788.

van Dyk, D. A. and Park, T. (2008). Partially collapsed Gibbs samplers:
theory and methods. *Journal of the American Statistical Association* 103,
790–796.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Bürkner, P.-C.
(2021). Rank-normalization, folding, and localization: an improved
$\hat R$ for assessing convergence of MCMC. *Bayesian Analysis* 16,
667–718.

Wand, M. P., Ormerod, J. T., Padoan, S. A. and Frühwirth, R. (2011). Mean
field variational Bayes for elaborate distributions. *Bayesian Analysis*
6, 847–900.

Wang, X. and Kockelman, K. M. (2009). Bayesian inference for ordered
response data with a dynamic spatial-ordered probit model. *Journal of
Regional Science* 49, 877–913.

Wilhelm, S. and de Matos, M. G. (2013). Estimating spatial probit models in
R. *The R Journal* 5(1), 130–143.

Wooldridge, J. M. (2005). Simple solutions to the initial conditions problem
in dynamic, nonlinear panel data models with unobserved heterogeneity.
*Journal of Applied Econometrics* 20, 39–54.

Wooldridge, J. M. (2010). *Econometric Analysis of Cross Section and Panel
Data*, 2nd ed. MIT Press.

LeSage, J. P., Pace, R. K., Lam, N., Campanella, R. and Liu, X. (2011).
New Orleans business recovery in the aftermath of Hurricane Katrina.
*Journal of the Royal Statistical Society A* 174, 1007–1027.

Martin, A. D., Quinn, K. M. and Park, J. H. (2011). MCMCpack: Markov chain
Monte Carlo in R. *Journal of Statistical Software* 42(9), 1–21.

Arbia, G., Billé, A. G. and Leorato, S. (2023). Concentrated partial ML
estimation for dynamic spatial panel data probit models with fixed effects.
SSRN working paper 4514415.

Billé, A. G. and Leorato, S. (2020). Partial ML estimation for spatial
autoregressive nonlinear probit models with autoregressive disturbances.
*Econometric Reviews* 39, 437–475.

Klier, T. and McMillen, D. P. (2008). Clustering of auto supplier plants in
the United States: generalized method of moments spatial logit for large
samples. *Journal of Business and Economic Statistics* 26, 460–471.

Piras, G. and Sarrias, M. (2023). GMM estimators for binary spatial models
in R. *Journal of Statistical Software* 107(8).

Spinelli, D. (2022). Fitting spatial autoregressive logit and probit models
using Stata: the spatbinary command. *Stata Journal* 22, 293–318.

Debarsy, N. and LeSage, J. P. (2018). Flexible dependence modeling using
convex combinations of different types of connectivity structures.
*Regional Science and Urban Economics* 69, 48–68.

Duan, L. L., Johndrow, J. E. and Dunson, D. B. (2018). Scaling up data
augmentation MCMC via calibration. *Journal of Machine Learning Research*
19(64), 1–34.

Yu, Y. and Meng, X.-L. (2011). To center or not to center: that is not the
question. An ancillarity–sufficiency interweaving strategy (ASIS) for
boosting MCMC efficiency. *Journal of Computational and Graphical
Statistics* 20, 531–570.

Geweke, J. (1991). Efficient simulation from the multivariate normal and
Student-$t$ distributions subject to linear constraints and the evaluation
of constraint probabilities. In E. M. Keramidas (ed.), *Computing Science
and Statistics: Proceedings of the 23rd Symposium on the Interface*,
571–578. Interface Foundation of North America. [TODO: check the pages.]

Geweke, J. (1993). Bayesian treatment of the independent Student-$t$
linear model. *Journal of Applied Econometrics* 8, S19–S40.

Hajivassiliou, V. A. and McFadden, D. L. (1998). The method of simulated
scores for the estimation of LDV models. *Econometrica* 66, 863–896.

Keane, M. P. (1994). A computationally practical simulation estimator for
panel data. *Econometrica* 62, 95–116.

Liu, C. (2004). Robit regression: a simple robust alternative to logistic
and probit regression. In A. Gelman and X.-L. Meng (eds), *Applied Bayesian
Modeling and Causal Inference from Incomplete-Data Perspectives*, 227–238.
Wiley.

Software:

R Core Team (2026). *R: A Language and Environment for Statistical
Computing*. R Foundation for Statistical Computing, Vienna.

Bates, D. and Eddelbuettel, D. (2013). Fast and elegant numerical linear
algebra using the RcppEigen package. *Journal of Statistical Software*
52(5), 1–24.

Bates, D., Maechler, M. and Jagan, M. (2026). *Matrix: Sparse and Dense
Matrix Classes and Methods*. R package version 1.7-5.

Bürkner, P.-C., Gabry, J., Kay, M. and Vehtari, A. (2026). *posterior:
Tools for Working with Posterior Distributions*. R package version 1.7.0.

Eddelbuettel, D. and François, R. (2011). Rcpp: seamless R and C++
integration. *Journal of Statistical Software* 40(8), 1–18.
