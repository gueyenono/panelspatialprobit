#' Dynamic specification
#'
#' Adds observed state dependence to [psprobit()] / [psoprobit()]: the
#' previous period's outcome \eqn{y_{i,t-1}} and, optionally, its spatial lag
#' \eqn{(W y_{t-1})_i} enter the latent equation as regressors. The first
#' period is used only as the initial condition, so estimation uses periods
#' \eqn{2, \dots, T}. For ordered outcomes the lag enters as one indicator per
#' category above the first.
#'
#' With `initial = "wooldridge"` (the default) the random effects are modelled
#' as \eqn{\alpha_i = \psi' y_{i1} + \bar x_i'\pi + a_i} (Wooldridge 2005): the
#' first-period outcome (columns `init_*`) and the unit means of all regressors
#' over the estimation periods (columns `mean_*`) are added to the model. This
#' handles the initial-conditions problem of dynamic random-effects models.
#'
#' Impacts of the regressors are short-run (within-period) impacts.
#'
#' With `latent = TRUE` the dynamics act on the latent outcome instead,
#' \eqn{z_t = \rho W z_t + \gamma z_{t-1} + X_t\beta + \alpha + u_t}
#' (Wang and Kockelman 2009). The coefficient \eqn{\gamma} is reported as
#' `lag_z`. Every period is modelled; the first has no lagged term, so its
#' latent outcome serves as the initial condition and shares the random effect
#' with the later periods (no separate initial-conditions device is needed).
#' Probabilities, impacts and predictions integrate over the latent history,
#' and [psp_response_path()] gives the propagation of a shock. Latent dynamics
#' need time-invariant weights and cannot be combined with lagged outcomes.
#' \eqn{\gamma} mixes slowly against \eqn{\rho} (it is drawn given the latent
#' values), so plan on longer chains than for static models and check the
#' convergence diagnostics.
#'
#' @param lag_y Logical; include \eqn{y_{i,t-1}} (columns `lag_*`).
#' @param lag_Wy Logical; include \eqn{(W y_{t-1})_i} (columns `W_lag_*`).
#' @param initial `"wooldridge"` or `"none"` (treat the first-period outcome as
#'   exogenous). Ignored with `latent = TRUE`.
#' @param latent Logical; latent AR(1) dynamics instead of lagged outcomes.
#' @return An object of class `psp_dynamic`.
#' @references Wang, X. and Kockelman, K. M. (2009). Bayesian inference for
#'   ordered response data with a dynamic spatial-ordered probit model.
#'   *Journal of Regional Science* 49, 877-913.
#'
#'   Wooldridge, J. M. (2005). Simple solutions to the initial
#'   conditions problem in dynamic, nonlinear panel data models with
#'   unobserved heterogeneity. *Journal of Applied Econometrics* 20, 39-54.
#' @examples
#' W <- psp_weights(rook_weights(5, 5), style = "W")
#' d <- sim_psprobit(W, T = 6, beta = c(-0.3, 1, -1), rho = 0.3, gamma = 0.8, seed = 1)
#' fit <- psprobit(y ~ x1 + x2, d, c("id", "time"), W, dynamic = psp_dynamic(),
#'                 impacts = FALSE,
#'                 control = psp_control(iter = 300, warmup = 200, chains = 1, seed = 1))
#' coef(fit)
#' @export
psp_dynamic <- function(lag_y = !latent, lag_Wy = FALSE, initial = c("wooldridge", "none"),
                        latent = FALSE) {
  initial <- match.arg(initial)
  if (isTRUE(latent)) {
    if (isTRUE(lag_y) || isTRUE(lag_Wy)) {
      stop("Latent dynamics cannot be combined with lagged outcomes.", call. = FALSE)
    }
    return(structure(list(lag_y = FALSE, lag_Wy = FALSE, initial = "none", latent = TRUE),
                     class = "psp_dynamic"))
  }
  if (!lag_y && !lag_Wy) stop("At least one of 'lag_y' and 'lag_Wy' must be TRUE.", call. = FALSE)
  structure(list(lag_y = isTRUE(lag_y), lag_Wy = isTRUE(lag_Wy), initial = initial,
                 latent = FALSE),
            class = "psp_dynamic")
}

# Does the fit / design have latent AR(1) dynamics?
.psp_latent_ar <- function(x) isTRUE(x$dynamic$latent)
