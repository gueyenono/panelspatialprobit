#' Student-t latent errors
#'
#' Replaces the normal errors \eqn{\varepsilon_{it} \sim N(0, 1)} of
#' [psprobit()] / [psoprobit()] by Student-t errors with `df` degrees of
#' freedom, which makes the model robust to outlying observations (for binary
#' outcomes this is the robit model of Liu 2004; `df = 7` approximates the
#' logit link). The errors are sampled as the scale mixture
#' \eqn{\varepsilon_{it} \sim N(0, v_{it})}, \eqn{v_{it} \sim IG(\nu/2,
#' \nu/2)}; the posterior means of \eqn{v_{it}} are returned in `fit$v_mean`
#' (large values flag observations the model fits poorly).
#'
#' The degrees of freedom are fixed rather than estimated: with discrete
#' outcomes the data carry very little information about them.
#'
#' With t errors the latent outcome is no longer normal, because the spatial
#' filter mixes the errors of several units. Probabilities and impacts are
#' therefore averaged over `control$mixture_draws` simulated sets of
#' \eqn{v} (fixed per fit, so results are reproducible): given \eqn{v} the
#' latent outcome is normal and the usual formulas apply. Above
#' `impacts_exact_max` units the reduced-form variance treats the own error
#' term exactly and the neighbours' terms at their mean.
#'
#' @param df Degrees of freedom, larger than 2.
#' @return An object of class `psp_errors`.
#' @references Liu, C. (2004). Robit regression: a simple robust alternative
#'   to logistic and probit regression. In Gelman, A. and Meng, X.-L. (eds.),
#'   *Applied Bayesian Modeling and Causal Inference from Incomplete-Data
#'   Perspectives*, 227-238. Wiley.
#' @examples
#' W <- psp_weights(rook_weights(5, 5), style = "W")
#' d <- sim_psprobit(W, T = 4, beta = c(0, 1, -1), rho = 0.4, df = 5, seed = 1)
#' fit <- psprobit(y ~ x1 + x2, d, c("id", "time"), W, errors = psp_t(df = 5),
#'                 impacts = FALSE,
#'                 control = psp_control(iter = 300, warmup = 200, chains = 1, seed = 1))
#' coef(fit)
#' @export
psp_t <- function(df = 7) {
  if (length(df) != 1 || !is.finite(df) || df <= 2) {
    stop("'df' must be a single number larger than 2.", call. = FALSE)
  }
  structure(list(type = "t", df = df), class = "psp_errors")
}

.psp_check_errors <- function(errors) {
  if (is.null(errors) || identical(errors, "normal")) {
    return(structure(list(type = "normal", df = Inf), class = "psp_errors"))
  }
  if (identical(errors, "t")) return(psp_t())
  if (!inherits(errors, "psp_errors")) {
    stop("'errors' must be \"normal\" or come from psp_t().", call. = FALSE)
  }
  errors
}

# Mixing variances v ~ IG(df/2, df/2) for t errors: N x mixture_draws, fixed per fit
.psp_t_mix <- function(fit) {
  if (is.null(fit$errors) || fit$errors$type == "normal") return(NULL)
  df <- fit$errors$df
  .psp_with_rng_restore({
    set.seed(fit$seed + 1L)
    matrix(1 / stats::rgamma(fit$N * fit$control$mixture_draws, df / 2, df / 2), fit$N)
  })
}

# Standard draws of the latent errors (normal or t)
.psp_rerr <- function(fit, n) {
  if (is.null(fit$errors) || fit$errors$type == "normal") stats::rnorm(n)
  else stats::rt(n, fit$errors$df)
}

# Mean over the mixture columns of f(s): s is a vector (normal errors) or a
# matrix with one column per mixture draw
.psp_mix_mean <- function(f, s) if (is.matrix(s)) rowMeans(f(s)) else f(s)
