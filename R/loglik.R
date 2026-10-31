#' Approximate log-likelihood and DIC (experimental)
#'
#' The likelihood of a spatial probit model is an \eqn{NT}-dimensional normal
#' rectangle probability with no closed form. `psp_loglik()` estimates it with
#' the GHK simulator (Geweke, Hajivassiliou and Keane), applied to the sparse
#' precision matrix of the latent outcomes so that it scales to large \eqn{N}.
#' For random-effects models the random effects are added to the latent
#' vector and integrated out with it (`type = "marginal"`), or held at their
#' posterior draws (`type = "conditional"`).
#'
#' The GHK estimate of the likelihood is unbiased, but its logarithm is biased
#' downwards when the replications vary a lot; increase `reps` until the
#' estimates settle, and compare models with the same `reps` and `seed`
#' (common random numbers). `attr(, "mcse")` gives a delta-method Monte Carlo
#' standard error for each value, and a warning is given when it exceeds 0.5.
#' Without random effects to integrate out (pooled models, or
#' `type = "conditional"`) the periods are independent and are simulated
#' separately, which keeps the error small even for large panels. The marginal
#' likelihood of a random-effects model is a single \eqn{NT}-dimensional
#' problem and becomes unreliable beyond a few hundred observations.
#'
#' `logLik()` evaluates the marginal log-likelihood at the posterior mean of
#' the parameters, so `AIC()` and `BIC()` work. `psp_dic()` computes the
#' deviance information criterion (Spiegelhalter et al. 2002) from the
#' marginal likelihood.
#'
#' The posterior draws are processed in parallel under the user's
#' [future::plan()] (see "Parallel post-processing" in [psp_control()]).
#'
#' @param fit,object A fitted `psprobit` object.
#' @param draws Number of posterior draws to evaluate.
#' @param reps Number of GHK replications per evaluation.
#' @param type `"marginal"` (random effects integrated out) or
#'   `"conditional"` (given the random-effect draws kept by the sampler).
#'   The two coincide for pooled models.
#' @param seed Seed for the GHK replications (the same random numbers are
#'   used for every draw).
#' @param ... Unused.
#' @return `psp_loglik()`: a vector of log-likelihood values, one per draw,
#'   with attribute `mcse`. `logLik()`: a `logLik` object. `psp_dic()`: a
#'   list with `dic`, `p_d` (effective number of parameters), `d_bar`
#'   (posterior mean deviance) and `d_hat` (deviance at the posterior mean).
#' @references
#' Geweke, J. (1991). Efficient simulation from the multivariate normal and
#' Student-t distributions subject to linear constraints. *Computer Science and
#' Statistics: Proceedings of the 23rd Symposium on the Interface*, 571-578.
#'
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A.
#' (2002). Bayesian measures of model complexity and fit. *JRSS B* 64,
#' 583-639.
#' @examples
#' W <- psp_weights(rook_weights(4, 4), style = "W")
#' d <- sim_psprobit(W, T = 3, beta = c(0, 1), rho = 0.4, seed = 2)
#' fit <- psprobit(y ~ x1, data = d, index = c("id", "time"), W = W, impacts = FALSE,
#'                 control = psp_control(iter = 200, warmup = 100, chains = 1, seed = 2))
#' logLik(fit, reps = 200)
#' psp_dic(fit, draws = 20, reps = 100)
#' @export
psp_loglik <- function(fit, draws = 100, reps = 200, type = c("marginal", "conditional"),
                       seed = 1) {
  if (!inherits(fit, "psprobit")) stop("'fit' must be a psprobit object.", call. = FALSE)
  type <- match.arg(type)
  .psp_ghk_normal_only(fit)
  conditional <- type == "conditional" && fit$effects != "pooled"
  if (conditional) {
    keep <- .psp_draw_subset_from(seq_along(fit$alpha_index), draws)
    idx <- fit$alpha_index[keep]
    alpha <- fit$alpha_draws[keep, , drop = FALSE]
  } else {
    idx <- .psp_draw_subset(fit, draws)
  }
  pm <- .psp_par_matrix(fit)[idx, , drop = FALSE]
  out <- .psp_ghk_loglik(fit, pm, reps, seed, alpha = if (conditional) alpha)
  .psp_ghk_check(out$mcse)
  structure(out$ll, mcse = out$mcse, type = type, reps = reps)
}

#' @rdname psp_loglik
#' @export
logLik.psprobit <- function(object, reps = 500, seed = 1, ...) {
  .psp_ghk_normal_only(object)
  par <- object$coefficients
  g <- .psp_ghk_loglik(object, matrix(par, 1, dimnames = list(NULL, names(par))), reps, seed)
  .psp_ghk_check(g$mcse)
  structure(g$ll, df = length(par), nobs = nobs.psprobit(object), mcse = g$mcse,
            class = "logLik")
}

#' @rdname psp_loglik
#' @export
psp_dic <- function(fit, draws = 100, reps = 200, seed = 1) {
  ll <- psp_loglik(fit, draws, reps, "marginal", seed)
  d_bar <- -2 * mean(ll)
  par <- fit$coefficients
  d_hat <- -2 * .psp_ghk_loglik(fit, matrix(par, 1, dimnames = list(NULL, names(par))), reps,
                                seed)$ll
  if (d_bar - d_hat < 0) {
    warning("Negative p_D: the GHK error dominates; increase 'reps'.", call. = FALSE)
  }
  p_d <- d_bar - d_hat
  list(dic = d_bar + p_d, p_d = p_d, d_bar = d_bar, d_hat = d_hat)
}

# GHK log-likelihood for each row of the draws x parameters matrix 'pm'
# (pspcore::ghk; list(ll, mcse)). The latent vector is (z_1, ..., z_T) and,
# for the marginal likelihood of a random-effects model, alpha. With
# A_t = R_t S_t, the centred kernel is
#   sum_t |A_t z_t - R_t alpha|^2 + |alpha|^2 / sigma^2,
# whose precision is sparse (with |B alpha|^2 for spatial random effects);
# given alpha (draws x N, conditional) only the z blocks remain. The
# unrestricted coordinates (alpha, missing outcomes) are ordered first, so
# GHK integrates them out exactly; independent periods are simulated as
# separate blocks. Every draw uses the same random numbers, from 'seed'.
.psp_ghk_loglik <- function(fit, pm, reps, seed, alpha = NULL) {
  parts <- .psp_post_map(psp_rs_ghk, .psp_rs_post_model(fit, NULL, NULL),
                         .psp_rs_draws(fit, pm, alpha), y = as.integer(fit$y),
                         reps = as.integer(reps), seed = .psp_rs_seeds(seed, 1))
  list(ll = unlist(lapply(parts, `[[`, "ll")), mcse = unlist(lapply(parts, `[[`, "mcse")))
}

.psp_ghk_normal_only <- function(fit) {
  if (!is.null(fit$errors) && fit$errors$type != "normal") {
    stop("The approximate likelihood is only available for normal errors.", call. = FALSE)
  }
}

.psp_ghk_check <- function(mcse) {
  if (any(mcse > 0.5)) {
    warning(sprintf(paste("The GHK estimates are unreliable (Monte Carlo se up to %.2f on the",
                          "log scale): increase 'reps'. GHK needs many replications for",
                          "large N x T."), max(mcse)), call. = FALSE)
  }
  invisible(mcse)
}
