#' Impacts (marginal effects) on outcome probabilities
#'
#' Computes average direct, indirect (spillover) and total impacts of each
#' regressor on the outcome probabilities, following LeSage and Pace (2009,
#' ch. 10) but accounting for the reduced-form standard deviation of the
#' latent outcome. Impacts are population-averaged: the random effects are
#' integrated out. For ordered models impacts are computed per category (they
#' sum to zero over categories).
#'
#' With `type = "unit"` the impacts are instead conditional on the posterior
#' draws of the random effects: \eqn{m_t = S^{-1}(X_t\beta + \alpha)} and
#' \eqn{s_i^2 = [S^{-1}S^{-\top}]_{ii}}. They describe the sample units
#' rather than the population and use the random-effect draws kept by the
#' sampler (about `control$impacts_draws` of them).
#'
#' For regressor \eqn{r} and period \eqn{t},
#' \deqn{\partial P(y_{it} = 1) / \partial x_{jt,r} =
#'   \phi(m_{it}/s_i)/s_i \, [S^{-1}(\beta_r I + \theta_r W)]_{ij},}
#' where \eqn{m_t = S^{-1}X_t\beta} and \eqn{s_i^2 = (1 + \sigma_\alpha^2)
#' [S^{-1}S^{-\top}]_{ii}}. Direct impacts average the \eqn{i = j} terms, total
#' impacts average row sums, and indirect = total - direct.
#'
#' The posterior draws are processed in parallel under the user's
#' [future::plan()] (see "Parallel post-processing" in [psp_control()]).
#'
#' @param fit A fitted `psprobit` object.
#' @param draws Number of posterior draws to use (default
#'   `control$impacts_draws` of the fit).
#' @param scale `"reduced_form"` (default) divides by \eqn{s_i} as above;
#'   `"lesage_pace"` uses \eqn{\phi(m_{it})} without it, as in the
#'   `spatialprobit` package, for comparison with published results.
#' @param type `"population"` (default, random effects integrated out) or
#'   `"unit"` (conditional on the random effects; needs `effects = "random"`).
#' @param by_period Logical; also report impacts averaged within each period
#'   (the default averages over all observations).
#' @param ... Unused.
#' @return An object of class `psp_impacts` with draws and a summary table.
#' @examples
#' W <- psp_weights(rook_weights(4, 4), style = "W")
#' d <- sim_psprobit(W, T = 3, beta = c(0, 1), rho = 0.4, seed = 2)
#' fit <- psprobit(y ~ x1, data = d, index = c("id", "time"), W = W, impacts = FALSE,
#'                 control = psp_control(iter = 200, warmup = 100, chains = 1, seed = 2))
#' psp_impacts(fit, draws = 50)
#' @export
psp_impacts <- function(fit, draws = NULL, scale = c("reduced_form", "lesage_pace"),
                        type = c("population", "unit"), by_period = FALSE) {
  if (!inherits(fit, "psprobit")) stop("'fit' must be a psprobit object.", call. = FALSE)
  scale <- match.arg(scale)
  type <- match.arg(type)
  vars <- names(fit$roles)[fit$roles == "x"]
  if (!length(vars)) stop("The model has no regressors for impacts.", call. = FALSE)
  if (type == "unit") {
    if (fit$effects == "pooled") {
      stop("type = \"unit\" needs a random-effects model.", call. = FALSE)
    }
    keep <- .psp_draw_subset_from(seq_along(fit$alpha_index), draws)
    idx <- fit$alpha_index[keep]
    alpha <- fit$alpha_draws[keep, , drop = FALSE]
  } else {
    idx <- .psp_draw_subset(fit, draws)
  }
  pm <- .psp_par_matrix(fit)[idx, , drop = FALSE]
  ordered <- inherits(fit, "psoprobit")
  ncat <- if (ordered) fit$J else 1L
  cats <- if (ordered) fit$levels else "1"
  periods <- c("all", if (by_period) as.character(fit$time_levels))
  theta_name <- if (length(fit$durbin_of)) stats::setNames(names(fit$durbin_of), fit$durbin_of) else character(0)
  xcols <- colnames(fit$X)
  vars_theta <- match(theta_name[vars], xcols) - 1L
  vars_theta[is.na(vars_theta)] <- -1L
  parts <- .psp_post_map(psp_rs_impacts, .psp_rs_post_model(fit),
                         .psp_rs_draws(fit, pm, if (type == "unit") alpha),
                         vars_beta = match(vars, xcols) - 1L, vars_theta = vars_theta,
                         by_period = by_period, lesage = scale == "lesage_pace")
  dims <- c(length(vars), ncat, 3, length(periods))
  res <- do.call(rbind, lapply(parts, matrix, ncol = prod(dims)))
  res <- array(res, c(length(idx), dims),
               dimnames = list(NULL, vars, cats, c("direct", "indirect", "total"), periods))
  smry <- .psp_impacts_summary(res)
  # draws: draws x variable x category x effect (x period when by_period)
  if (!by_period) res <- array(res, dim(res)[1:4], dimnames(res)[1:4])
  out <- list(draws = res, summary = smry, ordered = ordered,
              scale = scale, type = type, by_period = by_period, n_draws = length(idx))
  class(out) <- "psp_impacts"
  out
}

#' @rdname psp_impacts
#' @param obj A fitted `psprobit` object (method for `spatialreg::impacts`).
#' @exportS3Method spatialreg::impacts
impacts.psprobit <- function(obj, ...) psp_impacts(obj, ...)

.psp_impacts_summary <- function(res) {
  grid <- expand.grid(variable = dimnames(res)[[2]], category = dimnames(res)[[3]],
                      effect = dimnames(res)[[4]], period = dimnames(res)[[5]],
                      stringsAsFactors = FALSE)
  stats_ <- t(vapply(seq_len(nrow(grid)), function(i) {
    x <- res[, grid$variable[i], grid$category[i], grid$effect[i], grid$period[i]]
    c(mean = mean(x), sd = stats::sd(x), stats::quantile(x, c(0.025, 0.975), names = FALSE))
  }, numeric(4)))
  colnames(stats_) <- c("mean", "sd", "q2.5", "q97.5")
  cbind(grid, as.data.frame(stats_))
}

#' @export
print.psp_impacts <- function(x, digits = 4, ...) {
  cat(sprintf("Impacts on %s (posterior mean, sd, 95%% interval; %d draws)\n",
              if (x$ordered) "category probabilities" else "P(y = 1)", x$n_draws))
  if (x$scale == "lesage_pace") cat("Scale: LeSage-Pace (no reduced-form sd)\n")
  if (identical(x$type, "unit")) cat("Unit-specific (conditional on the random effects)\n")
  if (isTRUE(x$by_period)) {
    cat("Averaged over all periods; per-period impacts are in summary(<impacts>).\n")
  }
  s <- x$summary[x$summary$period == "all", ]
  for (cat_ in unique(s$category)) {
    if (x$ordered) cat(sprintf("\nCategory %s\n", cat_))
    for (eff in c("direct", "indirect", "total")) {
      sub <- s[s$category == cat_ & s$effect == eff, ]
      tab <- as.matrix(sub[, c("mean", "sd", "q2.5", "q97.5")])
      rownames(tab) <- sub$variable
      cat(sprintf("\n%s:\n", eff))
      print(round(tab, digits))
    }
  }
  invisible(x)
}

#' @export
summary.psp_impacts <- function(object, ...) {
  s <- object$summary
  if (!isTRUE(object$by_period)) s$period <- NULL
  s
}

# Posterior draws as a plain numeric matrix (draws x parameters)
.psp_par_matrix <- function(fit) {
  m <- posterior::as_draws_matrix(fit$draws)
  matrix(unclass(m), nrow(m), dimnames = list(NULL, colnames(m)))
}

# Evenly spaced subset of posterior draws
.psp_draw_subset <- function(fit, draws = NULL) {
  if (is.null(draws)) draws <- fit$control$impacts_draws
  .psp_draw_subset_from(seq_len(posterior::ndraws(fit$draws)), draws)
}

.psp_draw_subset_from <- function(pool, draws = NULL) {
  if (is.null(draws)) return(pool)
  pool[unique(round(seq(1, length(pool), length.out = min(draws, length(pool)))))]
}

# Rademacher probes, fixed per fit, for the stochastic diagonal estimator
.psp_probes <- function(fit) {
  if (fit$N <= fit$control$impacts_exact_max) return(NULL)
  P <- fit$control$impacts_probes
  if (.psp_latent_ar(fit)) P <- P * (fit$T + 1L)  # independent probes per period
  .psp_with_rng_restore({
    set.seed(fit$seed)
    matrix(sample(c(-1, 1), fit$N * P, replace = TRUE), fit$N, P)
  })
}

# Reduced-form quantities for one parameter draw, each of length n = N T
# (stacked time slow, unit fast), computed in Rust (pspcore::rf; the R
# original is .psp_ref_reduced_form in R/reference_postproc.R):
#   m  = S^{-1} X beta, s = reduced-form sd,
#   d1 = diag(S^{-1}), d2 = diag(S^{-1} W), t1 = S^{-1} 1, t2 = S^{-1} W 1.
# The latent covariance per period is sigma_alpha^2 S^{-1} V S^{-T} + K K' with
# K = S^{-1} R^{-1} (R = I - lambda M; K = S^{-1} without an error term) and
# V = I, or V = (B'B)^{-1} with B = I - phi W for spatial random effects.
# With 'alpha' (unit-specific), m = S^{-1}(X beta + alpha) and s drops sigma_alpha.
# For Student-t errors, 'mix' (N x Q mixing variances, .psp_t_mix) makes s an
# n x Q matrix: the sd given each simulated set of variances.
.psp_reduced_form <- function(fit, par, probes = NULL, need_impacts = FALSE, alpha = NULL,
                              mix = NULL) {
  pm <- matrix(par, 1, dimnames = list(NULL, names(par)))
  rf <- psp_rs_reduced_form(.psp_rs_post_model(fit, probes, mix),
                            .psp_rs_draws(fit, pm, if (!is.null(alpha)) matrix(alpha, 1)),
                            need_impacts)
  out <- list(m = rf$m, s = if (is.null(mix)) rf$s else matrix(rf$s, length(rf$m)))
  if (need_impacts) out[c("d1", "d2", "t1", "t2")] <- rf[c("d1", "d2", "t1", "t2")]
  out
}

# Full cut-point vector (-Inf, 0, free..., Inf) for one draw
.psp_cutpoints <- function(fit, par) {
  free <- par[grepl("^cut\\[", names(par))]
  c(-Inf, 0, unname(free), Inf)
}
