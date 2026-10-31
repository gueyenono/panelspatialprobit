#' Posterior predictive checks
#'
#' Simulates replicated outcomes from the posterior and compares summary
#' statistics of the replicates with those of the observed outcomes
#' (Gelman, Meng and Stern 1996). Replicates are drawn for the same units,
#' periods and regressors, conditional on the posterior draws of the random
#' effects kept by the sampler (pooled models have none). Missing outcomes are
#' left out of every statistic.
#'
#' The statistics are:
#' * `share`: the share of each outcome category (of \eqn{y = 1} for binary
#'   models) in each period;
#' * `moran`: Moran's I of the outcome (category codes for ordered models) in
#'   each period, which checks the spatial dependence;
#' * `persistence`: the share of units whose outcome equals that of the
#'   previous period, which checks the dependence over time (\eqn{T \ge 2}).
#'
#' The posterior draws are processed in parallel under the user's
#' [future::plan()] (see "Parallel post-processing" in [psp_control()]).
#'
#' @param fit A fitted `psprobit` object.
#' @param draws Number of posterior draws to replicate (default
#'   `control$impacts_draws`).
#' @param seed Optional seed for the replicates.
#' @return An object of class `psp_ppc`: a list with, for each statistic, the
#'   observed values (`observed`), the replicated values (`replicated`, draws x
#'   components) and posterior predictive p-values
#'   \eqn{P(T(y^{rep}) \ge T(y))} (`ppp`). Its `print()` method shows a summary.
#' @references Gelman, A., Meng, X.-L. and Stern, H. (1996). Posterior
#'   predictive assessment of model fitness via realized discrepancies.
#'   *Statistica Sinica* 6, 733-807.
#' @examples
#' W <- psp_weights(rook_weights(5, 5), style = "W")
#' d <- sim_psprobit(W, T = 4, beta = c(0, 1, -1), rho = 0.5, seed = 1)
#' fit <- psprobit(y ~ x1 + x2, d, c("id", "time"), W, impacts = FALSE,
#'                 control = psp_control(iter = 300, warmup = 200, chains = 1, seed = 1))
#' psp_ppc(fit, draws = 50, seed = 1)
#' @export
psp_ppc <- function(fit, draws = NULL, seed = NULL) {
  if (!inherits(fit, "psprobit")) stop("'fit' must be a psprobit object.", call. = FALSE)
  if (fit$effects != "pooled") {
    keep <- .psp_draw_subset_from(seq_along(fit$alpha_index),
                                  if (is.null(draws)) fit$control$impacts_draws else draws)
    idx <- fit$alpha_index[keep]
    alpha <- fit$alpha_draws[keep, , drop = FALSE]
  } else {
    idx <- .psp_draw_subset(fit, draws)
    alpha <- NULL
  }
  pm <- .psp_par_matrix(fit)[idx, , drop = FALSE]
  T <- fit$T
  obs <- fit$y  # sorted: time slow, unit fast; NA = missing
  # one random number stream per draw (Rust); seed = NULL draws the base seed
  # from the current RNG state, which is left unchanged
  if (is.null(seed)) seed <- .psp_with_rng_restore(sample.int(.Machine$integer.max, 1))
  dr <- .psp_rs_draws(fit, pm, alpha)
  dr$seeds <- .psp_rs_seeds(seed, length(idx))
  parts <- .psp_post_map(.psp_rs_ppc_call, .psp_rs_post_model(fit, NULL, NULL), dr)
  yrep <- do.call(cbind, lapply(parts, matrix, nrow = fit$N * T))
  yrep[is.na(obs), ] <- NA
  stats <- list(share = .psp_stat_share, moran = .psp_stat_moran)
  if (T >= 2) stats$persistence <- .psp_stat_persistence
  out <- lapply(stats, function(f) {
    o <- f(obs, fit)
    r <- t(apply(yrep, 2, f, fit = fit))
    if (length(o) == 1) r <- t(r)
    colnames(r) <- names(o)
    list(observed = o, replicated = r,
         ppp = colMeans(sweep(r, 2, o, ">=") * 1, na.rm = TRUE))
  })
  structure(c(out, list(n_draws = length(idx))), class = "psp_ppc")
}

.psp_stat_share <- function(y, fit) {
  per <- rep(fit$time_levels, each = fit$N)
  cats <- if (inherits(fit, "psoprobit")) seq_len(fit$J) else 1L
  labs <- if (inherits(fit, "psoprobit")) fit$levels else "1"
  out <- unlist(lapply(seq_along(cats), function(j) {
    tapply(y == cats[j], per, mean, na.rm = TRUE)
  }))
  names(out) <- paste0("P(y=", rep(labs, each = fit$T), ") t=", rep(fit$time_levels, length(cats)))
  out
}

.psp_stat_moran <- function(y, fit) {
  Y <- matrix(y, fit$N, fit$T)
  out <- vapply(seq_len(fit$T), function(t) {
    W <- .psp_Wt(fit$W, t)
    v <- Y[, t]
    ok <- !is.na(v)
    if (sum(ok) < 3 || stats::var(v[ok]) == 0) return(NA_real_)
    Wt <- W[ok, ok, drop = FALSE]
    s0 <- sum(Wt)
    if (s0 == 0) return(NA_real_)
    e <- v[ok] - mean(v[ok])
    sum(ok) / s0 * sum(e * as.vector(Wt %*% e)) / sum(e^2)
  }, numeric(1))
  names(out) <- paste0("Moran t=", fit$time_levels)
  out
}

.psp_stat_persistence <- function(y, fit) {
  Y <- matrix(y, fit$N, fit$T)
  out <- vapply(2:fit$T, function(t) mean(Y[, t] == Y[, t - 1], na.rm = TRUE), numeric(1))
  names(out) <- paste0("persist t=", fit$time_levels[-1])
  out
}

#' @export
print.psp_ppc <- function(x, digits = 3, ...) {
  cat(sprintf("Posterior predictive checks (%d replicates)\n", x$n_draws))
  cat("ppp = P(T(y_rep) >= T(y)); values near 0 or 1 signal misfit.\n")
  for (nm in setdiff(names(x), "n_draws")) {
    s <- x[[nm]]
    tab <- cbind(observed = s$observed, rep_mean = colMeans(s$replicated, na.rm = TRUE),
                 rep_q2.5 = apply(s$replicated, 2, stats::quantile, 0.025, na.rm = TRUE),
                 rep_q97.5 = apply(s$replicated, 2, stats::quantile, 0.975, na.rm = TRUE),
                 ppp = s$ppp)
    cat(sprintf("\n%s:\n", nm))
    print(round(tab, digits))
  }
  invisible(x)
}
