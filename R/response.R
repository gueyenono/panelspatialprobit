#' Dynamic response paths
#'
#' Traces how a change in a regressor propagates over time in a dynamic model
#' ([psp_dynamic()]): the shock changes the outcome probabilities in the period
#' it occurs (the short-run impact), which changes the outcomes, which enter
#' the next periods through the lagged outcomes (`lag_y`, `W_lag_y`), and so on.
#'
#' For each posterior draw, outcome histories are simulated forward from the
#' start period under the baseline (the observed regressors) and under the
#' shock, with the same random numbers in both scenarios. Outcomes before the
#' start period are the observed ones. In every period the outcome
#' probabilities given the simulated history are computed exactly, and the
#' response is their average difference between the two scenarios. The
#' simulation conditions on the posterior draws of the random effects kept by
#' the sampler (as [psp_ppc()] does), so it describes the sample units.
#' Correlated-random-effects columns (`mean_x`, `cre_x_*`) and the initial
#' conditions are held fixed.
#'
#' With latent dynamics (`psp_dynamic(latent = TRUE)`) the shock propagates
#' through the latent mean, \eqn{\Delta m_{t+h} = (\gamma S^{-1})^h
#' \Delta m_t}, and the responses are computed exactly from the reduced form
#' (no simulation; `nsim` is unused). Static models are accepted too; their
#' responses vanish after the shock ends.
#'
#' The posterior draws are processed in parallel under the user's
#' [future::plan()] (see "Parallel post-processing" in [psp_control()]).
#'
#' @param fit A fitted `psprobit` object.
#' @param variable Name of the regressor to shock.
#' @param shock Size of the change in `variable`.
#' @param start Period (a value of the time index) in which the shock starts;
#'   default the first modelled period.
#' @param units Identifiers of the shocked units (default: all units).
#' @param type `"temporary"` (the shock lasts one period) or `"permanent"` (it
#'   lasts from `start` to the last period).
#' @param draws Number of posterior draws.
#' @param nsim Number of simulated outcome histories per draw.
#' @param seed Optional seed for the simulation.
#' @return An object of class `psp_response`: a list with `summary`, a data
#'   frame with the posterior mean, sd and 95% interval of the response by
#'   `horizon` (0 = start period), `period`, `category` (ordered models) and
#'   `units` (`"all"` averages the change over every unit; when only some
#'   units are shocked, `"shocked"` averages over them and `"others"` over the
#'   remaining ones), and `draws`, the
#'   responses for each posterior draw. `print()` and `plot()` methods are
#'   available.
#' @examples
#' W <- psp_weights(rook_weights(5, 5), style = "W")
#' d <- sim_psprobit(W, T = 6, beta = c(-0.3, 1, -1), rho = 0.3, gamma = 0.8, seed = 1)
#' fit <- psprobit(y ~ x1 + x2, d, c("id", "time"), W, dynamic = psp_dynamic(),
#'                 impacts = FALSE,
#'                 control = psp_control(iter = 300, warmup = 200, chains = 1, seed = 1))
#' psp_response_path(fit, "x1", draws = 50, nsim = 5, seed = 1)
#' @export
psp_response_path <- function(fit, variable, shock = 1, start = NULL, units = NULL,
                              type = c("temporary", "permanent"), draws = 200, nsim = 20,
                              seed = NULL) {
  if (!inherits(fit, "psprobit")) stop("'fit' must be a psprobit object.", call. = FALSE)
  type <- match.arg(type)
  vars <- names(fit$roles)[fit$roles == "x"]
  if (length(variable) != 1 || !variable %in% vars) {
    stop("'variable' must be one of: ", paste(vars, collapse = ", "), call. = FALSE)
  }
  N <- fit$N
  T <- fit$T
  p0 <- if (is.null(start)) 1L else match(as.character(start), as.character(fit$time_levels))
  if (length(p0) != 1 || is.na(p0)) {
    stop("'start' must be one of the modelled periods: ",
         paste(fit$time_levels, collapse = ", "), call. = FALSE)
  }
  shocked <- if (is.null(units)) seq_len(N) else match(as.character(units), fit$unit_ids)
  if (anyNA(shocked)) stop("Unknown 'units'.", call. = FALSE)
  e <- numeric(N)
  e[shocked] <- shock
  periods <- seq.int(p0, T)

  if (fit$effects != "pooled") {
    keep <- .psp_draw_subset_from(seq_along(fit$alpha_index), draws)
    idx <- fit$alpha_index[keep]
    alpha <- fit$alpha_draws[keep, , drop = FALSE]
  } else {
    idx <- .psp_draw_subset(fit, draws)
    alpha <- NULL
  }
  pm <- .psp_par_matrix(fit)[idx, , drop = FALSE]
  ordered <- inherits(fit, "psoprobit")
  ncat <- if (ordered) fit$J else 1L
  cats <- if (ordered) fit$levels else "1"
  cn <- colnames(fit$X)
  ycols <- if (ordered) paste0("y", fit$levels[-1]) else "y"
  lag_cols <- paste0("lag_", ycols)
  wlag_cols <- paste0("W_lag_", ycols)
  dyn <- cn %in% c(lag_cols, wlag_cols)
  wvar <- paste0("W_", variable)
  probes <- .psp_probes(fit)
  mix <- .psp_t_mix(fit)
  groups <- list(all = seq_len(N))
  if (length(shocked) < N) {
    groups$shocked <- shocked
    groups$others <- setdiff(seq_len(N), shocked)
  }

  dims <- c(length(periods), ncat, length(groups))
  dn <- list(NULL, as.character(fit$time_levels[periods]), cats, names(groups))
  if (is.null(seed)) seed <- .psp_with_rng_restore(sample.int(.Machine$integer.max, 1))
  dr <- .psp_rs_draws(fit, pm, alpha)
  if (.psp_latent_ar(fit)) {
    # latent dynamics: exact reduced-form probabilities with and without the shock
    f1 <- fit
    for (p in if (type == "permanent") periods else p0) {
      rows <- (p - 1) * N + seq_len(N)
      f1$X[rows, variable] <- f1$X[rows, variable] + e
      if (wvar %in% cn) f1$X[rows, wvar] <- f1$X[rows, wvar] + as.vector(fit$W %*% e)
    }
    pr <- lapply(list(fit, f1), function(f) {
      parts <- .psp_post_map(psp_rs_probs, .psp_rs_post_model(f, probes, mix), dr)
      array(unlist(parts), c(N * T, ncat, length(idx)))
    })
    diff <- pr[[2]] - pr[[1]]
    res <- array(NA_real_, c(length(idx), dims), dimnames = dn)
    for (h in seq_along(periods)) {
      rows <- (periods[h] - 1) * N + seq_len(N)
      for (g in names(groups)) {
        res[, h, , g] <- t(colMeans(diff[rows, , , drop = FALSE][groups[[g]], , , drop = FALSE]))
      }
    }
  } else {
    col0 <- function(nm) {
      j <- match(nm, cn) - 1L
      ifelse(is.na(j), -1L, j)
    }
    spec <- list(dynamic = which(dyn) - 1L, lag = col0(lag_cols), wlag = col0(wlag_cols),
                 var = col0(variable), wvar = col0(wvar), shock = as.double(e),
                 start = as.integer(p0 - 1L), permanent = type == "permanent",
                 groups = lapply(groups, function(g) as.integer(g - 1L)),
                 nsim = as.integer(nsim))
    dr$seeds <- .psp_rs_seeds(seed, length(idx))
    parts <- .psp_post_map(.psp_rs_response_call, .psp_rs_post_model(fit, probes, mix), dr,
                           spec = spec)
    res <- array(do.call(rbind, lapply(parts, matrix, ncol = prod(dims))),
                 c(length(idx), dims), dimnames = dn)
  }
  grid <- expand.grid(period = dimnames(res)[[2]], category = cats, units = names(groups),
                      stringsAsFactors = FALSE)
  grid$horizon <- match(grid$period, dimnames(res)[[2]]) - 1L
  st <- t(vapply(seq_len(nrow(grid)), function(i) {
    x <- res[, grid$period[i], grid$category[i], grid$units[i]]
    c(mean = mean(x), sd = stats::sd(x), stats::quantile(x, c(0.025, 0.975), names = FALSE))
  }, numeric(4)))
  colnames(st) <- c("mean", "sd", "q2.5", "q97.5")
  smry <- cbind(grid[, c("horizon", "period", "category", "units")], as.data.frame(st))
  if (!ordered) smry$category <- NULL
  structure(list(summary = smry, draws = res, variable = variable, shock = shock, type = type,
                 ordered = ordered, dynamic = !is.null(fit$dynamic), n_draws = length(idx),
                 nsim = nsim),
            class = "psp_response")
}

#' @export
print.psp_response <- function(x, digits = 4, ...) {
  cat(sprintf("Response of %s to a %s shock of %s in %s (%d draws x %d histories)\n",
              if (x$ordered) "category probabilities" else "P(y = 1)", x$type,
              format(x$shock), x$variable, x$n_draws, x$nsim))
  if (!x$dynamic) cat("The model is static: responses end with the shock.\n")
  s <- x$summary
  for (u in unique(s$units)) {
    sub <- s[s$units == u, ]
    cat(sprintf("\nAverage change over %s units:\n",
                c(all = "all", shocked = "the shocked", others = "the other")[[u]]))
    tab <- sub[, intersect(c("horizon", "period", "category", "mean", "sd", "q2.5", "q97.5"),
                           names(sub))]
    num <- vapply(tab, is.numeric, logical(1)) & !names(tab) %in% "horizon"
    tab[num] <- lapply(tab[num], round, digits)
    print(tab, row.names = FALSE)
  }
  invisible(x)
}

#' @export
plot.psp_response <- function(x, units = "all", ...) {
  s <- x$summary[x$summary$units == units, ]
  cats <- if (x$ordered) unique(s$category) else "1"
  old <- graphics::par(mfrow = c(1, length(cats)), mar = c(4, 4, 2, 1))
  on.exit(graphics::par(old))
  for (k in cats) {
    sk <- if (x$ordered) s[s$category == k, ] else s
    graphics::plot(sk$horizon, sk$mean, type = "n", ylim = range(0, sk$q2.5, sk$q97.5),
                   xlab = "horizon", ylab = "change in probability",
                   main = if (x$ordered) paste("Category", k) else x$variable)
    graphics::polygon(c(sk$horizon, rev(sk$horizon)), c(sk$q2.5, rev(sk$q97.5)),
                      col = "grey85", border = NA)
    graphics::lines(sk$horizon, sk$mean, type = "b", pch = 19)
    graphics::abline(h = 0, lty = 2)
  }
  invisible(x)
}
