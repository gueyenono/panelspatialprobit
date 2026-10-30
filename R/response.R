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

  res <- array(NA_real_, c(length(idx), length(periods), ncat, length(groups)),
               dimnames = list(NULL, as.character(fit$time_levels[periods]), cats,
                               names(groups)))
  .psp_with_rng_restore({
    if (!is.null(seed)) set.seed(seed)
    for (d in seq_along(idx)) {
      par <- pm[d, ]
      beta <- par[cn]
      rho <- if ("rho" %in% names(par)) par[["rho"]] else 0
      lambda <- if ("lambda" %in% names(par)) par[["lambda"]] else 0
      a <- if (is.null(alpha)) 0 else alpha[d, ]
      kap <- if (ordered) .psp_cutpoints(fit, par) else c(-Inf, 0, Inf)
      if (.psp_latent_ar(fit)) {
        # latent dynamics: exact reduced-form probabilities with and without the shock
        f1 <- fit
        for (p in if (type == "permanent") periods else p0) {
          rows <- (p - 1) * N + seq_len(N)
          f1$X[rows, variable] <- f1$X[rows, variable] + e
          if (wvar %in% cn) f1$X[rows, wvar] <- f1$X[rows, wvar] + as.vector(fit$W %*% e)
        }
        a_d <- if (is.null(alpha)) NULL else alpha[d, ]
        pr <- lapply(list(fit, f1), function(f) {
          rf <- .psp_reduced_form(f, par, probes, alpha = a_d, mix = mix)
          .psp_category_probs(fit, par, rf$m, rf$s)
        })
        diff <- pr[[2]] - pr[[1]]
        for (h in seq_along(periods)) {
          rows <- (periods[h] - 1) * N + seq_len(N)
          for (g in names(groups)) {
            res[d, h, , g] <- colMeans(diff[rows, , drop = FALSE][groups[[g]], , drop = FALSE])
          }
        }
        next
      }
      Y <- list(base = NULL, shock = NULL)  # simulated outcomes of the previous period
      for (h in seq_along(periods)) {
        p <- periods[h]
        Wp <- .psp_Wt(fit$W, p)
        Mp <- .psp_Wt(fit$M, p)
        rows <- (p - 1) * N + seq_len(N)
        Xp <- fit$X[rows, , drop = FALSE]
        fixed <- as.vector(Xp[, !dyn, drop = FALSE] %*% beta[!dyn]) + a
        shift <- 0
        if (h == 1 || type == "permanent") {
          shift <- e * beta[[variable]]
          if (wvar %in% cn) shift <- shift + as.vector(Wp %*% e) * beta[[wvar]]
        }
        eta <- lapply(c(base = 0, shock = 1), function(sc) {
          lagged <- if (h == 1) {
            matrix(as.vector(Xp[, dyn, drop = FALSE] %*% beta[dyn]), N, nsim)
          } else {
            .psp_lag_effect(Y[[if (sc == 1) "shock" else "base"]], ordered, fit$J, Wp,
                            beta, lag_cols, wlag_cols, cn)
          }
          fixed + sc * shift + lagged
        })
        rf <- lapply(eta, function(x) {
          .psp_rf_block(Wp, Mp, x, 0, rho, lambda, probes, FALSE, mix = mix)
        })
        prob <- lapply(rf, function(r) {
          m <- matrix(r$m, N)
          # columns of r$s: mixture draws for t errors (one column otherwise)
          sm <- if (is.matrix(r$s)) r$s else matrix(r$s, N)
          vapply(seq_len(ncat), function(j) {
            rowMeans(vapply(seq_len(ncol(sm)), function(q) {
              s <- sm[, q]
              if (!ordered) return(rowMeans(stats::pnorm(m / s)))
              rowMeans(stats::pnorm((kap[j + 1] - m) / s) - stats::pnorm((kap[j] - m) / s))
            }, numeric(N)))
          }, numeric(N))
        })
        diff <- matrix(prob$shock - prob$base, N, ncat)
        for (g in names(groups)) {
          res[d, h, , g] <- colMeans(diff[groups[[g]], , drop = FALSE])
        }
        if (h < length(periods)) {
          # common latent errors: z = S^{-1}(eta + R^{-1} eps)
          eps <- matrix(.psp_rerr(fit, N * nsim), N, nsim)
          I <- Matrix::Diagonal(N)
          if (lambda != 0) eps <- as.matrix(Matrix::solve(I - lambda * Mp, eps))
          for (sc in names(eta)) {
            z <- eta[[sc]] + eps
            if (rho != 0) z <- as.matrix(Matrix::solve(I - rho * Wp, z))
            Y[[sc]] <- matrix(findInterval(z, kap[2:(length(kap) - 1)],
                                           left.open = TRUE) + if (ordered) 1L else 0L, N)
          }
        }
      }
    }
  })
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

# n x ncat matrix of outcome probabilities given the reduced form (m, s);
# s may hold mixture columns (t errors)
.psp_category_probs <- function(fit, par, m, s) {
  if (!inherits(fit, "psoprobit")) {
    return(matrix(.psp_mix_mean(function(s) stats::pnorm(m / s), s), ncol = 1))
  }
  kap <- .psp_cutpoints(fit, par)
  vapply(seq_len(fit$J), function(j) {
    .psp_mix_mean(function(s) stats::pnorm((kap[j + 1] - m) / s) - stats::pnorm((kap[j] - m) / s),
                  s)
  }, numeric(length(m)))
}

# Contribution of the lagged outcomes Y (N x nsim outcome codes) to the
# linear index: lag_* and W_lag_* columns rebuilt from Y.
.psp_lag_effect <- function(Y, ordered, J, Wp, beta, lag_cols, wlag_cols, cn) {
  out <- matrix(0, nrow(Y), ncol(Y))
  levels_ <- if (ordered) 2:J else 1L
  for (k in seq_along(levels_)) {
    ind <- (Y == levels_[k]) * 1
    if (lag_cols[k] %in% cn) out <- out + beta[[lag_cols[k]]] * ind
    if (wlag_cols[k] %in% cn) out <- out + beta[[wlag_cols[k]]] * as.matrix(Wp %*% ind)
  }
  out
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
