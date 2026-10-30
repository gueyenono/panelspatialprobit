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
  res <- array(NA_real_, c(length(idx), length(vars), ncat, 3, length(periods)),
               dimnames = list(NULL, vars, cats, c("direct", "indirect", "total"), periods))
  probes <- .psp_probes(fit)
  theta_name <- if (length(fit$durbin_of)) stats::setNames(names(fit$durbin_of), fit$durbin_of) else character(0)
  period <- rep(seq_len(fit$T), each = fit$N)
  mix <- .psp_t_mix(fit)
  for (d in seq_along(idx)) {
    rf <- .psp_reduced_form(fit, pm[d, ], probes, need_impacts = TRUE,
                            alpha = if (type == "unit") alpha[d, ], mix = mix)
    s <- if (scale == "lesage_pace") rep(1, length(rf$m)) else rf$s
    G <- .psp_density_weights(fit, pm[d, ], rf$m, s)  # n x ncat
    for (v in vars) {
      b <- pm[d, v]
      th <- if (v %in% names(theta_name)) pm[d, theta_name[[v]]] else 0
      dir_i <- b * rf$d1 + th * rf$d2
      tot_i <- b * rf$t1 + th * rf$t2
      dir <- colMeans(G * dir_i)
      tot <- colMeans(G * tot_i)
      res[d, v, , , 1] <- cbind(dir, tot - dir, tot)
      if (by_period) {
        dir_t <- rowsum(G * dir_i, period) / fit$N  # T x ncat
        tot_t <- rowsum(G * tot_i, period) / fit$N
        for (t in seq_len(fit$T)) {
          res[d, v, , , t + 1] <- cbind(dir_t[t, ], tot_t[t, ] - dir_t[t, ], tot_t[t, ])
        }
      }
    }
  }
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
# (stacked time slow, unit fast):
#   m  = S^{-1} X beta, s = reduced-form sd,
#   d1 = diag(S^{-1}), d2 = diag(S^{-1} W), t1 = S^{-1} 1, t2 = S^{-1} W 1.
# The latent covariance per period is sigma_alpha^2 S^{-1} V S^{-T} + K K' with
# K = S^{-1} R^{-1} (R = I - lambda M; K = S^{-1} without an error term) and
# V = I, or V = (B'B)^{-1} with B = I - phi W for spatial random effects.
# With 'alpha' (unit-specific), m = S^{-1}(X beta + alpha) and s drops sigma_alpha.
# Time-varying weights are handled period by period. For Student-t errors,
# 'mix' (N x Q mixing variances, .psp_t_mix) makes s an n x Q matrix: the sd
# given each simulated set of variances, KK' becoming K diag(v) K'.
.psp_reduced_form <- function(fit, par, probes = NULL, need_impacts = FALSE, alpha = NULL,
                              mix = NULL) {
  N <- fit$N
  T <- fit$T
  beta <- par[colnames(fit$X)]
  xb <- matrix(fit$X %*% beta, N, T)
  sig2 <- if (fit$effects != "pooled") par[["sigma_alpha"]]^2 else 0
  if (!is.null(alpha)) {
    xb <- xb + alpha
    sig2 <- 0
  }
  rho <- if ("rho" %in% names(par)) par[["rho"]] else 0
  lambda <- if ("lambda" %in% names(par)) par[["lambda"]] else 0
  phi <- if ("phi" %in% names(par)) par[["phi"]] else 0
  gamma <- if ("lag_z" %in% names(par)) par[["lag_z"]] else 0
  if (gamma != 0) {
    return(.psp_rf_latent_ar(fit$W, fit$M, xb, sig2, rho, lambda, phi, gamma, probes,
                             need_impacts, mix))
  }
  varies <- (rho != 0 && is.list(fit$W)) || (lambda != 0 && is.list(fit$M))
  if (!varies) {
    out <- .psp_rf_block(.psp_Wt(fit$W, 1), .psp_Wt(fit$M, 1), xb, sig2, rho, lambda,
                         probes, need_impacts, phi, mix)
    for (nm in setdiff(names(out), "m")) {
      out[[nm]] <- if (is.matrix(out[[nm]])) out[[nm]][rep(seq_len(N), T), , drop = FALSE]
                   else rep(out[[nm]], T)
    }
    return(out)
  }
  blocks <- lapply(seq_len(T), function(t) {
    .psp_rf_block(.psp_Wt(fit$W, t), .psp_Wt(fit$M, t), xb[, t, drop = FALSE], sig2, rho,
                  lambda, probes, need_impacts, phi, mix)
  })
  stats::setNames(lapply(names(blocks[[1]]), function(nm) {
    parts <- lapply(blocks, `[[`, nm)
    if (is.matrix(parts[[1]])) do.call(rbind, parts) else unlist(parts, use.names = FALSE)
  }), names(blocks[[1]]))
}

# Reduced form for periods sharing the weights W and M ('xb' is N x periods);
# m has one entry per observation, the other quantities one per unit. 'phi'
# is the parameter of spatial random effects (alpha = phi W alpha + a).
.psp_rf_block <- function(W, M, xb, sig2, rho, lambda, probes, need_impacts, phi = 0,
                          mix = NULL) {
  N <- nrow(xb)
  T <- ncol(xb)
  W1 <- Matrix::rowSums(W)
  if (sig2 == 0) phi <- 0
  if (rho == 0 && lambda == 0 && phi == 0) {
    return(list(m = as.vector(xb), s = if (is.null(mix)) rep(sqrt(1 + sig2), N)
                                       else sqrt(sig2 + mix),
                d1 = rep(1, N), d2 = rep(0, N), t1 = rep(1, N), t2 = W1))
  }
  I <- Matrix::Diagonal(N)
  S <- I - rho * W
  R <- if (lambda != 0) I - lambda * M
  B <- if (phi != 0) I - phi * W
  if (is.null(probes)) {
    Sinv <- as.matrix(Matrix::solve(S))
    m <- Sinv %*% xb
    K <- if (is.null(R)) Sinv else Sinv %*% as.matrix(Matrix::solve(R))
    G <- if (is.null(B)) Sinv else Sinv %*% as.matrix(Matrix::solve(B))
    out <- list(m = as.vector(m), s = if (is.null(mix)) sqrt(sig2 * rowSums(G^2) + rowSums(K^2))
                                      else sqrt(sig2 * rowSums(G^2) + K^2 %*% mix))
    if (need_impacts) {
      out$d1 <- diag(Sinv)
      out$d2 <- rowSums(Sinv * t(as.matrix(W)))
      out$t1 <- rowSums(Sinv)
      out$t2 <- as.vector(Sinv %*% W1)
    }
    return(out)
  }
  # stochastic diagonals: diag(G) ~ mean over probes v of v * (G v)
  P <- ncol(probes)
  WV <- as.matrix(W %*% probes)
  U <- as.matrix(Matrix::solve(Matrix::t(S), probes))           # S^{-T} v
  UK <- if (is.null(R)) U else
    as.matrix(Matrix::solve(R, Matrix::solve(Matrix::t(R), U)))  # R^{-1} R^{-T} S^{-T} v
  UG <- if (is.null(B)) U else
    as.matrix(Matrix::solve(B, Matrix::solve(Matrix::t(B), U)))  # B^{-1} B^{-T} S^{-T} v
  RV <- if (is.null(mix)) NULL else if (is.null(R)) probes else
    as.matrix(Matrix::solve(R, probes))                          # R^{-1} v, for diag(K)
  sol <- as.matrix(Matrix::solve(S, cbind(xb, 1, W1, probes, WV, UG, UK, RV)))
  cols <- function(start, len) sol[, start + seq_len(len) - 1L, drop = FALSE]
  m <- cols(1, T)
  dSS <- rowMeans(probes * cols(T + 3 + 2 * P, P))
  dKK <- rowMeans(probes * cols(T + 3 + 3 * P, P))
  s2 <- if (is.null(mix)) sig2 * dSS + dKK else {
    # t errors: own error term exact, the others at their mean variance
    dK <- rowMeans(probes * cols(T + 3 + 4 * P, P))
    Ev <- mean(mix)
    sig2 * dSS + dK^2 * mix + pmax(dKK - dK^2, 0) * Ev
  }
  out <- list(m = as.vector(m), s = sqrt(pmax(s2, 1e-12)))
  if (need_impacts) {
    out$t1 <- sol[, T + 1]
    out$t2 <- sol[, T + 2]
    out$d1 <- rowMeans(probes * cols(T + 3, P))
    out$d2 <- rowMeans(probes * cols(T + 3 + P, P))
  }
  out
}

# Reduced form with latent AR(1) dynamics, z_t = S^{-1}(X_t beta + alpha +
# gamma z_{t-1} + R^{-1} eps_t) with no lag in the first period (time-invariant
# weights). The mean follows m_t = S^{-1}(xb_t + gamma m_{t-1}); the variance
# is the contemporaneous one (.psp_rf_block, which also handles t errors)
# plus the lagged errors, sum_{h>=1} gamma^{2h} diag(S^{-h} K K' S^{-hT}),
# plus the random effects, sigma^2 diag(G_t V G_t') with
# G_t = S^{-1}(I + gamma G_{t-1}). The lagged errors enter t errors at their
# mean variance. Exact for small N; with probes, both sums are estimated by
# simulating the recursions with independent Rademacher vectors (the probe
# matrix then has P (T + 1) columns). Impact derivatives are within-period.
.psp_rf_latent_ar <- function(W, M, xb, sig2, rho, lambda, phi, gamma, probes, need_impacts,
                              mix) {
  N <- nrow(xb)
  T <- ncol(xb)
  P <- if (is.null(probes)) 0L else ncol(probes) %/% (T + 1L)
  pr0 <- if (P) probes[, seq_len(P), drop = FALSE]
  base <- .psp_rf_block(W, M, xb[, 1, drop = FALSE], 0, rho, lambda, pr0, need_impacts,
                        0, mix)
  I <- Matrix::Diagonal(N)
  S <- I - rho * W
  R <- if (lambda != 0) I - lambda * M
  B <- if (phi != 0 && sig2 > 0) I - phi * W
  m <- matrix(0, N, T)
  prev <- 0
  for (t in seq_len(T)) {
    m[, t] <- as.vector(Matrix::solve(S, xb[, t] + gamma * prev))
    prev <- m[, t]
  }
  lag_var <- matrix(0, N, T)  # lagged-error variance
  re_var <- matrix(0, N, T)   # random-effect variance (sigma^2 excluded)
  if (!P) {
    Ph <- as.matrix(Matrix::solve(S, if (is.null(R)) diag(N) else as.matrix(Matrix::solve(R))))
    Vh <- if (is.null(B)) diag(N) else as.matrix(Matrix::solve(B))
    acc <- 0
    for (t in seq_len(T)) {
      if (t > 1) {
        Ph <- as.matrix(Matrix::solve(S, Ph))
        acc <- acc + gamma^(2 * (t - 1)) * rowSums(Ph^2)
      }
      lag_var[, t] <- acc
    }
    if (sig2 > 0) {
      H <- 0
      for (t in seq_len(T)) {
        H <- as.matrix(Matrix::solve(S, Vh + gamma * H))
        re_var[, t] <- rowSums(H^2)
      }
    }
  } else {
    x <- matrix(0, N, P)
    yv <- 0
    c0 <- probes[, P + seq_len(P), drop = FALSE]
    if (!is.null(B)) c0 <- as.matrix(Matrix::solve(B, c0))
    for (t in seq_len(T)) {
      if (t > 1) {
        w <- probes[, (t + 1) * P - P + seq_len(P), drop = FALSE]
        Kw <- if (is.null(R)) w else as.matrix(Matrix::solve(R, w))
        x <- gamma * as.matrix(Matrix::solve(S, x + as.matrix(Matrix::solve(S, Kw))))
      }
      lag_var[, t] <- rowMeans(x^2)
      if (sig2 > 0) {
        yv <- as.matrix(Matrix::solve(S, c0 + gamma * yv))
        re_var[, t] <- rowMeans(yv^2)
      }
    }
  }
  Ev <- if (is.null(mix)) 1 else mean(mix)
  extra <- as.vector(Ev * lag_var + sig2 * re_var)
  s0 <- base$s
  s <- if (is.matrix(s0)) sqrt(s0[rep(seq_len(N), T), , drop = FALSE]^2 + extra)
       else sqrt(rep(s0, T)^2 + extra)
  out <- list(m = as.vector(m), s = unname(s))
  if (need_impacts) for (nm in c("d1", "d2", "t1", "t2")) out[[nm]] <- unname(rep(base[[nm]], T))
  out
}

# n x ncat matrix of d P(category) / d(index) for each observation
# (s may be a matrix of mixture columns for t errors; columns are averaged)
.psp_density_weights <- function(fit, par, m, s) {
  s_full <- if (is.matrix(s)) s else rep_len(s, length(m))
  if (!inherits(fit, "psoprobit")) {
    return(matrix(.psp_mix_mean(function(s) stats::dnorm(m / s) / s, s_full), ncol = 1))
  }
  kap <- .psp_cutpoints(fit, par)
  sapply(seq_len(fit$J), function(j) {
    .psp_mix_mean(function(s) {
      (stats::dnorm((kap[j] - m) / s) - stats::dnorm((kap[j + 1] - m) / s)) / s
    }, s_full)
  })
}

# Full cut-point vector (-Inf, 0, free..., Inf) for one draw
.psp_cutpoints <- function(fit, par) {
  free <- par[grepl("^cut\\[", names(par))]
  c(-Inf, 0, unname(free), Inf)
}
