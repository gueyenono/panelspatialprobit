# Test oracle for the Rust post-processing (plan/06 §6.5, §6.10): the R
# log-determinant, reduced form, impacts and predict code as it was before the
# port, kept verbatim (renamed .psp_ref_*). Used only by
# test-postproc-oracle.R; never called by the package itself.

# Log-determinant log|I - rho W| on a fine grid strictly inside (rmin, rmax).
.psp_ref_lndet <- function(W, rmin = -1, rmax = 1, method = c("auto", "eigen", "sparse"),
                       step = 5e-4) {
  method <- match.arg(method)
  N <- nrow(W)
  if (method == "auto") method <- if (N <= 3000) "eigen" else "sparse"
  grid <- seq(rmin + step / 2, rmax - step / 2, by = step)
  if (method == "eigen") {
    ev <- eigen(as.matrix(W), only.values = TRUE)$values
    re <- Re(ev)
    im <- Im(ev)
    lndet <- vapply(grid, function(r) {
      0.5 * sum(log((1 - r * re)^2 + (r * im)^2))
    }, numeric(1))
  } else {
    coarse <- unique(c(seq(rmin + step / 2, rmax - step / 2, length.out = 120),
                       rmax - c(0.02, 0.01, 0.005, 0.002) - step / 2))
    coarse <- sort(coarse[coarse > rmin & coarse < rmax])
    I <- Matrix::Diagonal(N)
    vals <- vapply(coarse, function(r) {
      as.numeric(Matrix::determinant(I - r * W, logarithm = TRUE)$modulus)
    }, numeric(1))
    lndet <- stats::splinefun(coarse, vals, method = "natural")(grid)
  }
  list(grid = grid, lndet = lndet, method = method)
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
.psp_ref_reduced_form <- function(fit, par, probes = NULL, need_impacts = FALSE, alpha = NULL,
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
    return(.psp_ref_rf_latent_ar(fit$W, fit$M, xb, sig2, rho, lambda, phi, gamma, probes,
                             need_impacts, mix))
  }
  varies <- (rho != 0 && is.list(fit$W)) || (lambda != 0 && is.list(fit$M))
  if (!varies) {
    out <- .psp_ref_rf_block(.psp_Wt(fit$W, 1), .psp_Wt(fit$M, 1), xb, sig2, rho, lambda,
                         probes, need_impacts, phi, mix)
    for (nm in setdiff(names(out), "m")) {
      out[[nm]] <- if (is.matrix(out[[nm]])) out[[nm]][rep(seq_len(N), T), , drop = FALSE]
                   else rep(out[[nm]], T)
    }
    return(out)
  }
  blocks <- lapply(seq_len(T), function(t) {
    .psp_ref_rf_block(.psp_Wt(fit$W, t), .psp_Wt(fit$M, t), xb[, t, drop = FALSE], sig2, rho,
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
.psp_ref_rf_block <- function(W, M, xb, sig2, rho, lambda, probes, need_impacts, phi = 0,
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
# is the contemporaneous one (.psp_ref_rf_block, which also handles t errors)
# plus the lagged errors, sum_{h>=1} gamma^{2h} diag(S^{-h} K K' S^{-hT}),
# plus the random effects, sigma^2 diag(G_t V G_t') with
# G_t = S^{-1}(I + gamma G_{t-1}). The lagged errors enter t errors at their
# mean variance. Exact for small N; with probes, both sums are estimated by
# simulating the recursions with independent Rademacher vectors (the probe
# matrix then has P (T + 1) columns). Impact derivatives are within-period.
.psp_ref_rf_latent_ar <- function(W, M, xb, sig2, rho, lambda, phi, gamma, probes, need_impacts,
                              mix) {
  N <- nrow(xb)
  T <- ncol(xb)
  P <- if (is.null(probes)) 0L else ncol(probes) %/% (T + 1L)
  pr0 <- if (P) probes[, seq_len(P), drop = FALSE]
  base <- .psp_ref_rf_block(W, M, xb[, 1, drop = FALSE], 0, rho, lambda, pr0, need_impacts,
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
.psp_ref_density_weights <- function(fit, par, m, s) {
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

# The impacts draws array of psp_impacts() (draws x variable x category x
# effect x period), without the summary.
.psp_ref_impacts <- function(fit, draws = NULL, scale = c("reduced_form", "lesage_pace"),
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
    rf <- .psp_ref_reduced_form(fit, pm[d, ], probes, need_impacts = TRUE,
                            alpha = if (type == "unit") alpha[d, ], mix = mix)
    s <- if (scale == "lesage_pace") rep(1, length(rf$m)) else rf$s
    G <- .psp_ref_density_weights(fit, pm[d, ], rf$m, s)  # n x ncat
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
  res
}

# The posterior means of the link and of the category probabilities of
# predict(), in the sorted (time-slow, unit-fast) order.
.psp_ref_predict <- function(object, draws = NULL) {
  idx <- .psp_draw_subset(object, draws)
  pm <- .psp_par_matrix(object)[idx, , drop = FALSE]
  probes <- .psp_probes(object)
  mix <- .psp_t_mix(object)
  ordered <- inherits(object, "psoprobit")
  n <- object$N * object$T
  acc_m <- numeric(n)
  acc_p <- matrix(0, n, if (ordered) object$J else 1)
  for (d in seq_along(idx)) {
    rf <- .psp_ref_reduced_form(object, pm[d, ], probes, mix = mix)
    s <- rf$s
    acc_m <- acc_m + rf$m
    if (ordered) {
      kap <- .psp_cutpoints(object, pm[d, ])
      for (j in seq_len(object$J)) {
        acc_p[, j] <- acc_p[, j] + .psp_mix_mean(function(s) {
          stats::pnorm((kap[j + 1] - rf$m) / s) - stats::pnorm((kap[j] - rf$m) / s)
        }, s)
      }
    } else {
      acc_p[, 1] <- acc_p[, 1] + .psp_mix_mean(function(s) stats::pnorm(rf$m / s), s)
    }
  }
  acc_m <- acc_m / length(idx)
  acc_p <- acc_p / length(idx)
  list(link = acc_m, prob = acc_p)
}
