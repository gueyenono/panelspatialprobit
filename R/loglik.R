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
  out <- lapply(seq_along(idx), function(d) {
    .psp_ghk_loglik(fit, pm[d, ], reps, seed, alpha = if (conditional) alpha[d, ])
  })
  mcse <- vapply(out, `[[`, numeric(1), "mcse")
  .psp_ghk_check(mcse)
  structure(vapply(out, `[[`, numeric(1), "ll"), mcse = mcse, type = type, reps = reps)
}

#' @rdname psp_loglik
#' @export
logLik.psprobit <- function(object, reps = 500, seed = 1, ...) {
  .psp_ghk_normal_only(object)
  par <- object$coefficients
  g <- .psp_ghk_loglik(object, par, reps, seed)
  .psp_ghk_check(g$mcse)
  structure(g$ll, df = length(par), nobs = nobs.psprobit(object), mcse = g$mcse,
            class = "logLik")
}

#' @rdname psp_loglik
#' @export
psp_dic <- function(fit, draws = 100, reps = 200, seed = 1) {
  ll <- psp_loglik(fit, draws, reps, "marginal", seed)
  d_bar <- -2 * mean(ll)
  d_hat <- -2 * .psp_ghk_loglik(fit, fit$coefficients, reps, seed)$ll
  if (d_bar - d_hat < 0) {
    warning("Negative p_D: the GHK error dominates; increase 'reps'.", call. = FALSE)
  }
  p_d <- d_bar - d_hat
  list(dic = d_bar + p_d, p_d = p_d, d_bar = d_bar, d_hat = d_hat)
}

# GHK log-likelihood for one parameter vector. The latent vector is
# (z_1, ..., z_T) and, for the marginal likelihood of a random-effects model,
# alpha. With A_t = R_t S_t, the centred kernel is
#   sum_t |A_t z_t - R_t alpha|^2 + |alpha|^2 / sigma^2,
# whose precision is sparse (with |B alpha|^2 for spatial random effects);
# given alpha only the z blocks remain.
.psp_ghk_loglik <- function(fit, par, reps, seed, alpha = NULL) {
  N <- fit$N
  T <- fit$T
  I <- Matrix::Diagonal(N)
  beta <- par[colnames(fit$X)]
  rho <- if ("rho" %in% names(par)) par[["rho"]] else 0
  lambda <- if ("lambda" %in% names(par)) par[["lambda"]] else 0
  marginal <- fit$effects != "pooled" && is.null(alpha)
  xb <- matrix(fit$X %*% beta, N, T) + (if (is.null(alpha)) 0 else alpha)
  S <- lapply(seq_len(T), function(t) I - rho * .psp_Wt(fit$W, t))
  R <- lapply(seq_len(T), function(t) {
    if (lambda == 0) I else I - lambda * .psp_Wt(fit$M, t)
  })
  # full operator: blocks A_t = R_t S_t, and -gamma R_t below the diagonal
  # for latent AR(1); the latent mean solves A m = (R_t X_t beta)_t
  A <- Matrix::bdiag(Map(function(r, s) r %*% s, R, S))
  gamma <- if ("lag_z" %in% names(par)) par[["lag_z"]] else 0
  if (gamma != 0) {
    sub <- Matrix::bdiag(c(lapply(R[-1], function(r) -gamma * r), list(Matrix::Matrix(0, N, N))))
    A <- A + rbind(Matrix::Matrix(0, N, N * T), sub[seq_len(N * (T - 1)), , drop = FALSE])
  }
  Rall <- do.call(rbind, R)
  m <- as.vector(Matrix::solve(A, as.vector(Matrix::bdiag(R) %*% as.vector(xb))))
  Q <- Matrix::crossprod(A)
  if (marginal) {
    C <- -Matrix::crossprod(A, Rall)
    B <- if ("phi" %in% names(par)) I - par[["phi"]] * .psp_Wt(fit$W, 1) else I
    Qaa <- Reduce(`+`, lapply(R, Matrix::crossprod)) +
      Matrix::crossprod(B) / par[["sigma_alpha"]]^2
    Q <- rbind(cbind(Q, C), cbind(Matrix::t(C), Qaa))
    m <- c(m, numeric(N))
  }
  # bounds of the latent outcomes; missing outcomes and alpha are unrestricted
  kap <- if (inherits(fit, "psoprobit")) .psp_cutpoints(fit, par) else c(-Inf, 0, Inf)
  code <- if (inherits(fit, "psoprobit")) fit$y else fit$y + 1L
  lo <- ifelse(is.na(code), -Inf, kap[pmax(code, 1L, na.rm = TRUE)])
  hi <- ifelse(is.na(code), Inf, kap[pmin(code + 1L, length(kap), na.rm = TRUE)])
  if (marginal) {
    lo <- c(lo, rep(-Inf, N))
    hi <- c(hi, rep(Inf, N))
  }
  # Fill-reducing order, except that the unrestricted coordinates (alpha,
  # missing outcomes) go first: GHK draws from the last coordinate to the
  # first, so they are drawn last and integrated out exactly instead of adding
  # simulation noise.
  # The free block is ordered by a fill-reducing permutation of its own
  # precision, the rest by one of its pattern after the free block is
  # eliminated (approximated by Q_rr + C C').
  Q <- methods::as(Matrix::forceSymmetric(Q), "CsparseMatrix")
  free <- which(is.infinite(lo) & is.infinite(hi))
  rest <- setdiff(seq_along(lo), free)
  amd <- function(P) {  # fill-reducing order of the pattern of P
    P <- methods::as(methods::as(Matrix::forceSymmetric(P), "generalMatrix"), "CsparseMatrix")
    P@x[] <- 1
    P <- P + Matrix::Diagonal(nrow(P), max(Matrix::rowSums(P)) + 1)  # diagonally dominant
    Matrix::Cholesky(P, perm = TRUE, LDL = FALSE, super = FALSE)@perm + 1L
  }
  p <- if (length(free)) {
    C <- Q[rest, free, drop = FALSE]
    pr <- Matrix::tcrossprod(abs(C))
    c(free[amd(Q[free, free, drop = FALSE])], rest[amd(abs(Q[rest, rest]) + pr)])
  } else {
    amd(Q)
  }
  ch <- Matrix::Cholesky(Q[p, p], perm = FALSE, LDL = FALSE, super = FALSE)
  L <- methods::as(methods::as(ch, "CsparseMatrix"), "generalMatrix")
  # without random effects to integrate out (or latent dynamics linking
  # them), the periods are independent
  block <- if (marginal || gamma != 0) rep(0L, length(lo)) else rep(seq_len(T) - 1L, each = N)
  nb <- max(block) + 1L
  logw <- .psp_with_rng_restore({
    set.seed(seed)
    psp_ghk_cpp(L, (lo - m)[p], (hi - m)[p], block[p], nb, reps)
  })
  ll <- 0
  rel2 <- 0
  for (b in seq_len(nb)) {
    mx <- max(logw[, b])
    w <- exp(logw[, b] - mx)
    ll <- ll + mx + log(mean(w))
    rel2 <- rel2 + (stats::sd(w) / mean(w))^2 / reps
  }
  list(ll = ll, mcse = sqrt(rel2))
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
