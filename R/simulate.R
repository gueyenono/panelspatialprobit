#' Simulate panel spatial probit and ordered probit data
#'
#' Draws data from the latent model
#' \deqn{z_t = \rho W z_t + X_t\beta + W X_t\theta + \alpha + u_t,
#'   \quad u_t = \lambda M u_t + \varepsilon_t}
#' with standard normal regressors. The outcome is binary (`cutpoints = NULL`)
#' or ordered with cut-points \eqn{(0, \kappa_2, \dots, \kappa_{J-1})}.
#'
#' @param W Spatial weights (`N x N`), used as given, or a list of `T` such
#'   matrices (time-varying weights, period by period).
#' @param T Number of periods.
#' @param beta Coefficients; the first is the intercept, the others belong to
#'   regressors `x1`, `x2`, ...
#' @param rho Spatial autoregressive parameter (0 for `"slx"`/`"none"`).
#' @param sigma_alpha Standard deviation of the random effects (0 for pooled).
#' @param theta Optional Durbin coefficients (one per regressor), giving the
#'   effects of \eqn{WX}.
#' @param gamma Optional coefficient(s) of the lagged outcome
#'   \eqn{y_{i,t-1}} (one per category above the first for ordered outcomes),
#'   as in [psp_dynamic()]. The first period is generated without lags.
#' @param eta Optional coefficient(s) of the spatially lagged outcome
#'   \eqn{(W y_{t-1})_i}.
#' @param lambda Spatial error parameter (0 for no spatial error term).
#' @param gamma_z Coefficient of latent AR(1) dynamics, \eqn{\gamma z_{t-1}}
#'   (as for `psp_dynamic(latent = TRUE)`; no lag in the first period).
#' @param df Optional degrees of freedom of Student-t latent errors (as for
#'   `errors = psp_t(df)`); normal errors by default.
#' @param phi Spatial parameter of the random effects, \eqn{\alpha = \phi W
#'   \alpha + a} with \eqn{a_i \sim N(0, \sigma_\alpha^2)} (as for
#'   `effects = "spatial_random"`; needs a single `W`).
#' @param M Error weights (default `W`); a matrix or a list of `T` matrices.
#' @param cutpoints Optional free cut-points \eqn{\kappa_2 < \dots <
#'   \kappa_{J-1}}, all positive (the first cut-point is fixed at 0). For
#'   example `cutpoints = 1` gives \eqn{J = 3} categories.
#' @param seed Optional seed.
#' @return A data frame with columns `id`, `time`, `y` and the regressors, with
#'   attributes `truth` (parameters, random effects and latent values) and `W`
#'   (the weights as passed, after [psp_weights()]).
#' @examples
#' W <- psp_weights(rook_weights(3, 3), style = "W")
#' d <- sim_psprobit(W, T = 2, beta = c(0, 1), rho = 0.3, seed = 1)
#' head(d)
#' @export
sim_psprobit <- function(W, T, beta = c(0, 1, -1), rho = 0.5, sigma_alpha = 1,
                         theta = NULL, cutpoints = NULL, seed = NULL, lambda = 0,
                         M = NULL, gamma = NULL, eta = NULL, phi = 0, df = NULL,
                         gamma_z = 0) {
  W <- if (.psp_is_weights_list(W)) lapply(W, psp_weights) else psp_weights(W)
  if (is.list(W) && length(W) != T) stop("A list 'W' needs one matrix per period.", call. = FALSE)
  if (is.null(M)) M <- W
  M <- if (.psp_is_weights_list(M)) lapply(M, psp_weights) else psp_weights(M)
  if (is.list(M) && length(M) != T) stop("A list 'M' needs one matrix per period.", call. = FALSE)
  W1 <- .psp_Wt(W, 1)
  N <- nrow(W1)
  ids <- if (is.null(rownames(W1))) seq_len(N) else rownames(W1)
  k <- length(beta) - 1L
  if (!is.null(theta) && length(theta) != k) {
    stop("'theta' needs one entry per regressor.", call. = FALSE)
  }
  if (!is.null(cutpoints) && (any(diff(c(0, cutpoints)) <= 0))) {
    stop("'cutpoints' must be positive and increasing.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(stats::rnorm(N * T * k), N * T, k, dimnames = list(NULL, paste0("x", seq_len(k))))
  alpha <- stats::rnorm(N, 0, sigma_alpha)
  if (phi != 0) {
    if (is.list(W)) stop("'phi' needs a single 'W'.", call. = FALSE)
    alpha <- as.vector(Matrix::solve(Matrix::Diagonal(N) - phi * W, alpha))
  }
  I <- Matrix::Diagonal(N)
  z <- numeric(N * T)
  y <- integer(N * T)
  ncat <- length(cutpoints) + 2L
  lag_cols <- function(yt) {
    if (is.null(cutpoints)) matrix(yt, ncol = 1) else outer(yt, 2:ncat, "==") * 1
  }
  for (t in seq_len(T)) {
    rows <- (t - 1) * N + seq_len(N)
    Xt <- X[rows, , drop = FALSE]
    Wt <- .psp_Wt(W, t)
    eps <- if (is.null(df)) stats::rnorm(N) else stats::rt(N, df)
    u <- if (lambda == 0) eps else as.vector(Matrix::solve(I - lambda * .psp_Wt(M, t), eps))
    mu <- beta[1] + Xt %*% beta[-1] + alpha + u
    if (!is.null(theta)) mu <- mu + as.matrix(Wt %*% Xt) %*% theta
    if (t > 1 && (!is.null(gamma) || !is.null(eta))) {
      Yl <- lag_cols(y[rows - N])
      if (!is.null(gamma)) mu <- mu + Yl %*% gamma
      if (!is.null(eta)) mu <- mu + as.matrix(Wt %*% Yl) %*% eta
    }
    if (t > 1 && gamma_z != 0) mu <- mu + gamma_z * z[rows - N]
    z[rows] <- as.vector(Matrix::solve(I - rho * Wt, mu))
    y[rows] <- if (is.null(cutpoints)) as.integer(z[rows] > 0) else
      as.integer(cut(z[rows], c(-Inf, 0, cutpoints, Inf), labels = FALSE))
  }
  d <- data.frame(id = rep(ids, T), time = rep(seq_len(T), each = N), y = y, X)
  attr(d, "truth") <- list(beta = beta, rho = rho, lambda = lambda, sigma_alpha = sigma_alpha,
                           phi = phi, df = df, gamma_z = gamma_z,
                           theta = theta, gamma = gamma, eta = eta,
                           cutpoints = cutpoints, alpha = alpha, z = z)
  attr(d, "W") <- W
  d
}
