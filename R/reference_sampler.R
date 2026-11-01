# Pure-R reference sampler (dense, slow): the test oracle for the Rust engine.
#
# Deliberately uses different algorithms where possible: the latent sweep
# conditions through the dense precision matrix and mean, and rho is drawn by
# random-walk Metropolis with exact log-determinants (beta not integrated out).
# The spatial error term (lambda, M) enters through A = R S, also by Metropolis.
# W and M may be lists with one matrix per period (time-varying weights).
# nu > 0 gives Student-t errors through precision weights w (scale mixture).
# lag_col (an index of X, whose column is ignored) adds latent AR(1) dynamics
# gamma z_{t-1}: the latent conditionals then use the full NT x NT operator and
# gamma is drawn by random-walk Metropolis.
# Same posterior, different route. Only for small problems.
.psp_reference_sampler <- function(y, X, W, N, T, ordered = FALSE, J = 2L, has_rho = TRUE,
                                   has_re = TRUE, has_lambda = FALSE, M = W,
                                   spatial_re = FALSE, nu = 0, lag_col = NULL,
                                   beta_sd = 10, rho_ab = c(1.01, 1.01),
                                   sigma_scale = 1, kappa_sd = 10, iter = 2000,
                                   warmup = 1000, rw_sd = 0.1) {
  tv <- is.list(W) || is.list(M)  # time-varying weights
  latent <- !is.null(lag_col)
  if (latent) X <- X[, -lag_col, drop = FALSE]
  gamma <- 0
  lagz <- function(z) c(rep(0, N), z[seq_len(N * (T - 1))])
  Wl <- lapply(if (is.list(W)) W else rep(list(W), T), as.matrix)
  Ml <- lapply(if (is.list(M)) M else rep(list(M), T), as.matrix)
  k <- ncol(X)
  n <- N * T
  I <- diag(N)
  Pb <- diag(1 / beta_sd^2, k)
  beta <- numeric(k)
  rho <- 0
  lambda <- 0
  phi <- 0
  sig2 <- 1
  xi <- 1
  alpha <- numeric(N)
  w <- rep(1, N * T)
  kappa <- if (ordered) c(-Inf, 0, seq_len(J - 2), Inf) else c(-Inf, 0, Inf)
  yb <- if (ordered) y else y + 1L
  z <- ifelse(yb == 1, -0.5, kappa[yb] + 0.5)
  Ss <- function(rho_) lapply(Wl, function(w) I - rho_ * w)
  Rs <- function(lambda_) lapply(Ml, function(m) I - lambda_ * m)
  # per-period map v -> Op_t v_t on stacked vectors
  pmul <- function(Ops, v) {
    V <- matrix(v, N)
    if (!tv) return(as.vector(Ops[[1]] %*% V))
    as.vector(vapply(seq_len(T), function(t) as.vector(Ops[[t]] %*% V[, t]), numeric(N)))
  }
  # log posterior of rho (or lambda) given everything else
  logpost <- function(r, rho_, lambda_) {
    if (r <= -1 || r >= 1) return(-Inf)
    S <- Ss(rho_)
    R <- Rs(lambda_)
    ld <- (if (tv) 1 else T) * sum(vapply(if (tv) seq_len(T) else 1L, function(t) {
      as.numeric(determinant(S[[t]], logarithm = TRUE)$modulus) +
        as.numeric(determinant(R[[t]], logarithm = TRUE)$modulus)
    }, numeric(1)))
    e <- pmul(R, pmul(S, z) - gamma * lagz(z) - X %*% beta - rep(alpha, T))
    ld - 0.5 * sum(w * e^2) + (rho_ab[1] - 1) * log(1 + r) + (rho_ab[2] - 1) * log(1 - r)
  }
  # spatial random effects: log density of alpha given phi (up to constants)
  Bmat <- function(phi_) I - phi_ * Wl[[1]]
  logpost_phi <- function(phi_) {
    if (phi_ <= -1 || phi_ >= 1) return(-Inf)
    B <- Bmat(phi_)
    as.numeric(determinant(B, logarithm = TRUE)$modulus) - 0.5 * sum((B %*% alpha)^2) / sig2 +
      (rho_ab[1] - 1) * log(1 + phi_) + (rho_ab[2] - 1) * log(1 - phi_)
  }
  logpost_gamma <- function(g) {
    e <- pmul(Rs(lambda), pmul(Ss(rho), z) - g * lagz(z) - X %*% beta - rep(alpha, T))
    -0.5 * sum(w * e^2) - 0.5 * g^2 / beta_sd^2
  }
  out <- matrix(NA_real_, iter, latent + k + has_rho + has_lambda + has_re + spatial_re +
                  (if (ordered) J - 2 else 0))
  for (it in seq_len(warmup + iter)) {
    S <- Ss(rho)
    R <- Rs(lambda)
    if (latent) {
      # full operator A (block lower-bidiagonal) and mean A^{-1} b
      Af <- matrix(0, N * T, N * T)
      for (t in seq_len(T)) {
        rows <- (t - 1) * N + seq_len(N)
        Af[rows, rows] <- R[[t]] %*% S[[t]]
        if (t > 1) Af[rows, rows - N] <- -gamma * R[[t]]
      }
      b <- pmul(R, X %*% beta + rep(alpha, T))
      muf <- solve(Af, b)
      H <- crossprod(sqrt(w) * Af)
      for (i in seq_len(N * T)) {
        cm <- muf[i] - sum(H[i, -i] * (z[-i] - muf[-i])) / H[i, i]
        lo <- kappa[yb[i]]
        hi <- kappa[yb[i] + 1]
        p <- stats::pnorm(c(lo, hi), cm, 1 / sqrt(H[i, i]))
        u <- stats::runif(1, p[1], p[2])
        z[i] <- stats::qnorm(min(max(u, 1e-300), 1 - 1e-16), cm, 1 / sqrt(H[i, i]))
      }
    }
    for (t in if (latent) integer(0) else seq_len(T)) {
      rows <- (t - 1) * N + seq_len(N)
      if (tv || t == 1 || nu > 0) {
        Si <- solve(S[[t]])
        H <- crossprod(sqrt(w[rows]) * (R[[t]] %*% S[[t]]))
      }
      mu <- Si %*% (X[rows, , drop = FALSE] %*% beta + alpha)
      zt <- z[rows]
      for (i in seq_len(N)) {
        cm <- mu[i] - sum(H[i, -i] * (zt[-i] - mu[-i])) / H[i, i]
        lo <- kappa[yb[rows[i]]]
        hi <- kappa[yb[rows[i]] + 1]
        p <- stats::pnorm(c(lo, hi), cm, 1 / sqrt(H[i, i]))
        u <- stats::runif(1, p[1], p[2])
        zt[i] <- stats::qnorm(min(max(u, 1e-300), 1 - 1e-16), cm, 1 / sqrt(H[i, i]))
      }
      z[rows] <- zt
    }
    # beta | rho, lambda, z, alpha
    ystar <- pmul(R, pmul(S, z) - gamma * lagz(z) - rep(alpha, T))
    Xs <- apply(X, 2, function(x) pmul(R, x))
    V <- solve(crossprod(Xs, w * Xs) + Pb)
    beta <- as.vector(V %*% crossprod(Xs, w * ystar) + t(chol(V)) %*% stats::rnorm(k))
    # rho and lambda | beta, z, alpha by random-walk Metropolis
    if (latent) {
      prop <- gamma + rw_sd * stats::rnorm(1)
      if (log(stats::runif(1)) < logpost_gamma(prop) - logpost_gamma(gamma)) gamma <- prop
    }
    if (has_rho) {
      prop <- rho + rw_sd * stats::rnorm(1)
      if (log(stats::runif(1)) < logpost(prop, prop, lambda) - logpost(rho, rho, lambda)) {
        rho <- prop
      }
    }
    if (has_lambda) {
      prop <- lambda + rw_sd * stats::rnorm(1)
      if (log(stats::runif(1)) < logpost(prop, rho, prop) - logpost(lambda, rho, lambda)) {
        lambda <- prop
      }
    }
    S <- Ss(rho)
    R <- Rs(lambda)
    if (has_re) {
      e <- matrix(pmul(S, z) - gamma * lagz(z) - X %*% beta, N, T)
      RtR <- lapply(seq_len(T), function(t) crossprod(R[[t]], w[(t - 1) * N + seq_len(N)] * R[[t]]))
      Pa <- Reduce(`+`, RtR) + crossprod(Bmat(phi)) / sig2
      Va <- solve(Pa)
      b <- Reduce(`+`, lapply(seq_len(T), function(t) RtR[[t]] %*% e[, t]))
      alpha <- as.vector(Va %*% b + t(chol(Va)) %*% stats::rnorm(N))
      if (spatial_re) {
        prop <- phi + rw_sd * stats::rnorm(1)
        if (log(stats::runif(1)) < logpost_phi(prop) - logpost_phi(phi)) phi <- prop
      }
      sig2 <- 1 / stats::rgamma(1, (N + 1) / 2, sum((Bmat(phi) %*% alpha)^2) / 2 + 1 / xi)
      xi <- 1 / stats::rgamma(1, 1, 1 / sig2 + 1 / sigma_scale^2)
    }
    if (nu > 0) {
      r <- pmul(R, pmul(S, z) - gamma * lagz(z) - X %*% beta - rep(alpha, T))
      w <- stats::rgamma(N * T, (nu + 1) / 2, rate = (nu + r^2) / 2)
    }
    if (ordered) {
      for (j in 2:(J - 1)) {
        lo <- max(c(z[yb == j], kappa[j]))
        hi <- min(c(z[yb == j + 1], kappa[j + 2]))
        p <- stats::pnorm(c(lo, hi), 0, kappa_sd)
        kappa[j + 1] <- stats::qnorm(stats::runif(1, p[1], p[2]), 0, kappa_sd)
      }
    }
    if (it > warmup) {
      bfull <- if (latent) append(beta, gamma, after = lag_col - 1L) else beta
      out[it - warmup, ] <- c(bfull, if (has_rho) rho, if (has_lambda) lambda,
                              if (has_re) sqrt(sig2), if (spatial_re) phi,
                              if (ordered) kappa[3:J])
    }
  }
  out
}
