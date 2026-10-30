fd_impacts <- function(fit, par, var, h = 1e-5, alpha = NULL) {
  lambda <- if ("lambda" %in% names(par)) par[["lambda"]] else 0
  N <- fit$N
  T <- fit$T
  Wt <- function(t) as.matrix(if (is.list(fit$W)) fit$W[[t]] else fit$W)
  Mt <- function(t) {
    if (is.null(fit$M)) return(Wt(t))
    as.matrix(if (is.list(fit$M)) fit$M[[t]] else fit$M)
  }
  beta <- par[colnames(fit$X)]
  rho <- if ("rho" %in% names(par)) par[["rho"]] else 0
  sig2 <- if ("sigma_alpha" %in% names(par)) par[["sigma_alpha"]]^2 else 0
  if (!is.null(alpha)) sig2 <- 0 else alpha <- numeric(N)
  phi <- if ("phi" %in% names(par)) par[["phi"]] else 0
  wcol <- paste0("W_", var)
  prob_t <- function(X, t) {
    Si <- solve(diag(N) - rho * Wt(t))
    K <- Si %*% solve(diag(N) - lambda * Mt(t))
    G <- Si %*% solve(diag(N) - phi * Wt(t))
    s <- sqrt(sig2 * rowSums(G^2) + rowSums(K^2))
    rows <- (t - 1) * N + seq_len(N)
    pnorm(as.vector(Si %*% (X[rows, ] %*% beta + alpha)) / s)
  }
  shift <- function(X, t, j, eps) {
    rows <- (t - 1) * N + seq_len(N)
    e <- if (is.null(j)) rep(eps, N) else replace(numeric(N), j, eps)
    X[rows, var] <- X[rows, var] + e
    if (wcol %in% colnames(X)) X[rows, wcol] <- X[rows, wcol] + as.vector(Wt(t) %*% e)
    X
  }
  direct <- total <- 0
  for (t in seq_len(T)) {
    for (j in seq_len(N)) {
      dp <- (prob_t(shift(fit$X, t, j, h), t) - prob_t(shift(fit$X, t, j, -h), t)) / (2 * h)
      direct <- direct + unname(dp[j])
    }
    dp <- (prob_t(shift(fit$X, t, NULL, h), t) - prob_t(shift(fit$X, t, NULL, -h), t)) / (2 * h)
    total <- total + sum(dp)
  }
  c(direct = direct / (N * T), total = total / (N * T))
}

test_that("impacts match finite differences of exact probabilities (SDM, RE)", {
  s <- small_data(theta = c(0.5, 0))
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, model = "sdm",
                                   impacts = FALSE, control = fast_control()))
  par <- panelspatialprobit:::.psp_par_matrix(fit)[1, ]
  imp <- psp_impacts(fit, draws = 1)
  for (v in c("x1", "x2")) {
    fd <- fd_impacts(fit, par, v)
    expect_equal(imp$draws[1, v, 1, "direct"], unname(fd["direct"]), tolerance = 1e-6)
    expect_equal(imp$draws[1, v, 1, "total"], unname(fd["total"]), tolerance = 1e-6)
  }
})

test_that("impacts for SLX and non-spatial models", {
  s <- small_data(rho = 0)
  for (mod in c("slx", "none")) {
    fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, model = mod,
                                     impacts = FALSE, control = fast_control()))
    par <- panelspatialprobit:::.psp_par_matrix(fit)[1, ]
    imp <- psp_impacts(fit, draws = 1)
    fd <- fd_impacts(fit, par, "x1")
    expect_equal(imp$draws[1, "x1", 1, "direct"], unname(fd["direct"]), tolerance = 1e-6)
    expect_equal(imp$draws[1, "x1", 1, "total"], unname(fd["total"]), tolerance = 1e-6)
  }
})

test_that("stochastic diagonal estimator approximates the exact computation", {
  s <- small_data()
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, impacts = FALSE,
                                   control = fast_control()))
  exact <- psp_impacts(fit, draws = 3)
  fit$control$impacts_exact_max <- 0
  fit$control$impacts_probes <- 4000L
  approx <- psp_impacts(fit, draws = 3)
  expect_equal(approx$draws, exact$draws, tolerance = 0.05)
})

test_that("ordered impacts sum to zero across categories", {
  s <- small_data(cutpoints = c(0.7, 1.5))
  fit <- suppressWarnings(psoprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, impacts = FALSE,
                                    control = fast_control()))
  imp <- psp_impacts(fit, draws = 5)
  sums <- apply(imp$draws, c(1, 2, 4), sum)
  expect_true(all(abs(sums) < 1e-10))
  expect_equal(dim(imp$draws)[3], 4)
})

test_that("unit-specific impacts match finite differences given alpha", {
  s <- small_data(seed = 5, theta = c(0.3, 0))
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, model = "sdm",
                                   impacts = FALSE, control = fast_control()))
  expect_equal(nrow(fit$alpha_draws), length(fit$alpha_index))
  imp <- psp_impacts(fit, draws = 1, type = "unit")
  d1 <- fit$alpha_index[1]
  par <- panelspatialprobit:::.psp_par_matrix(fit)[d1, ]
  fd <- fd_impacts(fit, par, "x1", alpha = fit$alpha_draws[1, ])
  expect_equal(imp$draws[1, "x1", 1, "direct"], unname(fd["direct"]), tolerance = 1e-6)
  expect_equal(imp$draws[1, "x1", 1, "total"], unname(fd["total"]), tolerance = 1e-6)
  pooled <- suppressWarnings(psprobit(y ~ x1, s$d, c("id", "time"), s$W, effects = "pooled",
                                      impacts = FALSE, control = fast_control()))
  expect_error(psp_impacts(pooled, type = "unit"), "random-effects")
})

test_that("per-period impacts average to the overall impacts", {
  s <- small_data(seed = 6)
  fit <- suppressWarnings(psoprobit(y ~ x1 + x2, small_data(seed = 6, cutpoints = 0.8)$d,
                                    c("id", "time"), s$W, impacts = FALSE,
                                    control = fast_control()))
  imp <- psp_impacts(fit, draws = 4, by_period = TRUE)
  expect_equal(dim(imp$draws), c(4, 2, 3, 3, 1 + fit$T))
  per <- apply(imp$draws[, , , , -1, drop = FALSE], 1:4, mean)
  expect_equal(per, imp$draws[, , , , 1])
  sm <- summary(imp)
  expect_true(all(c("period", "category") %in% names(sm)))
  expect_output(print(imp), "per-period")
})

test_that("impacts with a spatial error term match finite differences", {
  W <- psp_weights(rook_weights(4, 4), style = "W")
  d <- sim_psprobit(W, T = 3, beta = c(0, 1, -0.5), rho = 0.3, lambda = 0.4, theta = c(0.4, 0),
                    sigma_alpha = 0.5, seed = 4)
  for (mod in c("sem", "sdem", "sarar")) {
    fit <- suppressWarnings(psprobit(y ~ x1 + x2, d, c("id", "time"), W, model = mod,
                                     impacts = FALSE, control = fast_control()))
    par <- panelspatialprobit:::.psp_par_matrix(fit)[1, ]
    imp <- psp_impacts(fit, draws = 1)
    fd <- fd_impacts(fit, par, "x1")
    expect_equal(imp$draws[1, "x1", 1, "direct"], unname(fd["direct"]), tolerance = 1e-6)
    expect_equal(imp$draws[1, "x1", 1, "total"], unname(fd["total"]), tolerance = 1e-6)
    if (mod == "sem") expect_equal(imp$draws[1, "x1", 1, "indirect"], 0)
  }
  # stochastic diagonals with an error term
  fit$control$impacts_exact_max <- 0
  fit$control$impacts_probes <- 4000L
  approx <- psp_impacts(fit, draws = 1)
  expect_equal(approx$draws[1, "x1", 1, "total"], imp$draws[1, "x1", 1, "total"],
               tolerance = 0.03)
})

test_that("impacts with time-varying weights match finite differences", {
  s <- tv_data(seed = 8)
  for (mod in c("sdm", "sarar")) {
    fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, model = mod,
                                     impacts = FALSE, control = fast_control()))
    expect_true(is.list(fit$W))
    par <- panelspatialprobit:::.psp_par_matrix(fit)[1, ]
    imp <- psp_impacts(fit, draws = 1)
    fd <- fd_impacts(fit, par, "x1")
    expect_equal(imp$draws[1, "x1", 1, "direct"], unname(fd["direct"]), tolerance = 1e-6)
    expect_equal(imp$draws[1, "x1", 1, "total"], unname(fd["total"]), tolerance = 1e-6)
  }
})

test_that("impacts with spatial random effects match finite differences", {
  s <- small_data(seed = 12)
  for (mod in c("sdm", "slx")) {
    fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, model = mod,
                                     effects = "spatial_random", impacts = FALSE,
                                     control = fast_control()))
    par <- panelspatialprobit:::.psp_par_matrix(fit)[1, ]
    expect_true("phi" %in% names(par))
    imp <- psp_impacts(fit, draws = 1)
    fd <- fd_impacts(fit, par, "x1")
    expect_equal(imp$draws[1, "x1", 1, "direct"], unname(fd["direct"]), tolerance = 1e-6)
    expect_equal(imp$draws[1, "x1", 1, "total"], unname(fd["total"]), tolerance = 1e-6)
  }
  # stochastic diagonals agree with the exact ones
  fit$control$impacts_exact_max <- 1
  fit$control$impacts_probes <- 2000L
  par <- panelspatialprobit:::.psp_par_matrix(fit)[1, ]
  ex <- panelspatialprobit:::.psp_reduced_form(fit, par, NULL)
  st <- panelspatialprobit:::.psp_reduced_form(fit, par, panelspatialprobit:::.psp_probes(fit))
  expect_equal(st$s, ex$s, tolerance = 0.05)
})
