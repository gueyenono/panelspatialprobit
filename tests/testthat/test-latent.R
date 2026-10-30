latent_fit <- function(seed = 30, ordered = FALSE, model = "sar", effects = "random", ...) {
  W <- psp_weights(rook_weights(4, 4), style = "W")
  d <- sim_psprobit(W, T = 4, beta = c(0, 1, -0.5), rho = 0.3, sigma_alpha = 0.5,
                    gamma_z = 0.6, cutpoints = if (ordered) 0.8, seed = seed,
                    lambda = if (model == "sarar") 0.3 else 0)
  f <- if (ordered) psoprobit else psprobit
  suppressWarnings(f(y ~ x1 + x2, d, c("id", "time"), W, model = model, effects = effects,
                     dynamic = psp_dynamic(latent = TRUE), impacts = FALSE,
                     control = fast_control(...)))
}

# exact latent mean and variance by brute force on the stacked system
brute_rf <- function(fit, par) {
  N <- fit$N
  T <- fit$T
  W <- as.matrix(fit$W)
  M <- if (is.null(fit$M)) W else as.matrix(fit$M)
  rho <- if ("rho" %in% names(par)) par[["rho"]] else 0
  lambda <- if ("lambda" %in% names(par)) par[["lambda"]] else 0
  g <- par[["lag_z"]]
  S <- diag(N) - rho * W
  R <- diag(N) - lambda * M
  A <- kronecker(diag(T), R %*% S)
  for (t in 2:T) A[(t - 1) * N + 1:N, (t - 2) * N + 1:N] <- -g * R
  Ai <- solve(A)
  Rb <- kronecker(diag(T), R)
  m <- Ai %*% Rb %*% (fit$X %*% par[colnames(fit$X)])
  D <- Ai %*% Rb %*% kronecker(rep(1, T), diag(N))  # loading of alpha
  sig2 <- if ("sigma_alpha" %in% names(par)) par[["sigma_alpha"]]^2 else 0
  V <- tcrossprod(Ai) + sig2 * tcrossprod(D)
  list(m = as.vector(m), s = sqrt(diag(V)))
}

test_that("latent AR design, reduced form and simulation agree", {
  fit <- latent_fit()
  expect_true("lag_z" %in% colnames(fit$X))
  expect_equal(unname(fit$roles["lag_z"]), "latent")
  expect_false("lag_z" %in% names(fit$roles)[fit$roles == "x"])
  par <- panelspatialprobit:::.psp_par_matrix(fit)[1, ]
  ex <- panelspatialprobit:::.psp_reduced_form(fit, par)
  bf <- brute_rf(fit, par)
  expect_equal(ex$m, bf$m, tolerance = 1e-8)
  expect_equal(ex$s, bf$s, tolerance = 1e-8)
  fit$control$impacts_exact_max <- 1
  fit$control$impacts_probes <- 4000L
  st <- panelspatialprobit:::.psp_reduced_form(fit, par, panelspatialprobit:::.psp_probes(fit))
  expect_equal(st$s, bf$s, tolerance = 0.05)
  f2 <- latent_fit(31, ordered = TRUE, model = "sarar")
  par <- panelspatialprobit:::.psp_par_matrix(f2)[1, ]
  ex <- panelspatialprobit:::.psp_reduced_form(f2, par)
  bf <- brute_rf(f2, par)
  expect_equal(ex$s, bf$s, tolerance = 1e-8)
})

test_that("latent AR impacts and response paths", {
  fit <- latent_fit(32)
  fit$draws <- posterior::subset_draws(fit$draws, chain = 1, iteration = 1)
  imp <- psp_impacts(fit)
  h <- 1e-5
  shift <- function(eps) {
    f <- fit
    f$X[, "x1"] <- f$X[, "x1"] + eps
    mean(predict(f))
  }
  # a shift in every period also moves later periods through gamma, so compare
  # the within-period impact with a shift of one period at a time
  shift_t <- function(eps) {
    tot <- 0
    for (t in 1:4) {
      f <- fit
      rows <- (t - 1) * 16 + 1:16
      f$X[rows, "x1"] <- f$X[rows, "x1"] + eps
      p <- predict(f, draws = 1)
      p[fit$row_order] <- p
      tot <- tot + mean(p[rows])
    }
    tot / 4
  }
  expect_equal(imp$draws[1, "x1", 1, "total"], (shift_t(h) - shift_t(-h)) / (2 * h),
               tolerance = 1e-5)
  fit <- latent_fit(33)
  r <- psp_response_path(fit, "x1", draws = 20)
  expect_equal(r$summary$horizon, 0:3)
  expect_true(all(r$summary$mean[2:4] != 0))
  rp <- psp_response_path(fit, "x1", type = "permanent", start = 2, draws = 20)
  expect_equal(nrow(rp$summary), 3)
  expect_s3_class(psp_ppc(fit, draws = 10, seed = 1), "psp_ppc")
})

test_that("latent AR GHK agrees with mvtnorm", {
  skip_if_not_installed("mvtnorm")
  fit <- latent_fit(34, effects = "pooled")
  par <- fit$coefficients
  bf <- brute_rf(fit, par)
  N <- 16
  W <- as.matrix(fit$W)
  S <- diag(N) - par[["rho"]] * W
  A <- kronecker(diag(4), S)
  for (t in 2:4) A[(t - 1) * N + 1:N, (t - 2) * N + 1:N] <- -par[["lag_z"]] * diag(N)
  set.seed(1)
  ex <- log(mvtnorm::pmvnorm(lower = ifelse(fit$y == 1, 0, -Inf),
                             upper = ifelse(fit$y == 1, Inf, 0), mean = bf$m,
                             sigma = tcrossprod(solve(A)),
                             algorithm = mvtnorm::GenzBretz(maxpts = 5e5, abseps = 1e-8)))
  ll <- suppressWarnings(logLik(fit, reps = 5000))
  expect_lt(abs(as.numeric(ll) - as.numeric(ex)), 4 * attr(ll, "mcse") + 0.05)
})

test_that("latent AR argument checks", {
  expect_error(psp_dynamic(lag_y = TRUE, latent = TRUE), "cannot be combined")
  s <- tv_data(seed = 2, T = 3)
  expect_error(psprobit(y ~ x1, s$d, c("id", "time"), s$W,
                        dynamic = psp_dynamic(latent = TRUE)), "time-invariant")
})
