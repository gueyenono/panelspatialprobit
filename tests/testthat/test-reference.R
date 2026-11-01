# Distributional agreement between the Rust engine and the pure-R oracle,
# which uses a different latent update and a Metropolis step for rho.
compare_to_reference <- function(s, ordered = FALSE, J = 2L, model = "sar", effects = "random",
                                 df = NULL, latent = FALSE) {
  pri <- psp_priors(beta = prior_normal(0, 2), cutpoints = prior_normal(0, 2))
  f <- if (ordered) psoprobit else psprobit
  fit <- suppressWarnings(f(y ~ x1 + x2, s$d, c("id", "time"), s$W, priors = pri,
                            impacts = FALSE, model = model, effects = effects,
                            errors = if (is.null(df)) "normal" else psp_t(df),
                            dynamic = if (latent) psp_dynamic(latent = TRUE),
                            control = psp_control(iter = 4000, warmup = 1000, chains = 2,
                                                  seed = 21)))
  des <- suppressWarnings(panelspatialprobit:::.psp_design(
    y ~ x1 + x2, s$d, c("id", "time"), s$W, model, NULL, NULL, FALSE,
    if (ordered) "ordered" else "binary", NULL, if (latent) psp_dynamic(latent = TRUE)))
  set.seed(22)
  ref <- panelspatialprobit:::.psp_reference_sampler(des$y, des$X, des$W, des$N, des$T,
                                                     ordered = ordered, J = des$J, beta_sd = 2,
                                                     has_rho = model %in% c("sar", "sarar"),
                                                     has_lambda = model %in% c("sem", "sarar"),
                                                     M = des$W,
                                                     spatial_re = effects == "spatial_random",
                                                     nu = if (is.null(df)) 0 else df,
                                                     lag_col = if (latent) match("lag_z", colnames(des$X)),
                                                     kappa_sd = 2, iter = 6000, warmup = 1000)
  colnames(ref) <- posterior::variables(fit$draws)
  ref_d <- posterior::as_draws_array(array(ref, c(nrow(ref), 1, ncol(ref)),
                                           dimnames = list(NULL, NULL, colnames(ref))))
  m1 <- posterior::summarise_draws(fit$draws, "mean", "mcse_mean")
  m2 <- posterior::summarise_draws(ref_d, "mean", "mcse_mean")
  z <- (m1$mean - m2$mean) / sqrt(m1$mcse_mean^2 + m2$mcse_mean^2)
  expect_true(all(abs(z) < 4), info = paste(m1$variable, round(z, 2), collapse = "; "))
}

test_that("Sampler agrees with the R reference sampler (binary SAR, RE)", {
  skip_on_cran()
  compare_to_reference(small_data(seed = 4, T = 4))
})

test_that("Sampler agrees with the R reference sampler (ordered SAR, RE)", {
  skip_on_cran()
  compare_to_reference(small_data(seed = 5, T = 4, cutpoints = 0.8), ordered = TRUE)
})

test_that("Sampler agrees with the R reference sampler (binary SEM and SARAR, RE)", {
  skip_on_cran()
  W <- psp_weights(rook_weights(4, 4), style = "W")
  d <- sim_psprobit(W, T = 4, beta = c(0, 1, -0.5), rho = 0, lambda = 0.5,
                    sigma_alpha = 0.5, seed = 8)
  compare_to_reference(list(d = d, W = W), model = "sem")
  d <- sim_psprobit(W, T = 4, beta = c(0, 1, -0.5), rho = 0.4, lambda = -0.3,
                    sigma_alpha = 0.5, seed = 9)
  compare_to_reference(list(d = d, W = W), model = "sarar")
})

test_that("Sampler agrees with the R reference sampler (time-varying W)", {
  skip_on_cran()
  compare_to_reference(tv_data(seed = 10, T = 4))
  compare_to_reference(tv_data(seed = 11, T = 4, rho = 0.3, lambda = -0.4), model = "sarar")
})

test_that("Sampler agrees with the R reference sampler (spatial random effects)", {
  skip_on_cran()
  W <- psp_weights(rook_weights(4, 4), style = "W")
  d <- sim_psprobit(W, T = 4, beta = c(0, 1, -0.5), rho = 0.3, sigma_alpha = 0.8, phi = 0.5,
                    seed = 13)
  compare_to_reference(list(d = d, W = W), effects = "spatial_random")
  d <- sim_psprobit(W, T = 4, beta = c(0, 1, -0.5), rho = 0, lambda = 0.4, sigma_alpha = 0.8,
                    phi = 0.5, seed = 14)
  compare_to_reference(list(d = d, W = W), model = "sem", effects = "spatial_random")
})

test_that("Sampler agrees with the R reference sampler (Student-t errors)", {
  skip_on_cran()
  W <- psp_weights(rook_weights(4, 4), style = "W")
  d <- sim_psprobit(W, T = 4, beta = c(0, 1, -0.5), rho = 0.4, sigma_alpha = 0.5, df = 4,
                    seed = 15)
  compare_to_reference(list(d = d, W = W), df = 4)
  d <- sim_psprobit(W, T = 4, beta = c(0, 1, -0.5), rho = 0.3, lambda = 0.4,
                    sigma_alpha = 0.5, cutpoints = 0.8, df = 4, seed = 16)
  compare_to_reference(list(d = d, W = W), ordered = TRUE, model = "sarar", df = 4)
})

test_that("Sampler agrees with the R reference sampler (latent AR(1))", {
  skip_on_cran()
  W <- psp_weights(rook_weights(4, 4), style = "W")
  d <- sim_psprobit(W, T = 4, beta = c(0, 1, -0.5), rho = 0.3, sigma_alpha = 0.5,
                    gamma_z = 0.5, seed = 17)
  compare_to_reference(list(d = d, W = W), latent = TRUE)
  d <- sim_psprobit(W, T = 4, beta = c(0, 1, -0.5), rho = 0.2, lambda = 0.4, sigma_alpha = 0.5,
                    gamma_z = 0.5, cutpoints = 0.8, df = 5, seed = 18)
  compare_to_reference(list(d = d, W = W), ordered = TRUE, model = "sarar", df = 5,
                       latent = TRUE)
})
