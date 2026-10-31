# Model cases for the cross-build checks (plan/06 §6.10 layer 3). Each case is
# a function of nothing that fits one small model; both builds run the same
# code with the same seed. `exact` cases must give the same draws (to 1e-8);
# the others (sparse random-effects block, whose Cholesky ordering differs)
# are compared by posterior means in compare_sampler.R.
psp_case_data <- function(N_side = 6, T = 5, ...) {
  W <- psp_weights(rook_weights(N_side, N_side), style = "W")
  list(W = W, d = sim_psprobit(W, T = T, beta = c(0, 1, -1), seed = 3, ...))
}

psp_cases <- function(iter = 40, warmup = 20, chains = 2) {
  ctl <- function(...) psp_control(iter = iter, warmup = warmup, chains = chains, seed = 11,
                                   progress = FALSE, ...)
  f <- y ~ x1 + x2
  ix <- c("id", "time")
  bin <- function(model, effects = "random", ..., data = psp_case_data(rho = 0.4)) {
    function() psprobit(f, data$d, ix, data$W, model = model, effects = effects,
                        impacts = FALSE, control = ctl(), ...)
  }
  ord_data <- psp_case_data(rho = 0.4, cutpoints = c(0.7, 1.5))
  ord <- function(model, effects = "random", ..., control = ctl()) {
    function() psoprobit(f, ord_data$d, ix, ord_data$W, model = model, effects = effects,
                         impacts = FALSE, control = control, ...)
  }
  tv <- local({
    W1 <- psp_weights(rook_weights(6, 6), style = "W")
    W2 <- psp_weights(Matrix::t(rook_weights(6, 6)) + Matrix::Diagonal(36) * 0, style = "W")
    W2 <- psp_weights(Matrix::sparseMatrix(i = c(1:35, 2:36), j = c(2:36, 1:35), x = 1,
                                           dims = c(36, 36),
                                           dimnames = dimnames(W1)), style = "W")
    Wl <- list(W1, W2, W1, W2, W2)
    list(W = Wl, d = sim_psprobit(Wl, T = 5, beta = c(0, 1, -1), rho = 0.4, seed = 3))
  })
  miss <- local({
    x <- psp_case_data(rho = 0.4)
    x$d$y[c(3, 17, 40, 41, 100)] <- NA
    x
  })
  lat <- psp_case_data(rho = 0.3, gamma_z = 0.4)
  sem <- psp_case_data(rho = 0, lambda = 0.4)
  list(
    sar_binary_re = list(exact = TRUE, fit = bin("sar")),
    sar_binary_pooled = list(exact = TRUE, fit = bin("sar", "pooled")),
    slx_binary_re = list(exact = TRUE, fit = bin("slx")),
    none_binary_re = list(exact = TRUE, fit = bin("none")),
    sdm_binary_re = list(exact = TRUE, fit = bin("sdm")),
    sar_binary_missing = list(exact = TRUE, fit = bin("sar", data = miss)),
    sar_ordered_re = list(exact = TRUE, fit = ord("sar")),
    sar_ordered_ac = list(exact = TRUE,
                          fit = ord("sar", control = ctl(cutpoint_sampler = "albert_chib"))),
    sar_t_re = list(exact = TRUE, fit = bin("sar", errors = psp_t(5))),
    sar_ordered_t = list(exact = TRUE, fit = ord("sar", errors = psp_t(5))),
    sem_pooled = list(exact = TRUE, fit = bin("sem", "pooled", data = sem)),
    sarar_pooled = list(exact = TRUE, fit = bin("sarar", "pooled", data = sem)),
    sdem_pooled = list(exact = TRUE, fit = bin("sdem", "pooled", data = sem)),
    sar_timevarying_re = list(exact = TRUE, fit = bin("sar", data = tv)),
    sar_latent_ar_re = list(exact = TRUE, fit = bin("sar", data = lat,
                                                    dynamic = psp_dynamic(latent = TRUE))),
    sem_latent_ar_pooled = list(exact = TRUE, fit = bin("sem", "pooled", data = lat,
                                                        dynamic = psp_dynamic(latent = TRUE))),
    sar_lag_y_re = list(exact = TRUE, fit = bin("sar", dynamic = psp_dynamic())),
    sem_re = list(exact = FALSE, fit = bin("sem", data = sem)),
    sarar_re = list(exact = FALSE, fit = bin("sarar", data = sem)),
    sdem_re = list(exact = FALSE, fit = bin("sdem", data = sem)),
    sem_re_t = list(exact = FALSE, fit = bin("sem", data = sem, errors = psp_t(5))),
    sem_re_timevarying = list(exact = FALSE, fit = bin("sem", data = tv, M = tv$W)),
    sem_latent_ar_re = list(exact = FALSE, fit = bin("sem", data = lat,
                                                     dynamic = psp_dynamic(latent = TRUE))),
    sar_spatial_re = list(exact = FALSE, fit = bin("sar", "spatial_random")),
    sar_spatial_re_t = list(exact = FALSE, fit = bin("sar", "spatial_random",
                                                     errors = psp_t(5)))
  )
}
