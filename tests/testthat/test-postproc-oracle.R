# The Rust post-processing (pspcore::lndet, rf, post) against the R code it
# replaced, kept in R/reference_postproc.R (plan/06 §6.10 layer 5): the same
# numbers to 1e-10.

oracle_eq <- function(a, b, tol = 1e-10, info = NULL) {
  expect_equal(as.vector(a), as.vector(b), tolerance = tol, info = info)
}

rf_cases <- function() {
  W <- psp_weights(rook_weights(4, 5), style = "W")
  M <- psp_weights(rook_weights(5, 4), style = "W")
  dimnames(M) <- dimnames(W)
  set.seed(3)
  d <- sim_psprobit(W, T = 3, beta = c(0.2, 1, -0.5), rho = 0.4, sigma_alpha = 0.6, seed = 3)
  ctl <- function(...) psp_control(iter = 60, warmup = 40, chains = 1, seed = 4, ...)
  fit <- function(model = "sar", ..., ctl_args = list()) {
    suppressWarnings(psprobit(y ~ x1 + x2, d, c("id", "time"), W, model = model,
                              impacts = FALSE, control = do.call(ctl, ctl_args), ...))
  }
  Wt <- list(W, W, psp_weights(rook_weights(5, 4), style = "W"))
  dimnames(Wt[[3]]) <- dimnames(W)
  list(
    sar_re = fit("sar", effects = "random"),
    sdm = fit("sdm", effects = "pooled"),
    sem_re = fit("sem", effects = "random", M = M),
    sarar = fit("sarar", M = M, effects = "pooled"),
    spatial_re = fit("sar", effects = "spatial_random"),
    t_errors = fit("sar", errors = psp_t(5), effects = "pooled"),
    sarar_t = fit("sarar", M = M, errors = psp_t(5), effects = "random"),
    tv_weights = suppressWarnings(psprobit(y ~ x1 + x2, d, c("id", "time"), Wt, impacts = FALSE,
                                           control = ctl())),
    latent_ar = fit("sar", effects = "random", dynamic = psp_dynamic(latent = TRUE)),
    latent_ar_sem = fit("sem", M = M, dynamic = psp_dynamic(latent = TRUE), errors = psp_t(5))
  )
}

test_that("log-determinant grids match the R code", {
  W <- psp_weights(rook_weights(6, 7), style = "W")
  for (method in c("eigen", "sparse")) {
    a <- panelspatialprobit:::.psp_lndet(W, -1, 1, method)
    b <- panelspatialprobit:::.psp_ref_lndet(W, -1, 1, method)
    oracle_eq(a$lndet, b$lndet)
  }
})

test_that("reduced form, impacts and predict match the R code", {
  cases <- rf_cases()
  for (nm in names(cases)) {
    f <- cases[[nm]]
    pm <- panelspatialprobit:::.psp_par_matrix(f)
    mix <- panelspatialprobit:::.psp_t_mix(f)
    for (exact in c(TRUE, FALSE)) {
      probes <- if (exact) NULL else {
        P <- 7L * if (panelspatialprobit:::.psp_latent_ar(f)) f$T + 1L else 1L
        set.seed(9)
        matrix(sample(c(-1, 1), f$N * P, TRUE), f$N)
      }
      par <- pm[nrow(pm), ]
      a <- panelspatialprobit:::.psp_reduced_form(f, par, probes, need_impacts = TRUE, mix = mix)
      b <- panelspatialprobit:::.psp_ref_reduced_form(f, par, probes, need_impacts = TRUE, mix = mix)
      for (q in names(b)) {
        info <- sprintf("%s, %s, %s", nm, if (exact) "exact" else "probes", q)
        expect_identical(dim(a[[q]]), dim(b[[q]]), info = info)
        oracle_eq(a[[q]], b[[q]], info = info)
      }
      if (f$effects != "pooled") {
        al <- f$alpha_draws[1, ]
        a <- panelspatialprobit:::.psp_reduced_form(f, pm[f$alpha_index[1], ], probes, TRUE, alpha = al, mix = mix)
        b <- panelspatialprobit:::.psp_ref_reduced_form(f, pm[f$alpha_index[1], ], probes, TRUE, alpha = al, mix = mix)
        for (q in names(b)) oracle_eq(a[[q]], b[[q]], info = paste(nm, "alpha", q))
      }
    }
    f$control$impacts_exact_max <- if (nm == "sar_re") 0 else 1e4  # both paths
    oracle_eq(psp_impacts(f, draws = 5, by_period = TRUE)$draws,
              panelspatialprobit:::.psp_ref_impacts(f, draws = 5, by_period = TRUE))
    oracle_eq(psp_impacts(f, draws = 3, scale = "lesage_pace")$draws,
              panelspatialprobit:::.psp_ref_impacts(f, draws = 3, scale = "lesage_pace"))
    if (f$effects != "pooled") {
      oracle_eq(psp_impacts(f, draws = 3, type = "unit")$draws,
                panelspatialprobit:::.psp_ref_impacts(f, draws = 3, type = "unit"))
    }
    ref <- panelspatialprobit:::.psp_ref_predict(f, draws = 11)
    oracle_eq(predict(f, type = "link", draws = 11)[f$row_order], ref$link)
    oracle_eq(predict(f, type = "prob", draws = 11)[f$row_order], ref$prob[, 1])
  }
})

test_that("ordered models match the R code", {
  W <- psp_weights(rook_weights(4, 4), style = "W")
  d <- sim_psprobit(W, T = 3, beta = c(0.5, 1), rho = 0.3, cutpoints = c(0.7, 1.5), seed = 5)
  f <- suppressWarnings(psoprobit(y ~ x1, d, c("id", "time"), W, impacts = FALSE,
                                  effects = "random",
                                  control = psp_control(iter = 60, warmup = 40, chains = 1,
                                                        seed = 5)))
  oracle_eq(psp_impacts(f, draws = 5, by_period = TRUE)$draws,
            panelspatialprobit:::.psp_ref_impacts(f, draws = 5, by_period = TRUE))
  ref <- panelspatialprobit:::.psp_ref_predict(f, draws = 9)
  oracle_eq(predict(f, type = "prob", draws = 9)[f$row_order, ], ref$prob)
})
