# Results do not depend on how the work is split (plan/06 §6.10, R-D9):
# chains over threads (`cores`), and post-processing draws over the workers
# of a future::plan() (R-D11).

test_that("post-processing gives the same results under any future plan", {
  skip_on_cran()  # starts two background R sessions
  W <- psp_weights(rook_weights(5, 6), style = "W")
  d <- sim_psprobit(W, T = 3, beta = c(0.2, 1, -0.5), rho = 0.4, seed = 21)
  fit <- suppressWarnings(psprobit(
    y ~ x1 + x2, d, c("id", "time"), W, impacts = FALSE,
    control = psp_control(iter = 120, warmup = 60, chains = 2, seed = 3)))
  run_all <- function() {
    list(impacts = psp_impacts(fit, draws = 40, by_period = TRUE)$draws,
         unit = psp_impacts(fit, draws = 20, type = "unit")$draws,
         predict = predict(fit, draws = 37),
         ppc = psp_ppc(fit, draws = 30, seed = 4)$share$replicated,
         response = psp_response_path(fit, "x1", draws = 25, nsim = 3, seed = 5)$draws,
         loglik = suppressWarnings(psp_loglik(fit, draws = 19, reps = 20, type = "conditional")))
  }
  old <- future::plan("sequential")
  on.exit(future::plan(old), add = TRUE)
  seq_res <- run_all()
  future::plan("multisession", workers = 2)
  par_res <- run_all()
  for (nm in names(seq_res)) expect_identical(par_res[[nm]], seq_res[[nm]], info = nm)
})

test_that("chains give the same draws on one thread or several", {
  W <- psp_weights(rook_weights(4, 4), style = "W")
  d <- sim_psprobit(W, T = 3, beta = c(0, 1), rho = 0.4, seed = 22)
  fit <- function(cores) {
    suppressWarnings(psprobit(y ~ x1, d, c("id", "time"), W, impacts = FALSE,
                              control = psp_control(iter = 50, warmup = 50, chains = 2,
                                                    cores = cores, seed = 6)))$draws
  }
  expect_identical(fit(2), fit(1))
})
