test_that("binary fit returns a well-formed object and methods work", {
  s <- small_data()
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W,
                                   control = fast_control()))
  expect_s3_class(fit, "psprobit")
  expect_named(coef(fit), c("(Intercept)", "x1", "x2", "rho", "sigma_alpha"))
  expect_equal(dim(vcov(fit)), c(5, 5))
  expect_equal(dim(confint(fit)), c(5, 2))
  expect_equal(nobs(fit), 48)
  expect_equal(posterior::ndraws(posterior::as_draws(fit)), 400)
  expect_output(print(summary(fit)), "Posterior summary")
  tab <- panelspatialprobit:::.psp_round_table(summary(fit)$coefficients, 3)
  expect_equal(tab$ess_bulk, round(tab$ess_bulk))
  expect_equal(tab$mean, round(summary(fit)$coefficients$mean, 3))
  expect_output(print(fit), "Bayesian panel spatial probit")
  expect_s3_class(tidy(fit), "data.frame")
  expect_equal(nrow(glance(fit)), 1)
  p <- predict(fit, draws = 20)
  expect_length(p, 48)
  expect_true(all(p > 0 & p < 1))
  expect_true(all(predict(fit, type = "class", draws = 20) %in% c("0", "1")))
  expect_s3_class(fit$impacts, "psp_impacts")
  expect_output(print(fit$impacts), "direct")
})

test_that("predictions follow the row order of the data", {
  s <- small_data()
  perm <- sample(nrow(s$d))
  f1 <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, impacts = FALSE,
                                  control = fast_control()))
  f2 <- suppressWarnings(psprobit(y ~ x1 + x2, s$d[perm, ], c("id", "time"), s$W,
                                  impacts = FALSE, control = fast_control()))
  expect_equal(predict(f2, draws = 10), predict(f1, draws = 10)[perm])
})

test_that("results are reproducible and independent of cores", {
  skip_on_os("windows")
  s <- small_data()
  fit_fun <- function(cores) {
    suppressWarnings(psprobit(y ~ x1, s$d, c("id", "time"), s$W, impacts = FALSE,
                              control = psp_control(iter = 50, warmup = 20, chains = 2,
                                                    cores = cores, seed = 11)))
  }
  a <- fit_fun(1)
  b <- fit_fun(1)
  c <- fit_fun(2)
  expect_identical(unclass(a$draws), unclass(b$draws))
  expect_identical(unclass(a$draws), unclass(c$draws))
})

test_that("fitting does not disturb the user's RNG state beyond one draw", {
  s <- small_data()
  set.seed(5)
  invisible(suppressWarnings(psprobit(y ~ x1, s$d, c("id", "time"), s$W, impacts = FALSE,
                                      control = psp_control(iter = 20, warmup = 10, chains = 1,
                                                            seed = 1))))
  after_fit <- runif(1)
  set.seed(5)
  expect_identical(runif(1), after_fit)
})

test_that("ordered fit, pooled effects and model variants run", {
  s <- small_data(cutpoints = c(0.7, 1.5))
  fit <- suppressWarnings(psoprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, effects = "pooled",
                                    control = fast_control()))
  expect_s3_class(fit, "psoprobit")
  expect_true(all(c("cut[3|4]", "cut[2|3]") %in% names(coef(fit))))
  expect_false("sigma_alpha" %in% names(coef(fit)))
  pr <- predict(fit, draws = 10)
  expect_equal(dim(pr), c(48, 4))
  expect_equal(unname(rowSums(pr)), rep(1, 48), tolerance = 1e-10)
  for (mod in c("sdm", "slx", "none")) {
    f <- suppressWarnings(psprobit(y ~ x1 + x2, small_data()$d, c("id", "time"), s$W,
                                   model = mod, cre = ~ x1, time_effects = TRUE,
                                   impacts = FALSE, control = fast_control()))
    expect_equal("rho" %in% names(coef(f)), mod == "sdm")
  }
  fa <- suppressWarnings(psoprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, impacts = FALSE,
                                   priors = psp_priors(cutpoints = prior_flat(),
                                                       sigma_alpha = prior_inv_gamma(2, 1)),
                                   control = fast_control(cutpoint_sampler = "albert_chib")))
  expect_s3_class(fa, "psoprobit")
})

test_that("summary reports classification and polr-style thresholds", {
  s <- small_data(cutpoints = c(0.7, 1.5))
  fit <- suppressWarnings(psoprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, impacts = FALSE,
                                    control = fast_control()))
  sm <- summary(fit)
  expect_equal(rownames(sm$thresholds), c("1|2", "2|3", "3|4"))
  expect_equal(sm$thresholds["1|2", "mean"], -sm$coefficients["(Intercept)", "mean"])
  expect_equal(sm$thresholds["2|3", "mean"],
               sm$coefficients["cut[2|3]", "mean"] - sm$coefficients["(Intercept)", "mean"])
  expect_output(print(sm), "polr")
  cl <- psp_classification(fit, draws = 10)
  expect_true(cl$hit_rate >= 0 && cl$hit_rate <= 1)
  expect_true(is.na(cl$auc))
  expect_equal(sum(cl$confusion), 48)

  b <- small_data()
  fb <- suppressWarnings(psprobit(y ~ x1 + x2, b$d, c("id", "time"), b$W, impacts = FALSE,
                                  control = fast_control()))
  cb <- psp_classification(fb, draws = 10)
  expect_true(cb$auc > 0.5 && cb$auc <= 1)
  # AUC agrees with a direct pairwise computation
  p <- predict(fb, draws = 10)
  y <- b$d$y
  pairs <- outer(p[y == 1], p[y == 0], ">") + 0.5 * outer(p[y == 1], p[y == 0], "==")
  expect_equal(cb$auc, mean(pairs))
  expect_output(print(summary(fb)), "AUC")
})

test_that("Chamberlain CRE fits", {
  s <- small_data()
  f <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, impacts = FALSE,
                                 cre = list(type = "chamberlain", vars = ~ x1),
                                 control = fast_control()))
  expect_true(all(paste0("cre_x1_", 1:3) %in% names(coef(f))))
})

test_that("missing outcomes are imputed", {
  s <- small_data()
  d <- s$d
  miss <- c(3, 10, 20, 41)
  d$y[miss] <- NA
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, d, c("id", "time"), s$W, impacts = FALSE,
                                   control = fast_control()))
  expect_equal(nobs(fit), 44)
  expect_false(anyNA(posterior::as_draws_matrix(fit$draws)))
  p <- predict(fit, draws = 10)
  expect_length(p, 48)
  expect_false(anyNA(p))
  expect_equal(sum(psp_classification(fit, draws = 10)$confusion), 44)
  expect_output(print(summary(fit, classification = FALSE)), "<NA>")
  s2 <- small_data(cutpoints = c(0.7, 1.5))
  d2 <- s2$d
  d2$y[miss] <- NA
  fo <- suppressWarnings(psoprobit(y ~ x1 + x2, d2, c("id", "time"), s2$W, impacts = FALSE,
                                   control = fast_control()))
  expect_equal(nobs(fo), 44)
  expect_error(psprobit(y ~ x1, d, c("id", "time"), s$W, dynamic = psp_dynamic()),
               "dynamic")
  d$x1[1] <- NA
  expect_error(psprobit(y ~ x1, d, c("id", "time"), s$W), "regressors")
})
