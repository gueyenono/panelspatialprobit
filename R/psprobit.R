#' Bayesian panel spatial probit and ordered probit models
#'
#' Fits spatial probit (`psprobit()`) and spatial ordered probit
#' (`psoprobit()`) models to balanced panel data by Markov chain Monte Carlo.
#' The latent outcome of unit \eqn{i} in period \eqn{t} follows
#' \deqn{z_t = \rho W z_t + X_t\beta + W X_t\theta + \alpha + u_t,
#'   \quad u_t = \lambda M u_t + \varepsilon_t, \quad
#'   \varepsilon_t \sim N(0, I_N),}
#' with \eqn{y_{it} = 1\{z_{it} > 0\}} for binary outcomes and
#' \eqn{y_{it} = j} if \eqn{\kappa_{j-1} < z_{it} \le \kappa_j} for ordered
#' outcomes (\eqn{\kappa_1 = 0}).
#'
#' @param formula Model formula with an intercept. For `psprobit()` the
#'   response is 0/1, logical or a two-level factor; for `psoprobit()` an
#'   (ordered) factor or integer codes with at least three categories.
#' @param data A data frame with one row per unit and period (balanced panel).
#'   Outcomes may be missing (`NA`): their latent values are imputed during
#'   sampling, which keeps the spatial structure intact. Regressors may not.
#' @param index Character vector naming the unit and time columns of `data`.
#' @param W Spatial weights for the units (see [psp_weights()]); its row names
#'   identify the units. It is used as given, so row-standardise it beforehand
#'   (e.g. `psp_weights(W, style = "W")`) if desired. For weights that change
#'   over time, pass a list with one matrix per period, named by period (or in
#'   period order); \eqn{W_t} then replaces \eqn{W} in period \eqn{t},
#'   including in the Durbin terms and the spatially lagged outcome of
#'   [psp_dynamic()].
#' @param model Spatial specification:
#'   * `"sar"`: spatial lag of the latent outcome;
#'   * `"sdm"`: spatial Durbin, SAR plus \eqn{WX};
#'   * `"slx"`: \eqn{WX} only;
#'   * `"sem"`: spatial error, \eqn{(I - \lambda M)(z_t - X_t\beta - \alpha)
#'     = \varepsilon_t};
#'   * `"sdem"`: spatial Durbin error, SEM plus \eqn{WX};
#'   * `"sarar"`: spatial lag and spatial error;
#'   * `"none"`.
#' @param M Spatial weights for the error term of `"sem"`, `"sdem"` and
#'   `"sarar"` (default: `W`). Its row names must match those of `W`. Like
#'   `W`, it may be a list with one matrix per period.
#' @param effects Unit heterogeneity: `"random"` (\eqn{\alpha_i \sim N(0,
#'   \sigma_\alpha^2)}), `"pooled"` (\eqn{\alpha = 0}) or `"spatial_random"`
#'   (spatially autocorrelated effects \eqn{\alpha = \phi W\alpha + a},
#'   \eqn{a_i \sim N(0, \sigma_\alpha^2)}, as in Wang and Kockelman 2009;
#'   \eqn{W} must then be time-invariant and \eqn{\phi} uses `rho_range`).
#' @param cre Correlated random effects: a one-sided formula (or `TRUE` for all
#'   regressors) giving the regressors whose unit means are added to the model
#'   (Mundlak device, columns `mean_<x>`), or
#'   `list(type = "chamberlain", vars = ~ x1)` to add the values of the
#'   regressors in every period instead (Chamberlain device, columns
#'   `cre_<x>_<period>`). Chamberlain needs \eqn{T} extra columns per regressor.
#' @param durbin For `"sdm"`/`"slx"`: a one-sided formula selecting the
#'   regressors whose spatial lags enter the model (default: all).
#' @param time_effects Logical; add period dummies.
#' @param dynamic Optional dynamic specification from [psp_dynamic()]
#'   (lagged outcomes and initial conditions). The first period is then used
#'   only as the initial condition, and `predict()` returns values for the
#'   other periods.
#' @param errors Distribution of the latent errors: `"normal"` (probit) or
#'   Student-t from [psp_t()].
#' @param priors Priors, see [psp_priors()].
#' @param control MCMC settings, see [psp_control()].
#' @param impacts Logical; compute impacts (average direct, indirect and total
#'   marginal effects on probabilities) after sampling.
#' @section Convergence:
#' The sampler augments the data with the latent outcomes (Albert and Chib
#' 1993), which mixes well for balanced outcomes but slowly in some
#' situations:
#' * rare outcomes (roughly under 10% in one category);
#' * near-separation (outcomes that barely vary within units or periods);
#' * strong latent dynamics (`psp_dynamic(latent = TRUE)` with
#'   \eqn{|\gamma|} near or above 1);
#' * `"sarar"` with `M = W`.
#'
#' The intercept, \eqn{\rho} and \eqn{\sigma_\alpha} are then strongly
#' correlated a posteriori and their chains are autocorrelated. This slowness
#' is inherent to data augmentation with imbalanced outcomes (Johndrow et al.
#' 2019), and reparameterisations do not remove it. Run longer chains
#' (`psp_control(iter = , warmup = )`, with `thin` to save memory) until the
#' warning about R-hat or bulk ESS disappears.
#' @return An object of class `psprobit` (and `psoprobit` for ordered models).
#' @references
#' Albert, J. H. and Chib, S. (1993). Bayesian analysis of binary and
#' polychotomous response data. *JASA* 88, 669-679.
#'
#' Johndrow, J. E., Smith, A., Pillai, N. and Dunson, D. B. (2019). MCMC for
#' imbalanced categorical data. *JASA* 114, 1394-1403.
#'
#' Kakamu, K. and Wago, H. (2005). Bayesian spatial panel probit model with an
#' application to business cycle in Japan. *MODSIM 2005*, 856-862.
#'
#' LeSage, J. P. and Pace, R. K. (2009). *Introduction to Spatial
#' Econometrics*. CRC Press.
#'
#' Wang, X. and Kockelman, K. M. (2009). Bayesian inference for ordered
#' response data with a dynamic spatial-ordered probit model. *Journal of
#' Regional Science* 49, 877-913.
#' @examples
#' W <- psp_weights(rook_weights(5, 5), style = "W")
#' d <- sim_psprobit(W, T = 4, beta = c(0, 1, -1), rho = 0.5, sigma_alpha = 0.5,
#'                   seed = 1)
#' fit <- psprobit(y ~ x1 + x2, data = d, index = c("id", "time"), W = W,
#'                 control = psp_control(iter = 300, warmup = 200, chains = 2, seed = 1))
#' summary(fit)
#' @export
psprobit <- function(formula, data, index, W, model = c("sar", "sdm", "slx", "sem", "sdem", "sarar", "none"),
                     effects = c("random", "pooled", "spatial_random"), cre = NULL,
                     durbin = NULL,
                     time_effects = FALSE, M = NULL, dynamic = NULL, errors = "normal",
                     priors = psp_priors(), control = psp_control(), impacts = TRUE) {
  model <- .psp_match_model(model)
  effects <- match.arg(effects)
  design <- .psp_design(formula, data, index, W, model, cre, durbin, time_effects, "binary", M,
                        .psp_check_dynamic(dynamic))
  .psp_fit(design, "binary", model, effects, priors, control, impacts, match.call(),
           .psp_check_errors(errors))
}

#' @rdname psprobit
#' @export
psoprobit <- function(formula, data, index, W, model = c("sar", "sdm", "slx", "sem", "sdem", "sarar", "none"),
                      effects = c("random", "pooled", "spatial_random"), cre = NULL,
                      durbin = NULL,
                      time_effects = FALSE, M = NULL, dynamic = NULL, errors = "normal",
                      priors = psp_priors(), control = psp_control(), impacts = TRUE) {
  model <- .psp_match_model(model)
  effects <- match.arg(effects)
  design <- .psp_design(formula, data, index, W, model, cre, durbin, time_effects, "ordered", M,
                        .psp_check_dynamic(dynamic))
  .psp_fit(design, "ordered", model, effects, priors, control, impacts, match.call(),
           .psp_check_errors(errors))
}

.psp_check_dynamic <- function(dynamic) {
  if (is.null(dynamic) || isFALSE(dynamic)) return(NULL)
  if (isTRUE(dynamic)) return(psp_dynamic())
  if (!inherits(dynamic, "psp_dynamic")) {
    stop("'dynamic' must come from psp_dynamic().", call. = FALSE)
  }
  dynamic
}

.psp_match_model <- function(model) {
  match.arg(model, c("sar", "sdm", "slx", "sem", "sdem", "sarar", "none"))
}

.psp_fit <- function(design, family, model, effects, priors, control, impacts, call,
                     errors = .psp_check_errors("normal")) {
  if (!inherits(priors, "psp_priors")) stop("'priors' must come from psp_priors().", call. = FALSE)
  if (!inherits(control, "psp_control")) stop("'control' must come from psp_control().",
                                              call. = FALSE)
  t0 <- proc.time()
  has_rho <- model %in% c("sar", "sdm", "sarar")
  has_lambda <- model %in% c("sem", "sdem", "sarar")
  has_re <- effects != "pooled"
  spatial_re <- effects == "spatial_random"
  if (spatial_re && is.list(design$W)) {
    stop("effects = \"spatial_random\" needs time-invariant weights 'W'.", call. = FALSE)
  }
  ordered <- family == "ordered"
  X <- design$X
  k <- ncol(X)
  J <- design$J
  rr <- control$rho_range

  no_grid <- list(grid = c(-0.5, 0.5), lndet = c(0, 0), method = "none")
  ld <- if (has_rho || spatial_re) .psp_lndet_avg(design$W, rr[1], rr[2], control$lndet)
        else no_grid
  lr <- control$lambda_range
  ld_M <- no_grid
  if (has_lambda) {
    ld_M <- if (has_rho && identical(design$M, design$W) && identical(lr, rr)) ld
            else .psp_lndet_avg(design$M, lr[1], lr[2], control$lndet)
  }

  bm <- rep_len(priors$beta$mean, k)
  bs <- rep_len(priors$beta$sd, k)
  sp <- priors$sigma_alpha
  sigma_prior <- if (sp$family == "half_cauchy") 0L else 1L
  sigma_p <- if (sigma_prior == 0L) c(sp$scale, 0) else c(sp$shape, sp$scale)
  kp <- priors$cutpoints
  kappa_start <- numeric(0)
  # missing outcomes are coded 0 (ordered) or -1 (binary) for the sampler
  y_c <- design$y
  y_c[is.na(y_c)] <- if (ordered) 0L else -1L
  if (ordered) {
    cp <- stats::qnorm(cumsum(tabulate(design$y, J))[-J] / sum(!is.na(design$y)))
    kappa_start <- (cp - cp[1])[-1]
  }

  # keep about impacts_draws random-effect draws in total (unit-specific impacts)
  nsave <- control$iter %/% control$thin
  alpha_every <- if (has_re) max(1L, (nsave * control$chains) %/% control$impacts_draws) else 0L
  seed <- control$seed
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  streams <- .psp_streams(seed, control$chains)
  spec <- list(
    beta_mean = as.double(bm), beta_prec = as.double(diag(1 / bs^2, k)),
    rho = list(a = priors$rho$a, b = priors$rho$b, min = rr[1], max = rr[2],
               grid = as.double(ld$grid), lndet = as.double(ld$lndet)),
    lambda = list(a = priors$lambda$a, b = priors$lambda$b, min = lr[1], max = lr[2],
                  grid = as.double(ld_M$grid), lndet = as.double(ld_M$lndet)),
    phi_a = priors$phi$a, phi_b = priors$phi$b,
    nu = if (errors$type == "t") as.double(errors$df) else 0,
    sigma = list(family = sigma_prior, p1 = as.double(sigma_p[1]), p2 = as.double(sigma_p[2])),
    kappa_flat = kp$family == "flat",
    kappa_mean = if (kp$family == "flat") 0 else as.double(kp$mean),
    kappa_sd = if (kp$family == "flat") 1 else as.double(kp$sd),
    iter = as.integer(control$iter), warmup = as.integer(control$warmup),
    thin = as.integer(control$thin), z_sweeps = as.integer(control$z_sweeps),
    cut_sampler = if (control$cutpoint_sampler == "transform") 0L else 1L,
    scale_move = isTRUE(control$scale_move),
    beta_start = numeric(k), rho_start = 0, lambda_start = 0, sigma2_start = 1,
    kappa_start = as.double(kappa_start), alpha_every = as.integer(alpha_every))
  model_rs <- .psp_rs_model(design, y_c, ordered, has_rho, has_lambda, has_re, spatial_re)
  seeds <- vapply(streams, function(s) s[-1], integer(6))
  threads <- .psp_threads(control$cores, control$threads_per_chain)
  spec$threads_per_chain <- threads[["per_chain"]]
  out <- psp_rs_sample(model_rs, spec, seeds, threads[["cores"]], isTRUE(control$progress))

  pnames <- c(colnames(X), if (has_rho) "rho", if (has_lambda) "lambda",
              if (has_re) "sigma_alpha", if (spatial_re) "phi",
              if (ordered && J > 2) sprintf("cut[%s|%s]", design$levels[2:(J - 1)],
                                            design$levels[3:J]))
  arr <- simplify2array(lapply(out, function(o) o$draws))
  arr <- array(arr, dim = c(nrow(out[[1]]$draws), length(pnames), length(out)))
  arr <- aperm(arr, c(1, 3, 2))
  dimnames(arr) <- list(NULL, NULL, pnames)
  draws <- posterior::as_draws_array(arr)

  fit <- list(
    call = call, family = family, model = model, effects = effects, draws = draws,
    coefficients = colMeans(posterior::as_draws_matrix(draws)),
    y = design$y, levels = design$levels, J = J, X = X, roles = design$roles,
    durbin_of = design$durbin_of, W = design$W, M = if (has_lambda) design$M,
    N = design$N, T = design$T,
    unit_ids = design$unit_ids, time_levels = design$time_levels,
    row_order = design$row_order, index = design$index, terms = design$terms,
    dynamic = design$dynamic, spec = design$spec, errors = errors,
    lndet = ld, lndet_M = if (has_lambda) ld_M, priors = priors, control = control, seed = seed,
    z_mean = Reduce(`+`, lapply(out, `[[`, "z_mean")) / length(out),
    alpha_mean = if (has_re) Reduce(`+`, lapply(out, `[[`, "alpha_mean")) / length(out),
    alpha_draws = if (has_re) do.call(rbind, lapply(out, `[[`, "alpha_draws")),
    alpha_index = if (has_re) {
      unlist(lapply(seq_along(out), function(c) {
        (c - 1L) * nsave + seq(1L, by = alpha_every, length.out = nrow(out[[c]]$alpha_draws))
      }))
    },
    cut_accept = if (ordered) do.call(rbind, lapply(out, `[[`, "cut_accept")),
    v_mean = if (errors$type == "t") {
      v <- Reduce(`+`, lapply(out, `[[`, "v_mean")) / length(out)
      v_data <- v
      v_data[design$row_order] <- v  # data row order
      v_data
    }
  )
  class(fit) <- c(if (ordered) "psoprobit", "psprobit")
  fit$diagnostics <- .psp_diagnostics(fit)
  if (impacts && any(fit$roles == "x")) fit$impacts <- psp_impacts(fit)
  fit$time <- proc.time() - t0
  fit
}

# Independent L'Ecuyer-CMRG streams so results do not depend on 'cores'.
.psp_streams <- function(seed, n) {
  .psp_with_rng_restore({
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    s <- get(".Random.seed", envir = globalenv())
    streams <- vector("list", n)
    for (i in seq_len(n)) {
      streams[[i]] <- s
      s <- parallel::nextRNGStream(s)
    }
    streams
  })
}

.psp_with_rng_restore <- function(expr) {
  kind <- RNGkind()
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv())
  on.exit({
    RNGkind(kind[1], kind[2], kind[3])
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  expr
}

.psp_diagnostics <- function(fit) {
  s <- posterior::summarise_draws(fit$draws, "rhat", "ess_bulk", "ess_tail")
  nchain <- posterior::nchains(fit$draws)
  bad_rhat <- nchain > 1 & !is.na(s$rhat) & s$rhat > 1.01
  bad_ess <- !is.na(s$ess_bulk) & s$ess_bulk < 100 * nchain
  if (any(bad_rhat) || any(bad_ess)) {
    warning(sprintf(paste("Possible convergence problems (R-hat > 1.01 or bulk ESS < %d)",
                          "for: %s. Consider more iterations."), 100L * nchain,
                    paste(s$variable[bad_rhat | bad_ess], collapse = ", ")), call. = FALSE)
  }
  s
}
