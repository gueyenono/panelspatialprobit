#' Prior distributions
#'
#' Constructors for the prior distributions used by [psp_priors()].
#'
#' @param mean,sd Mean and standard deviation of a normal prior.
#' @param a,b Shape parameters of a beta prior, rescaled to the admissible
#'   interval of the spatial parameter.
#' @param scale Scale of a half-Cauchy prior on a standard deviation.
#' @param shape Shape of an inverse-gamma prior on a variance; `scale` is then
#'   its scale.
#' @return An object of class `psp_prior`.
#' @name priors
#' @examples
#' prior_normal(0, 5)
#' prior_half_cauchy(1)
NULL

#' @rdname priors
#' @export
prior_normal <- function(mean = 0, sd = 10) {
  stopifnot(is.numeric(mean), is.numeric(sd), all(sd > 0))
  structure(list(family = "normal", mean = mean, sd = sd), class = "psp_prior")
}

#' @rdname priors
#' @export
prior_beta <- function(a = 1.01, b = 1.01) {
  stopifnot(a > 0, b > 0)
  structure(list(family = "beta", a = a, b = b), class = "psp_prior")
}

#' @rdname priors
#' @export
prior_half_cauchy <- function(scale = 1) {
  stopifnot(scale > 0)
  structure(list(family = "half_cauchy", scale = scale), class = "psp_prior")
}

#' @rdname priors
#' @export
prior_inv_gamma <- function(shape = 1, scale = 1) {
  stopifnot(shape > 0, scale > 0)
  structure(list(family = "inv_gamma", shape = shape, scale = scale), class = "psp_prior")
}

#' @rdname priors
#' @export
prior_flat <- function() structure(list(family = "flat"), class = "psp_prior")

#' @export
print.psp_prior <- function(x, ...) {
  pars <- x[setdiff(names(x), "family")]
  cat(sprintf("%s(%s)\n", x$family,
              paste(names(pars), vapply(pars, function(p) paste(format(p), collapse = ","),
                                        character(1)), sep = " = ", collapse = ", ")))
  invisible(x)
}

#' Priors for panel spatial probit models
#'
#' @param beta Normal prior for all regression coefficients (including Durbin,
#'   correlated-random-effects and time-effect columns). `mean` and `sd` may be
#'   scalars or vectors with one entry per column of the design matrix.
#' @param rho Beta prior for the spatial autoregressive parameter, rescaled to
#'   the interval given by `rho_range` in [psp_control()].
#' @param lambda Beta prior for the spatial error parameter (models `"sem"`,
#'   `"sdem"` and `"sarar"`), rescaled to `lambda_range` in [psp_control()].
#' @param phi Beta prior for the spatial parameter of spatial random effects
#'   (`effects = "spatial_random"`), rescaled to `rho_range`.
#' @param sigma_alpha Prior for the standard deviation of the random effects:
#'   [prior_half_cauchy()] (default) or [prior_inv_gamma()] on the variance.
#' @param cutpoints Prior for the free cut-points of ordered models:
#'   [prior_normal()] (truncated to the ordering constraint) or [prior_flat()].
#' @return An object of class `psp_priors`.
#' @examples
#' psp_priors(beta = prior_normal(0, 5))
#' @export
psp_priors <- function(beta = prior_normal(0, 10), rho = prior_beta(1.01, 1.01),
                       lambda = prior_beta(1.01, 1.01), sigma_alpha = prior_half_cauchy(1),
                       cutpoints = prior_normal(0, 10), phi = prior_beta(1.01, 1.01)) {
  if (beta$family != "normal") stop("'beta' must be a normal prior.", call. = FALSE)
  if (rho$family != "beta") stop("'rho' must be a beta prior.", call. = FALSE)
  if (lambda$family != "beta") stop("'lambda' must be a beta prior.", call. = FALSE)
  if (phi$family != "beta") stop("'phi' must be a beta prior.", call. = FALSE)
  if (!sigma_alpha$family %in% c("half_cauchy", "inv_gamma")) {
    stop("'sigma_alpha' must be a half-Cauchy or inverse-gamma prior.", call. = FALSE)
  }
  if (!cutpoints$family %in% c("normal", "flat")) {
    stop("'cutpoints' must be a normal or flat prior.", call. = FALSE)
  }
  structure(list(beta = beta, rho = rho, lambda = lambda, sigma_alpha = sigma_alpha,
                 cutpoints = cutpoints, phi = phi),
            class = "psp_priors")
}

#' MCMC and computation settings
#'
#' @param iter Number of post-warm-up iterations per chain.
#' @param warmup Number of warm-up iterations per chain (discarded).
#' @param thin Thinning interval for stored draws.
#' @param chains Number of chains.
#' @param cores Number of chains run at once, on separate threads (on every
#'   operating system, Windows included). Results do not depend on `cores`.
#'   Under `R CMD check` (`_R_CHECK_LIMIT_CORES_`) at most 2 are used.
#' @param seed Optional integer seed. Given `seed`, results are reproducible
#'   on a given platform, whatever `cores` and the user's [RNGkind()]: each
#'   chain draws from its own L'Ecuyer-CMRG stream derived from `seed`, and
#'   the session's random number state is left unchanged. Without `seed`, one
#'   is drawn from the session's generator (`fit$seed`).
#' @param z_sweeps Number of Gibbs sweeps over the latent variables per
#'   iteration.
#' @param rho_range Admissible interval for the spatial autoregressive
#'   parameter.
#' @param lambda_range Admissible interval for the spatial error parameter.
#' @param lndet Method for log-determinants: `"auto"`, `"eigen"` or `"sparse"`.
#' @param cutpoint_sampler `"transform"` (Albert-Chib Gibbs step plus a joint
#'   transformation move of cut-points and latent values) or `"albert_chib"`.
#' @param scale_move Logical; add a generalised Gibbs move that rescales the
#'   latent values, coefficients, random effects and cut-points jointly
#'   (Liu and Sabatti, 2000). It leaves the posterior unchanged and greatly
#'   improves the mixing of the coefficients.
#' @param impacts_draws Number of posterior draws used to compute impacts and
#'   predictions.
#' @param impacts_exact_max Largest number of units for which diagonals of
#'   matrix inverses are computed exactly; above it, a stochastic estimator is
#'   used.
#' @param impacts_probes Number of probe vectors of the stochastic estimator.
#' @param mixture_draws Number of simulated sets of mixing variances over
#'   which probabilities and impacts are averaged for Student-t errors
#'   ([psp_t()]).
#' @param progress Logical; print sampling progress, for each chain.
#' @section Parallel post-processing:
#' [psp_impacts()], [predict()][predict.psprobit], [psp_ppc()],
#' [psp_response_path()], [psp_loglik()], `logLik()` and [psp_dic()] work
#' through the posterior draws in blocks and follow the [future::plan()] set
#' by the user: for example `future::plan("multisession", workers = 4)`, or
#' `future::plan("multicore", workers = 4)` on Unix-alikes. The default
#' sequential plan runs in the current R session. Results do not depend on
#' the plan or the number of workers.
#' @return An object of class `psp_control`.
#' @examples
#' psp_control(iter = 500, warmup = 250, chains = 2)
#' @export
psp_control <- function(iter = 2000, warmup = 1000, thin = 1, chains = 4, cores = 1,
                        seed = NULL, z_sweeps = 1, rho_range = c(-1, 1),
                        lambda_range = c(-1, 1),
                        lndet = c("auto", "eigen", "sparse"),
                        cutpoint_sampler = c("transform", "albert_chib"),
                        scale_move = TRUE, impacts_draws = 500, impacts_exact_max = 500,
                        impacts_probes = 50, mixture_draws = 50, progress = interactive()) {
  lndet <- match.arg(lndet)
  cutpoint_sampler <- match.arg(cutpoint_sampler)
  stopifnot(iter >= 1, warmup >= 0, thin >= 1, iter >= thin, chains >= 1, cores >= 1,
            z_sweeps >= 1, length(rho_range) == 2, rho_range[1] < rho_range[2],
            length(lambda_range) == 2, lambda_range[1] < lambda_range[2],
            impacts_draws >= 1)
  structure(list(iter = as.integer(iter), warmup = as.integer(warmup),
                 thin = as.integer(thin), chains = as.integer(chains),
                 cores = as.integer(cores), seed = seed, z_sweeps = as.integer(z_sweeps),
                 rho_range = rho_range, lambda_range = lambda_range, lndet = lndet, cutpoint_sampler = cutpoint_sampler,
                 scale_move = isTRUE(scale_move),
                 impacts_draws = as.integer(impacts_draws),
                 impacts_exact_max = impacts_exact_max,
                 impacts_probes = as.integer(impacts_probes),
                 mixture_draws = as.integer(mixture_draws),
                 progress = isTRUE(progress)),
            class = "psp_control")
}
