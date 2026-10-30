#' Methods for fitted panel spatial probit models
#'
#' @param x,object A fitted `psprobit` object.
#' @param digits Number of digits to print.
#' @param parm Parameters for which intervals are returned (default: all).
#' @param level Credible level.
#' @param type For `predict()`: `"prob"` (posterior mean probabilities; a
#'   matrix with one column per category for ordered models), `"link"`
#'   (posterior mean of the systematic part of the latent outcome) or
#'   `"class"` (most probable outcome).
#' @param draws Number of posterior draws used by `predict()`.
#' @param newdata Optional data frame for counterfactual predictions: the
#'   same units and periods as the fitted data, with new regressor values (the
#'   response column may be left out). Correlated-random-effects columns are
#'   kept at their fitted values, since they stand for fixed unit
#'   characteristics. Not available for models with lagged outcomes.
#' @param pars Parameters to plot (default: all).
#' @param ... Unused.
#' @name psprobit-methods
#' @return `coef()` a named vector, `vcov()` a matrix, `confint()` a matrix of
#'   interval bounds, `summary()` an object of class `summary.psprobit`,
#'   `predict()`/`fitted()` vectors or matrices in the row order of the data.
NULL

#' @rdname psprobit-methods
#' @export
print.psprobit <- function(x, digits = 3, ...) {
  cat(.psp_title(x), "\n")
  cat("Call: ", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  print(round(x$coefficients, digits))
  invisible(x)
}

.psp_title <- function(x) {
  sprintf("Bayesian panel spatial %sprobit (model = %s, effects = %s%s)",
          if (inherits(x, "psoprobit")) "ordered " else "", x$model, x$effects,
          if (!is.null(x$dynamic)) ", dynamic" else "")
}

#' @rdname psprobit-methods
#' @param classification Logical; report in-sample classification statistics
#'   (see [psp_classification()]) in `summary()`.
#' @export
summary.psprobit <- function(object, classification = TRUE, ...) {
  tab <- .psp_draws_table(object$draws)
  thresholds <- NULL
  if (inherits(object, "psoprobit")) {
    # MASS::polr parameterisation: no intercept, zeta_j = kappa_j - intercept
    m <- .psp_par_matrix(object)
    kap <- cbind(0, m[, grepl("^cut\\[", colnames(m)), drop = FALSE])
    zeta <- kap - m[, "(Intercept)"]
    colnames(zeta) <- sprintf("%s|%s", object$levels[-object$J], object$levels[-1])
    nch <- posterior::nchains(object$draws)
    arr <- array(zeta, c(nrow(zeta) / nch, nch, ncol(zeta)),
                 dimnames = list(NULL, NULL, colnames(zeta)))
    thresholds <- .psp_draws_table(posterior::as_draws_array(arr))
  }
  out <- list(title = .psp_title(object), call = object$call, coefficients = tab,
              thresholds = thresholds,
              classification = if (classification) psp_classification(object, draws = 100),
              N = object$N, T = object$T, chains = posterior::nchains(object$draws),
              draws = posterior::ndraws(object$draws),
              outcome = table(factor(object$y, labels = object$levels,
                                     levels = if (inherits(object, "psoprobit"))
                                       seq_len(object$J) else 0:1), useNA = "ifany"),
              impacts = object$impacts, time = object$time)
  class(out) <- "summary.psprobit"
  out
}

#' @export
print.summary.psprobit <- function(x, digits = 3, ...) {
  cat(x$title, "\n")
  cat("Call: ", paste(deparse(x$call), collapse = "\n"), "\n", sep = "")
  cat(sprintf("Units: %d, periods: %d, observations: %d\n", x$N, x$T, x$N * x$T))
  cat(sprintf("Chains: %d, posterior draws: %d\n", x$chains, x$draws))
  cat("Outcome distribution:\n")
  print(x$outcome)
  cat("\nPosterior summary:\n")
  print(.psp_round_table(x$coefficients, digits))
  if (!is.null(x$thresholds)) {
    cat("\nThresholds in the MASS::polr parameterisation (no intercept,",
        "zeta_j = cut_j - (Intercept)):\n")
    print(.psp_round_table(x$thresholds, digits))
  }
  if (!is.null(x$classification)) {
    cl <- x$classification
    cat(sprintf("\nIn-sample classification: hit rate %.3f%s\n", cl$hit_rate,
                if (!is.na(cl$auc)) sprintf(", AUC %.3f", cl$auc) else ""))
  }
  if (!is.null(x$impacts)) {
    cat("\n")
    print(x$impacts, digits = digits)
  }
  invisible(x)
}

# Round a posterior summary table for printing: effective sample sizes as
# whole numbers, everything else to 'digits' decimals (a data frame, so each
# column is formatted separately)
.psp_round_table <- function(tab, digits) {
  tab <- as.data.frame(tab)
  ess <- grepl("^ess_", names(tab))
  tab[!ess] <- lapply(tab[!ess], round, digits)
  tab[ess] <- lapply(tab[ess], round)
  tab
}

# Posterior summary table (one row per parameter)
.psp_draws_table <- function(draws) {
  d <- posterior::summarise_draws(draws, "mean", "sd",
                                  ~ posterior::quantile2(.x, probs = c(0.025, 0.975)),
                                  "rhat", "ess_bulk", "ess_tail")
  tab <- as.data.frame(d[, -1])
  rownames(tab) <- d$variable
  tab
}

#' In-sample classification statistics
#'
#' Compares posterior mean outcome probabilities (`predict(type = "prob")`)
#' with the observed outcomes (missing outcomes are left out).
#'
#' @param fit A fitted `psprobit` object.
#' @param draws Number of posterior draws used for the probabilities.
#' @return A list with `hit_rate` (share of observations whose most probable
#'   outcome is the observed one), `auc` (binary models: area under the ROC
#'   curve; `NA` for ordered models), `brier` (mean squared error of the
#'   probability vector) and `confusion` (observed by predicted table).
#' @examples
#' W <- psp_weights(rook_weights(4, 4), style = "W")
#' d <- sim_psprobit(W, T = 3, beta = c(0, 1), rho = 0.4, seed = 2)
#' fit <- psprobit(y ~ x1, data = d, index = c("id", "time"), W = W, impacts = FALSE,
#'                 control = psp_control(iter = 200, warmup = 100, chains = 1, seed = 2))
#' psp_classification(fit, draws = 50)
#' @export
psp_classification <- function(fit, draws = NULL) {
  if (!inherits(fit, "psprobit")) stop("'fit' must be a psprobit object.", call. = FALSE)
  p <- predict.psprobit(fit, type = "prob", draws = draws)
  y <- fit$y
  y[fit$row_order] <- fit$y  # data row order, like p
  obs <- !is.na(y)
  y <- y[obs]
  p <- if (is.matrix(p)) p[obs, , drop = FALSE] else p[obs]
  if (inherits(fit, "psoprobit")) {
    P <- as.matrix(p)
    pred <- max.col(P, ties.method = "first")
    onehot <- outer(y, seq_len(fit$J), "==")
    auc <- NA_real_
    brier <- mean(rowSums((P - onehot)^2))
  } else {
    pred <- as.integer(p > 0.5)
    r <- rank(p)
    n1 <- sum(y == 1)
    n0 <- sum(y == 0)
    auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)  # Mann-Whitney
    brier <- mean((p - y)^2)
  }
  lev <- if (inherits(fit, "psoprobit")) seq_len(fit$J) else 0:1
  confusion <- table(observed = factor(y, lev, fit$levels),
                     predicted = factor(pred, lev, fit$levels))
  list(hit_rate = mean(pred == y), auc = auc, brier = brier, confusion = confusion)
}

#' @rdname psprobit-methods
#' @export
coef.psprobit <- function(object, ...) object$coefficients

#' @rdname psprobit-methods
#' @export
vcov.psprobit <- function(object, ...) stats::cov(posterior::as_draws_matrix(object$draws))

#' @rdname psprobit-methods
#' @export
confint.psprobit <- function(object, parm, level = 0.95, ...) {
  m <- posterior::as_draws_matrix(object$draws)
  if (missing(parm)) parm <- colnames(m)
  a <- (1 - level) / 2
  ci <- t(apply(m[, parm, drop = FALSE], 2, stats::quantile, probs = c(a, 1 - a)))
  ci
}

#' @rdname psprobit-methods
#' @export
nobs.psprobit <- function(object, ...) sum(!is.na(object$y))

#' @rdname psprobit-methods
#' @export
predict.psprobit <- function(object, newdata = NULL, type = c("prob", "link", "class"),
                             draws = NULL, ...) {
  type <- match.arg(type)
  if (!is.null(newdata)) object <- .psp_with_newdata(object, newdata)
  idx <- .psp_draw_subset(object, draws)
  pm <- .psp_par_matrix(object)[idx, , drop = FALSE]
  probes <- .psp_probes(object)
  mix <- .psp_t_mix(object)
  ordered <- inherits(object, "psoprobit")
  n <- object$N * object$T
  acc_m <- numeric(n)
  acc_p <- matrix(0, n, if (ordered) object$J else 1)
  for (d in seq_along(idx)) {
    rf <- .psp_reduced_form(object, pm[d, ], probes, mix = mix)
    s <- rf$s
    acc_m <- acc_m + rf$m
    if (ordered) {
      kap <- .psp_cutpoints(object, pm[d, ])
      for (j in seq_len(object$J)) {
        acc_p[, j] <- acc_p[, j] + .psp_mix_mean(function(s) {
          stats::pnorm((kap[j + 1] - rf$m) / s) - stats::pnorm((kap[j] - rf$m) / s)
        }, s)
      }
    } else {
      acc_p[, 1] <- acc_p[, 1] + .psp_mix_mean(function(s) stats::pnorm(rf$m / s), s)
    }
  }
  acc_m <- acc_m / length(idx)
  acc_p <- acc_p / length(idx)
  back <- function(v) {
    if (is.matrix(v)) {
      out <- v
      out[object$row_order, ] <- v
      return(out)
    }
    out <- v
    out[object$row_order] <- v
    out
  }
  if (type == "link") return(back(acc_m))
  if (type == "prob") {
    if (ordered) {
      colnames(acc_p) <- object$levels
      return(back(acc_p))
    }
    return(back(acc_p[, 1]))
  }
  cls <- if (ordered) object$levels[max.col(acc_p, ties.method = "first")]
         else object$levels[(acc_p[, 1] > 0.5) + 1]
  back(cls)
}

# Swap in the design of 'newdata' (same units and periods, new regressors)
.psp_with_newdata <- function(object, newdata) {
  if (!is.null(object$dynamic) && !.psp_latent_ar(object)) {
    stop("predict(newdata = ) is not available for models with lagged outcomes.",
         call. = FALSE)
  }
  spec <- object$spec
  resp <- all.vars(spec$formula[[2]])
  if (length(resp) == 1 && !resp %in% names(newdata)) {
    newdata[[resp]] <- if (inherits(object, "psoprobit")) {
      factor(rep_len(object$levels, nrow(newdata)), levels = object$levels)
    } else {
      rep_len(c(0L, 1L), nrow(newdata))
    }
  }
  des <- suppressWarnings(.psp_design(spec$formula, newdata, object$index, object$W,
                                      object$model, spec$cre, spec$durbin, spec$time_effects,
                                      object$family, object$M, object$dynamic))
  if (!identical(as.character(des$time_levels), as.character(object$time_levels))) {
    stop("'newdata' must contain the periods of the fitted data.", call. = FALSE)
  }
  if (!identical(colnames(des$X), colnames(object$X))) {
    stop("The design of 'newdata' differs from the fitted one (columns: ",
         paste(setdiff(union(colnames(des$X), colnames(object$X)),
                       intersect(colnames(des$X), colnames(object$X))), collapse = ", "),
         ").", call. = FALSE)
  }
  keep <- object$roles %in% c("cre", "initial")
  des$X[, keep] <- object$X[, keep]
  object$X <- des$X
  object$row_order <- des$row_order
  object
}

#' @rdname psprobit-methods
#' @export
fitted.psprobit <- function(object, ...) predict.psprobit(object, type = "prob", ...)

#' @rdname psprobit-methods
#' @export
plot.psprobit <- function(x, pars = NULL, ...) {
  arr <- posterior::as_draws_array(x$draws)
  vars <- if (is.null(pars)) posterior::variables(arr) else pars
  old <- graphics::par(mfrow = c(min(length(vars), 4), 2), mar = c(4, 4, 2, 1))
  on.exit(graphics::par(old))
  for (v in vars) {
    a <- unclass(arr[, , v])[, , 1]
    graphics::matplot(a, type = "l", lty = 1, xlab = "iteration", ylab = v,
                      main = paste("Trace:", v))
    graphics::plot(stats::density(as.vector(a)), main = paste("Density:", v), xlab = v)
  }
  invisible(x)
}

#' @rdname psprobit-methods
#' @exportS3Method posterior::as_draws
as_draws.psprobit <- function(x, ...) x$draws

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @rdname psprobit-methods
#' @param conf.level Credible level for `tidy()`.
#' @export
tidy.psprobit <- function(x, conf.level = 0.95, ...) {
  ci <- confint.psprobit(x, level = conf.level)
  sds <- sqrt(diag(vcov.psprobit(x)))
  data.frame(term = names(x$coefficients), estimate = unname(x$coefficients),
             std.error = unname(sds[names(x$coefficients)]),
             conf.low = unname(ci[, 1]), conf.high = unname(ci[, 2]),
             rhat = x$diagnostics$rhat, ess_bulk = x$diagnostics$ess_bulk,
             row.names = NULL)
}

#' @rdname psprobit-methods
#' @export
glance.psprobit <- function(x, ...) {
  data.frame(family = x$family, model = x$model, effects = x$effects, nobs = x$N * x$T,
             n_units = x$N, n_periods = x$T, chains = posterior::nchains(x$draws),
             draws = posterior::ndraws(x$draws),
             max_rhat = suppressWarnings(max(x$diagnostics$rhat, na.rm = TRUE)),
             min_ess_bulk = min(x$diagnostics$ess_bulk, na.rm = TRUE))
}
