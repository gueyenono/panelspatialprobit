# R side of the Rust interface (plan/06 §6.5): R unpacks its objects into
# plain vectors and lists, and Rust copies them into its own structures once
# per call (R-D10).

# A sparse matrix as the CSC parts Rust reads (0-based row indices, as in a
# dgCMatrix).
.psp_csc <- function(M) {
  M <- methods::as(methods::as(methods::as(M, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  list(n_row = nrow(M), n_col = ncol(M), p = M@p, i = M@i, x = as.double(M@x))
}

# Weights for T periods: the distinct matrices and, per period, the 0-based
# index of the matrix it uses. W is one matrix (time-invariant) or a list with
# one matrix per period; identical matrices are stored once, as in
# .psp_lndet_avg().
.psp_rs_weights <- function(W, T) {
  if (!is.list(W)) return(list(mats = list(.psp_csc(W)), period = rep(0L, T)))
  if (length(W) != T) stop("'W' must have one matrix per period.", call. = FALSE)
  uniq <- list()
  grp <- integer(T)
  for (t in seq_len(T)) {
    j <- Position(function(u) identical(u, W[[t]]), uniq)
    if (is.na(j)) {
      uniq <- c(uniq, W[t])
      j <- length(uniq)
    }
    grp[t] <- j - 1L
  }
  list(mats = lapply(uniq, .psp_csc), period = grp)
}

# The data, weights and flags the Rust sampler reads (plan/06 §6.5). Missing
# outcomes are already coded in y (0 ordered, -1 binary).
.psp_rs_model <- function(design, y, ordered, has_rho, has_lambda, has_re, spatial_re) {
  X <- design$X
  list(y = as.integer(y), X = as.double(X), k = ncol(X), N = as.integer(design$N),
       T = as.integer(design$T), J = as.integer(design$J), ordered = ordered,
       has_rho = has_rho, has_lambda = has_lambda, has_re = has_re, spatial_re = spatial_re,
       W = .psp_rs_weights(design$W, design$T),
       M = .psp_rs_weights(if (has_lambda) design$M else design$W, design$T),
       lag_col = if ("lag_z" %in% colnames(X)) match("lag_z", colnames(X)) - 1L else -1L)
}

# Chains (or posterior draws) run at once: `cores`, capped at 2 under
# _R_CHECK_LIMIT_CORES_, as parallel::detectCores() users are on CRAN.
.psp_threads <- function(cores) {
  cores <- max(1L, as.integer(cores))
  limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(limit) && !identical(toupper(limit), "FALSE")) cores <- min(cores, 2L)
  cores
}

# A fit's data for the Rust post-processing (plan/06 §6.8): the design,
# weights, Rademacher probes (n_probes = 0: exact diagonals) and the N x Q
# mixing variances of t errors (n_mix = 0: normal errors).
.psp_rs_post_model <- function(fit, probes = .psp_probes(fit), mix = .psp_t_mix(fit)) {
  list(X = as.double(fit$X), k = ncol(fit$X), N = as.integer(fit$N), T = as.integer(fit$T),
       J = as.integer(if (inherits(fit, "psoprobit")) fit$J else 2L),
       ordered = inherits(fit, "psoprobit"),
       W = .psp_rs_weights(fit$W, fit$T),
       M = .psp_rs_weights(if (is.null(fit$M)) fit$W else fit$M, fit$T),
       probes = as.double(if (is.null(probes)) numeric(0) else probes),
       n_probes = if (is.null(probes)) 0L else ncol(probes),
       mix = as.double(if (is.null(mix)) numeric(0) else mix),
       n_mix = if (is.null(mix)) 0L else ncol(mix),
       df = if (is.null(fit$errors) || fit$errors$type == "normal") NA_real_
            else as.double(fit$errors$df))
}

# 6 x n integer matrix of L'Ecuyer-CMRG seeds for Rust: n consecutive
# streams from 'seed' (as .psp_streams()).
.psp_rs_seeds <- function(seed, n) {
  matrix(vapply(.psp_streams(seed, n), function(s) s[-1], integer(6)), 6)
}

# Calls whose per-draw seeds travel with the draws ('draws$seeds', 6 x draws),
# so that .psp_post_map() splits them with the draws.
.psp_rs_ppc_call <- function(model, draws) psp_rs_ppc(model, draws, draws$seeds)
.psp_rs_response_call <- function(model, draws, spec) {
  psp_rs_response(model, draws, draws$seeds, spec)
}

# Posterior draws for Rust, one column per draw: 'pm' is the draws x
# parameters matrix and 'alpha' (optional) the draws x N random effects.
.psp_rs_draws <- function(fit, pm, alpha = NULL) {
  D <- nrow(pm)
  col <- function(nm) if (nm %in% colnames(pm)) as.double(pm[, nm]) else numeric(D)
  cuts <- grepl("^cut\\[", colnames(pm))
  list(beta = t(pm[, colnames(fit$X), drop = FALSE]),
       sig2 = if (fit$effects != "pooled") col("sigma_alpha")^2 else numeric(D),
       rho = col("rho"), lambda = col("lambda"), phi = col("phi"), gamma = col("lag_z"),
       kappa = t(pm[, cuts, drop = FALSE]),
       alpha = if (is.null(alpha)) matrix(numeric(0), 0, D) else t(alpha))
}

.psp_rs_draws_subset <- function(draws, idx) {
  lapply(draws, function(x) if (is.matrix(x)) x[, idx, drop = FALSE] else x[idx])
}

# Runs fun(model, draws, ...) (a psp_rs_* post-processing function) over the
# posterior draws, spread over the user's future::plan() (R-D11): the draws
# are cut into blocks of 'block' draws, and each future gets a run of
# consecutive blocks. Results do not depend on the plan: per-draw results are
# the same wherever they run, and functions that sum over draws (predict:
# block = 8, pspcore's CHUNK) return block sums that the caller adds in block
# order. Returns fun's results, one per future, in draw order.
.psp_post_map <- function(fun, model, draws, ..., block = 1L) {
  n <- length(draws$sig2)
  blocks <- ceiling(n / block)
  workers <- if (inherits(future::plan(), "sequential")) 1L else future::nbrOfWorkers()
  groups <- max(1L, min(workers, blocks))
  if (groups == 1L) return(list(fun(model, draws, ...)))
  edges <- unique(pmin(round(seq(0, blocks, length.out = groups + 1)) * block, n))
  args <- list(...)
  futs <- lapply(seq_len(length(edges) - 1L), function(g) {
    sub <- .psp_rs_draws_subset(draws, seq.int(edges[g] + 1, edges[g + 1]))
    future::future(do.call(fun, c(list(model, sub), args)),
                   globals = list(fun = fun, model = model, sub = sub, args = args),
                   packages = "panelspatialprobit", seed = FALSE)
  })
  future::value(futs)
}
