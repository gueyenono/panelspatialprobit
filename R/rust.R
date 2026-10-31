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

# Which sampler engine runs a model: Rust (plan/06 R2), unless
# options(panelspatialprobit.engine = "cpp") asks for the C++ engine (for
# comparisons while both exist). R5 removes the C++ engine and this switch.
.psp_engine <- function() {
  forced <- getOption("panelspatialprobit.engine")
  if (is.null(forced)) "rust" else match.arg(forced, c("rust", "cpp"))
}
