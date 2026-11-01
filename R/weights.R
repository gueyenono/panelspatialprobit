#' Prepare a spatial weights matrix
#'
#' Coerces `W` to a sparse `dgCMatrix`, validates it and optionally
#' row-standardises it. Unit identifiers are taken from `ids`, or from the row
#' names of `W`, and are used by [psprobit()] to match rows of `W` to the units
#' in the data.
#'
#' @param W A square matrix, a `Matrix` sparse matrix, an `spdep` `listw` or
#'   `nb` object.
#' @param ids Optional vector of unit identifiers, one per row of `W`.
#' @param style `"asis"` (default) keeps the weights unchanged; `"W"`
#'   row-standardises them (rows of isolates stay zero).
#' @return A `dgCMatrix` with unit identifiers as dimnames.
#' @examples
#' W <- psp_weights(rook_weights(3, 3), style = "W")
#' Matrix::rowSums(W)
#' @export
psp_weights <- function(W, ids = NULL, style = c("asis", "W")) {
  style <- match.arg(style)
  if (inherits(W, "nb")) {
    if (!requireNamespace("spdep", quietly = TRUE)) {
      stop("Package 'spdep' is needed to convert 'nb' objects.", call. = FALSE)
    }
    W <- spdep::nb2listw(W, style = "B", zero.policy = TRUE)
  }
  if (inherits(W, "listw")) {
    if (!requireNamespace("spdep", quietly = TRUE)) {
      stop("Package 'spdep' is needed to convert 'listw' objects.", call. = FALSE)
    }
    nms <- attr(W$neighbours, "region.id")
    W <- spdep::listw2mat(W)
    if (is.null(ids) && !is.null(nms)) ids <- nms
  }
  if (is.null(ids)) ids <- rownames(W)
  W <- methods::as(methods::as(methods::as(Matrix::Matrix(W, sparse = TRUE), "dMatrix"),
                               "generalMatrix"), "CsparseMatrix")
  if (nrow(W) != ncol(W)) stop("'W' must be square.", call. = FALSE)
  if (anyNA(W@x) || any(!is.finite(W@x))) {
    stop("'W' contains missing or non-finite values.", call. = FALSE)
  }
  if (any(Matrix::diag(W) != 0)) stop("'W' must have a zero diagonal.", call. = FALSE)
  W <- Matrix::drop0(W)
  if (style == "W") {
    rs <- Matrix::rowSums(W)
    rs[rs == 0] <- 1
    W <- Matrix::Diagonal(x = 1 / rs) %*% W
    W <- methods::as(W, "CsparseMatrix")
  }
  if (!is.null(ids)) {
    if (length(ids) != nrow(W)) stop("'ids' must have one entry per row of 'W'.", call. = FALSE)
    if (anyDuplicated(ids)) stop("'ids' must be unique.", call. = FALSE)
    dimnames(W) <- list(as.character(ids), as.character(ids))
  }
  W
}

#' Rook contiguity weights on a regular grid
#'
#' Convenience constructor used in examples, tests and simulations.
#'
#' @param nrow,ncol Grid dimensions; units are numbered row by row.
#' @return A binary sparse `dgCMatrix` of dimension `nrow * ncol` with unit
#'   identifiers `"1"`, `"2"`, ... as dimnames.
#' @examples
#' rook_weights(2, 3)
#' @export
rook_weights <- function(nrow, ncol) {
  n <- nrow * ncol
  id <- matrix(seq_len(n), nrow, ncol, byrow = TRUE)
  pairs <- rbind(
    cbind(as.vector(id[, -ncol]), as.vector(id[, -1])),
    cbind(as.vector(id[-nrow, ]), as.vector(id[-1, ]))
  )
  W <- Matrix::sparseMatrix(i = c(pairs[, 1], pairs[, 2]), j = c(pairs[, 2], pairs[, 1]),
                            x = 1, dims = c(n, n),
                            dimnames = list(as.character(seq_len(n)), as.character(seq_len(n))))
  methods::as(W, "CsparseMatrix")
}

# Log-determinant log|I - rho W| on a fine grid strictly inside (rmin, rmax).
.psp_lndet <- function(W, rmin = -1, rmax = 1, method = c("auto", "eigen", "sparse"),
                       step = 5e-4) {
  method <- match.arg(method)
  N <- nrow(W)
  if (method == "auto") method <- if (N <= 3000) "eigen" else "sparse"
  grid <- seq(rmin + step / 2, rmax - step / 2, by = step)
  if (method == "eigen") {
    ev <- eigen(as.matrix(W), only.values = TRUE)$values
    lndet <- psp_rs_lndet_eigen(Re(ev), Im(ev), grid)
  } else {
    coarse <- unique(c(seq(rmin + step / 2, rmax - step / 2, length.out = 120),
                       rmax - c(0.02, 0.01, 0.005, 0.002) - step / 2))
    coarse <- sort(coarse[coarse > rmin & coarse < rmax])
    vals <- psp_rs_lndet_sparse(.psp_csc(W), coarse)  # one sparse LU per point
    lndet <- stats::splinefun(coarse, vals, method = "natural")(grid)
  }
  list(grid = grid, lndet = lndet, method = method)
}

# Time-varying weights are a list with one matrix per period; helpers below
# treat a single matrix and such a list alike.
.psp_is_weights_list <- function(W) is.list(W) && !inherits(W, c("nb", "listw"))

# Weights of period t (a single matrix serves every period)
.psp_Wt <- function(W, t) if (is.list(W)) W[[t]] else W

# Log-determinant grid averaged over periods, so that T * lndet is the sum of
# the per-period log-determinants; distinct matrices are computed once.
.psp_lndet_avg <- function(W, rmin, rmax, method) {
  if (!is.list(W)) return(.psp_lndet(W, rmin, rmax, method))
  uniq <- list()
  grp <- integer(length(W))
  for (t in seq_along(W)) {
    j <- Position(function(u) identical(u, W[[t]]), uniq)
    if (is.na(j)) {
      uniq <- c(uniq, W[t])
      j <- length(uniq)
    }
    grp[t] <- j
  }
  grids <- lapply(uniq, .psp_lndet, rmin = rmin, rmax = rmax, method = method)
  weight <- tabulate(grp, length(uniq)) / length(W)
  out <- grids[[1]]
  out$lndet <- Reduce(`+`, Map(function(g, w) w * g$lndet, grids, weight))
  out
}
