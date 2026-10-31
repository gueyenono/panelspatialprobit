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
