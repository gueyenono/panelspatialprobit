# Build the internal design: observations sorted by (time, unit in W order),
# regressor columns tagged by role. See plan/03 s3.2.
.psp_design <- function(formula, data, index, W, model, cre, durbin, time_effects, family,
                        M = NULL, dynamic = NULL) {
  call_spec <- list(formula = formula, cre = cre, durbin = durbin, time_effects = time_effects)
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  if (length(index) != 2 || !all(index %in% names(data))) {
    stop("'index' must name the unit and time columns of 'data'.", call. = FALSE)
  }
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass,
                           drop.unused.levels = FALSE)
  mt <- attr(mf, "terms")
  y_raw <- stats::model.response(mf)
  X <- stats::model.matrix(mt, mf)
  id <- data[[index[1]]]
  tm <- data[[index[2]]]
  if (anyNA(X) || anyNA(id) || anyNA(tm)) {
    stop("Missing values in the regressors or index are not supported ",
         "(missing outcomes are, and are imputed).", call. = FALSE)
  }
  if (!is.null(dynamic) && !isTRUE(dynamic$latent) && anyNA(y_raw)) {
    stop("Missing outcomes are not supported in dynamic models.", call. = FALSE)
  }
  if (attr(mt, "intercept") != 1) {
    stop("The model must include an intercept (the first cut-point is fixed at 0).",
         call. = FALSE)
  }

  # units: match data to the rows of W (the first matrix of time-varying W)
  W_in <- W
  W <- psp_weights(if (.psp_is_weights_list(W)) W[[1]] else W)
  id_chr <- as.character(id)
  if (is.null(rownames(W))) {
    ids_W <- as.character(sort(unique(id)))
    if (length(ids_W) != nrow(W)) {
      stop(sprintf("'W' has %d rows but the data have %d units.", nrow(W), length(ids_W)),
           call. = FALSE)
    }
    warning("'W' has no unit names; assuming its rows follow the sorted unit ids. ",
            "Name the rows of 'W' (or use psp_weights(ids = )) to silence this.",
            call. = FALSE)
    dimnames(W) <- list(ids_W, ids_W)
  } else {
    ids_W <- rownames(W)
    missing_units <- setdiff(unique(id_chr), ids_W)
    if (length(missing_units)) {
      stop("Units not found in the row names of 'W': ",
           paste(utils::head(missing_units, 5), collapse = ", "), call. = FALSE)
    }
  }
  N <- nrow(W)
  t_levels <- sort(unique(tm))
  T <- length(t_levels)
  if (.psp_is_weights_list(W_in)) W <- .psp_align_weights(W_in, ids_W, t_levels, "W")
  M <- .psp_error_weights(M, W, model, ids_W, t_levels)
  unit_idx <- match(id_chr, ids_W)
  time_idx <- match(tm, t_levels)
  if (anyDuplicated(cbind(unit_idx, time_idx))) {
    stop("Duplicated (unit, time) pairs in 'data'.", call. = FALSE)
  }
  if (nrow(data) != N * T) {
    stop(sprintf(paste("The panel must be balanced: %d units in 'W' x %d periods",
                       "= %d rows expected, %d found."), N, T, N * T, nrow(data)),
         call. = FALSE)
  }
  ord <- order(time_idx, unit_idx)
  X <- X[ord, , drop = FALSE]
  y_raw <- y_raw[ord]
  resp <- .psp_code_response(y_raw, family)
  roles <- ifelse(colnames(X) == "(Intercept)", "intercept", "x")
  names(roles) <- colnames(X)
  base_vars <- colnames(X)[roles == "x"]

  # dynamics: lagged outcomes and initial conditions, first period dropped
  if (!is.null(dynamic) && !isTRUE(dynamic$latent)) {
    dyn <- .psp_dynamic_columns(dynamic, resp, N, T, W)
    if (is.list(W)) W <- W[-1]
    if (is.list(M)) M <- M[-1]
    keep <- seq.int(N + 1, N * T)
    X <- cbind(X[keep, , drop = FALSE], dyn$X)
    roles <- c(roles, dyn$roles)
    resp$y <- resp$y[keep]
    if (length(unique(resp$y)) < resp$J) {  # no missing outcomes here
      stop("Every outcome category must occur after the first period.", call. = FALSE)
    }
    ord <- ord[keep]
    T <- T - 1L
    t_levels <- t_levels[-1]
    if (dynamic$initial == "wooldridge") {
      cre <- .psp_add_cre_vars(cre, base_vars)
    }
  }

  blocks <- split(seq_len(N * T), rep(seq_len(T), each = N))

  # Durbin (WX) columns
  durbin_of <- character(0)
  if (model %in% c("sdm", "slx", "sdem")) {
    dvars <- if (is.null(durbin)) base_vars else .psp_formula_vars(durbin, base_vars, "durbin")
    if (!length(dvars)) stop("No regressors available for Durbin terms.", call. = FALSE)
    WX <- matrix(0, N * T, length(dvars), dimnames = list(NULL, paste0("W_", dvars)))
    for (t in seq_len(T)) {
      b <- blocks[[t]]
      WX[b, ] <- as.matrix(.psp_Wt(W, t) %*% X[b, dvars, drop = FALSE])
    }
    X <- cbind(X, WX)
    roles <- c(roles, stats::setNames(rep("durbin", length(dvars)), colnames(WX)))
    durbin_of <- stats::setNames(dvars, colnames(WX))
  } else if (!is.null(durbin)) {
    warning("'durbin' is ignored unless model is \"sdm\", \"slx\" or \"sdem\".",
            call. = FALSE)
  }

  # correlated random effects: Mundlak unit means or Chamberlain period values
  if (!is.null(cre)) {
    spec <- .psp_cre_spec(cre)
    cvars <- .psp_formula_vars(spec$vars, base_vars, "cre")
    unit_sorted <- rep(seq_len(N), T)
    if (spec$type == "mundlak") {
      M <- vapply(cvars, function(v) {
        stats::ave(X[, v], unit_sorted)
      }, numeric(N * T))
      M <- matrix(M, ncol = length(cvars), dimnames = list(NULL, paste0("mean_", cvars)))
    } else {
      # column (v, s) holds x_{is,v} for every period of unit i
      M <- do.call(cbind, lapply(cvars, function(v) {
        xs <- matrix(X[, v], N, T)
        out <- xs[unit_sorted, , drop = FALSE]
        colnames(out) <- paste0("cre_", v, "_", t_levels)
        out
      }))
    }
    X <- cbind(X, M)
    roles <- c(roles, stats::setNames(rep("cre", ncol(M)), colnames(M)))
  }

  # time effects
  if (time_effects) {
    if (T < 2) stop("Time effects need at least two periods.", call. = FALSE)
    D <- matrix(0, N * T, T - 1, dimnames = list(NULL, paste0("time_", t_levels[-1])))
    for (t in 2:T) D[blocks[[t]], t - 1] <- 1
    X <- cbind(X, D)
    roles <- c(roles, stats::setNames(rep("time", ncol(D)), colnames(D)))
  }

  # aliased columns
  q <- qr(X)
  if (q$rank < ncol(X)) {
    drop <- colnames(X)[q$pivot[seq.int(q$rank + 1, ncol(X))]]
    warning("Dropping aliased columns: ", paste(drop, collapse = ", "), call. = FALSE)
    keep <- !colnames(X) %in% drop
    X <- X[, keep, drop = FALSE]
    roles <- roles[keep]
    durbin_of <- durbin_of[names(durbin_of) %in% colnames(X)]
  }

  # latent AR(1): placeholder column for z_{t-1}, filled by the sampler
  if (isTRUE(dynamic$latent)) {
    if (is.list(W) || is.list(M)) {
      stop("Latent dynamics need time-invariant weights.", call. = FALSE)
    }
    if (T < 2) stop("Latent dynamics need at least two periods.", call. = FALSE)
    X <- cbind(X, lag_z = 0)
    roles <- c(roles, lag_z = "latent")
  }

  list(y = resp$y, levels = resp$levels, J = resp$J, X = X, roles = roles,
       durbin_of = durbin_of, W = W, M = M, N = N, T = T, unit_ids = ids_W,
       time_levels = t_levels, row_order = as.integer(rank(ord)), data_rows = ord, terms = mt,
       index = index, dynamic = dynamic, spec = call_spec)
}

# Error weights M for sem/sdem/sarar, aligned with the units of W
.psp_error_weights <- function(M, W, model, ids, t_levels) {
  if (!model %in% c("sem", "sdem", "sarar")) {
    if (!is.null(M)) {
      warning("'M' is ignored unless model is \"sem\", \"sdem\" or \"sarar\".",
              call. = FALSE)
    }
    return(NULL)
  }
  M <- if (is.null(M)) W else .psp_align_weights(M, ids, t_levels, "M")
  if (model == "sarar" && isTRUE(all.equal(M, W))) {
    warning("With M = W, rho and lambda are identified only through the regressors ",
            "and may be weakly identified.", call. = FALSE)
  }
  M
}

# Align weights (one matrix, or a list with one matrix per period, named by
# period or in period order) with the units 'ids' and the periods 't_levels'.
# A list of identical matrices collapses to a single matrix.
.psp_align_weights <- function(W, ids, t_levels, what) {
  align <- function(w) {
    w <- psp_weights(w)
    if (nrow(w) != length(ids)) {
      stop(sprintf("'%s' has %d rows but 'W' has %d.", what, nrow(w), length(ids)),
           call. = FALSE)
    }
    if (is.null(rownames(w))) {
      dimnames(w) <- list(ids, ids)
    } else {
      if (!setequal(rownames(w), ids)) {
        stop(sprintf("The row names of '%s' must match those of 'W'.", what), call. = FALSE)
      }
      w <- w[ids, ids]
    }
    w
  }
  if (!.psp_is_weights_list(W)) return(align(W))
  T <- length(t_levels)
  if (length(W) != T) {
    stop(sprintf("'%s' is a list of %d matrices but the data have %d periods.", what,
                 length(W), T), call. = FALSE)
  }
  if (!is.null(names(W))) {
    pos <- match(as.character(t_levels), names(W))
    if (anyNA(pos)) {
      stop(sprintf("The names of the '%s' list must be the periods of the data.", what),
           call. = FALSE)
    }
    W <- W[pos]
  }
  W <- lapply(W, align)
  if (all(vapply(W, identical, logical(1), W[[1]]))) return(W[[1]])
  stats::setNames(W, as.character(t_levels))
}

# Lagged outcome columns (periods 2..T) for a dynamic model. Binary outcomes
# give lag_y (and W_lag_y); ordered outcomes give one indicator per category
# above the first. Wooldridge (2005) initial conditions add the first-period
# outcome (y0 columns); the unit means of the regressors are added as CRE.
.psp_dynamic_columns <- function(dynamic, resp, N, T, W) {
  if (T < 3) stop("Dynamic models need at least three periods.", call. = FALSE)
  if (resp$J == 2) {
    Y <- matrix(resp$y, ncol = 1, dimnames = list(NULL, "y"))
  } else {
    Y <- outer(resp$y, 2:resp$J, "==") * 1
    colnames(Y) <- paste0("y", resp$levels[-1])
  }
  lag_rows <- seq_len(N * (T - 1))  # period t-1 rows for periods 2..T
  out <- list()
  roles <- character(0)
  if (dynamic$lag_y) {
    L <- Y[lag_rows, , drop = FALSE]
    colnames(L) <- paste0("lag_", colnames(Y))
    out$lag <- L
    roles <- c(roles, stats::setNames(rep("dynamic", ncol(L)), colnames(L)))
  }
  if (dynamic$lag_Wy) {
    WL <- Y[lag_rows, , drop = FALSE]
    for (t in seq_len(T - 1)) {
      r <- (t - 1) * N + seq_len(N)
      WL[r, ] <- as.matrix(.psp_Wt(W, t + 1) %*% Y[r, , drop = FALSE])
    }
    colnames(WL) <- paste0("W_lag_", colnames(Y))
    out$wlag <- WL
    roles <- c(roles, stats::setNames(rep("dynamic", ncol(WL)), colnames(WL)))
  }
  if (dynamic$initial == "wooldridge") {
    Y0 <- Y[rep(seq_len(N), T - 1), , drop = FALSE]
    colnames(Y0) <- paste0("init_", colnames(Y))
    out$init <- Y0
    roles <- c(roles, stats::setNames(rep("initial", ncol(Y0)), colnames(Y0)))
  }
  list(X = do.call(cbind, unname(out)), roles = roles)
}

# Wooldridge initial conditions: Mundlak means of every regressor
.psp_add_cre_vars <- function(cre, base_vars) {
  if (!length(base_vars)) return(cre)
  if (is.null(cre)) return(TRUE)
  spec <- .psp_cre_spec(cre)
  if (spec$type != "mundlak") return(cre)
  TRUE
}

# 'cre' is a formula / TRUE (Mundlak) or list(type = "mundlak" | "chamberlain", vars = )
.psp_cre_spec <- function(cre) {
  if (is.list(cre) && !inherits(cre, "formula")) {
    type <- match.arg(cre$type, c("mundlak", "chamberlain"))
    vars <- if (is.null(cre$vars)) TRUE else cre$vars
    return(list(type = type, vars = vars))
  }
  list(type = "mundlak", vars = cre)
}

.psp_formula_vars <- function(f, available, what) {
  if (isTRUE(f)) return(available)
  if (!inherits(f, "formula")) stop(sprintf("'%s' must be a one-sided formula or TRUE.", what),
                                    call. = FALSE)
  vars <- attr(stats::terms(f), "term.labels")
  bad <- setdiff(vars, available)
  if (length(bad)) {
    stop(sprintf("'%s' terms must be numeric regressors of the model: %s", what,
                 paste(bad, collapse = ", ")), call. = FALSE)
  }
  vars
}

.psp_code_response <- function(y, family) {
  if (family == "binary") {
    if (is.logical(y)) y <- as.integer(y)
    if (is.factor(y)) {
      if (nlevels(y) != 2) stop("A factor response must have two levels.", call. = FALSE)
      lev <- levels(y)
      y <- as.integer(y) - 1L
    } else {
      if (!all(y[!is.na(y)] %in% c(0, 1))) {
        stop("The response must be 0/1, logical or a two-level factor.", call. = FALSE)
      }
      lev <- c("0", "1")
      y <- as.integer(y)
    }
    if (length(unique(y[!is.na(y)])) < 2) {
      stop("The response has only one outcome.", call. = FALSE)
    }
    return(list(y = y, levels = lev, J = 2L))
  }
  if (is.factor(y)) {
    lev <- levels(y)
    y <- as.integer(y)
  } else {
    if (!is.numeric(y) || any(y != round(y), na.rm = TRUE)) {
      stop("An ordered response must be an (ordered) factor or integer codes.", call. = FALSE)
    }
    lev <- as.character(sort(unique(y)))
    y <- match(y, sort(unique(y)))
  }
  J <- length(lev)
  counts <- tabulate(y, J)
  if (any(counts == 0)) {
    stop("Empty response categories: ", paste(lev[counts == 0], collapse = ", "),
         call. = FALSE)
  }
  if (J < 3) stop("An ordered response needs at least 3 categories; use psprobit() for binary outcomes.",
                  call. = FALSE)
  list(y = as.integer(y), levels = lev, J = as.integer(J))
}
