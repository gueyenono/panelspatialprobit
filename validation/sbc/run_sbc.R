# Simulation-based calibration (Talts et al. 2018; Modrak et al. 2025) for the
# sampler variants implemented so far. Usage:
#   Rscript validation/sbc/run_sbc.R [reps] [cores] [lib] [variant1,variant2,...]
# The environment variable SBC_SEED_OFFSET (default 0) shifts every seed, for
# an independent replication.
# Writes validation/sbc/results/sbc_<date>.rds and a text report.
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 200L
cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
lib <- if (length(args) >= 3) args[3] else NULL
library(panelspatialprobit, lib.loc = lib)

W <- psp_weights(rook_weights(4, 4), style = "W")
T <- 4
# time-varying weights: rook links, a different fifth of them removed per period
W_tv <- local({
  set.seed(424242)
  B <- as.matrix(rook_weights(4, 4))
  lapply(seq_len(T), function(t) {
    Bt <- B
    links <- which(upper.tri(Bt) & Bt > 0)
    Bt[sample(links, round(length(links) / 5))] <- 0
    Bt[lower.tri(Bt)] <- t(Bt)[lower.tri(Bt)]
    psp_weights(Bt, style = "W")
  })
})
priors <- psp_priors(beta = prior_normal(0, 1), rho = prior_beta(2, 2), lambda = prior_beta(2, 2),
                     sigma_alpha = prior_inv_gamma(3, 1), cutpoints = prior_normal(1, 1),
                     phi = prior_beta(2, 2))
ctrl <- function(seed, sampler = "transform", mult = 1) {
  psp_control(iter = 4000 * mult, warmup = 1000 * mult, thin = 20 * mult, chains = 2,
              seed = seed, cutpoint_sampler = sampler)
}

seed_offset <- as.integer(Sys.getenv("SBC_SEED_OFFSET", "0"))
only <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else NULL
variants <- list(
  sar_re_bin = list(model = "sar", effects = "random", J = 2),
  sar_pool_bin = list(model = "sar", effects = "pooled", J = 2),
  none_re_bin = list(model = "none", effects = "random", J = 2),
  sdm_re_bin = list(model = "sdm", effects = "random", J = 2),
  sar_re_ord3 = list(model = "sar", effects = "random", J = 3),
  sar_re_ord4 = list(model = "sar", effects = "random", J = 4),
  sar_re_ord4_ac = list(model = "sar", effects = "random", J = 4, sampler = "albert_chib"),
  slx_re_bin = list(model = "slx", effects = "random", J = 2),
  sar_mundlak_bin = list(model = "sar", effects = "random", J = 2, cre = ~ x1),
  sar_chamberlain_bin = list(model = "sar", effects = "random", J = 2,
                             cre = list(type = "chamberlain", vars = ~ x1)),
  sem_re_bin = list(model = "sem", effects = "random", J = 2),
  sem_pool_bin = list(model = "sem", effects = "pooled", J = 2),
  sdem_re_bin = list(model = "sdem", effects = "random", J = 2),
  sarar_re_bin = list(model = "sarar", effects = "random", J = 2),
  sem_re_ord3 = list(model = "sem", effects = "random", J = 3),
  sar_dyn_bin = list(model = "sar", effects = "random", J = 2, dynamic = TRUE),
  sar_re_bin_missing = list(model = "sar", effects = "random", J = 2, missing = 0.2),
  sar_re_ord3_missing = list(model = "sar", effects = "random", J = 3, missing = 0.2),
  sar_tv_re_bin = list(model = "sar", effects = "random", J = 2, tv = TRUE),
  sarar_tv_re_ord3 = list(model = "sarar", effects = "random", J = 3, tv = TRUE),
  sarar_re_ord3 = list(model = "sarar", effects = "random", J = 3),
  sar_sre_bin = list(model = "sar", effects = "spatial_random", J = 2),
  sem_sre_ord3 = list(model = "sem", effects = "spatial_random", J = 3),
  sar_t_re_bin = list(model = "sar", effects = "random", J = 2, df = 4),
  sarar_t_re_ord3 = list(model = "sarar", effects = "random", J = 3, df = 4),
  sem_t_re_ord3 = list(model = "sem", effects = "random", J = 3, df = 4),
  sarar_t_re_ord3_long = list(model = "sarar", effects = "random", J = 3, df = 4, mult = 2),
  sar_lat_re_bin = list(model = "sar", effects = "random", J = 2, latent = TRUE),
  sem_lat_re_ord3 = list(model = "sem", effects = "random", J = 3, latent = TRUE),
  sar_lat_re_bin_long = list(model = "sar", effects = "random", J = 2, latent = TRUE, mult = 3)
)

# Dynamic model with Wooldridge initial conditions, generated exactly as
# estimated: periods 2..T+1 are modelled given the period-1 outcome, with
# columns (Intercept), x1, x2, lag_y, init_y, mean_x1, mean_x2.
sim_dynamic <- function(th) {
  N <- nrow(W)
  Tt <- T + 1
  d <- sim_psprobit(W, Tt, beta = c(0, 0, 0), rho = 0, sigma_alpha = 0)
  X <- cbind(1, d$x1, d$x2)
  est <- d$time > 1
  y <- integer(N * Tt)
  y[!est] <- as.integer(rnorm(N) + d$x1[!est] > 0)
  mx <- cbind(tapply(d$x1[est], d$id[est], mean)[rownames(W)],
              tapply(d$x2[est], d$id[est], mean)[rownames(W)])
  beta <- rnorm(7)
  S <- Matrix::Diagonal(N) - th$rho * W
  alpha <- rnorm(N, 0, th$sigma_alpha)
  for (t in 2:Tt) {
    r <- (t - 1) * N + seq_len(N)
    mu <- X[r, ] %*% beta[1:3] + beta[4] * y[r - N] + beta[5] * y[seq_len(N)] +
      mx %*% beta[6:7] + alpha + rnorm(N)
    y[r] <- as.integer(as.vector(Matrix::solve(S, mu)) > 0)
  }
  d$y <- y
  if (length(unique(y[est])) < 2) return(NULL)  # redraw: one outcome only
  des <- panelspatialprobit:::.psp_design(y ~ x1 + x2, d, c("id", "time"), W, "sar", NULL, NULL,
                                          FALSE, "binary", NULL, psp_dynamic())
  stopifnot(identical(colnames(des$X), c("(Intercept)", "x1", "x2", "lag_y", "init_y",
                                         "mean_x1", "mean_x2")),
            isTRUE(all.equal(unname(des$X[, 6:7]), unname(mx[rep(seq_len(N), T), ]))))
  names(beta) <- colnames(des$X)
  list(d = d, beta = beta)
}

# Data from the model with the full design (Durbin and CRE columns built by
# the package's own design code); 'beta' covers every design column.
sim_design <- function(v, th) {
  d <- sim_psprobit(W, T, beta = c(0, 0, 0), rho = 0, sigma_alpha = 0)
  d$y <- rep(0:1, length.out = nrow(d))
  des <- panelspatialprobit:::.psp_design(y ~ x1 + x2, d, c("id", "time"), W, v$model,
                                          v$cre, NULL, FALSE, "binary")
  stopifnot(identical(des$row_order, seq_len(nrow(d))))
  beta <- rnorm(ncol(des$X))
  names(beta) <- colnames(des$X)
  N <- nrow(W)
  S <- Matrix::Diagonal(N) - th$rho * W
  alpha <- rnorm(N, 0, th$sigma_alpha)
  mu <- matrix(des$X %*% beta, N, T) + alpha + matrix(rnorm(N * T), N, T)
  z <- as.vector(as.matrix(Matrix::solve(S, mu)))
  d$y <- as.integer(z > 0)
  list(d = d, beta = beta)
}

draw_prior <- function(v) {
  beta <- rnorm(3)
  rho <- if (v$model %in% c("sar", "sdm", "sarar")) 2 * rbeta(1, 2, 2) - 1 else 0
  lambda <- if (v$model %in% c("sem", "sdem", "sarar")) 2 * rbeta(1, 2, 2) - 1 else 0
  sig <- if (v$effects != "pooled") sqrt(1 / rgamma(1, 3, 1)) else 0
  phi <- if (v$effects == "spatial_random") 2 * rbeta(1, 2, 2) - 1 else 0
  theta <- if (v$model %in% c("sdm", "slx", "sdem")) rnorm(2) else NULL
  cut <- NULL
  if (v$J > 2) {
    repeat {
      cut <- rnorm(v$J - 2, 1, 1)
      if (all(diff(c(0, cut)) > 0)) break
    }
  }
  gamma_z <- if (isTRUE(v$latent)) rnorm(1) else 0
  list(beta = beta, rho = rho, lambda = lambda, sigma_alpha = sig, phi = phi, theta = theta,
       cutpoints = cut, gamma_z = gamma_z)
}

one_rep <- function(rep, vname) {
  v <- variants[[vname]]
  Wv <- if (isTRUE(v$tv)) W_tv else W
  set.seed(1e5 * match(vname, names(variants)) + rep + 1e7 * seed_offset)
  repeat {
    th <- draw_prior(v)
    if (isTRUE(v$dynamic)) {
      sd_ <- sim_dynamic(th)
      if (is.null(sd_)) next
      d <- sd_$d
      th$beta <- sd_$beta
    } else if (!is.null(v$cre)) {
      sd_ <- sim_design(v, th)
      d <- sd_$d
      th$beta <- sd_$beta
    } else {
      d <- sim_psprobit(Wv, T, beta = th$beta, rho = th$rho, sigma_alpha = th$sigma_alpha,
                        theta = th$theta, cutpoints = th$cutpoints, lambda = th$lambda,
                        phi = th$phi, df = v$df, gamma_z = th$gamma_z)
    }
    if (!is.null(v$missing)) d$y[runif(nrow(d)) < v$missing] <- NA  # missing at random
    J_obs <- length(unique(na.omit(if (isTRUE(v$dynamic)) d$y[d$time > 1] else d$y)))
    if (J_obs == max(v$J, 2)) break
  }
  f <- if (v$J > 2) psoprobit else psprobit
  fit <- suppressWarnings(f(y ~ x1 + x2, d, c("id", "time"), Wv, model = v$model,
                            effects = v$effects, cre = v$cre, priors = priors, impacts = FALSE,
                            dynamic = if (isTRUE(v$dynamic)) psp_dynamic()
                                      else if (isTRUE(v$latent)) psp_dynamic(latent = TRUE),
                            errors = if (is.null(v$df)) "normal" else psp_t(v$df),
                            control = ctrl(rep, if (is.null(v$sampler)) "transform"
                                           else v$sampler,
                                           if (is.null(v$mult)) 1 else v$mult)))
  dm <- posterior::as_draws_matrix(fit$draws)
  truth <- c(th$beta, if (!is.null(th$theta)) th$theta, if (isTRUE(v$latent)) th$gamma_z,
             if (v$model %in% c("sar", "sdm", "sarar")) th$rho,
             if (v$model %in% c("sem", "sdem", "sarar")) th$lambda,
             if (v$effects != "pooled") th$sigma_alpha,
             if (v$effects == "spatial_random") th$phi, th$cutpoints)
  stopifnot(length(truth) == ncol(dm))
  ranks <- colSums(sweep(unclass(dm), 2, truth, "<"))
  c(ranks, max_rhat = max(fit$diagnostics$rhat, na.rm = TRUE))
}

res <- list()
for (vn in if (is.null(only)) names(variants) else only) {
  t0 <- Sys.time()
  out <- parallel::mclapply(seq_len(reps), one_rep, vname = vn, mc.cores = cores)
  ok <- !vapply(out, inherits, logical(1), "try-error")
  res[[vn]] <- do.call(rbind, out[ok])
  cat(sprintf("%-15s %d reps (%d failed) in %s\n", vn, sum(ok), sum(!ok),
              format(Sys.time() - t0, digits = 3)))
}

dir.create("validation/sbc/results", showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.Date())
tag <- paste0(if (is.null(only)) "" else paste0("_", paste(only, collapse = "-")),
              if (seed_offset != 0) paste0("_offset", seed_offset))
saveRDS(res, sprintf("validation/sbc/results/sbc_%s%s.rds", stamp, tag))

# Uniformity: chi-square on 10 rank bins (ranks are in 0..L, L = 400 draws)
L <- 400
report <- c(sprintf(paste("SBC report %s: %d reps per variant, %d posterior draws, N = 16,",
                          "T = %d, seed offset %d"), stamp, reps, L, T, seed_offset), "")
for (vn in names(res)) {
  r <- res[[vn]]
  par_cols <- setdiff(colnames(r), "max_rhat")
  report <- c(report, sprintf("[%s] median max R-hat %.3f", vn, median(r[, "max_rhat"])))
  for (p in par_cols) {
    bins <- cut(r[, p], breaks = seq(-0.5, L + 0.5, length.out = 11), labels = FALSE)
    pv <- suppressWarnings(chisq.test(tabulate(bins, 10))$p.value)
    report <- c(report, sprintf("  %-14s chi2 p = %.3f  mean rank/L = %.3f", p, pv,
                                mean(r[, p]) / L))
  }
}
writeLines(report, sprintf("validation/sbc/results/sbc_%s%s.txt", stamp, tag))
cat(report, sep = "\n")
