# Replication of the simulation study of Wang & Kockelman (2009; working
# paper version: literature/Wang_Kockelman_2008_dynamic_spatial_ordered_probit_TRB.pdf,
# section 5 and Table 1): a dynamic spatial ordered probit with latent AR(1)
# dynamics and spatially autocorrelated effects.
# Usage: Rscript validation/wang_kockelman/run_wk.R [reps] [cores] [lib]
# Writes validation/wang_kockelman/results/wk_<date>.{rds,txt}.
#
# Their design: U_t = lambda U_{t-1} + X_t beta + theta + eps_t, with S = 3
# categories, thresholds gamma_1 = 0 and gamma_2 = 2.1, no intercept,
# beta = (-1.7, 2, 1, 0.5), four U(0, 1) regressors, T = 8, regional effects
# theta = (I - rho W)^{-1} u with u ~ N(0, sigma^2 = 1), lambda in
# {0.1, 0.5, 0.9} x rho in {0.1, 0.6, 0.7, 0.9}, 50 data sets per cell,
# 2,000 draws with 1,000 burn-in.
#
# Differences, because panelspatialprobit has one random effect per unit and
# unit error variance:
#   * one individual per region (the case M = N their section 3.2 allows):
#     N = 300 units on a 15 x 20 queen-contiguity grid, row-standardised
#     (theirs: 30 regions of 10 individuals on a 5 x 6 grid);
#   * homoskedastic errors with variance 1 (theirs: variances 0.4-1.5 by
#     region, region 1 fixed for identification), a special case of their V;
#   * the first period has no lagged term (their U_0 is not specified);
#   * the fitted model has an intercept with the first cut-point fixed at 0;
#     their thresholds are gamma_j = cut_j - intercept, reported as gamma1
#     and gamma2 (truth 0 and 2.1);
#   * longer chains (2 chains of 15,000 after 5,000 warm-up, thinned by 5),
#     because lambda and rho mix slowly when both are large.
# Their rho is our phi and their lambda our lag_z; sigma2 = sigma_alpha^2.
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 50L
cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
lib <- if (length(args) >= 3) args[3] else NULL
library(panelspatialprobit, lib.loc = lib)

queen <- function(nr, nc) {
  g <- expand.grid(r = seq_len(nr), c = seq_len(nc))
  B <- outer(seq_len(nrow(g)), seq_len(nrow(g)), function(i, j) {
    pmax(abs(g$r[i] - g$r[j]), abs(g$c[i] - g$c[j])) == 1
  }) * 1
  dimnames(B) <- list(seq_len(nrow(g)), seq_len(nrow(g)))
  psp_weights(B, style = "W")
}
W <- queen(15, 20)
N <- nrow(W)
T <- 8
beta <- c(-1.7, 2, 1, 0.5)
truth_common <- c(x1 = -1.7, x2 = 2, x3 = 1, x4 = 0.5, sigma2 = 1, gamma1 = 0, gamma2 = 2.1)

sim_wk <- function(lambda, rho, seed) {
  set.seed(seed)
  alpha <- solve(diag(N) - rho * as.matrix(W), stats::rnorm(N))
  X <- matrix(stats::runif(N * T * 4), N * T, 4)
  z <- numeric(N * T)
  for (t in seq_len(T)) {
    r <- (t - 1) * N + seq_len(N)
    z[r] <- (if (t > 1) lambda * z[r - N] else 0) + X[r, ] %*% beta + alpha + stats::rnorm(N)
  }
  data.frame(id = rep(seq_len(N), T), time = rep(seq_len(T), each = N),
             y = cut(z, c(-Inf, 0, 2.1, Inf), labels = FALSE),
             x1 = X[, 1], x2 = X[, 2], x3 = X[, 3], x4 = X[, 4])
}

one_rep <- function(job) {
  d <- sim_wk(job$lambda, job$rho, job$seed)
  fit <- suppressWarnings(psoprobit(
    y ~ x1 + x2 + x3 + x4, d, c("id", "time"), W, model = "none",
    effects = "spatial_random", dynamic = psp_dynamic(latent = TRUE), impacts = FALSE,
    control = psp_control(iter = 15000, warmup = 5000, thin = 5, chains = 2, cores = 1,
                          seed = job$seed)))
  dm <- unclass(posterior::as_draws_matrix(fit$draws))
  draws <- cbind(dm[, c("x1", "x2", "x3", "x4", "lag_z", "phi")],
                 dm[, "sigma_alpha"]^2, -dm[, "(Intercept)"],
                 dm[, "cut[2|3]"] - dm[, "(Intercept)"])
  colnames(draws) <- c("x1", "x2", "x3", "x4", "lambda", "rho", "sigma2", "gamma1", "gamma2")
  truth <- c(truth_common, lambda = job$lambda, rho = job$rho)[colnames(draws)]
  data.frame(lambda_true = job$lambda, rho_true = job$rho, rep = job$rep,
             parameter = colnames(draws), true = truth, mean = colMeans(draws),
             lo90 = apply(draws, 2, stats::quantile, 0.05),
             hi90 = apply(draws, 2, stats::quantile, 0.95),
             max_rhat = max(fit$diagnostics$rhat, na.rm = TRUE),
             min_ess = min(fit$diagnostics$ess_bulk, na.rm = TRUE),
             share3 = mean(d$y == 3), row.names = NULL)
}

cells <- expand.grid(rho = c(0.1, 0.6, 0.7, 0.9), lambda = c(0.1, 0.5, 0.9))
jobs <- do.call(c, lapply(seq_len(nrow(cells)), function(k) {
  lapply(seq_len(reps), function(r) {
    list(lambda = cells$lambda[k], rho = cells$rho[k], rep = r, seed = 1e4 * k + r)
  })
}))
t0 <- Sys.time()
out <- parallel::mclapply(jobs, function(j) try(one_rep(j)), mc.cores = cores,
                          mc.preschedule = FALSE)
elapsed <- format(Sys.time() - t0, digits = 3)
ok <- !vapply(out, inherits, logical(1), "try-error")
res <- do.call(rbind, out[ok])

dir.create("validation/wang_kockelman/results", showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.Date())
saveRDS(res, sprintf("validation/wang_kockelman/results/wk_%s.rds", stamp))

# Wang & Kockelman Table 1 (averages of posterior means over 50 data sets),
# one row per cell in the order of 'cells' (rho fastest within lambda).
wk_table1 <- t(as.matrix(utils::read.csv("validation/wang_kockelman/wk_table1.csv",
                                         comment.char = "#")[, -(1:2)]))

report <- c(
  sprintf("Wang & Kockelman (2009) simulation replicated %s: %d data sets per cell, %d fits failed (%s)",
          stamp, reps, sum(!ok), elapsed),
  "N = 300 (queen grid), T = 8, S = 3; 2 chains x 15,000 draws after 5,000 warm-up, thinned by 5.",
  "Columns: true value; average posterior mean here (psp) and in W&K Table 1 (wk);",
  "RMSE of the posterior mean here; coverage of 90% intervals here.",
  "")
for (k in seq_len(nrow(cells))) {
  sub <- res[res$lambda_true == cells$lambda[k] & res$rho_true == cells$rho[k], ]
  first <- sub[!duplicated(sub$rep), ]
  report <- c(report, sprintf(
    "[lambda = %.1f, rho = %.1f] %d data sets; share y = 3: %.2f; max R-hat median %.3f (%d%% >= 1.05); min bulk ESS median %.0f",
    cells$lambda[k], cells$rho[k], nrow(first), mean(first$share3), stats::median(first$max_rhat),
    round(100 * mean(first$max_rhat >= 1.05)), stats::median(first$min_ess)))
  for (p in rownames(wk_table1)) {
    s <- sub[sub$parameter == p, ]
    err <- s$mean - s$true
    report <- c(report, sprintf("  %-7s true %6.3f  psp %7.3f  wk %7.3f  rmse %6.3f  cover90 %.2f",
                                p, s$true[1], mean(s$mean), wk_table1[p, k], sqrt(mean(err^2)),
                                mean(s$lo90 <= s$true & s$true <= s$hi90)))
  }
}
writeLines(report, sprintf("validation/wang_kockelman/results/wk_%s.txt", stamp))
cat(report, sep = "\n")
