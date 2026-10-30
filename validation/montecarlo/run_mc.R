# Monte Carlo study (plan/04, roadmap 3.6 and 5.4): frequentist properties of
# the posterior mean and 90% credible intervals at fixed true parameters.
# Usage: Rscript validation/montecarlo/run_mc.R [reps] [cores] [lib]
# Writes validation/montecarlo/results/mc_<date>.{rds,txt}.
#
# Designs:
#   sar_bin_N100, sar_bin_N400: binary SAR with random effects, T = 10;
#   sar_ord4_N100: ordered SAR (J = 4) with random effects, T = 10;
#   dyn_*: dynamic SAR data (state dependence gamma = 0.8, random effects
#     sigma_alpha = 0.8, T = 10 after the initial period) fitted three ways:
#     with Wooldridge initial conditions and random effects (dyn_wooldridge),
#     pooled without initial conditions (dyn_pooled, which confounds state
#     dependence with heterogeneity), and as a static random-effects model
#     (dyn_static, which omits the state dependence).
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 100L
cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
lib <- if (length(args) >= 3) args[3] else NULL
library(panelspatialprobit, lib.loc = lib)

ctrl <- function(seed) psp_control(iter = 2000, warmup = 1000, chains = 2, seed = seed,
                                   impacts_draws = 200)
grid <- function(n) psp_weights(rook_weights(n, n), style = "W")
truth_static <- list(beta = c(0, 1, -1), rho = 0.5, sigma_alpha = 0.7)

# true population-averaged total impact of x1 on P(y = 1) (binary) or on
# P(y = J) (ordered), evaluated on the fitted design
true_impact <- function(fit, th) {
  par <- c(stats::setNames(th$beta, colnames(fit$X)), rho = th$rho,
           sigma_alpha = th$sigma_alpha)
  if (!is.null(th$cutpoints)) {
    par <- c(par, stats::setNames(th$cutpoints, grep("^cut", posterior::variables(fit$draws),
                                                     value = TRUE)))
  }
  fit$draws <- posterior::as_draws_array(array(par, c(1, 1, length(par)),
                                               dimnames = list(NULL, NULL, names(par))))
  imp <- psp_impacts(fit, draws = 1)
  s <- summary(imp)
  s <- s[s$variable == "x1" & s$effect == "total", ]
  s$mean[nrow(s)]
}

est_impact <- function(fit) {
  s <- summary(fit$impacts)
  s <- s[s$variable == "x1" & s$effect == "total", ]
  ci <- quantile(fit$impacts$draws[, "x1", nrow(s), "total"], c(0.05, 0.95))
  c(mean = s$mean[nrow(s)], lo = ci[[1]], hi = ci[[2]])
}

summarise_fit <- function(fit, truth, pars) {
  dm <- posterior::as_draws_matrix(fit$draws)
  out <- lapply(pars, function(p) {
    x <- as.vector(dm[, p])
    q <- quantile(x, c(0.05, 0.95))
    c(true = truth[[p]], mean = mean(x), lo = q[[1]], hi = q[[2]])
  })
  names(out) <- pars
  out
}

static_rep <- function(rep, n, cutpoints = NULL) {
  set.seed(1000 * n + rep + if (is.null(cutpoints)) 0 else 5e5)
  W <- grid(n)
  d <- sim_psprobit(W, T = 10, beta = truth_static$beta, rho = truth_static$rho,
                    sigma_alpha = truth_static$sigma_alpha, cutpoints = cutpoints)
  f <- if (is.null(cutpoints)) psprobit else psoprobit
  fit <- suppressWarnings(f(y ~ x1 + x2, d, c("id", "time"), W, control = ctrl(rep)))
  th <- c(truth_static, list(cutpoints = cutpoints))
  truth <- c(list(x1 = 1, x2 = -1, rho = 0.5, sigma_alpha = 0.7),
             if (!is.null(cutpoints)) list(`cut[2|3]` = cutpoints[1], `cut[3|4]` = cutpoints[2]))
  res <- summarise_fit(fit, truth, names(truth))
  res$impact_x1 <- c(true = true_impact(fit, th), est_impact(fit))
  res
}

dyn_rep <- function(rep) {
  set.seed(9e5 + rep)
  W <- grid(10)
  d <- sim_psprobit(W, T = 11, beta = c(-0.3, 1, -1), rho = 0.4, sigma_alpha = 0.8,
                    gamma = 0.8)
  fits <- list(
    wooldridge = psprobit(y ~ x1 + x2, d, c("id", "time"), W, dynamic = psp_dynamic(),
                          impacts = FALSE, control = ctrl(rep)),
    pooled = psprobit(y ~ x1 + x2, d, c("id", "time"), W, effects = "pooled",
                      dynamic = psp_dynamic(initial = "none"), impacts = FALSE,
                      control = ctrl(rep)),
    static = psprobit(y ~ x1 + x2, d[d$time > 1, ], c("id", "time"), W, impacts = FALSE,
                      control = ctrl(rep))
  )
  fits <- lapply(fits, suppressWarnings)
  out <- list()
  for (nm in names(fits)) {
    truth <- list(x1 = 1, x2 = -1, rho = 0.4)
    if (nm != "static") truth$lag_y <- 0.8
    out[[nm]] <- summarise_fit(fits[[nm]], truth, names(truth))
  }
  out
}

run <- function(f, ...) parallel::mclapply(seq_len(reps), f, ..., mc.cores = cores)
t0 <- Sys.time()
res <- list(
  sar_bin_N100 = run(static_rep, n = 10),
  sar_bin_N400 = run(static_rep, n = 20),
  sar_ord4_N100 = run(static_rep, n = 10, cutpoints = c(0.8, 1.6))
)
dyn <- run(dyn_rep)
for (nm in c("wooldridge", "pooled", "static")) {
  res[[paste0("dyn_", nm)]] <- lapply(dyn, function(r) if (is.list(r)) r[[nm]] else r)
}
elapsed <- format(Sys.time() - t0, digits = 3)

dir.create("validation/montecarlo/results", showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.Date())
saveRDS(res, sprintf("validation/montecarlo/results/mc_%s.rds", stamp))

report <- c(sprintf("Monte Carlo study %s: %d replications per design (%s)", stamp, reps,
                    elapsed),
            "Posterior mean: bias and RMSE; 90% equal-tailed credible intervals: coverage and width.",
            "")
for (dn in names(res)) {
  ok <- Filter(function(r) is.list(r) && !inherits(r, "try-error"), res[[dn]])
  report <- c(report, sprintf("[%s] %d successful replications", dn, length(ok)))
  for (p in names(ok[[1]])) {
    m <- do.call(rbind, lapply(ok, `[[`, p))
    err <- m[, "mean"] - m[, "true"]
    report <- c(report, sprintf("  %-12s true %7.3f  bias %7.3f  rmse %6.3f  cover90 %.2f  width %6.3f",
                                p, mean(m[, "true"]), mean(err), sqrt(mean(err^2)),
                                mean(m[, "lo"] <= m[, "true"] & m[, "true"] <= m[, "hi"]),
                                mean(m[, "hi"] - m[, "lo"])))
  }
}
writeLines(report, sprintf("validation/montecarlo/results/mc_%s.txt", stamp))
cat(report, sep = "\n")
