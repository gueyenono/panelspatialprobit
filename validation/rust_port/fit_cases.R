# Child process of compare_sampler.R: fits the cases with one library and
# saves the draws and the per-chain summaries.
# Usage: Rscript fit_cases.R <lib> <out.rds> <iter> <warmup> <chains> [case ...]
a <- commandArgs(trailingOnly = TRUE)
library(panelspatialprobit, lib.loc = a[1])
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))),
                 "cases.R"))
cases <- psp_cases(as.integer(a[3]), as.integer(a[4]), as.integer(a[5]))
if (length(a) > 5) cases <- cases[a[-(1:5)]]
res <- lapply(cases, function(cs) {
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(suppressWarnings(cs$fit()), error = function(e) e)
  if (inherits(fit, "error")) return(list(error = conditionMessage(fit)))
  list(draws = posterior::as_draws_matrix(fit$draws), z_mean = fit$z_mean,
       alpha_mean = fit$alpha_mean, cut_accept = fit$cut_accept, v_mean = fit$v_mean,
       seconds = proc.time()[["elapsed"]] - t0)
})
saveRDS(res, a[2])
