# Threads, interrupts and fork() with the Rust sampler (plan/06 step 24;
# §6.10 layer 6). Unix only.
#   1. Ctrl-C (SIGINT) during psprobit() with 2 chains at once, each on 2
#      threads, stops it promptly with an R interrupt condition and leaves no
#      threads behind.
#   2. Progress is printed for every chain when chains run in parallel.
#   3. psprobit() runs inside parallel::mclapply() (forked children), and the
#      draws equal those of the same fit in the parent.
#   4. The draws do not depend on `cores`, or on `threads_per_chain`
#      (period-parallel latent sweep, R6), which items 1 and 3 also use.
# The fits use the sparse log-determinant: eigen() starts the BLAS thread
# pool, and on this machine (nix BLAS) a forked child that calls eigen() after
# the parent did deadlocks, with or without this package.
# Usage: Rscript validation/rust_port/check_fit_threads.R <lib>
library(panelspatialprobit, lib.loc = commandArgs(trailingOnly = TRUE)[1])
W <- psp_weights(rook_weights(30, 30), style = "W")
d <- sim_psprobit(W, T = 10, beta = c(0, 1, -1), rho = 0.5, seed = 1)
fit <- function(iter, cores, progress = FALSE, chains = 4, per_chain = 1) {
  psprobit(y ~ x1 + x2, d, c("id", "time"), W, impacts = FALSE,
           control = psp_control(iter = iter, warmup = iter, chains = chains, cores = cores,
                                 threads_per_chain = per_chain, seed = 5,
                                 progress = progress, lndet = "sparse"))
}

# 1. interrupt
nthreads <- function() length(list.files("/proc/self/task"))
job <- parallel::mcparallel({
  before <- nthreads()
  r <- tryCatch(fit(1e6, 2, per_chain = 2), interrupt = function(e) "interrupted")
  list(result = if (is.character(r)) r else "finished", threads = nthreads() - before)
})
Sys.sleep(3)
t_signal <- proc.time()[["elapsed"]]
tools::pskill(job$pid, tools::SIGINT)
res <- parallel::mccollect(job)[[1]]
cat(sprintf("1. interrupt: %s, %.3f s from SIGINT to the result, %d thread(s) more than before\n",
            res$result, proc.time()[["elapsed"]] - t_signal, res$threads))

# 2. progress with parallel chains
out <- capture.output(invisible(suppressWarnings(fit(200, 4, progress = TRUE))))
cat(sprintf("2. progress: %d lines; each chain reports 100%%: %s\n", length(out),
            all(sprintf("Chain %d: 100%% (sampling)", 1:4) %in% out)))

# 3. inside mclapply (fork) and 4. cores
ref <- suppressWarnings(posterior::as_draws_matrix(fit(100, 1)$draws))
kids <- parallel::mclapply(1:2, function(i) {
  suppressWarnings(posterior::as_draws_matrix(fit(100, 2, per_chain = 2)$draws))
}, mc.cores = 2)
cat(sprintf("3. mclapply: both children equal the parent's fit: %s\n",
            all(vapply(kids, identical, logical(1), ref))))
cat(sprintf("4. cores = 4 equals cores = 1: %s; threads_per_chain = 3 equals 1: %s\n",
            identical(suppressWarnings(posterior::as_draws_matrix(fit(100, 4)$draws)), ref),
            identical(suppressWarnings(posterior::as_draws_matrix(fit(100, 1, per_chain = 3)$draws)),
                      ref)))
