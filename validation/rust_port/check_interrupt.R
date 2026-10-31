# Interrupts during multi-threaded Rust work (plan/06 step 13; §6.10 layer 6).
# Forks a child R process (parallel::mcparallel, which also checks that the
# Rust code works after fork()) that runs a long dummy job on 4 threads,
# sends it SIGINT after 2 s, and reports how long the job took to stop,
# whether R saw an interrupt condition, how many threads the child had left,
# and that Rust code still runs afterwards. Unix only.
#
# (A child started with `&` from a non-interactive shell would not do:
# it inherits SIGINT as ignored, and R keeps that.)
# Usage: Rscript validation/rust_port/check_interrupt.R <lib>
library(panelspatialprobit, lib.loc = commandArgs(trailingOnly = TRUE)[1])
job <- parallel::mcparallel({
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(panelspatialprobit:::psp_rs_dummy_job(8L, 4L, 1000000L, 1:6, FALSE),
                interrupt = function(e) "interrupted")
  list(result = if (is.character(r)) r else "finished",
       stopped_at = proc.time()[["elapsed"]] - t0,
       threads = length(list.files("/proc/self/task")),
       again = length(panelspatialprobit:::psp_rs_dummy_job(2L, 2L, 5L, 1:6, FALSE)))
})
Sys.sleep(2)
t_signal <- proc.time()[["elapsed"]]
tools::pskill(job$pid, tools::SIGINT)
res <- parallel::mccollect(job)[[1]]
cat(sprintf("child: %s after %.2f s of work; %.3f s from SIGINT to the result\n",
            res$result, res$stopped_at, proc.time()[["elapsed"]] - t_signal))
cat(sprintf("threads left in the child: %d; a second job afterwards returned %d values\n",
            res$threads, res$again))
