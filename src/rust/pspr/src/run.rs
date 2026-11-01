//! Running pspcore tasks on worker threads from R (plan/06 §6.7): the main
//! thread checks for user interrupts every 100 ms and prints progress.

use crate::interrupt;
use pspcore::par::{self, Progress, RunError, TaskCtx};
use pspcore::PspError;
use std::time::Duration;

const POLL: Duration = Duration::from_millis(100);

/// Prints "<label> <task>: <pct>%" whenever a task crosses another 10% of
/// `total` progress units.
pub struct Printer {
    label: String,
    total: usize,
    printed: Vec<usize>,
    warmup: Option<usize>,
}

impl Printer {
    pub fn new(label: &str, tasks: usize, total: usize) -> Self {
        Self {
            label: label.to_string(),
            total: total.max(1),
            printed: vec![0; tasks],
            warmup: None,
        }
    }

    /// MCMC progress: "Chain c: x% (warm-up)" or "(sampling)", also for
    /// chains run in parallel.
    pub fn mcmc(chains: usize, warmup: usize, iter: usize) -> Self {
        let mut p = Self::new("Chain", chains, warmup + iter);
        p.warmup = Some(warmup);
        p
    }

    fn report(&mut self, p: &Progress) {
        for t in 0..p.len() {
            let tenths = (10 * p.get(t).min(self.total)) / self.total;
            if tenths > self.printed[t] {
                self.printed[t] = tenths;
                let msg = match self.warmup {
                    None => format!("{} {}: {}%", self.label, t + 1, 10 * tenths),
                    Some(w) => {
                        let it = tenths * (self.total / 10).max(1);
                        let phase = if it <= w { "warm-up" } else { "sampling" };
                        format!("{} {}: {:3}% ({phase})", self.label, t + 1, 10 * tenths)
                    }
                };
                savvy::io::r_print(&msg, true);
            }
        }
    }
}

/// Runs `task` for `0..n_tasks` on `threads` threads (see
/// `pspcore::par::run`), polling for interrupts and printing progress with
/// `printer` if given. An interrupt returns `savvy::Error::Aborted` after the
/// threads have stopped.
pub fn run_tasks<T, F>(
    n_tasks: usize,
    threads: usize,
    task: F,
    mut printer: Option<Printer>,
) -> savvy::Result<Vec<T>>
where
    T: Send,
    F: Fn(usize, &TaskCtx) -> Result<T, PspError> + Sync,
{
    let res = par::run(n_tasks, threads, POLL, task, |p| {
        if let Some(pr) = printer.as_mut() {
            pr.report(p);
        }
        interrupt::check()
    });
    match res {
        Ok(v) => Ok(v),
        Err(RunError::Task(e)) => Err(savvy::Error::new(e.to_string())),
        Err(RunError::Tick(e)) => Err(e),
    }
}
