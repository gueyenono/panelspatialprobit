//! Scoped worker threads (plan/06 §6.7, R-D4).
//!
//! Threads never outlive a call: [`run`] spawns them with
//! `std::thread::scope` and joins them before returning, so no thread pool
//! survives into a later `fork()` (e.g. `parallel::mclapply()`). Tasks are
//! taken from a shared queue; each task's result depends only on its index,
//! never on which thread ran it or on the thread count. The calling thread
//! does no task work: it polls `tick` (for R, the interrupt check and the
//! progress printing) until every task has finished.

use crate::PspError;
use std::ops::Range;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::mpsc::{self, RecvTimeoutError};
use std::time::Duration;

/// Draws per chunk for reductions over posterior draws (plan/06 §6.8):
/// partial results are combined in chunk order, so sums do not depend on the
/// thread count.
pub const CHUNK: usize = 8;

/// `0..n` split into consecutive ranges of `size` (the last may be shorter).
pub fn chunks(n: usize, size: usize) -> Vec<Range<usize>> {
    let size = size.max(1);
    (0..n.div_ceil(size))
        .map(|c| c * size..((c + 1) * size).min(n))
        .collect()
}

/// What a running task sees: the cancel flag and its own progress counter.
pub struct TaskCtx<'a> {
    cancel: &'a AtomicBool,
    progress: &'a AtomicUsize,
}

impl TaskCtx<'_> {
    pub fn cancelled(&self) -> bool {
        self.cancel.load(Ordering::Relaxed)
    }

    /// `Err(PspError::Cancelled)` once the run is cancelled. Long tasks call
    /// this regularly (the sampler: every 100 iterations).
    pub fn check(&self) -> Result<(), PspError> {
        if self.cancelled() {
            Err(PspError::Cancelled)
        } else {
            Ok(())
        }
    }

    /// Records how much of the task is done, in the task's own units.
    pub fn set_progress(&self, done: usize) {
        self.progress.store(done, Ordering::Relaxed);
    }
}

/// The progress counters of all tasks, read by the calling thread.
pub struct Progress<'a> {
    counters: &'a [AtomicUsize],
}

impl Progress<'_> {
    pub fn len(&self) -> usize {
        self.counters.len()
    }

    pub fn is_empty(&self) -> bool {
        self.counters.is_empty()
    }

    pub fn get(&self, task: usize) -> usize {
        self.counters[task].load(Ordering::Relaxed)
    }
}

/// Why [`run`] stopped early.
#[derive(Debug)]
pub enum RunError<E> {
    /// A task failed (the first failure; the others were cancelled).
    Task(PspError),
    /// `tick` asked to stop (for R: a user interrupt).
    Tick(E),
}

/// Runs tasks `0..n_tasks` on `threads` worker threads and returns their
/// results in task order.
///
/// `tick` runs on the calling thread at least every `every`, and after each
/// finished task. If it returns an error, or a task fails, the cancel flag is
/// set, the workers stop at their next check, and the error is returned once
/// all threads have been joined.
pub fn run<T, F, E>(
    n_tasks: usize,
    threads: usize,
    every: Duration,
    task: F,
    mut tick: impl FnMut(&Progress) -> Result<(), E>,
) -> Result<Vec<T>, RunError<E>>
where
    T: Send,
    F: Fn(usize, &TaskCtx) -> Result<T, PspError> + Sync,
{
    let cancel = AtomicBool::new(false);
    let next = AtomicUsize::new(0);
    let counters: Vec<AtomicUsize> = (0..n_tasks).map(|_| AtomicUsize::new(0)).collect();
    let workers = threads.clamp(1, n_tasks.max(1));
    let mut results: Vec<Option<T>> = (0..n_tasks).map(|_| None).collect();

    let outcome = std::thread::scope(|s| {
        let (tx, rx) = mpsc::channel::<(usize, Result<T, PspError>)>();
        for _ in 0..workers {
            let tx = tx.clone();
            let (cancel, next, counters, task) = (&cancel, &next, &counters, &task);
            s.spawn(move || loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= n_tasks || cancel.load(Ordering::Relaxed) {
                    break;
                }
                let ctx = TaskCtx {
                    cancel,
                    progress: &counters[i],
                };
                let res = task(i, &ctx);
                let failed = res.is_err();
                if tx.send((i, res)).is_err() || failed {
                    break;
                }
            });
        }
        drop(tx);

        let progress = Progress {
            counters: &counters,
        };
        let mut outcome: Result<(), RunError<E>> = Ok(());
        loop {
            match rx.recv_timeout(every) {
                Ok((i, Ok(v))) => results[i] = Some(v),
                Ok((_, Err(e))) => {
                    cancel.store(true, Ordering::Relaxed);
                    // Keep the first real failure, not the cancellations it causes.
                    if outcome.is_ok() && e != PspError::Cancelled {
                        outcome = Err(RunError::Task(e));
                    }
                }
                Err(RecvTimeoutError::Timeout) => {}
                Err(RecvTimeoutError::Disconnected) => break,
            }
            if outcome.is_ok() {
                if let Err(e) = tick(&progress) {
                    cancel.store(true, Ordering::Relaxed);
                    outcome = Err(RunError::Tick(e));
                }
            }
        }
        outcome
    });
    outcome?;
    results
        .into_iter()
        .map(|r| r.ok_or(RunError::Task(PspError::Cancelled)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::{RMrg32k3a, Rng};

    fn work(i: usize, ctx: &TaskCtx) -> Result<f64, PspError> {
        let mut rng = RMrg32k3a::from_seed([1, 2, 3, 4, 5, i as i32 + 6])?;
        let mut s = 0.0;
        for k in 0..20_000 {
            if k % 1000 == 0 {
                ctx.check()?;
                ctx.set_progress(k);
            }
            s += rng.norm();
        }
        Ok(s)
    }

    #[test]
    fn chunks_cover_the_range() {
        assert_eq!(chunks(0, 8), Vec::<Range<usize>>::new());
        assert_eq!(chunks(17, 8), vec![0..8, 8..16, 16..17]);
        assert_eq!(chunks(16, 8), vec![0..8, 8..16]);
    }

    #[test]
    fn results_do_not_depend_on_the_thread_count() {
        let tick = |_: &Progress| Ok::<(), ()>(());
        let every = Duration::from_millis(5);
        let one = run(13, 1, every, work, tick).unwrap();
        for threads in [2, 4, 32] {
            let many = run(13, threads, every, work, tick).unwrap();
            assert_eq!(
                one.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                many.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
            );
        }
        assert!(run(0, 4, every, work, tick).unwrap().is_empty());
    }

    #[test]
    fn a_failing_task_stops_the_run() {
        let task = |i: usize, ctx: &TaskCtx| {
            if i == 3 {
                return Err(PspError::InvalidInput("task 3".to_string()));
            }
            work(i, ctx)
        };
        let res = run(50, 4, Duration::from_millis(5), task, |_| Ok::<(), ()>(()));
        match res {
            Err(RunError::Task(PspError::InvalidInput(m))) => assert_eq!(m, "task 3"),
            other => panic!("unexpected: {other:?}"),
        }
    }

    #[test]
    fn tick_can_cancel() {
        let mut ticks = 0;
        let res = run(1000, 2, Duration::from_millis(1), work, |_| {
            ticks += 1;
            if ticks > 3 {
                Err("stop")
            } else {
                Ok(())
            }
        });
        assert!(matches!(res, Err(RunError::Tick("stop"))));
    }
}
