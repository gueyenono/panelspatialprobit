//! User interrupts (Ctrl-C, Esc) from the main thread (plan/06 §6.7).
//!
//! `R_CheckUserInterrupt()` longjmps when an interrupt is pending. Run inside
//! `savvy::unwind_protect`, the jump stops at R_UnwindProtect, and the result
//! is `Error::Aborted`: the caller then cancels and joins its threads and
//! returns the error, and savvy's C wrapper resumes R's unwind only after all
//! Rust frames have been dropped. Only the main thread calls this; worker
//! threads never touch R.

extern "C" {
    fn R_CheckUserInterrupt();
}

/// `Err(savvy::Error::Aborted)` if the user interrupted.
pub fn check() -> savvy::Result<()> {
    // SAFETY: called on R's main thread only; the closure calls a C function
    // of R's API that takes no arguments, and any longjmp it makes is caught
    // by R_UnwindProtect inside unwind_protect.
    unsafe {
        savvy::unwind_protect(|| {
            R_CheckUserInterrupt();
            std::ptr::null_mut()
        })
    }
    .map(|_| ())
}
