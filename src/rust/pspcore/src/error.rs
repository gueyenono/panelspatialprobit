use std::fmt;

/// Errors returned to the R binding, which turns them into R errors.
#[derive(Debug, Clone, PartialEq)]
pub enum PspError {
    /// An input failed validation; the message names it.
    InvalidInput(String),
    /// The run was cancelled (an interrupt, or another task failed).
    Cancelled,
}

impl fmt::Display for PspError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PspError::InvalidInput(msg) => write!(f, "{msg}"),
            PspError::Cancelled => write!(f, "cancelled"),
        }
    }
}

impl std::error::Error for PspError {}
