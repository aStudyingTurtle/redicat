use anyhow::Error;
use std::io;

/// Returns `true` if the error originated from a broken pipe.
#[inline]
pub fn is_broken_pipe(err: &Error) -> bool {
    err.root_cause()
        .downcast_ref::<io::Error>()
        .map(|io_err| io_err.kind() == io::ErrorKind::BrokenPipe)
        .unwrap_or(false)
}
