//! Backwards-compatible utility re-exports.
//!
//! The new architecture consolidates shared helpers under `crate::core`. This module
//! re-exports the previous `utils::*` API to avoid breaking existing downstream code
//! while steering new development toward `core::prelude`.

pub use crate::core::concurrency::{determine_allowed_cpus, set_rayon_global_pools_size};
pub use crate::core::errors::is_broken_pipe;
pub use crate::core::fs::{is_bgzipped, make_parent_dirs};
pub use crate::core::io::{get_reader, get_writer};
