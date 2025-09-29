pub mod concurrency;
pub mod error;
pub mod errors;
pub mod fs;
pub mod io;
pub mod read_filter;
pub mod sparse;

pub mod prelude {
    pub use super::concurrency::{determine_allowed_cpus, set_rayon_global_pools_size};
    pub use super::error::{RedicatError, Result};
    pub use super::errors::is_broken_pipe;
    pub use super::fs::{is_bgzipped, make_parent_dirs};
    pub use super::io::{get_reader, get_writer};
    pub use super::read_filter::{DefaultReadFilter, ReadFilter};
    pub use super::sparse::{SparseMatrixExt, SparseOps};
}
