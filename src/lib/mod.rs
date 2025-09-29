//! REDICAT: RNA Editing Cellular Assessment Toolkit
//!
//! REDICAT is a highly parallelized utility for analyzing RNA editing events in single-cell RNA-seq data.
//! The library is organized into three layers:
//!
//! - [`core`]: Shared concurrency, IO, error, read filtering, and sparse matrix utilities
//! - [`engine`]: High-performance processing primitives such as parallel genomic region schedulers
//!   and position data structures
//! - [`pipeline`]: End-user workflows, including BAM-to-matrix conversion and RNA editing analysis
//!
//! For backwards compatibility, the legacy `call` and `bam2mtx` modules are re-exported from the
//! new `pipeline` namespace.

pub mod core;
pub mod engine;
pub mod pipeline;
pub mod utils;

pub use pipeline::bam2mtx;
pub use pipeline::call;

/// Convenience exports for common helpers used across binaries and downstream crates.
pub mod prelude {
	pub use crate::core::prelude::*;
	pub use crate::engine::{par_granges::RegionProcessor, ParGranges};
	pub use crate::pipeline::call::{AnnDataContainer, EditingType, ReferenceGenome};
}

