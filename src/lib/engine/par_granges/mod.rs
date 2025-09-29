//! Parallel genomic range processing utilities.
//!
//! The [`ParGranges`] executor fans out genomic regions across a Rayon pool and
//! streams each region's results through a bounded crossbeam channel. Callers
//! implement [`RegionProcessor`] to define per-region work while benefitting from
//! shared scheduling, interval merging, and CRAM/BAM logistics.

mod intervals;
mod scheduler;
mod types;

pub use scheduler::ParGranges;
pub use types::{
    RegionProcessor, CHANNEL_SIZE_MODIFIER, CHANNEL_SIZE_MODIFIER_STR, CHUNKSIZE, CHUNKSIZE_STR,
};
