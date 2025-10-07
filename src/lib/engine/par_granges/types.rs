use lazy_static::lazy_static;
use serde::Serialize;

/// Number of bytes in a gigabyte.
pub const BYTES_IN_A_GIGABYTE: usize = 1024 * 1024 * 1024;

/// A modifier applied to the channel size formula: `(BYTES_IN_A_GIGABYTE * modifier) * threads / size_of(R::P)`.
/// A value of 0.5 corresponds to roughly 1_000_000 `PileupPosition` objects per worker with room to spare.
pub const CHANNEL_SIZE_MODIFIER: f64 = 0.25;

/// Ideal number of basepairs each worker receives. Total bp in memory at one time ≈ `threads * chunksize`.
pub const CHUNKSIZE: u32 = 500_000;

lazy_static! {
    /// [`CHANNEL_SIZE_MODIFIER`] as a string.
    pub static ref CHANNEL_SIZE_MODIFIER_STR: String = CHANNEL_SIZE_MODIFIER.to_string();
    /// [`CHUNKSIZE`] as a string.
    pub static ref CHUNKSIZE_STR: String = CHUNKSIZE.to_string();
}

/// Trait defining how genomic regions are processed.
pub trait RegionProcessor {
    /// The type returned when processing a region.
    type P: 'static + Send + Sync + Serialize;

    /// Process a genomic region defined by `tid`, `start`, and `stop` (0-based, half-open).
    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P>;
}
