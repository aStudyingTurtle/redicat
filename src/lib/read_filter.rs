//! A trait and default implementation of a read filter.
//!
//! This module provides functionality for filtering reads based on various criteria
//! such as SAM flags and mapping quality.
//!
//! The main trait is [`ReadFilter`], which defines the interface for read filtering.
//! The primary implementation is [`DefaultReadFilter`], which provides straightforward
//! filtering based on include/exclude flags and minimum mapping quality.
//!
//! # ReadFilter Trait
//!
//! The [`ReadFilter`] trait defines a single method `filter_read` that takes a
//! BAM record and optional alignment information and returns a boolean indicating
//! whether the read passes the filter (true) or fails (false).
//!
//! # Default Implementation
//!
//! [`DefaultReadFilter`] provides a standard implementation that filters reads based on:
//! - Required SAM flags (all specified flags must be present)
//! - Excluded SAM flags (none of the specified flags may be present)
//! - Minimum mapping quality
//!
//! # Example
//!
//! ```rust
//! use redicat_lib::read_filter::{ReadFilter, DefaultReadFilter};
//! use rust_htslib::bam::record::Record;
//!
//! // Create a filter that requires reads to be properly paired (flag 2)
//! // and excludes reads that are duplicates (flag 1024) or secondary alignments (flag 256)
//! // with a minimum mapping quality of 30
//! let filter = DefaultReadFilter::new(2, 1024 | 256, 30);
//!
//! // Use the filter on a record (assuming `record` is a BAM record)
//! // let passes = filter.filter_read(&record, None);
//! ```

use rust_htslib::bam::pileup::Alignment;
use rust_htslib::bam::record::Record;

/// A trait for filtering reads based on various criteria.
///
/// Implementors of this trait define how reads should be filtered. The `filter_read`
/// method should return `true` if the read passes the filter and `false` if it fails.
pub trait ReadFilter {
    /// Filter a read based on various criteria.
    ///
    /// # Arguments
    ///
    /// * `read` - The BAM record to filter
    /// * `alignment` - Optional alignment information (e.g., from pileup)
    ///
    /// # Returns
    ///
    /// `true` if the read passes the filter, `false` otherwise
    fn filter_read(&self, read: &Record, alignment: Option<&Alignment>) -> bool;
}

/// A straightforward read filter based on SAM flags and mapping quality.
///
/// This implementation filters reads based on required flags, excluded flags,
/// and minimum mapping quality.
pub struct DefaultReadFilter {
    /// SAM flags that must be present for a read to pass filtering.
    ///
    /// All bits set in this value must be present in the read's flags for it to pass.
    include_flags: u16,

    /// SAM flags that must be absent for a read to pass filtering.
    ///
    /// None of the bits set in this value may be present in the read's flags for it to pass.
    exclude_flags: u16,

    /// Minimum mapping quality for a read to pass filtering.
    ///
    /// The read's mapping quality must be greater than or equal to this value to pass.
    min_mapq: u8,
}

impl DefaultReadFilter {
    /// Create a new DefaultReadFilter with the specified criteria.
    ///
    /// # Arguments
    ///
    /// * `include_flags` - SAM flags that must be present
    /// * `exclude_flags` - SAM flags that must be absent
    /// * `min_mapq` - Minimum mapping quality
    ///
    /// # Returns
    ///
    /// A new DefaultReadFilter instance
    pub fn new(include_flags: u16, exclude_flags: u16, min_mapq: u8) -> Self {
        Self {
            include_flags,
            exclude_flags,
            min_mapq,
        }
    }
}

impl ReadFilter for DefaultReadFilter {
    /// Filter reads based on SAM flags and mapping quality.
    ///
    /// A read passes the filter if:
    /// 1. All flags in `include_flags` are present in the read's flags
    /// 2. No flags in `exclude_flags` are present in the read's flags
    /// 3. The read's mapping quality is >= `min_mapq`
    ///
    /// # Arguments
    ///
    /// * `read` - The BAM record to filter
    /// * `_alignment` - Optional alignment information (unused in this implementation)
    ///
    /// # Returns
    ///
    /// `true` if the read passes all filter criteria, `false` otherwise
    #[inline(always)]
    fn filter_read(&self, read: &Record, _alignment: Option<&Alignment>) -> bool {
        let flags = read.flags();
        (flags & self.include_flags) == self.include_flags
            && (flags & self.exclude_flags) == 0
            && read.mapq() >= self.min_mapq
    }
}
