//! A trait and default implementation of a read filter.
//!
//! This module provides functionality for filtering reads based on various criteria
//! primarily based on mapping quality.
//!
//! The main trait is [`ReadFilter`], which defines the interface for read filtering.
//! The primary implementation is [`DefaultReadFilter`], which provides straightforward
//! filtering based solely on minimum mapping quality.
//!
//! # ReadFilter Trait
//!
//! The [`ReadFilter`] trait defines a single method `filter_read` that takes a
//! BAM record and optional alignment information and returns a boolean indicating
//! whether the read passes the filter (true) or fails (false).
//!
//! # Default Implementation
//!
//! [`DefaultReadFilter`] provides a standard implementation that filters reads based on the
//! minimum mapping quality.
//!
//! # Example
//!
//! ```rust
//! use redicat_lib::read_filter::{ReadFilter, DefaultReadFilter};
//! use rust_htslib::bam::record::Record;
//!
//! // Create a filter that requires reads to meet a minimum mapping quality of 30
//! let filter = DefaultReadFilter::new(30);
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

/// A straightforward read filter based on mapping quality.
///
/// This implementation filters reads solely based on minimum mapping quality.
pub struct DefaultReadFilter {
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
    /// * `min_mapq` - Minimum mapping quality
    ///
    /// # Returns
    ///
    /// A new DefaultReadFilter instance
    pub fn new(min_mapq: u8) -> Self {
        Self { min_mapq }
    }
}

impl ReadFilter for DefaultReadFilter {
    /// Filter reads based on mapping quality.
    ///
    /// A read passes the filter if its mapping quality is greater than or equal to `min_mapq`.
    ///
    /// # Arguments
    ///
    /// * `read` - The BAM record to filter
    /// * `_alignment` - Optional alignment information (unused in this implementation)
    ///
    /// # Returns
    ///
    /// `true` if the read passes the filter, `false` otherwise
    #[inline(always)]
    fn filter_read(&self, read: &Record, _alignment: Option<&Alignment>) -> bool {
        read.mapq() >= self.min_mapq
    }
}
