//! Read filtering primitives used across REDICAT.
//!
//! This module exposes the [`ReadFilter`] trait along with a collection of
//! implementations, including the default mapping-quality based filter.

use rust_htslib::bam::pileup::Alignment;
use rust_htslib::bam::record::Record;

/// A trait for filtering reads based on various criteria.
///
/// Implementors define how reads should be filtered. Implementations should return
/// `true` if the read passes the filter and `false` otherwise.
pub trait ReadFilter {
    /// Filter a read based on various criteria.
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
    /// Create a new [`DefaultReadFilter`] with the specified criteria.
    pub fn new(min_mapq: u8) -> Self {
        Self { min_mapq }
    }
}

impl ReadFilter for DefaultReadFilter {
    /// Filter reads based on mapping quality.
    #[inline(always)]
    fn filter_read(&self, read: &Record, _alignment: Option<&Alignment>) -> bool {
        read.mapq() >= self.min_mapq
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record_with_mapq(mapq: u8) -> Record {
        let mut record = Record::new();
        record.set_mapq(mapq);
        record
    }

    #[test]
    fn rejects_low_quality_reads() {
        let filter = DefaultReadFilter::new(30);
        let record = record_with_mapq(10);
        assert!(!filter.filter_read(&record, None));
    }

    #[test]
    fn accepts_high_quality_reads() {
        let filter = DefaultReadFilter::new(20);
        let record = record_with_mapq(25);
        assert!(filter.filter_read(&record, None));
    }
}
