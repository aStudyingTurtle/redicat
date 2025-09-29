//! An implementation of `Position` for dealing with pileups.
//!
//! This module provides the [`PileupPosition`] struct, which holds detailed information
//! about a genomic position derived from a pileup. This includes:
//!
//! - Reference sequence name
//! - Position
//! - Reference base
//! - Depth information
//! - Nucleotide counts (A, C, G, T, N)
//! - Indel counts (insertions, deletions)
//! - Reference skip counts
//! - Failed read counts
//! - Near max depth flag
//!
//! The struct provides methods for creating instances from pileup data and updating counts.
//!
//! # PileupPosition Structure
//!
//! The [`PileupPosition`] struct contains all the information needed to represent
//! a genomic position in a pileup analysis. It's designed to be serializable for
//! easy output and implements the [`Position`] trait for compatibility with the
//! par_granges framework.
//!
//! # Usage
//!
//! Typically, [`PileupPosition`] instances are created from htslib pileup objects:
//!
//! ```rust
//! // This is typically done within a RegionProcessor implementation
//! // let position = PileupPosition::from_pileup(pileup, header, read_filter, base_filter);
//! ```

use crate::core::read_filter::ReadFilter;
use crate::engine::position::Position;
use rust_htslib::bam::{
    self,
    pileup::{Alignment, Pileup},
    record::Record,
    HeaderView,
};
use serde::Serialize;
use smartstring::{alias::String, LazyCompact, SmartString};
use std::default;

/// Hold all information about a position.
// NB: The max depth that htslib will return is i32::MAX, and the type of pos for htlib is u32
// There is no reason to go bigger, for now at least
#[derive(Debug, Serialize, Default)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub struct PileupPosition {
    /// Reference sequence name.
    #[serde(rename = "CHR")]
    pub ref_seq: SmartString<LazyCompact>,
    /// 1-based position in the sequence.
    pub pos: u32,
    /// The reference base at this position.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ref_base: Option<char>,
    /// Total depth at this position.
    pub depth: u32,
    /// Number of A bases at this position.
    pub a: u32,
    /// Number of C bases at this position.
    pub c: u32,
    /// Number of G bases at this position.
    pub g: u32,
    /// Number of T bases at this position.
    pub t: u32,
    /// Number of N bases at this position. Any unrecognized base will be counted as an N.
    pub n: u32,
    /// Number of insertions that start to the right of this position.
    /// Does not count toward depth.
    pub ins: u32,
    /// Number of deletions at this position.
    pub del: u32,
    /// Number of refskips at this position. Does not count toward depth.
    pub ref_skip: u32,
    /// Number of reads failing filters at this position.
    pub fail: u32,
    /// Depth is within 1% of max_depth
    pub near_max_depth: bool,
}

impl Position for PileupPosition {
    /// Create a new position for the given ref_seq name.
    fn new(ref_seq: String, pos: u32) -> Self {
        PileupPosition {
            ref_seq: SmartString::from(ref_seq.as_str()),
            pos,
            ..default::Default::default()
        }
    }
}

impl PileupPosition {
    /// Given a record, update the counts at this position
    #[inline(always)]
    fn update<F: ReadFilter>(
        &mut self,
        alignment: &Alignment,
        record: &Record,
        read_filter: &F,
        base_filter: Option<u8>,
    ) {
        // 1. 过滤不满足条件的 read
        if !read_filter.filter_read(record, Some(alignment)) {
            self.depth -= 1;
            self.fail += 1;
            return;
        }

        // 2. 优先处理 refskip 和 deletion
        // 注意：这里的顺序很重要，is_refskip() 为 true 时 is_del() 也可能为 true
        if alignment.is_refskip() {
            self.ref_skip += 1;
            self.depth -= 1; // Refskip 不计入 depth
        } else if alignment.is_del() {
            self.del += 1;
        } else {
            // 3. 处理匹配的碱基
            // .expect() 在调试时能提供更清晰的错误信息，发布版中与 unwrap() 无性能差异
            let qpos = alignment
                .qpos()
                .expect("Pileup alignment should have a query position");

            // 使用 map_or 优雅地处理 Option，判断碱基质量是否过低
            let is_low_qual = base_filter.map_or(false, |cutoff| record.qual()[qpos] < cutoff);

            if is_low_qual {
                self.n += 1;
            } else {
                // 直接对 u8 (字节) 进行匹配，避免了 u8 -> char 的转换
                match record.seq()[qpos].to_ascii_uppercase() {
                    b'A' => self.a += 1,
                    b'C' => self.c += 1,
                    b'G' => self.g += 1,
                    b'T' => self.t += 1,
                    _ => self.n += 1,
                }
            }

            // 4. 检查在此位置之后开始的 insertion
            if let bam::pileup::Indel::Ins(_len) = alignment.indel() {
                self.ins += 1;
            }
        }
    }

    /// Convert a pileup into a `Position`.
    ///
    /// This will walk over each of the alignments and count the number each nucleotide it finds.
    /// It will also count the number of Ins/Dels/Skips that are at each position.
    ///
    /// # Arguments
    ///
    /// * `pileup` - a pileup at a genomic position
    /// * `header` - a headerview for the bam file being read, to get the sequence name
    /// * `read_filter` - a function to filter out reads, returning false will cause a read to be filtered
    /// * `base_filter` - an optional base quality score. If Some(number) if the base quality is not >= that number the base is treated as an `N`
    #[inline]
    pub fn from_pileup<F: ReadFilter>(
        pileup: Pileup,
        header: &bam::HeaderView,
        read_filter: &F,
        base_filter: Option<u8>,
    ) -> Self {
        let name = Self::compact_refseq(header, pileup.tid());
        // make output 1-based
        let mut pos = Self::new(name, pileup.pos());
        pos.depth = pileup.depth();

        for alignment in pileup.alignments() {
            let record = alignment.record();
            Self::update(&mut pos, &alignment, &record, read_filter, base_filter);
        }
        pos
    }

    /// Convert a tid to a [`SmartString<LazyCompact>`].
    #[inline]
    pub fn compact_refseq(header: &HeaderView, tid: u32) -> SmartString<LazyCompact> {
        let name = std::str::from_utf8(header.tid2name(tid)).unwrap();
        String::from(name)
    }
}

#[cfg(test)]
mod tests {

    use crate::core::read_filter::DefaultReadFilter;

    use super::*;
    use rust_htslib::bam::{IndexedReader, Read};
    use std::path::Path;
    #[test]
    fn test_pileup_position_from_bam() {
        let bam_path = Path::new("test/chr22.bam");
        // test for chr22:50783283
        let site_pos = 50783283;
        // Only run this test if  the file exists
        if bam_path.exists() {
            // Open the BAM file
            let mut bam = IndexedReader::from_path(bam_path).expect("Failed to open BAM file");

            // Set the reference (not needed for this specific BAM)

            bam.fetch(("chr22", site_pos - 1, site_pos))
                .expect("Failed to fetch region");

            // Get the header
            let header = bam.header().to_owned();
            let read_filter = DefaultReadFilter::new(255);

            // Process the pileup
            for pileup_result in bam.pileup() {
                let pileup = pileup_result.expect("Failed to get pileup");
                // Convert to PileupPosition
                let position = PileupPosition::from_pileup(
                    pileup,
                    &header,
                    &read_filter,
                    Some(30), // No base quality filtering
                );
                if !(position.pos == site_pos - 1) {
                    continue;
                }
                // Print the position information
                eprintln!("Reference: {}", position.ref_seq);
                eprintln!("Position: {}", position.pos + 1);
                eprintln!("Reference base: {:?}", position.ref_base);
                eprintln!("Depth: {}", position.depth);
                eprintln!(
                    "A: {}, C: {}, G: {}, T: {}, N: {} Insertions: {}, Deletions: {}, RefSkips: {}, Fail: {}",
                    position.a, position.c, position.g, position.t, position.n, position.ins, position.del, position.ref_skip, position.fail
                );
                assert_eq!(position.c, 2208);
                assert_eq!(position.t, 3);
                eprintln!("Near max depth: {}", position.near_max_depth);

                // We only want the first position
                // break;
            }
        } else {
            eprintln!("Test BAM file not found, skipping test");
            // Panic to see if the test is running
            panic!("Test BAM file not found");
        }
    }
}
