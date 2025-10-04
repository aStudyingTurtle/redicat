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
//! Typically, [`PileupPosition`] instances are constructed by aggregating
//! coverage across [`noodles::bam`] records when traversing a genomic region.

use crate::engine::position::Position;
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
        /// Convenience constructor mirroring the [`Position`] trait implementation.
        pub fn create(ref_seq: SmartString<LazyCompact>, pos: u32) -> Self {
            <Self as Position>::new(ref_seq, pos)
        }

    /// Flag the position when observed depth is within 1% of the configured max depth.
    pub fn mark_near_max_depth(&mut self, max_depth: u32) {
        if max_depth == 0 {
            self.near_max_depth = false;
            return;
        }

        let max_depth_u32 = u32::from(max_depth);
        let tolerance = (max_depth_u32 + 99) / 100;

        self.near_max_depth =
            max_depth_u32.saturating_sub(self.depth) <= tolerance;
    }

}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn mark_near_max_depth_accounts_for_failures() {
        let mut position = PileupPosition {
            depth: 98_500,
            fail: 900,
            ..Default::default()
        };

        position.mark_near_max_depth(100_000);
        assert!(!position.near_max_depth);

        position.depth = 80_000;
        position.fail = 0;
        position.near_max_depth = true;
        position.mark_near_max_depth(100_000);
        assert!(!position.near_max_depth);

        position.depth = 10;
        position.fail = 0;
        position.mark_near_max_depth(0);
        assert!(!position.near_max_depth);
    }
}
