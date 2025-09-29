//! A set of implementations of `Position` for different use cases
//!
//! This module provides data structures for storing information about genomic positions.
//! The main trait is [`Position`], which defines the interface for position data structures.
//!
//! The primary implementation is [`pileup_position::PileupPosition`], which holds detailed
//! information about a genomic position derived from a pileup, including depth, nucleotide
//! counts, and indel counts.
//!
//! # Position Trait
//!
//! The [`Position`] trait defines the basic interface that all position data structures
//! must implement. It requires:
//!
//! - `Default` implementation for creating empty positions
//! - `Serialize` implementation for output
//! - `new` method for creating positions with a reference sequence and position
//!
//! # Implementations
//!
//! - [`pileup_position::PileupPosition`]: Detailed per-base information from pileup analysis
//! - [`range_positions::RangePosition`]: Position information for genomic ranges

pub mod pileup_position;
pub mod range_positions;

use serde::Serialize;
use smartstring::alias::String;

/// A serializable object meant to hold all information about a position.
///
/// This trait defines the interface for position data structures. All implementations
/// must be serializable and provide a way to create a new position with a reference
/// sequence and position.
pub trait Position: Default + Serialize {
    /// Create a new position with the given reference sequence name and position.
    ///
    /// # Arguments
    ///
    /// * `ref_seq` - The reference sequence name
    /// * `pos` - The position on the reference sequence (typically 1-based)
    ///
    /// # Returns
    ///
    /// A new instance of the position struct
    fn new(ref_seq: String, pos: u32) -> Self;
}
