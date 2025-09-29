//! Genomic position primitives used across REDICAT pipelines.
//!
//! This module provides data structures for storing information about genomic positions.
//! The main trait is [`Position`], which defines the interface for position data structures.
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
    fn new(ref_seq: String, pos: u32) -> Self;
}
