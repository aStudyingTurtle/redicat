//! REDICAT: RNA Editing Cellular Assessment Toolkit
//!
//! REDICAT is a highly parallelized utility for analyzing RNA editing events in single-cell RNA-seq data.
//! The library provides functionality for:
//! 1. Per-base analysis for indel detection in reduced representation sequencing
//! 2. Single-cell matrix generation from BAM files
//! 3. Parallel processing of genomic regions
//! 4. RNA editing detection and analysis in single-cell data
//!
//! # Modules
//!
//! The main modules are:
//! - [`bam2mtx`]: BAM to matrix conversion functionality for single-cell analysis
//! - [`par_granges`]: Parallel processing of genomic regions
//! - [`position`]: Data structures for genomic positions
//! - [`read_filter`]: Filtering of reads based on various criteria
//! - [`utils`]: Utility functions used throughout the library
//! - [`call`]: RNA editing detection and analysis functionality

pub mod bam2mtx;
pub mod call;
pub mod par_granges;
pub mod position;
pub mod read_filter;
pub mod utils;
