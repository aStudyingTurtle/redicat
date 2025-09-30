//! BAM to matrix conversion functionality for single-cell data
//!
//! This module provides functionality to process BAM files and convert them
//! to sparse matrices suitable for single-cell analysis. It supports:
//! - Processing genomic positions from TSV files
//! - Parallel processing using the par_granges framework
//! - UMI deduplication
//! - Cell barcode validation
//! - Output in AnnData (H5AD) format
//!
//! # Key Components
//!
//! - [`anndata_output`]: AnnData output functionality with performance optimizations
//! - [`barcode`]: Cell barcode processing functionality
//! - [`processor`]: Core BAM file processing logic
//! - [`region_processor`]: Region processor implementation for par_granges
//! - [`utils`]: Utility functions for bam2mtx processing

pub mod anndata_output;
pub mod barcode;
pub mod processor;

pub use anndata_output::AnnDataConverter;
pub use barcode::BarcodeProcessor;
pub use processor::BamProcessor;
