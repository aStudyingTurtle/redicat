//! RNA editing detection and analysis module
//!
//! This module provides the core functionality for RNA editing analysis in
//! single-cell RNA-seq data. It includes components for:
//! - AnnData file I/O operations
//! - Base matrix operations
//! - Editing analysis pipelines
//! - Reference genome operations
//! - Sparse matrix operations
//! - Error handling
//!
//! # Key Components
//!
//! - [`anndata_ops`]: Reading and writing AnnData files in H5AD format
//! - [`base_matrix`]: Operations for working with base count matrices
//! - [`editing_analysis`]: Core RNA editing analysis pipeline functions
//! - [`error`]: Error types specific to the call module
//! - [`reference_genome`]: Reference genome operations
//! - [`sparse_ops`]: Sparse matrix operations
//! - [`editing`]: Editing type definitions and related functionality
//! - [`validation`]: Input validation functions

pub mod anndata_ops;
pub mod base_matrix;
pub mod editing;
pub mod editing_analysis;
pub mod error;
pub mod reference_genome;
pub mod sparse_ops;
pub mod validation;
// Re-export main types
pub use anndata_ops::AnnDataContainer;
pub use editing::{load_rediportal_parallel, EditingType};
pub use editing_analysis::{
    annotate_variants_pipeline, calculate_cei, calculate_ref_alt_matrices,
    calculate_site_mismatch_stats,
};
pub use error::{RedicatError, Result};
pub use reference_genome::ReferenceGenome;
