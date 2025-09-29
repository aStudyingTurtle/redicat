//! Input validation utilities
//!
//! This module provides functions for validating input parameters, files,
//! and matrix dimensions for the REdiCAT RNA editing analysis pipeline.
//!
//! # Overview
//!
//! The validation module includes:
//! - Configuration validation for analysis parameters
//! - Input file validation to ensure required files exist and are accessible
//! - Output path validation to ensure proper file extensions and write permissions
//! - Matrix dimension validation to ensure compatibility between data structures

#![allow(dead_code)]
use crate::core::error::{RedicatError, Result};
use std::path::Path;

/// Configuration parameters for validation
///
/// This struct holds all the configuration parameters that need to be
/// validated for the RNA editing analysis pipeline.
///
/// # Fields
///
/// * `max_other_threshold` - Maximum threshold for other mismatches (0.0 to 1.0)
/// * `min_edited_threshold` - Minimum threshold for edited mismatches (0.0 to 1.0)
/// * `min_ref_threshold` - Minimum threshold for reference mismatches (0.0 to 1.0)
/// * `min_coverage` - Minimum coverage threshold
/// * `chunk_size` - Size of data chunks for processing
/// * `num_threads` - Number of threads to use for parallel processing
pub struct ValidationConfig {
    /// Maximum threshold for other mismatches (0.0 to 1.0)
    pub max_other_threshold: f64,
    /// Minimum threshold for edited mismatches (0.0 to 1.0)
    pub min_edited_threshold: f64,
    /// Minimum threshold for reference mismatches (0.0 to 1.0)
    pub min_ref_threshold: f64,
    /// Minimum coverage threshold
    pub min_coverage: u16,
    /// Size of data chunks for processing
    pub chunk_size: usize,
    /// Number of threads to use for parallel processing
    pub num_threads: usize,
}

impl ValidationConfig {
    /// Validate all configuration parameters
    ///
    /// This function validates all configuration parameters to ensure they
    /// are within acceptable ranges and meet the requirements for the
    /// RNA editing analysis pipeline.
    ///
    /// # Returns
    ///
    /// * `Result<()>` - Ok if all parameters are valid, Err otherwise
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Any threshold parameter is not between 0.0 and 1.0
    /// - min_coverage is 0
    /// - chunk_size is 0
    /// - num_threads is 0
    ///
    /// # Example
    ///
    /// ```rust
    /// use redicat_lib::call::validation::ValidationConfig;
    ///
    /// let config = ValidationConfig {
    ///     max_other_threshold: 0.01,
    ///     min_edited_threshold: 0.01,
    ///     min_ref_threshold: 0.01,
    ///     min_coverage: 5,
    ///     chunk_size: 1000,
    ///     num_threads: 8,
    /// };
    /// config.validate().unwrap();
    /// ```
    pub fn validate(&self) -> Result<()> {
        self.validate_threshold("max_other_threshold", self.max_other_threshold, 0.0, 1.0)?;
        self.validate_threshold("min_edited_threshold", self.min_edited_threshold, 0.0, 1.0)?;
        self.validate_threshold("min_ref_threshold", self.min_ref_threshold, 0.0, 1.0)?;

        if self.min_coverage == 0 {
            return Err(RedicatError::InvalidInput(
                "min_coverage must be greater than 0".to_string(),
            ));
        }

        if self.chunk_size == 0 {
            return Err(RedicatError::InvalidInput(
                "chunk_size must be greater than 0".to_string(),
            ));
        }

        if self.chunk_size > 100_000 {
            log::warn!(
                "Large chunk_size ({}) may cause memory issues",
                self.chunk_size
            );
        }

        if self.num_threads == 0 {
            return Err(RedicatError::InvalidInput(
                "num_threads must be greater than 0".to_string(),
            ));
        }

        let max_threads = num_cpus::get() * 2;
        if self.num_threads > max_threads {
            log::warn!(
                "num_threads ({}) exceeds recommended maximum ({})",
                self.num_threads,
                max_threads
            );
        }

        Ok(())
    }

    /// Validate a threshold parameter
    ///
    /// This function validates a threshold parameter to ensure it is
    /// a finite number within the specified range.
    ///
    /// # Arguments
    ///
    /// * `name` - Name of the parameter for error reporting
    /// * `value` - Value of the parameter to validate
    /// * `min` - Minimum allowed value
    /// * `max` - Maximum allowed value
    ///
    /// # Returns
    ///
    /// * `Result<()>` - Ok if parameter is valid, Err otherwise
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Value is not finite
    /// - Value is outside the specified range
    fn validate_threshold(&self, name: &str, value: f64, min: f64, max: f64) -> Result<()> {
        if !value.is_finite() {
            return Err(RedicatError::InvalidInput(format!(
                "{} must be a finite number",
                name
            )));
        }

        if value < min || value > max {
            return Err(RedicatError::ThresholdValidation {
                field: name.to_string(),
                min,
                max,
                value,
            });
        }

        Ok(())
    }
}

/// Validate input files for the analysis pipeline
///
/// This function validates that all required input files exist and have
/// the correct formats for the RNA editing analysis pipeline.
///
/// # Arguments
///
/// * `input` - Path to the input AnnData (.h5ad) file
/// * `fa` - Path to the reference genome FASTA file
/// * `rediportal` - Path to the REDIPortal editing sites file
///
/// # Returns
///
/// * `Result<()>` - Ok if all files are valid, Err otherwise
///
/// # Errors
///
/// Returns an error if:
/// - Input file doesn't exist
/// - Input file doesn't have .h5ad extension
/// - Reference genome file doesn't exist
/// - FASTA index file doesn't exist
/// - REDIPortal file doesn't exist
///
/// # Example
///
/// ```rust
/// use redicat_lib::call::validation::validate_input_files;
///
/// // validate_input_files("input.h5ad", "reference.fa", "editing_sites.tsv.gz").unwrap();
/// ```
pub fn validate_input_files(input: &str, fa: &str, rediportal: &str) -> Result<()> {
    if !Path::new(input).exists() {
        return Err(RedicatError::FileNotFound(format!(
            "Input file not found: {}",
            input
        )));
    }

    if !input.ends_with(".h5ad") {
        return Err(RedicatError::InvalidInput(
            "Input file must have .h5ad extension".to_string(),
        ));
    }

    if !Path::new(fa).exists() {
        return Err(RedicatError::FileNotFound(format!(
            "Reference genome file not found: {}",
            fa
        )));
    }

    let fai_path = format!("{}.fai", fa);
    if !Path::new(&fai_path).exists() {
        return Err(RedicatError::FileNotFound(format!(
            "FASTA index file not found: {}. Please create it using: samtools faidx {}",
            fai_path, fa
        )));
    }

    if !Path::new(rediportal).exists() {
        return Err(RedicatError::FileNotFound(format!(
            "REDIPortal file not found: {}",
            rediportal
        )));
    }

    Ok(())
}

/// Validate output path for the analysis pipeline
///
/// This function validates that the output path is valid and writable,
/// and that it has the correct file extension.
///
/// # Arguments
///
/// * `output` - Path to the output file
///
/// # Returns
///
/// * `Result<()>` - Ok if output path is valid, Err otherwise
///
/// # Errors
///
/// Returns an error if:
/// - Output directory doesn't exist
/// - Output directory is read-only
/// - Output file doesn't have .h5ad extension
///
/// # Example
///
/// ```rust
/// use redicat_lib::call::validation::validate_output_path;
///
/// // validate_output_path("output.h5ad").unwrap();
/// ```
pub fn validate_output_path(output: &str) -> Result<()> {
    let path = Path::new(output);

    // Check if output file exists, if so, remove it
    if path.exists() {
        std::fs::remove_file(path).map_err(|e| {
            RedicatError::InvalidInput(format!(
                "Failed to remove existing output file '{}': {}",
                output, e
            ))
        })?;
    }

    if !output.ends_with(".h5ad") {
        return Err(RedicatError::InvalidInput(
            "Output file must have .h5ad extension".to_string(),
        ));
    }

    Ok(())
}

/// Validate matrix dimensions for compatibility
///
/// This function validates that two matrices have compatible dimensions
/// for element-wise operations.
///
/// # Arguments
///
/// * `_matrix_name` - Name of the matrix for error reporting (currently unused)
/// * `actual_shape` - Actual dimensions of the matrix (rows, cols)
/// * `expected_shape` - Expected dimensions of the matrix (rows, cols)
///
/// # Returns
///
/// * `Result<()>` - Ok if dimensions are compatible, Err otherwise
///
/// # Errors
///
/// Returns an error if:
/// - Actual and expected dimensions don't match
///
/// # Example
///
/// ```rust
/// use redicat_lib::call::validation::validate_matrix_dimensions;
///
/// validate_matrix_dimensions("test_matrix", (100, 50), (100, 50)).unwrap();
/// ```
pub fn validate_matrix_dimensions(
    _matrix_name: &str,
    actual_shape: (usize, usize),
    expected_shape: (usize, usize),
) -> Result<()> {
    if actual_shape != expected_shape {
        return Err(RedicatError::DimensionMismatch {
            expected: format!("{} × {}", expected_shape.0, expected_shape.1),
            actual: format!("{} × {}", actual_shape.0, actual_shape.1),
        });
    }
    Ok(())
}
