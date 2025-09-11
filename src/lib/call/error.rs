//! Error types for the REdiCAT library

use thiserror::Error;

#[derive(Error, Debug)]
pub enum RedicatError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("AnnData error: {0}")]
    AnnData(#[from] anyhow::Error),

    #[error("Polars error: {0}")]
    Polars(#[from] polars::error::PolarsError),

    #[error("Sparse matrix error: {0}")]
    SparseMatrix(String),

    #[error("Invalid input: {0}")]
    InvalidInput(String),

    #[error("File not found: {0}")]
    FileNotFound(String),

    #[error("Parse error: {0}")]
    Parse(String),

    #[error("Reference genome error: {0}")]
    ReferenceGenome(String),

    #[error("Data processing error: {0}")]
    DataProcessing(String),

    #[error("Threshold validation error: {field} must be between {min} and {max}, got {value}")]
    ThresholdValidation {
        field: String,
        min: f64,
        max: f64,
        value: f64,
    },

    #[error("Dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch { expected: String, actual: String },

    #[error("Configuration error: {0}")]
    Config(String), // 添加这个变体

    #[error("Empty data: {0}")]
    EmptyData(String),
}

pub type Result<T> = std::result::Result<T, RedicatError>;

impl From<nalgebra_sparse::SparseFormatError> for RedicatError {
    fn from(err: nalgebra_sparse::SparseFormatError) -> Self {
        RedicatError::SparseMatrix(format!("Sparse format error: {:?}", err))
    }
}
