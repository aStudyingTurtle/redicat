pub mod bam2mtx;
pub mod base_depth;
pub mod preprocess;
pub mod call;

// Re-export Bulk struct with new name
pub use base_depth::Bulk;
