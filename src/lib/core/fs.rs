use anyhow::Result;
use std::ffi::OsStr;
use std::fs;
use std::path::Path;

/// Create parent directories for a path when missing.
pub fn make_parent_dirs<P: AsRef<Path>>(path: P) -> Result<()> {
    if let Some(parent) = path.as_ref().parent() {
        fs::create_dir_all(parent)?;
    }
    Ok(())
}

/// Detect whether a path uses a BGZF-compatible extension.
pub fn is_bgzipped<P: AsRef<Path>>(path: P) -> bool {
    matches!(
        path.as_ref().extension().unwrap_or_else(|| OsStr::new("")),
        ext if ext == "gz" || ext == "gzip" || ext == "bgzf"
    )
}
