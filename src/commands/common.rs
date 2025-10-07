use anyhow::{anyhow, Context, Result};
use once_cell::sync::OnceCell;
use rayon::ThreadPoolBuilder;
use redicat_lib::utils;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rustc_hash::FxHashSet;
use std::path::{Path, PathBuf};

static GLOBAL_RAYON_THREADS: OnceCell<usize> = OnceCell::new();

/// Ensure an output path ends with a gzip-compatible extension.
///
/// If the provided path doesn't already end with `.gz`, `.gzip`, or `.bgzf`,
/// a `.gz` suffix is appended to the filename while preserving the original
/// parent directory.
pub fn ensure_gz_path(path: &Path) -> PathBuf {
    if utils::is_bgzipped(path) {
        return path.to_path_buf();
    }

    let mut adjusted = path.to_path_buf();
    if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
        adjusted.set_file_name(format!("{}.gz", name));
    } else {
        adjusted.set_extension("gz");
    }
    adjusted
}

/// Configure the global Rayon thread pool exactly once, returning the active
/// worker count. Subsequent calls reuse the first configured pool and emit a
/// warning when the requested thread count differs from the established size.
pub fn configure_global_thread_pool(threads: usize) -> Result<usize> {
    let requested = utils::determine_allowed_cpus(threads)?;

    if let Some(active) = GLOBAL_RAYON_THREADS.get() {
        if *active != requested {
            log::warn!(
                "Rayon global thread pool already initialised with {} threads; ignoring request for {}",
                active,
                requested
            );
        }
        return Ok(*active);
    }

    match ThreadPoolBuilder::new()
        .num_threads(requested)
        .stack_size(2 * 1024 * 1024) // 2 MB stack per thread to reduce memory overhead
        .build_global()
    {
        Ok(_) => {
            GLOBAL_RAYON_THREADS
                .set(requested)
                .map_err(|_| anyhow!("Failed to record global Rayon thread count"))?;
            Ok(requested)
        }
        Err(err) => {
            // The pool was likely initialised elsewhere; fall back to the current size.
            log::debug!("Global Rayon thread pool initialisation skipped: {}", err);
            let fallback = rayon::current_num_threads();
            if fallback != requested {
                log::warn!(
                    "Using existing Rayon pool with {} threads instead of requested {}",
                    fallback,
                    requested
                );
            }
            GLOBAL_RAYON_THREADS.set(fallback).ok();
            Ok(fallback)
        }
    }
}

/// Collect the numerical target identifiers (TIDs) from an indexed BAM/CRAM
/// header whose contig names satisfy the supplied predicate.
///
/// The predicate receives contig names as UTF-8 strings and should return
/// `true` when the contig should be retained.
pub fn collect_tids_with_filter<P, F>(reads: P, mut predicate: F) -> Result<FxHashSet<u32>>
where
    P: AsRef<Path>,
    F: FnMut(&str) -> bool,
{
    let reader = bam::IndexedReader::from_path(reads.as_ref())
        .with_context(|| format!("Failed to open {}", reads.as_ref().display()))?;
    let header = reader.header().to_owned();

    let mut allowed = FxHashSet::default();
    for tid in 0..header.target_count() {
        let name = std::str::from_utf8(header.tid2name(tid))
            .with_context(|| format!("Invalid contig name at TID {}", tid))?;
        if predicate(name) {
            allowed.insert(tid);
        }
    }

    Ok(allowed)
}
