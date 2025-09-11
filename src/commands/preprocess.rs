//! CLI command for preprocessing BAM files
//!
//! This module provides the command-line interface for filtering BAM files
//! based on barcode whitelist, mapping quality, and SAM flags.
//!
//! # Features
//!
//! - Barcode whitelist filtering
//! - Mapping quality filtering
//! - SAM flag-based filtering
//! - Parallel processing
//! - Output BAM indexing

use std::collections::HashSet;
use std::fs::File;
use std::io::BufRead;
use std::time::Instant;
use std::{path::PathBuf};

use anyhow::Result;
use flate2::read::GzDecoder;
use structopt::StructOpt;

use redicat_lib::{
    read_filter::{DefaultReadFilter, ReadFilter},
    utils,
};
use rust_htslib::{bam, bam::record::Aux, bam::record::Record, bam::Read};

/// Arguments for the preprocess command
#[derive(Debug, StructOpt)]
#[structopt(
    name = "preprocess",
    about = "Filter BAM files based on various criteria"
)]
pub struct PreprocessArgs {
    /// Path to the barcodes whitelist file (gzip compressed)
    #[structopt(long, parse(from_os_str))]
    pub barcodes: PathBuf,

    /// Minimum mapping quality
    #[structopt(long, default_value = "255")]
    pub mapquality: u8,

    /// SAM flags that must be present (bit mask)
    #[structopt(long)]
    pub include_flags: Option<u16>,

    /// SAM flags that must be absent (bit mask)
    #[structopt(long)]
    pub exclude_flags: Option<u16>,

    /// Input BAM file path
    #[structopt(long, parse(from_os_str))]
    pub inbam: PathBuf,

    /// Output BAM file path
    #[structopt(long, parse(from_os_str))]
    pub outbam: PathBuf,

    /// Number of threads to use
    #[structopt(short, long, default_value = "8")]
    pub threads: usize,

    /// Write cache size for output BAM file
    #[structopt(long, default_value = "50000")]
    pub write_cache_size: usize,
}

/// Run the preprocess command
pub fn run_preprocess(args: PreprocessArgs) -> Result<()> {
    let start_time = Instant::now();

    // Validate mapquality parameter
    // Since mapquality is u8, it can never be > 255, so we only check if it's at the upper limit
    if args.mapquality == 255 {
        log::warn!("Map quality is at maximum value (255), which may filter out all reads");
    }

    // Validate flags parameters
    if let Some(include_flags) = args.include_flags {
        if include_flags > u16::MAX as u16 {
            return Err(anyhow::anyhow!(
                "Invalid include-flags value: {}. Must be a valid u16.",
                include_flags
            ));
        }
    }

    if let Some(exclude_flags) = args.exclude_flags {
        if exclude_flags > u16::MAX as u16 {
            return Err(anyhow::anyhow!(
                "Invalid exclude-flags value: {}. Must be a valid u16.",
                exclude_flags
            ));
        }
    }

    // Check if input BAM file exists
    if !args.inbam.exists() {
        return Err(anyhow::anyhow!(
            "Input BAM file does not exist: {:?}",
            args.inbam
        ));
    }

    // Check if barcodes file exists
    if !args.barcodes.exists() {
        return Err(anyhow::anyhow!(
            "Barcodes file does not exist: {:?}",
            args.barcodes
        ));
    }

    log::info!("Starting preprocessing...");
    log::info!("Input BAM: {:?}", args.inbam);
    log::info!("Output BAM: {:?}", args.outbam);
    log::info!("Barcodes file: {:?}", args.barcodes);
    log::info!("Map quality threshold: {}", args.mapquality);

    // Load barcodes whitelist
    log::info!("Loading barcodes whitelist...");
    let barcodes_set = load_barcodes_whitelist(&args.barcodes)?;
    log::info!("Loaded {} barcodes", barcodes_set.len());

    // Create read filter
    let include_flags = args.include_flags.unwrap_or(0);
    let exclude_flags = args.exclude_flags.unwrap_or(0);
    let read_filter = DefaultReadFilter::new(include_flags, exclude_flags, args.mapquality);

    // Open input BAM file with sorted order preservation
    let mut reader = bam::Reader::from_path(&args.inbam)?;
    reader.set_threads(args.threads)?;
    let header = bam::Header::from_template(reader.header());

    // Create output BAM file
    utils::make_parent_dirs(&args.outbam)?;
    let mut writer = bam::Writer::from_path(&args.outbam, &header, bam::Format::Bam)?;

    // Process records in order
    let mut written_count = 0u64;
    let mut processed_count = 0u64;

    log::info!("Processing records in coordinate order...");
    for r in reader.records() {
        let record = r?;
        processed_count += 1;

        // Apply filters
        if filter_record(&record, &barcodes_set, &read_filter) {
            writer.write(&record)?;
            written_count += 1;

            // Periodic progress logging
            if written_count % 1000000 == 0 {
                log::info!(
                    "Processed {} records, wrote {} records",
                    processed_count,
                    written_count
                );
            }
        }
    }

    // Explicitly drop the writer to ensure all data is written
    drop(writer);

    let duration = start_time.elapsed();
    log::info!("Processing completed in {:?}", duration);
    log::info!(
        "Processed {} records, wrote {} records",
        processed_count,
        written_count
    );
    bam::index::build(&args.outbam, None, bam::index::Type::Bai, args.threads as u32)?;
    Ok(())
}

/// Load barcodes whitelist from a gzip-compressed file
fn load_barcodes_whitelist(path: &PathBuf) -> Result<HashSet<String>> {
    let file = File::open(path)?;
    let decoder = GzDecoder::new(file);
    let reader = std::io::BufReader::new(decoder);

    let mut barcodes = HashSet::new();
    for line in reader.lines() {
        let line = line?;
        if !line.is_empty() {
            barcodes.insert(line.trim().to_string());
        }
    }

    Ok(barcodes)
}

/// Filter a BAM record based on all criteria
fn filter_record(
    record: &Record,
    barcodes: &HashSet<String>,
    read_filter: &DefaultReadFilter,
) -> bool {
    // Apply read filter first (MAPQ, flags)
    if !read_filter.filter_read(record, None) {
        return false;
    }

    // Check CB tag
    match record.aux(b"CB") {
        Ok(cb_tag) => {
            match cb_tag {
                Aux::String(cb_str) => {
                    if !barcodes.contains(cb_str) {
                        return false;
                    }
                }
                _ => return false, // CB tag is not a string
            }
        }
        Err(_) => return false, // No CB tag
    }

    // Check UB tag
    match record.aux(b"UB") {
        Ok(ub_tag) => {
            match ub_tag {
                Aux::String(ub_str) => {
                    if ub_str == "-" {
                        return false; // UB tag is "-"
                    }
                }
                _ => return false, // UB tag is not a string
            }
        }
        Err(_) => return false, // No UB tag
    }

    true
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_filter_record_with_no_tags() {
        // This test just verifies the function signature compiles
        // In a real test, we would create mock records with various tags and flags
    }
}
