//! # Base Depth Analysis
//!
//! Performs a single pass over a BAM/CRAM file to calculate depth at each position
//! as well as depth per nucleotide. Additionally counts the number of
//! insertions / deletions at each position.
//!
//! This implementation is optimized for detecting indels in reduced representation
//! sequencing data, which is the primary purpose of REDICAT.
//!
//! # Features
//!
//! - Parallel processing using the par_granges framework
//! - Configurable filtering based on base quality and depth thresholds
//! - Support for both 0-based and 1-based coordinate output
//! - Automatic gzip compression for output
//! - Support for BED and BCF/VCF region restrictions
use anyhow::{Context, Result};
use log::*;
use parking_lot::Mutex;
use redicat_lib::{
    par_granges::{self, RegionProcessor},
    position::pileup_position::PileupPosition,
    read_filter::DefaultReadFilter,
    utils,
};
use rust_htslib::{bam, bam::Read};
use std::cmp::max;
use std::sync::Arc;
use std::{convert::TryInto, path::PathBuf};
use structopt::StructOpt;

use crate::commands::is_standard_contig;
use rustc_hash::FxHashSet;

/// Calculate the depth at each base, per-nucleotide.
#[derive(StructOpt)]
#[structopt(author, name = "bulk")]
pub struct Bulk {
    /// Input indexed BAM/CRAM to analyze.
    reads: PathBuf,

    /// Output path (required, will automatically append .gz if not present).
    #[structopt(long, short = "o")]
    output: PathBuf,

    /// The number of threads to use.
    #[structopt(long, short = "t", default_value = "10")]
    threads: usize,

    /// The ideal number of basepairs each worker receives. Total bp in memory at one time is (threads - 2) * chunksize.
    #[structopt(long, short = "c", default_value=par_granges::CHUNKSIZE_STR.as_str())]
    chunksize: u32,

    /// Minium base quality for a base to be counted toward [A, C, T, G]. If the base is less than the specified
    /// quality score it will instead be counted as an `N`. If nothing is set for this to 30.
    #[structopt(long, short = "Q")]
    min_baseq: Option<u8>,

    /// Minimum mapping quality for reads to be counted (default matches legacy preprocess).
    #[structopt(long, default_value = "255", short = "q")]
    mapquality: u8,

    /// Output positions as 0-based instead of 1-based.
    #[structopt(long, short = "z")]
    zero_base: bool,

    /// Set the max depth for a pileup. If a positions depth is within 1% of max-depth the `NEAR_MAX_DEPTH`
    /// output field will be set to true and that position should be viewed as suspect.
    #[structopt(long, short = "D", default_value = "8000")]
    max_depth: u32,
    /// The minimum valid depth to report on. If a position has a depth less than this it will not be reported.
    #[structopt(long, short = "d", default_value = "10")]
    min_depth: u32,
    /// The number of N at a position must be less than or equal to this fraction of the total depth to be reported.
    #[structopt(long, short = "n", default_value = "20")]
    max_n_fraction: u32,
    /// Report all positions, not just edited ones. When set, applies less stringent filtering.
    #[structopt(long, short = "a")]
    all: bool,

    /// Editing threshold for valid value calculation (pos.depth / editing_threshold).
    #[structopt(long = "editing-threshold", short = "et", default_value = "1000")]
    editing_threshold: u32,

    /// When set, traverse every contig present in the BAM header instead of restricting to canonical chroms.
    #[structopt(long = "allcontigs", short = "A")]
    all_contigs: bool,
}

impl Bulk {
    pub fn run(self) -> Result<()> {
        info!("Running redicat-bulk on: {:?}", self.reads);
        let cpus = utils::determine_allowed_cpus(self.threads)?;

        // Ensure output file has .gz extension
        let output_path = if self.output.extension().and_then(|ext| ext.to_str()) == Some("gz") {
            self.output.clone()
        } else {
            let mut new_path = self.output.clone();
            let new_name = format!("{}.gz", self.output.file_name().unwrap().to_string_lossy());
            new_path.set_file_name(new_name);
            new_path
        };

        let mut writer = utils::get_writer(
            &Some(output_path),
            true, // Always use compression
            true,
            1, // Single-threaded compression
            6, // Default compression level
        )?;

        let read_filter = Arc::new(DefaultReadFilter::new(self.mapquality));

        let allowed_tids = (!self.all_contigs)
            .then(|| Self::resolve_standard_tids(&self.reads))
            .transpose()?;

        let base_processor = BaseProcessor::new(
            self.reads.clone(),
            if self.zero_base { 0 } else { 1 },
            self.max_depth,
            self.min_depth,
            self.max_n_fraction,
            self.min_baseq,
            !self.all, // Use edited filtering logic based on command line parameter (inverted: --all means less stringent)
            self.editing_threshold,
            read_filter,
            allowed_tids,
        );

        let par_granges_runner = par_granges::ParGranges::new(
            self.reads.clone(),
            None,
            None,
            None,
            false,
            Some(cpus),
            Some(self.chunksize),
            Some(0.5), // Fixed channel_size_modifier
            base_processor,
        );

        let receiver = par_granges_runner.process()?;

        for pos in receiver.into_iter() {
            writer.serialize(pos)?
        }

        writer.flush()?;
        Ok(())
    }

    /// Collect the numeric IDs of canonical chromosomes to prune work during pileup iteration.
    fn resolve_standard_tids(reads: &PathBuf) -> Result<FxHashSet<u32>> {
        let header = bam::IndexedReader::from_path(reads)
            .with_context(|| format!("Failed to open {}", reads.display()))?
            .header()
            .to_owned();
        let mut allowed = FxHashSet::default();
        for tid in 0..header.target_count() {
            let name = std::str::from_utf8(header.tid2name(tid)).unwrap_or("");
            if is_standard_contig(name) {
                allowed.insert(tid);
            }
        }
        Ok(allowed)
    }
}

/// Holds the info needed for [par_granges::RegionProcessor] implementation
struct BaseProcessor {
    /// path to indexed BAM/CRAM
    reads: PathBuf,
    /// 0-based or 1-based coordiante output
    coord_base: u32,
    /// max depth to pass to htslib pileup engine, max value is MAX(i32)
    max_depth: u32,
    /// min depth to report on
    min_depth: u32,
    /// max percentage of N to report on
    max_n_fraction: u32,
    /// an optional base quality score. If Some(number) if the base quality is not >= that number the base is treated as an `N`
    min_baseq: Option<u8>,
    /// whether to use edited filtering logic
    edited: bool,
    /// editing threshold for valid value calculation
    editing_threshold: u32,
    /// composite read filter (MAPQ-only)
    read_filter: Arc<DefaultReadFilter>,
    /// pool of reusable BAM readers to avoid reopen cost
    reader_pool: Mutex<Vec<bam::IndexedReader>>,
    /// Optional whitelist of target IDs corresponding to canonical contigs
    allowed_tids: Option<FxHashSet<u32>>,
}

impl BaseProcessor {
    /// Create a new BaseProcessor
    fn new(
        reads: PathBuf,
        coord_base: u32,
        max_depth: u32,
        min_depth: u32,
        max_n_fraction: u32,
        min_baseq: Option<u8>,
        edited: bool,
        editing_threshold: u32,
        read_filter: Arc<DefaultReadFilter>,
        allowed_tids: Option<FxHashSet<u32>>,
    ) -> Self {
        let min_baseq = min_baseq.or(Some(30u8));
        Self {
            reads,
            coord_base,
            max_depth,
            min_depth,
            max_n_fraction,
            min_baseq,
            edited,
            editing_threshold,
            read_filter,
            reader_pool: Mutex::new(Vec::new()),
            allowed_tids,
        }
    }
}

/// Implement [par_granges::RegionProcessor] for [BaseProcessor]
impl RegionProcessor for BaseProcessor {
    /// Objects of [PileupPosition] will be returned by each call to [BaseProcessor::process_region]
    type P = PileupPosition;

    /// Process a region by fetching it from a BAM/CRAM, getting a pileup, and then
    /// walking the pileup (checking bounds) to create Position objects according to
    /// the defined filters
    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<PileupPosition> {
        if let Some(allowed) = &self.allowed_tids {
            if !allowed.contains(&tid) {
                return Vec::new();
            }
        }

        // Create a reader
        let mut reader = {
            let mut pool = self.reader_pool.lock();
            pool.pop()
        }
        .unwrap_or_else(|| {
            bam::IndexedReader::from_path(&self.reads)
                .unwrap_or_else(|e| panic!("Failed to open {}: {}", self.reads.display(), e))
        });

        let header = reader.header().to_owned();
        // fetch the region of interest
        reader.fetch((tid, start, stop)).expect("Fetched a region");
        // Walk over pileups
        let mut pileup: bam::pileup::Pileups<'_, bam::IndexedReader> = reader.pileup();
        pileup.set_max_depth(std::cmp::min(
            i32::max_value().try_into().unwrap(),
            self.max_depth,
        ));

        let mut result = Vec::new();
        for p in pileup {
            let pileup = p.expect("Extracted a pileup");
            // Verify that we are within the bounds of the chunk we are iterating on
            if pileup.pos() >= start && pileup.pos() < stop {
                let mut pos = PileupPosition::from_pileup(
                    pileup,
                    &header,
                    self.read_filter.as_ref(),
                    self.min_baseq,
                );
                pos.pos += self.coord_base;

                let should_include = if self.edited {
                    let valid_value = max(pos.depth / self.editing_threshold, 2);
                    let mut count = 0;
                    if pos.a > valid_value {
                        count += 1;
                    }
                    if pos.t > valid_value {
                        count += 1;
                    }
                    if pos.g > valid_value {
                        count += 1;
                    }
                    if pos.c > valid_value {
                        count += 1;
                    }

                    ((pos.depth - pos.n) >= self.min_depth)
                        && (count >= 2)
                        && (pos.n <= max(pos.depth / self.max_n_fraction, 2))
                } else {
                    // Original simple filtering: just check if coverage is above min_depth
                    (pos.depth - pos.n) >= self.min_depth
                };

                if should_include {
                    result.push(pos);
                }
            }
        }
        {
            let mut pool = self.reader_pool.lock();
            pool.push(reader);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use std::sync::Arc;

    #[test]
    fn test_process_region_chr22() {
        // Check if the test file exists
        let test_bam_path = PathBuf::from("test/chr22.bam");
        if !test_bam_path.exists() {
            info!(
                "Test BAM file not found at {:?}, skipping test",
                test_bam_path
            );
            return;
        }

        // Create a BaseProcessor
        let read_filter = Arc::new(DefaultReadFilter::new(255));
        let processor = BaseProcessor::new(
            test_bam_path,
            1,        // 1-based coordinates
            10000,    // max_depth
            10,       // min_depth
            20,       // max_n_fraction
            Some(30), // min_baseq
            false,    // edited
            1000,     // editing_threshold
            read_filter,
            None,
        );

        // Process the region chr22:11252003-11252024 (0-based coordinates would be 11252002-11252024)
        // For 1-based coordinates, we pass 11252003-11252024 directly
        let results = processor.process_region(14, 22901237, 22901238); // tid=0 for chr22, 0-based start/end

        // Print results using log::info! as requested
        info!("Found {} positions", results.len());
        for pos in &results {
            info!(
                "Position: {}, Depth: {}, A: {}, T: {}, C: {}, G: {}, N: {} Ins: {}, Del: {}, RefSkip: {}, Fail: {}",
                pos.pos, pos.depth, pos.a, pos.t, pos.c, pos.g, pos.n, pos.ins, pos.del, pos.ref_skip, pos.fail
            );
        }

        // Just print the results without asserting - useful for debugging
        // In a real test with actual data, you would add appropriate assertions
    }
}
