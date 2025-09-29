//! Region processor for bam2mtx functionality
//!
//! This module implements the RegionProcessor trait for processing genomic regions
//! in the context of single-cell BAM to matrix conversion.
//!
//! The [`Bam2MtxProcessor`] struct implements the [`RegionProcessor`] trait to enable
//! parallel processing of genomic regions using the par_granges framework. It processes
//! BAM files and generates count data for single-cell analysis.
//!
//! # Key Features
//!
//! - Parallel processing of genomic regions using par_granges
//! - Cell barcode validation using BarcodeProcessor
//! - UMI deduplication
//! - Base quality filtering
//! - Stranded and unstranded processing modes
//! - Support for custom BAM tags for UMIs and cell barcodes

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::Result;
use parking_lot::Mutex;
use rust_htslib::bam::{self, pileup::Alignment, Read};
use rustc_hash::FxHashMap;
use smartstring::alias::String as CompactString;

use crate::pipeline::bam2mtx::{
    barcode::BarcodeProcessor,
    processor::{self, BamProcessorConfig, PositionData, StrandBaseCounts, UMI_CONFLICT_CODE},
};
use crate::engine::par_granges::RegionProcessor;

/// Region processor for bam2mtx that implements the RegionProcessor trait
pub struct Bam2MtxProcessor {
    /// Path to the BAM file
    bam_path: PathBuf,
    /// Configuration for BAM processing
    config: BamProcessorConfig,
    /// Barcode processor for filtering valid cell barcodes
    barcode_processor: Arc<BarcodeProcessor>,
    /// Pool of reusable BAM readers to amortize IO setup cost
    reader_pool: Mutex<Vec<bam::IndexedReader>>,
}

impl Bam2MtxProcessor {
    /// Create a new Bam2MtxProcessor
    pub fn new(
        bam_path: PathBuf,
        config: BamProcessorConfig,
        barcode_processor: Arc<BarcodeProcessor>,
    ) -> Self {
        Self {
            bam_path,
            config,
            barcode_processor,
            reader_pool: Mutex::new(Vec::new()),
        }
    }

    /// Process a single genomic position
    fn process_position(
        &self,
        reader: &mut bam::IndexedReader,
        tid: u32,
        pos: u32,
        header: &bam::HeaderView,
    ) -> Result<PositionData> {
        // Convert tid to chromosome name
        let chrom = std::str::from_utf8(header.tid2name(tid))?.to_string();

        // Fetch the position (using 0-based coordinate for rust-htslib)
        reader.fetch((tid, pos, pos + 1))?;

        let mut counts: FxHashMap<CompactString, StrandBaseCounts> = FxHashMap::default();
        let mut umi_consensus: FxHashMap<(CompactString, CompactString), u8> = FxHashMap::default();

        let barcode_capacity = self.barcode_processor.len().max(64).min(4096);
        counts.reserve(barcode_capacity);
        umi_consensus.reserve(barcode_capacity * 4);

        // Process pileup
        for pileup in reader.pileup() {
            let pileup = pileup?;
            if pileup.pos() != pos {
                continue;
            }

            let mut processed = 0u32;
            let mut truncated = false;

            for alignment in pileup.alignments() {
                if !self.should_process_alignment(&alignment) {
                    continue;
                }

                processed = processed.saturating_add(1);
                if processed > self.config.max_depth {
                    truncated = true;
                    break;
                }

                let record = alignment.record();
                let cell_barcode = match self.get_cell_barcode(&record)? {
                    Some(barcode) => barcode,
                    None => continue,
                };

                if !self.barcode_processor.is_valid(cell_barcode.as_str()) {
                    continue;
                }

                let umi = match self.get_umi(&record)? {
                    Some(umi) => umi,
                    None => continue,
                };

                let base = self.get_base_at_position(&record, alignment.qpos())?;
                if base == 'N' {
                    continue;
                }
                if let Some(encoded) =
                    processor::encode_call(self.config.stranded, base, record.is_reverse())
                {
                    umi_consensus
                        .entry((cell_barcode, umi))
                        .and_modify(|existing| {
                            if *existing != encoded {
                                *existing = UMI_CONFLICT_CODE;
                            }
                        })
                        .or_insert(encoded);
                }
            }

            if truncated && self.config.skip_max_depth {
                return Ok(PositionData {
                    chrom,
                    pos: (pos + 1) as u64,
                    counts: HashMap::new(),
                });
            }
        }

        // Aggregate counts by cell barcode
        let mut position_counts = HashMap::with_capacity(counts.len());
        for ((cell_barcode, _umi), encoded) in umi_consensus.drain() {
            if encoded == UMI_CONFLICT_CODE {
                continue;
            }

            let counts_entry = counts
                .entry(cell_barcode)
                .or_insert_with(StrandBaseCounts::default);

            processor::apply_encoded_call(self.config.stranded, encoded, counts_entry);
        }

        for (barcode, strand_counts) in counts.drain() {
            position_counts.insert(barcode.into(), strand_counts);
        }

        Ok(PositionData {
            chrom,
            pos: (pos + 1) as u64, // Convert back to 1-based position
            counts: position_counts,
        })
    }

    /// Check if an alignment should be processed
    fn should_process_alignment(&self, alignment: &Alignment) -> bool {
        if alignment.is_del() || alignment.is_refskip() {
            return false;
        }

        let record = alignment.record();

        // Check mapping quality
        if record.mapq() < self.config.min_mapping_quality {
            return false;
        }

        // Check base quality
        if let Some(qpos) = alignment.qpos() {
            if let Some(qual) = record.qual().get(qpos) {
                if *qual < self.config.min_base_quality {
                    return false;
                }
            }
        }

        true
    }

    /// Get cell barcode from record
    fn get_cell_barcode(&self, record: &bam::record::Record) -> Result<Option<CompactString>> {
        let tag = self.config.cell_barcode_tag.as_bytes();

        match record.aux(tag) {
            Ok(value) => {
                match value {
                    bam::record::Aux::String(s) => {
                        let barcode = s;
                        // Handle barcodes with suffixes like "-1"
                        let clean = barcode.split('-').next().unwrap_or(barcode);
                        if clean.is_empty() || clean == "-" {
                            return Ok(None);
                        }
                        Ok(Some(CompactString::from(clean)))
                    }
                    bam::record::Aux::ArrayU8(arr) => {
                        let u8_vec: Vec<u8> = arr.iter().collect();
                        let barcode_str = std::str::from_utf8(&u8_vec)?;
                        let clean = barcode_str.split('-').next().unwrap_or(barcode_str);
                        if clean.is_empty() || clean == "-" {
                            return Ok(None);
                        }
                        Ok(Some(CompactString::from(clean)))
                    }
                    _ => Ok(None),
                }
            }
            Err(_) => Ok(None),
        }
    }

    /// Get UMI from record
    fn get_umi(&self, record: &bam::record::Record) -> Result<Option<CompactString>> {
        let tag = self.config.umi_tag.as_bytes();

        match record.aux(tag) {
            Ok(value) => match value {
                bam::record::Aux::String(s) => {
                    if s.is_empty() || s == "-" {
                        Ok(None)
                    } else {
                        Ok(Some(CompactString::from(s)))
                    }
                }
                bam::record::Aux::ArrayU8(arr) => {
                    let u8_vec: Vec<u8> = arr.iter().collect();
                    let umi_str = std::str::from_utf8(&u8_vec)?;
                    if umi_str.is_empty() || umi_str == "-" {
                        Ok(None)
                    } else {
                        Ok(Some(CompactString::from(umi_str)))
                    }
                }
                _ => Ok(None),
            },
            Err(_) => Ok(None),
        }
    }

    /// Get base at specific position
    fn get_base_at_position(
        &self,
        record: &bam::record::Record,
        qpos: Option<usize>,
    ) -> Result<char> {
        let qpos = qpos.ok_or_else(|| anyhow::anyhow!("Invalid query position"))?;

        let seq = record.seq();
        let base = seq.as_bytes()[qpos];

        match base {
            b'A' | b'a' => Ok('A'),
            b'T' | b't' => Ok('T'),
            b'G' | b'g' => Ok('G'),
            b'C' | b'c' => Ok('C'),
            _ => Ok('N'),
        }
    }
}

impl RegionProcessor for Bam2MtxProcessor {
    /// The output type is a vector of PositionData
    type P = PositionData;

    /// Process a genomic region and return position data for each position
    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P> {
        let mut reader = {
            let mut pool = self.reader_pool.lock();
            pool.pop()
        }
        .unwrap_or_else(|| {
            bam::IndexedReader::from_path(&self.bam_path).unwrap_or_else(|e| {
                panic!(
                    "Failed to open BAM {} for region processing: {}",
                    self.bam_path.display(),
                    e
                )
            })
        });

        let header = reader.header().to_owned();
        let mut results = Vec::new();

        // Process each position in the region
        for pos in start..stop {
            match self.process_position(&mut reader, tid, pos, &header) {
                Ok(position_data) => {
                    // Only include positions that have data
                    if !position_data.counts.is_empty() {
                        results.push(position_data);
                    }
                }
                Err(e) => {
                    log::warn!(
                        "Failed to process position {}:{}",
                        std::str::from_utf8(header.tid2name(tid)).unwrap_or("unknown"),
                        pos + 1
                    ); // 1-based for logging
                    log::debug!("Error: {}", e);
                }
            }
        }

        {
            let mut pool = self.reader_pool.lock();
            pool.push(reader);
        }

        results
    }
}
