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
use rust_htslib::bam::{self, pileup::Alignment, Read};

use crate::bam2mtx::{
    barcode::BarcodeProcessor,
    processor::{BamProcessorConfig, PositionData, StrandBaseCounts},
};
use crate::par_granges::RegionProcessor;

/// Region processor for bam2mtx that implements the RegionProcessor trait
pub struct Bam2MtxProcessor {
    /// Path to the BAM file
    bam_path: PathBuf,
    /// Configuration for BAM processing
    config: BamProcessorConfig,
    /// Barcode processor for filtering valid cell barcodes
    barcode_processor: Arc<BarcodeProcessor>,
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

        let mut counts: HashMap<String, StrandBaseCounts> = HashMap::new();
        let mut umi_consensus: HashMap<String, String> = HashMap::new();

        // Process pileup
        for pileup in reader.pileup() {
            let pileup = pileup?;
            if pileup.pos() != pos {
                continue;
            }

            for alignment in pileup.alignments() {
                if !self.should_process_alignment(&alignment) {
                    continue;
                }

                let record = alignment.record();
                let cell_barcode = self.get_cell_barcode(&record)?;

                if !self.barcode_processor.is_valid(&cell_barcode) {
                    continue;
                }

                let umi = self.get_umi(&record)?;
                let base = self.get_base_at_position(&record, alignment.qpos())?;
                let strand = if record.is_reverse() { '-' } else { '+' };

                // UMI consensus logic
                let key = format!("{}:{}", cell_barcode, umi);
                let base_strand = if self.config.stranded {
                    format!("{}{}", base, strand)
                } else {
                    base.to_string()
                };

                if let Some(existing) = umi_consensus.get(&key) {
                    if existing != &base_strand {
                        umi_consensus.insert(key, "-".to_string());
                    }
                } else {
                    umi_consensus.insert(key, base_strand);
                }
            }
        }

        // Aggregate counts by cell barcode
        for (key, base_strand) in umi_consensus {
            if base_strand == "-" {
                continue; // Skip conflicting UMIs
            }

            let parts: Vec<_> = key.split(':').collect();
            let cell_barcode = parts[0].to_string();

            let counts_entry = counts
                .entry(cell_barcode)
                .or_insert_with(StrandBaseCounts::default);

            if self.config.stranded {
                match base_strand.as_str() {
                    "A+" => counts_entry.forward.a += 1,
                    "T+" => counts_entry.forward.t += 1,
                    "G+" => counts_entry.forward.g += 1,
                    "C+" => counts_entry.forward.c += 1,
                    "A-" => counts_entry.reverse.a += 1,
                    "T-" => counts_entry.reverse.t += 1,
                    "G-" => counts_entry.reverse.g += 1,
                    "C-" => counts_entry.reverse.c += 1,
                    _ => {}
                }
            } else {
                match base_strand.as_str() {
                    "A" => counts_entry.forward.a += 1,
                    "T" => counts_entry.forward.t += 1,
                    "G" => counts_entry.forward.g += 1,
                    "C" => counts_entry.forward.c += 1,
                    _ => {}
                }
            }
        }

        Ok(PositionData {
            chrom,
            pos: (pos + 1) as u64, // Convert back to 1-based position
            counts,
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
    fn get_cell_barcode(&self, record: &bam::record::Record) -> Result<String> {
        let tag = self.config.cell_barcode_tag.as_bytes();

        match record.aux(tag) {
            Ok(value) => {
                match value {
                    bam::record::Aux::String(s) => {
                        let barcode = s;
                        // Handle barcodes with suffixes like "-1"
                        Ok(barcode.split('-').next().unwrap_or(barcode).to_string())
                    }
                    bam::record::Aux::ArrayU8(arr) => {
                        let u8_vec: Vec<u8> = arr.iter().collect();
                        let barcode_str = std::str::from_utf8(&u8_vec)?;
                        Ok(barcode_str
                            .split('-')
                            .next()
                            .unwrap_or(barcode_str)
                            .to_string())
                    }
                    _ => Ok("-".to_string()),
                }
            }
            Err(_) => Ok("-".to_string()),
        }
    }

    /// Get UMI from record
    fn get_umi(&self, record: &bam::record::Record) -> Result<String> {
        let tag = self.config.umi_tag.as_bytes();

        match record.aux(tag) {
            Ok(value) => match value {
                bam::record::Aux::String(s) => Ok(s.to_string()),
                bam::record::Aux::ArrayU8(arr) => {
                    let u8_vec: Vec<u8> = arr.iter().collect();
                    Ok(std::str::from_utf8(&u8_vec)?.to_string())
                }
                _ => Ok("-".to_string()),
            },
            Err(_) => Ok("-".to_string()),
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
        let mut reader = match bam::IndexedReader::from_path(&self.bam_path) {
            Ok(reader) => reader,
            Err(e) => {
                log::error!("Failed to create BAM reader: {}", e);
                return vec![];
            }
        };

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

        results
    }
}
