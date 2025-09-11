//! BAM file processing for single-cell data

use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use anyhow::Result;
use rust_htslib::bam::{self, pileup::Alignment, record::Record, Read};

use crate::bam2mtx::barcode::BarcodeProcessor;

/// Base counts for a specific position
#[derive(Debug, Clone, Default, serde::Serialize)]
pub struct BaseCounts {
    /// Count of adenine (A) nucleotides
    pub a: u32,
    /// Count of thymine (T) nucleotides
    pub t: u32,
    /// Count of guanine (G) nucleotides
    pub g: u32,
    /// Count of cytosine (C) nucleotides
    pub c: u32,
}

/// Strand-specific base counts
#[derive(Debug, Clone, Default, serde::Serialize)]
pub struct StrandBaseCounts {
    /// Base counts for the forward strand
    pub forward: BaseCounts,
    /// Base counts for the reverse strand
    pub reverse: BaseCounts,
}

/// Processed data for a specific genomic position
#[derive(Debug, Clone, serde::Serialize)]
pub struct PositionData {
    /// Chromosome name
    pub chrom: String,
    /// 1-based genomic position
    pub pos: u64,
    /// Counts per cell barcode
    pub counts: HashMap<String, StrandBaseCounts>,
}

/// Configuration for BAM processing
#[derive(Debug, Clone)]
pub struct BamProcessorConfig {
    /// Minimum mapping quality for a read to be considered
    pub min_mapping_quality: u8,
    /// Minimum base quality for a base to be counted
    pub min_base_quality: u8,
    /// Whether the data is stranded (true) or unstranded (false)
    pub stranded: bool,
    /// Tag name for UMI (Unique Molecular Identifier)
    pub umi_tag: String,
    /// Tag name for cell barcode
    pub cell_barcode_tag: String,
}

impl Default for BamProcessorConfig {
    fn default() -> Self {
        Self {
            min_mapping_quality: 255,
            min_base_quality: 30,
            stranded: true,
            umi_tag: "UB".to_string(),
            cell_barcode_tag: "CB".to_string(),
        }
    }
}

/// Main processor for BAM files
pub struct BamProcessor {
    /// Configuration for BAM processing
    config: BamProcessorConfig,
    /// Processor for validating cell barcodes
    barcode_processor: Arc<BarcodeProcessor>,
}

impl BamProcessor {
    /// Create a new BamProcessor
    pub fn new(config: BamProcessorConfig, barcode_processor: Arc<BarcodeProcessor>) -> Self {
        Self {
            config,
            barcode_processor,
        }
    }

    /// Process a single genomic position
    pub fn process_position(&self, bam_path: &Path, chrom: &str, pos: u64) -> Result<PositionData> {
        let mut reader = bam::IndexedReader::from_path(bam_path)?;

        // Convert to 0-based position for rust-htslib
        let start_pos = (pos - 1) as u32;
        let end_pos = pos as u32;

        // Get chromosome ID
        let header = reader.header().to_owned();
        let tid = header
            .tid(chrom.as_bytes())
            .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' not found", chrom))?;

        // Fetch the region
        reader.fetch((tid, start_pos, end_pos))?;
        let mut pileups: bam::pileup::Pileups<'_, bam::IndexedReader> = reader.pileup();
        pileups.set_max_depth(u32::max_value());
        let mut counts = HashMap::new();
        let mut umi_consensus = HashMap::new();

        // Process pileup
        for pileup in pileups {
            let pileup = pileup?;
            if pileup.pos() != start_pos {
                continue;
            }

            for read in pileup.alignments() {
                if !self.should_process_read(&read) {
                    continue;
                }

                let record = read.record();
                let cell_barcode = self.get_cell_barcode(&record)?;

                if !self.barcode_processor.is_valid(&cell_barcode) {
                    continue;
                }

                let umi = self.get_umi(&record)?;
                let base = self.get_base_at_position(&record, read.qpos())?;
                let strand = if record.is_reverse() { '-' } else { '+' };

                // UMI consensus logic
                let key = format!("{}: {}", cell_barcode, umi);
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
            let cell_barcode = parts[0].trim().to_string();

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
            chrom: chrom.to_string(),
            pos,
            counts,
        })
    }

    /// Check if a read should be processed
    fn should_process_read(&self, read: &Alignment) -> bool {
        if read.is_del() || read.is_refskip() {
            return false;
        }

        let record = read.record();

        // Check mapping quality
        if record.mapq() < self.config.min_mapping_quality {
            return false;
        }

        // Check base quality
        if let Some(qpos) = read.qpos() {
            if let Some(qual) = record.qual().get(qpos) {
                if *qual < self.config.min_base_quality {
                    return false;
                }
            }
        }

        true
    }

    /// Get cell barcode from record
    fn get_cell_barcode(&self, record: &Record) -> Result<String> {
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
    fn get_umi(&self, record: &Record) -> Result<String> {
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
    fn get_base_at_position(&self, record: &Record, qpos: Option<usize>) -> Result<char> {
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

/// Process multiple positions in parallel
pub fn process_positions_parallel(
    bam_path: &Path,
    positions: &[(String, u64)],
    config: &BamProcessorConfig,
    barcode_processor: Arc<BarcodeProcessor>,
) -> Result<Vec<PositionData>> {
    use rayon::prelude::*;

    let processor = BamProcessor::new(config.clone(), barcode_processor);

    let results: Result<Vec<_>> = positions
        .par_iter()
        .map(|(chrom, pos)| processor.process_position(bam_path, chrom, *pos))
        .collect();

    results
}

/// Optimized position processing to reduce memory allocations
pub fn process_positions_optimized(
    bam_path: &Path,
    positions: &[(String, u64)],
    config: &BamProcessorConfig,
    barcode_processor: Arc<BarcodeProcessor>,
) -> Result<Vec<PositionData>> {
    use rayon::prelude::*;

    use std::sync::Arc;
    let processor = Arc::new(BamProcessor::new(config.clone(), barcode_processor));

    // Process positions by chromosome to minimize BAM file seeking
    let mut positions_by_chrom: HashMap<String, Vec<u64>> = HashMap::new();
    for (chrom, pos) in positions {
        positions_by_chrom
            .entry(chrom.clone())
            .or_insert_with(Vec::new)
            .push(*pos);
    }

    // Sort positions within each chromosome
    for positions in positions_by_chrom.values_mut() {
        positions.sort_unstable();
    }

    let results: Result<Vec<_>> = positions_by_chrom
        .into_iter()
        .par_bridge()
        .flat_map(|(chrom, positions)| {
            let processor = processor.clone();
            positions
                .into_iter()
                .map(move |pos| processor.process_position(bam_path, &chrom, pos))
                .collect::<Vec<_>>()
        })
        .collect();

    results
}
