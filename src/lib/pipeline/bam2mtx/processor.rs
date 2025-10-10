//! BAM file processing for single-cell data
use std::path::Path;
use std::sync::Arc;

use anyhow::{anyhow, Result};
use rust_htslib::bam::{self, pileup::Alignment, record::Record, Read};
use rustc_hash::FxHashMap;

use crate::pipeline::bam2mtx::barcode::BarcodeProcessor;

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
    /// Numeric contig identifier (matches the BAM header TID)
    pub contig_id: u32,
    /// 1-based genomic position
    pub pos: u64,
    /// Counts per cell barcode (indexed by whitelist order)
    pub counts: FxHashMap<u32, StrandBaseCounts>,
}

/// Consensus code used to indicate conflicting UMI calls
pub const UMI_CONFLICT_CODE: u8 = u8::MAX;

fn clean_tag_value(raw: &str) -> Option<String> {
    let clean = raw.split('-').next().unwrap_or(raw).trim();
    if clean.is_empty() || clean == "-" {
        None
    } else {
        Some(clean.to_string())
    }
}

/// Extract and normalize a cell barcode from the requested BAM tag.
pub fn decode_cell_barcode(record: &Record, tag: &[u8]) -> Result<Option<String>> {
    match record.aux(tag) {
        Ok(bam::record::Aux::String(s)) => Ok(clean_tag_value(s)),
        Ok(bam::record::Aux::ArrayU8(arr)) => {
            let bytes: Vec<u8> = arr.iter().collect();
            let raw = std::str::from_utf8(&bytes)?;
            Ok(clean_tag_value(raw))
        }
        Ok(_) => Ok(None),
        Err(_) => Ok(None),
    }
}

/// Extract and normalize a UMI from the requested BAM tag.
pub fn decode_umi(record: &Record, tag: &[u8]) -> Result<Option<String>> {
    match record.aux(tag) {
        Ok(bam::record::Aux::String(s)) => Ok(clean_tag_value(s)),
        Ok(bam::record::Aux::ArrayU8(arr)) => {
            let bytes: Vec<u8> = arr.iter().collect();
            let raw = std::str::from_utf8(&bytes)?;
            Ok(clean_tag_value(raw))
        }
        Ok(_) => Ok(None),
        Err(_) => Ok(None),
    }
}

/// Retrieve the canonical base at the requested query position.
pub fn decode_base(record: &Record, qpos: Option<usize>) -> Result<char> {
    let qpos = qpos.ok_or_else(|| anyhow!("Invalid query position"))?;
    let seq = record.seq();
    let base = seq.as_bytes()[qpos];

    Ok(match base {
        b'A' | b'a' => 'A',
        b'T' | b't' => 'T',
        b'G' | b'g' => 'G',
        b'C' | b'c' => 'C',
        _ => 'N',
    })
}

#[inline]
pub fn encode_call(stranded: bool, base: char, is_reverse: bool) -> Option<u8> {
    let base_code = match base {
        'A' => 0,
        'T' => 1,
        'G' => 2,
        'C' => 3,
        _ => return None,
    };

    if stranded {
        let strand_bit = if is_reverse { 1 } else { 0 };
        Some((base_code << 1) | strand_bit)
    } else {
        Some(base_code)
    }
}

#[inline]
pub fn apply_encoded_call(stranded: bool, code: u8, counts_entry: &mut StrandBaseCounts) {
    if stranded {
        let strand_bit = code & 1;
        let base_code = code >> 1;
        let target = if strand_bit == 1 {
            &mut counts_entry.reverse
        } else {
            &mut counts_entry.forward
        };

        match base_code {
            0 => target.a += 1,
            1 => target.t += 1,
            2 => target.g += 1,
            3 => target.c += 1,
            _ => {}
        }
    } else {
        match code {
            0 => counts_entry.forward.a += 1,
            1 => counts_entry.forward.t += 1,
            2 => counts_entry.forward.g += 1,
            3 => counts_entry.forward.c += 1,
            _ => {}
        }
    }
}

/// Configuration for BAM processing
#[derive(Debug, Clone)]
pub struct BamProcessorConfig {
    /// Minimum mapping quality for a read to be considered
    pub min_mapping_quality: u8,
    /// Minimum base quality for a base to be counted
    pub min_base_quality: u8,
    /// Minimum depth (excluding Ns) required to keep a position
    pub min_depth: u32,
    /// Maximum allowed N fraction denominator (depth / max_n_fraction)
    pub max_n_fraction: u32,
    /// Editing threshold used to require multi-base support
    pub editing_threshold: u32,
    /// Whether the data is stranded (true) or unstranded (false)
    pub stranded: bool,
    /// Maximum pileup depth to examine per genomic position
    pub max_depth: u32,
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
            min_depth: 10,
            max_n_fraction: 20,
            editing_threshold: 1000,
            stranded: true,
            max_depth: 65_536,
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
        pileups.set_max_depth(self.config.max_depth.min(i32::MAX as u32));
        let mut counts: FxHashMap<u32, StrandBaseCounts> = FxHashMap::default();
        let mut umi_consensus: FxHashMap<(u32, String), u8> = FxHashMap::default();

        // Process pileup
        for pileup in pileups {
            let pileup = pileup?;
            if pileup.pos() != start_pos {
                continue;
            }

            if (pileup.depth() as u32) >= self.config.max_depth {
                continue;
            }

            // let mut processed = 0u32;

            for read in pileup.alignments() {
                if !self.should_process_read(&read) {
                    continue;
                }

                // processed = processed.saturating_add(1);
                // if processed > self.config.max_depth {
                //     break;
                // }

                let record = read.record();
                let cell_id =
                    match decode_cell_barcode(&record, self.config.cell_barcode_tag.as_bytes())? {
                        Some(barcode) => match self.barcode_processor.id_of(&barcode) {
                            Some(id) => id,
                            None => continue,
                        },
                        None => continue,
                    };

                let umi = match decode_umi(&record, self.config.umi_tag.as_bytes())? {
                    Some(umi) => umi,
                    None => continue,
                };

                let base = decode_base(&record, read.qpos())?;
                if let Some(encoded) = encode_call(self.config.stranded, base, record.is_reverse())
                {
                    umi_consensus
                        .entry((cell_id, umi))
                        .and_modify(|existing| {
                            if *existing != encoded {
                                *existing = UMI_CONFLICT_CODE;
                            }
                        })
                        .or_insert(encoded);
                }
            }
        }

        // Aggregate counts by cell barcode
        for ((cell_id, _umi), encoded) in umi_consensus.drain() {
            if encoded == UMI_CONFLICT_CODE {
                continue;
            }

            let counts_entry = counts.entry(cell_id).or_default();

            apply_encoded_call(self.config.stranded, encoded, counts_entry);
        }

        Ok(PositionData {
            contig_id: tid,
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
}
