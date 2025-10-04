//! BAM file processing for single-cell data

use std::collections::HashMap;
use std::convert::{TryFrom, TryInto};
use std::io;
use std::path::Path;
use std::sync::Arc;

use anyhow::{anyhow, Context, Result};
use noodles::{
    bam,
    core::{self, Position, Region},
    sam,
};
use noodles::bam::Record;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Value;
use noodles::sam::alignment::record::data::field::value::Array;
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
    /// Chromosome name
    pub chrom: String,
    /// 1-based genomic position
    pub pos: u64,
    /// Counts per cell barcode
    pub counts: HashMap<String, StrandBaseCounts>,
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
    let tag: [u8; 2] = tag
        .try_into()
        .map_err(|_| anyhow!("Tag must be 2 bytes long"))?;

    match record.data().get(&tag) {
        Some(Ok(Value::String(value))) => {
            let raw = std::str::from_utf8(value.as_ref())?;
            Ok(clean_tag_value(raw))
        }
        Some(Ok(Value::Array(array))) => match array {
            Array::UInt8(values) => {
                let bytes: Vec<u8> = values.iter().collect::<io::Result<Vec<_>>>()?;
                let raw = std::str::from_utf8(&bytes)?;
                Ok(clean_tag_value(raw))
            }
            _ => Ok(None),
        },
        Some(Ok(_)) => Ok(None),
        Some(Err(_)) | None => Ok(None),
    }
}

/// Extract and normalize a UMI from the requested BAM tag.
pub fn decode_umi(record: &Record, tag: &[u8]) -> Result<Option<String>> {
    let tag: [u8; 2] = tag
        .try_into()
        .map_err(|_| anyhow!("Tag must be 2 bytes long"))?;

    match record.data().get(&tag) {
        Some(Ok(Value::String(value))) => {
            let raw = std::str::from_utf8(value.as_ref())?;
            Ok(clean_tag_value(raw))
        }
        Some(Ok(Value::Array(array))) => match array {
            Array::UInt8(values) => {
                let bytes: Vec<u8> = values.iter().collect::<io::Result<Vec<_>>>()?;
                let raw = std::str::from_utf8(&bytes)?;
                Ok(clean_tag_value(raw))
            }
            _ => Ok(None),
        },
        Some(Ok(_)) => Ok(None),
        Some(Err(_)) | None => Ok(None),
    }
}

/// Retrieve the canonical base at the requested query position.
pub fn decode_base(record: &Record, qpos: usize) -> Result<char> {
    let base = record
        .sequence()
        .get(qpos)
        .ok_or_else(|| anyhow!("Invalid query position"))?;

    Ok(match base.to_ascii_uppercase() {
        b'A' => 'A',
        b'T' => 'T',
        b'G' => 'G',
        b'C' => 'C',
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

fn query_offset_at(record: &Record, target_ref: i64) -> Option<usize> {
    let start = record.alignment_start().transpose().ok().flatten()?;
    let mut ref_pos = usize::from(start) as i64 - 1;
    let mut read_pos = 0usize;

    for op_result in record.cigar().iter() {
        let op = match op_result {
            Ok(op) => op,
            Err(_) => return None,
        };

        let len = op.len() as usize;

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let span_end = ref_pos + len as i64;
                if target_ref >= ref_pos && target_ref < span_end {
                    let offset = (target_ref - ref_pos) as usize;
                    return Some(read_pos + offset);
                }
                ref_pos = span_end;
                read_pos += len;
            }
            Kind::Insertion => {
                read_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                let span_end = ref_pos + len as i64;
                if target_ref >= ref_pos && target_ref < span_end {
                    return None;
                }
                ref_pos = span_end;
            }
            Kind::SoftClip => {
                read_pos += len;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    None
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
    /// Threshold for skipping sites when depth is excessive (handled upstream)
    pub skip_max_depth: u32,
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
            max_depth: u32::MAX,
            skip_max_depth: u32::MAX,
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
        if pos == 0 {
            return Err(anyhow!("Positions must be 1-based"));
        }

        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .with_context(|| format!("Failed to open {}", bam_path.display()))?;
        let header: sam::Header = reader
            .read_header()
            .context("Failed to read BAM header")?;

        let reference_sequences = header.reference_sequences();
        let tid = reference_sequences
            .get_index_of(chrom.as_bytes())
            .ok_or_else(|| anyhow!("Chromosome '{}' not found", chrom))?;

        let position = Position::try_from(pos as usize)
            .map_err(|_| anyhow!("Invalid genomic coordinate {}", pos))?;
        let interval = core::region::Interval::from(position..=position);
        let region = Region::new(chrom.to_string(), interval);

        let mut query = reader
            .query(&header, &region)
            .with_context(|| format!("Failed to query region {}:{}", chrom, pos))?;

        let target_ref = (pos - 1) as i64;
        let mut counts: FxHashMap<String, StrandBaseCounts> = FxHashMap::default();
        let mut umi_consensus: FxHashMap<(String, String), u8> = FxHashMap::default();
        let mut processed = 0u32;

        while let Some(record_result) = query.next() {
            let record = record_result?;

            if record.flags().is_unmapped() {
                continue;
            }

            let read_tid = match record.reference_sequence_id().transpose()? {
                Some(id) => id,
                None => continue,
            };

            if read_tid != tid {
                continue;
            }

            let mapq = record.mapping_quality().map(u8::from).unwrap_or(0);
            if mapq < self.config.min_mapping_quality {
                continue;
            }

            let Some(qpos) = query_offset_at(&record, target_ref) else {
                continue;
            };

            let base_qual = record
                .quality_scores()
                .as_ref()
                .get(qpos)
                .copied()
                .unwrap_or(0);
            if base_qual < self.config.min_base_quality {
                continue;
            }

            let base = decode_base(&record, qpos)?;
            if base == 'N' {
                continue;
            }

            let cell_barcode = match decode_cell_barcode(&record, self.config.cell_barcode_tag.as_bytes())? {
                Some(barcode) => barcode,
                None => continue,
            };

            if !self.barcode_processor.is_valid(&cell_barcode) {
                continue;
            }

            let umi = match decode_umi(&record, self.config.umi_tag.as_bytes())? {
                Some(umi) => umi,
                None => continue,
            };

            let encoded = match encode_call(
                self.config.stranded,
                base,
                record.flags().is_reverse_complemented(),
            ) {
                Some(value) => value,
                None => continue,
            };

            umi_consensus
                .entry((cell_barcode, umi))
                .and_modify(|existing| {
                    if *existing != encoded {
                        *existing = UMI_CONFLICT_CODE;
                    }
                })
                .or_insert(encoded);

            processed = processed.saturating_add(1);
            if processed >= self.config.max_depth {
                break;
            }
        }

        let mut position_counts = HashMap::new();
        for ((cell_barcode, _umi), encoded) in umi_consensus.drain() {
            if encoded == UMI_CONFLICT_CODE {
                continue;
            }

            let entry = counts.entry(cell_barcode).or_default();
            apply_encoded_call(self.config.stranded, encoded, entry);
        }

        for (barcode, strand_counts) in counts.drain() {
            position_counts.insert(barcode, strand_counts);
        }

        Ok(PositionData {
            chrom: chrom.to_string(),
            pos,
            counts: position_counts,
        })
    }
}
