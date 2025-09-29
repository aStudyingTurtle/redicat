//! CLI command for bam2mtx functionality with streaming aggregation
//!
//! This module provides the command-line interface for converting BAM files to
//! single-cell count matrices. It supports various optimizations for handling
//! large datasets efficiently including parallel processing, per-site
//! filtering, and configurable chunk sizes.
//!
//! # Features
//!
//! - Parallel processing of genomic positions
//! - UMI deduplication
//! - Cell barcode validation
//! - Support for stranded and unstranded data
//! - Streaming aggregation to keep memory usage bounded
//! - Configurable chunk sizes and matrix density estimation

use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use crate::commands::{base_depth::Bulk as BulkCommand, is_standard_contig};
use anyhow::{anyhow, Result};
use parking_lot::Mutex;
use rayon::prelude::*;
use rust_htslib::bam::Read;
use structopt::StructOpt;

use redicat_lib::bam2mtx::{
    anndata_output::{AnnDataConfig, AnnDataConverter},
    barcode::BarcodeProcessor,
    processor::{
        apply_encoded_call, encode_call, BamProcessorConfig, PositionData, StrandBaseCounts,
        UMI_CONFLICT_CODE,
    },
};
use redicat_lib::utils;
use rustc_hash::FxHashMap;

/// Arguments for the bam2mtx command
#[derive(Debug, StructOpt)]
#[structopt(name = "bam2mtx", about = "Convert BAM files to single-cell matrices")]
pub struct Bam2MtxArgs {
    /// Path to the input BAM file (must be indexed)
    #[structopt(short, long, parse(from_os_str))]
    pub bam: PathBuf,

    /// Optional TSV file with genomic positions (CHR/POS). When omitted in two-pass mode, a list is generated automatically.
    #[structopt(long, parse(from_os_str))]
    pub tsv: Option<PathBuf>,

    /// Path to the cell barcode whitelist file.
    #[structopt(long, parse(from_os_str))]
    pub barcodes: PathBuf,

    /// Enable the two-pass workflow: run `bulk` to build a site list before matrix generation.
    #[structopt(long)]
    pub two_pass: bool,

    /// Output path for the H5AD file
    #[structopt(short, long, parse(from_os_str))]
    pub output: PathBuf,

    /// Number of threads to use (default: 10)
    #[structopt(short, long, default_value = "10")]
    pub threads: usize,

    /// Minimum mapping quality
    #[structopt(long, default_value = "255", short = "q")]
    pub min_mapq: u8,

    /// Minimum base quality
    #[structopt(long, default_value = "30", short = "Q")]
    pub min_baseq: u8,

    /// Minimum effective depth (excluding Ns) required to keep a site
    #[structopt(long = "min-depth", default_value = "10", short = "d")]
    pub min_depth: u32,

    /// Maximum allowed fraction of ambiguous bases (denominator form, like bulk)
    #[structopt(long = "max-n-fraction", default_value = "20", short = "n")]
    pub max_n_fraction: u32,

    /// Editing threshold used to detect multi-base support
    #[structopt(long = "editing-threshold", default_value = "1000", short = "et")]
    pub editing_threshold: u32,

    /// Whether data is stranded (default: unstranded)
    #[structopt(long, short = "S")]
    pub stranded: bool,

    /// Maximum pileup depth to examine per site (reads beyond this limit are ignored)
    #[structopt(long = "max-depth", default_value = "20000", short = "D")]
    pub max_depth: u32,

    /// Skip sites whose observed depth exceeds the configured max-depth (alias: -sd)
    #[structopt(long = "skip-max-depth", short = "s", visible_alias = "sd")]
    pub skip_max_depth: bool,

    /// UMI tag name
    #[structopt(long, default_value = "UB")]
    pub umi_tag: String,

    /// Cell barcode tag name
    #[structopt(long, default_value = "CB")]
    pub cb_tag: String,

    /// Path to reference FASTA file (required for CRAM files)
    #[structopt(long, parse(from_os_str), short = "r")]
    pub reference: Option<PathBuf>,

    /// Chunk size for parallel processing (number of positions per chunk)
    #[structopt(long, default_value = "3000", short = "c")]
    pub chunksize: u32,

    /// Matrix density estimate used to pre-size sparse buffers.
    /// Typical values: sparse matrix 0.001-0.01, medium density 0.01-0.1, dense matrix 0.1+
    #[structopt(long, default_value = "0.005")]
    pub matrix_density: f64,

    /// Include all contigs from the BAM header instead of restricting to canonical chromosomes.
    #[structopt(long = "allcontigs", short = "A")]
    pub all_contigs: bool,
}

/// Represents a genomic position from TSV file
#[derive(Debug, Clone)]
pub struct GenomicPosition {
    pub chrom: String,
    pub pos: u64,
}

/// Chunk of genomic positions for parallel processing
#[derive(Debug, Clone)] // 添加 Clone trait
pub struct PositionChunk {
    pub positions: Vec<GenomicPosition>,
}

/// Read TSV file with streaming parser for better memory efficiency
fn read_tsv_file<P: AsRef<Path>>(tsv_path: P) -> Result<Vec<GenomicPosition>> {
    use flate2::read::GzDecoder;
    use std::io::{BufRead, BufReader};

    let path = tsv_path.as_ref();
    let file = File::open(path)?;
    let reader: Box<dyn BufRead> = if path.extension().map_or(false, |ext| ext == "gz") {
        let decoder = GzDecoder::new(file);
        Box::new(BufReader::with_capacity(256 * 1024, decoder)) // Larger buffer
    } else {
        Box::new(BufReader::with_capacity(256 * 1024, file))
    };

    let mut positions = Vec::new();
    let mut line_buffer = String::with_capacity(128); // Reuse buffer

    for line_result in reader.lines() {
        line_buffer.clear();
        line_buffer = line_result?;
        let line_trimmed = line_buffer.trim_end();

        // Split by tab, use only the first two columns
        let mut cols = line_trimmed.split('\t');
        if let (Some(chrom), Some(pos_str)) = (cols.next(), cols.next()) {
            if let Ok(pos) = pos_str.parse::<u64>() {
                positions.push(GenomicPosition {
                    chrom: chrom.to_string(),
                    pos,
                });
            }
        }
    }

    // Use unstable sort for better performance
    positions.par_sort_unstable_by(|a, b| {
        let chrom_cmp = a.chrom.cmp(&b.chrom);
        if chrom_cmp == std::cmp::Ordering::Equal {
            a.pos.cmp(&b.pos)
        } else {
            chrom_cmp
        }
    });

    positions.shrink_to_fit();
    Ok(positions)
}

/// Split positions into chunks with improved locality
fn chunk_positions(positions: Vec<GenomicPosition>, chunk_size: usize) -> Vec<PositionChunk> {
    if positions.is_empty() {
        return Vec::new();
    }

    let mut chunks = Vec::with_capacity((positions.len() + chunk_size - 1) / chunk_size);
    let mut current_chunk = Vec::with_capacity(chunk_size);
    let mut current_chrom: Option<String> = None;

    for pos in positions {
        if current_chrom
            .as_ref()
            .map_or(false, |chr| chr != &pos.chrom)
            && !current_chunk.is_empty()
        {
            chunks.push(PositionChunk {
                positions: std::mem::replace(&mut current_chunk, Vec::with_capacity(chunk_size)),
            });
            current_chrom = None;
        }

        if current_chunk.is_empty() {
            current_chrom = Some(pos.chrom.clone());
        }

        current_chunk.push(pos);

        if current_chunk.len() >= chunk_size {
            chunks.push(PositionChunk {
                positions: std::mem::replace(&mut current_chunk, Vec::with_capacity(chunk_size)),
            });
            current_chrom = None;
        }
    }

    if !current_chunk.is_empty() {
        chunks.push(PositionChunk {
            positions: current_chunk,
        });
    }

    chunks.shrink_to_fit();
    chunks
}

// Enhanced processor with memory pooling
struct OptimizedChunkProcessor {
    bam_path: PathBuf,
    config: BamProcessorConfig,
    barcode_processor: Arc<BarcodeProcessor>,
    tid_lookup: FxHashMap<String, u32>,
    count_capacity_hint: usize,
    umi_capacity_hint: usize,
    reader_pool: Mutex<Vec<rust_htslib::bam::IndexedReader>>,
}

impl OptimizedChunkProcessor {
    pub fn new(
        bam_path: PathBuf,
        config: BamProcessorConfig,
        barcode_processor: Arc<BarcodeProcessor>,
    ) -> Result<Self> {
        let header = rust_htslib::bam::IndexedReader::from_path(&bam_path)?
            .header()
            .to_owned();

        let mut tid_lookup =
            FxHashMap::with_capacity_and_hasher(header.target_count() as usize, Default::default());

        for tid in 0..header.target_count() {
            if let Ok(name) = std::str::from_utf8(header.tid2name(tid)) {
                tid_lookup.insert(name.to_string(), tid);
            }
        }

        let barcode_count = barcode_processor.len();
        let count_capacity_hint = barcode_count.clamp(64, 4096);
        let umi_capacity_hint = count_capacity_hint.saturating_mul(8);

        Ok(Self {
            bam_path,
            config,
            barcode_processor,
            tid_lookup,
            count_capacity_hint,
            umi_capacity_hint,
            reader_pool: Mutex::new(Vec::new()),
        })
    }

    pub fn process_chunk(&self, chunk: &PositionChunk) -> Result<Vec<PositionData>> {
        if chunk.positions.is_empty() {
            return Ok(Vec::new());
        }

        let mut reader = {
            let mut pool = self.reader_pool.lock();
            pool.pop()
        }
        .unwrap_or_else(|| {
            rust_htslib::bam::IndexedReader::from_path(&self.bam_path)
                .unwrap_or_else(|e| panic!("Failed to open BAM {}: {}", self.bam_path.display(), e))
        });
        debug_assert!(chunk
            .positions
            .iter()
            .all(|p| p.chrom == chunk.positions[0].chrom));

        let chrom = &chunk.positions[0].chrom;
        let tid = match self.tid_lookup.get(chrom) {
            Some(tid) => *tid,
            None => {
                let mut pool = self.reader_pool.lock();
                pool.push(reader);
                return Ok(Vec::new());
            }
        };

        let fetch_start = chunk
            .positions
            .first()
            .map(|p| p.pos.saturating_sub(1) as u32)
            .unwrap_or(0);
        let fetch_end = chunk
            .positions
            .last()
            .map(|p| p.pos as u32)
            .unwrap_or(fetch_start);

        reader.fetch((tid, fetch_start, fetch_end))?;

        let target_positions: Vec<u32> = chunk
            .positions
            .iter()
            .map(|p| p.pos.saturating_sub(1) as u32)
            .collect();

        let mut chunk_results = Vec::with_capacity(chunk.positions.len());
        let mut position_index = 0usize;

        let mut counts: FxHashMap<String, StrandBaseCounts> = FxHashMap::default();
        counts.reserve(self.count_capacity_hint);
        let mut umi_consensus: FxHashMap<(String, String), u8> = FxHashMap::default();
        umi_consensus.reserve(self.umi_capacity_hint);

        for pileup in reader.pileup() {
            let pileup = pileup?;
            let pile_pos = pileup.pos();

            while position_index < target_positions.len()
                && target_positions[position_index] < pile_pos
            {
                position_index += 1;
            }

            if position_index >= target_positions.len() {
                break;
            }

            if target_positions[position_index] != pile_pos {
                continue;
            }

            counts.clear();
            umi_consensus.clear();
            let mut n_count: u32 = 0;
            let mut processed: u32 = 0;
            let mut truncated = false;

            for alignment in pileup.alignments() {
                let record = alignment.record();

                if record.mapq() < self.config.min_mapping_quality {
                    continue;
                }

                let qpos = match alignment.qpos() {
                    Some(q) => q,
                    None => continue,
                };

                processed = processed.saturating_add(1);
                if processed > self.config.max_depth {
                    truncated = true;
                    break;
                }

                let base_qual = record.qual().get(qpos).copied().unwrap_or(0);
                if base_qual < self.config.min_base_quality {
                    n_count = n_count.saturating_add(1);
                    continue;
                }

                let base = self.get_base_at_position(&record, Some(qpos))?;
                if base == 'N' {
                    n_count = n_count.saturating_add(1);
                    continue;
                }

                let cell_barcode = self.get_cell_barcode(&record)?;
                if !self.barcode_processor.is_valid(&cell_barcode) {
                    continue;
                }

                let umi = self.get_umi(&record)?;
                if umi == "-" {
                    continue;
                }

                if let Some(encoded) = encode_call(self.config.stranded, base, record.is_reverse())
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

            let current_index = position_index;
            position_index += 1;

            if truncated && self.config.skip_max_depth {
                continue;
            }

            for ((cell_barcode, _umi), encoded) in umi_consensus.drain() {
                if encoded == UMI_CONFLICT_CODE {
                    continue;
                }

                let counts_entry = counts
                    .entry(cell_barcode)
                    .or_insert_with(StrandBaseCounts::default);

                apply_encoded_call(self.config.stranded, encoded, counts_entry);
            }

            if counts.is_empty() {
                continue;
            }

            let mut position_counts: HashMap<String, StrandBaseCounts> =
                HashMap::with_capacity(counts.len());
            position_counts.extend(counts.drain());

            let mut a_total = 0u32;
            let mut t_total = 0u32;
            let mut g_total = 0u32;
            let mut c_total = 0u32;

            for strand_counts in position_counts.values() {
                a_total = a_total
                    .saturating_add(strand_counts.forward.a)
                    .saturating_add(strand_counts.reverse.a);
                t_total = t_total
                    .saturating_add(strand_counts.forward.t)
                    .saturating_add(strand_counts.reverse.t);
                g_total = g_total
                    .saturating_add(strand_counts.forward.g)
                    .saturating_add(strand_counts.reverse.g);
                c_total = c_total
                    .saturating_add(strand_counts.forward.c)
                    .saturating_add(strand_counts.reverse.c);
            }

            let effective_depth = a_total + t_total + g_total + c_total;
            if effective_depth < self.config.min_depth {
                continue;
            }

            let total_depth = effective_depth + n_count;
            let n_limit = if self.config.max_n_fraction == 0 {
                u32::MAX
            } else {
                std::cmp::max(total_depth / self.config.max_n_fraction, 2)
            };

            if n_count > n_limit {
                continue;
            }

            let editing_denom = std::cmp::max(self.config.editing_threshold, 1);
            let valid_value = std::cmp::max(total_depth / editing_denom, 2);
            let mut support = 0;
            if a_total > valid_value {
                support += 1;
            }
            if t_total > valid_value {
                support += 1;
            }
            if g_total > valid_value {
                support += 1;
            }
            if c_total > valid_value {
                support += 1;
            }

            if support < 2 {
                continue;
            }

            let position_meta = &chunk.positions[current_index];
            chunk_results.push(PositionData {
                chrom: position_meta.chrom.clone(),
                pos: position_meta.pos,
                counts: position_counts,
            });
        }

        {
            let mut pool = self.reader_pool.lock();
            pool.push(reader);
        }

        Ok(chunk_results)
    }

    fn get_cell_barcode(&self, record: &rust_htslib::bam::record::Record) -> Result<String> {
        let tag = self.config.cell_barcode_tag.as_bytes();
        match record.aux(tag) {
            Ok(rust_htslib::bam::record::Aux::String(s)) => {
                let clean_barcode = s.split('-').next().unwrap_or(s);
                Ok(clean_barcode.to_string())
            }
            Ok(rust_htslib::bam::record::Aux::ArrayU8(arr)) => {
                let u8_vec: Vec<u8> = arr.iter().collect();
                match std::str::from_utf8(&u8_vec) {
                    Ok(barcode_str) => {
                        let clean_barcode = barcode_str.split('-').next().unwrap_or(barcode_str);
                        Ok(clean_barcode.to_string())
                    }
                    Err(_) => Ok("-".to_string()),
                }
            }
            _ => Ok("-".to_string()),
        }
    }

    fn get_umi(&self, record: &rust_htslib::bam::record::Record) -> Result<String> {
        let tag = self.config.umi_tag.as_bytes();
        match record.aux(tag) {
            Ok(rust_htslib::bam::record::Aux::String(s)) => Ok(s.to_string()),
            Ok(rust_htslib::bam::record::Aux::ArrayU8(arr)) => {
                let u8_vec: Vec<u8> = arr.iter().collect();
                match std::str::from_utf8(&u8_vec) {
                    Ok(s) => Ok(s.to_string()),
                    Err(_) => Ok("-".to_string()),
                }
            }
            _ => Ok("-".to_string()),
        }
    }

    fn get_base_at_position(
        &self,
        record: &rust_htslib::bam::record::Record,
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

/// Run the bam2mtx command with intelligent configuration selection
pub fn run_bam2mtx(args: Bam2MtxArgs) -> Result<()> {
    let start_time = Instant::now();

    log::info!("Starting optimized bam2mtx processing...");
    log::info!("BAM file: {:?}", args.bam);
    log::info!("Barcode whitelist: {:?}", args.barcodes);
    log::info!("Output file: {:?}", args.output);
    if args.two_pass {
        log::info!("Two-pass mode enabled: running bulk first pass to build site list");
    }

    let positions_path = prepare_positions_file(&args)?;
    log::info!("Using position list: {:?}", positions_path);

    if let Err(e) = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
    {
        log::debug!("Thread pool already initialized: {}", e);
    }

    log::info!("Loading cell barcodes...");
    let barcode_processor = Arc::new(BarcodeProcessor::from_file(&args.barcodes)?);
    let estimated_cells = barcode_processor.len();
    log::info!("Loaded {} valid cell barcodes", estimated_cells);

    log::info!("Reading TSV file...");
    let mut positions = read_tsv_file(&positions_path)?;
    if !args.all_contigs {
        let original = positions.len();
        positions = filter_to_standard_contigs(positions);
        log::info!(
            "Filtered positions to canonical contigs: {} → {} entries",
            original,
            positions.len()
        );
    }
    let estimated_positions = positions.len();
    log::info!("Loaded {} positions from TSV", estimated_positions);

    if estimated_positions == 0 {
        return Err(anyhow!(
            "Position list {:?} contains no usable entries",
            positions_path
        ));
    }

    let chunk_size = usize::max(args.chunksize as usize, 1);

    let adata_config = AnnDataConfig {
        stranded: args.stranded,
        compression: Some("gzip".to_string()),
        threads: args.threads,
        chunk_size,
        matrix_density: args.matrix_density,
        batch_size: chunk_size,
    };

    let config = BamProcessorConfig {
        min_mapping_quality: args.min_mapq,
        min_base_quality: args.min_baseq,
        min_depth: args.min_depth,
        max_n_fraction: args.max_n_fraction,
        editing_threshold: args.editing_threshold,
        stranded: args.stranded,
        max_depth: args.max_depth,
        skip_max_depth: args.skip_max_depth,
        umi_tag: args.umi_tag.clone(),
        cell_barcode_tag: args.cb_tag.clone(),
    };

    let chunks = chunk_positions(positions, adata_config.chunk_size);
    let total_chunks = chunks.len();
    log::info!("Split into {} chunks for parallel processing", total_chunks);

    let processor = OptimizedChunkProcessor::new(args.bam.clone(), config, barcode_processor)?;

    log::info!("Processing {} chunks...", total_chunks);
    let processing_start = Instant::now();

    let all_position_data: Vec<PositionData> = chunks
        .into_par_iter()
        .enumerate()
        .try_fold(Vec::new, |mut acc, (idx, chunk)| {
            if idx % 50 == 0 {
                log::debug!("Processing chunk {} of {}", idx + 1, total_chunks);
            }
            let mut data = processor.process_chunk(&chunk)?;
            acc.append(&mut data);
            Ok::<Vec<PositionData>, anyhow::Error>(acc)
        })
        .try_reduce(Vec::new, |mut acc, mut part| {
            acc.append(&mut part);
            Ok::<Vec<PositionData>, anyhow::Error>(acc)
        })?;

    let processing_duration = processing_start.elapsed();
    log::info!(
        "Chunk processing completed in {:?} ({} positions)",
        processing_duration,
        all_position_data.len()
    );

    log::info!("Converting to AnnData format...");
    let converter = AnnDataConverter::new(adata_config);
    let adata = converter.convert(&all_position_data, &args.output)?;

    log::info!("Writing output file...");
    converter.write_to_file(&adata, &args.output)?;

    let duration = start_time.elapsed();
    log::info!("Processing completed in {:?}", duration);
    log::info!("Output written to: {:?}", args.output);

    Ok(())
}

fn prepare_positions_file(args: &Bam2MtxArgs) -> Result<PathBuf> {
    if args.two_pass {
        let target_path = args.tsv.clone().unwrap_or_else(|| {
            args.output
                .parent()
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| PathBuf::from("."))
                .join("1pass.tsv.gz")
        });
        run_bulk_first_pass(args, &target_path)?;
        Ok(target_path)
    } else {
        args.tsv
            .clone()
            .ok_or_else(|| anyhow!("--tsv must be provided unless --two-pass is enabled"))
    }
}

/// Execute the `bulk` subcommand as the first pass of the two-pass workflow.
fn run_bulk_first_pass(args: &Bam2MtxArgs, target: &Path) -> Result<()> {
    log::info!("Running bulk first pass to generate {:?}", target);
    utils::make_parent_dirs(target)?;

    let cli = vec![
        "bulk".to_string(),
        args.bam.to_string_lossy().into_owned(),
        "-o".to_string(),
        target.to_string_lossy().into_owned(),
        "--mapquality".to_string(),
        args.min_mapq.to_string(),
        "-t".to_string(),
        args.threads.to_string(),
        "-c".to_string(),
        args.chunksize.to_string(),
        "-Q".to_string(),
        args.min_baseq.to_string(),
        "--min-depth".to_string(),
        args.min_depth.to_string(),
        "--max-n-fraction".to_string(),
        args.max_n_fraction.to_string(),
        "-D".to_string(),
        "8000".to_string(),
        "--editing-threshold".to_string(),
        args.editing_threshold.to_string(),
    ];

    let cli_refs: Vec<&str> = cli.iter().map(|s| s.as_str()).collect();
    let bulk_args = BulkCommand::from_iter_safe(&cli_refs)?;
    bulk_args.run()
}

/// Discard genomic positions that do not belong to the canonical chromosome list.
fn filter_to_standard_contigs(positions: Vec<GenomicPosition>) -> Vec<GenomicPosition> {
    positions
        .into_iter()
        .filter(|pos| is_standard_contig(&pos.chrom))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bam2mtx_args() {
        let args = Bam2MtxArgs::from_iter_safe(&[
            "test",
            "--bam",
            "test.bam",
            "--tsv",
            "test.tsv",
            "--barcodes",
            "barcodes.tsv",
            "--output",
            "output.h5ad",
        ])
        .unwrap();

        assert_eq!(args.bam, PathBuf::from("test.bam"));
        assert_eq!(args.tsv, Some(PathBuf::from("test.tsv")));
        assert_eq!(args.barcodes, PathBuf::from("barcodes.tsv"));
        assert_eq!(args.output, PathBuf::from("output.h5ad"));
        assert_eq!(args.matrix_density, 0.005);
        assert_eq!(args.min_depth, 10);
        assert_eq!(args.max_n_fraction, 20);
        assert_eq!(args.editing_threshold, 1000);
        assert!(!args.two_pass);
        assert_eq!(args.max_depth, 20000);
        assert!(!args.skip_max_depth);
    }
}
