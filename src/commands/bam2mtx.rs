//! CLI command for bam2mtx functionality with optimized memory management
//!
//! This module provides the command-line interface for converting BAM files to
//! single-cell count matrices. It supports various optimizations for handling
//! large datasets efficiently including parallel processing, memory mapping,
//! and configurable chunk sizes.
//!
//! # Features
//!
//! - Parallel processing of genomic positions
//! - UMI deduplication
//! - Cell barcode validation
//! - Support for stranded and unstranded data
//! - Memory-efficient processing for large datasets
//! - Configurable chunk sizes and matrix density estimation

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use crate::commands::base_depth::Bulk as BulkCommand;
use anyhow::{anyhow, Result};
use rayon::prelude::*;
use rust_htslib::bam::Read;
use structopt::StructOpt;

use redicat_lib::bam2mtx::{
    anndata_output::{AnnDataConfig, AnnDataConverter},
    barcode::BarcodeProcessor,
    processor::{BamProcessorConfig, PositionData, StrandBaseCounts},
};
use redicat_lib::utils;

/// Arguments for the bam2mtx command
#[derive(Debug, StructOpt)]
#[structopt(name = "bam2mtx", about = "Convert BAM files to single-cell matrices")]
pub struct Bam2MtxArgs {
    /// Path to the input BAM file (must be indexed)
    #[structopt(short, long, parse(from_os_str))]
    pub bam: PathBuf,

    /// Optional TSV file with genomic positions (CHR/POS). 若未提供且启用 two-pass，将自动生成。
    #[structopt(long, parse(from_os_str))]
    pub tsv: Option<PathBuf>,

    /// Path to the cell barcodes file
    #[structopt(long, parse(from_os_str))]
    pub barcodes: PathBuf,

    /// Enable two-pass模式：先运行 `bulk` 生成 1pass.tsv.gz，再执行 bam2mtx。
    #[structopt(long)]
    pub two_pass: bool,

    /// Output path for the H5AD file
    #[structopt(short, long, parse(from_os_str))]
    pub output: PathBuf,

    /// Number of threads to use (default: 10)
    #[structopt(short, long, default_value = "10")]
    pub threads: usize,

    /// Minimum mapping quality
    #[structopt(long, default_value = "255")]
    pub min_mapq: u8,

    /// Minimum base quality
    #[structopt(long, default_value = "30")]
    pub min_baseq: u8,

    /// Whether data is stranded (default: unstranded)
    #[structopt(long)]
    pub stranded: bool,

    /// UMI tag name
    #[structopt(long, default_value = "UB")]
    pub umi_tag: String,

    /// Cell barcode tag name
    #[structopt(long, default_value = "CB")]
    pub cb_tag: String,

    /// Path to reference FASTA file (required for CRAM files)
    #[structopt(long, parse(from_os_str))]
    pub reference: Option<PathBuf>,

    /// Chunk size for parallel processing
    /// Controls the number of genomic positions processed in each parallel chunk.
    /// Larger values reduce scheduling overhead but increase memory usage.
    /// Smaller values improve load balancing but may increase scheduling overhead.
    #[structopt(long, default_value = "2500")]
    pub chunksize: u32,

    /// Matrix density estimation for memory优化
    /// Typical values: sparse matrix 0.001-0.01, medium density 0.01-0.1, dense matrix 0.1+
    #[structopt(long, default_value = "0.005")]
    pub matrix_density: f64,
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
    use std::fs::File;
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
    let mut current_chrom = positions[0].chrom.clone();

    for pos in positions {
        // Adaptive chunking: smaller chunks within chromosomes for better cache locality
        let effective_chunk_size = if pos.chrom == current_chrom {
            std::cmp::min(chunk_size, 5000) // Smaller chunks within chromosome
        } else {
            current_chrom = pos.chrom.clone();
            chunk_size
        };

        if current_chunk.len() >= effective_chunk_size {
            chunks.push(PositionChunk {
                positions: std::mem::replace(&mut current_chunk, Vec::with_capacity(chunk_size)),
            });
        }
        current_chunk.push(pos);
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
}

impl OptimizedChunkProcessor {
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

    pub fn process_chunk(&self, chunk: &PositionChunk) -> Result<Vec<PositionData>> {
        use rust_htslib::bam::{self, Read};

        let mut reader = bam::IndexedReader::from_path(&self.bam_path)?;
        let mut chunk_results = Vec::with_capacity(chunk.positions.len());

        // Pre-allocate reusable containers to reduce allocations
        let mut counts = std::collections::HashMap::with_capacity(1000);
        let mut umi_consensus = std::collections::HashMap::with_capacity(10000);

        for pos_data in &chunk.positions {
            counts.clear();
            umi_consensus.clear();

            let chrom = pos_data.chrom.as_str();
            let pos = pos_data.pos;

            let tid = match reader.header().tid(chrom.as_bytes()) {
                Some(tid) => tid,
                None => continue,
            };

            if let Ok(position_data) = self.process_single_position_optimized(
                &mut reader,
                tid,
                pos - 1,
                chrom,
                &mut counts,
                &mut umi_consensus,
            ) {
                chunk_results.push(position_data);
            }
        }

        Ok(chunk_results)
    }

    fn process_single_position_optimized(
        &self,
        reader: &mut rust_htslib::bam::IndexedReader,
        tid: u32,
        pos: u64,
        chrom: &str,
        counts: &mut std::collections::HashMap<String, StrandBaseCounts>,
        umi_consensus: &mut std::collections::HashMap<String, String>,
    ) -> Result<PositionData> {
        reader.fetch((tid, pos as u32, pos as u32 + 1))?;

        for pileup in reader.pileup() {
            let pileup = pileup?;
            if pileup.pos() != pos as u32 {
                continue;
            }

            for alignment in pileup.alignments() {
                let record = alignment.record();

                if !self.should_process_alignment(&record, alignment.qpos())? {
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
                let base = self.get_base_at_position(&record, alignment.qpos())?;
                let strand = if record.is_reverse() { "-" } else { "+" };

                // Use string interning to reduce memory allocations
                let key = format!("{}:{}", cell_barcode, umi);
                let base_strand = if self.config.stranded {
                    format!("{}{}", base, strand)
                } else {
                    base.to_string()
                };

                umi_consensus
                    .entry(key.clone())
                    .and_modify(|existing| {
                        if *existing != base_strand {
                            *existing = "-".to_string();
                        }
                    })
                    .or_insert(base_strand);
            }
        }

        // Process UMI consensus with optimized aggregation
        for (key, base_strand) in umi_consensus.drain() {
            if base_strand == "-" {
                continue;
            }

            let separator_pos = key.find(':').unwrap_or(key.len());
            let cell_barcode = &key[..separator_pos];

            let counts_entry = counts
                .entry(cell_barcode.to_string())
                .or_insert_with(StrandBaseCounts::default);

            // Optimized base counting with match on bytes for better performance
            if self.config.stranded {
                match base_strand.as_bytes() {
                    [b'A', b'+'] => counts_entry.forward.a += 1,
                    [b'T', b'+'] => counts_entry.forward.t += 1,
                    [b'G', b'+'] => counts_entry.forward.g += 1,
                    [b'C', b'+'] => counts_entry.forward.c += 1,
                    [b'A', b'-'] => counts_entry.reverse.a += 1,
                    [b'T', b'-'] => counts_entry.reverse.t += 1,
                    [b'G', b'-'] => counts_entry.reverse.g += 1,
                    [b'C', b'-'] => counts_entry.reverse.c += 1,
                    _ => {}
                }
            } else {
                match base_strand.as_bytes()[0] {
                    b'A' => counts_entry.forward.a += 1,
                    b'T' => counts_entry.forward.t += 1,
                    b'G' => counts_entry.forward.g += 1,
                    b'C' => counts_entry.forward.c += 1,
                    _ => {}
                }
            }
        }

        Ok(PositionData {
            chrom: chrom.to_string(),
            pos: pos + 1,
            counts: std::mem::take(counts),
        })
    }

    fn should_process_alignment(
        &self,
        record: &rust_htslib::bam::record::Record,
        qpos: Option<usize>,
    ) -> Result<bool> {
        if let Some(qpos) = qpos {
            if let Some(qual) = record.qual().get(qpos) {
                if *qual < self.config.min_base_quality {
                    return Ok(false);
                }
            }
        } else {
            return Ok(false);
        }

        if record.mapq() < self.config.min_mapping_quality {
            return Ok(false);
        }

        Ok(true)
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
    let positions = read_tsv_file(&positions_path)?;
    let estimated_positions = positions.len();
    log::info!("Loaded {} positions from TSV", estimated_positions);

    if estimated_positions == 0 {
        return Err(anyhow!(
            "Position list {:?} contains no usable entries",
            positions_path
        ));
    }

    let auto_use_mmap = estimated_positions >= 1_000_000 || estimated_cells >= 200_000;
    let chunk_size = if auto_use_mmap {
        std::cmp::max(args.chunksize as usize, 10_000)
    } else {
        args.chunksize as usize
    };

    let adata_config = AnnDataConfig {
        stranded: args.stranded,
        compression: Some("gzip".to_string()),
        threads: args.threads,
        chunk_size,
        matrix_density: args.matrix_density,
        use_mmap: auto_use_mmap,
        batch_size: chunk_size,
    };

    log::info!("Configuration summary:");
    log::info!("  - Threads: {}", adata_config.threads);
    log::info!("  - Chunk size: {}", adata_config.chunk_size);
    log::info!("  - Matrix density: {:.4}", adata_config.matrix_density);
    log::info!("  - Use mmap: {}", adata_config.use_mmap);
    log::info!("  - Batch size: {}", adata_config.batch_size);
    log::info!("  - Stranded: {}", adata_config.stranded);

    let config = BamProcessorConfig {
        min_mapping_quality: args.min_mapq,
        min_base_quality: args.min_baseq,
        stranded: args.stranded,
        umi_tag: args.umi_tag.clone(),
        cell_barcode_tag: args.cb_tag.clone(),
    };

    let chunks = chunk_positions(positions, adata_config.chunk_size);
    let total_chunks = chunks.len();
    log::info!("Split into {} chunks for parallel processing", total_chunks);

    let processor = OptimizedChunkProcessor::new(args.bam.clone(), config, barcode_processor);

    log::info!("Processing {} chunks...", total_chunks);
    let processing_start = Instant::now();

    let all_results: Vec<Vec<PositionData>> = chunks
        .into_par_iter()
        .enumerate()
        .map(|(idx, chunk)| {
            if idx % 50 == 0 {
                log::debug!("Processing chunk {} of {}", idx + 1, total_chunks);
            }
            processor.process_chunk(&chunk).unwrap_or_default()
        })
        .collect();

    let processing_duration = processing_start.elapsed();
    log::info!("Chunk processing completed in {:?}", processing_duration);

    let total_length: usize = all_results.iter().map(|v| v.len()).sum();
    let mut all_position_data = Vec::with_capacity(total_length);
    for chunk_result in all_results {
        all_position_data.extend(chunk_result);
    }
    log::info!("Processed {} positions total", all_position_data.len());

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

fn run_bulk_first_pass(args: &Bam2MtxArgs, target: &Path) -> Result<()> {
    log::info!("Running bulk first pass to generate {:?}", target);
    utils::make_parent_dirs(target)?;

    let mut cli = vec![
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
        "--all".to_string(),
    ];

    cli.push("--barcodes".to_string());
    cli.push(args.barcodes.to_string_lossy().into_owned());

    let cli_refs: Vec<&str> = cli.iter().map(|s| s.as_str()).collect();
    let bulk_args = BulkCommand::from_iter_safe(&cli_refs)?;
    bulk_args.run()
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
        assert!(!args.two_pass);
    }
}
