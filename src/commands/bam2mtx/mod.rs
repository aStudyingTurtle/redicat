mod args;
mod engine;
mod input;
mod workflow;

use anyhow::{anyhow, Result};
use log::info;
use rayon::prelude::*;
use redicat_lib::bam2mtx::anndata_output::{AnnDataConfig, AnnDataConverter};
use redicat_lib::bam2mtx::barcode::BarcodeProcessor;
use redicat_lib::bam2mtx::processor::{BamProcessorConfig, PositionData};
use redicat_lib::utils;
use std::fs::File;
use std::io::Write;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

use crate::commands::{common, is_standard_contig};

pub use args::Bam2MtxArgs;
use engine::{OptimizedChunkProcessor, SkippedSite};
use input::{chunk_positions, filter_positions_by, read_positions};
use workflow::prepare_positions_file;

struct ChunkOutcome {
    positions: Vec<PositionData>,
    skipped: Vec<SkippedSite>,
}

/// Entry point for the `bam2mtx` command.
pub fn run_bam2mtx(args: Bam2MtxArgs) -> Result<()> {
    let start_time = Instant::now();

    info!("Starting bam2mtx processing for {:?}", args.bam);
    if args.two_pass {
        info!("Two-pass mode enabled – running bulk first");
    }

    let positions_path = prepare_positions_file(&args)?;
    info!("Using position list: {:?}", positions_path);

    let active_threads = common::configure_global_thread_pool(args.threads)?;

    utils::make_parent_dirs(&args.output)?;

    info!("Loading {} threads into Rayon pool", active_threads);
    info!("Loading cell barcodes from {:?}", args.barcodes);
    let barcode_processor = Arc::new(BarcodeProcessor::from_file(&args.barcodes)?);
    info!("Loaded {} valid barcodes", barcode_processor.len());

    info!("Reading positions...");
    let mut positions = read_positions(&positions_path)?;
    if !args.all_contigs {
        let before = positions.len();
        positions = filter_positions_by(positions, is_standard_contig);
        info!(
            "Filtered positions to canonical contigs: {} → {} entries",
            before,
            positions.len()
        );
    }

    if positions.is_empty() {
        return Err(anyhow!(
            "Position list {:?} contains no usable entries",
            positions_path
        ));
    }

    let manifest_positions = positions.len();

    let chunk_size = args.chunk_size();
    let chunk_size_max_depth = args.chunk_size_max_depth();
    let adata_config = AnnDataConfig {
        stranded: args.stranded,
        compression: Some("gzip".to_string()),
        threads: active_threads,
        chunk_size,
        matrix_density: args.matrix_density,
        batch_size: chunk_size,
        total_positions: manifest_positions,
        triplet_spill_nnz: (chunk_size as usize * 16).max(500_000),
    };

    let processor_config = BamProcessorConfig {
        min_mapping_quality: args.min_mapq,
        min_base_quality: args.min_baseq,
        min_depth: args.min_depth,
        max_n_fraction: args.max_n_fraction,
        editing_threshold: args.editing_threshold,
        stranded: args.stranded,
        max_depth: args.max_depth,
        umi_tag: args.umi_tag.clone(),
        cell_barcode_tag: args.cb_tag.clone(),
    };

    info!(
        "Chunking positions with chunk_size={} and chunk_size_max_depth={}",
        chunk_size, chunk_size_max_depth
    );

    let chunks = chunk_positions(positions, chunk_size, chunk_size_max_depth);
    let total_near_max_positions: usize = chunks
        .iter()
        .map(|chunk| chunk.near_max_depth_count())
        .sum();
    info!(
        "Chunked positions into {} batches ({} near-max-depth positions)",
        chunks.len(),
        total_near_max_positions
    );

    let processor = OptimizedChunkProcessor::new(
        args.bam.clone(),
        processor_config,
        Arc::clone(&barcode_processor),
    )?;
    let contig_names = processor.contig_names();

    info!("Processing chunks with parallel aggregation...");
    let processing_start = Instant::now();
    let total_chunks = chunks.len();
    let processed_chunks = AtomicUsize::new(0);
    let total_positions = AtomicUsize::new(0);
    let near_max_remaining = AtomicUsize::new(total_near_max_positions);
    let long_report_interval = Duration::from_secs(1800);
    let last_timed_report_secs = AtomicU64::new(0);
    let last_percent_bucket = AtomicUsize::new(0);

    let chunk_results = chunks
        .into_par_iter()
        .enumerate()
        .try_fold(
            || Vec::new(),
            |mut local, (_idx, chunk)| -> Result<Vec<ChunkOutcome>> {
                let chunk_near = chunk.near_max_depth_count();
                let (data, skipped) = processor.process_chunk(&chunk)?;

                if !data.is_empty() {
                    total_positions.fetch_add(data.len(), Ordering::Relaxed);
                }

                let near_remaining = if chunk_near == 0 {
                    near_max_remaining.load(Ordering::Relaxed)
                } else {
                    near_max_remaining
                        .fetch_sub(chunk_near, Ordering::Relaxed)
                        .saturating_sub(chunk_near)
                };

                let completed = processed_chunks.fetch_add(1, Ordering::Relaxed) + 1;
                let elapsed = processing_start.elapsed();
                let elapsed_secs = elapsed.as_secs();

                if total_chunks > 0 {
                    let percent = (completed as f64 * 100.0)
                        / total_chunks.max(1) as f64;
                    let bucket = ((percent / 5.0).floor() as usize).min(20);

                    let mut percent_log_due = false;
                    if completed < total_chunks {
                        let previous = last_percent_bucket.load(Ordering::Relaxed);
                        if bucket > previous
                            && last_percent_bucket
                                .compare_exchange(
                                    previous,
                                    bucket,
                                    Ordering::Relaxed,
                                    Ordering::Relaxed,
                                )
                                .is_ok()
                        {
                            percent_log_due = true;
                        }
                    }

                    let mut timed_log_due = false;
                    if completed < total_chunks {
                        let last_marker = last_timed_report_secs.load(Ordering::Relaxed);
                        if elapsed_secs
                            .saturating_sub(last_marker)
                            >= long_report_interval.as_secs()
                            && last_timed_report_secs
                                .compare_exchange(
                                    last_marker,
                                    elapsed_secs,
                                    Ordering::Relaxed,
                                    Ordering::Relaxed,
                                )
                                .is_ok()
                        {
                            timed_log_due = true;
                        }
                    }

                    if percent_log_due || timed_log_due {
                        info!(
                            "Processed {:.1}% ({} / {} chunks, {} near-max-depth positions remaining)",
                            percent,
                            completed,
                            total_chunks,
                            near_remaining
                        );
                    }
                }

                local.push(ChunkOutcome { positions: data, skipped });
                Ok(local)
            },
        )
        .try_reduce(
            || Vec::new(),
            |mut acc, mut local| {
                acc.append(&mut local);
                Ok(acc)
            },
        )?;

    info!("Processed 100.0%");

    info!(
        "Chunk traversal finished in {:?} ({} positions)",
        processing_start.elapsed(),
        total_positions.load(Ordering::Relaxed)
    );

    let mut skipped_records: Vec<SkippedSite> = Vec::new();
    let mut position_batches: Vec<Vec<PositionData>> = Vec::with_capacity(chunk_results.len());
    for outcome in chunk_results {
        let ChunkOutcome { positions, skipped } = outcome;
        if !skipped.is_empty() {
            skipped_records.extend(skipped);
        }
        if !positions.is_empty() {
            position_batches.push(positions);
        }
    }

    let output_path = args.output.clone();
    let converter = AnnDataConverter::new(adata_config, Arc::clone(&barcode_processor), contig_names);
    info!(
        "Assembling sparse matrices from {} result batches...",
        position_batches.len()
    );
    let parallel_assembly_start = Instant::now();
    let adata = converter.convert_parallel_chunks(position_batches, &output_path)?;
    info!(
        "Sparse matrix assembly completed in {:?}",
        parallel_assembly_start.elapsed()
    );
    converter.write_to_file(&adata, &output_path)?;
    info!("AnnData write finished -> {:?}", output_path);

    let skipped_count = skipped_records.len();
    let mut skipped_path = output_path.clone();
    let stem = skipped_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("output");
    skipped_path.set_file_name(format!("{}_skiped_sites.txt", stem));
    let mut file = File::create(&skipped_path)?;
    writeln!(file, "CHR\tPOS\tDEPTH")?;
    for site in skipped_records {
        writeln!(file, "{}\t{}\t{}", site.contig, site.pos, site.depth)?;
    }
    info!(
        "Recorded {} depth-capped sites at {:?}",
        skipped_count, skipped_path
    );

    info!(
        "Chunk processing completed in {:?} ({} positions)",
        processing_start.elapsed(),
        total_positions.load(Ordering::Relaxed)
    );

    info!("bam2mtx workflow finished in {:?}", start_time.elapsed());

    Ok(())
}
