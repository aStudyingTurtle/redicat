mod args;
mod engine;
mod input;
mod workflow;

use anyhow::{anyhow, Result};
use crossbeam::channel::bounded;
use log::info;
use parking_lot::Mutex;
use rayon::prelude::*;
use redicat_lib::bam2mtx::anndata_output::{AnnDataConfig, AnnDataConverter};
use redicat_lib::bam2mtx::barcode::BarcodeProcessor;
use redicat_lib::bam2mtx::processor::BamProcessorConfig;
use redicat_lib::utils;
use std::fs::File;
use std::io::Write;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};

use crate::commands::{common, is_standard_contig};

pub use args::Bam2MtxArgs;
use engine::{OptimizedChunkProcessor, SkippedSite};
use input::{chunk_positions, filter_positions_by, read_positions};
use workflow::prepare_positions_file;

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

    info!("Processing chunks in parallel with streaming AnnData conversion...");
    let processing_start = Instant::now();
    let total_chunks = chunks.len();
    let processed_chunks = AtomicUsize::new(0);
    let total_positions = AtomicUsize::new(0);
    let log_step = usize::max(1, total_chunks.max(1) / 20);
    let near_max_remaining = AtomicUsize::new(total_near_max_positions);
    let long_report_interval = Duration::from_secs(1800);
    let last_long_report = AtomicU64::new(0);

    let channel_capacity = usize::max(active_threads.saturating_mul(2), 1);
    info!(
        "Creating bounded channel with capacity {} (threads={}, multiplier=2)",
        channel_capacity, active_threads
    );
    let (sender, receiver) = bounded(channel_capacity);
    let output_path = args.output.clone();
    let converter =
        AnnDataConverter::new(adata_config, Arc::clone(&barcode_processor), contig_names);
    let output_path_writer = output_path.clone();
    let skipped_sites: Arc<Mutex<Vec<SkippedSite>>> = Arc::new(Mutex::new(Vec::new()));

    let writer_handle = thread::spawn(move || -> Result<()> {
        info!(
            "Streaming AnnData writer listening with buffer capacity {} entries",
            channel_capacity
        );
        let adata = converter.convert_streaming(
            receiver
                .into_iter()
                .flat_map(|chunk: Vec<_>| chunk.into_iter()),
            &output_path_writer,
        )?;
        converter.write_to_file(&adata, &output_path_writer)?;
        info!(
            "AnnData streaming write finished -> {:?}",
            output_path_writer
        );
        Ok(())
    });

    // CRITICAL FIX: Ensure all sender clones are dropped before waiting for writer thread
    // Use a scoped block to guarantee cleanup of all sender clones from Rayon threads
    {
        let sender_clone = sender.clone();
        let skipped_sites_clone = Arc::clone(&skipped_sites);
        chunks.into_par_iter().enumerate().try_for_each_init(
            || (sender_clone.clone(), Arc::clone(&skipped_sites_clone)),
            |(tx, skipped_bucket), (idx, chunk)| -> Result<()> {
                let chunk_near = chunk.near_max_depth_count();
                let (data, skipped) = processor.process_chunk(&chunk)?;
                if !skipped.is_empty() {
                    let mut guard = skipped_bucket.lock();
                    guard.extend(skipped);
                }
                let chunk_positions = data.len();
                if chunk_positions > 0 {
                    tx.send(data)
                        .map_err(|err| anyhow!("failed to send chunk {idx}: {err}"))?;
                    total_positions.fetch_add(chunk_positions, Ordering::Relaxed);
                }

                let near_remaining = if chunk_near == 0 {
                    near_max_remaining.load(Ordering::Relaxed)
                } else {
                    near_max_remaining
                        .fetch_sub(chunk_near, Ordering::Relaxed)
                        .saturating_sub(chunk_near)
                };

                let completed = processed_chunks.fetch_add(1, Ordering::Relaxed) + 1;
                let percent = if completed < total_chunks {
                    let scaled = (completed as f64 * 1000.0)
                        / total_chunks.max(1) as f64;
                    (scaled.floor()) / 10.0
                } else {
                    100.0
                };
                if completed == total_chunks || completed.is_multiple_of(log_step) {
                    info!(
                        "Processed {:.1}% ({} / {} chunks, {} near-max-depth positions remaining)",
                        percent,
                        completed,
                        total_chunks,
                        near_remaining
                    );
                }

                let elapsed = processing_start.elapsed();
                let elapsed_secs = elapsed.as_secs();
                let last_marker = last_long_report.load(Ordering::Relaxed);
                if elapsed_secs.saturating_sub(last_marker) >= long_report_interval.as_secs() {
                    if last_long_report
                        .compare_exchange(
                            last_marker,
                            elapsed_secs,
                            Ordering::Relaxed,
                            Ordering::Relaxed,
                        )
                        .is_ok()
                    {
                        let remaining_chunks = total_chunks.saturating_sub(completed);
                        info!(
                            "Long-running status: {} chunks remaining, {} near-max-depth positions pending (elapsed {:?})",
                            remaining_chunks,
                            near_remaining,
                            elapsed
                        );
                    }
                }

                Ok(())
            },
        )?;
        // sender_clone is dropped here when the scope ends
    }
    
    // Now drop the original sender - this ensures ALL senders are dropped
    drop(sender);
    info!("All chunk senders dropped, waiting for writer thread to complete...");

    let writer_result = writer_handle.join().map_err(|err| {
        let panic_msg = if let Some(msg) = err.downcast_ref::<&str>() {
            *msg
        } else if let Some(msg) = err.downcast_ref::<String>() {
            msg.as_str()
        } else {
            "unknown panic"
        };
        anyhow!("AnnData writer thread panicked: {panic_msg}")
    })?;
    writer_result?;

    let skipped_records = {
        let mut guard = skipped_sites.lock();
        std::mem::take(&mut *guard)
    };
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
