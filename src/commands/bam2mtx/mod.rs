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
use input::{chunk_positions, filter_positions_by, read_positions, PositionChunk};
use workflow::prepare_positions_file;

/// Calculate adaptive channel capacity based on workload characteristics.
/// 
/// This function estimates optimal channel capacity considering:
/// - Thread count (parallelism level)
/// - Chunk sizes and weights (memory footprint)
/// - Available system memory
/// - Total number of chunks (queue depth)
///
/// Returns a capacity that balances throughput and memory usage.
fn estimate_channel_capacity(
    threads: usize,
    chunks: &[PositionChunk],
    available_memory_gb: u64,
) -> usize {
    if chunks.is_empty() {
        return threads * 2;  // Fallback to simple heuristic
    }

    // Calculate average chunk weight (positions × depth)
    let total_weight: u64 = chunks.iter().map(|c| c.total_weight()).sum();
    let avg_weight = total_weight / chunks.len() as u64;

    // Estimate memory per buffered chunk (rough approximation)
    // Each position typically uses ~100-500 bytes depending on depth
    let bytes_per_weight_unit = 200_u64;  // Conservative estimate
    let avg_chunk_memory_mb = (avg_weight * bytes_per_weight_unit) / (1024 * 1024);

    // Base capacity: threads * 2 (producer-consumer pattern)
    let base_capacity = threads.saturating_mul(2);

    // Maximum by memory: don't use more than 10% of available memory for buffering
    let max_by_memory = if avg_chunk_memory_mb > 0 {
        ((available_memory_gb * 1024 * 10) / 100 / avg_chunk_memory_mb) as usize
    } else {
        base_capacity * 4
    };

    // Maximum by chunk count: don't buffer more than 25% of total chunks
    let max_by_chunks = chunks.len() / 4;

    // Choose the minimum of all constraints, but at least base_capacity
    let capacity = base_capacity
        .max(8)  // Minimum 8 to avoid thrashing
        .min(max_by_memory)
        .min(max_by_chunks)
        .min(1024);  // Hard cap at 1024 to prevent excessive memory

    capacity
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

    info!("Processing chunks in parallel with streaming AnnData conversion...");
    let processing_start = Instant::now();
    let total_chunks = chunks.len();
    let processed_chunks = AtomicUsize::new(0);
    let total_positions = AtomicUsize::new(0);
    let log_step = usize::max(1, total_chunks.max(1) / 20);
    let near_max_remaining = AtomicUsize::new(total_near_max_positions);
    let long_report_interval = Duration::from_secs(1800);
    let last_long_report = AtomicU64::new(0);

    // Adaptive channel capacity based on workload characteristics
    let available_memory_gb = {
        // Try to get system memory, fallback to conservative estimate
        #[cfg(target_os = "linux")]
        {
            std::fs::read_to_string("/proc/meminfo")
                .ok()
                .and_then(|content| {
                    content.lines()
                        .find(|line| line.starts_with("MemAvailable:"))
                        .and_then(|line| {
                            line.split_whitespace()
                                .nth(1)
                                .and_then(|s| s.parse::<u64>().ok())
                                .map(|kb| kb / (1024 * 1024))  // KB to GB
                        })
                })
                .unwrap_or(8)  // Default to 8GB if detection fails
        }
        #[cfg(not(target_os = "linux"))]
        {
            8  // Conservative default for non-Linux systems
        }
    };

    let channel_capacity = estimate_channel_capacity(
        active_threads,
        &chunks,
        available_memory_gb,
    );
    
    info!(
        "Creating adaptive bounded channel: capacity={} (threads={}, chunks={}, memory={}GB)",
        channel_capacity, active_threads, chunks.len(), available_memory_gb
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
