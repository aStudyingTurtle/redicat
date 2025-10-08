mod args;
mod engine;
mod input;
mod workflow;

use anyhow::{anyhow, Result};
use crossbeam::channel::bounded;
use log::info;
use rayon::prelude::*;
use redicat_lib::bam2mtx::anndata_output::{AnnDataConfig, AnnDataConverter};
use redicat_lib::bam2mtx::barcode::BarcodeProcessor;
use redicat_lib::bam2mtx::processor::BamProcessorConfig;
use redicat_lib::utils;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::Instant;

use crate::commands::{common, is_standard_contig};

pub use args::Bam2MtxArgs;
use engine::OptimizedChunkProcessor;
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
        triplet_spill_nnz: (chunk_size as usize * 32).max(500_000),
    };

    let processor_config = BamProcessorConfig {
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

    info!(
        "Chunking positions with chunk_size={} and chunk_size_max_depth={}",
        chunk_size, chunk_size_max_depth
    );

    let chunks = chunk_positions(positions, chunk_size, chunk_size_max_depth);
    info!("Chunked positions into {} batches", chunks.len());

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
    let log_step = usize::max(1, total_chunks.max(1) / 10);

    let channel_capacity = usize::max(active_threads.saturating_mul(2), 1);
    let (sender, receiver) = bounded(channel_capacity);
    let output_path = args.output.clone();
    let converter = AnnDataConverter::new(
        adata_config,
        Arc::clone(&barcode_processor),
        contig_names,
    );

    let writer_handle = thread::spawn(move || -> Result<()> {
        info!(
            "Streaming AnnData writer listening with buffer capacity {} entries",
            channel_capacity
        );
        let adata = converter.convert_streaming(
            receiver
                .into_iter()
                .flat_map(|chunk: Vec<_>| chunk.into_iter()),
            &output_path,
        )?;
        converter.write_to_file(&adata, &output_path)?;
        info!("AnnData streaming write finished -> {:?}", output_path);
        Ok(())
    });

    let sender_clone = sender.clone();
    chunks.into_par_iter().enumerate().try_for_each_init(
        || sender_clone.clone(),
        |tx, (idx, chunk)| -> Result<()> {
                let data = processor.process_chunk(&chunk)?;
                let chunk_positions = data.len();
                if chunk_positions > 0 {
                    tx.send(data)
                    .map_err(|err| anyhow!("failed to send chunk {idx}: {err}"))?;
                total_positions.fetch_add(chunk_positions, Ordering::Relaxed);
            }

            let completed = processed_chunks.fetch_add(1, Ordering::Relaxed) + 1;
            if completed == total_chunks || completed.is_multiple_of(log_step) {
                let percent = (completed as f64 / total_chunks.max(1) as f64) * 100.0;
                info!(
                    "Processed {:.1}% ({} / {} chunks)",
                    percent, completed, total_chunks
                );
            }

            Ok(())
        },
    )?;

    drop(sender);

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

    info!(
        "Chunk processing completed in {:?} ({} positions)",
        processing_start.elapsed(),
        total_positions.load(Ordering::Relaxed)
    );

    info!("bam2mtx workflow finished in {:?}", start_time.elapsed());

    Ok(())
}
