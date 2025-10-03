mod args;
mod engine;
mod input;
mod workflow;

use anyhow::{anyhow, Result};
use log::info;
use rayon::prelude::*;
use redicat_lib::bam2mtx::anndata_output::{AnnDataConfig, AnnDataConverter};
use redicat_lib::bam2mtx::barcode::BarcodeProcessor;
use redicat_lib::bam2mtx::processor::BamProcessorConfig;
use redicat_lib::utils;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
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

    let chunk_size = args.chunk_size();
    let adata_config = AnnDataConfig {
        stranded: args.stranded,
        compression: Some("gzip".to_string()),
        threads: active_threads,
        chunk_size,
        matrix_density: args.matrix_density,
        batch_size: chunk_size,
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

    let chunks = chunk_positions(positions, chunk_size);
    info!("Chunked positions into {} batches", chunks.len());

    let processor =
        OptimizedChunkProcessor::new(args.bam.clone(), processor_config, barcode_processor)?;

    info!("Processing chunks in parallel...");
    let processing_start = Instant::now();
    let total_chunks = chunks.len();
    let processed_chunks = AtomicUsize::new(0);
    let log_step = usize::max(1, total_chunks / 10);
    let all_position_data: Vec<_> = chunks
        .into_par_iter()
        .enumerate()
        .try_fold(Vec::new, |mut acc, (_idx, chunk)| {
            let mut data = processor.process_chunk(&chunk)?;
            let completed = processed_chunks.fetch_add(1, Ordering::Relaxed) + 1;
            if completed == total_chunks || completed.is_multiple_of(log_step) {
                let percent = (completed as f64 / total_chunks as f64) * 100.0;
                info!(
                    "Processed {:.1}% ({} / {} chunks)",
                    percent,
                    completed,
                    total_chunks
                );
            }
            acc.append(&mut data);
            Ok::<Vec<_>, anyhow::Error>(acc)
        })
        .try_reduce(Vec::new, |mut acc, mut part| {
            acc.append(&mut part);
            Ok::<Vec<_>, anyhow::Error>(acc)
        })?;

    info!(
        "Chunk processing completed in {:?} ({} positions)",
        processing_start.elapsed(),
        all_position_data.len()
    );

    info!("Converting to AnnData and writing {:?}", args.output);
    let converter = AnnDataConverter::new(adata_config);
    let adata = converter.convert(&all_position_data, &args.output)?;
    converter.write_to_file(&adata, &args.output)?;

    info!("bam2mtx workflow finished in {:?}", start_time.elapsed());

    Ok(())
}
