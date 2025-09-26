//! RNA Editing Analysis Pipeline
//!
//! This module provides the command-line interface for the RNA editing detection
//! and analysis pipeline. It processes single-cell RNA-seq data to identify RNA
//! editing events and calculate editing indices in a strand-aware manner.
//!
//! # Pipeline Steps
//!
//! 1. Variant annotation using REDIPortal database
//! 2. Base matrix processing and filtering
//! 3. Reference base assignment
//! 4. Mismatch filtering
//! 5. Ref/Alt matrix calculation with strand-aware processing
//! 6. Editing index calculation
//! 7. Site-level mismatch analysis
//!
//! # Strand-Aware Processing
//!
//! The analysis pipeline performs strand-aware processing to correctly handle
//! RNA editing events on both positive and negative DNA strands. This is important
//! because the same editing event can be observed from either strand due to
//! complementary base pairing:
//! - A>G editing on the positive strand appears as T>C editing on the negative strand
//! - A>C editing on the positive strand appears as T>G editing on the negative strand
//! - And so on for all editing types
//!
//! This ensures comprehensive detection of RNA editing events regardless of the
//! DNA strand they originate from.

use anyhow::Result;
use log::{error, info};
use num_cpus;
use std::fs;
use std::path::Path;
use std::sync::Arc;
use structopt::StructOpt;

use redicat_lib::call::anndata_ops::{read_anndata_h5ad, write_anndata_h5ad};
use redicat_lib::call::annotate_variants_pipeline;
use redicat_lib::call::calculate_cei;
use redicat_lib::call::calculate_ref_alt_matrices;
use redicat_lib::call::calculate_site_mismatch_stats;
use redicat_lib::call::editing::load_rediportal_parallel;
use redicat_lib::call::editing::EditingType;
use redicat_lib::call::error::RedicatError;
use redicat_lib::call::reference_genome::ReferenceGenome;
use redicat_lib::call::validation::{validate_input_files, validate_output_path, ValidationConfig};

#[derive(StructOpt, Debug)]
#[structopt(name = "call")]
#[structopt(about = "REDICAT analysis pipeline - Rust implementation")]
#[structopt(version = "0.1.0")]
pub struct CallArgs {
    #[structopt(long, help = "Input AnnData (.h5ad) file path")]
    pub input: String,

    #[structopt(long, help = "Output AnnData (.h5ad) file path")]
    pub output: String,

    #[structopt(long, help = "Reference genome FASTA file path")]
    pub fa: String,

    #[structopt(
        long,
        help = "Site white list TSV file path (at least containing CHR and POS columns as the first two columns)"
    )]
    pub site_white_list: String,

    #[structopt(
        long,
        default_value = "ag",
        help = "Editing type, one of ag, ac, at, ca, cg, ct"
    )]
    pub editingtype: EditingType,

    #[structopt(
        long,
        default_value = "0.01",
        help = "Maximum threshold for other mismatches"
    )]
    pub max_other_threshold: f64,

    #[structopt(
        long,
        default_value = "0.01",
        help = "Minimum threshold for edited mismatches"
    )]
    pub min_edited_threshold: f64,

    #[structopt(
        long,
        default_value = "0.01",
        help = "Minimum threshold for reference mismatches"
    )]
    pub min_ref_threshold: f64,

    #[structopt(
        long,
        short = "c",
        default_value = "100000",
        help = "Chunk size for parallel processing"
    )]
    pub chunksize: usize,

    #[structopt(short, long, help = "Number of threads to use (default: 2)")]
    pub threads: Option<usize>,

    #[structopt(long, default_value = "5", help = "Minimum coverage threshold")]
    pub min_coverage: u16,

    #[structopt(long, short = "v", help = "Verbose output")]
    pub verbose: bool,

    #[structopt(long, help = "Dry run - validate inputs without processing")]
    pub dry_run: bool,
}

impl CallArgs {
    /// Validate command line arguments and input files
    ///
    /// This function performs several validation steps:
    /// 1. Validates configuration parameters
    /// 2. Validates input files exist and are readable
    /// 3. Creates output directory if it doesn't exist
    /// 4. Validates output path is writable
    pub fn validate(&self) -> Result<(), RedicatError> {
        // Validate configuration parameters
        let config = ValidationConfig {
            max_other_threshold: self.max_other_threshold,
            min_edited_threshold: self.min_edited_threshold,
            min_ref_threshold: self.min_ref_threshold,
            min_coverage: self.min_coverage,
            chunk_size: self.chunksize,
            num_threads: self.threads.unwrap_or_else(num_cpus::get),
        };

        config.validate()?;
        info!("Configuration validation passed");

        // Validate input files exist and are readable
        validate_input_files(&self.input, &self.fa, &self.site_white_list)?;

        // Create output directory if it doesn't exist
        if let Some(parent) = Path::new(&self.output).parent() {
            fs::create_dir_all(parent).map_err(|e| {
                RedicatError::InvalidInput(format!("Failed to create output directory: {}", e))
            })?;
        }

        validate_output_path(&self.output)?;
        info!("File validation passed");

        Ok(())
    }
}

/// Main entry point for the RNA editing analysis pipeline
///
/// This function orchestrates the entire REDICAT analysis workflow:
/// 1. Validates command line arguments
/// 2. Configures the thread pool for parallel processing
/// 3. Executes the main analysis pipeline
/// 4. Handles errors and logs results
///
/// # Arguments
///
/// * `args` - Command line arguments parsed into CallArgs struct
///
/// # Returns
///
/// * `Result<()>` - Ok if successful, Err with error details if failed
pub fn run_call(args: CallArgs) -> Result<()> {
    info!("Starting REDICAT analysis pipeline");
    info!("Arguments: {:?}", args);

    // Validate command line arguments
    args.validate()?;

    if args.dry_run {
        info!("Dry run completed successfully - all validations passed");
        return Ok(());
    }

    // Configure number of threads (default to 2 for stability)
    let num_threads = args.threads.unwrap_or_else(|| 2);
    // Alternative: let num_threads = args.threads.unwrap_or_else(num_cpus::get);

    // Set up thread pool for parallel processing
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .map_err(|e| RedicatError::Config(format!("Failed to build thread pool: {}", e)))?;

    info!("Thread pool configured with {} threads", num_threads);

    // Execute main analysis pipeline
    let result = run_analysis(&args);

    match result {
        Ok(_) => {
            info!("Analysis completed successfully!");
            info!("Output saved to: {}", args.output);
        }
        Err(e) => {
            error!("Analysis failed: {}", e);
            return Err(e.into());
        }
    }

    Ok(())
}

/// Execute the main analysis pipeline
///
/// This function performs the core RNA editing analysis steps with strand-aware processing:
/// 1. Reads the input AnnData file
/// 2. Loads reference genome and editing site data
/// 3. Annotates variants using the annotation pipeline
/// 4. Calculates reference/alternate matrices for the specified editing type using strand-aware logic
/// 5. Computes the Cell Editing Index (CEI)
/// 6. Calculates site-level mismatch statistics
/// 7. Writes the results to an output file
///
/// The analysis is performed in a strand-aware manner to correctly handle RNA editing
/// events on both positive and negative DNA strands. This ensures comprehensive
/// detection of editing events regardless of the DNA strand they originate from.
///
/// # Arguments
///
/// * `args` - Reference to CallArgs containing pipeline parameters
///
/// # Returns
///
/// * `Result<(), RedicatError>` - Ok if successful, Err with error details if failed
fn run_analysis(args: &CallArgs) -> Result<(), RedicatError> {
    info!("Reading input AnnData file: {}", args.input);
    let mut adata = read_anndata_h5ad(&args.input)
        .map_err(|e| RedicatError::DataProcessing(format!("Failed to read input file: {}", e)))?;

    info!(
        "Loaded AnnData with shape: {} × {}",
        adata.n_obs, adata.n_vars
    );

    info!("Loading reference genome: {}", args.fa);
    let reference = Arc::new(ReferenceGenome::new(&args.fa)?);

    info!("Loading site white list data: {}", args.site_white_list);
    let editing_sites = Arc::new(load_rediportal_parallel(&args.site_white_list)?);

    info!("Annotating variants...");
    adata = annotate_variants_pipeline(
        adata,
        editing_sites.clone(),
        reference.clone(),
        args.max_other_threshold as f32,
        args.min_edited_threshold as f32,
        args.min_ref_threshold as f32,
        args.min_coverage as u32,
    )?;

    info!(
        "Calculating ref/alt matrices for editing type: {:?} (strand-aware processing)",
        args.editingtype
    );
    adata = calculate_ref_alt_matrices(adata, &args.editingtype)?;

    info!("Calculating CEI...");
    adata = calculate_cei(adata)?;

    info!("Calculating site-level mismatch statistics...");
    let (ref_base, alt_base) = args.editingtype.to_bases();
    adata = calculate_site_mismatch_stats(adata, ref_base, alt_base)?;

    info!("Writing output: {}", args.output);
    write_anndata_h5ad(&adata, &args.output)
        .map_err(|e| RedicatError::DataProcessing(format!("Failed to write output file: {}", e)))?;

    Ok(())
}
