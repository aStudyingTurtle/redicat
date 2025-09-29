use anyhow::Result;
use log::{error, info};
use std::sync::Arc;

use crate::commands::common;
use redicat_lib::call::anndata_ops::{read_anndata_h5ad, write_anndata_h5ad};
use redicat_lib::call::annotate_variants_pipeline;
use redicat_lib::call::calculate_cei;
use redicat_lib::call::calculate_ref_alt_matrices;
use redicat_lib::call::calculate_site_mismatch_stats;
use redicat_lib::call::editing::load_rediportal_parallel;
use redicat_lib::call::reference_genome::ReferenceGenome;
use redicat_lib::call::RedicatError;

use super::args::CallArgs;

pub fn run_call(args: CallArgs) -> Result<()> {
    info!("Starting REDICAT analysis pipeline");
    info!("Arguments: {:?}", args);

    args.validate()?;

    if args.dry_run {
        info!("Dry run completed successfully - all validations passed");
        return Ok(());
    }

    let requested_threads = args.effective_threads();
    let active_threads = common::configure_global_thread_pool(requested_threads)?;
    info!(
        "Rayon thread pool configured with {} threads",
        active_threads
    );

    match execute_pipeline(&args) {
        Ok(_) => {
            info!("Analysis completed successfully");
            info!("Output written to: {}", args.output);
            Ok(())
        }
        Err(err) => {
            error!("Analysis failed: {}", err);
            Err(err.into())
        }
    }
}

fn execute_pipeline(args: &CallArgs) -> Result<(), RedicatError> {
    info!("Reading input AnnData file: {}", args.input);
    let mut adata = read_anndata_h5ad(&args.input)
        .map_err(|e| RedicatError::DataProcessing(format!("Failed to read input file: {}", e)))?;

    info!(
        "Loaded AnnData matrix with shape: {} × {}",
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
        "Calculating ref/alt matrices for editing type: {:?}",
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
