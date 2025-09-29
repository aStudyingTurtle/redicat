use anyhow::{anyhow, Result};
use std::path::{Path, PathBuf};

use crate::commands::bulk::{run_bulk, BulkArgs};
use redicat_lib::utils;

use super::args::Bam2MtxArgs;

pub fn prepare_positions_file(args: &Bam2MtxArgs) -> Result<PathBuf> {
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

    let bulk_args = BulkArgs {
        reads: args.bam.clone(),
        output: target.to_path_buf(),
        threads: args.threads,
        chunksize: args.chunksize,
        min_baseq: Some(args.min_baseq),
        mapquality: args.min_mapq,
        zero_base: false,
        max_depth: 8000,
        min_depth: args.min_depth,
        max_n_fraction: args.max_n_fraction,
        all: false,
        editing_threshold: args.editing_threshold,
        all_contigs: args.all_contigs,
    };

    run_bulk(bulk_args)
}
