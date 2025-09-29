//! REDICAT - RNA Editing Cellular Assessment Toolkit
//!
//! REDICAT (RNA Editing Cellular Assessment Toolkit) is a highly parallelized utility
//! for analyzing RNA editing events in single-cell RNA-seq data. Originally designed for
//! detecting indels in reduced representation sequencing data, REDICAT has been extended
//! to include powerful functionality for comprehensive RNA editing analysis.
//!
//! # Tools
//!
//! REDICAT provides several subcommands for different analysis types:
//!
//! - `bulk`: Calculate depth and nucleotide counts at each base position
//! - `bam2mtx`: Convert BAM files to single-cell matrices
//! - `call`: RNA editing detection and analysis pipeline
//!
//! # Usage
//!
//! ```bash
//! # Analyze base depth at each position
//! redicat bulk input.bam
//!
//! # Convert BAM to single-cell matrix
//! redicat bam2mtx --bam input.bam --tsv positions.tsv --barcodes barcodes.tsv --output matrix.h5ad
//!
//! # Run RNA editing analysis pipeline
//! redicat call --input input.h5ad --output output.h5ad --fa reference.fa --site-white-list editing_sites.tsv.gz
//!
//! # Run bam2mtx with automatic site discovery
//! redicat bam2mtx --bam input.bam --barcodes barcodes.tsv --output matrix.h5ad --two-pass
//! ```
//!
//! For more detailed usage information, see the documentation for each subcommand.

extern crate redicat_lib;
pub mod commands;
use anyhow::Result;
use env_logger::Env;
use log::*;
use redicat_lib::utils;
use structopt::StructOpt;

#[derive(StructOpt)]
#[structopt(rename_all = "kebab-case", author, about)]
/// Commands for generating per-base analysis with REDICAT
struct Args {
    #[structopt(subcommand)]
    subcommand: Subcommand,
}

#[derive(StructOpt)]
enum Subcommand {
    /// Calculate the depth at each base, per-nucleotide
    Bulk(commands::BulkArgs),
    /// Convert BAM files to single-cell matrices
    Bam2mtx(commands::Bam2MtxArgs),
    /// RNA editing detection and analysis pipeline
    Call(commands::CallArgs),
}

impl Subcommand {
    fn run(self) -> Result<()> {
        match self {
            Subcommand::Bulk(args) => commands::run_bulk(args)?,
            Subcommand::Bam2mtx(args) => commands::run_bam2mtx(args)?,
            Subcommand::Call(args) => commands::run_call(args)?,
        }
        Ok(())
    }
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    if let Err(err) = Args::from_args().subcommand.run() {
        if utils::is_broken_pipe(&err) {
            std::process::exit(0);
        }
        error!("{}", err);
        std::process::exit(1);
    }
    Ok(())
}
