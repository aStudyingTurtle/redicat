use anyhow::Result;
use log::info;
use num_cpus;
use redicat_lib::call::editing::EditingType;
use redicat_lib::call::validation::{validate_input_files, validate_output_path, ValidationConfig};
use redicat_lib::call::RedicatError;
use std::fs;
use std::path::Path;
use structopt::StructOpt;

#[derive(StructOpt, Debug, Clone)]
#[structopt(
    name = "call",
    about = "REDICAT analysis pipeline - Rust implementation"
)]
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
    pub fn validate(&self) -> Result<(), RedicatError> {
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

        validate_input_files(&self.input, &self.fa, &self.site_white_list)?;

        if let Some(parent) = Path::new(&self.output).parent() {
            fs::create_dir_all(parent).map_err(|e| {
                RedicatError::InvalidInput(format!("Failed to create output directory: {}", e))
            })?;
        }

        validate_output_path(&self.output)?;
        info!("File validation passed");

        Ok(())
    }

    #[inline]
    pub fn effective_threads(&self) -> usize {
        self.threads.unwrap_or(2)
    }
}
