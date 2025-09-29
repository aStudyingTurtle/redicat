use redicat_lib::engine::par_granges;
use std::path::PathBuf;
use structopt::StructOpt;

/// CLI arguments for the `bulk` subcommand.
#[derive(Debug, Clone, StructOpt)]
#[structopt(author, name = "bulk")]
pub struct BulkArgs {
    /// Input indexed BAM/CRAM to analyze.
    pub reads: PathBuf,

    /// Output path (required, `.gz` will be appended when missing).
    #[structopt(long, short = "o")]
    pub output: PathBuf,

    /// Number of worker threads to use.
    #[structopt(long, short = "t", default_value = "10")]
    pub threads: usize,

    /// Ideal base pairs per worker. Total bp resident ≈ (`threads` * `chunksize`).
    #[structopt(long, short = "c", default_value = par_granges::CHUNKSIZE_STR.as_str())]
    pub chunksize: u32,

    /// Minimum base quality to treat a base as [A, C, G, T]. Lower bases count as `N`.
    #[structopt(long, short = "Q")]
    pub min_baseq: Option<u8>,

    /// Minimum mapping quality for reads to be counted.
    #[structopt(long, default_value = "255", short = "q")]
    pub mapquality: u8,

    /// Output positions as 0-based instead of 1-based.
    #[structopt(long, short = "z")]
    pub zero_base: bool,

    /// Maximum pileup depth inspected per site. Positions within 1% of the cap are flagged.
    #[structopt(long, short = "D", default_value = "8000")]
    pub max_depth: u32,

    /// Minimum non-`N` depth required to report a site.
    #[structopt(long, short = "d", default_value = "10")]
    pub min_depth: u32,

    /// Maximum tolerated ambiguous depth fraction (depth / value).
    #[structopt(long, short = "n", default_value = "20")]
    pub max_n_fraction: u32,

    /// Report every position instead of the editing-enriched subset.
    #[structopt(long, short = "a")]
    pub all: bool,

    /// Editing threshold for valid value calculation (`depth / editing_threshold`).
    #[structopt(long = "editing-threshold", short = "et", default_value = "1000")]
    pub editing_threshold: u32,

    /// Visit every contig in the BAM header instead of the canonical subset.
    #[structopt(long = "allcontigs", short = "A")]
    pub all_contigs: bool,
}

/// Normalised configuration derived from [`BulkArgs`].
#[derive(Debug, Clone)]
pub struct BulkConfig {
    pub reads: PathBuf,
    pub output: PathBuf,
    pub threads: usize,
    pub chunksize: u32,
    pub min_baseq: Option<u8>,
    pub mapquality: u8,
    pub coord_offset: u32,
    pub max_depth: u32,
    pub min_depth: u32,
    pub max_n_fraction: u32,
    pub report_all: bool,
    pub editing_threshold: u32,
    pub all_contigs: bool,
}

impl From<BulkArgs> for BulkConfig {
    fn from(args: BulkArgs) -> BulkConfig {
        BulkConfig {
            reads: args.reads,
            output: args.output,
            threads: args.threads,
            chunksize: args.chunksize,
            min_baseq: args.min_baseq,
            mapquality: args.mapquality,
            coord_offset: if args.zero_base { 0 } else { 1 },
            max_depth: args.max_depth,
            min_depth: args.min_depth,
            max_n_fraction: args.max_n_fraction,
            report_all: args.all,
            editing_threshold: args.editing_threshold,
            all_contigs: args.all_contigs,
        }
    }
}

impl BulkConfig {
    #[inline]
    pub fn editing_mode(&self) -> bool {
        !self.report_all
    }
}
