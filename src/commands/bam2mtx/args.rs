use std::path::PathBuf;
use structopt::StructOpt;

/// Arguments for the `bam2mtx` command.
#[derive(Debug, Clone, StructOpt)]
#[structopt(name = "bam2mtx", about = "Convert BAM files to single-cell matrices")]
pub struct Bam2MtxArgs {
    /// Path to the input BAM file (must be indexed).
    #[structopt(short, long, parse(from_os_str))]
    pub bam: PathBuf,

    /// Optional TSV file with genomic positions (CHR/POS). Can be auto-generated via `--two-pass`.
    #[structopt(long, parse(from_os_str))]
    pub tsv: Option<PathBuf>,

    /// Path to the cell barcode whitelist file.
    #[structopt(long, parse(from_os_str))]
    pub barcodes: PathBuf,

    /// Enable the two-pass workflow: run `bulk` to build a site list before matrix generation.
    #[structopt(long)]
    pub two_pass: bool,

    /// Output path for the `.h5ad` file.
    #[structopt(short, long, parse(from_os_str))]
    pub output: PathBuf,

    /// Number of threads to use (default: 10).
    #[structopt(short, long, default_value = "10")]
    pub threads: usize,

    /// Minimum mapping quality.
    #[structopt(long, default_value = "255", short = "q")]
    pub min_mapq: u8,

    /// Minimum base quality.
    #[structopt(long, default_value = "30", short = "Q")]
    pub min_baseq: u8,

    /// Minimum effective depth (excluding Ns) required to keep a site.
    #[structopt(long = "min-depth", default_value = "10", short = "d")]
    pub min_depth: u32,

    /// Maximum allowed fraction of ambiguous bases (denominator form, like `bulk`).
    #[structopt(long = "max-n-fraction", default_value = "20", short = "n")]
    pub max_n_fraction: u32,

    /// Editing threshold used to detect multi-base support.
    #[structopt(long = "editing-threshold", default_value = "1000", short = "et")]
    pub editing_threshold: u32,

    /// Whether the library is stranded (default: unstranded).
    #[structopt(long, short = "S")]
    pub stranded: bool,

    /// Maximum pileup depth to examine per site.
    #[structopt(long = "max-depth", default_value = "655360", short = "D")]
    pub max_depth: u32,

    /// UMI tag name.
    #[structopt(long, default_value = "UB")]
    pub umi_tag: String,

    /// Cell barcode tag name.
    #[structopt(long, default_value = "CB")]
    pub cb_tag: String,

    /// Path to reference FASTA file (required for CRAM files).
    #[structopt(long, parse(from_os_str), short = "r")]
    pub reference: Option<PathBuf>,

    /// Chunk size for parallel processing (weight budget for normal positions). Default: 100_000
    #[structopt(long, default_value = "100000", short = "c")]
    pub chunksize: u32,

    /// Chunk size applied to high-depth loci marked by `NEAR_MAX_DEPTH` in the TSV.
    /// Default: 2
    #[structopt(long = "chunk-size-max-depth", default_value = "2")]
    pub chunk_size_max_depth: u32,

    /// Matrix density estimate used to pre-size sparse buffers.
    #[structopt(long, default_value = "0.005")]
    pub matrix_density: f64,

    /// Include all contigs from the BAM header instead of restricting to canonical chromosomes.
    #[structopt(long = "allcontigs", short = "A")]
    pub all_contigs: bool,
}

impl Bam2MtxArgs {
    #[inline]
    pub fn chunk_size(&self) -> usize {
        usize::max(self.chunksize as usize, 1)
    }

    #[inline]
    pub fn chunk_size_max_depth(&self) -> usize {
        usize::max(self.chunk_size_max_depth as usize, 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn parses_minimal_arguments() {
        let args = Bam2MtxArgs::from_iter_safe(&[
            "bam2mtx",
            "--bam",
            "test.bam",
            "--tsv",
            "test.tsv",
            "--barcodes",
            "barcodes.tsv",
            "--output",
            "output.h5ad",
        ])
        .unwrap();

        assert_eq!(args.bam, PathBuf::from("test.bam"));
        assert_eq!(args.tsv, Some(PathBuf::from("test.tsv")));
        assert_eq!(args.barcodes, PathBuf::from("barcodes.tsv"));
        assert_eq!(args.output, PathBuf::from("output.h5ad"));
        assert_eq!(args.max_depth, 655_360);
        assert_eq!(args.chunk_size_max_depth, 2);
        assert!(!args.two_pass);
    }
}
