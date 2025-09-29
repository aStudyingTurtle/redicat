pub mod bam2mtx;
pub mod base_depth;
pub mod call;

/// Canonical human contigs processed by default for BAM-driven analytics.
///
/// The list follows UCSC naming (chr-prefixed autosomes, sex chromosomes, and mitochondrial DNA)
/// to match the majority of single-cell alignments. Use the `--allcontigs` flag in the CLI to opt-in
/// to alternative contigs, decoys, or spike-in references when required.
pub const STANDARD_CONTIGS: &[&str] = &[
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22", "chrX", "chrY", "chrM",
];

/// Returns `true` when a contig name matches one of the canonical autosomes, sex chromosomes,
/// or mitochondrial chromosome listed in [`STANDARD_CONTIGS`].
#[inline]
pub fn is_standard_contig(name: &str) -> bool {
    STANDARD_CONTIGS
        .iter()
        .any(|contig| contig.eq_ignore_ascii_case(name))
}

// Re-export Bulk struct with new name
pub use base_depth::Bulk;
