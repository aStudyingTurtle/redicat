# REdiCAT API Reference

This document provides detailed documentation for the REdiCAT Rust library API.

## Library Structure

The REdiCAT library is organized into several modules:

- `call`: RNA editing analysis pipeline
- `bam2mtx`: BAM to matrix conversion
- `par_granges`: Parallel genomic region processing
- `position`: Genomic position data structures
- `read_filter`: Read filtering functionality
- `utils`: Utility functions

## call Module

The `call` module is the core of the RNA editing analysis functionality.

### Key Components

1. `anndata_ops`: Reading and writing AnnData files
2. `base_matrix`: Working with base count matrices
3. `editing_analysis`: Core RNA editing analysis pipeline
4. `error`: Error types specific to the call module
5. `reference_genome`: Reference genome operations
6. `sparse_ops`: Sparse matrix operations
7. `validation`: Input validation functions
8. `editing`: Editing type definitions and related functionality

### Main Types and Functions

- `AnnDataContainer`: Main data structure for holding AnnData information
- `ReferenceGenome`: Thread-safe reference genome access
- `EditingType`: Enum representing different RNA editing types
- `annotate_variants_pipeline`: Main variant annotation pipeline
- `calculate_ref_alt_matrices`: Calculate reference/alternate matrices
- `calculate_cei`: Calculate Cell Editing Index
- `calculate_site_mismatch_stats`: Analyze site-level mismatch patterns

### anndata_ops

Handles reading and writing of AnnData files in H5AD format.

#### AnnDataContainer

Main data structure holding all AnnData information:

```rust
pub struct AnnDataContainer {
    pub obs: DataFrame,                     // Cell-level annotations
    pub var: DataFrame,                     // Site-level annotations
    pub x: Option<CsrMatrix<f64>>,          // Main data matrix
    pub layers: HashMap<String, CsrMatrix<u32>>, // Additional data matrices (u32 sparse matrices)
    pub n_obs: usize,                       // Number of observations (cells)
    pub n_vars: usize,                      // Number of variables (sites)
    pub var_names: Vec<String>,             // Variable names
    pub obs_names: Vec<String>,             // Observation names
}
```

#### Functions

- `read_anndata_h5ad(path: &str) -> Result<AnnDataContainer>`: Read AnnData from H5AD file
- `write_anndata_h5ad(adata: &AnnDataContainer, path: &str) -> Result<()>`: Write AnnData to H5AD file
- `compute_layer_row_sums(&self, layer_name: &str) -> Option<Vec<u32>>`: Compute row sums for a layer
- `compute_layer_col_sums(&self, layer_name: &str) -> Option<Vec<u32>>`: Compute column sums for a layer
- `compute_total_coverage(&self) -> Vec<u32>`: Get total coverage across all base layers

### base_matrix

Operations for working with base count matrices.

#### Functions

- `calculate_coverage(adata: &mut AnnDataContainer) -> Result<()>`: Calculate coverage matrix
- `filter_sites_by_coverage(adata: AnnDataContainer, min_coverage: u32) -> Result<AnnDataContainer>`: Filter base matrix by coverage
- `count_base_levels(adata: AnnDataContainer) -> Result<AnnDataContainer>`: Count base levels at each site
- `apply_site_filter(adata: AnnDataContainer, filter_mask: &[bool]) -> Result<AnnDataContainer>`: Apply site filter

### editing_analysis

Core RNA editing analysis pipeline functions.

#### Functions

- `annotate_variants_pipeline(...)`: Main variant annotation pipeline
- `calculate_ref_alt_matrices(adata: AnnDataContainer, editing_type: &EditingType) -> Result<AnnDataContainer>`: Calculate ref/alt matrices
- `calculate_cei(adata: AnnDataContainer) -> Result<AnnDataContainer>`: Calculate editing index
- `calculate_site_mismatch_stats(adata: AnnDataContainer, ref_base: char, alt_base: char) -> Result<AnnDataContainer>`: Analyze site-level mismatches

### error

Error types for the call module.

#### RedicatError

Enumeration of all possible errors in the call module:

```rust
pub enum RedicatError {
    Io(std::io::Error),
    FileNotFound(String),
    DataProcessing(String),
    InvalidInput(String),
    Parse(String),
    DimensionMismatch { expected: String, actual: String },
    SparseMatrix(String),
    EmptyData(String),
    ReferenceGenome(String),
    Config(String),
}
```

### reference_genome

Reference genome operations.

#### ReferenceGenome

Thread-safe reference genome accessor:

```rust
pub struct ReferenceGenome {
    reader: parking_lot::Mutex<IndexedReader<std::fs::File>>,
    sequences: Vec<String>,
}
```

#### Functions

- `new(fasta_path: &str) -> Result<Self>`: Create new ReferenceGenome
- `get_ref_of_pos(genomic_pos: &str) -> Result<char>`: Get reference base at position
- `get_multiple_refs_chunked(positions: &[String], chunk_size: usize) -> Result<Vec<char>>`: Get multiple reference bases
- `validate_position(genomic_pos: &str) -> bool`: Validate position

### sparse_ops

Sparse matrix operations with performance optimizations.

#### SparseOps

Utility struct for sparse matrix operations with optimizations for large-scale data:

```rust
pub struct SparseOps;
```

The sparse matrix operations are optimized for performance with large single-cell datasets:
- Uses native `nalgebra_sparse` operations for maximum efficiency
- Implements parallel processing for computationally intensive operations
- Employs efficient memory management strategies to handle large matrices
- Uses optimized algorithms for common operations like matrix addition and filtering

#### Functions

- `from_triplets_u32(nrows: usize, ncols: usize, triplets: Vec<(usize, usize, u32)>) -> Result<CsrMatrix<u32>>`: Create CSR matrix from u32 triplets
- `add_matrices(a: &CsrMatrix<u32>, b: &CsrMatrix<u32>) -> Result<CsrMatrix<u32>>`: Add two u32 CSR matrices
- `filter_columns_u32(matrix: &CsrMatrix<u32>, keep_indices: &[usize]) -> Result<CsrMatrix<u32>>`: Filter matrix columns for u32 matrices
- `compute_row_sums(matrix: &CsrMatrix<u32>) -> Vec<u32>`: Compute row sums for u32 matrix
- `compute_col_sums(matrix: &CsrMatrix<u32>) -> Vec<u32>`: Compute column sums for u32 matrix

### validation

Input validation functions.

#### Functions

- `validate_input_files(input: &str, fa: &str, site_white_list: &str) -> Result<()>`: Validate input files
- `validate_output_path(output: &str) -> Result<()>`: Validate output path
- `validate(&self) -> Result<()>`: Validate configuration parameters
- `validate_matrix_dimensions(matrix: &CsrMatrix<u32>, expected_rows: usize, expected_cols: usize) -> Result<()>`: Validate matrix dimensions

### editing

Editing type definitions and related functionality with strand-aware logic.

#### EditingType

Enumeration of RNA editing types with strand-aware processing:

```rust
pub enum EditingType {
    AG,  // A to G editing
    AC,  // A to C editing
    AT,  // A to T editing
    CA,  // C to A editing
    CG,  // C to G editing
    CT,  // C to T editing
}
```

The strand-aware logic properly handles editing events on both positive and negative strands by:
- Determining valid reference bases for each strand
- Computing the corresponding alternate bases considering strand orientation
- Ensuring proper complementarity for negative strand editing events

#### Functions

- `load_rediportal_parallel(path: &str) -> Result<HashMap<String, u8>>`: Load REDIPortal editing sites
- `to_bases(&self) -> (char, char)`: Get the reference and alternate bases for this editing type
- `get_strand_aware_ref_bases(&self) -> [char; 2]`: Get allowed reference bases for both strands
- `get_alt_base_for_ref(&self, ref_base: char) -> char`: Get the corresponding alt base for a given ref base considering strand

## bam2mtx Module

BAM to matrix conversion functionality.

### Key Components

1. `anndata_output`: AnnData output functionality
2. `barcode`: Cell barcode processing
3. `processor`: BAM file processing
4. `region_processor`: Region processing for par_granges
5. `utils`: Utility functions

### anndata_output

AnnData output functionality with performance optimizations.

#### AnnDataConfig

Configuration for AnnData output:

```rust
pub struct AnnDataConfig {
    pub stranded: bool,
    pub compression: Option<String>,
    pub threads: usize,
    pub chunk_size: usize,
    pub matrix_density: f64,
    pub use_mmap: bool,
    pub batch_size: usize,
}
```

#### AnnDataConverter

Converter for transforming processed position data into AnnData format:

```rust
pub struct AnnDataConverter {
    config: AnnDataConfig,
}
```

#### Functions

- `new(config: AnnDataConfig) -> Self`: Create new converter
- `convert(position_data: &[PositionData], output_path: &Path) -> Result<AnnData<H5>>`: Convert to AnnData
- `write_to_file(adata: &AnnData<H5>, output_path: &Path) -> Result<()>`: Write to file

### barcode

Cell barcode processing functionality.

#### BarcodeProcessor

Processor for cell barcodes:

```rust
pub struct BarcodeProcessor {
    valid_barcodes: Arc<HashSet<String>>,
}
```

#### Functions

- `from_file<P: AsRef<Path>>(path: P) -> Result<Self>`: Create from file
- `is_valid(&self, barcode: &str) -> bool`: Check if barcode is valid
- `len(&self) -> usize`: Get number of valid barcodes

### processor

BAM file processing for single-cell data.

#### BaseCounts

Base counts for a specific position:

```rust
pub struct BaseCounts {
    pub a: u32,
    pub t: u32,
    pub g: u32,
    pub c: u32,
}
```

#### StrandBaseCounts

Strand-specific base counts:

```rust
pub struct StrandBaseCounts {
    pub forward: BaseCounts,
    pub reverse: BaseCounts,
}
```

#### PositionData

Processed data for a specific genomic position:

```rust
pub struct PositionData {
    pub chrom: String,
    pub pos: u64,
    pub counts: HashMap<String, StrandBaseCounts>,
}
```

#### BamProcessorConfig

Configuration for BAM processing:

```rust
pub struct BamProcessorConfig {
    pub min_mapping_quality: u8,
    pub min_base_quality: u8,
    pub stranded: bool,
    pub umi_tag: String,
    pub cell_barcode_tag: String,
}
```

#### BamProcessor

Main processor for BAM files:

```rust
pub struct BamProcessor {
    config: BamProcessorConfig,
    barcode_processor: Arc<BarcodeProcessor>,
}
```

#### Functions

- `new(config: BamProcessorConfig, barcode_processor: Arc<BarcodeProcessor>) -> Self`: Create new processor
- `process_position(&self, bam_path: &Path, chrom: &str, pos: u64) -> Result<PositionData>`: Process single position
- `process_positions_parallel(...) -> Result<Vec<PositionData>>`: Process positions in parallel
- `process_positions_optimized(...) -> Result<Vec<PositionData>>`: Optimized position processing

## par_granges Module

Parallel processing of genomic regions.

### RegionProcessor

Trait defining methods for processing genomic regions:

```rust
pub trait RegionProcessor {
    type P: 'static + Send + Sync + Serialize;
    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P>;
}
```

### ParGranges

Main struct for parallel genomic region processing:

```rust
pub struct ParGranges<R: 'static + RegionProcessor + Send + Sync> {
    reads: PathBuf,
    ref_fasta: Option<PathBuf>,
    regions_bed: Option<PathBuf>,
    regions_bcf: Option<PathBuf>,
    merge_regions: bool,
    threads: usize,
    chunksize: u32,
    channel_size_modifier: f64,
    pool: rayon::ThreadPool,
    processor: R,
}
```

#### Functions

- `new(...) -> Self`: Create new ParGranges
- `process(self) -> Result<Receiver<R::P>>`: Process regions

## position Module

Data structures for genomic positions.

### PileupPosition

Implementation of `Position` for dealing with pileups:

```rust
pub struct PileupPosition {
    pub ref_seq: SmartString<LazyCompact>,
    pub pos: u32,
    pub ref_base: Option<char>,
    pub depth: u32,
    pub a: u32,
    pub c: u32,
    pub g: u32,
    pub t: u32,
    pub n: u32,
    pub ins: u32,
    pub del: u32,
    pub ref_skip: u32,
    pub fail: u32,
    pub near_max_depth: bool,
}
```

#### Functions

- `from_pileup(...) -> Self`: Convert pileup to PileupPosition
- `compact_refseq(header: &HeaderView, tid: u32) -> SmartString<LazyCompact>`: Convert tid to reference sequence name

## read_filter Module

Read filtering functionality.

### ReadFilter

Trait for filtering reads based on various criteria:

```rust
pub trait ReadFilter {
    fn filter_read(&self, read: &Record, alignment: Option<&Alignment>) -> bool;
}
```

### DefaultReadFilter

Straightforward read filter based on SAM flags and mapping quality:

```rust
pub struct DefaultReadFilter {
    include_flags: u16,
    exclude_flags: u16,
    min_mapq: u8,
}
```

#### Functions

- `new(include_flags: u16, exclude_flags: u16, min_mapq: u8) -> Self`: Create new filter
- `filter_read(&self, read: &Record, _alignment: Option<&Alignment>) -> bool`: Filter read

## utils Module

Utility functions used throughout the library.

### Functions

- `get_writer(...) -> Result<Writer>`: Get appropriate writer for output
- `determine_allowed_cpus(requested: usize) -> Result<usize>`: Determine allowed CPUs
- `is_broken_pipe(err: &anyhow::Error) -> bool`: Check if error is broken pipe