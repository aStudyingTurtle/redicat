# REdiCAT Examples

This document provides practical examples of how to use REdiCAT for various analysis tasks.

## Table of Contents

1. [RNA Editing Analysis](#rna-editing-analysis)
2. [Single-cell Matrix Generation](#single-cell-matrix-generation)
3. [Bulk Analysis](#bulk-analysis)
4. [Preprocessing](#preprocessing)
5. [Library Usage Examples](#library-usage-examples)

## RNA Editing Analysis

### Basic RNA Editing Pipeline

This example shows how to run the complete RNA editing analysis pipeline:

```bash
redicat call \
  --input input.h5ad \
  --output output.h5ad \
  --fa reference.fa \
  --site-white-list editing_sites.tsv.gz \
  --editingtype ag \
  --num-threads 16 \
  --min-coverage 10
```

Note: The default number of threads is 2. In this example, we explicitly set it to 16 for better performance on multi-core systems. The analysis is performed in a strand-aware manner, correctly handling complementary bases on the negative strand.

### Advanced RNA Editing Analysis

For more detailed analysis with custom thresholds:

```bash
redicat call \
  --input input.h5ad \
  --output output.h5ad \
  --fa reference.fa \
  --site-white-list editing_sites.tsv.gz \
  --editingtype ag \
  --max-other-threshold 0.05 \
  --min-edited-threshold 0.02 \
  --min-ref-threshold 0.02 \
  --min-coverage 5 \
  --num-threads 32
```

This analysis uses strand-aware logic to properly handle editing events on both positive and negative strands. The strand-aware processing considers that the same editing event can be observed from either DNA strand due to complementary base pairing:
- A>G editing on the positive strand appears as T>C editing on the negative strand
- A>C editing on the positive strand appears as T>G editing on the negative strand
- And so on for all editing types

This ensures comprehensive detection of RNA editing events regardless of the DNA strand they originate from.

## Single-cell Matrix Generation

### Basic BAM to Matrix Conversion

Generate a count matrix from BAM files:

```bash
redicat bam2mtx \
  --bam input.bam \
  --tsv positions.tsv \
  --barcodes barcodes.tsv \
  --output output.h5ad \
  --num-threads 8
```

Note: The default number of threads is 10. In this example, we explicitly set it to 8.

### Stranded BAM to Matrix Conversion

For stranded data:

```bash
redicat bam2mtx \
  --bam input.bam \
  --tsv positions.tsv \
  --barcodes barcodes.tsv \
  --output output.h5ad \
  --stranded \
  --num-threads 8
```

### Large Dataset Optimization

For large datasets (>10M cells or positions):

```bash
redicat bam2mtx \
  --bam input.bam \
  --tsv positions.tsv \
  --barcodes barcodes.tsv \
  --output output.h5ad \
  --large-dataset \
  --num-threads 16
```

### Memory-Efficient Processing

For resource-constrained environments:

```bash
redicat bam2mtx \
  --bam input.bam \
  --tsv positions.tsv \
  --barcodes barcodes.tsv \
  --output output.h5ad \
  --memory-efficient \
  --num-threads 4
```

## Bulk Analysis

### Basic Base Depth Analysis

Analyze base depth across the entire genome:

```bash
redicat bulk input.bam > output.tsv
```

Note: Uses default parameters including 10 threads and minimum depth of 10.

### Region-specific Analysis

Analyze only specific genomic regions:

```bash
redicat bulk \
  --bed-file regions.bed \
  --min-mapq 30 \
  --exclude-flags 3848 \
  input.bam > output.tsv
```

### Output with Reference Bases

Include reference bases in the output:

```bash
redicat bulk \
  --ref-fasta reference.fa \
  input.bam > output.tsv
```

### High-depth Analysis

For high-depth sequencing data:

```bash
redicat bulk \
  --max-depth 50000 \
  --min-depth 100 \
  input.bam > output.tsv
```

### Edited Mode Analysis

Use edited filtering logic for more stringent filtering:

```bash
redicat bulk \
  --edited \
  --min-base-quality-score 20 \
  input.bam > output.tsv
```

## Preprocessing

### Basic BAM Filtering

Filter BAM files based on quality criteria:

```bash
redicat preprocess \
  --barcodes whitelist.tsv \
  --inbam input.bam \
  --outbam filtered.bam
```

### Advanced BAM Filtering

Filter with custom parameters:

```bash
redicat preprocess \
  --barcodes whitelist.tsv \
  --inbam input.bam \
  --outbam filtered.bam \
  --mapquality 30 \
  --exclude-flags 1796 \
  --include-flags 2 \
  --threads 16
```

Note: The default number of threads is 8. In this example, we explicitly set it to 16 for better performance on multi-core systems.

## Library Usage Examples

### Reading and Writing AnnData Files

```rust
use redicat_lib::call::anndata_ops::{read_anndata_h5ad, write_anndata_h5ad};

fn process_anndata() -> Result<(), Box<dyn std::error::Error>> {
    // Read an AnnData file
    let mut adata = read_anndata_h5ad("input.h5ad")?;
    
    // Process the data (example: add a new column)
    // ... your processing code here ...
    
    // Write the processed data back to disk
    write_anndata_h5ad(&adata, "output.h5ad")?;
    
    Ok(())
}
```

### Working with Reference Genomes

```rust
use redicat_lib::call::reference_genome::ReferenceGenome;

fn get_reference_bases() -> Result<(), Box<dyn std::error::Error>> {
    // Load reference genome (requires .fai index file)
    let reference = ReferenceGenome::new("reference.fa")?;
    
    // Get reference base at specific positions
    let positions = vec!["chr1:1000", "chr1:2000", "chr2:1500"];
    for pos in positions {
        let base = reference.get_ref_of_pos(pos)?;
        println!("Reference base at {}: {}", pos, base);
    }
    
    Ok(())
}
```

### RNA Editing Type Analysis

```rust
use redicat_lib::call::editing::EditingType;
use redicat_lib::call::editing_analysis::calculate_ref_alt_matrices;

fn analyze_editing_types() -> Result<(), Box<dyn std::error::Error>> {
    // Read input data
    let adata = read_anndata_h5ad("input.h5ad")?;
    
    // Analyze different editing types
    let editing_types = vec![EditingType::AG, EditingType::CT, EditingType::AC];
    
    for editing_type in editing_types {
        let result = calculate_ref_alt_matrices(adata.clone(), &editing_type)?;
        let output_file = format!("output_{}.h5ad", editing_type);
        write_anndata_h5ad(&result, &output_file)?;
        println!("Processed editing type: {}", editing_type);
    }
    
    Ok(())
}
```

### Custom Analysis Pipeline

```rust
use redicat_lib::call::editing::EditingType;
use redicat_lib::call::editing::load_rediportal_parallel;
use redicat_lib::call::reference_genome::ReferenceGenome;
use redicat_lib::call::editing_analysis::*;
use std::sync::Arc;

fn custom_analysis_pipeline() -> Result<(), Box<dyn std::error::Error>> {
    // Load input data
    let mut adata = read_anndata_h5ad("input.h5ad")?;
    
    // Load reference data
    let reference = Arc::new(ReferenceGenome::new("reference.fa")?);
    let editing_sites = Arc::new(load_rediportal_parallel("editing_sites.tsv.gz")?);
    
    // Run variant annotation pipeline
    adata = annotate_variants_pipeline(
        adata,
        editing_sites,
        reference,
        0.01,  // max_other_threshold
        0.01,  // min_edited_threshold
        0.01,  // min_ref_threshold
        5,     // min_coverage
    )?;
    
    // Process for A-to-G editing
    let editing_type = EditingType::AG;
    adata = calculate_ref_alt_matrices(adata, &editing_type)?;
    
    // Calculate editing index
    adata = calculate_cei(adata)?;
    
    // Save results
    write_anndata_h5ad(&adata, "output.h5ad")?;
    
    Ok(())
}
```

### Parallel Processing Example

```rust
use redicat_lib::call::editing_analysis::calculate_site_mismatch_stats;
use rayon::prelude::*;

fn parallel_mismatch_analysis() -> Result<(), Box<dyn std::error::Error>> {
    let adata = read_anndata_h5ad("input.h5ad")?;
    
    // Define editing types to analyze in parallel
    let editing_types = vec![('A', 'G'), ('C', 'T'), ('A', 'C')];
    
    // Process all editing types in parallel
    let results: Vec<_> = editing_types
        .par_iter()
        .map(|(ref_base, alt_base)| {
            calculate_site_mismatch_stats(
                adata.clone(),
                *ref_base,
                *alt_base
            )
        })
        .collect();
    
    // Handle results
    for (i, result) in results.into_iter().enumerate() {
        let (ref_base, alt_base) = editing_types[i];
        let filename = format!("mismatch_{}{}.h5ad", ref_base, alt_base);
        write_anndata_h5ad(&result?, &filename)?;
    }
    
    Ok(())
}
```

## Input File Formats

### Position Files (.tsv)

Position files should be tab-separated with at least two columns:

```
# Chromosome and position (1-based)
chr1	1000
chr1	1050
chr2	2000
chrX	3000
```

### Barcode Files (.tsv)

Barcode files should contain one barcode per line:

```
AAACCTGAGAAGGCAC
AAACCTGAGATAGGAT
AAACCTGCATGCATGC
AAACCTGGTAGAAGGA
```

### REDIPortal Files

REDIPortal files should be tab-separated with chromosome in the first column and position in the second:

```
chromosome	position	...
chr1	1000	...
chr1	1050	...
chr2	2000	...
```

## Output File Formats

### Bulk Analysis Output

The bulk analysis produces a tab-separated file like [perbase](https://github.com/sstadick/perbase) with the following columns:

| Column | Description |
|--------|-------------|
| CHR | Reference sequence name |
| POS | Position on the reference sequence |
| REF_BASE | Reference base at the position |
| DEPTH | Total depth at the position |
| A,C,G,T,N | Counts for each nucleotide |
| INS | Insertion counts |
| DEL | Deletion counts |
| REF_SKIP | Reference skip counts |
| FAIL | Failed read counts |
| NEAR_MAX_DEPTH | Flag for positions near max depth |

### AnnData Output

The RNA editing analysis produces AnnData files (.h5ad) containing:

1. **obs**: Cell-level annotations and metrics
2. **var**: Site-level annotations and editing information
3. **layers**: 
   - `ref`: Reference allele counts
   - `alt`: Alternate allele counts
   - `others`: Other nucleotide counts
   - Base-specific layers (A1, T1, G1, C1, etc.)
4. **Additional metrics**:
   - CEI (Cell Editing Index)
   - Site-level mismatch counts
   - Filtering information

## Performance Tips

1. **Use appropriate thread counts**: Set `--num-threads` based on your system's CPU cores
2. **Adjust chunk sizes**: For large datasets, tune `--chunk-size` parameter
3. **Filter low-quality data**: Use appropriate thresholds for MAPQ and coverage
4. **Use compressed inputs**: Gzipped input files can reduce I/O time
5. **Parallel processing**: Split large analyses into smaller parallel jobs when possible
6. **Memory management**: Use `--memory-efficient` flag for resource-constrained environments
7. **Large dataset optimization**: Use `--large-dataset` flag for datasets with >10M cells or positions

## Troubleshooting

### Common Issues

1. **Missing reference index**: Ensure your FASTA file has a corresponding `.fai` index
   ```bash
   samtools faidx reference.fa
   ```

2. **Memory issues**: Reduce `--num-threads` or `--chunk-size` for large datasets

3. **File not found errors**: Verify all input file paths are correct

4. **Permission errors**: Ensure you have read/write permissions for all files

5. **Invalid file formats**: Check that input files are in the correct format and not corrupted

### Debugging

Use the `--verbose` flag to get more detailed output:

```bash
redicat call --verbose --input input.h5ad --output output.h5ad --fa reference.fa --site-white-list editing_sites.tsv.gz
```

For dry runs to validate inputs without processing:

```bash
redicat call --dry-run --input input.h5ad --output output.h5ad --fa reference.fa --site-white-list editing_sites.tsv.gz
```

### Performance Profiling

To profile performance, you can use tools like `htop` or `perf`:

```bash
# Run with time to see execution time
time redicat call --input input.h5ad --output output.h5ad --fa reference.fa --site-white-list editing_sites.tsv.gz

# For more detailed profiling, use perf (Linux)
perf record -g redicat call --input input.h5ad --output output.h5ad --fa reference.fa --site-white-list editing_sites.tsv.gz
perf report
```