# REDICAT Documentation

REDICAT (RNA Editing Cellular Assessment Toolkit) is a comprehensive toolkit for analyzing RNA editing events in single-cell RNA-seq data.

## Table of Contents

- [REDICAT Documentation](#redicat-documentation)
  - [Table of Contents](#table-of-contents)
  - [Project Overview](#project-overview)
  - [Installation](#installation)
    - [Prerequisites](#prerequisites)
    - [Building from Source](#building-from-source)
  - [Quick Start](#quick-start)
    - [Basic Usage](#basic-usage)
  - [Tools](#tools)
    - [bulk](#bulk)
      - [Usage](#usage)
    - [bam2mtx](#bam2mtx)
      - [Usage](#usage-1)
    - [preprocess](#preprocess)
      - [Usage](#usage-2)
    - [call](#call)
      - [Usage](#usage-3)
  - [Library API](#library-api)
  - [Examples](#examples)
    - [RNA Editing Analysis Pipeline](#rna-editing-analysis-pipeline)
    - [Single-cell Matrix Generation](#single-cell-matrix-generation)
  - [Contributing](#contributing)

## Project Overview

REdiCAT is designed to provide a complete solution for RNA editing analysis in single-cell RNA-seq data. It combines efficient parallel processing with comprehensive analysis pipelines to identify and quantify RNA editing events at single-cell resolution.

Key features:
- High-performance parallel processing
- Support for various file formats (BAM/CRAM, H5AD, TSV)
- Comprehensive RNA editing analysis pipeline
- Single-cell matrix generation
- Quality control and filtering capabilities

## Installation

### Prerequisites

- Rust (latest stable version)
- Cargo
- samtools (for creating FASTA index files)

### Building from Source

```bash
# Clone the repository
git clone https://github.com/aStudyingTurtle/perbase_redicat.git
cd perbase_redicat

# Build the project
cargo build --release

# The executable will be located at target/release/redicat
```

## Quick Start

### Basic Usage

```bash
# Analyze base depth at each position
redicat bulk input.bam

# Convert BAM to single-cell matrix
redicat bam2mtx --bam input.bam --tsv positions.tsv --barcodes barcodes.tsv --output matrix.h5ad

# Filter BAM file
redicat preprocess --barcodes whitelist.tsv --inbam input.bam --outbam filtered.bam

# Run RNA editing analysis pipeline
redicat call --input input.h5ad --output output.h5ad --fa reference.fa --site-white-list editing_sites.tsv.gz
```

## Tools

### bulk

The `bulk` tool walks over every position in the BAM/CRAM file and calculates the depth, as well as the number of each nucleotide at the given position. Additionally, it counts the numbers of Ins/Dels at each position.

#### Usage

```bash
redicat bulk [FLAGS] [OPTIONS] <reads>
```

#### Default Parameters

- `--threads`: 10
- `--compression-threads`: 2
- `--compression-level`: 2
- `--chunksize`: 500000
- `--channel-size-modifier`: 0.5
- `--min-baseq`: 30
- `--zero-base`: false (1-based coordinates by default)
- `--max-depth`: 8000
- `--min-depth`: 10
- `--max-n-fraction`: 20
- `--edited`: false

### bam2mtx

The `bam2mtx` tool converts BAM files to single-cell matrices. This tool is designed for single-cell RNA-seq or ATAC-seq analysis where you want to generate count matrices from BAM files.

#### Usage

```bash
redicat bam2mtx [FLAGS] [OPTIONS] --bam <bam> --barcodes <barcodes> --output <output> --tsv <tsv>
```

#### Default Parameters

- `--threads`: 10
- `--min-mapq`: 255
- `--min-baseq`: 30
- `--cb-tag`: "CB"
- `--umi-tag`: "UB"
- `--chunksize`: 2500
- `--matrix-density`: 0.005
- `--use-mmap`: false
- `--large-dataset`: false
- `--memory-efficient`: false

### preprocess

The `preprocess` tool filters BAM files based on various criteria such as mapping quality, flags, and regions of interest.

#### Usage

```bash
redicat preprocess [FLAGS] [OPTIONS] --barcodes <barcodes> --inbam <inbam> --outbam <outbam>
```

#### Default Parameters

- `--mapquality`: 255
- `--include-flags`: None
- `--exclude-flags`: None
- `--threads`: 8
- `--write-cache-size`: 50000

### call

The `call` tool is the RNA editing detection and analysis pipeline. It processes single-cell RNA-seq data to identify RNA editing events and calculate editing indices.

#### Usage

```bash
redicat call [FLAGS] [OPTIONS] --fa <fa> --input <input> --output <output> --site-white-list <site-white-list>
```

#### Default Parameters

- `--editingtype`: ag
- `--max-other-threshold`: 0.01
- `--min-edited-threshold`: 0.01
- `--min-ref-threshold`: 0.01
- `--chunk-size`: 100000
- `--num-threads`: 2
- `--min-coverage`: 5
- `--verbose`: false
- `--dry-run`: false

## Library API

REdiCAT can also be used as a Rust library. The main modules are:

- `bam2mtx`: BAM to matrix conversion functionality
- `par_granges`: Parallel processing of genomic regions
- `position`: Data structures for genomic positions
- `read_filter`: Filtering of reads based on various criteria
- `utils`: Utility functions
- `call`: RNA editing detection and analysis functionality

For detailed API documentation, see the [API Reference](api.md).

## Examples

### RNA Editing Analysis Pipeline

1. Prepare input data (AnnData file, reference genome, editing site database)
2. Run the analysis pipeline:

```bash
redicat call \
  --input input.h5ad \
  --output output.h5ad \
  --fa reference.fa \
  --site-white-list editing_sites.tsv.gz \
  --editingtype ag \
  --num-threads 16
```

### Single-cell Matrix Generation

1. Prepare input files (BAM, position list, barcode list)
2. Generate matrix:

```bash
redicat bam2mtx \
  --bam input.bam \
  --tsv positions.tsv \
  --barcodes barcodes.tsv \
  --output output.h5ad \
  --num-threads 8
```

## Contributing

We welcome contributions to REdiCAT! Please see our [contributing guidelines](CONTRIBUTING.md) for more information.