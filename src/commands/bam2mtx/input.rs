use anyhow::Result;
use flate2::read::GzDecoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Represents a genomic position parsed from the TSV manifest.
#[derive(Debug, Clone)]
pub struct GenomicPosition {
    pub chrom: String,
    pub pos: u64,
}

/// Chunk of genomic positions for parallel processing.
#[derive(Debug, Clone)]
pub struct PositionChunk {
    pub positions: Vec<GenomicPosition>,
}

/// Read a TSV/TSV.GZ file into a sorted list of [`GenomicPosition`].
pub fn read_positions<P: AsRef<Path>>(tsv_path: P) -> Result<Vec<GenomicPosition>> {
    let path = tsv_path.as_ref();
    let file = File::open(path)?;

    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|ext| ext == "gz") {
        Box::new(BufReader::with_capacity(256 * 1024, GzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(256 * 1024, file))
    };

    let mut positions = Vec::new();
    let mut buffer = String::with_capacity(128);

    for line in reader.lines() {
        buffer.clear();
        buffer = line?;
        let trimmed = buffer.trim_end();

        let mut columns = trimmed.split('\t');
        if let (Some(chrom), Some(pos)) = (columns.next(), columns.next()) {
            if let Ok(pos) = pos.parse::<u64>() {
                positions.push(GenomicPosition {
                    chrom: chrom.to_string(),
                    pos,
                });
            }
        }
    }

    positions.par_sort_unstable_by(|a, b| match a.chrom.cmp(&b.chrom) {
        std::cmp::Ordering::Equal => a.pos.cmp(&b.pos),
        other => other,
    });

    positions.shrink_to_fit();
    Ok(positions)
}

/// Split positions into contig-aware chunks to preserve locality.
pub fn chunk_positions(
    mut positions: Vec<GenomicPosition>,
    chunk_size: usize,
) -> Vec<PositionChunk> {
    if positions.is_empty() {
        return Vec::new();
    }

    let mut chunks = Vec::with_capacity(positions.len().div_ceil(chunk_size));
    let mut current_chunk = Vec::with_capacity(chunk_size);
    let mut current_chrom: Option<String> = None;

    for position in positions.drain(..) {
        if current_chrom
            .as_ref()
            .is_some_and(|chromosome| chromosome != &position.chrom)
            && !current_chunk.is_empty()
        {
            chunks.push(PositionChunk {
                positions: std::mem::take(&mut current_chunk),
            });
            current_chrom = None;
        }

        if current_chunk.is_empty() {
            current_chrom = Some(position.chrom.clone());
        }

        current_chunk.push(position);

        if current_chunk.len() >= chunk_size {
            chunks.push(PositionChunk {
                positions: std::mem::take(&mut current_chunk),
            });
            current_chrom = None;
        }
    }

    if !current_chunk.is_empty() {
        chunks.push(PositionChunk {
            positions: current_chunk,
        });
    }

    chunks.shrink_to_fit();
    chunks
}

/// Retain only positions that satisfy the supplied predicate.
pub fn filter_positions_by<F>(positions: Vec<GenomicPosition>, mut keep: F) -> Vec<GenomicPosition>
where
    F: FnMut(&str) -> bool,
{
    positions
        .into_iter()
        .filter(|pos| keep(&pos.chrom))
        .collect()
}
