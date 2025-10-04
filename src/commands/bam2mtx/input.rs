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
    pub near_max_depth: bool,
}

/// Chunk of genomic positions for parallel processing.
#[derive(Debug, Clone)]
pub struct PositionChunk {
    pub positions: Vec<GenomicPosition>,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum ChunkKind {
    Normal,
    NearMax,
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
        if let (Some(chrom), Some(pos_str)) = (columns.next(), columns.next()) {
            if let Ok(pos) = pos_str.parse::<u64>() {
                let near_max_depth = columns
                    .last()
                    .map(|value| {
                        let trimmed = value.trim();
                        trimmed.eq_ignore_ascii_case("true") || trimmed == "1"
                    })
                    .unwrap_or(false);

                positions.push(GenomicPosition {
                    chrom: chrom.to_string(),
                    pos,
                    near_max_depth,
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
    chunk_size_max_depth: usize,
) -> Vec<PositionChunk> {
    if positions.is_empty() {
        return Vec::new();
    }

    let normal_capacity = chunk_size.max(1);
    let near_capacity = chunk_size_max_depth.max(1);
    let base_capacity = normal_capacity.max(near_capacity);

    let mut chunks = Vec::with_capacity(positions.len().div_ceil(normal_capacity));
    let mut current_chunk = Vec::with_capacity(base_capacity);
    let mut current_chrom: Option<String> = None;
    let mut current_kind = ChunkKind::Normal;
    let mut current_capacity = normal_capacity;

    for position in positions.drain(..) {
        let position_kind = if position.near_max_depth {
            ChunkKind::NearMax
        } else {
            ChunkKind::Normal
        };

        let contig_changed = current_chrom
            .as_ref()
            .is_some_and(|chromosome| chromosome != &position.chrom);

        if (!current_chunk.is_empty()) && (contig_changed || current_kind != position_kind) {
            chunks.push(PositionChunk {
                positions: std::mem::take(&mut current_chunk),
            });
            current_chunk.reserve(base_capacity);
            current_chrom = None;
        }

        if current_chunk.is_empty() {
            current_kind = position_kind;
            current_capacity = match current_kind {
                ChunkKind::Normal => normal_capacity,
                ChunkKind::NearMax => near_capacity,
            };
            current_chrom = Some(position.chrom.clone());
        }

        current_chunk.push(position);

        if current_chunk.len() >= current_capacity {
            chunks.push(PositionChunk {
                positions: std::mem::take(&mut current_chunk),
            });
            current_chunk.reserve(base_capacity);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chunk_positions_respects_near_max_chunk_size() {
        let positions = vec![
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 1,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 2,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 3,
                near_max_depth: true,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 4,
                near_max_depth: true,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 5,
                near_max_depth: true,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 6,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 7,
                near_max_depth: false,
            },
        ];

        let chunks = chunk_positions(positions, 4, 2);
        assert_eq!(chunks.len(), 4);

        for chunk in &chunks {
            if chunk.positions.iter().any(|p| p.near_max_depth) {
                assert!(chunk.positions.len() <= 2);
                assert!(chunk
                    .positions
                    .iter()
                    .all(|position| position.near_max_depth));
            } else {
                assert!(chunk.positions.len() <= 4);
            }
        }
    }
}
