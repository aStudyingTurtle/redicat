use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use polars::prelude::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{Cursor, Read};
use std::path::Path;

/// Represents a genomic position parsed from the TSV manifest.
#[derive(Debug, Clone)]
pub struct GenomicPosition {
    pub chrom: String,
    pub pos: u64,
    pub depth: u32,
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

    let mut raw = Vec::new();
    if path.extension().is_some_and(|ext| ext == "gz") {
        let mut decoder = GzDecoder::new(File::open(path)?);
        decoder.read_to_end(&mut raw)?;
    } else {
        let mut file = File::open(path)?;
        file.read_to_end(&mut raw)?;
    }

    let mut processed = Vec::with_capacity(raw.len());
    let mut start = 0;
    while start < raw.len() {
        let end = raw[start..]
            .iter()
            .position(|&b| b == b'\n')
            .map(|offset| start + offset)
            .unwrap_or(raw.len());
        let line = &raw[start..end];

        if !line.is_empty() && line[0] != b'#' {
            for &byte in line {
                if byte == b'\r' {
                    continue;
                }
                processed.push(if byte == b'\t' { b',' } else { byte });
            }
            processed.push(b'\n');
        }

        start = end.saturating_add(1);
    }

    let dataframe = CsvReader::new(Cursor::new(processed))
        .finish()
        .with_context(|| format!("Failed to parse positions from {:?}", path))?;

    let chrom_series = dataframe
        .column("CHR")
        .context("positions file is missing CHR column")?
        .str()
        .context("CHR column must be UTF-8")?;

    let pos_series_owned = dataframe
        .column("POS")
        .context("positions file is missing POS column")?
        .cast(&DataType::UInt64)
        .context("unable to cast POS to UInt64")?
        .clone();
    let pos_series = pos_series_owned
        .u64()
        .context("POS column must be unsigned integer")?;

    let depth_series_owned = dataframe
        .column("DEPTH")
        .context("positions file is missing DEPTH column")?
        .cast(&DataType::UInt32)
        .context("unable to cast DEPTH to UInt32")?
        .clone();
    let depth_series = depth_series_owned
        .u32()
        .context("DEPTH column must be unsigned integer")?;

    let near_series_owned = dataframe
        .column("NEAR_MAX_DEPTH")
        .context("positions file is missing NEAR_MAX_DEPTH column")?
        .cast(&DataType::Boolean)
        .context("unable to cast NEAR_MAX_DEPTH to Boolean")?
        .clone();
    let near_series = near_series_owned
        .bool()
        .context("NEAR_MAX_DEPTH column must be boolean")?;

    let mut positions = Vec::with_capacity(dataframe.height());

    for idx in 0..dataframe.height() {
        let Some(chrom) = chrom_series.get(idx) else {
            continue;
        };
        let Some(pos) = pos_series.get(idx) else {
            continue;
        };
        let Some(depth) = depth_series.get(idx) else {
            continue;
        };
        let near = near_series.get(idx).unwrap_or(false);

        positions.push(GenomicPosition {
            chrom: chrom.to_string(),
            pos,
            depth,
            near_max_depth: near,
        });
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
    positions: Vec<GenomicPosition>,
    chunk_size: usize,
    chunk_size_max_depth: usize,
) -> Vec<PositionChunk> {
    if positions.is_empty() {
        return Vec::new();
    }

    let normal_limit = chunk_size.max(1);
    let hotspot_limit = chunk_size_max_depth.max(1);
    let base_capacity = normal_limit.max(hotspot_limit);

    let estimated_chunks = if hotspot_limit == 0 {
        positions.len()
    } else {
        positions.len().saturating_add(hotspot_limit - 1) / hotspot_limit
    };
    let mut chunks = Vec::with_capacity(estimated_chunks.max(1));
    let mut current_chunk: Vec<GenomicPosition> = Vec::with_capacity(base_capacity);
    let mut current_chrom: Option<String> = None;
    let mut current_kind = ChunkKind::Normal;
    let mut current_weight: usize = 0;

    let flush_current = |chunks: &mut Vec<PositionChunk>,
                         chunk: &mut Vec<GenomicPosition>,
                         weight: &mut usize,
                         chrom: &mut Option<String>| {
        if chunk.is_empty() {
            return;
        }
        chunks.push(PositionChunk {
            positions: std::mem::take(chunk),
        });
        chunk.reserve(base_capacity);
        *weight = 0;
        *chrom = None;
    };

    for mut position in positions.into_iter() {
        let position_kind = if position.near_max_depth {
            ChunkKind::NearMax
        } else {
            ChunkKind::Normal
        };

        let contig_changed = current_chrom
            .as_ref()
            .is_some_and(|chrom| chrom != &position.chrom);

        if !current_chunk.is_empty() && (contig_changed || current_kind != position_kind) {
            flush_current(
                &mut chunks,
                &mut current_chunk,
                &mut current_weight,
                &mut current_chrom,
            );
        }

        if current_chunk.is_empty() {
            current_kind = position_kind;
            current_chrom = Some(position.chrom.clone());
        }

        if current_kind == ChunkKind::Normal && position.depth == 0 {
            position.depth = 1;
        }

        let increment = match current_kind {
            ChunkKind::Normal => usize::max(position.depth as usize, 1),
            ChunkKind::NearMax => 1,
        };

        current_chunk.push(position);

        match current_kind {
            ChunkKind::Normal => {
                current_weight = current_weight.saturating_add(increment);
                if current_weight >= normal_limit {
                    flush_current(
                        &mut chunks,
                        &mut current_chunk,
                        &mut current_weight,
                        &mut current_chrom,
                    );
                }
            }
            ChunkKind::NearMax => {
                if current_chunk.len() >= hotspot_limit {
                    flush_current(
                        &mut chunks,
                        &mut current_chunk,
                        &mut current_weight,
                        &mut current_chrom,
                    );
                }
            }
        }
    }

    if !current_chunk.is_empty() {
        flush_current(
            &mut chunks,
            &mut current_chunk,
            &mut current_weight,
            &mut current_chrom,
        );
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
                depth: 2,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 2,
                depth: 2,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 3,
                depth: 120,
                near_max_depth: true,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 4,
                depth: 95,
                near_max_depth: true,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 5,
                depth: 88,
                near_max_depth: true,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 6,
                depth: 1,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 7,
                depth: 2,
                near_max_depth: false,
            },
        ];

        let chunks = chunk_positions(positions, 4, 2);

        for chunk in &chunks {
            if chunk.positions.iter().any(|p| p.near_max_depth) {
                assert!(chunk.positions.len() <= 2);
                assert!(chunk
                    .positions
                    .iter()
                    .all(|position| position.near_max_depth));
            } else {
                let total_weight: usize = chunk.positions.iter().map(|p| p.depth as usize).sum();
                assert!(total_weight <= 4);
            }
        }
    }

    #[test]
    fn chunk_positions_accumulates_depth_weight() {
        let positions = vec![
            GenomicPosition {
                chrom: "chr2".to_string(),
                pos: 10,
                depth: 5,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr2".to_string(),
                pos: 11,
                depth: 3,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr2".to_string(),
                pos: 12,
                depth: 4,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr2".to_string(),
                pos: 13,
                depth: 6,
                near_max_depth: false,
            },
        ];

        let chunks = chunk_positions(positions, 10, 2);
        assert_eq!(chunks.len(), 2);
        assert_eq!(chunks[0].positions.len(), 3);
        assert_eq!(chunks[1].positions.len(), 1);

        let first_weights: usize = chunks[0].positions.iter().map(|p| p.depth as usize).sum();
        assert!(first_weights >= 10);
        let second_weights: usize = chunks[1].positions.iter().map(|p| p.depth as usize).sum();
        assert!(second_weights < 10);
    }
}
