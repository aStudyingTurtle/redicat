use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use polars::prelude::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::{Cursor, Read};
use std::path::Path;

/// Minimum buffer capacity used when chunking positions to avoid frequent reallocations.
const MIN_CHUNK_BUFFER_CAPACITY: usize = 512;

/// Represents a genomic position parsed from the TSV manifest.
#[derive(Debug, Clone)]
pub struct GenomicPosition {
    pub chrom: String,
    pub pos: u64,
    pub depth: u32,
    pub ins: u32,
    pub del: u32,
    pub ref_skip: u32,
    pub fail: u32,
    pub near_max_depth: bool,
}

/// Chunk of genomic positions for parallel processing.
#[derive(Debug, Clone)]
pub struct PositionChunk {
    pub positions: Vec<GenomicPosition>,
    pub near_max_depth_count: usize,
}

impl PositionChunk {
    #[inline]
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    #[inline]
    pub fn near_max_depth_count(&self) -> usize {
        self.near_max_depth_count
    }

    /// Calculate total weight (positions × depth) for capacity estimation.
    /// This helps estimate memory footprint for adaptive channel sizing.
    #[inline]
    pub fn total_weight(&self) -> u64 {
        self.positions
            .iter()
            .map(|p| p.depth as u64)
            .sum::<u64>()
            .saturating_mul(self.positions.len() as u64)
    }
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

    let ins_series_owned = dataframe
        .column("INS")
        .context("positions file is missing INS column")?
        .cast(&DataType::UInt32)
        .context("unable to cast INS to UInt32")?
        .clone();
    let ins_series = ins_series_owned
        .u32()
        .context("INS column must be unsigned integer")?;

    let del_series_owned = dataframe
        .column("DEL")
        .context("positions file is missing DEL column")?
        .cast(&DataType::UInt32)
        .context("unable to cast DEL to UInt32")?
        .clone();
    let del_series = del_series_owned
        .u32()
        .context("DEL column must be unsigned integer")?;

    let ref_skip_series_owned = dataframe
        .column("REF_SKIP")
        .context("positions file is missing REF_SKIP column")?
        .cast(&DataType::UInt32)
        .context("unable to cast REF_SKIP to UInt32")?
        .clone();
    let ref_skip_series = ref_skip_series_owned
        .u32()
        .context("REF_SKIP column must be unsigned integer")?;

    let fail_series_owned = dataframe
        .column("FAIL")
        .context("positions file is missing FAIL column")?
        .cast(&DataType::UInt32)
        .context("unable to cast FAIL to UInt32")?
        .clone();
    let fail_series = fail_series_owned
        .u32()
        .context("FAIL column must be unsigned integer")?;

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
        let Some(ins) = ins_series.get(idx) else {
            continue;
        };
        let Some(del) = del_series.get(idx) else {
            continue;
        };
        let Some(ref_skip) = ref_skip_series.get(idx) else {
            continue;
        };
        let Some(fail) = fail_series.get(idx) else {
            continue;
        };
        let near = near_series.get(idx).unwrap_or(false);

        positions.push(GenomicPosition {
            chrom: chrom.to_string(),
            pos,
            depth,
            ins,
            del,
            ref_skip,
            fail,
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
    let base_capacity = normal_limit
        .max(hotspot_limit)
        .max(MIN_CHUNK_BUFFER_CAPACITY);

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
                         chrom: &mut Option<String>,
                         near_count: &mut usize| {
        if chunk.is_empty() {
            return;
        }
        let near_total = *near_count;
        chunks.push(PositionChunk {
            positions: std::mem::take(chunk),
            near_max_depth_count: near_total,
        });
        chunk.reserve(base_capacity);
        *weight = 0;
        *chrom = None;
        *near_count = 0;
    };

    let mut current_near: usize = 0;

    for position in positions.into_iter() {
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
                &mut current_near,
            );
        }

        if current_chunk.is_empty() {
            current_kind = position_kind;
            current_chrom = Some(position.chrom.clone());
        }

        if position.near_max_depth {
            current_near = current_near.saturating_add(1);
        }

        // Calculate the weight for this position: DEPTH + INS + DEL + REF_SKIP + FAIL
        let position_weight = if current_kind == ChunkKind::Normal {
            (position.depth as usize)
                .saturating_add(position.ins as usize)
                .saturating_add(position.del as usize)
                .saturating_add(position.ref_skip as usize)
                .saturating_add(position.fail as usize)
                .max(1)
        } else {
            1
        };

        current_chunk.push(position);

        match current_kind {
            ChunkKind::Normal => {
                current_weight = current_weight.saturating_add(position_weight);
                if current_weight >= normal_limit {
                    flush_current(
                        &mut chunks,
                        &mut current_chunk,
                        &mut current_weight,
                        &mut current_chrom,
                        &mut current_near,
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
                        &mut current_near,
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
            &mut current_near,
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
                ins: 0,
                del: 0,
                ref_skip: 0,
                fail: 0,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 2,
                depth: 2,
                ins: 0,
                del: 0,
                ref_skip: 0,
                fail: 0,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 3,
                depth: 120,
                ins: 5,
                del: 3,
                ref_skip: 2,
                fail: 1,
                near_max_depth: true,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 4,
                depth: 95,
                ins: 4,
                del: 2,
                ref_skip: 1,
                fail: 2,
                near_max_depth: true,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 5,
                depth: 88,
                ins: 3,
                del: 1,
                ref_skip: 0,
                fail: 1,
                near_max_depth: true,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 6,
                depth: 1,
                ins: 0,
                del: 0,
                ref_skip: 0,
                fail: 0,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr1".to_string(),
                pos: 7,
                depth: 2,
                ins: 0,
                del: 0,
                ref_skip: 0,
                fail: 0,
                near_max_depth: false,
            },
        ];

        let chunks = chunk_positions(positions, 4, 2);

        for chunk in &chunks {
            let computed_near = chunk.positions.iter().filter(|p| p.near_max_depth).count();
            assert_eq!(chunk.near_max_depth_count(), computed_near);
            if chunk.positions.iter().any(|p| p.near_max_depth) {
                assert!(chunk.positions.len() <= 2);
                assert!(chunk
                    .positions
                    .iter()
                    .all(|position| position.near_max_depth));
                assert_eq!(chunk.near_max_depth_count(), chunk.positions.len());
            } else {
                let total_weight: usize = chunk
                    .positions
                    .iter()
                    .map(|p| (p.depth + p.ins + p.del + p.ref_skip + p.fail) as usize)
                    .sum();
                assert!(total_weight <= 4);
                assert_eq!(chunk.near_max_depth_count(), 0);
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
                ins: 1,
                del: 0,
                ref_skip: 0,
                fail: 0,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr2".to_string(),
                pos: 11,
                depth: 3,
                ins: 0,
                del: 1,
                ref_skip: 0,
                fail: 0,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr2".to_string(),
                pos: 12,
                depth: 4,
                ins: 0,
                del: 0,
                ref_skip: 1,
                fail: 0,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr2".to_string(),
                pos: 13,
                depth: 6,
                ins: 0,
                del: 0,
                ref_skip: 0,
                fail: 1,
                near_max_depth: false,
            },
        ];

        let chunks = chunk_positions(positions, 10, 2);
        assert_eq!(chunks.len(), 2);
        assert_eq!(chunks[0].positions.len(), 2); // First two positions: weight 6 + 4 = 10
        assert_eq!(chunks[1].positions.len(), 2); // Last two positions: weight 5 + 7 = 12
        assert_eq!(chunks[0].near_max_depth_count(), 0);
        assert_eq!(chunks[1].near_max_depth_count(), 0);

        let first_weights: usize = chunks[0]
            .positions
            .iter()
            .map(|p| (p.depth + p.ins + p.del + p.ref_skip + p.fail) as usize)
            .sum();
        assert!(first_weights >= 10);
        let second_weights: usize = chunks[1]
            .positions
            .iter()
            .map(|p| (p.depth + p.ins + p.del + p.ref_skip + p.fail) as usize)
            .sum();
        assert!(second_weights < 20);
    }

    #[test]
    fn chunk_positions_uses_comprehensive_weight() {
        // Test that the weight calculation includes all components
        let positions = vec![
            GenomicPosition {
                chrom: "chr3".to_string(),
                pos: 100,
                depth: 10,
                ins: 2,
                del: 3,
                ref_skip: 1,
                fail: 4,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr3".to_string(),
                pos: 101,
                depth: 15,
                ins: 1,
                del: 1,
                ref_skip: 2,
                fail: 1,
                near_max_depth: false,
            },
            GenomicPosition {
                chrom: "chr3".to_string(),
                pos: 102,
                depth: 5,
                ins: 0,
                del: 0,
                ref_skip: 0,
                fail: 0,
                near_max_depth: false,
            },
        ];

        // First position has weight: 10 + 2 + 3 + 1 + 4 = 20
        // Second position has weight: 15 + 1 + 1 + 2 + 1 = 20
        // Third position has weight: 5
        // With chunksize=25:
        //   - First position added: weight = 20 < 25, continue
        //   - Second position added: weight = 40 >= 25, flush (both in chunk 1)
        //   - Third position starts chunk 2: weight = 5
        let chunks = chunk_positions(positions, 25, 10);
        assert_eq!(chunks.len(), 2, "Should create 2 chunks");
        assert_eq!(
            chunks[0].positions.len(),
            2,
            "First chunk should have 2 positions"
        );
        assert_eq!(
            chunks[1].positions.len(),
            1,
            "Second chunk should have 1 position"
        );

        // Verify first chunk weight includes all components
        let first_weight: usize = chunks[0]
            .positions
            .iter()
            .map(|p| (p.depth + p.ins + p.del + p.ref_skip + p.fail) as usize)
            .sum();
        assert_eq!(
            first_weight, 40,
            "First chunk weight should be 40 (20 + 20)"
        );

        // Verify second chunk weight
        let second_weight: usize = chunks[1]
            .positions
            .iter()
            .map(|p| (p.depth + p.ins + p.del + p.ref_skip + p.fail) as usize)
            .sum();
        assert_eq!(second_weight, 5, "Second chunk weight should be 5");
    }
}
