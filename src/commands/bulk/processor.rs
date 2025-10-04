use anyhow::{Context, Result};
use noodles::{
    bam,
    bgzf,
    core::{self, Position, Region},
    sam,
};
use noodles::sam::alignment::record::cigar::op::Kind;
use redicat_lib::core::read_filter::{DefaultReadFilter, ReadFilter};
use redicat_lib::engine::position::pileup_position::PileupPosition;
use redicat_lib::engine::RegionProcessor;
use rustc_hash::{FxHashMap, FxHashSet};
use smartstring::SmartString;
use std::cmp::max;
use std::convert::TryFrom;
use std::fs::File;
use std::path::PathBuf;
use std::sync::Arc;

use super::args::BulkConfig;

type IndexedBamReader = bam::io::IndexedReader<bgzf::io::Reader<File>>;

/// Implements [`RegionProcessor`] for bulk pileup traversal.
pub struct BaseProcessor {
    reads: PathBuf,
    coord_base: u32,
    max_depth: u32,
    skip_max_depth: u32,
    min_depth: u32,
    max_n_fraction: u32,
    min_baseq: Option<u8>,
    edited_mode: bool,
    editing_threshold: u32,
    read_filter: Arc<DefaultReadFilter>,
    header: Arc<sam::Header>,
    allowed_tids: Option<FxHashSet<u32>>,
}

impl BaseProcessor {
    pub fn from_config(
        config: &BulkConfig,
        read_filter: Arc<DefaultReadFilter>,
        allowed_tids: Option<FxHashSet<u32>>,
    ) -> Result<Self> {
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&config.reads)
            .with_context(|| format!("Failed to open {}", config.reads.display()))?;
        let header = Arc::new(reader.read_header().context("Failed to read BAM header")?);
        drop(reader);

        let adjusted_max_depth = if config.skip_max_depth < 2_000_000_000 {
            config.skip_max_depth + 1000
        } else {
            config.max_depth
        };

        Ok(Self {
            reads: config.reads.clone(),
            coord_base: config.coord_offset,
            max_depth: adjusted_max_depth,
            skip_max_depth: config.skip_max_depth,
            min_depth: config.min_depth,
            max_n_fraction: config.max_n_fraction,
            min_baseq: config.min_baseq.or(Some(30)),
            edited_mode: config.editing_mode(),
            editing_threshold: config.editing_threshold,
            read_filter,
            header,
            allowed_tids,
        })
    }

    fn contig_name(&self, tid: usize) -> Option<String> {
        self
            .header
            .reference_sequences()
            .get_index(tid)
            .and_then(|(name, _)| std::str::from_utf8(name.as_ref()).ok())
            .map(|s| s.to_string())
    }

    fn process_region_inner(
        &self,
        reader: &mut IndexedBamReader,
        tid: u32,
        start: u32,
        stop: u32,
    ) -> Vec<PileupPosition> {
        if start >= stop {
            return Vec::new();
        }

        let Some(chrom_name) = self.contig_name(tid as usize) else {
            return Vec::new();
        };

        let start_pos = match Position::try_from(start as usize + 1) {
            Ok(pos) => pos,
            Err(_) => return Vec::new(),
        };

        let end_pos = match Position::try_from(stop as usize) {
            Ok(pos) => pos,
            Err(_) => return Vec::new(),
        };

        let interval = core::region::Interval::from(start_pos..=end_pos);

        let region = Region::new(chrom_name.clone(), interval);

        let mut query = match reader.query(&self.header, &region) {
            Ok(query) => query,
            Err(err) => {
                log::warn!("Failed to query region {chrom_name}:{start}-{stop}: {err}");
                return Vec::new();
            }
        };

        let mut positions: FxHashMap<u32, PileupPosition> = FxHashMap::default();
        let chrom_label = SmartString::from(chrom_name.as_str());
        let region_start = start as i64;
        let region_stop = stop as i64;

        while let Some(record_result) = query.next() {
            let record = match record_result {
                Ok(record) => record,
                Err(err) => {
                    log::debug!("Skipping malformed record: {err}");
                    continue;
                }
            };

            let is_target_tid = match record.reference_sequence_id().transpose() {
                Ok(Some(id)) => id == tid as usize,
                _ => false,
            };

            if !is_target_tid {
                continue;
            }

            let alignment_start = match record.alignment_start().transpose() {
                Ok(Some(pos)) => usize::from(pos).saturating_sub(1),
                _ => continue,
            } as i64;

            let passes_filter = self.read_filter.filter_read(&record);
            let sequence = record.sequence();
            let qualities = record.quality_scores();

            let mut ref_pos = alignment_start;
            let mut read_pos = 0usize;

            for op_result in record.cigar().iter() {
                let op = match op_result {
                    Ok(op) => op,
                    Err(err) => {
                        log::debug!("Skipping invalid CIGAR op: {err}");
                        continue;
                    }
                };

                let len = op.len() as usize;

                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        for i in 0..len {
                            let current_ref = ref_pos + i as i64;
                            if current_ref < region_start || current_ref >= region_stop || current_ref < 0 {
                                continue;
                            }

                            let Ok(pos_u32) = u32::try_from(current_ref) else {
                                continue;
                            };

                            let entry = positions.entry(pos_u32).or_insert_with(|| {
                                PileupPosition::create(chrom_label.clone(), pos_u32)
                            });

                            if !passes_filter {
                                entry.fail = entry.fail.saturating_add(1);
                                continue;
                            }

                            if entry.depth >= self.max_depth {
                                continue;
                            }

                            entry.depth = entry.depth.saturating_add(1);

                            let qpos = read_pos + i;
                            let base_qual = qualities.as_ref().get(qpos).copied().unwrap_or(0);
                            if let Some(cutoff) = self.min_baseq {
                                if base_qual < cutoff {
                                    entry.n = entry.n.saturating_add(1);
                                    continue;
                                }
                            }

                            match sequence.get(qpos).map(|b| b.to_ascii_uppercase()) {
                                Some(b'A') => entry.a += 1,
                                Some(b'C') => entry.c += 1,
                                Some(b'G') => entry.g += 1,
                                Some(b'T') => entry.t += 1,
                                _ => entry.n += 1,
                            }
                        }

                        ref_pos += len as i64;
                        read_pos += len;
                    }
                    Kind::Deletion => {
                        for i in 0..len {
                            let current_ref = ref_pos + i as i64;
                            if current_ref < region_start || current_ref >= region_stop || current_ref < 0 {
                                continue;
                            }

                            let Ok(pos_u32) = u32::try_from(current_ref) else {
                                continue;
                            };

                            let entry = positions.entry(pos_u32).or_insert_with(|| {
                                PileupPosition::create(chrom_label.clone(), pos_u32)
                            });

                            if !passes_filter {
                                entry.fail = entry.fail.saturating_add(1);
                                continue;
                            }

                            if entry.depth >= self.max_depth {
                                continue;
                            }

                            entry.depth = entry.depth.saturating_add(1);
                            entry.del += 1;
                        }

                        ref_pos += len as i64;
                    }
                    Kind::Skip => {
                        for i in 0..len {
                            let current_ref = ref_pos + i as i64;
                            if current_ref < region_start || current_ref >= region_stop || current_ref < 0 {
                                continue;
                            }

                            let Ok(pos_u32) = u32::try_from(current_ref) else {
                                continue;
                            };

                            let entry = positions.entry(pos_u32).or_insert_with(|| {
                                PileupPosition::create(chrom_label.clone(), pos_u32)
                            });

                            if !passes_filter {
                                entry.fail = entry.fail.saturating_add(1);
                            } else {
                                entry.ref_skip = entry.ref_skip.saturating_add(1);
                            }
                        }

                        ref_pos += len as i64;
                    }
                    Kind::Insertion => {
                        if passes_filter {
                            let insertion_pos = ref_pos - 1;
                            if insertion_pos >= region_start
                                && insertion_pos < region_stop
                                && insertion_pos >= 0
                            {
                                if let Ok(pos_u32) = u32::try_from(insertion_pos) {
                                    let entry = positions.entry(pos_u32).or_insert_with(|| {
                                        PileupPosition::create(chrom_label.clone(), pos_u32)
                                    });
                                    entry.ins = entry.ins.saturating_add(1);
                                }
                            }
                        }
                        read_pos += len;
                    }
                    Kind::SoftClip => {
                        read_pos += len;
                    }
                    Kind::HardClip | Kind::Pad => {}
                }
            }
        }

        let mut aggregate: Vec<_> = positions
            .into_iter()
            .map(|(_, mut pos)| {
                pos.pos = pos.pos.saturating_add(self.coord_base);
                pos
            })
            .collect();

        aggregate.sort_by_key(|p| p.pos);

        let mut filtered = Vec::with_capacity(aggregate.len());
        for mut pos in aggregate {
            pos.mark_near_max_depth(self.max_depth);

            if pos.depth > self.skip_max_depth {
                continue;
            }

            let should_include = if self.edited_mode {
                let valid_value = max(pos.depth / self.editing_threshold, 2);
                let mut count = 0;
                if pos.a > valid_value {
                    count += 1;
                }
                if pos.t > valid_value {
                    count += 1;
                }
                if pos.g > valid_value {
                    count += 1;
                }
                if pos.c > valid_value {
                    count += 1;
                }

                (pos.depth >= self.min_depth)
                    && (count >= 2)
                    && (pos.n <= max(pos.depth / self.max_n_fraction, 2))
            } else {
                pos.depth >= self.min_depth
            };

            if should_include {
                filtered.push(pos);
            }
        }

        filtered
    }
}

impl RegionProcessor for BaseProcessor {
    type P = PileupPosition;

    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P> {
        if let Some(allowed) = &self.allowed_tids {
            if !allowed.contains(&tid) {
                return Vec::new();
            }
        }

        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&self.reads)
            .unwrap_or_else(|e| panic!("Failed to open {}: {}", self.reads.display(), e));

        reader
            .read_header()
            .unwrap_or_else(|e| panic!("Failed to read header for {}: {}", self.reads.display(), e));

        self.process_region_inner(&mut reader, tid, start, stop)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use redicat_lib::core::read_filter::DefaultReadFilter;
    use std::path::PathBuf;
    use std::sync::Arc;

    #[test]
    fn test_process_region_chr22_smoke() {
        let bam = PathBuf::from("test/chr22.bam");
        if !bam.exists() {
            log::info!("Skipping chr22 smoke test – fixture missing");
            return;
        }

        let config = BulkConfig {
            reads: bam.clone(),
            output: PathBuf::from("test.tsv.gz"),
            threads: 2,
            chunksize: 10_000,
            min_baseq: Some(30),
            mapquality: 255,
            coord_offset: 1,
            max_depth: 10_000,
            skip_max_depth: u32::MAX,
            min_depth: 10,
            max_n_fraction: 20,
            report_all: true,
            editing_threshold: 1000,
            all_contigs: true,
        };

        let processor = BaseProcessor::from_config(
            &config,
            Arc::new(DefaultReadFilter::new(config.mapquality)),
            None,
        )
        .expect("processor initialisation");

        let results = processor.process_region(14, 22_901_237, 22_901_238);
        log::info!("Smoke test processed {} positions", results.len());
    }
}
