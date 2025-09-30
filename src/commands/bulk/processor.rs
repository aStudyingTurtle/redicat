use parking_lot::Mutex;
use redicat_lib::core::read_filter::DefaultReadFilter;
use redicat_lib::engine::position::pileup_position::PileupPosition;
use redicat_lib::engine::RegionProcessor;
use rust_htslib::{bam, bam::Read};
use rustc_hash::FxHashSet;
use std::cmp::max;
use std::convert::TryInto;
use std::path::PathBuf;
use std::sync::Arc;

use super::args::BulkConfig;

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
    reader_pool: Mutex<Vec<bam::IndexedReader>>,
    allowed_tids: Option<FxHashSet<u32>>,
}

impl BaseProcessor {
    pub fn from_config(
        config: &BulkConfig,
        read_filter: Arc<DefaultReadFilter>,
        allowed_tids: Option<FxHashSet<u32>>,
    ) -> Self {
        Self {
            reads: config.reads.clone(),
            coord_base: config.coord_offset,
            max_depth: config.max_depth,
            skip_max_depth: config.skip_max_depth,
            min_depth: config.min_depth,
            max_n_fraction: config.max_n_fraction,
            min_baseq: config.min_baseq.or(Some(30)),
            edited_mode: config.editing_mode(),
            editing_threshold: config.editing_threshold,
            read_filter,
            reader_pool: Mutex::new(Vec::new()),
            allowed_tids,
        }
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

        let mut reader = {
            let mut pool = self.reader_pool.lock();
            pool.pop()
        }
        .unwrap_or_else(|| {
            bam::IndexedReader::from_path(&self.reads)
                .unwrap_or_else(|e| panic!("Failed to open {}: {}", self.reads.display(), e))
        });

        let header = reader.header().to_owned();
        reader
            .fetch((tid, start, stop))
            .expect("Fetched requested region");
        let mut pileup = reader.pileup();
        pileup.set_max_depth(std::cmp::min(i32::MAX.try_into().unwrap(), self.max_depth));

        let estimated = (stop.saturating_sub(start) / 10).max(16) as usize;
        let mut result = Vec::with_capacity(estimated);
        for p in pileup {
            let pileup = match p {
                Ok(p) => p,
                Err(_) => continue,
            };

            if pileup.pos() < start || pileup.pos() >= stop {
                continue;
            }

            let mut pos = PileupPosition::from_pileup(
                pileup,
                &header,
                self.read_filter.as_ref(),
                self.min_baseq,
            );
            pos.pos += self.coord_base;

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

                ((pos.depth - pos.n) >= self.min_depth)
                    && (count >= 2)
                    && (pos.n <= max(pos.depth / self.max_n_fraction, 2))
            } else {
                (pos.depth - pos.n) >= self.min_depth
            };

            if should_include {
                result.push(pos);
            }
        }

        {
            let mut pool = self.reader_pool.lock();
            pool.push(reader);
        }

        result
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
        );

        let results = processor.process_region(14, 22_901_237, 22_901_238);
        log::info!("Smoke test processed {} positions", results.len());
    }
}
