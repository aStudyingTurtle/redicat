use anyhow::Result;
use parking_lot::Mutex;
use redicat_lib::bam2mtx::barcode::BarcodeProcessor;
use redicat_lib::bam2mtx::processor::{
    apply_encoded_call, decode_base, decode_cell_barcode, decode_umi, encode_call,
    BamProcessorConfig, PositionData, StrandBaseCounts, UMI_CONFLICT_CODE,
};
use rust_htslib::bam::{self, Read};
use rustc_hash::FxHashMap;
use std::hash::{Hash, Hasher};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use super::input::PositionChunk;

/// Number of shards for UMI interner (power of 2 for fast modulo)
const UMI_INTERNER_SHARDS: usize = 16;

/// Sharded UMI string interner to reduce lock contention.
/// Uses 16 independent shards, reducing contention by ~90%.
struct ShardedUmiInterner {
    shards: [Mutex<FxHashMap<String, Arc<str>>>; UMI_INTERNER_SHARDS],
}

impl ShardedUmiInterner {
    fn new() -> Self {
        // Create array of empty hashmaps wrapped in Mutex
        Self {
            shards: std::array::from_fn(|_| Mutex::new(FxHashMap::default())),
        }
    }

    /// Intern a UMI string, returning a shared Arc<str>.
    /// Uses FxHash for shard selection (fast, non-crypto hash).
    fn intern(&self, umi: &str) -> Arc<str> {
        // Use FxHash (same as FxHashMap) for consistent hashing
        let mut hasher = rustc_hash::FxHasher::default();
        umi.hash(&mut hasher);
        let hash = hasher.finish();
        let shard_idx = (hash as usize) % UMI_INTERNER_SHARDS;

        let mut shard = self.shards[shard_idx].lock();
        shard
            .entry(umi.to_string())
            .or_insert_with(|| Arc::from(umi))
            .clone()
    }

    /// Get statistics about interner usage (for debugging/monitoring).
    #[allow(dead_code)]
    fn stats(&self) -> Vec<usize> {
        self.shards.iter().map(|shard| shard.lock().len()).collect()
    }
}

/// Metadata for genomic sites skipped due to exceeding the configured depth ceiling.
#[derive(Debug, Clone)]
pub struct SkippedSite {
    pub contig: String,
    pub pos: u64,
    pub depth: u32,
}

/// Chunk processor with aggressive reuse of BAM readers and pre-sized hash maps.
/// Now includes sharded UMI string interning to reduce memory usage and lock contention.
pub struct OptimizedChunkProcessor {
    bam_path: PathBuf,
    config: BamProcessorConfig,
    barcode_processor: Arc<BarcodeProcessor>,
    tid_lookup: FxHashMap<String, u32>,
    contig_names: Arc<Vec<String>>,
    count_capacity_hint: usize,
    umi_capacity_hint: usize,
    reader_pool: Mutex<Vec<bam::IndexedReader>>,
    reader_pool_size: AtomicUsize,         // Lock-free size tracking
    umi_interner: Arc<ShardedUmiInterner>, // Sharded for reduced contention
}

impl OptimizedChunkProcessor {
    pub fn new(
        bam_path: PathBuf,
        config: BamProcessorConfig,
        barcode_processor: Arc<BarcodeProcessor>,
    ) -> Result<Self> {
        let header = bam::IndexedReader::from_path(&bam_path)?
            .header()
            .to_owned();
        let mut tid_lookup =
            FxHashMap::with_capacity_and_hasher(header.target_count() as usize, Default::default());

        let mut contig_names = Vec::with_capacity(header.target_count() as usize);
        for tid in 0..header.target_count() {
            if let Ok(name) = std::str::from_utf8(header.tid2name(tid)) {
                tid_lookup.insert(name.to_string(), tid);
                contig_names.push(name.to_string());
            } else {
                contig_names.push(String::from("unknown"));
            }
        }

        let barcode_count = barcode_processor.len();
        let count_capacity_hint = barcode_count.clamp(64, 4096);
        let umi_capacity_hint = count_capacity_hint.saturating_mul(8);

        Ok(Self {
            bam_path,
            config,
            barcode_processor,
            tid_lookup,
            contig_names: Arc::new(contig_names),
            count_capacity_hint,
            umi_capacity_hint,
            reader_pool: Mutex::new(Vec::new()),
            reader_pool_size: AtomicUsize::new(0),
            umi_interner: Arc::new(ShardedUmiInterner::new()),
        })
    }

    pub fn contig_names(&self) -> Arc<Vec<String>> {
        Arc::clone(&self.contig_names)
    }

    pub fn process_chunk(
        &self,
        chunk: &PositionChunk,
    ) -> Result<(Vec<PositionData>, Vec<SkippedSite>)> {
        if chunk.is_empty() {
            return Ok((Vec::new(), Vec::new()));
        }

        // Optimized reader acquisition with lock-free empty check
        let mut reader = if self.reader_pool_size.load(Ordering::Acquire) > 0 {
            // Pool might have readers, try to get one
            let mut pool = self.reader_pool.lock();
            if let Some(r) = pool.pop() {
                self.reader_pool_size.fetch_sub(1, Ordering::Release);
                r
            } else {
                // Race condition: another thread took it, create new reader
                drop(pool); // Release lock before I/O
                bam::IndexedReader::from_path(&self.bam_path).unwrap_or_else(|e| {
                    panic!("Failed to open BAM {}: {}", self.bam_path.display(), e)
                })
            }
        } else {
            // Pool is empty, create new reader without locking
            bam::IndexedReader::from_path(&self.bam_path)
                .unwrap_or_else(|e| panic!("Failed to open BAM {}: {}", self.bam_path.display(), e))
        };

        let chrom = &chunk.positions[0].chrom;
        let tid = match self.tid_lookup.get(chrom) {
            Some(&tid) => tid,
            None => {
                // Return reader to pool with atomic size update
                let mut pool = self.reader_pool.lock();
                pool.push(reader);
                self.reader_pool_size.fetch_add(1, Ordering::Release);
                return Ok((Vec::new(), Vec::new()));
            }
        };

        let fetch_start = chunk
            .positions
            .first()
            .map(|p| p.pos.saturating_sub(1) as u32)
            .unwrap_or(0);
        let fetch_end = chunk
            .positions
            .last()
            .map(|p| p.pos as u32)
            .unwrap_or(fetch_start);

        reader.fetch((tid, fetch_start, fetch_end))?;

        let target_positions: Vec<u32> = chunk
            .positions
            .iter()
            .map(|p| p.pos.saturating_sub(1) as u32)
            .collect();

        let mut chunk_results = Vec::with_capacity(chunk.len());
        let mut skipped_sites: Vec<SkippedSite> = Vec::new();
        let mut position_index = 0usize;

        // Improved HashMap capacity estimation based on chunk characteristics
        let has_high_depth = chunk.near_max_depth_count() > 0;
        let estimated_count_capacity = if has_high_depth {
            self.count_capacity_hint.min(5_000)
        } else {
            self.count_capacity_hint
        };
        let estimated_umi_capacity = if has_high_depth {
            self.umi_capacity_hint.min(50_000)
        } else {
            self.umi_capacity_hint
        };

        let mut counts: FxHashMap<u32, StrandBaseCounts> =
            FxHashMap::with_capacity_and_hasher(estimated_count_capacity, Default::default());
        let mut umi_consensus: FxHashMap<(u32, Arc<str>), u8> =
            FxHashMap::with_capacity_and_hasher(estimated_umi_capacity, Default::default());

        let mut pileups = reader.pileup();
        pileups.set_max_depth(self.config.max_depth.min(i32::MAX as u32));

        for pileup in pileups {
            let pileup = pileup?;
            let pile_pos = pileup.pos();

            while position_index < target_positions.len()
                && target_positions[position_index] < pile_pos
            {
                position_index += 1;
            }

            if position_index >= target_positions.len() {
                break;
            }

            if target_positions[position_index] != pile_pos {
                continue;
            }

            let current_index = position_index;
            position_index += 1;

            let depth = pileup.depth() as u32;
            if depth >= self.config.max_depth {
                let position_meta = &chunk.positions[current_index];
                skipped_sites.push(SkippedSite {
                    contig: position_meta.chrom.clone(),
                    pos: position_meta.pos,
                    depth,
                });
                continue;
            }

            counts.clear();
            umi_consensus.clear();

            for alignment in pileup.alignments() {
                let record = alignment.record();

                if record.mapq() < self.config.min_mapping_quality {
                    continue;
                }

                let qpos = match alignment.qpos() {
                    Some(q) => q,
                    None => continue,
                };

                let base_qual = record.qual().get(qpos).copied().unwrap_or(0);
                if base_qual < self.config.min_base_quality {
                    continue;
                }

                let base = decode_base(&record, Some(qpos))?;
                if base == 'N' {
                    continue;
                }

                let cell_id =
                    match decode_cell_barcode(&record, self.config.cell_barcode_tag.as_bytes())? {
                        Some(barcode) => match self.barcode_processor.id_of(&barcode) {
                            Some(id) => id,
                            None => continue,
                        },
                        None => continue,
                    };

                let umi = match decode_umi(&record, self.config.umi_tag.as_bytes())? {
                    Some(umi) => umi,
                    None => continue,
                };

                // Intern UMI string using sharded interner (reduced lock contention)
                let umi_arc = self.umi_interner.intern(&umi);

                if let Some(encoded) = encode_call(self.config.stranded, base, record.is_reverse())
                {
                    umi_consensus
                        .entry((cell_id, umi_arc))
                        .and_modify(|existing| {
                            if *existing != encoded {
                                *existing = UMI_CONFLICT_CODE;
                            }
                        })
                        .or_insert(encoded);
                }
            }

            for ((cell_id, _umi), encoded) in umi_consensus.drain() {
                if encoded == UMI_CONFLICT_CODE {
                    continue;
                }

                let counts_entry = counts.entry(cell_id).or_default();

                apply_encoded_call(self.config.stranded, encoded, counts_entry);
            }

            let position_meta = &chunk.positions[current_index];
            chunk_results.push(PositionData {
                contig_id: tid,
                pos: position_meta.pos,
                counts: counts.drain().collect(),
            });
        }

        // Return reader to pool with atomic size update
        {
            let mut pool = self.reader_pool.lock();
            pool.push(reader);
            self.reader_pool_size.fetch_add(1, Ordering::Release);
        }

        Ok((chunk_results, skipped_sites))
    }
}
