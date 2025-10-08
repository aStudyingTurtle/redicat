use anyhow::Result;
use parking_lot::Mutex;
use redicat_lib::bam2mtx::barcode::BarcodeProcessor;
use redicat_lib::bam2mtx::processor::{
    apply_encoded_call, decode_base, decode_cell_barcode, decode_umi, encode_call,
    BamProcessorConfig, PositionData, StrandBaseCounts, UMI_CONFLICT_CODE,
};
use rust_htslib::bam::{self, Read};
use rustc_hash::FxHashMap;
use std::path::PathBuf;
use std::sync::Arc;

use super::input::PositionChunk;

/// Chunk processor with aggressive reuse of BAM readers and pre-sized hash maps.
pub struct OptimizedChunkProcessor {
    bam_path: PathBuf,
    config: BamProcessorConfig,
    barcode_processor: Arc<BarcodeProcessor>,
    tid_lookup: FxHashMap<String, u32>,
    contig_names: Arc<Vec<String>>,
    count_capacity_hint: usize,
    umi_capacity_hint: usize,
    reader_pool: Mutex<Vec<bam::IndexedReader>>,
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
        })
    }

    pub fn contig_names(&self) -> Arc<Vec<String>> {
        Arc::clone(&self.contig_names)
    }

    pub fn process_chunk(&self, chunk: &PositionChunk) -> Result<Vec<PositionData>> {
        if chunk.positions.is_empty() {
            return Ok(Vec::new());
        }

        let mut reader = {
            let mut pool = self.reader_pool.lock();
            pool.pop()
        }
        .unwrap_or_else(|| {
            bam::IndexedReader::from_path(&self.bam_path)
                .unwrap_or_else(|e| panic!("Failed to open BAM {}: {}", self.bam_path.display(), e))
        });

        let chrom = &chunk.positions[0].chrom;
        let tid = match self.tid_lookup.get(chrom) {
            Some(&tid) => tid,
            None => {
                let mut pool = self.reader_pool.lock();
                pool.push(reader);
                return Ok(Vec::new());
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

        let mut chunk_results = Vec::with_capacity(chunk.positions.len());
        let mut position_index = 0usize;

    let mut counts: FxHashMap<u32, StrandBaseCounts> = FxHashMap::default();
        counts.reserve(self.count_capacity_hint);
    let mut umi_consensus: FxHashMap<(u32, String), u8> = FxHashMap::default();
        umi_consensus.reserve(self.umi_capacity_hint);

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

                let cell_id = match decode_cell_barcode(
                    &record,
                    self.config.cell_barcode_tag.as_bytes(),
                )? {
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

                if let Some(encoded) = encode_call(self.config.stranded, base, record.is_reverse())
                {
                    umi_consensus
                        .entry((cell_id, umi))
                        .and_modify(|existing| {
                            if *existing != encoded {
                                *existing = UMI_CONFLICT_CODE;
                            }
                        })
                        .or_insert(encoded);
                }
            }

            let current_index = position_index;
            position_index += 1;

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

        {
            let mut pool = self.reader_pool.lock();
            pool.push(reader);
        }

        Ok(chunk_results)
    }
}
