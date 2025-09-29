use anyhow::{anyhow, Result};
use parking_lot::Mutex;
use redicat_lib::bam2mtx::barcode::BarcodeProcessor;
use redicat_lib::bam2mtx::processor::{
    apply_encoded_call, encode_call, BamProcessorConfig, PositionData, StrandBaseCounts,
    UMI_CONFLICT_CODE,
};
use rust_htslib::bam::{self, Read};
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use super::input::PositionChunk;

/// Chunk processor with aggressive reuse of BAM readers and pre-sized hash maps.
pub struct OptimizedChunkProcessor {
    bam_path: PathBuf,
    config: BamProcessorConfig,
    barcode_processor: Arc<BarcodeProcessor>,
    tid_lookup: FxHashMap<String, u32>,
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

        for tid in 0..header.target_count() {
            if let Ok(name) = std::str::from_utf8(header.tid2name(tid)) {
                tid_lookup.insert(name.to_string(), tid);
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
            count_capacity_hint,
            umi_capacity_hint,
            reader_pool: Mutex::new(Vec::new()),
        })
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

        let mut counts: FxHashMap<String, StrandBaseCounts> = FxHashMap::default();
        counts.reserve(self.count_capacity_hint);
        let mut umi_consensus: FxHashMap<(String, String), u8> = FxHashMap::default();
        umi_consensus.reserve(self.umi_capacity_hint);

        for pileup in reader.pileup() {
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
            let mut processed: u32 = 0;
            let mut truncated = false;

            for alignment in pileup.alignments() {
                let record = alignment.record();

                if record.mapq() < self.config.min_mapping_quality {
                    continue;
                }

                let qpos = match alignment.qpos() {
                    Some(q) => q,
                    None => continue,
                };

                processed = processed.saturating_add(1);
                if processed > self.config.max_depth {
                    truncated = true;
                    break;
                }

                let base_qual = record.qual().get(qpos).copied().unwrap_or(0);
                if base_qual < self.config.min_base_quality {
                    continue;
                }

                let base = self.get_base_at_position(&record, Some(qpos))?;
                if base == 'N' {
                    continue;
                }

                let cell_barcode = self.get_cell_barcode(&record)?;
                if !self.barcode_processor.is_valid(&cell_barcode) {
                    continue;
                }

                let umi = self.get_umi(&record)?;
                if umi == "-" {
                    continue;
                }

                if let Some(encoded) = encode_call(self.config.stranded, base, record.is_reverse())
                {
                    umi_consensus
                        .entry((cell_barcode, umi))
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

            if truncated && self.config.skip_max_depth {
                continue;
            }

            for ((cell_barcode, _umi), encoded) in umi_consensus.drain() {
                if encoded == UMI_CONFLICT_CODE {
                    continue;
                }

                let counts_entry = counts
                    .entry(cell_barcode)
                    .or_insert_with(StrandBaseCounts::default);

                apply_encoded_call(self.config.stranded, encoded, counts_entry);
            }

            let mut position_counts: HashMap<String, StrandBaseCounts> =
                HashMap::with_capacity(counts.len());
            position_counts.extend(counts.drain());

            let position_meta = &chunk.positions[current_index];
            chunk_results.push(PositionData {
                chrom: position_meta.chrom.clone(),
                pos: position_meta.pos,
                counts: position_counts,
            });
        }

        {
            let mut pool = self.reader_pool.lock();
            pool.push(reader);
        }

        Ok(chunk_results)
    }

    fn get_cell_barcode(&self, record: &bam::Record) -> Result<String> {
        let tag = self.config.cell_barcode_tag.as_bytes();
        match record.aux(tag) {
            Ok(bam::record::Aux::String(s)) => {
                let clean_barcode = s.split('-').next().unwrap_or(s);
                Ok(clean_barcode.to_string())
            }
            Ok(bam::record::Aux::ArrayU8(arr)) => {
                let data: Vec<u8> = arr.iter().collect();
                match std::str::from_utf8(&data) {
                    Ok(barcode) => {
                        let clean_barcode = barcode.split('-').next().unwrap_or(barcode);
                        Ok(clean_barcode.to_string())
                    }
                    Err(_) => Ok("-".to_string()),
                }
            }
            _ => Ok("-".to_string()),
        }
    }

    fn get_umi(&self, record: &bam::Record) -> Result<String> {
        let tag = self.config.umi_tag.as_bytes();
        match record.aux(tag) {
            Ok(bam::record::Aux::String(s)) => Ok(s.to_string()),
            Ok(bam::record::Aux::ArrayU8(arr)) => {
                let bytes: Vec<u8> = arr.iter().collect();
                match std::str::from_utf8(&bytes) {
                    Ok(umi) => Ok(umi.to_string()),
                    Err(_) => Ok("-".to_string()),
                }
            }
            _ => Ok("-".to_string()),
        }
    }

    fn get_base_at_position(&self, record: &bam::Record, qpos: Option<usize>) -> Result<char> {
        let qpos = qpos.ok_or_else(|| anyhow!("Invalid query position"))?;
        let seq = record.seq();
        let base = seq.as_bytes()[qpos];

        let base = match base {
            b'A' | b'a' => 'A',
            b'T' | b't' => 'T',
            b'G' | b'g' => 'G',
            b'C' | b'c' => 'C',
            _ => 'N',
        };
        Ok(base)
    }
}
