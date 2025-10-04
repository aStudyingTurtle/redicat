use anyhow::{Context, Result};
use noodles::{
    bam,
    bgzf,
    core::{self, Position, Region},
    sam,
};
use noodles::sam::alignment::record::cigar::op::Kind;
use parking_lot::Mutex;
use redicat_lib::bam2mtx::barcode::BarcodeProcessor;
use redicat_lib::bam2mtx::processor::{
    apply_encoded_call, decode_base, decode_cell_barcode, decode_umi, encode_call,
    BamProcessorConfig, PositionData, StrandBaseCounts, UMI_CONFLICT_CODE,
};
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::convert::TryFrom;
use std::fs::File;
use std::path::PathBuf;
use std::sync::Arc;

use super::input::PositionChunk;

type IndexedBamReader = bam::io::IndexedReader<bgzf::io::Reader<File>>;

struct ThreadSafeIndexedReader(IndexedBamReader);

impl ThreadSafeIndexedReader {
    fn new(reader: IndexedBamReader) -> Self {
        Self(reader)
    }

    fn into_inner(self) -> IndexedBamReader {
        self.0
    }
}

unsafe impl Send for ThreadSafeIndexedReader {}

/// Chunk processor with aggressive reuse of BAM readers and pre-sized hash maps.
pub struct OptimizedChunkProcessor {
    bam_path: PathBuf,
    config: BamProcessorConfig,
    barcode_processor: Arc<BarcodeProcessor>,
    tid_lookup: FxHashMap<String, usize>,
    header: Arc<sam::Header>,
    count_capacity_hint: usize,
    umi_capacity_hint: usize,
    reader_pool: Mutex<Vec<ThreadSafeIndexedReader>>,
}

impl OptimizedChunkProcessor {
    pub fn new(
        bam_path: PathBuf,
        config: BamProcessorConfig,
        barcode_processor: Arc<BarcodeProcessor>,
    ) -> Result<Self> {
        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&bam_path)
            .with_context(|| format!("Failed to open indexed BAM {}", bam_path.display()))?;
        let header: sam::Header = reader
            .read_header()
            .with_context(|| format!("Failed to read header from {}", bam_path.display()))?;

        let reference_sequences = header.reference_sequences();
        let mut tid_lookup = FxHashMap::with_capacity_and_hasher(reference_sequences.len(), Default::default());
        for (tid, (name, _)) in reference_sequences.iter().enumerate() {
            let contig = String::from_utf8_lossy(name.as_ref()).into_owned();
            tid_lookup.insert(contig, tid);
        }

        let header = Arc::new(header);

        let barcode_count = barcode_processor.len();
        let count_capacity_hint = barcode_count.clamp(64, 4096);
        let umi_capacity_hint = count_capacity_hint.saturating_mul(8);

    let mut reader_pool = Vec::new();
    reader_pool.push(ThreadSafeIndexedReader::new(reader));

        Ok(Self {
            bam_path,
            config,
            barcode_processor,
            tid_lookup,
            header,
            count_capacity_hint,
            umi_capacity_hint,
            reader_pool: Mutex::new(reader_pool),
        })
    }

    pub fn process_chunk(&self, chunk: &PositionChunk) -> Result<Vec<PositionData>> {
        if chunk.positions.is_empty() {
            return Ok(Vec::new());
        }

        let header = self.header.as_ref();

        let reader_opt = {
            let mut pool = self.reader_pool.lock();
            pool.pop()
        };

        let mut reader = match reader_opt {
            Some(reader) => reader.into_inner(),
            None => {
                let mut reader = bam::io::indexed_reader::Builder::default()
                    .build_from_path(&self.bam_path)
                    .with_context(|| {
                        format!("Failed to open indexed BAM {}", self.bam_path.display())
                    })?;
                reader
                    .read_header()
                    .with_context(|| {
                        format!("Failed to read header from {}", self.bam_path.display())
                    })?;
                reader
            }
        };

        let chrom = &chunk.positions[0].chrom;
        let tid = match self.tid_lookup.get(chrom) {
            Some(&tid) => tid,
            None => {
                let mut pool = self.reader_pool.lock();
                pool.push(ThreadSafeIndexedReader::new(reader));
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

        let start_pos = Position::try_from((fetch_start as usize) + 1)
            .with_context(|| format!("Invalid fetch start coordinate {}", fetch_start + 1))?;
        let end_pos = Position::try_from(fetch_end as usize)
            .with_context(|| format!("Invalid fetch end coordinate {}", fetch_end))?;
        let interval = core::region::Interval::from(start_pos..=end_pos);
        let region = Region::new(chrom.clone(), interval);

        let mut query = reader
            .query(header, &region)
            .with_context(|| {
                format!(
                    "Failed to query region {}:{}-{}",
                    chrom,
                    fetch_start + 1,
                    fetch_end
                )
            })?;

        let target_positions: Vec<u32> = chunk
            .positions
            .iter()
            .map(|p| p.pos.saturating_sub(1) as u32)
            .collect();

        let mut position_to_index =
            FxHashMap::with_capacity_and_hasher(target_positions.len(), Default::default());
        for (idx, pos) in target_positions.iter().enumerate() {
            position_to_index.insert(*pos, idx);
        }

        let min_target = i64::from(*target_positions.first().unwrap());
        let max_target = i64::from(*target_positions.last().unwrap());

        let mut position_umis: Vec<Option<FxHashMap<(String, String), u8>>> =
            vec![None; target_positions.len()];
        let mut visited = vec![false; target_positions.len()];
        let mut processed = 0u32;

        let cell_barcode_tag = self.config.cell_barcode_tag.as_bytes();
        let umi_tag = self.config.umi_tag.as_bytes();

        'records: while let Some(record_result) = query.next() {
            let record = record_result?;
            if record.flags().is_unmapped() {
                continue;
            }

            let read_tid = match record.reference_sequence_id().transpose()? {
                Some(id) => id,
                None => continue,
            };

            if read_tid != tid {
                continue;
            }

            let mapq = record.mapping_quality().map(u8::from).unwrap_or(0);
            let mapq_pass = mapq >= self.config.min_mapping_quality;

            let alignment_start = match record.alignment_start().transpose()? {
                Some(pos) => usize::from(pos).saturating_sub(1),
                None => continue,
            };

            let mut ref_pos = alignment_start as i64;
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
                            let genomic = ref_pos + i as i64;
                            if genomic < min_target || genomic > max_target {
                                continue;
                            }
                            if genomic < 0 {
                                continue;
                            }
                            let Ok(genomic_u32) = u32::try_from(genomic) else {
                                continue;
                            };
                            let Some(&pos_index) = position_to_index.get(&genomic_u32) else {
                                continue;
                            };

                            visited[pos_index] = true;

                            if !mapq_pass {
                                continue;
                            }

                            let qpos = read_pos + i;

                            let base_qual = record
                                .quality_scores()
                                .as_ref()
                                .get(qpos)
                                .copied()
                                .unwrap_or(0);
                            if base_qual < self.config.min_base_quality {
                                continue;
                            }

                            let base = decode_base(&record, qpos)?;
                            if base == 'N' {
                                continue;
                            }

                            let cell_barcode =
                                match decode_cell_barcode(&record, cell_barcode_tag)? {
                                    Some(barcode) => barcode,
                                    None => continue,
                                };
                            if !self.barcode_processor.is_valid(&cell_barcode) {
                                continue;
                            }

                            let umi = match decode_umi(&record, umi_tag)? {
                                Some(umi) => umi,
                                None => continue,
                            };

                            if let Some(encoded) = encode_call(
                                self.config.stranded,
                                base,
                                record.flags().is_reverse_complemented(),
                            ) {
                                let umi_map = position_umis[pos_index].get_or_insert_with(|| {
                                    let mut map = FxHashMap::default();
                                    map.reserve(self.umi_capacity_hint);
                                    map
                                });

                                umi_map
                                    .entry((cell_barcode, umi))
                                    .and_modify(|existing| {
                                        if *existing != encoded {
                                            *existing = UMI_CONFLICT_CODE;
                                        }
                                    })
                                    .or_insert(encoded);

                                processed = processed.saturating_add(1);
                                if processed >= self.config.max_depth {
                                    break 'records;
                                }
                            }
                        }

                        ref_pos += len as i64;
                        read_pos += len;
                    }
                    Kind::Deletion | Kind::Skip => {
                        for i in 0..len {
                            let genomic = ref_pos + i as i64;
                            if genomic < min_target || genomic > max_target {
                                continue;
                            }
                            if genomic < 0 {
                                continue;
                            }
                            if let Ok(genomic_u32) = u32::try_from(genomic) {
                                if let Some(&pos_index) = position_to_index.get(&genomic_u32) {
                                    visited[pos_index] = true;
                                }
                            }
                        }
                        ref_pos += len as i64;
                    }
                    Kind::Insertion => {
                        read_pos += len;
                    }
                    Kind::SoftClip => {
                        read_pos += len;
                    }
                    Kind::HardClip | Kind::Pad => {}
                }
            }
        }

        let mut chunk_results = Vec::with_capacity(chunk.positions.len());
        for (idx, position_meta) in chunk.positions.iter().enumerate() {
            if !visited[idx] {
                continue;
            }

            let consensus_map = position_umis[idx]
                .take()
                .unwrap_or_else(FxHashMap::default);

            let mut position_counts: HashMap<String, StrandBaseCounts> =
                HashMap::with_capacity(consensus_map.len().min(self.count_capacity_hint));

            for ((cell_barcode, _umi), encoded) in consensus_map {
                if encoded == UMI_CONFLICT_CODE {
                    continue;
                }

                let entry = position_counts
                    .entry(cell_barcode)
                    .or_insert_with(StrandBaseCounts::default);
                apply_encoded_call(self.config.stranded, encoded, entry);
            }

            chunk_results.push(PositionData {
                chrom: position_meta.chrom.clone(),
                pos: position_meta.pos,
                counts: position_counts,
            });
        }

        {
            let mut pool = self.reader_pool.lock();
            pool.push(ThreadSafeIndexedReader::new(reader));
        }

        Ok(chunk_results)
    }
}
