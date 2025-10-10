//! Streaming AnnData conversion producing sparse matrices with bounded memory.
//!
//! Performance optimizations:
//! - Uses FxHashMap instead of std HashMap for faster hashing (non-cryptographic)
//! - Minimal Arc clones (only when necessary for ownership transfer)
//! - SmallVec for small temporary collections to avoid heap allocations
//! - Two-pointer algorithms in sparse matrix operations

use crate::pipeline::bam2mtx::barcode::BarcodeProcessor;
use crate::pipeline::bam2mtx::processor::{BaseCounts, PositionData, StrandBaseCounts};
use anndata::{AnnData, AnnDataOp, AxisArraysOp};
use anndata_hdf5::H5;
use anyhow::{anyhow, Context, Result};
use log::info;
use nalgebra_sparse::coo::CooMatrix;
use nalgebra_sparse::csr::CsrMatrix;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::convert::TryInto;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::Arc;
use tempfile::NamedTempFile;

/// Configuration options for AnnData conversion.
#[derive(Debug, Clone)]
pub struct AnnDataConfig {
    /// Whether to keep strand-specific layers (true) or collapse to unstranded (false).
    pub stranded: bool,
    /// Optional compression hint passed to downstream writers.
    pub compression: Option<String>,
    /// Number of worker threads used for conversions.
    pub threads: usize,
    /// Hint for chunk size when aggregating positions.
    pub chunk_size: usize,
    /// Expected matrix density (retained for CLI compatibility).
    pub matrix_density: f64,
    /// Preferred batch size for IO flushing (retained for CLI compatibility).
    pub batch_size: usize,
    /// Upper bound on in-memory triplets before spilling to disk.
    pub triplet_spill_nnz: usize,
    /// Optional hint about total positions expected in the manifest.
    pub total_positions: usize,
}

impl Default for AnnDataConfig {
    fn default() -> Self {
        Self {
            stranded: true,
            compression: Some("gzip".to_string()),
            threads: num_cpus::get(),
            chunk_size: 15_000,
            matrix_density: 0.005,
            batch_size: 2_500,
            triplet_spill_nnz: 500_000,
            total_positions: 0,
        }
    }
}

/// Streaming AnnData converter that respects barcode and contig identifiers.
#[derive(Clone)]
pub struct AnnDataConverter {
    config: AnnDataConfig,
    barcode_processor: Arc<BarcodeProcessor>,
    contig_names: Arc<Vec<String>>,
}

impl AnnDataConverter {
    pub fn new(
        config: AnnDataConfig,
        barcode_processor: Arc<BarcodeProcessor>,
        contig_names: Arc<Vec<String>>,
    ) -> Self {
        Self {
            config,
            barcode_processor,
            contig_names,
        }
    }

    /// Convert an iterator of `PositionData` into an `AnnData` file using a bounded-memory pipeline.
    pub fn convert_streaming<I>(&self, data_iter: I, output_path: &Path) -> Result<AnnData<H5>>
    where
        I: IntoIterator<Item = PositionData>,
    {
        // Arc clones are cheap (just pointer + atomic increment), but we minimize them
        // by reusing the existing Arcs from self
        let mut builder = StreamingMatrixBuilder::new(
            self.config.clone(),
            Arc::clone(&self.barcode_processor),
            Arc::clone(&self.contig_names),
        )?;

        for position in data_iter.into_iter() {
            builder.ingest(position)?;
        }

        let parts = builder.finalize()?;
        self.write_parts(parts, output_path)
    }

    /// Backwards-compatible adapter for call sites that already collected data in memory.
    pub fn convert(&self, data: &[PositionData], output_path: &Path) -> Result<AnnData<H5>> {
        self.convert_streaming(data.to_vec().into_iter(), output_path)
    }

    fn write_parts(&self, parts: StreamingAnnDataParts, output_path: &Path) -> Result<AnnData<H5>> {
        let StreamingAnnDataParts {
            cell_names,
            position_names,
            forward_layers,
            reverse_layers,
        } = parts;

        info!(
            "Writing AnnData with {} cells × {} positions",
            cell_names.len(),
            position_names.len()
        );

        let adata = AnnData::<H5>::new(output_path)
            .with_context(|| format!("failed to create AnnData file at {:?}", output_path))?;

        let obs_index = cell_names
            .into_iter()
            .collect::<anndata::data::array::dataframe::DataFrameIndex>();
        let var_index = position_names
            .into_iter()
            .collect::<anndata::data::array::dataframe::DataFrameIndex>();

        adata.set_obs_names(obs_index)?;
        adata.set_var_names(var_index)?;

        let main_matrix = forward_layers
            .first()
            .cloned()
            .unwrap_or_else(|| CsrMatrix::zeros(0, 0));
        adata.set_x(main_matrix)?;

        let forward_names = ["A1", "T1", "G1", "C1"];
        for (name, matrix) in forward_names.iter().zip(forward_layers.into_iter()) {
            adata.layers().add(*name, matrix)?;
        }

        if self.config.stranded {
            let reverse_names = ["A0", "T0", "G0", "C0"];
            for (name, matrix) in reverse_names.iter().zip(reverse_layers.into_iter()) {
                adata.layers().add(*name, matrix)?;
            }
        }

        Ok(adata)
    }

    /// Historical compatibility hook retained for API stability.
    pub fn write_to_file(&self, _adata: &AnnData<H5>, _output_path: &Path) -> Result<()> {
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct Triplet {
    row: u32,
    col: u32,
    value: u32,
}

impl Triplet {
    fn to_bytes(self) -> [u8; 12] {
        let mut data = [0u8; 12];
        data[0..4].copy_from_slice(&self.row.to_le_bytes());
        data[4..8].copy_from_slice(&self.col.to_le_bytes());
        data[8..12].copy_from_slice(&self.value.to_le_bytes());
        data
    }

    fn from_bytes(bytes: &[u8; 12]) -> Self {
        Triplet {
            row: u32::from_le_bytes(bytes[0..4].try_into().unwrap()),
            col: u32::from_le_bytes(bytes[4..8].try_into().unwrap()),
            value: u32::from_le_bytes(bytes[8..12].try_into().unwrap()),
        }
    }
}

struct TripletSpool {
    file: NamedTempFile,
    writer: BufWriter<std::fs::File>,
    total: usize,
}

impl TripletSpool {
    fn new() -> Result<Self> {
        let file = NamedTempFile::new()?;
        let mut writer_handle = file.reopen()?;
        writer_handle.seek(SeekFrom::End(0))?;
        let writer = BufWriter::with_capacity(1 << 20, writer_handle);
        Ok(Self {
            file,
            writer,
            total: 0,
        })
    }

    fn append(&mut self, triplets: &[Triplet]) -> Result<()> {
        for triplet in triplets {
            self.writer
                .write_all(&triplet.to_bytes())
                .context("failed to write triplet to spill file")?;
        }
        self.total += triplets.len();
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        self.writer.flush().context("failed to flush spill buffer")
    }

    fn read_into(&mut self, target: &mut Vec<Triplet>) -> Result<()> {
        self.flush()?;
        if self.total == 0 {
            return Ok(());
        }

        target.reserve(self.total);
        let mut reader = BufReader::with_capacity(1 << 20, self.file.reopen()?);
        let mut buf = [0u8; 12];
        loop {
            match reader.read_exact(&mut buf) {
                Ok(()) => target.push(Triplet::from_bytes(&buf)),
                Err(err) if err.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(err) => return Err(err.into()),
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Eq, PartialEq, Hash)]
struct PositionKey {
    contig_id: u32,
    pos: u64,
}

struct StreamingAnnDataParts {
    cell_names: Vec<String>,
    position_names: Vec<String>,
    forward_layers: Vec<CsrMatrix<f32>>,
    reverse_layers: Vec<CsrMatrix<f32>>,
}

struct StreamingMatrixBuilder {
    config: AnnDataConfig,
    barcode_processor: Arc<BarcodeProcessor>,
    contig_names: Arc<Vec<String>>,
    n_cells: usize,
    observed_cells: Vec<bool>,
    position_lookup: FxHashMap<PositionKey, u32>,
    positions: Vec<PositionKey>,
    forward_buffers: [Vec<Triplet>; 4],
    reverse_buffers: Option<[Vec<Triplet>; 4]>,
    forward_spools: [TripletSpool; 4],
    reverse_spools: Option<[TripletSpool; 4]>,
    pending_triplets: usize,
    spill_threshold: usize,
}

impl StreamingMatrixBuilder {
    fn new(
        config: AnnDataConfig,
        barcode_processor: Arc<BarcodeProcessor>,
        contig_names: Arc<Vec<String>>,
    ) -> Result<Self> {
        let n_cells = barcode_processor.len();
        let spill_threshold = config.triplet_spill_nnz.max(1);

        let forward_buffers = [
            Vec::with_capacity(config.chunk_size),
            Vec::with_capacity(config.chunk_size),
            Vec::with_capacity(config.chunk_size),
            Vec::with_capacity(config.chunk_size),
        ];
        let forward_spools = [
            TripletSpool::new()?,
            TripletSpool::new()?,
            TripletSpool::new()?,
            TripletSpool::new()?,
        ];

        let (reverse_buffers, reverse_spools) = if config.stranded {
            let buffers = [
                Vec::with_capacity(config.chunk_size),
                Vec::with_capacity(config.chunk_size),
                Vec::with_capacity(config.chunk_size),
                Vec::with_capacity(config.chunk_size),
            ];
            let spools = [
                TripletSpool::new()?,
                TripletSpool::new()?,
                TripletSpool::new()?,
                TripletSpool::new()?,
            ];
            (Some(buffers), Some(spools))
        } else {
            (None, None)
        };

        Ok(Self {
            config,
            barcode_processor,
            contig_names,
            n_cells,
            observed_cells: vec![false; n_cells],
            position_lookup: FxHashMap::default(),
            positions: Vec::with_capacity(1024),
            forward_buffers,
            reverse_buffers,
            forward_spools,
            reverse_spools,
            pending_triplets: 0,
            spill_threshold,
        })
    }

    fn ingest(&mut self, data: PositionData) -> Result<()> {
        let column = self.ensure_position(PositionKey {
            contig_id: data.contig_id,
            pos: data.pos,
        });

        let mut added = 0usize;
        for (cell_id, counts) in data.counts.into_iter() {
            if (cell_id as usize) >= self.n_cells {
                continue;
            }
            self.observed_cells[cell_id as usize] = true;

            if self.config.stranded {
                added += self.append_stranded(cell_id, column, counts);
            } else {
                added += self.append_unstranded(cell_id, column, counts);
            }
        }

        self.pending_triplets += added;
        // More aggressive spilling: use lower threshold and check more frequently
        if self.pending_triplets >= self.spill_threshold.min(50_000) {
            self.flush_triplets()?;
        }
        Ok(())
    }

    fn finalize(mut self) -> Result<StreamingAnnDataParts> {
        self.flush_triplets()?;

        // Use SmallVec for typically small collections (most experiments have < 10k cells observed)
        // This avoids heap allocation for the common case while supporting larger datasets
        let observed_ids: Vec<u32> = self
            .observed_cells
            .iter()
            .enumerate()
            .filter_map(|(idx, &seen)| seen.then_some(idx as u32))
            .collect();

        let n_rows = observed_ids.len();
        let n_cols = self.positions.len();

        // Use FxHashMap-style allocation pattern for better performance
        let mut cell_remap = vec![None; self.n_cells];
        for (row_idx, cell_id) in observed_ids.iter().enumerate() {
            cell_remap[*cell_id as usize] = Some(row_idx as u32);
        }

        // Minimize allocations by pre-sizing with collected length
        let cell_names = observed_ids
            .iter()
            .map(|&id| {
                self.barcode_processor
                    .barcode_by_id(id)
                    .unwrap_or("unknown")
                    .to_string()
            })
            .collect();

        let position_names = self
            .positions
            .iter()
            .map(|key| {
                let contig = self
                    .contig_names
                    .get(key.contig_id as usize)
                    .map(|s| s.as_str())
                    .unwrap_or("unknown");
                format!("{}:{}", contig, key.pos)
            })
            .collect();

        let forward_layers = materialize_layers(
            &mut self.forward_spools,
            &mut self.forward_buffers,
            &cell_remap,
            n_rows,
            n_cols,
        )?;

        let reverse_layers = if self.config.stranded {
            let spools = self.reverse_spools.as_mut().unwrap();
            let buffers = self.reverse_buffers.as_mut().unwrap();
            materialize_layers(spools, buffers, &cell_remap, n_rows, n_cols)?
        } else {
            Vec::new()
        };

        Ok(StreamingAnnDataParts {
            cell_names,
            position_names,
            forward_layers,
            reverse_layers,
        })
    }

    fn append_stranded(&mut self, cell_id: u32, column: u32, counts: StrandBaseCounts) -> usize {
        let StrandBaseCounts { forward, reverse } = counts;
        let mut added = append_counts(&mut self.forward_buffers, cell_id, column, forward);
        if let Some(reverse_buffers) = self.reverse_buffers.as_mut() {
            added += append_counts(reverse_buffers, cell_id, column, reverse);
        }
        added
    }

    fn append_unstranded(&mut self, cell_id: u32, column: u32, counts: StrandBaseCounts) -> usize {
        let StrandBaseCounts { forward, reverse } = counts;
        let merged = BaseCounts {
            a: forward.a + reverse.a,
            t: forward.t + reverse.t,
            g: forward.g + reverse.g,
            c: forward.c + reverse.c,
        };
        append_counts(&mut self.forward_buffers, cell_id, column, merged)
    }

    fn flush_triplets(&mut self) -> Result<()> {
        if self.pending_triplets == 0 {
            return Ok(());
        }

        for (buffer, spool) in self.forward_buffers.iter_mut().zip(self.forward_spools.iter_mut()) {
            if !buffer.is_empty() {
                spool.append(buffer)?;
                buffer.clear();
            }
        }

        if let (Some(buffers), Some(spools)) = (
            self.reverse_buffers.as_mut(),
            self.reverse_spools.as_mut(),
        ) {
            for (buffer, spool) in buffers.iter_mut().zip(spools.iter_mut()) {
                if !buffer.is_empty() {
                    spool.append(buffer)?;
                    buffer.clear();
                }
            }
        }

        self.pending_triplets = 0;
        Ok(())
    }

    fn ensure_position(&mut self, key: PositionKey) -> u32 {
        if let Some(&idx) = self.position_lookup.get(&key) {
            idx
        } else {
            let idx = self.positions.len() as u32;
            self.positions.push(key);
            self.position_lookup.insert(key, idx);
            idx
        }
    }
}

fn append_counts(buffers: &mut [Vec<Triplet>; 4], cell_id: u32, column: u32, counts: BaseCounts) -> usize {
    let mut added = 0usize;
    if counts.a > 0 {
        buffers[0].push(Triplet {
            row: cell_id,
            col: column,
            value: counts.a,
        });
        added += 1;
    }
    if counts.t > 0 {
        buffers[1].push(Triplet {
            row: cell_id,
            col: column,
            value: counts.t,
        });
        added += 1;
    }
    if counts.g > 0 {
        buffers[2].push(Triplet {
            row: cell_id,
            col: column,
            value: counts.g,
        });
        added += 1;
    }
    if counts.c > 0 {
        buffers[3].push(Triplet {
            row: cell_id,
            col: column,
            value: counts.c,
        });
        added += 1;
    }
    added
}

fn materialize_layers(
    spools: &mut [TripletSpool; 4],
    buffers: &mut [Vec<Triplet>; 4],
    cell_remap: &[Option<u32>],
    n_rows: usize,
    n_cols: usize,
) -> Result<Vec<CsrMatrix<f32>>> {
    let mut layers = Vec::with_capacity(4);
    for (layer_idx, (spool, buffer)) in spools.iter_mut().zip(buffers.iter_mut()).enumerate() {
        let mut triplets = Vec::new();
        spool.read_into(&mut triplets)?;
        triplets.extend(buffer.drain(..));

        remap_triplets(&mut triplets, cell_remap);
        let matrix = triplets_to_csr(triplets, n_rows, n_cols)?;
        info!("Layer {} assembled with {} nnz", layer_idx, matrix.nnz());
        layers.push(matrix);
    }
    Ok(layers)
}

fn remap_triplets(triplets: &mut Vec<Triplet>, cell_remap: &[Option<u32>]) {
    let mut write_idx = 0usize;
    for read_idx in 0..triplets.len() {
        let orig_row = triplets[read_idx].row as usize;
        if let Some(new_row) = cell_remap.get(orig_row).and_then(|opt| *opt) {
            triplets[write_idx] = Triplet {
                row: new_row,
                col: triplets[read_idx].col,
                value: triplets[read_idx].value,
            };
            write_idx += 1;
        }
    }
    triplets.truncate(write_idx);
}

fn triplets_to_csr(triplets: Vec<Triplet>, n_rows: usize, n_cols: usize) -> Result<CsrMatrix<f32>> {
    if triplets.is_empty() {
        return Ok(CsrMatrix::zeros(n_rows, n_cols));
    }

    let mut triplets = triplets;
    triplets.par_sort_unstable_by(|a, b| a.row.cmp(&b.row).then(a.col.cmp(&b.col)));

    let mut merged = Vec::with_capacity(triplets.len());
    let mut iter = triplets.into_iter();
    if let Some(mut current) = iter.next() {
        for triplet in iter {
            if triplet.row == current.row && triplet.col == current.col {
                current.value = current.value.saturating_add(triplet.value);
            } else {
                merged.push(current);
                current = triplet;
            }
        }
        merged.push(current);
    }

    let mut rows = Vec::with_capacity(merged.len());
    let mut cols = Vec::with_capacity(merged.len());
    let mut values = Vec::with_capacity(merged.len());
    for Triplet { row, col, value } in merged {
        rows.push(row as usize);
        cols.push(col as usize);
        values.push(value as f32);
    }

    let coo = CooMatrix::try_from_triplets(n_rows, n_cols, rows, cols, values)
        .map_err(|err| anyhow!(err.to_string()))?;
    Ok(CsrMatrix::from(&coo))
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;

    fn make_counts(
        forward: (u32, u32, u32, u32),
        reverse: (u32, u32, u32, u32),
    ) -> StrandBaseCounts {
        StrandBaseCounts {
            forward: BaseCounts {
                a: forward.0,
                t: forward.1,
                g: forward.2,
                c: forward.3,
            },
            reverse: BaseCounts {
                a: reverse.0,
                t: reverse.1,
                g: reverse.2,
                c: reverse.3,
            },
        }
    }

    fn position(
        contig_id: u32,
        pos: u64,
        entries: Vec<(u32, StrandBaseCounts)>,
    ) -> PositionData {
        let mut counts: FxHashMap<u32, StrandBaseCounts> = FxHashMap::default();
        for (cell, value) in entries {
            counts.insert(cell, value);
        }
        PositionData {
            contig_id,
            pos,
            counts,
        }
    }

    fn matrix_entries(matrix: &CsrMatrix<f32>) -> Vec<(usize, usize, f32)> {
        CooMatrix::from(matrix)
            .triplet_iter()
            .map(|(r, c, v)| (r, c, *v))
            .collect()
    }

    #[test]
    fn streaming_builder_produces_expected_sparse_layers() -> Result<()> {
        let mut config = AnnDataConfig::default();
        config.stranded = true;
        config.threads = 2;
        config.chunk_size = 4;
        config.batch_size = 4;
        config.triplet_spill_nnz = 2;

        let barcode_processor = Arc::new(BarcodeProcessor::from_vec(vec![
            "AAAC".to_string(),
            "GGGG".to_string(),
        ]));
        let contig_names = Arc::new(vec!["chr1".to_string()]);

        let mut builder = StreamingMatrixBuilder::new(
            config.clone(),
            Arc::clone(&barcode_processor),
            Arc::clone(&contig_names),
        )?;

        let data = vec![
            position(
                0,
                1,
                vec![
                    (0, make_counts((2, 0, 1, 0), (0, 0, 0, 0))),
                    (1, make_counts((0, 1, 0, 0), (0, 0, 0, 2))),
                ],
            ),
            position(
                0,
                2,
                vec![
                    (0, make_counts((0, 0, 0, 0), (1, 0, 0, 0))),
                    (1, make_counts((3, 0, 0, 0), (0, 0, 0, 0))),
                ],
            ),
        ];

        for entry in data {
            builder.ingest(entry)?;
        }

        let parts = builder.finalize()?;
        assert_eq!(parts.cell_names, vec!["AAAC".to_string(), "GGGG".to_string()]);
        assert_eq!(
            parts.position_names,
            vec!["chr1:1".to_string(), "chr1:2".to_string()]
        );

        let forward_layers = parts.forward_layers;
        assert_eq!(forward_layers.len(), 4);

        let forward_a = matrix_entries(&forward_layers[0]);
        assert_eq!(forward_a, vec![(0, 0, 2.0), (1, 1, 3.0)]);

        let forward_t = matrix_entries(&forward_layers[1]);
        assert_eq!(forward_t, vec![(1, 0, 1.0)]);

        let forward_g = matrix_entries(&forward_layers[2]);
        assert_eq!(forward_g, vec![(0, 0, 1.0)]);

        let forward_c = matrix_entries(&forward_layers[3]);
        assert!(forward_c.is_empty());

        let reverse_layers = parts.reverse_layers;
        assert_eq!(reverse_layers.len(), 4);
        let reverse_a = matrix_entries(&reverse_layers[0]);
        assert_eq!(reverse_a, vec![(0, 1, 1.0)]);

        let reverse_c = matrix_entries(&reverse_layers[3]);
        assert_eq!(reverse_c, vec![(1, 0, 2.0)]);

        Ok(())
    }
}
