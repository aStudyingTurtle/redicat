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
use nalgebra_sparse::csr::CsrMatrix;
use rayon::{prelude::*, ThreadPool, ThreadPoolBuilder};
use rustc_hash::{FxHashMap, FxHashSet};
use std::cmp::Ordering;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::convert::TryInto;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::Arc;
use tempfile::NamedTempFile;

#[cfg(unix)]
use std::os::unix::fs::FileExt;
#[cfg(windows)]
use std::os::windows::fs::FileExt;

/// Lower bound for per-layer buffer allocations to keep ingestion efficient.
const MIN_STREAM_BUFFER_CAPACITY: usize = 512;

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

    /// Aggregate chunk-local results that were captured per worker thread and assemble the
    /// final sparse matrices in parallel before writing the `.h5ad` file.
    pub fn convert_parallel_chunks(
        &self,
        mut chunks: Vec<Vec<PositionData>>,
        output_path: &Path,
    ) -> Result<AnnData<H5>> {
        let total_positions: usize = chunks.iter().map(|c| c.len()).sum();
        let mut flattened = Vec::with_capacity(total_positions);
        for mut chunk in chunks.drain(..) {
            flattened.append(&mut chunk);
        }
        self.convert_parallel(&flattened, output_path)
    }

    /// Convert an in-memory collection of [`PositionData`] using the parallel sparse assembler.
    pub fn convert_parallel(
        &self,
        data: &[PositionData],
        output_path: &Path,
    ) -> Result<AnnData<H5>> {
        let parts = assemble_parallel_parts(
            &self.config,
            &self.barcode_processor,
            &self.contig_names,
            data,
        )?;
        self.write_parts(parts, output_path)
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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Triplet {
    row: u32,
    col: u32,
    value: u32,
}

impl Ord for Triplet {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.row.cmp(&other.row) {
            Ordering::Equal => self.col.cmp(&other.col),
            ord => ord,
        }
    }
}

impl PartialOrd for Triplet {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
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
    _file: NamedTempFile,
    reader: Arc<std::fs::File>,
    writer: BufWriter<std::fs::File>,
    runs: Vec<RunMeta>,
    bytes_written: u64,
    total: usize,
    sort_pool: Arc<ThreadPool>,
}

struct RunMeta {
    offset: u64,
    len: usize,
}

impl TripletSpool {
    fn new(sort_pool: Arc<ThreadPool>) -> Result<Self> {
        let file = NamedTempFile::new()?;
        let reader_handle = Arc::new(file.reopen()?);
        let mut writer_handle = file.reopen()?;
        writer_handle.seek(SeekFrom::End(0))?;
        let writer = BufWriter::with_capacity(1 << 20, writer_handle);
        Ok(Self {
            _file: file,
            reader: reader_handle,
            writer,
            runs: Vec::new(),
            bytes_written: 0,
            total: 0,
            sort_pool,
        })
    }

    fn append_from_buffer(&mut self, buffer: &mut Vec<Triplet>) -> Result<()> {
        if buffer.is_empty() {
            return Ok(());
        }

        const PARALLEL_SORT_THRESHOLD: usize = 32_768;
        if buffer.len() > 1 {
            if buffer.len() >= PARALLEL_SORT_THRESHOLD {
                self.sort_pool.install(|| buffer.par_sort_unstable());
            } else {
                buffer.sort_unstable();
            }
        }

        let mut write_len = 0usize;
        for idx in 0..buffer.len() {
            if write_len == 0 {
                buffer[write_len] = buffer[idx];
                write_len = 1;
                continue;
            }

            let prev = buffer[write_len - 1];
            let current = buffer[idx];
            if prev.row == current.row && prev.col == current.col {
                let merged = Triplet {
                    row: prev.row,
                    col: prev.col,
                    value: prev.value.saturating_add(current.value),
                };
                buffer[write_len - 1] = merged;
            } else {
                buffer[write_len] = current;
                write_len += 1;
            }
        }

        if write_len == 0 {
            buffer.clear();
            return Ok(());
        }

        let offset = self.bytes_written;
        let bytes = (write_len as u64) * 12;

        for triplet in buffer.iter().take(write_len) {
            self.writer
                .write_all(&triplet.to_bytes())
                .context("failed to write triplet to spill file")?;
        }

        self.runs.push(RunMeta {
            offset,
            len: write_len,
        });
        self.bytes_written += bytes;
        self.total += write_len;
        buffer.clear();
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        self.writer.flush().context("failed to flush spill buffer")
    }

    fn stream(&mut self) -> Result<TripletStream> {
        self.flush()?;
        TripletStream::new(self)
    }
}

struct TripletStream {
    runs: Vec<RunState>,
    heap: BinaryHeap<Reverse<HeapEntry>>,
}

struct RunState {
    file: Arc<std::fs::File>,
    next_offset: u64,
    remaining: usize,
    buffer: Vec<u8>,
    buffer_pos: usize,
    buffer_len: usize,
}

impl RunState {
    fn read_next(&mut self) -> Result<Option<Triplet>> {
        if self.buffer_pos >= self.buffer_len {
            self.fill_buffer()?;
            if self.buffer_len == 0 {
                return Ok(None);
            }
        }

        let chunk = &self.buffer[self.buffer_pos..self.buffer_pos + 12];
        self.buffer_pos += 12;

        let bytes: [u8; 12] = chunk.try_into().expect("triplet buffer slice length");
        Ok(Some(Triplet::from_bytes(&bytes)))
    }

    fn fill_buffer(&mut self) -> Result<()> {
        if self.remaining == 0 {
            self.buffer_len = 0;
            self.buffer_pos = 0;
            return Ok(());
        }

        const TRIPLETS_PER_READ: usize = 8192;
        const BYTES_PER_TRIPLET: usize = 12;

        let to_read = self.remaining.min(TRIPLETS_PER_READ);
        let byte_len = to_read * BYTES_PER_TRIPLET;

        if self.buffer.len() < byte_len {
            self.buffer.resize(byte_len, 0);
        }

        read_exact_shared(&self.file, &mut self.buffer[..byte_len], self.next_offset)
            .context("failed to read spill run data")?;

        self.next_offset += byte_len as u64;
        self.remaining -= to_read;
        self.buffer_pos = 0;
        self.buffer_len = byte_len;
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct HeapEntry {
    run_idx: usize,
    triplet: Triplet,
}

impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.triplet == other.triplet && self.run_idx == other.run_idx
    }
}

impl Eq for HeapEntry {}

impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.triplet.cmp(&other.triplet) {
            Ordering::Equal => self.run_idx.cmp(&other.run_idx),
            ord => ord,
        }
    }
}

impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl TripletStream {
    fn new(spool: &TripletSpool) -> Result<Self> {
        let mut runs = Vec::with_capacity(spool.runs.len());
        let mut heap = BinaryHeap::new();
        let file = Arc::clone(&spool.reader);

        for meta in spool.runs.iter() {
            if meta.len == 0 {
                continue;
            }

            let mut run_state = RunState {
                file: Arc::clone(&file),
                next_offset: meta.offset,
                remaining: meta.len,
                buffer: Vec::new(),
                buffer_pos: 0,
                buffer_len: 0,
            };

            if let Some(first) = run_state.read_next()? {
                let run_idx = runs.len();
                heap.push(Reverse(HeapEntry {
                    run_idx,
                    triplet: first,
                }));
                runs.push(run_state);
            }
        }

        Ok(Self { runs, heap })
    }

    fn next(&mut self) -> Result<Option<Triplet>> {
        let Some(Reverse(entry)) = self.heap.pop() else {
            return Ok(None);
        };

        let triplet = entry.triplet;
        if let Some(next_triplet) = self.runs[entry.run_idx].read_next()? {
            self.heap.push(Reverse(HeapEntry {
                run_idx: entry.run_idx,
                triplet: next_triplet,
            }));
        }

        Ok(Some(triplet))
    }
}

fn read_exact_shared(file: &std::fs::File, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
    #[cfg(unix)]
    {
        file.read_exact_at(buf, offset)
    }

    #[cfg(windows)]
    {
        let mut read = 0usize;
        while read < buf.len() {
            let slice = &mut buf[read..];
            let count = file.seek_read(slice, offset + read as u64)?;
            if count == 0 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    "unexpected EOF while reading spill file",
                ));
            }
            read += count;
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Eq, PartialEq, Hash, Debug)]
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

        let sort_threads = std::cmp::max(1, config.threads / 4).min(8);
        let sort_pool = Arc::new(
            ThreadPoolBuilder::new()
                .num_threads(sort_threads)
                .thread_name(|idx| format!("anndata-sort-{idx}"))
                .build()
                .context("failed to build AnnData sorting thread pool")?,
        );

        let layer_capacity = config.chunk_size.max(MIN_STREAM_BUFFER_CAPACITY);
        let forward_buffers = [
            Vec::with_capacity(layer_capacity),
            Vec::with_capacity(layer_capacity),
            Vec::with_capacity(layer_capacity),
            Vec::with_capacity(layer_capacity),
        ];
        let forward_spools = [
            TripletSpool::new(Arc::clone(&sort_pool))?,
            TripletSpool::new(Arc::clone(&sort_pool))?,
            TripletSpool::new(Arc::clone(&sort_pool))?,
            TripletSpool::new(Arc::clone(&sort_pool))?,
        ];

        let (reverse_buffers, reverse_spools) = if config.stranded {
            let buffers = [
                Vec::with_capacity(layer_capacity),
                Vec::with_capacity(layer_capacity),
                Vec::with_capacity(layer_capacity),
                Vec::with_capacity(layer_capacity),
            ];
            let spools = [
                TripletSpool::new(Arc::clone(&sort_pool))?,
                TripletSpool::new(Arc::clone(&sort_pool))?,
                TripletSpool::new(Arc::clone(&sort_pool))?,
                TripletSpool::new(Arc::clone(&sort_pool))?,
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

        for (buffer, spool) in self
            .forward_buffers
            .iter_mut()
            .zip(self.forward_spools.iter_mut())
        {
            if !buffer.is_empty() {
                spool.append_from_buffer(buffer)?;
            }
        }

        if let (Some(buffers), Some(spools)) =
            (self.reverse_buffers.as_mut(), self.reverse_spools.as_mut())
        {
            for (buffer, spool) in buffers.iter_mut().zip(spools.iter_mut()) {
                if !buffer.is_empty() {
                    spool.append_from_buffer(buffer)?;
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

fn append_counts(
    buffers: &mut [Vec<Triplet>; 4],
    cell_id: u32,
    column: u32,
    counts: BaseCounts,
) -> usize {
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
        if !buffer.is_empty() {
            spool.append_from_buffer(buffer)?;
        }

        let mut stream = spool.stream()?;
        let matrix = stream_to_csr(&mut stream, cell_remap, n_rows, n_cols)?;
        info!("Layer {} assembled with {} nnz", layer_idx, matrix.nnz());
        layers.push(matrix);
    }
    Ok(layers)
}

fn stream_to_csr(
    stream: &mut TripletStream,
    cell_remap: &[Option<u32>],
    n_rows: usize,
    n_cols: usize,
) -> Result<CsrMatrix<f32>> {
    if n_rows == 0 || n_cols == 0 {
        return Ok(CsrMatrix::zeros(n_rows, n_cols));
    }

    let mut column_indices: Vec<usize> = Vec::new();
    let mut values: Vec<f32> = Vec::new();
    let mut row_counts = vec![0usize; n_rows];

    let mut pending_key: Option<(usize, usize)> = None;
    let mut pending_value: u32 = 0;

    while let Some(triplet) = stream.next()? {
        let Some(remapped_row) = cell_remap.get(triplet.row as usize).and_then(|opt| *opt) else {
            continue;
        };

        let row = remapped_row as usize;
        let col = triplet.col as usize;
        if row >= n_rows || col >= n_cols {
            continue;
        }

        if let Some((prev_row, prev_col)) = pending_key {
            if prev_row == row && prev_col == col {
                pending_value = pending_value.saturating_add(triplet.value);
                continue;
            } else {
                row_counts[prev_row] = row_counts[prev_row].saturating_add(1);
                column_indices.push(prev_col);
                values.push(pending_value as f32);
            }
        }

        pending_key = Some((row, col));
        pending_value = triplet.value;
    }

    if let Some((row, col)) = pending_key {
        row_counts[row] = row_counts[row].saturating_add(1);
        column_indices.push(col);
        values.push(pending_value as f32);
    }

    let mut row_offsets = Vec::with_capacity(n_rows + 1);
    row_offsets.push(0);
    let mut nnz = 0usize;
    for count in row_counts.into_iter() {
        nnz = nnz.saturating_add(count);
        row_offsets.push(nnz);
    }

    debug_assert_eq!(nnz, column_indices.len());
    debug_assert_eq!(column_indices.len(), values.len());

    CsrMatrix::try_from_csr_data(n_rows, n_cols, row_offsets, column_indices, values)
        .map_err(|err| anyhow!(err.to_string()))
}

fn assemble_parallel_parts(
    config: &AnnDataConfig,
    barcode_processor: &Arc<BarcodeProcessor>,
    contig_names: &Arc<Vec<String>>,
    data: &[PositionData],
) -> Result<StreamingAnnDataParts> {
    if data.is_empty() {
        let mut forward_layers = Vec::with_capacity(4);
        for _ in 0..4 {
            forward_layers.push(CsrMatrix::zeros(0, 0));
        }
        let reverse_layers = if config.stranded {
            let mut layers = Vec::with_capacity(4);
            for _ in 0..4 {
                layers.push(CsrMatrix::zeros(0, 0));
            }
            layers
        } else {
            Vec::new()
        };

        return Ok(StreamingAnnDataParts {
            cell_names: Vec::new(),
            position_names: Vec::new(),
            forward_layers,
            reverse_layers,
        });
    }

    let observed_cells: FxHashSet<u32> = data
        .par_iter()
        .fold(
            || FxHashSet::default(),
            |mut acc, record| {
                acc.extend(record.counts.keys().copied());
                acc
            },
        )
        .reduce(|| FxHashSet::default(), |mut left, right| {
            left.extend(right.into_iter());
            left
        });

    let mut observed: Vec<u32> = observed_cells.into_iter().collect();
    observed.sort_unstable();

    let n_cells_total = barcode_processor.len();
    let mut cell_remap = vec![None; n_cells_total];
    for (row_idx, cell_id) in observed.iter().enumerate() {
        if let Some(slot) = cell_remap.get_mut(*cell_id as usize) {
            *slot = Some(row_idx as u32);
        }
    }

    let cell_names = observed
        .iter()
        .map(|&id| {
            barcode_processor
                .barcode_by_id(id)
                .unwrap_or("unknown")
                .to_string()
        })
        .collect();

    let mut position_keys: Vec<PositionKey> = data
        .iter()
        .map(|p| PositionKey {
            contig_id: p.contig_id,
            pos: p.pos,
        })
        .collect();

    position_keys.par_sort_unstable_by(|a, b| match a.contig_id.cmp(&b.contig_id) {
        Ordering::Equal => a.pos.cmp(&b.pos),
        ord => ord,
    });
    position_keys.dedup();

    let position_names = position_keys
        .iter()
        .map(|key| {
            let contig = contig_names
                .get(key.contig_id as usize)
                .map(|s| s.as_str())
                .unwrap_or("unknown");
            format!("{}:{}", contig, key.pos)
        })
        .collect();

    let mut position_lookup = FxHashMap::with_capacity_and_hasher(position_keys.len(), Default::default());
    for (idx, key) in position_keys.iter().enumerate() {
        position_lookup.insert(*key, idx as u32);
    }

    let n_rows = observed.len();
    let n_cols = position_keys.len();
    let capacity_hint = config.chunk_size.max(MIN_STREAM_BUFFER_CAPACITY);

    let accumulator = data
        .par_iter()
        .try_fold(
            || ParallelLayerAccumulator::new(config.stranded, capacity_hint),
            |mut acc, position| -> Result<ParallelLayerAccumulator> {
                let key = PositionKey {
                    contig_id: position.contig_id,
                    pos: position.pos,
                };
                let column = *position_lookup
                    .get(&key)
                    .ok_or_else(|| anyhow!("position {:?} missing from lookup", key))?;
                acc.ingest(position, column, &cell_remap, config.stranded);
                Ok(acc)
            },
        )
        .try_reduce(
            || ParallelLayerAccumulator::new(config.stranded, capacity_hint),
            |mut left, right| -> Result<ParallelLayerAccumulator> {
                left.merge(right);
                Ok(left)
            },
        )?;

    let (forward_layers, reverse_layers) =
        accumulator.finish(n_rows, n_cols, config.stranded)?;

    Ok(StreamingAnnDataParts {
        cell_names,
        position_names,
        forward_layers,
        reverse_layers,
    })
}

struct ParallelLayerAccumulator {
    forward: [Vec<Triplet>; 4],
    reverse: Option<[Vec<Triplet>; 4]>,
}

impl ParallelLayerAccumulator {
    fn new(stranded: bool, capacity_hint: usize) -> Self {
        let forward = [
            Vec::with_capacity(capacity_hint),
            Vec::with_capacity(capacity_hint),
            Vec::with_capacity(capacity_hint),
            Vec::with_capacity(capacity_hint),
        ];
        let reverse = if stranded {
            Some([
                Vec::with_capacity(capacity_hint),
                Vec::with_capacity(capacity_hint),
                Vec::with_capacity(capacity_hint),
                Vec::with_capacity(capacity_hint),
            ])
        } else {
            None
        };

        Self { forward, reverse }
    }

    fn ingest(
        &mut self,
        position: &PositionData,
        column: u32,
        cell_remap: &[Option<u32>],
        stranded: bool,
    ) {
        for (&cell_id, counts) in position.counts.iter() {
            let Some(row_opt) = cell_remap.get(cell_id as usize) else {
                continue;
            };
            let Some(row) = row_opt else {
                continue;
            };

            if stranded {
                self.push_stranded(*row, column, counts);
            } else {
                self.push_unstranded(*row, column, counts);
            }
        }
    }

    fn push_stranded(&mut self, row: u32, column: u32, counts: &StrandBaseCounts) {
        push_counts(&mut self.forward, row, column, &counts.forward);
        if let Some(reverse) = self.reverse.as_mut() {
            push_counts(reverse, row, column, &counts.reverse);
        }
    }

    fn push_unstranded(&mut self, row: u32, column: u32, counts: &StrandBaseCounts) {
        let merged = BaseCounts {
            a: counts.forward.a + counts.reverse.a,
            t: counts.forward.t + counts.reverse.t,
            g: counts.forward.g + counts.reverse.g,
            c: counts.forward.c + counts.reverse.c,
        };
        push_counts(&mut self.forward, row, column, &merged);
    }

    fn merge(&mut self, other: ParallelLayerAccumulator) {
        let [o_a, o_t, o_g, o_c] = other.forward;
        self.forward[0].extend(o_a);
        self.forward[1].extend(o_t);
        self.forward[2].extend(o_g);
        self.forward[3].extend(o_c);

        if let (Some(dst), Some([r_a, r_t, r_g, r_c])) = (self.reverse.as_mut(), other.reverse) {
            dst[0].extend(r_a);
            dst[1].extend(r_t);
            dst[2].extend(r_g);
            dst[3].extend(r_c);
        }
    }

    fn finish(
        self,
        n_rows: usize,
        n_cols: usize,
        stranded: bool,
    ) -> Result<(Vec<CsrMatrix<f32>>, Vec<CsrMatrix<f32>>)> {
        let [f_a, f_t, f_g, f_c] = self.forward;
        let forward_layers = vec![
            triplets_to_csr(f_a, n_rows, n_cols)?,
            triplets_to_csr(f_t, n_rows, n_cols)?,
            triplets_to_csr(f_g, n_rows, n_cols)?,
            triplets_to_csr(f_c, n_rows, n_cols)?,
        ];

        let reverse_layers = if stranded {
            if let Some([r_a, r_t, r_g, r_c]) = self.reverse {
                vec![
                    triplets_to_csr(r_a, n_rows, n_cols)?,
                    triplets_to_csr(r_t, n_rows, n_cols)?,
                    triplets_to_csr(r_g, n_rows, n_cols)?,
                    triplets_to_csr(r_c, n_rows, n_cols)?,
                ]
            } else {
                vec![
                    CsrMatrix::zeros(n_rows, n_cols),
                    CsrMatrix::zeros(n_rows, n_cols),
                    CsrMatrix::zeros(n_rows, n_cols),
                    CsrMatrix::zeros(n_rows, n_cols),
                ]
            }
        } else {
            Vec::new()
        };

        Ok((forward_layers, reverse_layers))
    }
}

fn push_counts(buffers: &mut [Vec<Triplet>; 4], row: u32, column: u32, counts: &BaseCounts) {
    if counts.a > 0 {
        buffers[0].push(Triplet {
            row,
            col: column,
            value: counts.a,
        });
    }
    if counts.t > 0 {
        buffers[1].push(Triplet {
            row,
            col: column,
            value: counts.t,
        });
    }
    if counts.g > 0 {
        buffers[2].push(Triplet {
            row,
            col: column,
            value: counts.g,
        });
    }
    if counts.c > 0 {
        buffers[3].push(Triplet {
            row,
            col: column,
            value: counts.c,
        });
    }
}

fn triplets_to_csr(
    mut triplets: Vec<Triplet>,
    n_rows: usize,
    n_cols: usize,
) -> Result<CsrMatrix<f32>> {
    if n_rows == 0 || n_cols == 0 {
        return Ok(CsrMatrix::zeros(n_rows, n_cols));
    }

    if triplets.is_empty() {
        return Ok(CsrMatrix::zeros(n_rows, n_cols));
    }

    if triplets.len() > 1 {
        triplets.par_sort_unstable();
    }

    let mut merged: Vec<Triplet> = Vec::with_capacity(triplets.len());
    let mut iter = triplets.into_iter();
    if let Some(mut current) = iter.next() {
        for entry in iter {
            if current.row == entry.row && current.col == entry.col {
                current.value = current.value.saturating_add(entry.value);
            } else {
                merged.push(current);
                current = entry;
            }
        }
        merged.push(current);
    }

    let mut column_indices: Vec<usize> = Vec::with_capacity(merged.len());
    let mut values: Vec<f32> = Vec::with_capacity(merged.len());
    let mut row_counts = vec![0usize; n_rows];

    for triplet in merged.into_iter() {
        let row = triplet.row as usize;
        let col = triplet.col as usize;
        if row >= n_rows || col >= n_cols {
            continue;
        }
        row_counts[row] = row_counts[row].saturating_add(1);
        column_indices.push(col);
        values.push(triplet.value as f32);
    }

    let mut row_offsets = Vec::with_capacity(n_rows + 1);
    row_offsets.push(0);
    let mut nnz = 0usize;
    for count in row_counts.into_iter() {
        nnz = nnz.saturating_add(count);
        row_offsets.push(nnz);
    }

    CsrMatrix::try_from_csr_data(n_rows, n_cols, row_offsets, column_indices, values)
        .map_err(|err| anyhow!(err.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Context;
    use anyhow::Result;
    use std::sync::Arc;

    fn test_sort_pool() -> Result<Arc<ThreadPool>> {
        Ok(Arc::new(
            ThreadPoolBuilder::new()
                .num_threads(1)
                .build()
                .context("failed to build test sorting pool")?,
        ))
    }

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

    fn position(contig_id: u32, pos: u64, entries: Vec<(u32, StrandBaseCounts)>) -> PositionData {
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
        let mut entries = Vec::new();
        for row in 0..matrix.nrows() {
            let view = matrix.row(row);
            for (&col, &value) in view.col_indices().iter().zip(view.values()) {
                entries.push((row, col, value));
            }
        }
        entries
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
        assert_eq!(
            parts.cell_names,
            vec!["AAAC".to_string(), "GGGG".to_string()]
        );
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

    #[test]
    fn parallel_assembler_matches_streaming_builder() -> Result<()> {
        let mut config = AnnDataConfig::default();
        config.stranded = true;
        config.threads = 2;
        config.chunk_size = 4;
        config.batch_size = 4;
        config.triplet_spill_nnz = 8;

        let barcode_processor = Arc::new(BarcodeProcessor::from_vec(vec![
            "AAAC".to_string(),
            "GGGG".to_string(),
        ]));
        let contig_names = Arc::new(vec!["chr1".to_string()]);

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

        let mut builder = StreamingMatrixBuilder::new(
            config.clone(),
            Arc::clone(&barcode_processor),
            Arc::clone(&contig_names),
        )?;
        for entry in data.clone() {
            builder.ingest(entry)?;
        }
        let streaming_parts = builder.finalize()?;

        let parallel_parts = assemble_parallel_parts(
            &config,
            &barcode_processor,
            &contig_names,
            &data,
        )?;

        assert_eq!(streaming_parts.cell_names, parallel_parts.cell_names);
        assert_eq!(streaming_parts.position_names, parallel_parts.position_names);

        let streaming_forward: Vec<Vec<(usize, usize, f32)>> = streaming_parts
            .forward_layers
            .iter()
            .map(matrix_entries)
            .collect();
        let parallel_forward: Vec<Vec<(usize, usize, f32)>> = parallel_parts
            .forward_layers
            .iter()
            .map(matrix_entries)
            .collect();
        assert_eq!(streaming_forward, parallel_forward);

        let streaming_reverse: Vec<Vec<(usize, usize, f32)>> = streaming_parts
            .reverse_layers
            .iter()
            .map(matrix_entries)
            .collect();
        let parallel_reverse: Vec<Vec<(usize, usize, f32)>> = parallel_parts
            .reverse_layers
            .iter()
            .map(matrix_entries)
            .collect();
        assert_eq!(streaming_reverse, parallel_reverse);

        Ok(())
    }

    #[test]
    fn triplet_spool_merges_runs_in_order() -> Result<()> {
        let sort_pool = test_sort_pool()?;
        let mut spool = TripletSpool::new(sort_pool)?;

        let mut run_a = vec![
            Triplet {
                row: 2,
                col: 5,
                value: 1,
            },
            Triplet {
                row: 1,
                col: 3,
                value: 2,
            },
            Triplet {
                row: 2,
                col: 5,
                value: 4,
            },
        ];
        spool.append_from_buffer(&mut run_a)?;

        let mut run_b = vec![
            Triplet {
                row: 0,
                col: 0,
                value: 7,
            },
            Triplet {
                row: 1,
                col: 3,
                value: 5,
            },
            Triplet {
                row: 3,
                col: 1,
                value: 6,
            },
        ];
        spool.append_from_buffer(&mut run_b)?;

        let mut stream = spool.stream()?;
        let identity_remap = (0..4).map(|idx| Some(idx as u32)).collect::<Vec<_>>();
        let matrix = stream_to_csr(&mut stream, &identity_remap, 4, 6)?;
        let mut observed = matrix_entries(&matrix);
        observed.sort_unstable_by_key(|&(row, col, _)| (row, col));

        assert_eq!(
            observed,
            vec![(0, 0, 7.0), (1, 3, 7.0), (2, 5, 5.0), (3, 1, 6.0),]
        );

        Ok(())
    }
}
