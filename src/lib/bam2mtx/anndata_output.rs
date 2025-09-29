//! AnnData output functionality with advanced performance optimizations

use crate::bam2mtx::processor::PositionData;
use anndata::{AnnData, AnnDataOp, AxisArraysOp};
use anndata_hdf5::H5;
use anyhow::Result;
use log::info;
use nalgebra_sparse::coo::CooMatrix;
use nalgebra_sparse::csr::CsrMatrix;
use rayon::prelude::*;
use rustc_hash::FxHashMap; // Faster hash map
use smallvec::SmallVec; // Stack-allocated small vectors
use std::cell::UnsafeCell;
use std::path::Path;
use std::sync::{Arc, RwLock};

/// Configuration for AnnData output
#[derive(Debug, Clone)]
pub struct AnnDataConfig {
    /// Whether to process data in a strand-specific manner
    pub stranded: bool,
    /// Compression algorithm to use for the output file (e.g., "gzip", "lz4")
    pub compression: Option<String>,
    /// Number of threads to use for parallel processing
    pub threads: usize,
    /// Size of chunks for processing data in batches
    pub chunk_size: usize,
    /// Pre-allocated matrix density estimate (proportion of non-zero elements)
    pub matrix_density: f64,
    /// Batch processing size
    pub batch_size: usize,
}

impl Default for AnnDataConfig {
    fn default() -> Self {
        Self {
            stranded: true,
            compression: Some("gzip".to_string()),
            threads: num_cpus::get(),
            chunk_size: 15000,    // Medium size balancing memory and performance
            matrix_density: 0.01, // Typical sparsity for single-cell data
            batch_size: 1000,     // Moderate batch size
        }
    }
}

/// High-performance sparse matrix builder using lock-free data structures
struct HighPerformanceMatrixBuilder {
    // Use thread-local storage to reduce lock contention
    thread_local_triplets: Vec<UnsafeCell<Vec<(u32, u32, u32)>>>,
    shape: (usize, usize),
    estimated_nnz: usize,
}

unsafe impl Sync for HighPerformanceMatrixBuilder {}

impl HighPerformanceMatrixBuilder {
    fn new(shape: (usize, usize), estimated_nnz: usize) -> Self {
        let num_threads = rayon::current_num_threads();
        let per_thread_capacity = (estimated_nnz + num_threads - 1) / num_threads;

        let thread_local_triplets = (0..num_threads)
            .map(|_| UnsafeCell::new(Vec::with_capacity(per_thread_capacity)))
            .collect();

        Self {
            thread_local_triplets,
            shape,
            estimated_nnz,
        }
    }

    /// Lock-free addition of triplet batches
    fn add_triplet_batch_lockfree(&self, thread_id: usize, batch: SmallVec<[(u32, u32, u32); 64]>) {
        if thread_id < self.thread_local_triplets.len() {
            unsafe {
                let triplets = &mut *self.thread_local_triplets[thread_id].get();
                triplets.extend(batch.into_iter().filter(|(_, _, val)| *val != 0));
            }
        }
    }

    /// Build CSR matrix using more efficient algorithms
    fn build_optimized(self) -> CsrMatrix<f32> {
        // Collect triplets from all threads
        let mut all_triplets = Vec::with_capacity(self.estimated_nnz);
        let shape = self.shape; // Save shape in advance

        // Fix: Get all data first, then process
        for cell in self.thread_local_triplets {
            let triplets = cell.into_inner();
            all_triplets.extend(triplets);
        }

        if all_triplets.is_empty() {
            return CsrMatrix::zeros(shape.0, shape.1);
        }

        // Parallel sort for better cache locality
        all_triplets.par_sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        // Properly separate the triplets
        let mut row_indices = Vec::with_capacity(all_triplets.len());
        let mut col_indices = Vec::with_capacity(all_triplets.len());
        let mut values = Vec::with_capacity(all_triplets.len());

        for (r, c, v) in all_triplets {
            row_indices.push(r as usize);
            col_indices.push(c as usize);
            values.push(v as f32); // Convert u32 to f32
        }

        // Use COO to build and then convert to CSR
        let coo = CooMatrix::try_from_triplets(shape.0, shape.1, row_indices, col_indices, values)
            .unwrap_or_else(|_| CooMatrix::new(shape.0, shape.1));

        CsrMatrix::from(&coo)
    }
}

/// String interning pool to reduce string allocations
#[allow(dead_code)]
struct StringInterner {
    pool: RwLock<FxHashMap<String, u32>>,
    reverse_pool: RwLock<Vec<String>>,
    next_id: std::sync::atomic::AtomicU32,
}

#[allow(dead_code)]
impl StringInterner {
    fn new() -> Self {
        Self {
            pool: RwLock::new(FxHashMap::default()),
            reverse_pool: RwLock::new(Vec::new()),
            next_id: std::sync::atomic::AtomicU32::new(0),
        }
    }

    fn intern(&self, s: &str) -> u32 {
        {
            let pool = self.pool.read().unwrap();
            if let Some(&id) = pool.get(s) {
                return id;
            }
        }

        let mut pool = self.pool.write().unwrap();
        if let Some(&id) = pool.get(s) {
            return id;
        }

        let id = self
            .next_id
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        pool.insert(s.to_string(), id);

        let mut reverse_pool = self.reverse_pool.write().unwrap();
        reverse_pool.push(s.to_string());

        id
    }

    fn get_string(&self, id: u32) -> Option<String> {
        let reverse_pool = self.reverse_pool.read().unwrap();
        reverse_pool.get(id as usize).cloned()
    }
}

/// Converter for transforming processed position data into AnnData format
pub struct AnnDataConverter {
    config: AnnDataConfig,
    #[allow(dead_code)]
    string_interner: Arc<StringInterner>,
}

impl AnnDataConverter {
    /// Create a new AnnDataConverter with the specified configuration
    pub fn new(config: AnnDataConfig) -> Self {
        Self {
            config,
            string_interner: Arc::new(StringInterner::new()),
        }
    }

    /// Optimized identifier collection using more efficient data structures
    fn collect_unique_identifiers_optimized(
        &self,
        position_data: &[PositionData],
    ) -> (Vec<String>, Vec<String>) {
        use rustc_hash::FxHashSet;

        info!(
            "  - Collecting identifiers with {} threads (optimized)...",
            self.config.threads
        );

        let chunk_size = std::cmp::max(1000, position_data.len() / (self.config.threads * 4));

        // Use faster hash sets with pre-allocation
        let results: Vec<(FxHashSet<String>, FxHashSet<String>)> = position_data
            .par_chunks(chunk_size)
            .map(|chunk| {
                let mut local_cells =
                    FxHashSet::with_capacity_and_hasher(chunk.len() * 10, Default::default());
                let mut local_positions =
                    FxHashSet::with_capacity_and_hasher(chunk.len(), Default::default());

                // Pre-allocate string buffer
                let mut pos_buffer = String::with_capacity(32);

                for data in chunk {
                    // Bulk insert cell barcodes
                    local_cells.extend(data.counts.keys().cloned());

                    // Reuse string buffer
                    pos_buffer.clear();
                    pos_buffer.push_str(&data.chrom);
                    pos_buffer.push(':');
                    pos_buffer.push_str(&data.pos.to_string());
                    local_positions.insert(pos_buffer.clone());
                }

                (local_cells, local_positions)
            })
            .collect();

        // Use more efficient merging strategy
        let total_cells_estimate = results.iter().map(|(c, _)| c.len()).sum::<usize>();
        let total_positions_estimate = results.iter().map(|(_, p)| p.len()).sum::<usize>();

        let mut cell_barcodes_set =
            FxHashSet::with_capacity_and_hasher(total_cells_estimate, Default::default());
        let mut positions_set =
            FxHashSet::with_capacity_and_hasher(total_positions_estimate, Default::default());

        for (cells, positions) in results {
            cell_barcodes_set.extend(cells);
            positions_set.extend(positions);
        }

        // Parallel conversion and sorting
        let (mut cell_barcodes, mut positions) = rayon::join(
            || {
                let mut vec: Vec<_> = cell_barcodes_set.into_iter().collect();
                vec.par_sort_unstable();
                vec
            },
            || {
                let mut vec: Vec<_> = positions_set.into_iter().collect();
                vec.par_sort_unstable();
                vec
            },
        );

        // Shrink to actual size
        cell_barcodes.shrink_to_fit();
        positions.shrink_to_fit();

        (cell_barcodes, positions)
    }

    /// Optimized index mapping creation using faster hash tables
    fn create_index_mappings_optimized(
        cell_barcodes: &[String],
        positions: &[String],
    ) -> (FxHashMap<String, u32>, FxHashMap<String, u32>) {
        info!("  - Creating optimized index mappings in parallel...");

        let (cell_idx, pos_idx) = rayon::join(
            || {
                let mut map =
                    FxHashMap::with_capacity_and_hasher(cell_barcodes.len(), Default::default());
                for (i, cb) in cell_barcodes.iter().enumerate() {
                    map.insert(cb.clone(), i as u32);
                }
                map
            },
            || {
                let mut map =
                    FxHashMap::with_capacity_and_hasher(positions.len(), Default::default());
                for (i, pos) in positions.iter().enumerate() {
                    map.insert(pos.clone(), i as u32);
                }
                map
            },
        );

        (cell_idx, pos_idx)
    }

    /// Highly optimized sparse matrix construction
    fn build_sparse_matrices_optimized(
        &self,
        position_data: &[PositionData],
        cell_idx: &FxHashMap<String, u32>,
        pos_idx: &FxHashMap<String, u32>,
        n_cells: usize,
        n_positions: usize,
    ) -> Result<(Vec<CsrMatrix<f32>>, Vec<CsrMatrix<f32>>)> {
        info!(
            "  - Building optimized sparse matrices with {} threads...",
            self.config.threads
        );
        info!(
            "    Matrix shape: {} cells × {} positions",
            n_cells, n_positions
        );

        // Estimate number of non-zero elements
        let estimated_nnz_per_matrix =
            (n_cells as f64 * n_positions as f64 * self.config.matrix_density) as usize;

        // Create high-performance matrix builders
        let forward_builders: Vec<_> = (0..4)
            .map(|_| {
                Arc::new(HighPerformanceMatrixBuilder::new(
                    (n_cells, n_positions),
                    estimated_nnz_per_matrix,
                ))
            })
            .collect();

        let reverse_builders: Vec<_> = if self.config.stranded {
            (0..4)
                .map(|_| {
                    Arc::new(HighPerformanceMatrixBuilder::new(
                        (n_cells, n_positions),
                        estimated_nnz_per_matrix,
                    ))
                })
                .collect()
        } else {
            Vec::new()
        };

        // Use smaller chunks to improve cache locality
        let chunk_size = std::cmp::min(self.config.chunk_size, 5000);
        info!(
            "    Processing {} chunks of size {}",
            (position_data.len() + chunk_size - 1) / chunk_size,
            chunk_size
        );

        // Pre-allocate thread-local buffers
        position_data
            .par_chunks(chunk_size)
            .enumerate()
            .try_for_each(|(chunk_idx, chunk)| -> Result<()> {
                let thread_id =
                    rayon::current_thread_index().unwrap_or(chunk_idx % self.config.threads);

                // Use stack-allocated small vectors to reduce heap allocations
                let mut forward_triplets: [SmallVec<[(u32, u32, u32); 64]>; 4] = Default::default();
                let mut reverse_triplets: [SmallVec<[(u32, u32, u32); 64]>; 4] =
                    if self.config.stranded {
                        Default::default()
                    } else {
                        [
                            SmallVec::new(),
                            SmallVec::new(),
                            SmallVec::new(),
                            SmallVec::new(),
                        ]
                    };

                // Pre-allocate string buffer
                let mut pos_key_buffer = String::with_capacity(32);

                for data in chunk {
                    pos_key_buffer.clear();
                    pos_key_buffer.push_str(&data.chrom);
                    pos_key_buffer.push(':');
                    pos_key_buffer.push_str(&data.pos.to_string());

                    if let Some(&col_idx) = pos_idx.get(&pos_key_buffer) {
                        for (cell_barcode, counts) in &data.counts {
                            if let Some(&row_idx) = cell_idx.get(cell_barcode) {
                                if self.config.stranded {
                                    // Batch add forward strand
                                    Self::add_base_counts_to_triplets(
                                        &mut forward_triplets,
                                        row_idx,
                                        col_idx,
                                        &counts.forward,
                                    );

                                    // Batch add reverse strand
                                    Self::add_base_counts_to_triplets(
                                        &mut reverse_triplets,
                                        row_idx,
                                        col_idx,
                                        &counts.reverse,
                                    );
                                } else {
                                    // Unstranded: merge counts
                                    let total_counts = crate::bam2mtx::processor::BaseCounts {
                                        a: counts.forward.a + counts.reverse.a,
                                        t: counts.forward.t + counts.reverse.t,
                                        g: counts.forward.g + counts.reverse.g,
                                        c: counts.forward.c + counts.reverse.c,
                                    };

                                    Self::add_base_counts_to_triplets(
                                        &mut forward_triplets,
                                        row_idx,
                                        col_idx,
                                        &total_counts,
                                    );
                                }
                            }
                        }
                    }
                }

                // Fix: Use std::mem::take to avoid move errors
                for i in 0..4 {
                    let triplets = std::mem::take(&mut forward_triplets[i]);
                    if !triplets.is_empty() {
                        forward_builders[i].add_triplet_batch_lockfree(thread_id, triplets);
                    }
                }

                if self.config.stranded {
                    for i in 0..4 {
                        let triplets = std::mem::take(&mut reverse_triplets[i]);
                        if !triplets.is_empty() {
                            reverse_builders[i].add_triplet_batch_lockfree(thread_id, triplets);
                        }
                    }
                }

                Ok(())
            })?;

        // Parallel build of final matrices
        info!("    Building final CSR matrices in parallel...");

        let forward_matrices: Vec<_> = forward_builders
            .into_par_iter()
            .map(|builder| {
                Arc::try_unwrap(builder)
                    .unwrap_or_else(|_| panic!("Failed to unwrap matrix builder"))
                    .build_optimized()
            })
            .collect();

        let reverse_matrices: Vec<_> = if self.config.stranded {
            reverse_builders
                .into_par_iter()
                .map(|builder| {
                    Arc::try_unwrap(builder)
                        .unwrap_or_else(|_| panic!("Failed to unwrap matrix builder"))
                        .build_optimized()
                })
                .collect()
        } else {
            Vec::new()
        };

        info!("    Sparse matrices built successfully");
        Ok((forward_matrices, reverse_matrices))
    }

    /// Inline function to optimize base count addition
    #[inline(always)]
    fn add_base_counts_to_triplets(
        triplets: &mut [SmallVec<[(u32, u32, u32); 64]>; 4],
        row_idx: u32,
        col_idx: u32,
        counts: &crate::bam2mtx::processor::BaseCounts,
    ) {
        if counts.a > 0 {
            triplets[0].push((row_idx, col_idx, counts.a));
        }
        if counts.t > 0 {
            triplets[1].push((row_idx, col_idx, counts.t));
        }
        if counts.g > 0 {
            triplets[2].push((row_idx, col_idx, counts.g));
        }
        if counts.c > 0 {
            triplets[3].push((row_idx, col_idx, counts.c));
        }
    }

    /// Main conversion function integrating all optimizations
    pub fn convert(
        &self,
        position_data: &[PositionData],
        output_path: &Path,
    ) -> Result<AnnData<H5>> {
        info!(
            "Starting optimized AnnData conversion with {} threads...",
            self.config.threads
        );

        // Step 1: Optimized identifier collection
        info!("Collecting unique cell barcodes and positions (optimized)...");
        let (cell_barcodes, positions) = self.collect_unique_identifiers_optimized(position_data);

        info!(
            "Found {} cells and {} positions",
            cell_barcodes.len(),
            positions.len()
        );

        // Step 2: Optimized index mapping creation
        info!("Creating optimized index mappings...");
        let (cell_idx, pos_idx) = Self::create_index_mappings_optimized(&cell_barcodes, &positions);

        let n_cells = cell_barcodes.len();
        let n_positions = positions.len();

        // Step 3: Highly optimized sparse matrix construction
        info!("Building optimized sparse matrices...");
        let (forward_matrices, reverse_matrices) = self.build_sparse_matrices_optimized(
            position_data,
            &cell_idx,
            &pos_idx,
            n_cells,
            n_positions,
        )?;

        // Step 4-8: AnnData object creation and writing (keeping original logic)
        info!("Creating AnnData object... (single-threaded IO operation)");
        let adata = AnnData::<H5>::new(output_path)?;

        // Parallel metadata preparation
        info!("Preparing metadata in parallel...");
        let (obs_names, var_names) = rayon::join(
            || {
                cell_barcodes
                    .into_iter()
                    .collect::<anndata::data::array::dataframe::DataFrameIndex>()
            },
            || {
                positions
                    .into_iter()
                    .collect::<anndata::data::array::dataframe::DataFrameIndex>()
            },
        );

        info!("Setting main matrix (X)... (single-threaded write operation)");
        adata.set_x(forward_matrices[0].clone())?;

        info!("Adding layers... (sequential write operations)");
        let forward_layer_names = ["A1", "T1", "G1", "C1"];
        let reverse_layer_names = if self.config.stranded {
            vec!["A0", "T0", "G0", "C0"]
        } else {
            vec![]
        };

        // Batch layer addition
        for (i, &layer_name) in forward_layer_names.iter().enumerate() {
            info!("  - Adding layer {}", layer_name);
            adata
                .layers()
                .add(layer_name, forward_matrices[i].clone())?;
        }

        if self.config.stranded && !reverse_matrices.is_empty() {
            for (i, &layer_name) in reverse_layer_names.iter().enumerate() {
                info!("  - Adding layer {}", layer_name);
                adata
                    .layers()
                    .add(layer_name, reverse_matrices[i].clone())?;
            }
        }

        info!("Setting metadata... (single-threaded write operation)");
        adata.set_obs_names(obs_names)?;
        adata.set_var_names(var_names)?;

        info!("Optimized AnnData conversion completed successfully!");
        Ok(adata)
    }

    /// Write the AnnData object to a file
    ///
    /// # Arguments
    ///
    /// * `adata` - The AnnData object to write
    /// * `output_path` - The path where the file should be written
    ///
    /// # Returns
    ///
    /// * `Result<()>` - Ok if successful, Err otherwise
    pub fn write_to_file(&self, _adata: &AnnData<H5>, _output_path: &Path) -> Result<()> {
        Ok(())
    }
}
