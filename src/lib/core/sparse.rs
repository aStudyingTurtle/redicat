//! High-performance sparse matrix utilities shared across REDICAT

use crate::core::error::{RedicatError, Result};
use itertools::Itertools;
use nalgebra_sparse::ops::serial::spadd_csr_prealloc;
use nalgebra_sparse::ops::Op;
use nalgebra_sparse::{CooMatrix, CsrMatrix};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use smallvec::SmallVec;

pub struct SparseOps;

impl SparseOps {
    /// Create CSR matrix from COO format using nalgebra_sparse native conversion
    pub fn from_triplets_u32(
        nrows: usize,
        ncols: usize,
        triplets: Vec<(usize, usize, u32)>,
    ) -> Result<CsrMatrix<u32>> {
        if nrows == 0 || ncols == 0 {
            return Ok(CsrMatrix::zeros(nrows, ncols));
        }

        if triplets.is_empty() {
            return Ok(CsrMatrix::zeros(nrows, ncols));
        }

        // Validate indices
        for &(row, col, _) in &triplets {
            if row >= nrows || col >= ncols {
                return Err(RedicatError::InvalidInput(format!(
                    "Index ({}, {}) exceeds matrix dimensions ({}, {})",
                    row, col, nrows, ncols
                )));
            }
        }

        // Use COO format first, then convert to CSR using native nalgebra_sparse
        let (row_indices, col_indices, values): (Vec<_>, Vec<_>, Vec<_>) =
            triplets.into_iter().multiunzip();

        let coo = CooMatrix::try_from_triplets(nrows, ncols, row_indices, col_indices, values)
            .map_err(|e| RedicatError::SparseMatrix(format!("COO creation failed: {:?}", e)))?;

        // Use native conversion from COO to CSR
        let csr = CsrMatrix::from(&coo);
        Ok(csr)
    }

    /// Create CSR matrix from triplets with u8 values
    pub fn from_triplets(
        nrows: usize,
        ncols: usize,
        triplets: Vec<(usize, usize, u8)>,
    ) -> Result<CsrMatrix<u8>> {
        if nrows == 0 || ncols == 0 {
            return Ok(CsrMatrix::zeros(nrows, ncols));
        }

        if triplets.is_empty() {
            return Ok(CsrMatrix::zeros(nrows, ncols));
        }

        let (row_indices, col_indices, values): (Vec<_>, Vec<_>, Vec<_>) =
            triplets.into_iter().multiunzip();

        let coo = CooMatrix::try_from_triplets(nrows, ncols, row_indices, col_indices, values)
            .map_err(|e| RedicatError::SparseMatrix(format!("COO creation failed: {:?}", e)))?;

        Ok(CsrMatrix::from(&coo))
    }

    /// Highly optimized matrix addition using nalgebra_sparse's spadd operation
    pub fn add_matrices(a: &CsrMatrix<u32>, b: &CsrMatrix<u32>) -> Result<CsrMatrix<u32>> {
        if a.nrows() != b.nrows() || a.ncols() != b.ncols() {
            return Err(RedicatError::DimensionMismatch {
                expected: format!("{}×{}", a.nrows(), a.ncols()),
                actual: format!("{}×{}", b.nrows(), b.ncols()),
            });
        }

        // Use nalgebra_sparse's native sparse addition pattern computation
        let pattern = nalgebra_sparse::ops::serial::spadd_pattern(a.pattern(), b.pattern());

        // Pre-allocate result matrix with computed pattern
        let mut result =
            CsrMatrix::try_from_pattern_and_values(pattern.clone(), vec![0u32; pattern.nnz()])
                .map_err(|e| {
                    RedicatError::SparseMatrix(format!("Failed to create result matrix: {:?}", e))
                })?;

        // Use native sparse addition operation with correct API
        // API signature: spadd_csr_prealloc(beta, C, alpha, Op<A>)
        spadd_csr_prealloc(1u32, &mut result, 1u32, Op::NoOp(a))
            .map_err(|e| RedicatError::SparseMatrix(format!("Sparse addition failed: {:?}", e)))?;

        spadd_csr_prealloc(1u32, &mut result, 1u32, Op::NoOp(b))
            .map_err(|e| RedicatError::SparseMatrix(format!("Sparse addition failed: {:?}", e)))?;

        Ok(result)
    }

    /// Optimized column filtering using sparsity pattern operations
    pub fn filter_columns_u32(
        matrix: &CsrMatrix<u32>,
        keep_indices: &[usize],
    ) -> Result<CsrMatrix<u32>> {
        let nrows = matrix.nrows();
        let new_ncols = keep_indices.len();

        if new_ncols == 0 {
            return Ok(CsrMatrix::zeros(nrows, 0));
        }

        // Create efficient column mapping using FxHashMap for better performance
        let col_map: FxHashMap<usize, usize> = keep_indices
            .iter()
            .enumerate()
            .map(|(new_idx, &old_idx)| (old_idx, new_idx))
            .collect();

        // Build new sparsity pattern efficiently
        let mut new_row_offsets = Vec::with_capacity(nrows + 1);
        let mut new_col_indices = Vec::new();
        let mut new_values = Vec::new();

        new_row_offsets.push(0);

        for row_idx in 0..nrows {
            let row = matrix.row(row_idx);

            for (&old_col, &val) in row.col_indices().iter().zip(row.values()) {
                if let Some(&new_col) = col_map.get(&old_col) {
                    new_col_indices.push(new_col);
                    new_values.push(val);
                }
            }

            new_row_offsets.push(new_col_indices.len());
        }

        // Create CSR matrix directly using nalgebra_sparse
        CsrMatrix::try_from_csr_data(
            nrows,
            new_ncols,
            new_row_offsets,
            new_col_indices,
            new_values,
        )
        .map_err(|e| {
            RedicatError::SparseMatrix(format!("Failed to create filtered matrix: {:?}", e))
        })
    }

    /// Optimized row sums using native CSR structure access
    pub fn compute_row_sums(matrix: &CsrMatrix<u32>) -> Vec<u32> {
        (0..matrix.nrows())
            .into_par_iter()
            .map(|row_idx| {
                let row = matrix.row(row_idx);
                row.values()
                    .iter()
                    .fold(0u64, |acc, &val| acc.saturating_add(val as u64))
                    .min(u32::MAX as u64) as u32
            })
            .collect()
    }

    /// Optimized column sums using parallel reduction over CSR structure
    pub fn compute_col_sums(matrix: &CsrMatrix<u32>) -> Vec<u32> {
        let ncols = matrix.ncols();

        // Use parallel reduction with chunked processing
        let chunk_size = std::cmp::max(1, matrix.nrows() / rayon::current_num_threads());

        (0..matrix.nrows())
            .into_par_iter()
            .chunks(chunk_size)
            .map(|chunk| {
                let mut local_sums = vec![0u64; ncols];
                for row_idx in chunk {
                    let row = matrix.row(row_idx);
                    for (&col_idx, &val) in row.col_indices().iter().zip(row.values()) {
                        local_sums[col_idx] = local_sums[col_idx].saturating_add(val as u64);
                    }
                }
                local_sums
            })
            .reduce(
                || vec![0u64; ncols],
                |mut acc, local| {
                    for (i, val) in local.into_iter().enumerate() {
                        acc[i] = acc[i].saturating_add(val);
                    }
                    acc
                },
            )
            .into_iter()
            .map(|sum| (sum.min(u32::MAX as u64)) as u32)
            .collect()
    }

    /// Element-wise multiplication using optimized two-pointer algorithm for sparse matrices
    /// 
    /// This implementation uses a sorted two-pointer merge algorithm which is significantly
    /// faster than HashMap-based intersection for sparse matrices. CSR format guarantees
    /// column indices are already sorted within each row, allowing O(nnz_a + nnz_b) complexity
    /// instead of O(nnz_a * log(nnz_b)) with HashMap lookups.
    pub fn element_wise_multiply(a: &CsrMatrix<u32>, b: &CsrMatrix<u8>) -> Result<CsrMatrix<u32>> {
        if a.nrows() != b.nrows() || a.ncols() != b.ncols() {
            return Err(RedicatError::DimensionMismatch {
                expected: format!("{}×{}", a.nrows(), a.ncols()),
                actual: format!("{}×{}", b.nrows(), b.ncols()),
            });
        }

        // Use parallel processing with two-pointer algorithm for element-wise multiplication
        // SmallVec optimizes for typical sparse row sizes (most rows have < 32 non-zero elements)
        let triplets: Vec<(usize, usize, u32)> = (0..a.nrows())
            .into_par_iter()
            .flat_map(|row_idx| {
                let a_row = a.row(row_idx);
                let b_row = b.row(row_idx);

                let a_cols = a_row.col_indices();
                let a_vals = a_row.values();
                let b_cols = b_row.col_indices();
                let b_vals = b_row.values();

                // Two-pointer algorithm for sorted sparse row intersection
                // CSR format guarantees column indices are sorted, so we can use merge-like algorithm
                let mut result: SmallVec<[(usize, usize, u32); 32]> = SmallVec::new();
                let mut a_idx = 0;
                let mut b_idx = 0;

                while a_idx < a_cols.len() && b_idx < b_cols.len() {
                    let a_col = a_cols[a_idx];
                    let b_col = b_cols[b_idx];

                    match a_col.cmp(&b_col) {
                        std::cmp::Ordering::Equal => {
                            // Column indices match - this is an intersection point
                            if b_vals[b_idx] > 0 {
                                result.push((row_idx, a_col, a_vals[a_idx]));
                            }
                            a_idx += 1;
                            b_idx += 1;
                        }
                        std::cmp::Ordering::Less => {
                            // a has a column that b doesn't have yet
                            a_idx += 1;
                        }
                        std::cmp::Ordering::Greater => {
                            // b has a column that a doesn't have yet
                            b_idx += 1;
                        }
                    }
                }

                // Convert SmallVec to Vec for flat_map compatibility
                result.into_vec()
            })
            .collect();

        Self::from_triplets_u32(a.nrows(), a.ncols(), triplets)
    }

    /// Transpose operation using nalgebra_sparse native transpose
    pub fn transpose_u32(matrix: &CsrMatrix<u32>) -> CsrMatrix<u32> {
        matrix.transpose()
    }

    /// Matrix-vector multiplication using native spmv operation
    pub fn matrix_vector_multiply(matrix: &CsrMatrix<u32>, vector: &[u32]) -> Result<Vec<u32>> {
        if matrix.ncols() != vector.len() {
            return Err(RedicatError::DimensionMismatch {
                expected: format!("vector length = {}", matrix.ncols()),
                actual: format!("vector length = {}", vector.len()),
            });
        }

        let mut result = vec![0u64; matrix.nrows()];

        // Use parallel processing for matrix-vector multiplication
        result
            .par_iter_mut()
            .enumerate()
            .for_each(|(row_idx, result_val)| {
                let row = matrix.row(row_idx);
                *result_val = row.col_indices().iter().zip(row.values()).fold(
                    0u64,
                    |acc, (&col_idx, &mat_val)| {
                        acc.saturating_add((mat_val as u64) * (vector[col_idx] as u64))
                    },
                );
            });

        Ok(result
            .into_iter()
            .map(|val| (val.min(u32::MAX as u64)) as u32)
            .collect())
    }

    /// Get matrix density statistics
    pub fn get_density_stats(matrix: &CsrMatrix<u32>) -> (f64, usize, usize) {
        let total_elements = matrix.nrows() * matrix.ncols();
        let nnz = matrix.nnz();
        let density = if total_elements > 0 {
            nnz as f64 / total_elements as f64
        } else {
            0.0
        };
        (density, nnz, total_elements)
    }
}

/// Trait extension for additional sparse matrix operations
pub trait SparseMatrixExt<T> {
    fn apply_threshold(&self, threshold: T) -> CsrMatrix<T>
    where
        T: Copy + PartialOrd + Default + nalgebra::Scalar;
}

impl SparseMatrixExt<u32> for CsrMatrix<u32> {
    /// Apply threshold to sparse matrix values
    fn apply_threshold(&self, threshold: u32) -> CsrMatrix<u32> {
        let triplets: Vec<(usize, usize, u32)> = self
            .triplet_iter()
            .filter_map(|(row, col, &val)| {
                if val >= threshold {
                    Some((row, col, val))
                } else {
                    None
                }
            })
            .collect();

        SparseOps::from_triplets_u32(self.nrows(), self.ncols(), triplets)
            .unwrap_or_else(|_| CsrMatrix::zeros(self.nrows(), self.ncols()))
    }
}
