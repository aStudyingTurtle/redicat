//! Optimized AnnData operations with better memory management and error handling

use crate::core::error::{RedicatError, Result};
use crate::core::sparse::SparseOps;
use anndata::data::array::dataframe::DataFrameIndex;
use anndata::{
    data::*,
    traits::{AnnDataOp, AxisArraysOp},
    AnnData, Backend,
};
use anndata_hdf5::H5;
use log::{debug, info, warn};
use nalgebra_sparse::CsrMatrix;
use polars::prelude::*;
use std::collections::HashMap;
use std::convert::TryFrom;
use std::path::Path;

/// Optimized AnnData container with better memory management
#[derive(Debug, Clone)]
pub struct AnnDataContainer {
    pub obs: DataFrame,
    pub var: DataFrame,
    pub x: Option<CsrMatrix<f64>>,
    pub layers: HashMap<String, CsrMatrix<u32>>,
    pub n_obs: usize,
    pub n_vars: usize,
    pub var_names: Vec<String>,
    pub obs_names: Vec<String>,
}

impl AnnDataContainer {
    /// Create a new empty AnnDataContainer
    pub fn new(n_obs: usize, n_vars: usize) -> Self {
        let obs_names: Vec<String> = (0..n_obs).map(|i| format!("cell_{}", i)).collect();
        let var_names: Vec<String> = (0..n_vars).map(|i| format!("gene_{}", i)).collect();

        let obs = DataFrame::new(vec![
            Series::new("obs_names".into(), obs_names.clone()).into()
        ])
        .unwrap();

        let var = DataFrame::new(vec![
            Series::new("var_names".into(), var_names.clone()).into()
        ])
        .unwrap();

        Self {
            obs,
            var,
            x: None,
            layers: HashMap::new(),
            n_obs,
            n_vars,
            var_names,
            obs_names,
        }
    }

    /// Compute row sums for a layer with error handling
    pub fn compute_layer_row_sums(&self, layer_name: &str) -> Option<Vec<u32>> {
        match self.layers.get(layer_name) {
            Some(matrix) => {
                debug!("Computing row sums for layer: {}", layer_name);
                Some(SparseOps::compute_row_sums(matrix))
            }
            None => {
                warn!(
                    "Layer '{}' not found. Available layers: {:?}",
                    layer_name,
                    self.layers.keys().collect::<Vec<_>>()
                );
                None
            }
        }
    }

    /// Compute column sums for a layer with error handling
    pub fn compute_layer_col_sums(&self, layer_name: &str) -> Option<Vec<u32>> {
        match self.layers.get(layer_name) {
            Some(matrix) => {
                debug!("Computing column sums for layer: {}", layer_name);
                Some(SparseOps::compute_col_sums(matrix))
            }
            None => {
                warn!(
                    "Layer '{}' not found. Available layers: {:?}",
                    layer_name,
                    self.layers.keys().collect::<Vec<_>>()
                );
                None
            }
        }
    }

    /// Get total coverage across all base layers with optimized computation
    pub fn compute_total_coverage(&self) -> Vec<u32> {
        // Determine layer strategy based on available layers
        let layer_names = if self.layers.contains_key("A0") {
            // Stranded data
            vec!["A0", "T0", "G0", "C0", "A1", "T1", "G1", "C1"]
        } else {
            // Unstranded data
            vec!["A1", "T1", "G1", "C1"]
        };

        debug!("Computing total coverage using layers: {:?}", layer_names);

        // Collect existing matrices
        let matrices: Vec<&CsrMatrix<u32>> = layer_names
            .iter()
            .filter_map(|&name| self.layers.get(name))
            .collect();

        if matrices.is_empty() {
            warn!("No base matrices found for coverage calculation");
            return vec![0; self.n_obs];
        }

        // Sum matrices efficiently using nalgebra_sparse operations
        let total_matrix =
            matrices
                .into_iter()
                .fold(None, |acc: Option<CsrMatrix<u32>>, matrix| match acc {
                    None => Some(matrix.clone()),
                    Some(existing) => SparseOps::add_matrices(&existing, matrix)
                        .map_err(|e| warn!("Failed to add matrices: {}", e))
                        .ok()
                        .or(Some(existing)),
                });

        match total_matrix {
            Some(matrix) => SparseOps::compute_row_sums(&matrix),
            None => {
                warn!("Failed to compute total coverage matrix");
                vec![0; self.n_obs]
            }
        }
    }

    /// Validate matrix dimensions with flexible obs/var handling
    pub fn validate_dimensions(&self) -> Result<()> {
        // Check that obs_names length matches n_obs
        if self.obs_names.len() != self.n_obs {
            return Err(RedicatError::DimensionMismatch {
                expected: format!("obs_names length = {}", self.n_obs),
                actual: format!("obs_names length = {}", self.obs_names.len()),
            });
        }

        // Check that var_names length matches n_vars
        if self.var_names.len() != self.n_vars {
            return Err(RedicatError::DimensionMismatch {
                expected: format!("var_names length = {}", self.n_vars),
                actual: format!("var_names length = {}", self.var_names.len()),
            });
        }

        // For obs and var DataFrames, allow them to be empty or have different heights
        // as long as they can be reconstructed from the names
        if !self.obs.is_empty() && self.obs.height() != self.n_obs {
            warn!(
                "obs DataFrame height ({}) doesn't match n_obs ({}), will use obs_names",
                self.obs.height(),
                self.n_obs
            );
        }

        if !self.var.is_empty() && self.var.height() != self.n_vars {
            warn!(
                "var DataFrame height ({}) doesn't match n_vars ({}), will use var_names",
                self.var.height(),
                self.n_vars
            );
        }

        // Check X matrix dimensions if present
        if let Some(ref x_matrix) = self.x {
            if x_matrix.nrows() != self.n_obs || x_matrix.ncols() != self.n_vars {
                return Err(RedicatError::DimensionMismatch {
                    expected: format!("X matrix {}×{}", self.n_obs, self.n_vars),
                    actual: format!("X matrix {}×{}", x_matrix.nrows(), x_matrix.ncols()),
                });
            }
        }

        // Check layer dimensions
        for (layer_name, matrix) in &self.layers {
            if matrix.nrows() != self.n_obs || matrix.ncols() != self.n_vars {
                return Err(RedicatError::DimensionMismatch {
                    expected: format!("Layer '{}' {}×{}", layer_name, self.n_obs, self.n_vars),
                    actual: format!(
                        "Layer '{}' {}×{}",
                        layer_name,
                        matrix.nrows(),
                        matrix.ncols()
                    ),
                });
            }
        }

        Ok(())
    }

    /// Ensure obs and var DataFrames have correct dimensions
    pub fn fix_dataframe_dimensions(&mut self) -> Result<()> {
        // Fix obs DataFrame if needed
        if self.obs.is_empty() || self.obs.height() != self.n_obs {
            info!(
                "Reconstructing obs DataFrame with {} observations",
                self.n_obs
            );
            self.obs = DataFrame::new(vec![Series::new(
                "obs_names".into(),
                self.obs_names.clone(),
            )
            .into()])?;
        }

        // Fix var DataFrame if needed
        if self.var.is_empty() || self.var.height() != self.n_vars {
            info!(
                "Reconstructing var DataFrame with {} variables",
                self.n_vars
            );
            self.var = DataFrame::new(vec![Series::new(
                "var_names".into(),
                self.var_names.clone(),
            )
            .into()])?;
        }

        Ok(())
    }

    /// Get memory usage statistics
    pub fn get_memory_stats(&self) -> HashMap<String, usize> {
        let mut stats = HashMap::new();

        // Estimate obs DataFrame size
        stats.insert("obs_bytes".to_string(), estimate_dataframe_size(&self.obs));
        stats.insert("var_bytes".to_string(), estimate_dataframe_size(&self.var));

        // X matrix size
        if let Some(ref x_matrix) = self.x {
            stats.insert(
                "x_bytes".to_string(),
                estimate_csr_matrix_size_f64(x_matrix),
            );
        }

        // Layer sizes
        let mut total_layer_bytes = 0;
        for (layer_name, matrix) in &self.layers {
            let size = estimate_csr_matrix_size_u32(matrix);
            stats.insert(format!("layer_{}_bytes", layer_name), size);
            total_layer_bytes += size;
        }
        stats.insert("total_layer_bytes".to_string(), total_layer_bytes);

        stats
    }

    /// Optimize memory usage by removing empty layers
    pub fn optimize_memory(&mut self) {
        let mut layers_to_remove = Vec::new();

        for (layer_name, matrix) in &self.layers {
            if matrix.nnz() == 0 {
                layers_to_remove.push(layer_name.clone());
            }
        }

        for layer_name in layers_to_remove {
            info!("Removing empty layer: {}", layer_name);
            self.layers.remove(&layer_name);
        }
    }
}

/// Optimized AnnData writing with compression and validation
pub fn write_anndata_h5ad(adata: &AnnDataContainer, path: &str) -> Result<()> {
    info!("Writing AnnData to: {}", path);

    // Create a mutable copy for fixing dimensions
    let mut adata_copy = adata.clone();
    adata_copy.fix_dataframe_dimensions()?;
    adata_copy.validate_dimensions()?;

    // Log memory usage
    let stats = adata_copy.get_memory_stats();
    info!(
        "Memory usage: obs={} KB, var={} KB, layers={} KB",
        stats.get("obs_bytes").unwrap_or(&0) / 1024,
        stats.get("var_bytes").unwrap_or(&0) / 1024,
        stats.get("total_layer_bytes").unwrap_or(&0) / 1024
    );

    let h5_adata = AnnData::<H5>::new(Path::new(path))?;

    // Set observation and variable indices
    let obs_index: DataFrameIndex = adata_copy.obs_names.iter().cloned().collect();
    let var_index: DataFrameIndex = adata_copy.var_names.iter().cloned().collect();
    h5_adata.set_obs_names(obs_index)?;
    h5_adata.set_var_names(var_index)?;

    // Set main matrix X with compression
    if let Some(ref x_matrix) = adata_copy.x {
        let x_f32 = convert_f64_to_f32_csr(x_matrix)?;
        h5_adata.set_x(x_f32)?;
        info!(
            "  - Written X matrix: {}×{} with {} non-zeros",
            x_matrix.nrows(),
            x_matrix.ncols(),
            x_matrix.nnz()
        );
    } else {
        let zero_matrix = CsrMatrix::<f32>::zeros(adata_copy.n_obs, adata_copy.n_vars);
        h5_adata.set_x(zero_matrix)?;
        info!(
            "  - Written empty X matrix: {}×{}",
            adata_copy.n_obs, adata_copy.n_vars
        );
    }

    // Set layers with priority order for important layers
    let priority_layers = ["ref", "alt", "others", "coverage"];
    let mut written_layers = 0;

    // Write priority layers first
    for layer_name in &priority_layers {
        if let Some(layer_matrix) = adata_copy.layers.get(*layer_name) {
            info!(
                "  - Writing layer: {} ({}×{}, {} non-zeros)",
                layer_name,
                layer_matrix.nrows(),
                layer_matrix.ncols(),
                layer_matrix.nnz()
            );
            let f32_matrix = convert_u32_to_f32_csr(layer_matrix)?;
            h5_adata.layers().add(layer_name, f32_matrix)?;
            written_layers += 1;
        }
    }

    // // Write remaining layers
    // for (layer_name, layer_matrix) in &adata_copy.layers {
    //     if !priority_layers.contains(&layer_name.as_str()) {
    //         info!("  - Writing layer: {} ({}×{}, {} non-zeros)",
    //               layer_name, layer_matrix.nrows(), layer_matrix.ncols(), layer_matrix.nnz());
    //         let f32_matrix = convert_u32_to_f32_csr(layer_matrix)?;
    //         h5_adata.layers().add(layer_name, f32_matrix)?;
    //         written_layers += 1;
    //     }
    // }

    // Set annotations
    if !adata_copy.obs.is_empty() {
        h5_adata.set_obs(adata_copy.obs.clone())?;
        info!(
            "  - Written obs annotations: {} rows, {} columns",
            adata_copy.obs.height(),
            adata_copy.obs.width()
        );
    }

    if !adata_copy.var.is_empty() {
        h5_adata.set_var(adata_copy.var.clone())?;
        info!(
            "  - Written var annotations: {} rows, {} columns",
            adata_copy.var.height(),
            adata_copy.var.width()
        );
    }

    h5_adata.set_n_obs(adata_copy.n_obs)?;
    h5_adata.set_n_vars(adata_copy.n_vars)?;

    info!(
        "Successfully wrote AnnData with shape: {} × {}, {} layers",
        adata_copy.n_obs, adata_copy.n_vars, written_layers
    );
    Ok(())
}

/// Optimized AnnData reading with better error handling and dimension fixing
pub fn read_anndata_h5ad(path: &str) -> Result<AnnDataContainer> {
    info!("Reading H5AD file: {}", path);

    if !std::path::Path::new(path).exists() {
        return Err(RedicatError::FileNotFound(format!(
            "File not found: {}",
            path
        )));
    }

    let adata =
        AnnData::<H5>::open(H5::open(path).map_err(|e| {
            RedicatError::DataProcessing(format!("Failed to open H5 file: {:?}", e))
        })?)?;

    let n_obs = adata.n_obs();
    let n_vars = adata.n_vars();
    info!("AnnData shape: {} obs × {} vars", n_obs, n_vars);

    if n_obs == 0 || n_vars == 0 {
        return Err(RedicatError::EmptyData(format!(
            "Empty AnnData: {} obs × {} vars",
            n_obs, n_vars
        )));
    }

    let obs_names = read_names(&adata.obs_names())?;
    let var_names = read_names(&adata.var_names())?;

    // Validate that names match dimensions
    if obs_names.len() != n_obs {
        return Err(RedicatError::DimensionMismatch {
            expected: format!("obs_names length = {}", n_obs),
            actual: format!("obs_names length = {}", obs_names.len()),
        });
    }

    if var_names.len() != n_vars {
        return Err(RedicatError::DimensionMismatch {
            expected: format!("var_names length = {}", n_vars),
            actual: format!("var_names length = {}", var_names.len()),
        });
    }

    let obs = read_obs_dataframe(&adata, &obs_names, n_obs)?;
    let var = read_var_dataframe(&adata, &var_names, n_vars)?;
    let x = read_x_matrix(&adata)?;
    let layers = read_layers_as_u32(&adata)?;

    info!(
        "Successfully loaded AnnData with {} layers: {:?}",
        layers.len(),
        layers.keys().collect::<Vec<_>>()
    );

    let mut container = AnnDataContainer {
        obs,
        var,
        x,
        layers,
        n_obs,
        n_vars,
        var_names,
        obs_names,
    };

    // Fix any dimension mismatches
    container.fix_dataframe_dimensions()?;

    // Validate the fixed data
    container.validate_dimensions()?;

    Ok(container)
}

// Helper functions with better error handling

fn read_names(index: &DataFrameIndex) -> Result<Vec<String>> {
    Ok(index.clone().into_vec())
}

fn read_obs_dataframe(
    adata: &AnnData<H5>,
    obs_names: &[String],
    n_obs: usize,
) -> Result<DataFrame> {
    match adata.read_obs() {
        Ok(obs_df) => {
            debug!(
                "Read obs DataFrame: {} rows, {} columns",
                obs_df.height(),
                obs_df.width()
            );
            if obs_df.height() == n_obs {
                Ok(obs_df)
            } else {
                warn!(
                    "obs DataFrame height ({}) doesn't match n_obs ({}), creating from names",
                    obs_df.height(),
                    n_obs
                );
                DataFrame::new(vec![
                    Series::new("obs_names".into(), obs_names.to_vec()).into()
                ])
                .map_err(|e| {
                    RedicatError::DataProcessing(format!("Failed to create obs DataFrame: {}", e))
                })
            }
        }
        Err(e) => {
            warn!("Failed to read obs DataFrame: {:?}, creating from names", e);
            DataFrame::new(vec![
                Series::new("obs_names".into(), obs_names.to_vec()).into()
            ])
            .map_err(|e| {
                RedicatError::DataProcessing(format!("Failed to create obs DataFrame: {}", e))
            })
        }
    }
}

fn read_var_dataframe(
    adata: &AnnData<H5>,
    var_names: &[String],
    n_vars: usize,
) -> Result<DataFrame> {
    match adata.read_var() {
        Ok(var_df) => {
            debug!(
                "Read var DataFrame: {} rows, {} columns",
                var_df.height(),
                var_df.width()
            );
            if var_df.height() == n_vars {
                Ok(var_df)
            } else {
                warn!(
                    "var DataFrame height ({}) doesn't match n_vars ({}), creating from names",
                    var_df.height(),
                    n_vars
                );
                DataFrame::new(vec![
                    Series::new("var_names".into(), var_names.to_vec()).into()
                ])
                .map_err(|e| {
                    RedicatError::DataProcessing(format!("Failed to create var DataFrame: {}", e))
                })
            }
        }
        Err(e) => {
            warn!("Failed to read var DataFrame: {:?}, creating from names", e);
            DataFrame::new(vec![
                Series::new("var_names".into(), var_names.to_vec()).into()
            ])
            .map_err(|e| {
                RedicatError::DataProcessing(format!("Failed to create var DataFrame: {}", e))
            })
        }
    }
}

fn read_x_matrix(adata: &AnnData<H5>) -> Result<Option<CsrMatrix<f64>>> {
    let mut x_elem = match adata.x().extract() {
        Some(elem) => elem,
        None => {
            debug!("No X matrix found");
            return Ok(None);
        }
    };

    let shape = x_elem.shape();
    if shape.ndim() == 0 || shape.as_ref().contains(&0) {
        debug!("Empty X matrix shape: {:?}", shape.as_ref());
        return Ok(None);
    }

    match x_elem.data() {
        Ok(array_data) => match convert_array_to_csr_f64(array_data) {
            Ok(matrix) => {
                info!(
                    "Read X matrix: {}×{} with {} non-zeros",
                    matrix.nrows(),
                    matrix.ncols(),
                    matrix.nnz()
                );
                Ok(Some(matrix))
            }
            Err(e) => {
                warn!("Failed to convert X matrix: {}", e);
                Ok(None)
            }
        },
        Err(e) => {
            warn!("Failed to extract X matrix data: {:?}", e);
            Ok(None)
        }
    }
}

// Fixed layer reading function
fn read_layers_as_u32(adata: &AnnData<H5>) -> Result<HashMap<String, CsrMatrix<u32>>> {
    let mut layers: HashMap<String, CsrMatrix<u32>> = HashMap::new();
    let layers_ref = adata.layers();

    // Common layer names to try
    let common_layer_names = vec![
        "A0", "T0", "G0", "C0", "A1", "T1", "G1", "C1", "ref", "alt", "others", "coverage",
    ];

    info!("Attempting to load common layers: {:?}", common_layer_names);

    for layer_name in common_layer_names {
        match layers_ref.get_item::<ArrayData>(layer_name) {
            Ok(Some(array_data)) => match convert_array_to_csr_u32(array_data) {
                Ok(matrix) => {
                    info!(
                        "  - Loaded layer '{}': {}×{} with {} non-zeros",
                        layer_name,
                        matrix.nrows(),
                        matrix.ncols(),
                        matrix.nnz()
                    );
                    layers.insert(layer_name.to_string(), matrix);
                }
                Err(e) => {
                    warn!("  - Failed to convert layer '{}': {}", layer_name, e);
                }
            },
            Ok(None) => {
                debug!("  - Layer '{}' not found (normal)", layer_name);
            }
            Err(_) => {
                debug!(
                    "  - Could not access layer '{}' (normal if it doesn't exist)",
                    layer_name
                );
            }
        }
    }

    info!("Successfully loaded {} layers", layers.len());
    Ok(layers)
}

// Conversion functions with better error handling

fn convert_f64_to_f32_csr(matrix: &CsrMatrix<f64>) -> Result<CsrMatrix<f32>> {
    let (row_offsets, col_indices, values) = matrix.csr_data();
    let values_f32: Vec<f32> = values.iter().map(|&x| x as f32).collect();

    CsrMatrix::try_from_csr_data(
        matrix.nrows(),
        matrix.ncols(),
        row_offsets.to_vec(),
        col_indices.to_vec(),
        values_f32,
    )
    .map_err(|e| RedicatError::DataProcessing(format!("Failed to convert f64 to f32: {:?}", e)))
}

fn convert_f32_to_f64_csr(matrix: &CsrMatrix<f32>) -> Result<CsrMatrix<f64>> {
    let (row_offsets, col_indices, values) = matrix.csr_data();
    let values_f64: Vec<f64> = values.iter().map(|&x| x as f64).collect();

    CsrMatrix::try_from_csr_data(
        matrix.nrows(),
        matrix.ncols(),
        row_offsets.to_vec(),
        col_indices.to_vec(),
        values_f64,
    )
    .map_err(|e| RedicatError::DataProcessing(format!("Failed to convert f32 to f64: {:?}", e)))
}

fn convert_u32_to_f32_csr(matrix: &CsrMatrix<u32>) -> Result<CsrMatrix<f32>> {
    let (row_offsets, col_indices, values) = matrix.csr_data();
    let values_f32: Vec<f32> = values.iter().map(|&x| x as f32).collect();

    CsrMatrix::try_from_csr_data(
        matrix.nrows(),
        matrix.ncols(),
        row_offsets.to_vec(),
        col_indices.to_vec(),
        values_f32,
    )
    .map_err(|e| RedicatError::DataProcessing(format!("Failed to convert u32 to f32: {:?}", e)))
}

fn convert_u32_to_f64_csr(matrix: &CsrMatrix<u32>) -> Result<CsrMatrix<f64>> {
    let (row_offsets, col_indices, values) = matrix.csr_data();
    let values_f64: Vec<f64> = values.iter().map(|&x| x as f64).collect();

    CsrMatrix::try_from_csr_data(
        matrix.nrows(),
        matrix.ncols(),
        row_offsets.to_vec(),
        col_indices.to_vec(),
        values_f64,
    )
    .map_err(|e| RedicatError::DataProcessing(format!("Failed to convert u32 to f64: {:?}", e)))
}

fn convert_array_to_csr_f64(array_data: ArrayData) -> Result<CsrMatrix<f64>> {
    if let Ok(matrix) = CsrMatrix::<f64>::try_from(array_data.clone()) {
        return Ok(matrix);
    }

    if let Ok(matrix_f32) = CsrMatrix::<f32>::try_from(array_data.clone()) {
        return convert_f32_to_f64_csr(&matrix_f32);
    }

    if let Ok(matrix_u32) = CsrMatrix::<u32>::try_from(array_data.clone()) {
        return convert_u32_to_f64_csr(&matrix_u32);
    }

    Err(RedicatError::DataProcessing(format!(
        "Unsupported array data type for X matrix: {:?}",
        array_data.data_type()
    )))
}

fn convert_array_to_csr_u32(array_data: ArrayData) -> Result<CsrMatrix<u32>> {
    if let Ok(matrix) = CsrMatrix::<u32>::try_from(array_data.clone()) {
        return Ok(matrix);
    }

    if let Ok(matrix_f32) = CsrMatrix::<f32>::try_from(array_data.clone()) {
        let (row_offsets, col_indices, values) = matrix_f32.csr_data();
        let values_u32: Vec<u32> = values.iter().map(|&x| x as u32).collect();
        return CsrMatrix::try_from_csr_data(
            matrix_f32.nrows(),
            matrix_f32.ncols(),
            row_offsets.to_vec(),
            col_indices.to_vec(),
            values_u32,
        )
        .map_err(|e| {
            RedicatError::DataProcessing(format!("Failed to convert f32 to u32: {:?}", e))
        });
    }

    if let Ok(matrix_f64) = CsrMatrix::<f64>::try_from(array_data.clone()) {
        let (row_offsets, col_indices, values) = matrix_f64.csr_data();
        let values_u32: Vec<u32> = values.iter().map(|&x| x as u32).collect();
        return CsrMatrix::try_from_csr_data(
            matrix_f64.nrows(),
            matrix_f64.ncols(),
            row_offsets.to_vec(),
            col_indices.to_vec(),
            values_u32,
        )
        .map_err(|e| {
            RedicatError::DataProcessing(format!("Failed to convert f64 to u32: {:?}", e))
        });
    }

    Err(RedicatError::DataProcessing(format!(
        "Unsupported array data type for layer: {:?}",
        array_data.data_type()
    )))
}

fn estimate_dataframe_size(df: &DataFrame) -> usize {
    df.get_columns()
        .iter()
        .map(|column| column.as_materialized_series().estimated_size())
        .sum()
}

fn estimate_csr_matrix_size_f64(matrix: &CsrMatrix<f64>) -> usize {
    let (row_offsets, col_indices, values) = matrix.csr_data();
    std::mem::size_of_val(row_offsets)
        + std::mem::size_of_val(col_indices)
        + std::mem::size_of_val(values)
}

fn estimate_csr_matrix_size_u32(matrix: &CsrMatrix<u32>) -> usize {
    let (row_offsets, col_indices, values) = matrix.csr_data();
    std::mem::size_of_val(row_offsets)
        + std::mem::size_of_val(col_indices)
        + std::mem::size_of_val(values)
}

/// Estimate total memory footprint in bytes for the AnnDataContainer
pub fn estimate_anndata_memory_usage(adata: &AnnDataContainer) -> usize {
    let mut total = 0;
    total += estimate_dataframe_size(&adata.obs);
    total += estimate_dataframe_size(&adata.var);
    if let Some(ref x) = adata.x {
        total += estimate_csr_matrix_size_f64(x);
    }
    for matrix in adata.layers.values() {
        total += estimate_csr_matrix_size_u32(matrix);
    }
    total
}
