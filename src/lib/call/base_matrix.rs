//! Optimized base matrix operations with efficient memory usage

use crate::call::anndata_ops::AnnDataContainer;
use crate::call::error::Result;
use crate::call::sparse_ops::SparseOps;
use log::info;
use nalgebra_sparse::CsrMatrix;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;

/// Calculate coverage matrix using optimized sparse matrix addition
pub fn calculate_coverage(adata: &mut AnnDataContainer) -> Result<()> {
    info!("Calculating coverage matrix...");

    // Determine which layers to use based on availability
    let layer_names = if adata.layers.contains_key("A0") {
        // Stranded data - use both strands
        vec!["A0", "T0", "G0", "C0", "A1", "T1", "G1", "C1"]
    } else {
        // Unstranded data - use only positive strand
        vec!["A1", "T1", "G1", "C1"]
    };

    // Collect all existing matrices first to avoid borrowing issues
    let matrices: Vec<CsrMatrix<u32>> = layer_names
        .iter()
        .filter_map(|&name| adata.layers.get(name))
        .cloned()
        .collect();

    if matrices.is_empty() {
        let coverage = CsrMatrix::<u32>::zeros(adata.n_obs, adata.n_vars);
        adata.layers.insert("coverage".to_string(), coverage);
        return Ok(());
    }

    // Efficiently sum all matrices using reduce operation
    let coverage = matrices
        .into_iter()
        .reduce(|acc, matrix| {
            SparseOps::add_matrices(&acc, &matrix).unwrap_or(acc)
        })
        .unwrap_or_else(|| CsrMatrix::<u32>::zeros(adata.n_obs, adata.n_vars));

    adata.layers.insert("coverage".to_string(), coverage);
    info!("Coverage matrix calculated with {} non-zero elements", 
          adata.layers.get("coverage").map(|m| m.nnz()).unwrap_or(0));
    Ok(())
}

/// Optimized site filtering with parallel processing
pub fn filter_sites_by_coverage(
    mut adata: AnnDataContainer,
    min_coverage: u32,
) -> Result<AnnDataContainer> {
    info!("Filtering sites with min_coverage: {}", min_coverage);

    // Ensure coverage is calculated
    calculate_coverage(&mut adata)?;

    // Compute site coverage in parallel
    let site_coverage = adata.compute_layer_col_sums("coverage")
        .unwrap_or_else(|| vec![0; adata.n_vars]);

    // Create filter mask using parallel iterator
    let filter_mask: Vec<bool> = adata.var_names
        .par_iter()
        .zip(site_coverage.par_iter())
        .map(|(name, &cov)| {
            name.starts_with("chr") && cov >= min_coverage
        })
        .collect();

    let kept_sites = filter_mask.iter().filter(|&&x| x).count();
    info!("Keeping {} out of {} sites ({}% retained)", 
          kept_sites, adata.n_vars, 
          (kept_sites as f64 / adata.n_vars as f64 * 100.0) as u32);

    if kept_sites == 0 {
        return Err(crate::call::error::RedicatError::EmptyData(
            "No sites passed coverage filter".to_string(),
        ));
    }

    apply_site_filter(adata, &filter_mask)
}

/// Optimized site filtering with better memory usage
pub fn apply_site_filter(
    mut adata: AnnDataContainer,
    filter_mask: &[bool],
) -> Result<AnnDataContainer> {
    // Pre-compute selected indices once
    let selected_indices: Vec<usize> = filter_mask
        .par_iter()
        .enumerate()
        .filter_map(|(i, &keep)| if keep { Some(i) } else { None })
        .collect();

    if selected_indices.is_empty() {
        return Err(crate::call::error::RedicatError::EmptyData(
            "No sites selected after filtering".to_string(),
        ));
    }

    // Process layers in parallel to improve performance
    let layer_names: Vec<String> = adata.layers.keys().cloned().collect();
    let filtered_layers: Result<HashMap<String, CsrMatrix<u32>>> = layer_names
        .into_par_iter()
        .map(|name| {
            let matrix = adata.layers.get(&name).unwrap();
            SparseOps::filter_columns_u32(matrix, &selected_indices)
                .map(|filtered| (name, filtered))
        })
        .collect();

    // Filter var DataFrame
    let filtered_var = filter_dataframe_by_indices(&adata.var, &selected_indices)?;
    
    // Filter var_names in parallel
    let filtered_var_names: Vec<String> = selected_indices
        .par_iter()
        .map(|&i| adata.var_names[i].clone())
        .collect();

    // Update adata with filtered data
    adata.layers = filtered_layers?;
    adata.var = filtered_var;
    adata.n_vars = selected_indices.len();
    adata.var_names = filtered_var_names;

    info!("Site filtering completed: {} sites retained", adata.n_vars);
    Ok(adata)
}

/// Optimized base level counting with parallel processing
pub fn count_base_levels(mut adata: AnnDataContainer) -> Result<AnnDataContainer> {
    info!("Counting base levels at each site...");

    // Process all bases in parallel
    let base_counts: Vec<(char, Vec<u32>)> = ['A', 'T', 'G', 'C']
        .par_iter()
        .map(|&base| {
            let layer_name = format!("{}1", base);
            let counts = adata.compute_layer_col_sums(&layer_name)
                .unwrap_or_else(|| vec![0u32; adata.n_vars]);
            (base, counts)
        })
        .collect();

    // Add all base count columns - fix the mutation issue
    for (base, counts) in base_counts {
        let series = Series::new(base.to_string().into(), counts);
        adata.var = adata.var.with_column(series).unwrap().clone(); // Fixed: assign back to adata.var
    }

    // Calculate total coverage using vectorized operations
    let coverage_expr = col("A") + col("T") + col("G") + col("C");
    adata.var = adata.var
        .lazy()
        .with_columns([coverage_expr.alias("Coverage")])
        .collect()?;

    info!("Base levels counted successfully");
    Ok(adata)
}

/// Optimized DataFrame filtering using vectorized operations
fn filter_dataframe_by_indices(df: &DataFrame, indices: &[usize]) -> Result<DataFrame> {
    if indices.is_empty() {
        return Err(crate::call::error::RedicatError::EmptyData(
            "No indices provided for filtering".to_string(),
        ));
    }

    // Create boolean mask more efficiently
    let mask = (0..df.height())
        .into_par_iter()
        .map(|i| indices.binary_search(&i).is_ok())
        .collect::<Vec<bool>>();

    let mask_chunked = BooleanChunked::from_slice("mask".into(), &mask);

    df.filter(&mask_chunked).map_err(|e| {
        crate::call::error::RedicatError::DataProcessing(format!(
            "Failed to filter DataFrame: {}", e
        ))
    })
}