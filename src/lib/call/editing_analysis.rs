//! RNA editing analysis optimized with nalgebra_sparse native operations
//!
//! This module provides the core analysis functions for RNA editing detection
//! and quantification in single-cell data. The analysis is performed in a
//! strand-aware manner to properly handle editing events on both positive
//! and negative DNA strands.

use crate::call::anndata_ops::AnnDataContainer;
use crate::call::base_matrix::*;
use crate::call::error::Result;
use crate::call::sparse_ops::SparseOps;
use crate::call::{EditingType, ReferenceGenome};
use itertools::Itertools;
use log::info;
use nalgebra_sparse::CsrMatrix;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

/// Optimized ref/alt matrix calculation using vectorized sparse operations
///
/// This function calculates reference and alternate allele count matrices for
/// the specified RNA editing type in a strand-aware manner. The strand-aware
/// processing ensures that editing events are correctly identified and
/// quantified regardless of which DNA strand they originate from.
///
/// The function performs the following steps:
/// 1. Filters genomic sites to only those relevant for the specified editing type
/// 2. Extracts reference base information for all sites
/// 3. Computes reference, alternate, and other base count matrices using
///    strand-aware logic that considers both positive and negative strand
///    editing events
/// 4. Sets the computed matrices as layers in the AnnData container
/// 5. Calculates observation-level statistics
///
/// # Arguments
///
/// * `adata` - The AnnDataContainer with base count matrices
/// * `editing_type` - The editing type to analyze (e.g., AG, CT, etc.)
///
/// # Returns
///
/// Updated AnnDataContainer with ref/alt matrices and observation statistics
pub fn calculate_ref_alt_matrices(
    mut adata: AnnDataContainer,
    editing_type: &EditingType,
) -> Result<AnnDataContainer> {
    info!(
        "Calculating strand-aware ref/alt matrices for editing type: {:?}",
        editing_type
    );

    // Filter sites by editing type using strand-aware logic
    adata = filter_by_editing_type_strand_aware(adata, editing_type)?;

    if adata.n_vars == 0 {
        return Err(crate::call::error::RedicatError::EmptyData(
            "No sites remain after strand-aware editing type filtering".to_string(),
        ));
    }

    // Extract reference sequence information
    let ref_bases = extract_reference_bases(&adata.var)?;

    // Calculate editing matrices using highly optimized vectorized operations
    // with strand-aware base assignment
    let (ref_matrix, alt_matrix, others_matrix) =
        compute_editing_matrices_vectorized(&adata, &ref_bases, editing_type)?;

    // Set matrices
    adata.layers.insert("ref".to_string(), ref_matrix);
    adata.layers.insert("alt".to_string(), alt_matrix.clone());
    adata.layers.insert("others".to_string(), others_matrix);

    // Set main matrix X as f64 version of alt matrix
    adata.x = Some(convert_u32_to_f64_csr(&alt_matrix));

    // Calculate observation-level statistics using vectorized operations
    adata = calculate_observation_sums_vectorized(adata)?;

    info!("Strand-aware ref/alt matrices calculated using vectorized operations");
    Ok(adata)
}

/// Highly optimized computation using nalgebra_sparse vectorized operations
fn compute_editing_matrices_vectorized(
    adata: &AnnDataContainer,
    ref_bases: &[char],
    editing_type: &EditingType,
) -> Result<(CsrMatrix<u32>, CsrMatrix<u32>, CsrMatrix<u32>)> {
    info!("Computing editing matrices with vectorized sparse operations");

    // Pre-compute base mappings for vectorized operations
    let (ref_site_masks, alt_site_masks) = create_vectorized_base_masks(ref_bases, editing_type);

    // Process all base layers in parallel using native sparse operations
    let bases = ['A', 'T', 'G', 'C'];
    let base_matrices: Vec<(char, &CsrMatrix<u32>)> = bases
        .iter()
        .filter_map(|&base| {
            let layer_name = format!("{}1", base);
            adata.layers.get(&layer_name).map(|matrix| (base, matrix))
        })
        .collect();

    if base_matrices.is_empty() {
        return Err(crate::call::error::RedicatError::DataProcessing(
            "No base matrices found for editing calculation".to_string(),
        ));
    }

    // Use parallel reduction to compute final matrices
    let (ref_matrix, alt_matrix, others_matrix) = base_matrices
        .into_par_iter()
        .map(|(base, matrix)| {
            compute_base_contributions_vectorized(matrix, &ref_site_masks, &alt_site_masks, base)
        })
        .reduce(
            || {
                Ok((
                    CsrMatrix::zeros(adata.n_obs, adata.n_vars),
                    CsrMatrix::zeros(adata.n_obs, adata.n_vars),
                    CsrMatrix::zeros(adata.n_obs, adata.n_vars),
                ))
            },
            |acc, curr| match (acc, curr) {
                (Ok((acc_ref, acc_alt, acc_others)), Ok((curr_ref, curr_alt, curr_others))) => {
                    Ok((
                        SparseOps::add_matrices(&acc_ref, &curr_ref)?,
                        SparseOps::add_matrices(&acc_alt, &curr_alt)?,
                        SparseOps::add_matrices(&acc_others, &curr_others)?,
                    ))
                }
                (Err(e), _) | (_, Err(e)) => Err(e),
            },
        )?;

    info!("Completed vectorized editing matrix computation");
    Ok((ref_matrix, alt_matrix, others_matrix))
}

/// Create vectorized base masks for efficient sparse operations
fn create_vectorized_base_masks(
    ref_bases: &[char],
    editing_type: &EditingType,
) -> (HashMap<char, Vec<bool>>, HashMap<char, Vec<bool>>) {
    let bases = ['A', 'T', 'G', 'C'];
    let n_sites = ref_bases.len();

    let mut ref_site_masks = HashMap::new();
    let mut alt_site_masks = HashMap::new();

    // Initialize masks for all bases
    for &base in &bases {
        ref_site_masks.insert(base, vec![false; n_sites]);
        alt_site_masks.insert(base, vec![false; n_sites]);
    }

    // Vectorized mask computation using immutable operations
    let ref_updates: Vec<(usize, char)> = ref_bases
        .par_iter()
        .enumerate()
        .map(|(site_idx, &ref_base)| (site_idx, ref_base))
        .collect();

    let alt_updates: Vec<(usize, char)> = ref_bases
        .par_iter()
        .enumerate()
        .map(|(site_idx, &ref_base)| (site_idx, editing_type.get_alt_base_for_ref(ref_base)))
        .collect();

    // Apply updates to masks
    for (site_idx, ref_base) in ref_updates {
        if let Some(ref_mask) = ref_site_masks.get_mut(&ref_base) {
            ref_mask[site_idx] = true;
        }
    }

    for (site_idx, alt_base) in alt_updates {
        if alt_base != 'N' {
            if let Some(alt_mask) = alt_site_masks.get_mut(&alt_base) {
                alt_mask[site_idx] = true;
            }
        }
    }

    (ref_site_masks, alt_site_masks)
}

/// Compute base contributions using highly optimized sparse operations
fn compute_base_contributions_vectorized(
    base_matrix: &CsrMatrix<u32>,
    ref_site_masks: &HashMap<char, Vec<bool>>,
    alt_site_masks: &HashMap<char, Vec<bool>>,
    base: char,
) -> Result<(CsrMatrix<u32>, CsrMatrix<u32>, CsrMatrix<u32>)> {
    // Get masks for this base
    let ref_mask = ref_site_masks.get(&base).unwrap();
    let alt_mask = alt_site_masks.get(&base).unwrap();

    // Partition triplets without intermediate allocation
    let nnz = base_matrix.nnz();
    let mask_len = ref_mask.len();

    let mut ref_triplets = Vec::with_capacity(nnz / 3 + 1);
    let mut alt_triplets = Vec::with_capacity(nnz / 3 + 1);
    let mut others_triplets = Vec::with_capacity(nnz / 3 + 1);

    for (row_idx, col_idx, &val) in base_matrix.triplet_iter() {
        if val == 0 || col_idx >= mask_len {
            continue;
        }

        if ref_mask[col_idx] {
            ref_triplets.push((row_idx, col_idx, val));
        } else if alt_mask[col_idx] {
            alt_triplets.push((row_idx, col_idx, val));
        } else {
            others_triplets.push((row_idx, col_idx, val));
        }
    }

    // Create matrices using native nalgebra_sparse operations
    let ref_matrix =
        SparseOps::from_triplets_u32(base_matrix.nrows(), base_matrix.ncols(), ref_triplets)?;

    let alt_matrix =
        SparseOps::from_triplets_u32(base_matrix.nrows(), base_matrix.ncols(), alt_triplets)?;

    let others_matrix =
        SparseOps::from_triplets_u32(base_matrix.nrows(), base_matrix.ncols(), others_triplets)?;

    Ok((ref_matrix, alt_matrix, others_matrix))
}

/// Vectorized observation-level sum calculation
fn calculate_observation_sums_vectorized(mut adata: AnnDataContainer) -> Result<AnnDataContainer> {
    info!("Calculating observation-level sums using vectorized operations");

    // Use parallel processing for all layer sum calculations
    let layer_sums: Vec<(String, Vec<u32>)> = ["ref", "alt", "others"]
        .par_iter()
        .filter_map(|&layer_name| {
            adata
                .compute_layer_row_sums(layer_name)
                .map(|sums| (layer_name.to_string(), sums))
        })
        .collect();

    if !layer_sums.is_empty() {
        let columns: Vec<Column> = layer_sums
            .into_iter()
            .map(|(layer_name, sums)| Series::new(layer_name.into(), sums).into_column())
            .collect();
        adata.obs.hstack_mut(&columns)?;
    }

    info!("Vectorized observation-level sums calculated");
    Ok(adata)
}

// Keep the other helper functions with DataFrame mutation fixes

pub fn annotate_variants_pipeline(
    adata: AnnDataContainer,
    editing_sites: Arc<HashMap<String, u8>>,
    reference: Arc<ReferenceGenome>,
    max_other_threshold: f32,
    min_edited_threshold: f32,
    min_ref_threshold: f32,
    min_coverage: u32,
) -> Result<AnnDataContainer> {
    info!("Starting variant annotation pipeline...");

    let mut adata = adata;
    adata = mark_editing_sites(adata, &editing_sites)?;
    adata = filter_sites_by_coverage(adata, min_coverage)?;
    adata = count_base_levels(adata)?;
    adata = add_reference_bases(adata, reference)?;
    adata = apply_mismatch_filtering(
        adata,
        max_other_threshold,
        min_edited_threshold,
        min_ref_threshold,
    )?;

    info!("Variant annotation completed");
    Ok(adata)
}

pub fn calculate_cei(mut adata: AnnDataContainer) -> Result<AnnDataContainer> {
    info!("Calculating Cell Editing Index (CEI)...");

    let cei_expr = col("alt").cast(DataType::Float32)
        / (col("ref").cast(DataType::Float32) + col("alt").cast(DataType::Float32));

    let cei_series = adata
        .obs
        .clone()
        .lazy()
        .with_columns([cei_expr.fill_null(0.0).alias("CEI")])
        .collect()?
        .column("CEI")?
        .clone();

    // Fixed DataFrame mutation
    adata.obs.hstack_mut(&[cei_series.into_column()])?;
    info!("CEI calculated");
    Ok(adata)
}

pub fn calculate_site_mismatch_stats(
    mut adata: AnnDataContainer,
    ref_base: char,
    alt_base: char,
) -> Result<AnnDataContainer> {
    info!(
        "Calculating site-level mismatch stats for {}>{}",
        ref_base, alt_base
    );

    let results: Vec<(u32, u32, u32)> = (0..adata.n_vars)
        .into_par_iter()
        .map(|site_idx| calculate_site_stats(&adata.var, site_idx, ref_base, alt_base))
        .collect();

    let (ref_counts, alt_counts, others_counts): (Vec<u32>, Vec<u32>, Vec<u32>) =
        results.into_iter().multiunzip();

    // Add columns to var DataFrame - Fixed DataFrame mutations
    let ref_col_name = format!("{}{}_ref", ref_base, alt_base);
    let alt_col_name = format!("{}{}_alt", ref_base, alt_base);
    let others_col_name = format!("{}{}_others", ref_base, alt_base);

    let mismatch_columns: Vec<Column> = vec![
        Series::new(ref_col_name.into(), ref_counts).into_column(),
        Series::new(alt_col_name.into(), alt_counts).into_column(),
        Series::new(others_col_name.into(), others_counts).into_column(),
    ];
    adata.var.hstack_mut(&mismatch_columns)?;

    info!("Site-level mismatch stats calculated");
    Ok(adata)
}

// Helper functions with DataFrame mutation fixes

/// Filter sites by editing type using strand-aware logic
///
/// This function filters genomic sites to only those that are relevant for the
/// specified editing type, taking into account strand-aware processing. The
/// strand-aware filtering considers that the same editing event can be observed
/// from either DNA strand due to complementary base pairing.
///
/// For example, A>G editing on the positive strand appears as T>C editing on the
/// negative strand. This function uses the editing type's strand-aware reference
/// base definitions to identify all potentially relevant sites.
///
/// # Arguments
///
/// * `adata` - The AnnDataContainer containing genomic site information
/// * `editing_type` - The editing type to filter for
///
/// # Returns
///
/// Filtered AnnDataContainer containing only sites relevant for the editing type
fn filter_by_editing_type_strand_aware(
    adata: AnnDataContainer,
    editing_type: &EditingType,
) -> Result<AnnDataContainer> {
    info!(
        "Filtering sites by strand-aware editing type: {:?}",
        editing_type
    );

    // Get the set of reference bases that are valid for this editing type
    // on either DNA strand
    let allowed_ref_bases = editing_type.get_strand_aware_ref_bases();

    let ref_col = adata.var.column("ref")?;
    let filter_mask: Vec<bool> = ref_col
        .str()?
        .par_iter()
        .map(|opt_str| {
            opt_str
                .and_then(|s| s.chars().next())
                .map(|c| allowed_ref_bases.contains(&c))
                .unwrap_or(false)
        })
        .collect();

    let kept_count = filter_mask.par_iter().filter(|&&x| x).count();
    info!(
        "Keeping {} sites after strand-aware editing type filtering",
        kept_count
    );

    apply_site_filter(adata, &filter_mask)
}

fn extract_reference_bases(var_df: &DataFrame) -> Result<Vec<char>> {
    let ref_col = var_df.column("ref")?;
    let ref_bases: Vec<char> = ref_col
        .str()?
        .par_iter()
        .map(|opt_str| opt_str.and_then(|s| s.chars().next()).unwrap_or('N'))
        .collect();

    Ok(ref_bases)
}

fn convert_u32_to_f64_csr(matrix: &CsrMatrix<u32>) -> CsrMatrix<f64> {
    let (row_offsets, col_indices, values) = matrix.csr_data();
    let values_f64: Vec<f64> = values.par_iter().map(|&x| x as f64).collect();

    CsrMatrix::try_from_csr_data(
        matrix.nrows(),
        matrix.ncols(),
        row_offsets.to_vec(),
        col_indices.to_vec(),
        values_f64,
    )
    .expect("Failed to convert u32 to f64 CSR matrix")
}

fn mark_editing_sites(
    mut adata: AnnDataContainer,
    editing_sites: &HashMap<String, u8>,
) -> Result<AnnDataContainer> {
    info!("Marking known editing sites...");

    let is_editing_site: Vec<bool> = adata
        .var_names
        .par_iter()
        .map(|name| editing_sites.contains_key(name))
        .collect();

    let marked_count = is_editing_site.par_iter().filter(|&&x| x).count();
    info!(
        "Marked {} editing sites out of {}",
        marked_count, adata.n_vars
    );

    let filter_column = Series::new("is_editing_site".into(), is_editing_site).into_column();
    adata.var.hstack_mut(&[filter_column])?;

    Ok(adata)
}

fn add_reference_bases(
    mut adata: AnnDataContainer,
    reference: Arc<ReferenceGenome>,
) -> Result<AnnDataContainer> {
    info!("Adding reference bases...");

    let ref_bases: Vec<String> = adata
        .var_names
        .par_iter()
        .map(|name| reference.get_ref_of_pos(name).unwrap_or('N').to_string())
        .collect();

    let n_count = ref_bases.par_iter().filter(|&x| x == "N").count();
    info!(
        "Retrieved {} valid reference bases, {} unknown",
        ref_bases.len() - n_count,
        n_count
    );

    let ref_column = Series::new("ref".into(), ref_bases).into_column();
    adata.var.hstack_mut(&[ref_column])?;

    Ok(adata)
}

fn apply_mismatch_filtering(
    mut adata: AnnDataContainer,
    max_other_threshold: f32,
    min_edited_threshold: f32,
    min_ref_threshold: f32,
) -> Result<AnnDataContainer> {
    info!("Applying mismatch filtering...");

    let mismatch_results: Vec<String> = (0..adata.n_vars)
        .into_par_iter()
        .map(|site_idx| {
            classify_mismatch(
                &adata.var,
                site_idx,
                max_other_threshold,
                min_edited_threshold,
                min_ref_threshold,
            )
        })
        .collect();

    let valid_count = mismatch_results.par_iter().filter(|&x| x != "-").count();
    info!(
        "Found {} valid mismatches out of {} sites",
        valid_count, adata.n_vars
    );

    let mismatch_column = Series::new("Mismatch".into(), mismatch_results).into_column();
    adata.var.hstack_mut(&[mismatch_column])?;

    Ok(adata)
}

fn classify_mismatch(
    var_df: &DataFrame,
    site_idx: usize,
    max_other_threshold: f32,
    min_edited_threshold: f32,
    min_ref_threshold: f32,
) -> String {
    let coverage = get_var_u32(var_df, site_idx, "Coverage").unwrap_or(0) as f32;
    let ref_base_str = get_var_string(var_df, site_idx, "ref").unwrap_or("N".to_string());

    if ref_base_str == "N" || coverage < 1.0 {
        return "-".to_string();
    }

    let ref_char = ref_base_str.chars().next().unwrap();

    // Calculate thresholds
    let other_max = (max_other_threshold * coverage).ceil() as u32;
    let edited_min = (min_edited_threshold * coverage).ceil().max(1.0) as u32;
    let ref_min = (min_ref_threshold * coverage).ceil().max(1.0) as u32;

    // Get counts for each base
    let mut base_counts = HashMap::new();
    for &base in &['A', 'T', 'G', 'C'] {
        let count = get_var_u32(var_df, site_idx, &base.to_string()).unwrap_or(0);
        base_counts.insert(base, count);
    }

    let ref_count = base_counts.get(&ref_char).copied().unwrap_or(0);
    if ref_count < ref_min {
        return "-".to_string();
    }

    // Find the highest non-reference base
    let mut non_ref_bases: Vec<(char, u32)> = base_counts
        .iter()
        .filter(|(&base, _)| base != ref_char)
        .map(|(&base, &count)| (base, count))
        .collect();

    non_ref_bases.sort_by(|a, b| b.1.cmp(&a.1));

    if non_ref_bases.is_empty() {
        return "-".to_string();
    }

    let edited_base = non_ref_bases[0].0;
    let edited_count = non_ref_bases[0].1;
    let others_count: u32 = non_ref_bases[1..].iter().map(|(_, count)| count).sum();

    if others_count > other_max || edited_count < edited_min {
        return "-".to_string();
    }

    format!("{}{}", ref_char, edited_base)
}

fn calculate_site_stats(
    var_df: &DataFrame,
    site_idx: usize,
    ref_base: char,
    alt_base: char,
) -> (u32, u32, u32) {
    let mut base_counts = HashMap::new();

    for &base in &['A', 'T', 'G', 'C'] {
        let count = get_var_u32(var_df, site_idx, &base.to_string()).unwrap_or(0);
        base_counts.insert(base, count);
    }

    let ref_count = base_counts.get(&ref_base).copied().unwrap_or(0);
    let alt_count = base_counts.get(&alt_base).copied().unwrap_or(0);
    let total_count: u32 = base_counts.values().sum();
    let others_count = total_count
        .saturating_sub(ref_count)
        .saturating_sub(alt_count);

    (ref_count, alt_count, others_count)
}

fn get_var_u32(var_df: &DataFrame, row: usize, col: &str) -> Option<u32> {
    var_df.column(col).ok().and_then(|column| {
        column.get(row).ok().and_then(|value| match value {
            AnyValue::UInt32(v) => Some(v),
            AnyValue::Int32(v) => Some(v.max(0) as u32),
            AnyValue::UInt64(v) => Some(v.min(u32::MAX as u64) as u32),
            AnyValue::Int64(v) => Some(v.max(0).min(u32::MAX as i64) as u32),
            _ => None,
        })
    })
}

fn get_var_string(var_df: &DataFrame, row: usize, col: &str) -> Option<String> {
    var_df.column(col).ok().and_then(|column| {
        column.get(row).ok().and_then(|value| match value {
            AnyValue::String(s) => Some(s.to_string()),
            AnyValue::StringOwned(s) => Some(s.to_string()),
            _ => None,
        })
    })
}
