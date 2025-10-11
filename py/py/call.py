"""Utilities for computing RNA editing metrics from base-level AnnData matrices."""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Literal, Optional, Tuple

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from sklearn.preprocessing import OneHotEncoder


def calculate_coverage(base_matrix: anndata.AnnData) -> scipy.sparse.spmatrix:
    """Sum all strand-specific nucleotide layers into a coverage matrix."""
    coverage_matrix = (
        base_matrix.layers['A0']
        + base_matrix.layers['A1']
        + base_matrix.layers['C0']
        + base_matrix.layers['C1']
        + base_matrix.layers['G0']
        + base_matrix.layers['G1']
        + base_matrix.layers['T0']
        + base_matrix.layers['T1']
    )
    return coverage_matrix


def filter_base_matrix(
    base_matrix: anndata.AnnData,
    min_coverage: int = 5,
) -> anndata.AnnData:
    """Filter positions by chromosome prefix and minimum coverage."""

    base_matrix.layers['coverage'] = calculate_coverage(base_matrix)
    coverage_totals = base_matrix.layers['coverage'].sum(axis=0)
    base_matrix.var['Coverage'] = np.asarray(coverage_totals).reshape(-1, 1)

    selected_pos = base_matrix.var_names[
        base_matrix.var.index.str.contains('^chr')
        & (base_matrix.var.Coverage >= min_coverage)
    ]
    return base_matrix[:, selected_pos]


def calculate_sparsity(coverage_matrix: scipy.sparse.spmatrix) -> None:
    """Print sparsity statistics for the provided coverage matrix."""

    print('Shape:', coverage_matrix.shape)
    non_zero_values = coverage_matrix.data
    quartiles = np.percentile(non_zero_values, [90, 95, 99, 99.9, 99.99, 99.999])
    print('90% 95% 99% 99.9% 99.99% 99.999% quartiles:', quartiles)
    total_elements = coverage_matrix.shape[0] * coverage_matrix.shape[1]
    non_zero_elements = coverage_matrix.count_nonzero()
    sparsity = non_zero_elements / total_elements
    print('稀疏度：', sparsity)


def count_base_site_levels(base_matrix: anndata.AnnData) -> anndata.AnnData:
    """Aggregate strand-specific base counts into total layers."""

    base_matrix.var['A_P'] = base_matrix.layers['A1'].sum(axis=0).reshape(-1, 1)
    base_matrix.var['G_P'] = base_matrix.layers['G1'].sum(axis=0).reshape(-1, 1)
    base_matrix.var['T_P'] = base_matrix.layers['T1'].sum(axis=0).reshape(-1, 1)
    base_matrix.var['C_P'] = base_matrix.layers['C1'].sum(axis=0).reshape(-1, 1)
    base_matrix.var['A_N'] = base_matrix.layers['A0'].sum(axis=0).reshape(-1, 1)
    base_matrix.var['G_N'] = base_matrix.layers['G0'].sum(axis=0).reshape(-1, 1)
    base_matrix.var['T_N'] = base_matrix.layers['T0'].sum(axis=0).reshape(-1, 1)
    base_matrix.var['C_N'] = base_matrix.layers['C0'].sum(axis=0).reshape(-1, 1)

    base_matrix.var['A'] = base_matrix.var['A_P'] + base_matrix.var['A_N']
    base_matrix.var['G'] = base_matrix.var['G_P'] + base_matrix.var['G_N']
    base_matrix.var['T'] = base_matrix.var['T_P'] + base_matrix.var['T_N']
    base_matrix.var['C'] = base_matrix.var['C_P'] + base_matrix.var['C_N']
    base_matrix.var['CoverageP'] = (
        base_matrix.var['A_P']
        + base_matrix.var['G_P']
        + base_matrix.var['T_P']
        + base_matrix.var['C_P']
    )
    base_matrix.var['CoverageN'] = (
        base_matrix.var['A_N']
        + base_matrix.var['G_N']
        + base_matrix.var['T_N']
        + base_matrix.var['C_N']
    )

    return base_matrix


def test_strand_specific_matrix(base_matrix: anndata.AnnData) -> Dict[str, int]:
    """Summarize strand-specific behavior of the matrix."""

    return {
        'Shape': base_matrix.shape[1],
        'Strict Strand Specific': base_matrix.var[
            (base_matrix.var.CoverageP == 0) | (base_matrix.var.CoverageN == 0)
        ].shape[0],
        'Strand Specific 10X': base_matrix.var[
            (base_matrix.var.CoverageP // base_matrix.var.CoverageN >= 10)
            | (base_matrix.var.CoverageN // base_matrix.var.CoverageP >= 10)
        ].shape[0],
        'Strand Specific 5X': base_matrix.var[
            (base_matrix.var.CoverageP // base_matrix.var.CoverageN >= 5)
            | (base_matrix.var.CoverageN // base_matrix.var.CoverageP >= 5)
        ].shape[0],
        'Strand Specific 2X': base_matrix.var[
            (base_matrix.var.CoverageP // base_matrix.var.CoverageN >= 2)
            | (base_matrix.var.CoverageN // base_matrix.var.CoverageP >= 2)
        ].shape[0],
    }


def downsample_by_group(
    adata: anndata.AnnData,
    group_col: str = 'gene',
    n_samples: int = 50,
    random_state: int = 2025,
) -> anndata.AnnData:
    """Downsample cells per group without materializing dense arrays."""

    if group_col not in adata.obs.columns:
        raise ValueError(f"'{group_col}' not found in adata.obs")

    rng = np.random.default_rng(random_state)
    groups = adata.obs[group_col].unique()
    indices_to_keep: list[int] = []

    for group in groups:
        group_mask = adata.obs[group_col] == group
        group_indices = np.where(group_mask)[0]

        if len(group_indices) <= n_samples:
            indices_to_keep.extend(group_indices)
        else:
            sampled_indices = rng.choice(group_indices, size=n_samples, replace=False)
            indices_to_keep.extend(sampled_indices)

    indices_to_keep = np.sort(indices_to_keep)
    return adata[indices_to_keep].copy()


def _normalize_metadata_entry(entry: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if entry is None:
        return {}
    normalized = {str(k).lower(): v for k, v in entry.items()}
    return normalized


def infer_ref_from_var_name(var_name: str) -> Optional[str]:
    """Infer the reference base from a position identifier when possible."""

    tokens = [token.strip() for token in var_name.replace('|', ':').replace('-', ':').split(':') if token]
    for token in reversed(tokens):
        upper = token.upper()
        if upper in {'A', 'C', 'G', 'T'}:
            return upper
        if len(upper) == 2 and upper[0] in {'A', 'C', 'G', 'T'}:
            return upper[0]
    return None


def load_editing_site_dict(path: Optional[str]) -> Dict[str, Dict[str, Any]]:
    """Load per-site annotations from a JSON or tabular file."""

    if path is None:
        return {}

    site_path = Path(path)
    if not site_path.exists():
        raise FileNotFoundError(f'Editing-site annotation file not found: {path}')

    suffix = site_path.suffix.lower()
    if suffix in {'.json', '.jsonl'}:
        with site_path.open('r', encoding='utf-8') as handle:
            data = json.load(handle)
        if isinstance(data, dict):
            return {str(key): (value if isinstance(value, dict) else {}) for key, value in data.items()}
        if isinstance(data, list):
            result: Dict[str, Dict[str, Any]] = {}
            for item in data:
                if isinstance(item, dict):
                    key = item.get('site') or item.get('id') or item.get('name')
                    if key is not None:
                        result[str(key)] = {k: v for k, v in item.items() if k not in {'site', 'id', 'name'}}
            return result
        raise ValueError('Unsupported JSON structure for editing annotations.')

    if suffix in {'.tsv', '.csv', '.txt'}:
        sep = ',' if suffix == '.csv' else '\t'
        df = pd.read_csv(site_path, sep=sep)
    elif suffix in {'.parquet', '.pq'}:
        df = pd.read_parquet(site_path)
    else:
        raise ValueError(f'Unsupported annotation format: {suffix}')

    if df.index.name is None or df.index.name == '':
        for candidate in ('site', 'id', 'name', 'position', 'pos'):
            if candidate in df.columns:
                df = df.set_index(candidate)
                break
    if df.index.name is None or df.index.name == '':
        raise ValueError('Annotation table must contain a column identifying each site.')

    df = df.applymap(lambda x: x if not (isinstance(x, float) and np.isnan(x)) else None)
    return {str(idx): row.dropna().to_dict() for idx, row in df.iterrows()}


def annotate_var(
    base_matrix: anndata.AnnData,
    editing_site_dict: Dict[str, Dict[str, Any]],
    min_coverage: int = 5,
) -> anndata.AnnData:
    """Annotate variant metadata and attach per-site filters."""

    base_matrix.var['filter'] = base_matrix.var.index.map(lambda x: x in editing_site_dict)
    base_matrix = filter_base_matrix(base_matrix, min_coverage=min_coverage)
    base_matrix = count_base_site_levels(base_matrix)

    editing_metadata = pd.DataFrame.from_dict(editing_site_dict, orient='index') if editing_site_dict else pd.DataFrame()
    if not editing_metadata.empty:
        editing_metadata.index = editing_metadata.index.astype(str)
        base_matrix.var = base_matrix.var.join(editing_metadata, how='left')

    def resolve_ref(var_name: str) -> str:
        entry = _normalize_metadata_entry(editing_site_dict.get(var_name))
        ref_candidate = entry.get('ref') or entry.get('reference') or entry.get('ref_base')
        if ref_candidate is None:
            ref_candidate = infer_ref_from_var_name(var_name)
        if ref_candidate is None:
            raise ValueError(f'Unable to determine reference base for {var_name}.')
        return str(ref_candidate).upper()

    base_matrix.var['ref'] = base_matrix.var.index.map(resolve_ref)

    def build_filter_tuple(var_name: str) -> Tuple[bool, bool, int, int]:
        entry = _normalize_metadata_entry(editing_site_dict.get(var_name))
        stranded = bool(entry.get('stranded', entry.get('strand_specific', False)))
        mismatch = bool(entry.get('mismatch', entry.get('is_mismatch', False)))
        mismatch_p = int(entry.get('mismatch_p', entry.get('mismatch_plus', 0) or 0))
        mismatch_n = int(entry.get('mismatch_n', entry.get('mismatch_minus', 0) or 0))
        return stranded, mismatch, mismatch_p, mismatch_n

    filter_values = pd.DataFrame(
        [build_filter_tuple(name) for name in base_matrix.var_names],
        index=base_matrix.var_names,
        columns=['Stranded', 'Mismatch', 'Mismatch_P', 'Mismatch_N'],
    )
    base_matrix.var[['Stranded', 'Mismatch', 'Mismatch_P', 'Mismatch_N']] = filter_values

    return base_matrix


def filter_by_gene(
    base_matrix: anndata.AnnData,
    gene: str = 'DDX39B',
    nt: str = 'non-targeting',
    umi_threshold: int = 5000,
    gene_col: str = 'gene',
    umi_col: str = 'UMI_count',
    kd_threshold: float = -0.5,
    nt_threshold: float = 0.5,
) -> anndata.AnnData:
    """Filter cells by gene-level perturbation status."""

    umi_mask = base_matrix.obs[umi_col] >= umi_threshold
    gene_mask = (base_matrix.obs[gene_col] == gene) & (base_matrix.obs[gene] < kd_threshold)
    nt_mask = (base_matrix.obs[gene_col] == nt) & (base_matrix.obs[gene] > nt_threshold)

    combined_mask = umi_mask & (gene_mask | nt_mask)
    return base_matrix[combined_mask]


def get_base_mask_optimized(categories: Iterable[str], ref_array: np.ndarray) -> scipy.sparse.csr_matrix:
    """Build a one-hot encoded mask for reference bases."""

    try:
        encoder = OneHotEncoder(categories=[list(categories)], dtype=np.uint8, sparse_output=True)
    except TypeError:
        encoder = OneHotEncoder(categories=[list(categories)], dtype=np.uint8, sparse=True)

    mask = encoder.fit_transform(ref_array)
    if not scipy.sparse.issparse(mask):
        mask = scipy.sparse.csr_matrix(mask)
    return mask.tocsr()


def get_matrix_with_mask_optimized(read_matrix: scipy.sparse.spmatrix, mask: scipy.sparse.spmatrix) -> scipy.sparse.csr_matrix:
    """Apply a mask to a read matrix while retaining sparsity."""

    if scipy.sparse.issparse(read_matrix):
        result = read_matrix.multiply(mask)
        return scipy.sparse.csr_matrix(result.sum(axis=2), dtype=np.uint16)

    scaled = mask.multiply(read_matrix)
    return scipy.sparse.csr_matrix(np.sum(scaled, axis=2), dtype=np.uint16)


def calculate_cei(sc_adata: anndata.AnnData) -> anndata.AnnData:
    """Compute the composite editing index per observation."""

    ref_alt_sum = sc_adata.obs['ref'] + sc_adata.obs['alt']
    sc_adata.obs['CEI'] = np.divide(
        sc_adata.obs['alt'],
        ref_alt_sum,
        out=np.zeros_like(sc_adata.obs['alt'], dtype=float),
        where=ref_alt_sum != 0,
    )
    return sc_adata


def get_ref_alt_matrix_optimized(
    sc_adata: anndata.AnnData,
    editing_type: Literal['AG', 'CU'],
) -> anndata.AnnData:
    """Project strand-specific matrices into ref/alt space using sparse ops."""

    if editing_type == 'AG':
        pos_to_keep = (sc_adata.var.ref == 'A') | (sc_adata.var.ref == 'T')
    elif editing_type == 'CU':
        pos_to_keep = (sc_adata.var.ref == 'C') | (sc_adata.var.ref == 'G')
    else:
        raise ValueError("editing_type must be 'AG' or 'CU'")

    sc_adata = sc_adata[:, pos_to_keep].copy()
    ref_array = sc_adata.var.ref.values.reshape(-1, 1)

    ref_mask = get_base_mask_optimized(['A', 'C', 'G', 'T'], ref_array)
    alt_mask = get_base_mask_optimized(['G', 'T', 'A', 'C'], ref_array)
    ones = scipy.sparse.csr_matrix(np.ones(ref_mask.shape, dtype=np.uint8))
    others_mask = ones - ref_mask - alt_mask

    sc_adata.layers['ref'] = scipy.sparse.csr_matrix(sc_adata.shape, dtype=np.uint16)
    sc_adata.layers['alt'] = scipy.sparse.csr_matrix(sc_adata.shape, dtype=np.uint16)
    sc_adata.layers['others'] = scipy.sparse.csr_matrix(sc_adata.shape, dtype=np.uint16)

    bases = ['A', 'C', 'G', 'T']
    for idx, base in enumerate(bases):
        print(f'Processing {base}...')
        base_matrix = (sc_adata.layers[f'{base}0'] + sc_adata.layers[f'{base}1']).tocsr()

        del sc_adata.layers[f'{base}0']
        del sc_adata.layers[f'{base}1']

        ref_positions = ref_mask[:, idx].toarray().ravel()
        alt_positions = alt_mask[:, idx].toarray().ravel()
        others_positions = others_mask[:, idx].toarray().ravel()

        ref_contribution = base_matrix.multiply(ref_positions.reshape(1, -1))
        alt_contribution = base_matrix.multiply(alt_positions.reshape(1, -1))
        others_contribution = base_matrix.multiply(others_positions.reshape(1, -1))

        sc_adata.layers['ref'] += ref_contribution
        sc_adata.layers['alt'] += alt_contribution
        sc_adata.layers['others'] += others_contribution

        del base_matrix, ref_contribution, alt_contribution, others_contribution
        del ref_positions, alt_positions, others_positions
        gc.collect()

        print(f'Finished processing {base}')

    sc_adata.obs['ref'] = np.asarray(sc_adata.layers['ref'].sum(axis=1)).ravel()
    sc_adata.obs['alt'] = np.asarray(sc_adata.layers['alt'].sum(axis=1)).ravel()
    sc_adata.obs['others'] = np.asarray(sc_adata.layers['others'].sum(axis=1)).ravel()

    sc_adata.X = sc_adata.layers['alt']

    del ref_mask, alt_mask, others_mask, ref_array
    gc.collect()

    return sc_adata


def get_ref_alt_matrix_chunked(
    sc_adata: anndata.AnnData,
    editing_type: Literal['AG', 'CU'],
    chunk_size: int = 1000,
) -> anndata.AnnData:
    """Chunked version of :func:`get_ref_alt_matrix_optimized` for large data."""

    if editing_type == 'AG':
        pos_to_keep = (sc_adata.var.ref == 'A') | (sc_adata.var.ref == 'T')
    elif editing_type == 'CU':
        pos_to_keep = (sc_adata.var.ref == 'C') | (sc_adata.var.ref == 'G')
    else:
        raise ValueError("editing_type must be 'AG' or 'CU'")

    sc_adata = sc_adata[:, pos_to_keep].copy()
    n_obs = sc_adata.n_obs

    ref_sums = np.zeros(n_obs, dtype=np.uint32)
    alt_sums = np.zeros(n_obs, dtype=np.uint32)
    others_sums = np.zeros(n_obs, dtype=np.uint32)

    for start_idx in range(0, n_obs, chunk_size):
        end_idx = min(start_idx + chunk_size, n_obs)
        chunk_adata = sc_adata[start_idx:end_idx].copy()
        chunk_result = get_ref_alt_matrix_optimized(chunk_adata, editing_type)

        ref_sums[start_idx:end_idx] = chunk_result.obs['ref']
        alt_sums[start_idx:end_idx] = chunk_result.obs['alt']
        others_sums[start_idx:end_idx] = chunk_result.obs['others']

        del chunk_adata, chunk_result
        gc.collect()

        print(
            f'Processed chunk {start_idx // chunk_size + 1}/{(n_obs - 1) // chunk_size + 1}'
        )

    sc_adata.obs['ref'] = ref_sums
    sc_adata.obs['alt'] = alt_sums
    sc_adata.obs['others'] = others_sums

    return sc_adata


def get_memory_usage() -> float:
    """Return current process memory usage in MB."""

    import psutil

    process = psutil.Process()
    return process.memory_info().rss / 1024 / 1024


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Compute ref/alt matrices and CEI scores.')
    parser.add_argument('--input-h5ad', required=True, help='Input AnnData file with base layers.')
    parser.add_argument('--editing-sites', help='Annotation file (JSON/CSV/TSV/Parquet) describing editing sites.')
    parser.add_argument('--output-h5ad', required=True, help='Destination AnnData file for results.')
    parser.add_argument(
        '--editing-type',
        choices=['AG', 'CU'],
        default='AG',
        help='Editing type to retain when constructing ref/alt matrices.',
    )
    parser.add_argument(
        '--min-coverage',
        type=int,
        default=5,
        help='Minimum total coverage required to retain a genomic position.',
    )
    parser.add_argument(
        '--chunk-size',
        type=int,
        default=0,
        help='Optional chunk size for processing very large matrices (0 disables chunked mode).',
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    editing_site_dict = load_editing_site_dict(args.editing_sites)
    matrix_adata = sc.read(args.input_h5ad)

    matrix_adata = annotate_var(matrix_adata, editing_site_dict, min_coverage=args.min_coverage)

    if args.chunk_size and args.chunk_size > 0:
        matrix_adata = get_ref_alt_matrix_chunked(matrix_adata, args.editing_type, chunk_size=args.chunk_size)
    else:
        matrix_adata = get_ref_alt_matrix_optimized(matrix_adata, args.editing_type)

    matrix_adata = calculate_cei(matrix_adata)
    matrix_adata.write(args.output_h5ad)


if __name__ == '__main__':
    main()
