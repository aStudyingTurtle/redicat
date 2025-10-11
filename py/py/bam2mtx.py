"""Convert pileup-derived single nucleotide counts into an AnnData matrix."""

from __future__ import annotations

import argparse
import gzip
import multiprocessing as mp
from collections import Counter
from functools import partial
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import anndata
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

try:
    import pysam  # type: ignore[import]
except ImportError as exc:  # pragma: no cover - environment dependent
    pysam = None  # type: ignore[assignment]
    _PYSAM_IMPORT_ERROR: Optional[ImportError] = exc
else:
    _PYSAM_IMPORT_ERROR = None


def require_pysam() -> Any:
    """Ensure pysam is available before performing BAM operations."""

    if pysam is None:
        raise ImportError(
            'pysam is required for BAM handling. Install it with "pip install pysam".'
        ) from _PYSAM_IMPORT_ERROR
    return pysam


def process_pos(
    chromosome: str,
    position: int,
    bam_file: Any,
    cell_barcodes: Dict[str, bool],
    stranded: bool = True,
) -> Dict[str, List[int]]:
    """Accumulate counts per cell barcode and UMI for a single genomic position."""

    observed: Dict[str, Dict[str, str]] = {}

    for pileup_column in bam_file.pileup(
        contig=chromosome,
        start=position - 1,
        end=position,
        truncate=True,
        min_mapping_quality=255,
        min_base_quality=30,
    ):
        for pileup_read in pileup_column.pileups:
            if pileup_read.is_del or pileup_read.is_refskip:
                continue

            read = pileup_read.alignment
            base_at_pos = read.query_sequence[pileup_read.query_position].upper()
            umi = read.get_tag('UB') if read.has_tag('UB') else '-'
            cb_tag = read.get_tag('CB') if read.has_tag('CB') else '-'
            cb_tag = cb_tag.split('-')[0]

            if cb_tag not in cell_barcodes:
                continue

            observed.setdefault(cb_tag, {})

            if stranded:
                strand = '+' if not read.is_reverse else '-'
                base_descriptor = f'{base_at_pos}{strand}'
            else:
                base_descriptor = base_at_pos

            previous = observed[cb_tag].get(umi)
            if previous is None:
                observed[cb_tag][umi] = base_descriptor
            elif previous != base_descriptor:
                observed[cb_tag][umi] = '-'

    result: Dict[str, List[int]] = {}

    for cb_tag, umi_dict in observed.items():
        counts = Counter(value.upper() for value in umi_dict.values())
        if stranded:
            result[cb_tag] = [
                counts.get('A+', 0),
                counts.get('T+', 0),
                counts.get('G+', 0),
                counts.get('C+', 0),
                counts.get('A-', 0),
                counts.get('T-', 0),
                counts.get('G-', 0),
                counts.get('C-', 0),
            ]
        else:
            result[cb_tag] = [
                counts.get('A', 0),
                counts.get('T', 0),
                counts.get('G', 0),
                counts.get('C', 0),
            ]

    return result


def get_cell_barcodes(path: str) -> Dict[str, bool]:
    """Load valid cell barcodes, supporting optional gzip compression."""

    barcode_path = Path(path)
    if not barcode_path.exists():
        raise FileNotFoundError(f'Barcode file not found: {path}')

    if barcode_path.suffix == '.gz':
        with gzip.open(barcode_path, 'rt', encoding='utf-8') as handle:
            entries = handle.read().splitlines()
    else:
        entries = barcode_path.read_text(encoding='utf-8').splitlines()

    return {entry.split('-')[0]: True for entry in entries if entry}


def split_dataframe(df: pd.DataFrame, chunk_size: int = 2000) -> List[pd.DataFrame]:
    """Partition a DataFrame into equally sized chunks."""

    total_rows = df.shape[0]
    num_chunks = (total_rows + chunk_size - 1) // chunk_size
    return [df.iloc[i * chunk_size : (i + 1) * chunk_size] for i in range(num_chunks)]


def get_sites_to_process(path: str, chunk_size: int = 2000) -> List[pd.DataFrame]:
    """Read a TSV/CSV file describing loci to process and split into chunks."""

    site_path = Path(path)
    if not site_path.exists():
        raise FileNotFoundError(f'Site definition file not found: {path}')

    df = pd.read_csv(site_path, sep='\t' if site_path.suffix == '.tsv' else ',', usecols=['REF', 'POS'])
    return split_dataframe(df, chunk_size)


def bam2mtx_worker(
    chunk_df: pd.DataFrame,
    bam_file_path: str,
    cell_barcodes: Dict[str, bool],
    stranded: bool = True,
) -> Dict[str, Dict[str, List[int]]]:
    """Worker function executed in parallel to accumulate counts for a chunk."""

    pysam_module = require_pysam()
    bam_file = pysam_module.AlignmentFile(bam_file_path, 'rb')
    result: Dict[str, Dict[str, List[int]]] = {}

    chrom_list = chunk_df['REF'].tolist()
    pos_list = chunk_df['POS'].tolist()

    for chrom, pos in zip(chrom_list, pos_list):
        chrom_pos = f'{chrom}:{pos}'
        result[chrom_pos] = process_pos(chrom, int(pos), bam_file, cell_barcodes, stranded)

    bam_file.close()
    return result


def res_to_anndata(result_dict: Dict[str, Dict[str, List[int]]], stranded: bool = True) -> anndata.AnnData:
    """Convert nested dictionaries of counts into an AnnData object."""

    A1: List[int] = []
    T1: List[int] = []
    G1: List[int] = []
    C1: List[int] = []
    A0: List[int] = []
    T0: List[int] = []
    G0: List[int] = []
    C0: List[int] = []

    obs_index: Dict[str, int] = {}
    var_index: Dict[str, int] = {}
    obs_coords: List[int] = []
    var_coords: List[int] = []

    for pos, barcode_dict in result_dict.items():
        col_idx = var_index.setdefault(pos, len(var_index))
        for cell, counts in barcode_dict.items():
            row_idx = obs_index.setdefault(cell, len(obs_index))
            obs_coords.append(row_idx)
            var_coords.append(col_idx)

            A1.append(counts[0])
            T1.append(counts[1])
            G1.append(counts[2])
            C1.append(counts[3])

            if stranded:
                A0.append(counts[4])
                T0.append(counts[5])
                G0.append(counts[6])
                C0.append(counts[7])

    shape = (len(obs_index), len(var_index))
    sparse_A1 = csc_matrix((A1, (obs_coords, var_coords)), shape=shape, dtype=np.uint16)
    sparse_T1 = csc_matrix((T1, (obs_coords, var_coords)), shape=shape, dtype=np.uint16)
    sparse_G1 = csc_matrix((G1, (obs_coords, var_coords)), shape=shape, dtype=np.uint16)
    sparse_C1 = csc_matrix((C1, (obs_coords, var_coords)), shape=shape, dtype=np.uint16)

    obs_names = [None] * len(obs_index)
    for name, idx in obs_index.items():
        obs_names[idx] = name

    var_names = [None] * len(var_index)
    for name, idx in var_index.items():
        var_names[idx] = name

    obs_df = pd.DataFrame(index=pd.Index(obs_names, name='cell'))
    var_df = pd.DataFrame(index=pd.Index(var_names, name='site'))

    adata = anndata.AnnData(X=sparse_A1, obs=obs_df, var=var_df)
    adata.layers['A1'] = sparse_A1
    adata.layers['T1'] = sparse_T1
    adata.layers['G1'] = sparse_G1
    adata.layers['C1'] = sparse_C1

    if stranded:
        sparse_A0 = csc_matrix((A0, (obs_coords, var_coords)), shape=shape, dtype=np.uint16)
        sparse_T0 = csc_matrix((T0, (obs_coords, var_coords)), shape=shape, dtype=np.uint16)
        sparse_G0 = csc_matrix((G0, (obs_coords, var_coords)), shape=shape, dtype=np.uint16)
        sparse_C0 = csc_matrix((C0, (obs_coords, var_coords)), shape=shape, dtype=np.uint16)

        adata.layers['A0'] = sparse_A0
        adata.layers['T0'] = sparse_T0
        adata.layers['G0'] = sparse_G0
        adata.layers['C0'] = sparse_C0

    return adata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Convert BAM pileups into strand-aware count matrices.')
    parser.add_argument('--threads', type=int, default=1, help='Number of worker processes to spawn.')
    parser.add_argument('--barcodes', required=True, help='Path to the whitelist of cell barcodes (optionally gzipped).')
    parser.add_argument('--bam', required=True, help='Aligned BAM file produced by the sequencing experiment.')
    parser.add_argument('--sites', required=True, help='TSV/CSV file with columns REF and POS defining loci to extract.')
    parser.add_argument('--output-h5ad', required=True, help='Destination path for the generated AnnData object.')
    parser.add_argument('--chunk-size', type=int, default=2000, help='Number of loci to process per worker chunk.')
    parser.add_argument('--unstranded', action='store_true', help='Treat reads as unstranded when aggregating counts.')
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if pysam is None:
        raise ImportError(
            'pysam could not be imported. Please install pysam to run bam2mtx.'
        ) from _PYSAM_IMPORT_ERROR

    barcodes = get_cell_barcodes(args.barcodes)
    chunks = get_sites_to_process(args.sites, chunk_size=args.chunk_size)

    worker = partial(
        bam2mtx_worker,
        bam_file_path=args.bam,
        cell_barcodes=barcodes,
        stranded=not args.unstranded,
    )

    if args.threads > 1:
        with mp.Pool(args.threads) as pool:
            partial_results = pool.map(worker, chunks)
    else:
        partial_results = [worker(chunk) for chunk in chunks]

    merged: Dict[str, Dict[str, List[int]]] = {}
    for chunk_dict in partial_results:
        merged.update(chunk_dict)

    adata = res_to_anndata(merged, stranded=not args.unstranded)
    adata.write_h5ad(args.output_h5ad, compression='gzip')


if __name__ == '__main__':
    main()
