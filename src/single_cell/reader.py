"""Readers for single-cell and single-nucleus RNA-seq assets (H5AD, 10x Matrix Market, RDS rejection)."""

import gzip
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
from scipy import io, sparse

from src.qc.result import Finding, Severity
from src.single_cell.models import InvalidRawCountMatrixException, SingleCellMatrix


class UnsupportedSingleCellFormatException(Exception):
    """Raised when an unsupported binary format (e.g. .rds R serialization) is encountered."""

    pass


def _decode_h5_value(val: Any) -> Any:
    """Decode bytes or arrays of bytes from HDF5."""
    if isinstance(val, (bytes, np.bytes_)):
        return val.decode("utf-8", errors="replace")
    if isinstance(val, np.ndarray) and val.dtype.kind in ("S", "O"):
        return [_decode_h5_value(v) for v in val]
    return val


def read_h5ad(file_path: Path | str) -> tuple[SingleCellMatrix, pd.DataFrame, pd.DataFrame]:
    """Read an AnnData H5AD file via h5py without requiring the full anndata package.

    Extracts:
    - X (expression matrix; transposed to genes x cells if stored as standard AnnData cells x genes)
    - obs (cell annotations dataframe)
    - var (gene annotations dataframe)
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"H5AD file not found: {path}")

    with h5py.File(path, "r") as h5:
        # 1. Read obs (cell metadata)
        obs_dict: dict[str, Any] = {}
        obs_index = None

        if "obs" in h5:
            obs_group = h5["obs"]
            if isinstance(obs_group, h5py.Group):
                for key in obs_group.keys():
                    item = obs_group[key]
                    if isinstance(item, h5py.Dataset):
                        data = item[()]
                        if isinstance(data, np.ndarray) and data.dtype.kind in ("S", "O", "U"):
                            obs_dict[key] = [
                                v.decode("utf-8", errors="replace") if isinstance(v, bytes) else str(v)
                                for v in data
                            ]
                        else:
                            obs_dict[key] = data
                    elif isinstance(item, h5py.Group) and "categories" in item and "codes" in item:
                        # Categorical column in AnnData
                        cats = [
                            c.decode("utf-8", errors="replace") if isinstance(c, bytes) else str(c)
                            for c in item["categories"][()]
                        ]
                        codes = item["codes"][()]
                        obs_dict[key] = [cats[c] if 0 <= c < len(cats) else None for c in codes]

                # Check index
                if "_index" in obs_group.attrs:
                    idx_name = obs_group.attrs["_index"]
                    if idx_name in obs_dict:
                        obs_index = obs_dict.pop(idx_name)
                elif "_index" in obs_dict:
                    obs_index = obs_dict.pop("_index")
                elif "index" in obs_dict:
                    obs_index = obs_dict.pop("index")
            elif isinstance(obs_group, h5py.Dataset):
                # Compound dataset
                arr = obs_group[()]
                for name in arr.dtype.names:
                    obs_dict[name] = arr[name]

        # 2. Read var (gene metadata)
        var_dict: dict[str, Any] = {}
        var_index = None

        if "var" in h5:
            var_group = h5["var"]
            if isinstance(var_group, h5py.Group):
                for key in var_group.keys():
                    item = var_group[key]
                    if isinstance(item, h5py.Dataset):
                        data = item[()]
                        if isinstance(data, np.ndarray) and data.dtype.kind in ("S", "O", "U"):
                            var_dict[key] = [
                                v.decode("utf-8", errors="replace") if isinstance(v, bytes) else str(v)
                                for v in data
                            ]
                        else:
                            var_dict[key] = data
                    elif isinstance(item, h5py.Group) and "categories" in item and "codes" in item:
                        cats = [
                            c.decode("utf-8", errors="replace") if isinstance(c, bytes) else str(c)
                            for c in item["categories"][()]
                        ]
                        codes = item["codes"][()]
                        var_dict[key] = [cats[c] if 0 <= c < len(cats) else None for c in codes]

                if "_index" in var_group.attrs:
                    idx_name = var_group.attrs["_index"]
                    if idx_name in var_dict:
                        var_index = var_dict.pop(idx_name)
                elif "_index" in var_dict:
                    var_index = var_dict.pop("_index")
                elif "index" in var_dict:
                    var_index = var_dict.pop("index")
            elif isinstance(var_group, h5py.Dataset):
                arr = var_group[()]
                for name in arr.dtype.names:
                    var_dict[name] = arr[name]

        # 3. Read X matrix
        if "X" not in h5:
            raise ValueError(f"H5AD file missing root 'X' matrix: {path}")

        x_item = h5["X"]
        if isinstance(x_item, h5py.Group):
            # Sparse matrix: csr_matrix or csc_matrix
            data = x_item["data"][()]
            indices = x_item["indices"][()]
            indptr = x_item["indptr"][()]
            shape = tuple(x_item.attrs.get("shape", ()))
            if not shape:
                shape = (len(indptr) - 1, int(np.max(indices)) + 1 if len(indices) else 0)

            encoding = x_item.attrs.get("encoding-type", "csr_matrix")
            if isinstance(encoding, bytes):
                encoding = encoding.decode("utf-8")

            if "csc" in encoding.lower():
                mat = sparse.csc_matrix((data, indices, indptr), shape=shape)
            else:
                mat = sparse.csr_matrix((data, indices, indptr), shape=shape)
        else:
            mat = x_item[()]

        # Determine dimensions and orientations
        n_obs = len(obs_index) if obs_index is not None else (len(next(iter(obs_dict.values()))) if obs_dict else mat.shape[0])
        n_vars = len(var_index) if var_index is not None else (len(next(iter(var_dict.values()))) if var_dict else mat.shape[1])

        # In standard AnnData, X is shape (n_obs, n_vars): cells x genes.
        # Our SingleCellMatrix requires shape (n_genes, n_cells): genes x cells.
        if mat.shape == (n_obs, n_vars):
            mat_genes_by_cells = mat.T
        elif mat.shape == (n_vars, n_obs):
            mat_genes_by_cells = mat
        else:
            raise ValueError(
                f"Matrix shape {mat.shape} incompatible with obs ({n_obs}) and var ({n_vars})."
            )

        if obs_index is None:
            obs_index = [f"cell_{i+1}" for i in range(mat_genes_by_cells.shape[1])]
        if var_index is None:
            var_index = [f"gene_{i+1}" for i in range(mat_genes_by_cells.shape[0])]

        obs_df = pd.DataFrame(obs_dict, index=pd.Index([str(x) for x in obs_index], name="cell_id"))
        var_df = pd.DataFrame(var_dict, index=pd.Index([str(x) for x in var_index], name="gene_id"))

        single_cell_mat = SingleCellMatrix(
            gene_ids=tuple(var_df.index),
            cell_ids=tuple(obs_df.index),
            matrix=mat_genes_by_cells,
        )

        return single_cell_mat, obs_df, var_df


def read_10x_matrix_market(dir_or_mtx_path: Path | str) -> tuple[SingleCellMatrix, pd.DataFrame, pd.DataFrame]:
    """Read 10x Genomics Matrix Market files (matrix.mtx[.gz], features.tsv[.gz], barcodes.tsv[.gz]).

    Returns SingleCellMatrix (genes x cells), obs_df (barcodes), var_df (features).
    """
    path = Path(dir_or_mtx_path)
    if path.is_file():
        directory = path.parent
        mtx_file = path
    else:
        directory = path
        candidates = list(directory.glob("matrix.mtx*"))
        if not candidates:
            raise FileNotFoundError(f"No matrix.mtx or matrix.mtx.gz found in {directory}")
        mtx_file = candidates[0]

    # Find barcodes file
    barcode_candidates = list(directory.glob("barcodes.tsv*"))
    if not barcode_candidates:
        raise FileNotFoundError(f"No barcodes.tsv or barcodes.tsv.gz found in {directory}")
    barcodes_file = barcode_candidates[0]

    # Find features / genes file
    feature_candidates = list(directory.glob("features.tsv*")) or list(directory.glob("genes.tsv*"))
    if not feature_candidates:
        raise FileNotFoundError(f"No features.tsv or genes.tsv found in {directory}")
    features_file = feature_candidates[0]

    # 1. Read barcodes
    if barcodes_file.suffix == ".gz":
        with gzip.open(barcodes_file, "rt", encoding="utf-8") as f:
            barcodes = [line.strip() for line in f if line.strip()]
    else:
        with barcodes_file.open("r", encoding="utf-8") as f:
            barcodes = [line.strip() for line in f if line.strip()]

    # 2. Read features
    feature_rows: list[list[str]] = []
    if features_file.suffix == ".gz":
        with gzip.open(features_file, "rt", encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    feature_rows.append(line.strip().split("\t"))
    else:
        with features_file.open("r", encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    feature_rows.append(line.strip().split("\t"))

    gene_ids = [r[0] for r in feature_rows]
    gene_symbols = [r[1] if len(r) > 1 else r[0] for r in feature_rows]
    feature_types = [r[2] if len(r) > 2 else "Gene Expression" for r in feature_rows]

    var_df = pd.DataFrame(
        {"gene_symbol": gene_symbols, "feature_type": feature_types},
        index=pd.Index(gene_ids, name="gene_id"),
    )
    obs_df = pd.DataFrame(index=pd.Index(barcodes, name="cell_id"))

    # 3. Read MTX (scipy.io.mmread handles gzip or raw file)
    mat = io.mmread(str(mtx_file))
    if not sparse.issparse(mat):
        mat = sparse.csr_matrix(mat)
    else:
        mat = mat.tocsr()

    # 10x standard: matrix is (n_genes, n_cells)
    if mat.shape != (len(gene_ids), len(barcodes)):
        # Check if transposed
        if mat.shape == (len(barcodes), len(gene_ids)):
            mat = mat.T.tocsr()
        else:
            raise ValueError(
                f"Matrix shape {mat.shape} does not match "
                f"({len(gene_ids)} genes, {len(barcodes)} barcodes)."
            )

    single_cell_mat = SingleCellMatrix(
        gene_ids=tuple(gene_ids),
        cell_ids=tuple(barcodes),
        matrix=mat,
    )

    return single_cell_mat, obs_df, var_df


def detect_and_reject_rds(file_path: Path | str) -> None:
    """Inspect file for R serialization signature (.rds) and emit a clear rejection message."""
    p = Path(file_path)
    if p.suffix.lower() == ".rds":
        raise UnsupportedSingleCellFormatException(
            f"R serialized object detected at '{p}'. Direct ingestion of R serialized objects (.rds, "
            "Seurat, SingleCellExperiment) is unsupported. Please export the dataset to AnnData H5AD (.h5ad) "
            "or 10x Genomics Matrix Market format (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz) "
            "using SeuratDisk::SaveH5Seurat or Seurat::write10xCounts before ingestion."
        )

    # Check magic bytes for R serialization (gzip containing RDS or raw RDX)
    try:
        with p.open("rb") as f:
            header = f.read(8)
            if header.startswith(b"RDX2") or header.startswith(b"RDS2") or header.startswith(b"X\n"):
                raise UnsupportedSingleCellFormatException(
                    f"R serialization header detected in '{p}'. Please export to AnnData H5AD or 10x Matrix Market."
                )
            if header.startswith(b"\x1f\x8b"):
                # Gzipped; check first decompressed bytes
                with gzip.open(p, "rb") as gz:
                    inner_header = gz.read(8)
                    if inner_header.startswith(b"RDX2") or inner_header.startswith(b"RDS2") or inner_header.startswith(b"X\n"):
                        raise UnsupportedSingleCellFormatException(
                            f"Gzipped R serialization detected in '{p}'. Please export to AnnData H5AD or 10x Matrix Market."
                        )
    except (OSError, UnicodeError):
        pass
