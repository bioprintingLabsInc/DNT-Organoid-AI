"""Single-cell and single-nucleus RNA-seq ingestion, quality control, and pseudobulk module."""

from src.single_cell.models import (
    CellMetadata,
    InvalidRawCountMatrixException,
    PseudobulkConfig,
    PseudobulkDataset,
    PseudobulkSample,
    SingleCellMatrix,
    SingleCellQCConfig,
    SingleCellQCResult,
)
from src.single_cell.pipeline import run_single_cell_ingestion_pipeline
from src.single_cell.pseudobulk import aggregate_pseudobulk
from src.single_cell.qc import evaluate_single_cell_qc
from src.single_cell.reader import (
    UnsupportedSingleCellFormatException,
    detect_and_reject_rds,
    read_10x_matrix_market,
    read_h5ad,
)

__all__ = [
    "CellMetadata",
    "SingleCellMatrix",
    "SingleCellQCConfig",
    "SingleCellQCResult",
    "PseudobulkConfig",
    "PseudobulkSample",
    "PseudobulkDataset",
    "read_h5ad",
    "read_10x_matrix_market",
    "detect_and_reject_rds",
    "UnsupportedSingleCellFormatException",
    "InvalidRawCountMatrixException",
    "evaluate_single_cell_qc",
    "aggregate_pseudobulk",
    "run_single_cell_ingestion_pipeline",
]
