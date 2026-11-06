"""Canonical data models for single-cell and single-nucleus RNA-seq ingestion, QC, and pseudobulk."""

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy import sparse

from src.qc.matrix import CountMatrix
from src.qc.result import Finding, Severity, Status


@dataclass(frozen=True, slots=True)
class CellMetadata:
    """Metadata and quality metrics for a single cell / nucleus barcode."""

    cell_id: str
    donor_id: str
    cell_type: str
    condition_id: str
    total_counts: int
    detected_genes: int
    mitochondrial_fraction: float
    is_pass_qc: bool
    qc_filter_reason: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class SingleCellMatrix:
    """Memory-efficient representation of a sparse or dense single-cell UMI count matrix."""

    def __init__(
        self,
        gene_ids: Sequence[str],
        cell_ids: Sequence[str],
        matrix: sparse.spmatrix | np.ndarray,
    ) -> None:
        self._gene_ids = tuple(str(g) for g in gene_ids)
        self._cell_ids = tuple(str(c) for c in cell_ids)
        self._gene_index = {g: i for i, g in enumerate(self._gene_ids)}
        self._cell_index = {c: i for i, c in enumerate(self._cell_ids)}

        if sparse.issparse(matrix):
            self._matrix = matrix.tocsr()
        else:
            self._matrix = np.asarray(matrix, dtype=np.int64)

        if self._matrix.shape != (len(self._gene_ids), len(self._cell_ids)):
            raise ValueError(
                f"Matrix shape {self._matrix.shape} does not match "
                f"({len(self._gene_ids)} genes, {len(self._cell_ids)} cells)."
            )

    @property
    def gene_ids(self) -> tuple[str, ...]:
        return self._gene_ids

    @property
    def cell_ids(self) -> tuple[str, ...]:
        return self._cell_ids

    @property
    def number_of_genes(self) -> int:
        return len(self._gene_ids)

    @property
    def number_of_cells(self) -> int:
        return len(self._cell_ids)

    @property
    def raw_matrix(self) -> sparse.csr_matrix | np.ndarray:
        return self._matrix

    def get_cell_column(self, cell_id: str) -> np.ndarray:
        idx = self._cell_index[cell_id]
        if sparse.issparse(self._matrix):
            return np.asarray(self._matrix[:, idx].toarray()).flatten()
        return self._matrix[:, idx]

    def get_gene_row(self, gene_id: str) -> np.ndarray:
        idx = self._gene_index[gene_id]
        if sparse.issparse(self._matrix):
            return np.asarray(self._matrix[idx, :].toarray()).flatten()
        return self._matrix[idx, :]

    def get_value(self, gene_id: str, cell_id: str) -> int:
        g_idx = self._gene_index[gene_id]
        c_idx = self._cell_index[cell_id]
        val = self._matrix[g_idx, c_idx]
        return int(val)

    def total_counts_per_cell(self) -> np.ndarray:
        if sparse.issparse(self._matrix):
            return np.asarray(self._matrix.sum(axis=0)).flatten().astype(np.int64)
        return np.asarray(self._matrix.sum(axis=0), dtype=np.int64)

    def detected_genes_per_cell(self) -> np.ndarray:
        if sparse.issparse(self._matrix):
            # Non-zero entries per column
            bin_mat = self._matrix.copy()
            bin_mat.data = np.ones_like(bin_mat.data)
            return np.asarray(bin_mat.sum(axis=0)).flatten().astype(np.int64)
        return np.asarray((self._matrix > 0).sum(axis=0), dtype=np.int64)

    def subset_cells(self, cell_ids: Sequence[str]) -> "SingleCellMatrix":
        indices = [self._cell_index[c] for c in cell_ids if c in self._cell_index]
        sub_cells = [self._cell_ids[i] for i in indices]
        sub_mat = self._matrix[:, indices]
        return SingleCellMatrix(self._gene_ids, sub_cells, sub_mat)


@dataclass(frozen=True, slots=True)
class SingleCellQCConfig:
    """Configurable quality thresholds for single-cell / single-nucleus filtering."""

    min_genes: int = 200
    min_counts: int = 500
    max_mito_fraction: float = 0.20
    mito_prefixes: tuple[str, ...] = ("MT-", "mt-", "Mt-")


@dataclass(frozen=True, slots=True)
class SingleCellQCResult:
    """Outcome of cell-level quality control evaluation."""

    total_cells: int
    passed_cells: int
    filtered_cells: int
    cell_metrics: tuple[CellMetadata, ...]
    findings: tuple[Finding, ...]
    status: Status

    def to_dict(self) -> dict[str, Any]:
        return {
            "total_cells": self.total_cells,
            "passed_cells": self.passed_cells,
            "filtered_cells": self.filtered_cells,
            "status": str(self.status),
            "findings": [asdict(f) for f in self.findings],
            "passed_cell_ids": [c.cell_id for c in self.cell_metrics if c.is_pass_qc],
        }


@dataclass(frozen=True, slots=True)
class PseudobulkConfig:
    """Parameters governing pseudobulk aggregation."""

    min_cells_per_sample: int = 20


@dataclass(frozen=True, slots=True)
class PseudobulkSample:
    """A biological sample pseudobulk aggregation unit for a specific cell type and donor."""

    sample_id: str
    donor_id: str
    cell_type: str
    condition_id: str
    cell_count: int
    total_umi_count: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class PseudobulkDataset:
    """Result of aggregating quality-controlled single-cell data into canonical pseudobulk."""

    dataset_id: str
    cell_types: tuple[str, ...]
    count_matrix: CountMatrix
    per_cell_type_matrices: dict[str, CountMatrix]
    pseudobulk_samples: tuple[PseudobulkSample, ...]
    metadata: dict[str, Any]
    provenance: dict[str, Any]
    findings: tuple[Finding, ...]
    status: Status
    emitted_matrix_path: Path | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "dataset_id": self.dataset_id,
            "cell_types": list(self.cell_types),
            "number_of_genes": len(self.count_matrix.gene_ids),
            "number_of_pseudobulk_samples": len(self.pseudobulk_samples),
            "pseudobulk_samples": [s.to_dict() for s in self.pseudobulk_samples],
            "status": str(self.status),
            "provenance": self.provenance,
            "emitted_matrix_path": str(self.emitted_matrix_path) if self.emitted_matrix_path else None,
            "findings": [asdict(f) for f in self.findings],
        }
