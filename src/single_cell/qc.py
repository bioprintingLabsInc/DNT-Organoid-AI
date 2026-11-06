"""Quality control metrics computation and cell filtering for single-cell/single-nucleus RNA-seq."""

from typing import Any, Sequence
import numpy as np
import pandas as pd

from src.qc.result import Finding, Severity, Status
from src.single_cell.models import (
    CellMetadata,
    SingleCellMatrix,
    SingleCellQCConfig,
    SingleCellQCResult,
)


def _is_mitochondrial_gene(gene_id: str, gene_symbol: str | None, prefixes: tuple[str, ...]) -> bool:
    """Check whether a gene is mitochondrial based on prefix or Ensembl mitochondrial ID."""
    g_id = gene_id.strip()
    sym = (gene_symbol or "").strip()

    for p in prefixes:
        if g_id.startswith(p) or sym.startswith(p):
            return True

    # Common human mitochondrial gene Ensembl prefix or known MT genes
    if g_id.upper().startswith("MT-") or sym.upper().startswith("MT-"):
        return True

    return False


def evaluate_single_cell_qc(
    matrix: SingleCellMatrix,
    obs_df: pd.DataFrame,
    var_df: pd.DataFrame | None = None,
    config: SingleCellQCConfig | None = None,
    metadata: dict[str, Any] | None = None,
) -> SingleCellQCResult:
    """Compute cell-level quality control metrics and filter low-quality barcodes.

    Metrics:
    - Total UMI counts per cell
    - Number of detected genes (count > 0) per cell
    - Mitochondrial read fraction per cell
    """
    if config is None:
        config = SingleCellQCConfig()

    findings: list[Finding] = []

    # 1. Identify mitochondrial genes
    mito_indices: list[int] = []
    gene_symbols = {}
    if var_df is not None and "gene_symbol" in var_df.columns:
        gene_symbols = var_df["gene_symbol"].to_dict()

    for i, g_id in enumerate(matrix.gene_ids):
        sym = gene_symbols.get(g_id)
        if _is_mitochondrial_gene(g_id, sym, config.mito_prefixes):
            mito_indices.append(i)

    # 2. Vectorized metric computation
    total_counts = matrix.total_counts_per_cell()
    detected_genes = matrix.detected_genes_per_cell()

    if mito_indices:
        # Sum rows corresponding to mito genes
        if hasattr(matrix.raw_matrix, "tocsr"):
            sub_mat = matrix.raw_matrix[mito_indices, :]
            mito_counts = np.asarray(sub_mat.sum(axis=0)).flatten().astype(np.int64)
        else:
            mito_counts = np.asarray(matrix.raw_matrix[mito_indices, :].sum(axis=0), dtype=np.int64)
        mito_fraction = np.where(total_counts > 0, mito_counts / total_counts, 0.0)
    else:
        mito_counts = np.zeros(matrix.number_of_cells, dtype=np.int64)
        mito_fraction = np.zeros(matrix.number_of_cells, dtype=np.float64)

    # 3. Resolve cell metadata columns
    # Find donor / sample ID column
    donor_col = next((c for c in ("donor_id", "donor", "sample_id", "sample", "patient", "orig.ident") if c in obs_df.columns), None)
    # Find cell type / cluster column
    ct_col = next((c for c in ("cell_type", "celltype", "cluster", "cell_annotation", "seurat_clusters", "predicted_celltype") if c in obs_df.columns), None)
    # Find condition / treatment column
    cond_col = next((c for c in ("condition_id", "condition", "treatment", "group", "status") if c in obs_df.columns), None)

    cell_metrics: list[CellMetadata] = []
    passed_count = 0
    filtered_count = 0

    for idx, cell_id in enumerate(matrix.cell_ids):
        t_counts = int(total_counts[idx])
        d_genes = int(detected_genes[idx])
        m_frac = float(mito_fraction[idx])

        # Cell metadata attributes
        donor = str(obs_df.iloc[idx][donor_col]) if donor_col else "donor_1"
        cell_type = str(obs_df.iloc[idx][ct_col]) if ct_col else "unclassified"
        condition = str(obs_df.iloc[idx][cond_col]) if cond_col else "control"

        # Quality filtering
        reasons: list[str] = []
        if d_genes < config.min_genes:
            reasons.append(f"low_detected_genes ({d_genes} < {config.min_genes})")
        if t_counts < config.min_counts:
            reasons.append(f"low_total_counts ({t_counts} < {config.min_counts})")
        if m_frac > config.max_mito_fraction:
            reasons.append(f"high_mitochondrial_fraction ({m_frac:.2%} > {config.max_mito_fraction:.2%})")

        is_pass = len(reasons) == 0
        if is_pass:
            passed_count += 1
            reason_str = None
        else:
            filtered_count += 1
            reason_str = "; ".join(reasons)

        cell_metrics.append(
            CellMetadata(
                cell_id=cell_id,
                donor_id=donor,
                cell_type=cell_type,
                condition_id=condition,
                total_counts=t_counts,
                detected_genes=d_genes,
                mitochondrial_fraction=m_frac,
                is_pass_qc=is_pass,
                qc_filter_reason=reason_str,
            )
        )

    # 4. Global QC findings
    if passed_count == 0:
        findings.append(
            Finding(
                Severity.ERROR,
                "all_cells_filtered",
                "cell_qc",
                None,
                "cell_filtering",
                f"All {matrix.number_of_cells} cells failed single-cell quality thresholds (min_genes={config.min_genes}, min_counts={config.min_counts}, max_mito={config.max_mito_fraction:.0%}).",
            )
        )
        status = Status.FAIL
    elif filtered_count > 0.5 * matrix.number_of_cells:
        findings.append(
            Finding(
                Severity.REVIEW,
                "high_cell_filtration_rate",
                "cell_qc",
                None,
                "cell_filtering",
                f"{filtered_count}/{matrix.number_of_cells} cells ({filtered_count/matrix.number_of_cells:.1%}) filtered out during QC.",
            )
        )
        status = Status.NEEDS_REVIEW
    else:
        status = Status.PASS

    return SingleCellQCResult(
        total_cells=matrix.number_of_cells,
        passed_cells=passed_count,
        filtered_cells=filtered_count,
        cell_metrics=tuple(cell_metrics),
        findings=tuple(findings),
        status=status,
    )
