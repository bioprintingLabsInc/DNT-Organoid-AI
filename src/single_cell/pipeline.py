"""Single-cell and single-nucleus RNA-seq end-to-end ingestion and pseudobulk pipeline."""

from pathlib import Path
from typing import Any

from src.qc.result import Finding, Severity, Status
from src.single_cell.models import (
    PseudobulkConfig,
    PseudobulkDataset,
    SingleCellQCConfig,
    SingleCellQCResult,
)
from src.single_cell.pseudobulk import aggregate_pseudobulk
from src.single_cell.qc import evaluate_single_cell_qc
from src.single_cell.reader import (
    UnsupportedSingleCellFormatException,
    detect_and_reject_rds,
    read_10x_matrix_market,
    read_h5ad,
)


def run_single_cell_ingestion_pipeline(
    asset_path: Path | str,
    metadata: dict[str, Any] | None = None,
    qc_config: SingleCellQCConfig | None = None,
    pseudobulk_config: PseudobulkConfig | None = None,
    output_tsv_path: Path | str | None = None,
    dataset_id: str = "single_cell_dataset",
) -> PseudobulkDataset:
    """Execute end-to-end scRNA-seq / snRNA-seq ingestion:

    1. Read AnnData H5AD or 10x Matrix Market (rejecting .rds R serialization).
    2. Compute cell-level QC (total UMIs, detected genes, mitochondrial fraction).
    3. Filter low-quality cells against QC thresholds.
    4. Group cells by (donor, cell_type, condition) and strictly sum raw integer UMI counts.
    5. Enforce donor biological independence (no pseudo-replicate inflation).
    6. Emit canonical bulk raw count matrix and aligned metadata for downstream DESeq2.
    """
    path = Path(asset_path)
    meta = metadata or {}

    # 1. Format check & RDS rejection
    detect_and_reject_rds(path)

    # 2. Intake: H5AD or 10x Matrix Market
    if path.is_file() and path.suffix.lower() == ".h5ad":
        matrix, obs_df, var_df = read_h5ad(path)
    elif path.is_dir() or path.name.lower().startswith("matrix.mtx"):
        matrix, obs_df, var_df = read_10x_matrix_market(path)
    else:
        # Fallback inspection
        if path.is_file():
            # Check if HDF5
            with path.open("rb") as f:
                magic = f.read(8)
                if magic.startswith(b"\x89HDF\r\n\x1a\n"):
                    matrix, obs_df, var_df = read_h5ad(path)
                else:
                    raise UnsupportedSingleCellFormatException(
                        f"Unrecognized single-cell file format for '{path}'. Supported: .h5ad, 10x Matrix Market directory."
                    )
        else:
            raise UnsupportedSingleCellFormatException(
                f"Single-cell asset path '{path}' does not exist or is unrecognized."
            )

    # 3. Cell QC
    qc_res = evaluate_single_cell_qc(
        matrix=matrix,
        obs_df=obs_df,
        var_df=var_df,
        config=qc_config,
        metadata=meta,
    )

    if qc_res.status == Status.FAIL:
        from src.qc.matrix import CountMatrix
        empty_mat = CountMatrix(matrix.gene_ids, (), ())
        return PseudobulkDataset(
            dataset_id=dataset_id,
            cell_types=(),
            count_matrix=empty_mat,
            per_cell_type_matrices={},
            pseudobulk_samples=(),
            metadata={},
            provenance={"status": "FAILED_QC"},
            findings=qc_res.findings,
            status=Status.FAIL,
        )

    # 4. Pseudobulk Aggregation
    pb_dataset = aggregate_pseudobulk(
        matrix=matrix,
        qc_result=qc_res,
        config=pseudobulk_config,
        output_tsv_path=output_tsv_path,
        dataset_id=dataset_id,
        metadata_template=meta,
    )

    return pb_dataset
