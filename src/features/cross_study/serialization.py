"""Deterministic serialization services for Stage 8A.1 Condition-Level Gene Response.

Exports registry records to long-format pandas DataFrame and CSV files
with stable column ordering and deterministic row sorting.
"""

import json
from pathlib import Path
from typing import Any
import pandas as pd

from src.features.cross_study.models import ConditionGeneResponseRegistry


CANONICAL_LONG_COLUMNS = [
    "study_id",
    "study_family_id",
    "experiment_id",
    "contrast_id",
    "compound_id",
    "compound_name",
    "dose_value",
    "dose_unit",
    "dose_original",
    "developmental_stage",
    "exposure_duration",
    "model_type",
    "canonical_ensembl_gene_id",
    "gene_symbol",
    "base_mean",
    "log2_fold_change",
    "log2_fold_change_se",
    "wald_statistic",
    "p_value",
    "adjusted_p_value",
    "response_direction",
    "estimate_status",
    "missing_reason",
    "independent_study_id",
    "cross_study_status",
    "provenance",
]


def export_registry_to_dataframe(registry: ConditionGeneResponseRegistry) -> pd.DataFrame:
    """Convert ConditionGeneResponseRegistry to a canonical long-format DataFrame."""
    rows: list[dict[str, Any]] = []
    for r in registry.records:
        prov_summary = json.dumps(
            {
                "source_stage": r.source_stage,
                "source_contrast_id": r.source_contrast_id,
                "de_version": r.de_version,
                "registry_version": r.registry_version,
            },
            sort_keys=True,
        )

        row = {
            "study_id": r.study_id,
            "study_family_id": r.study_family_id,
            "experiment_id": r.experiment_id,
            "contrast_id": r.contrast_id,
            "compound_id": r.compound_id,
            "compound_name": r.compound_name,
            "dose_value": r.dose_value,
            "dose_unit": r.dose_unit,
            "dose_original": r.dose_original,
            "developmental_stage": r.developmental_stage,
            "exposure_duration": r.exposure_duration,
            "model_type": r.model_type,
            "canonical_ensembl_gene_id": r.canonical_ensembl_gene_id,
            "gene_symbol": r.gene_symbol,
            "base_mean": r.base_mean,
            "log2_fold_change": r.log2_fold_change,
            "log2_fold_change_se": r.log2_fold_change_se,
            "wald_statistic": r.wald_statistic,
            "p_value": r.p_value,
            "adjusted_p_value": r.adjusted_p_value,
            "response_direction": str(r.response_direction),
            "estimate_status": str(r.estimate_status),
            "missing_reason": r.missing_reason,
            "independent_study_id": r.independent_study_id,
            "cross_study_status": str(r.cross_study_status),
            "provenance": prov_summary,
        }
        rows.append(row)

    df = pd.DataFrame(rows, columns=CANONICAL_LONG_COLUMNS)
    if not df.empty:
        df = df.sort_values(
            by=["study_id", "contrast_id", "canonical_ensembl_gene_id"],
            ascending=[True, True, True],
        ).reset_index(drop=True)
    return df


def export_registry_to_csv(registry: ConditionGeneResponseRegistry, output_path: str | Path) -> Path:
    """Export the registry to a deterministic CSV file."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df = export_registry_to_dataframe(registry)
    df.to_csv(path, index=False)
    return path
