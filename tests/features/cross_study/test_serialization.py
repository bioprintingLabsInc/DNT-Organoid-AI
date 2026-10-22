"""Unit tests for Stage 8A.1 serialization services."""

import tempfile
from pathlib import Path
import pandas as pd
import pytest

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
)
from src.features.cross_study.registry import (
    build_condition_gene_response_registry,
)
from src.features.cross_study.serialization import (
    CANONICAL_LONG_COLUMNS,
    export_registry_to_csv,
    export_registry_to_dataframe,
)


def test_dataframe_columns_match_governed_schema(
    synthetic_contrast_vpa_200: DifferentialExpressionContrastResult,
) -> None:
    registry = build_condition_gene_response_registry([synthetic_contrast_vpa_200])
    df = export_registry_to_dataframe(registry)

    # Verify all governed columns exist in exact canonical sequence
    assert list(df.columns) == CANONICAL_LONG_COLUMNS
    assert len(df) == 8

    # Check key columns
    assert "canonical_ensembl_gene_id" in df.columns
    assert "log2_fold_change" in df.columns
    assert "log2_fold_change_se" in df.columns
    assert "wald_statistic" in df.columns
    assert "p_value" in df.columns
    assert "adjusted_p_value" in df.columns
    assert "response_direction" in df.columns
    assert "estimate_status" in df.columns
    assert "independent_study_id" in df.columns


def test_deterministic_csv_export(
    synthetic_contrast_vpa_200: DifferentialExpressionContrastResult,
) -> None:
    registry = build_condition_gene_response_registry([synthetic_contrast_vpa_200])

    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path1 = Path(tmpdir) / "registry_1.csv"
        csv_path2 = Path(tmpdir) / "registry_2.csv"

        export_registry_to_csv(registry, csv_path1)
        export_registry_to_csv(registry, csv_path2)

        # 27. Output serialization is deterministic: identical file bytes
        content1 = csv_path1.read_bytes()
        content2 = csv_path2.read_bytes()
        assert content1 == content2

        # Verify CSV readability and content
        df_read = pd.read_csv(csv_path1)
        assert len(df_read) == 8
        assert list(df_read.columns) == CANONICAL_LONG_COLUMNS
        assert df_read["compound_id"].iloc[0] == "CID:3121"
        assert df_read["independent_study_id"].iloc[0] == "GSE232218"
