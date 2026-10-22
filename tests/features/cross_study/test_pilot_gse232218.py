"""Integration verification tests for GSE232218 VPA 200 uM Stage 8A.1 Pilot."""

import json
from pathlib import Path
import pandas as pd
import pytest

from src.features.cross_study.serialization import CANONICAL_LONG_COLUMNS


@pytest.fixture
def pilot_artifacts() -> tuple[Path, Path]:
    csv_path = Path("data/processed/GSE232218_VPA_200uM_condition_gene_response_registry.csv")
    qc_path = Path("data/processed/GSE232218_VPA_200uM_stage8a_qc.json")
    if not csv_path.exists() or not qc_path.exists():
        pytest.skip("Pilot artifacts not present. Run scripts/run_vpa_pilot_stage8a.py first.")
    return csv_path, qc_path


def test_pilot_qc_json(pilot_artifacts: tuple[Path, Path]) -> None:
    _, qc_path = pilot_artifacts
    with open(qc_path, encoding="utf-8") as f:
        qc_data = json.load(f)

    assert qc_data["is_valid"] is True
    assert qc_data["total_records_evaluated"] == 60207
    assert qc_data["numerical_mismatches_count"] == 0
    assert qc_data["gene_identity_mismatches_count"] == 0
    assert qc_data["direction_mismatches_count"] == 0
    assert qc_data["context_mismatches_count"] == 0
    assert qc_data["leakage_violations_count"] == 0

    acc = qc_data["reconciled_accounting"]
    assert acc["stage7_canonical_evaluated_genes"] == 60207
    assert acc["stage8a_records_produced"] == 60207
    assert acc["valid_estimates_ok"] == 32140
    assert acc["zero_counts"] == 27980
    assert acc["cooks_outliers"] == 87
    assert acc["single_study_only_count"] == 60207


def test_pilot_csv_schema_and_contents(pilot_artifacts: tuple[Path, Path]) -> None:
    csv_path, _ = pilot_artifacts
    df = pd.read_csv(csv_path)

    # 1. Total row count matches 60,207 evaluated genes
    assert len(df) == 60207

    # 2. Columns match canonical governed long schema
    assert list(df.columns) == CANONICAL_LONG_COLUMNS

    # 3. Independent study identity is GSE232218 across all rows
    assert (df["independent_study_id"] == "GSE232218").all()

    # 4. Cross-study status is SINGLE_STUDY_ONLY across all rows
    assert (df["cross_study_status"] == "SINGLE_STUDY_ONLY").all()

    # 5. Compound identity is CID:3121 / Valproic Acid across all rows
    assert (df["compound_id"] == "CID:3121").all()
    assert (df["compound_name"] == "Valproic Acid").all()

    # 6. Dose is 200.0 uM across all rows
    assert (df["dose_value"] == 200.0).all()
    assert (df["dose_unit"] == "uM").all()

    # 7. Exact count of OK estimates matches 32,140
    assert (df["estimate_status"] == "OK").sum() == 32140
    assert (df["estimate_status"] == "ZERO_COUNTS").sum() == 27980
    assert (df["estimate_status"] == "COOKS_OUTLIER").sum() == 87


def test_numerical_identity_with_stage8c_universe(pilot_artifacts: tuple[Path, Path]) -> None:
    csv_path, _ = pilot_artifacts
    stage8c_univ_path = Path("data/processed/GSE232218_VPA_200uM_statistically_eligible_universe.csv")
    if not stage8c_univ_path.exists():
        pytest.skip("Stage 8C universe file not found.")

    df_registry = pd.read_csv(csv_path)
    df_univ = pd.read_csv(stage8c_univ_path)

    # 32,140 statistically eligible genes in Stage 8C universe
    assert len(df_univ) == 32140

    # Merge on canonical_gene_id
    merged = pd.merge(
        df_univ,
        df_registry,
        left_on="canonical_gene_id",
        right_on="canonical_ensembl_gene_id",
        suffixes=("_univ", "_reg"),
    )
    assert len(merged) == 32140

    # Verify numerical preservation across all 32,140 genes
    diff_lfc = (merged["log2_fold_change_univ"] - merged["log2_fold_change_reg"]).abs()
    assert (diff_lfc < 1e-5).all(), f"Max log2FC diff: {diff_lfc.max()}"

    diff_se = (merged["lfc_standard_error"] - merged["log2_fold_change_se"]).abs()
    assert (diff_se < 1e-5).all(), f"Max SE diff: {diff_se.max()}"

    diff_wald = (merged["wald_statistic_univ"] - merged["wald_statistic_reg"]).abs()
    assert (diff_wald < 1e-5).all(), f"Max Wald diff: {diff_wald.max()}"
