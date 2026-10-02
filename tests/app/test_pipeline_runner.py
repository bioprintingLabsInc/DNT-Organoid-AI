"""Unit and integration tests for the Streamlit application pipeline runner."""

import io
import json
import pytest
import pandas as pd

from src.app.pipeline_runner import (
    ChemicalParams,
    PipelineRunResult,
    build_canonical_metadata,
    generate_demo_data,
    parse_count_matrix,
    parse_sample_info,
    run_full_analysis_pipeline,
)
from src.response_builder.models import ContrastStatus


def test_parse_count_matrix_valid_df() -> None:
    df_raw = pd.DataFrame(
        {
            "gene": ["ENSG00000100292", "ENSG00000128272"],
            "sample_1": [100, 200],
            "sample_2": [110, 210],
        }
    )
    df, sids = parse_count_matrix(df_raw)
    assert sids == ["sample_1", "sample_2"]
    assert "gene_id" in df.columns
    assert len(df) == 2
    assert df["sample_1"].iloc[0] == 100


def test_parse_count_matrix_csv_and_tsv_strings() -> None:
    csv_str = "gene_id,s1,s2\nENSG00000100292,50,60\nENSG00000128272,70,80\n"
    df_csv, sids_csv = parse_count_matrix(csv_str)
    assert sids_csv == ["s1", "s2"]
    assert len(df_csv) == 2

    tsv_str = "gene_id\ts1\ts2\nENSG00000100292\t50\t60\nENSG00000128272\t70\t80\n"
    df_tsv, sids_tsv = parse_count_matrix(tsv_str)
    assert sids_tsv == ["s1", "s2"]
    assert len(df_tsv) == 2


def test_parse_count_matrix_invalid() -> None:
    # Less than 2 columns
    with pytest.raises(ValueError, match="at least 2 sample columns"):
        parse_count_matrix("gene_id\nENSG00000100292\n")

    # Non-numeric counts
    with pytest.raises(ValueError, match="non-numeric"):
        parse_count_matrix("gene_id,s1,s2\nENSG00000100292,abc,100\n")

    # Negative counts
    with pytest.raises(ValueError, match="negative"):
        parse_count_matrix("gene_id,s1,s2\nENSG00000100292,-5,100\n")


def test_parse_sample_info_explicit() -> None:
    sample_df_raw = pd.DataFrame(
        {
            "sample_id": ["ctrl_1", "trt_1"],
            "condition_type": ["control", "treatment"],
            "biological_replicate_id": ["b_c1", "b_t1"],
        }
    )
    res_df = parse_sample_info(sample_df_raw, ["ctrl_1", "trt_1"])
    assert list(res_df["condition_type"]) == ["control", "treatment"]
    assert list(res_df["biological_replicate_id"]) == ["b_c1", "b_t1"]


def test_parse_sample_info_auto_inference() -> None:
    sids = ["dmso_rep1", "dmso_rep2", "rotenone_rep1", "rotenone_rep2"]
    res_df = parse_sample_info(None, sids)
    assert len(res_df) == 4
    assert res_df.loc[res_df["sample_id"] == "dmso_rep1", "condition_type"].iloc[0] == "control"
    assert res_df.loc[res_df["sample_id"] == "rotenone_rep1", "condition_type"].iloc[0] == "treatment"


def test_build_canonical_metadata() -> None:
    samples_df = pd.DataFrame(
        {
            "sample_id": ["c1", "c2", "t1", "t2"],
            "condition_type": ["control", "control", "treatment", "treatment"],
            "biological_replicate_id": ["b_c1", "b_c2", "b_t1", "b_t2"],
        }
    )
    params = ChemicalParams(agent_name="Rotenone", concentration_or_dose=0.5, developmental_age="day_35")
    meta = build_canonical_metadata(samples_df, params)

    assert "studies" in meta
    assert "experiments" in meta
    assert "conditions" in meta
    assert "samples" in meta
    assert len(meta["samples"]) == 4
    assert meta["conditions"][0]["developmental_age_or_stage"] == "day_35"
    assert meta["exposures"][0]["agent_name"] == "Rotenone"
    assert meta["treatment_control_relationships"][0]["relationship_id"] == "rel_app_01"


def test_generate_demo_data() -> None:
    counts_df, samples_df, params = generate_demo_data()
    assert len(counts_df) == 40
    assert len(samples_df) == 6
    assert sum(samples_df["condition_type"] == "control") == 3
    assert sum(samples_df["condition_type"] == "treatment") == 3
    assert params.agent_name == "Rotenone"


def test_run_pipeline_insufficient_replicates() -> None:
    counts_df = pd.DataFrame(
        {
            "gene_id": ["ENSG00000100292"],
            "ctrl_1": [100],
            "trt_1": [200],
        }
    )
    samples_df = pd.DataFrame(
        {
            "sample_id": ["ctrl_1", "trt_1"],
            "condition_type": ["control", "treatment"],
            "biological_replicate_id": ["rep_1", "rep_2"],
        }
    )
    params = ChemicalParams()
    res = run_full_analysis_pipeline(counts_df, samples_df, params)
    assert not res.is_success
    assert "at least 2 biological replicates" in res.error_message


def test_run_pipeline_demo_data_live() -> None:
    """Live full integration test running all 5 backend stages with real DESeq2 in R."""
    counts_df, samples_df, params = generate_demo_data()
    res = run_full_analysis_pipeline(counts_df, samples_df, params)

    assert res.is_success, f"Pipeline failed: {res.error_message}"
    assert res.contrast_result is not None
    assert res.contrast_result.status == ContrastStatus.ELIGIBLE
    assert res.contrast_result.is_inferentially_eligible

    # Gene conservation
    assert res.summary_metrics["is_gene_universe_conserved"] is True
    assert res.summary_metrics["total_input_genes"] == 40
    assert res.summary_metrics["evaluated_canonical_genes"] == 40
    assert res.summary_metrics["excluded_genes_count"] == 0

    # DE statistical outputs
    assert len(res.de_genes_df) == 40
    assert res.summary_metrics["significant_upregulated"] >= 5
    assert res.summary_metrics["significant_downregulated"] >= 5

    # Check key columns
    expected_cols = [
        "canonical_gene_id",
        "approved_symbol",
        "base_mean",
        "log2_fold_change",
        "lfc_standard_error",
        "wald_statistic",
        "p_value",
        "adjusted_p_value_bh",
        "significance",
        "direction",
    ]
    for col in expected_cols:
        assert col in res.de_genes_df.columns

    # JSON serialization
    summary_dict = res.to_summary_dict()
    json_str = json.dumps(summary_dict)
    assert len(json_str) > 0
