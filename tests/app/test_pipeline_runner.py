"""Unit and integration tests for the Streamlit application and pipeline runner.

Verifies:
1. Synthetic demo mode is no longer exposed in the application.
2. Treatment/control status and biological replicates cannot be inferred from filenames or sample names.
3. Missing or invalid explicit sample metadata fails safely with clear validation errors.
4. Valid explicit user-provided sample metadata executes through the full locked backend.
"""

from pathlib import Path
import random
import pandas as pd
import pytest

from src.app.pipeline_runner import (
    ChemicalParams,
    PipelineRunResult,
    build_canonical_metadata,
    parse_count_matrix,
    parse_sample_info,
    run_full_analysis_pipeline,
)
from src.response_builder.models import ContrastStatus


def _create_test_dataset(n_genes: int = 40) -> tuple[pd.DataFrame, pd.DataFrame, ChemicalParams]:
    """Helper to assemble a test expression matrix and explicit metadata for testing."""
    sids = ["donor1_ctrl_rep1", "donor1_ctrl_rep2", "donor1_ctrl_rep3", "donor1_trt_rep1", "donor1_trt_rep2", "donor1_trt_rep3"]

    genes = [
        "ENSG00000100292", "ENSG00000128272", "ENSG00000175197", "ENSG00000001084",
        "ENSG00000181019", "ENSG00000044574", "ENSG00000116717", "ENSG00000087088",
        "ENSG00000164305", "ENSG00000124762", "ENSG00000161011", "ENSG00000291237",
        "ENSG00000090266", "ENSG00000167792", "ENSG00000078018", "ENSG00000258947",
        "ENSG00000008056", "ENSG00000132639", "ENSG00000167281", "ENSG00000102003",
        "ENSG00000106089", "ENSG00000176884", "ENSG00000155511", "ENSG00000070808",
        "ENSG00000001617", "ENSG00000001631", "ENSG00000111640", "ENSG00000075624",
        "ENSG00000166710", "ENSG00000089157", "ENSG00000150991", "ENSG00000196262",
        "ENSG00000112592", "ENSG00000168488", "ENSG00000102144", "ENSG00000164924",
        "ENSG00000073578", "ENSG00000144713", "ENSG00000000003", "ENSG00000000419",
    ][:n_genes]

    mean_levels = [40, 80, 150, 300, 600, 1200, 2500, 5000]
    rng = random.Random(2026)
    rows = []

    for i, gid in enumerate(genes):
        base_mu = mean_levels[i % len(mean_levels)]
        if i < 12:
            fc = 3.5
        elif i < 26:
            fc = 0.28
        else:
            fc = 1.0

        row = [gid]
        for sid in sids:
            is_trt = "trt" in sid
            mu = base_mu * (fc if is_trt else 1.0)
            val = max(10, int(round(rng.gauss(mu, max(3.0, mu * 0.12)))))
            row.append(val)
        rows.append(row)

    counts_df = pd.DataFrame(rows, columns=["gene_id"] + sids)

    # Explicit sample metadata
    samples_records = []
    for sid in sids:
        ctype = "treatment" if "trt" in sid else "control"
        rep_id = f"biorep_{sid}"
        samples_records.append(
            {
                "sample_id": sid,
                "condition_type": ctype,
                "biological_replicate_id": rep_id,
            }
        )
    samples_df = pd.DataFrame(samples_records)

    params = ChemicalParams(
        agent_name="Bisphenol_A",
        agent_identifier="CID:6623",
        concentration_or_dose=10.0,
        concentration_or_dose_unit="uM",
        developmental_age="day_35",
        exposure_duration=24.0,
        exposure_duration_unit="h",
        vehicle="0.1% DMSO",
        organoid_context_id="ctx_cerebral_organoid",
        organoid_type="cerebral_organoid",
        brain_region="cerebral_cortex",
        biological_source_id="ipsc_line_donor_1",
        species="Homo sapiens",
    )

    return counts_df, samples_df, params


# -----------------------------------------------------------------------------
# 1. Verification that Demo Mode is Not Exposed
# -----------------------------------------------------------------------------
def test_demo_mode_not_exposed_in_app() -> None:
    """Verify that app.py does not expose synthetic demo mode or generate_demo_data."""
    app_path = Path(__file__).resolve().parents[2] / "app.py"
    assert app_path.exists()
    content = app_path.read_text(encoding="utf-8")

    assert "generate_demo_data" not in content, "generate_demo_data found in app.py"
    assert "Load Demo Data" not in content, "Demo Data option found in app.py"
    assert "Rotenone" not in content, "Rotenone hardcoded example found in app.py"


def test_generate_demo_data_removed_from_pipeline_runner() -> None:
    """Verify that pipeline_runner module no longer exports generate_demo_data."""
    import src.app.pipeline_runner as pr
    assert not hasattr(pr, "generate_demo_data"), "generate_demo_data must be removed from pipeline_runner"


# -----------------------------------------------------------------------------
# 2. Verification that Automatic Inference is Prohibited
# -----------------------------------------------------------------------------
def test_treatment_control_cannot_be_inferred_from_sample_names() -> None:
    """Attempting to parse sample info without explicit metadata must fail."""
    sids = ["ctrl_sample_1", "ctrl_sample_2", "treatment_sample_1", "treatment_sample_2"]
    with pytest.raises(ValueError, match="Explicit sample metadata is required"):
        parse_sample_info(None, sids)


def test_missing_condition_column_fails_safely() -> None:
    """Sample metadata without an explicit condition column must raise ValueError."""
    sids = ["sample_1", "sample_2"]
    invalid_df = pd.DataFrame(
        {
            "sample_id": sids,
            "biological_replicate_id": ["rep_1", "rep_2"],
        }
    )
    with pytest.raises(ValueError, match="explicit condition column"):
        parse_sample_info(invalid_df, sids)


def test_ambiguous_condition_values_fail_safely() -> None:
    """Sample metadata with ambiguous condition values must raise ValueError."""
    sids = ["s1", "s2"]
    ambiguous_df = pd.DataFrame(
        {
            "sample_id": sids,
            "condition_type": ["treated_maybe", "unknown_group"],
            "biological_replicate_id": ["rep_1", "rep_2"],
        }
    )
    with pytest.raises(ValueError, match="Invalid condition value"):
        parse_sample_info(ambiguous_df, sids)


def test_missing_biological_replicate_column_fails_safely() -> None:
    """Sample metadata without an explicit biological replicate column must raise ValueError."""
    sids = ["s1", "s2"]
    no_rep_df = pd.DataFrame(
        {
            "sample_id": sids,
            "condition_type": ["control", "treatment"],
        }
    )
    with pytest.raises(ValueError, match="explicit biological replicate column"):
        parse_sample_info(no_rep_df, sids)


def test_empty_biological_replicate_value_fails_safely() -> None:
    """Sample metadata with empty biological replicate identifier must raise ValueError."""
    sids = ["s1", "s2"]
    empty_rep_df = pd.DataFrame(
        {
            "sample_id": sids,
            "condition_type": ["control", "treatment"],
            "biological_replicate_id": ["rep_1", ""],
        }
    )
    with pytest.raises(ValueError, match="Missing biological replicate identifier"):
        parse_sample_info(empty_rep_df, sids)


def test_missing_matrix_sample_in_metadata_fails_safely() -> None:
    """Count matrix sample missing from metadata must raise ValueError."""
    matrix_sids = ["s1", "s2", "s3"]
    partial_df = pd.DataFrame(
        {
            "sample_id": ["s1", "s2"],
            "condition_type": ["control", "treatment"],
            "biological_replicate_id": ["rep_1", "rep_2"],
        }
    )
    with pytest.raises(ValueError, match="missing explicit entries for matrix sample"):
        parse_sample_info(partial_df, matrix_sids)


# -----------------------------------------------------------------------------
# 3. Pipeline Gating & Safe Failure
# -----------------------------------------------------------------------------
def test_pipeline_rejects_missing_sample_metadata() -> None:
    """Pipeline runner must return is_success=False when samples_input is None."""
    counts_df, _, params = _create_test_dataset()
    res = run_full_analysis_pipeline(counts_df, None, params)
    assert not res.is_success
    assert "Explicit sample metadata is required" in (res.error_message or "")


def test_pipeline_rejects_missing_chemical_name() -> None:
    """Pipeline runner must reject empty chemical agent name."""
    counts_df, samples_df, params = _create_test_dataset()
    empty_agent_params = ChemicalParams(agent_name="", developmental_age="day_35")
    res = run_full_analysis_pipeline(counts_df, samples_df, empty_agent_params)
    assert not res.is_success
    assert "Chemical / Agent Name is required" in (res.error_message or "")


def test_pipeline_rejects_missing_developmental_age() -> None:
    """Pipeline runner must reject empty developmental age."""
    counts_df, samples_df, params = _create_test_dataset()
    empty_age_params = ChemicalParams(agent_name="Compound_X", developmental_age="")
    res = run_full_analysis_pipeline(counts_df, samples_df, empty_age_params)
    assert not res.is_success
    assert "Developmental Age / Stage is required" in (res.error_message or "")


def test_pipeline_rejects_under_replicated_cohort() -> None:
    """Pipeline runner must safely block cohorts with < 2 biological replicates."""
    counts_df, samples_df, params = _create_test_dataset()
    # Filter to 1 treatment and 1 control
    under_replicated_samples = samples_df.iloc[[0, 3]].copy()
    sub_counts = counts_df[["gene_id"] + list(under_replicated_samples["sample_id"])]

    res = run_full_analysis_pipeline(sub_counts, under_replicated_samples, params)
    assert not res.is_success
    assert "at least 2 biological replicates" in (res.error_message or "")


# -----------------------------------------------------------------------------
# 4. End-to-End Pipeline Execution with Valid Explicit Metadata
# -----------------------------------------------------------------------------
def test_pipeline_runs_with_valid_explicit_metadata() -> None:
    """Live integration test running all 5 backend stages with real DESeq2 using explicit metadata."""
    counts_df, samples_df, params = _create_test_dataset(n_genes=40)
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
