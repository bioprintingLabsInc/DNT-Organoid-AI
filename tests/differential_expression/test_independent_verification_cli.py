"""Unit and regression tests for independent DESeq2 verification.

Validates that:
1. Fabricated biological replicates or missing replicate metadata are rejected.
2. Fabricated chemical or experimental parameters are rejected.
3. Simple condition designs (~ condition) are verified with independent R.
4. Paired biological source designs (~ biological_source + condition) are verified.
5. Batch-adjusted designs (~ batch + condition) are verified.
6. Combined source and batch designs (~ biological_source + batch + condition) are verified.
7. Confounded statistical designs fail before execution.
8. Missing-value mismatches between pipeline and direct R are caught and rejected.
9. Missing genes in either pipeline or direct R are caught and rejected.
10. Design formula or rank mismatches between pipeline and direct R are rejected.
11. Treatment-vs-control directionality mismatches are caught and rejected.
12. Arbitrary custom valid datasets pass verification cleanly without source modifications.
13. Missing biological source identifiers when required by design are strictly rejected.
14. Missing batch identifiers when required by design are strictly rejected.
15. Conflicting metadata identifiers across fields are strictly rejected.
"""

import dataclasses
import json
from pathlib import Path
import subprocess
import tempfile
from typing import Any

import pandas as pd
import pytest

from scripts.independent_deseq2_verification import run_independent_verification
from src.app.pipeline_runner import ChemicalParams

# 15 canonical Ensembl GRCh38 genes for synthetic fixtures
SYNTHETIC_CANONICAL_GENES = [
    "ENSG00000000003",
    "ENSG00000000005",
    "ENSG00000000419",
    "ENSG00000000457",
    "ENSG00000000460",
    "ENSG00000000938",
    "ENSG00000000971",
    "ENSG00000001036",
    "ENSG00000001084",
    "ENSG00000001167",
    "ENSG00000001460",
    "ENSG00000001461",
    "ENSG00000001497",
    "ENSG00000001540",
    "ENSG00000001541",
]


def _make_counts_df(sample_ids: list[str], trt_sids: list[str]) -> pd.DataFrame:
    """Create deterministic counts with varied baselines and differential fold changes."""
    data: dict[str, Any] = {"gene_id": SYNTHETIC_CANONICAL_GENES}
    trt_set = set(trt_sids)
    for sid in sample_ids:
        counts = []
        is_trt = sid in trt_set
        s_idx = sample_ids.index(sid)
        for g_idx in range(len(SYNTHETIC_CANONICAL_GENES)):
            base = 100 * (g_idx + 1)
            noise = ((s_idx + 1) * (g_idx + 3) * 17) % 31 - 15
            if is_trt:
                fc = 2.5 if (g_idx % 2 == 1) else 0.4
                val = max(15, int(base * fc + noise))
            else:
                val = max(15, int(base + noise))
            counts.append(val)
        data[sid] = counts
    return pd.DataFrame(data)


def _make_valid_chem_params() -> ChemicalParams:
    """Create valid explicit ChemicalParams for synthetic testing."""
    return ChemicalParams(
        agent_name="Test_Compound_A",
        agent_identifier="CID:12345",
        concentration_or_dose=25.0,
        concentration_or_dose_unit="uM",
        developmental_age="Day 35",
        exposure_duration=48.0,
        exposure_duration_unit="h",
        vehicle="DMSO",
        organoid_context_id="ctx_synthetic_01",
        organoid_type="cerebral_organoid",
        brain_region="cortex",
        biological_source_id="donor_1",
        species="Homo sapiens",
    )


# 1. Missing Biological Replicate Metadata Fails
def test_missing_biological_replicate_metadata_fails():
    """Sample metadata lacking explicit biological replicate IDs must raise ValueError."""
    samples = ["s1", "s2", "s3", "s4"]
    counts_df = _make_counts_df(samples, trt_sids=["s3", "s4"])
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control", "control", "treatment", "treatment"],
            "biological_replicate_id": ["mock_rep_1", "", "mock_rep_3", "mock_rep_4"],  # s2 blank
        }
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)
        with pytest.raises(ValueError, match="biological replicate"):
            run_independent_verification(
                counts_path=counts_csv,
                samples_df=samples_df,
                chemical_params=_make_valid_chem_params(),
            )


# 2. Missing Required Experimental Metadata Fails
def test_missing_required_experimental_metadata_fails():
    """Missing ChemicalParams or missing required fields must raise ValueError."""
    samples = ["s1", "s2", "s3", "s4"]
    counts_df = _make_counts_df(samples, trt_sids=["s3", "s4"])
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control", "control", "treatment", "treatment"],
            "biological_replicate_id": ["mock_rep_1", "mock_rep_2", "mock_rep_3", "mock_rep_4"],
        }
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)

        # Neither chemical_params nor fallback
        with pytest.raises(ValueError, match="Chemical exposure and experimental parameters"):
            run_independent_verification(
                counts_path=counts_csv,
                samples_df=samples_df,
                chemical_params=None,
            )

        # Missing agent_name
        invalid_chem = ChemicalParams(
            agent_name="",
            concentration_or_dose=10.0,
            developmental_age="Day 30",
            exposure_duration=24.0,
        )
        with pytest.raises(ValueError, match="agent_name"):
            run_independent_verification(
                counts_path=counts_csv,
                samples_df=samples_df,
                chemical_params=invalid_chem,
            )


# 3. Valid Simple Treatment-Control Design (~ condition)
def test_valid_simple_treatment_control_design():
    """Verify that ~ condition design runs and matches 100% with direct R."""
    samples = ["c1", "c2", "t1", "t2"]
    counts_df = _make_counts_df(samples, trt_sids=["t1", "t2"])
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control", "control", "treatment", "treatment"],
            "biological_replicate_id": ["mock_rep_c1", "mock_rep_c2", "mock_rep_t1", "mock_rep_t2"],
            "biological_source_id": ["donor_1", "donor_1", "donor_1", "donor_1"],
        }
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)
        rep = run_independent_verification(
            counts_path=counts_csv,
            samples_df=samples_df,
            chemical_params=_make_valid_chem_params(),
        )
        assert rep["status"] == "PASS"
        assert rep["design_formula"] == "~ condition"
        assert rep["design_matrix_rank"] == 2
        assert rep["residual_degrees_of_freedom"] == 2
        assert rep["design_matches"] is True
        assert rep["rank_matches"] is True
        assert rep["df_matches"] is True
        assert rep["direction_agreement_pct"] == 100.0


# 4. Valid Paired Biological Source Design (~ biological_source + condition)
def test_valid_paired_biological_source_design():
    """Verify that ~ biological_source + condition runs and matches direct R."""
    samples = ["c_A1", "c_A2", "c_B1", "c_B2", "t_A1", "t_A2", "t_B1", "t_B2"]
    counts_df = _make_counts_df(samples, trt_sids=["t_A1", "t_A2", "t_B1", "t_B2"])
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control"] * 4 + ["treatment"] * 4,
            "biological_replicate_id": [f"mock_rep_{i}" for i in range(1, 9)],
            "biological_source_id": ["donor_A", "donor_A", "donor_B", "donor_B"] * 2,
        }
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)
        rep = run_independent_verification(
            counts_path=counts_csv,
            samples_df=samples_df,
            chemical_params=_make_valid_chem_params(),
        )
        assert rep["status"] == "PASS"
        assert rep["design_formula"] == "~ biological_source + condition"
        assert rep["design_matrix_rank"] == 3
        assert rep["residual_degrees_of_freedom"] == 5
        assert rep["design_matches"] is True
        assert rep["rank_matches"] is True
        assert rep["df_matches"] is True


# 5. Valid Batch-Adjusted Design (~ batch + condition)
def test_valid_batch_adjusted_design():
    """Verify that ~ batch + condition runs and matches direct R."""
    samples = ["c_b1_1", "c_b1_2", "c_b2_1", "c_b2_2", "t_b1_1", "t_b1_2", "t_b2_1", "t_b2_2"]
    counts_df = _make_counts_df(samples, trt_sids=["t_b1_1", "t_b1_2", "t_b2_1", "t_b2_2"])
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control"] * 4 + ["treatment"] * 4,
            "biological_replicate_id": [f"mock_rep_{i}" for i in range(1, 9)],
            "biological_source_id": ["donor_1"] * 8,
            "experimental_batch_id": ["batch_1", "batch_1", "batch_2", "batch_2"] * 2,
        }
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)
        rep = run_independent_verification(
            counts_path=counts_csv,
            samples_df=samples_df,
            chemical_params=_make_valid_chem_params(),
        )
        assert rep["status"] == "PASS"
        assert rep["design_formula"] == "~ batch + condition"
        assert rep["design_matrix_rank"] == 3
        assert rep["residual_degrees_of_freedom"] == 5
        assert rep["design_matches"] is True
        assert rep["rank_matches"] is True
        assert rep["df_matches"] is True


# 6. Valid Combined Source and Batch Design (~ biological_source + batch + condition)
def test_valid_combined_source_and_batch_design():
    """Verify that ~ biological_source + batch + condition runs and matches direct R."""
    samples = ["c_A1", "c_A2", "c_B1", "c_B2", "t_A1", "t_A2", "t_B1", "t_B2"]
    counts_df = _make_counts_df(samples, trt_sids=["t_A1", "t_A2", "t_B1", "t_B2"])
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control"] * 4 + ["treatment"] * 4,
            "biological_replicate_id": [f"mock_rep_{i}" for i in range(1, 9)],
            "biological_source_id": ["donor_A", "donor_A", "donor_B", "donor_B"] * 2,
            "experimental_batch_id": ["batch_1", "batch_2", "batch_1", "batch_2"] * 2,
        }
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)
        rep = run_independent_verification(
            counts_path=counts_csv,
            samples_df=samples_df,
            chemical_params=_make_valid_chem_params(),
        )
        assert rep["status"] == "PASS"
        assert rep["design_formula"] == "~ biological_source + batch + condition"
        assert rep["design_matrix_rank"] == 4
        assert rep["residual_degrees_of_freedom"] == 4
        assert rep["design_matches"] is True
        assert rep["rank_matches"] is True
        assert rep["df_matches"] is True


# 7. Confounded Statistical Design Fails
def test_confounded_statistical_design_fails():
    """Treatment completely confounded with biological source must block execution."""
    samples = ["c1", "c2", "t1", "t2"]
    counts_df = _make_counts_df(samples, trt_sids=["t1", "t2"])
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control", "control", "treatment", "treatment"],
            "biological_replicate_id": ["mock_rep_1", "mock_rep_2", "mock_rep_3", "mock_rep_4"],
            "biological_source_id": ["donor_A", "donor_A", "donor_B", "donor_B"],  # completely confounded
        }
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)
        with pytest.raises(RuntimeError, match="confounded"):
            run_independent_verification(
                counts_path=counts_csv,
                samples_df=samples_df,
                chemical_params=_make_valid_chem_params(),
            )


# 8. Numerical vs Missing Value Mismatch Fails
def test_numerical_vs_missing_value_mismatch_fails(monkeypatch):
    """Tampering direct R results so a numeric p-value becomes missing must raise ValueError."""
    samples = ["c1", "c2", "t1", "t2"]
    counts_df = _make_counts_df(samples, trt_sids=["t1", "t2"])
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control", "control", "treatment", "treatment"],
            "biological_replicate_id": ["mock_rep_1", "mock_rep_2", "mock_rep_3", "mock_rep_4"],
            "biological_source_id": ["donor_1"] * 4,
        }
    )
    real_run = subprocess.run

    def mock_subprocess_run(cmd, *args, **kwargs):
        res = real_run(cmd, *args, **kwargs)
        if any("independent_run.R" in str(arg) for arg in cmd):
            script_path = Path(cmd[-1])
            csv_path = script_path.parent / "direct_r_results.csv"
            if csv_path.exists():
                df = pd.read_csv(csv_path)
                df.loc[0, "p_value"] = None
                df.to_csv(csv_path, index=False)
        return res

    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)

    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)
        with pytest.raises(ValueError, match="Missing-value mismatch"):
            run_independent_verification(
                counts_path=counts_csv,
                samples_df=samples_df,
                chemical_params=_make_valid_chem_params(),
            )


# 9. Missing Gene Mismatch Fails
def test_missing_gene_mismatch_fails(monkeypatch):
    """Omitting a gene from direct R results must raise ValueError with gene universe mismatch."""
    samples = ["c1", "c2", "t1", "t2"]
    counts_df = _make_counts_df(samples, trt_sids=["t1", "t2"])
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control", "control", "treatment", "treatment"],
            "biological_replicate_id": ["mock_rep_1", "mock_rep_2", "mock_rep_3", "mock_rep_4"],
            "biological_source_id": ["donor_1"] * 4,
        }
    )
    real_run = subprocess.run

    def mock_subprocess_run(cmd, *args, **kwargs):
        res = real_run(cmd, *args, **kwargs)
        if any("independent_run.R" in str(arg) for arg in cmd):
            script_path = Path(cmd[-1])
            csv_path = script_path.parent / "direct_r_results.csv"
            if csv_path.exists():
                df = pd.read_csv(csv_path)
                df = df.iloc[1:]  # drop first gene
                df.to_csv(csv_path, index=False)
        return res

    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)

    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)
        with pytest.raises(ValueError, match="Gene universe mismatch"):
            run_independent_verification(
                counts_path=counts_csv,
                samples_df=samples_df,
                chemical_params=_make_valid_chem_params(),
            )


# 10. Incorrect Statistical Design Fails
def test_incorrect_statistical_design_fails(monkeypatch):
    """Direct R diagnostics reporting wrong design formula must raise ValueError."""
    samples = ["c1", "c2", "t1", "t2"]
    counts_df = _make_counts_df(samples, trt_sids=["t1", "t2"])
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control", "control", "treatment", "treatment"],
            "biological_replicate_id": ["mock_rep_1", "mock_rep_2", "mock_rep_3", "mock_rep_4"],
            "biological_source_id": ["donor_1"] * 4,
        }
    )
    real_run = subprocess.run

    def mock_subprocess_run(cmd, *args, **kwargs):
        res = real_run(cmd, *args, **kwargs)
        if any("independent_run.R" in str(arg) for arg in cmd):
            script_path = Path(cmd[-1])
            diag_path = script_path.parent / "direct_r_diagnostics.json"
            if diag_path.exists():
                with diag_path.open("r", encoding="utf-8") as f:
                    diag = json.load(f)
                diag["design_formula"] = "~ batch + condition"  # tamper formula
                with diag_path.open("w", encoding="utf-8") as f:
                    json.dump(diag, f)
        return res

    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)

    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)
        with pytest.raises(ValueError, match="Statistical design mismatch"):
            run_independent_verification(
                counts_path=counts_csv,
                samples_df=samples_df,
                chemical_params=_make_valid_chem_params(),
            )


# 11. Treatment-Control Direction Mismatch Fails
def test_treatment_control_direction_mismatch_fails(monkeypatch):
    """Inverting log2FC sign between pipeline and direct R must raise ValueError."""
    samples = ["c1", "c2", "t1", "t2"]
    counts_df = _make_counts_df(samples, trt_sids=["t1", "t2"])
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control", "control", "treatment", "treatment"],
            "biological_replicate_id": ["mock_rep_1", "mock_rep_2", "mock_rep_3", "mock_rep_4"],
            "biological_source_id": ["donor_1"] * 4,
        }
    )
    real_run = subprocess.run

    def mock_subprocess_run(cmd, *args, **kwargs):
        res = real_run(cmd, *args, **kwargs)
        if any("independent_run.R" in str(arg) for arg in cmd):
            script_path = Path(cmd[-1])
            csv_path = script_path.parent / "direct_r_results.csv"
            if csv_path.exists():
                df = pd.read_csv(csv_path)
                df.loc[0, "log2_fold_change"] = -1.0 * float(df.loc[0, "log2_fold_change"])
                df.to_csv(csv_path, index=False)
        return res

    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)

    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)
        with pytest.raises(ValueError, match="Treatment-vs-control direction mismatch"):
            run_independent_verification(
                counts_path=counts_csv,
                samples_df=samples_df,
                chemical_params=_make_valid_chem_params(),
            )


# 12. Custom Valid Dataset Passes Cleanly
def test_custom_valid_dataset_passes():
    """Verify that a completely arbitrary custom valid dataset runs without code changes."""
    samples = ["sample_A1", "sample_A2", "sample_B1", "sample_B2"]
    counts_df = _make_counts_df(samples, trt_sids=["sample_B1", "sample_B2"])
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control", "control", "treatment", "treatment"],
            "biological_replicate_id": ["mock_rep_A1", "mock_rep_A2", "mock_rep_B1", "mock_rep_B2"],
            "biological_source_id": ["donor_X", "donor_X", "donor_X", "donor_X"],
        }
    )
    chem_params = ChemicalParams(
        agent_name="Novel_Compound_Z",
        agent_identifier="CID:77777",
        concentration_or_dose=100.0,
        concentration_or_dose_unit="nM",
        developmental_age="Day 60",
        exposure_duration=96.0,
        exposure_duration_unit="h",
        vehicle="Ethanol",
        organoid_context_id="ctx_novel_01",
        organoid_type="cortical_spheroid",
        brain_region="forebrain",
        biological_source_id="donor_X",
        species="Homo sapiens",
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        counts_csv = tmp / "custom_counts.csv"
        out_dir = tmp / "verification_output"
        counts_df.to_csv(counts_csv, index=False)

        report = run_independent_verification(
            counts_path=counts_csv,
            samples_df=samples_df,
            chemical_params=chem_params,
            output_dir=out_dir,
            tolerance=1e-6,
        )

        assert report["status"] == "PASS"
        assert report["genes_compared"] == report["genes_evaluated"]
        assert report["genes_compared"] > 0
        assert report["missing_genes_count"] == 0
        assert report["max_base_mean_diff"] <= 1e-6
        assert report["max_log2fc_diff"] <= 1e-6
        assert report["max_lfc_se_diff"] <= 1e-6
        assert report["max_wald_stat_diff"] <= 1e-6
        assert report["max_pvalue_diff"] <= 1e-6
        assert report["max_padj_diff"] <= 1e-6
        assert report["direction_agreement_pct"] == 100.0

        assert (out_dir / "direct_r_results.csv").exists()
        assert (out_dir / "pipeline_results.csv").exists()
        assert (out_dir / "verification_report.json").exists()


# 13. Missing Biological Source Identifier Rejected When Required
def test_missing_biological_source_identifier_rejected_when_required():
    """Missing biological source identifier when design requires biological_source must raise ValueError."""
    samples = ["c_A1", "c_A2", "c_B1", "c_B2", "t_A1", "t_A2", "t_B1", "t_B2"]
    counts_df = _make_counts_df(samples, trt_sids=["t_A1", "t_A2", "t_B1", "t_B2"])
    sources = ["donor_A", "donor_A", "donor_B", "donor_B"] * 2
    sources[0] = ""  # sample c_A1 lacks explicit biological source
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control"] * 4 + ["treatment"] * 4,
            "biological_replicate_id": [f"mock_rep_{i}" for i in range(1, 9)],
            "biological_source_id": sources,
        }
    )
    chem_params = dataclasses.replace(_make_valid_chem_params(), biological_source_id="donor_A")

    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)
        with pytest.raises(ValueError, match="requires 'biological_source'"):
            run_independent_verification(
                counts_path=counts_csv,
                samples_df=samples_df,
                chemical_params=chem_params,
            )


# 14. Missing Batch Identifier Rejected When Required
def test_missing_batch_identifier_rejected_when_required(monkeypatch):
    """Missing batch identifier when design requires batch must raise ValueError."""
    samples = ["c_b1_1", "c_b1_2", "c_b2_1", "c_b2_2", "t_b1_1", "t_b1_2", "t_b2_1", "t_b2_2"]
    counts_df = _make_counts_df(samples, trt_sids=["t_b1_1", "t_b1_2", "t_b2_1", "t_b2_2"])
    batches = ["batch_1", "batch_1", "batch_2", "batch_2"] * 2
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control"] * 4 + ["treatment"] * 4,
            "biological_replicate_id": [f"mock_rep_{i}" for i in range(1, 9)],
            "biological_source_id": ["donor_1"] * 8,
            "experimental_batch_id": batches,
        }
    )
    import scripts.independent_deseq2_verification as idv
    real_run_pipeline = idv.run_full_analysis_pipeline

    def mock_run_pipeline(*args, **kwargs):
        kwargs_copy = dict(kwargs)
        kwargs_copy["samples_input"] = samples_df.copy()
        return real_run_pipeline(*args, **kwargs_copy)

    monkeypatch.setattr(idv, "run_full_analysis_pipeline", mock_run_pipeline)

    samples_df_missing = samples_df.copy()
    samples_df_missing.loc[0, "experimental_batch_id"] = ""

    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)
        with pytest.raises(ValueError, match="requires 'batch'"):
            run_independent_verification(
                counts_path=counts_csv,
                samples_df=samples_df_missing,
                chemical_params=_make_valid_chem_params(),
            )


# 15. Conflicting Metadata Identifiers Rejected
def test_conflicting_metadata_identifiers_rejected():
    """Conflicting biological source or batch identifiers in sample metadata must raise ValueError."""
    samples = ["c_A1", "c_A2", "c_B1", "c_B2", "t_A1", "t_A2", "t_B1", "t_B2"]
    counts_df = _make_counts_df(samples, trt_sids=["t_A1", "t_A2", "t_B1", "t_B2"])
    sources = ["donor_A", "donor_A", "donor_B", "donor_B"] * 2
    cell_lines = ["line_conflict"] + sources[1:]
    samples_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition_type": ["control"] * 4 + ["treatment"] * 4,
            "biological_replicate_id": [f"mock_rep_{i}" for i in range(1, 9)],
            "biological_source_id": sources,
            "cell_line_id": cell_lines,
        }
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        counts_csv = Path(tmpdir) / "counts.csv"
        counts_df.to_csv(counts_csv, index=False)
        with pytest.raises(ValueError, match="conflicting biological source"):
            run_independent_verification(
                counts_path=counts_csv,
                samples_df=samples_df,
                chemical_params=_make_valid_chem_params(),
            )

