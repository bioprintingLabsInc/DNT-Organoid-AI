"""Unit tests for the independent DESeq2 verification CLI and engine."""

from pathlib import Path
import tempfile
import pytest

import pandas as pd

from scripts.independent_deseq2_verification import run_independent_verification
from src.app.pipeline_runner import ChemicalParams


def test_independent_verification_missing_counts_file():
    """Missing count matrix file raises FileNotFoundError with clear message."""
    with pytest.raises(FileNotFoundError, match="Raw count matrix file not found"):
        run_independent_verification(
            counts_path="non_existent_matrix_file.csv",
            treatment_sids=["s1", "s2"],
            control_sids=["s3", "s4"],
        )


def test_independent_verification_missing_metadata():
    """Missing sample metadata and missing treatment/control IDs raises ValueError."""
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        f.write("gene_id\ts1\ts2\nENSG00000000003\t10\t20\n")
        temp_path = f.name

    try:
        with pytest.raises(ValueError, match="Either samples_path, samples_df, or both"):
            run_independent_verification(counts_path=temp_path)
    finally:
        Path(temp_path).unlink(missing_ok=True)


def test_independent_verification_insufficient_replicates():
    """Fewer than 2 biological replicates raises ValueError."""
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        f.write("gene_id\ts1\ts2\ts3\nENSG00000000003\t10\t20\t30\n")
        temp_path = f.name

    try:
        with pytest.raises(ValueError, match="requires >= 2 replicates per arm"):
            run_independent_verification(
                counts_path=temp_path,
                treatment_sids=["s1"],  # only 1 replicate
                control_sids=["s2", "s3"],
            )
    finally:
        Path(temp_path).unlink(missing_ok=True)


def test_independent_verification_different_valid_dataset():
    """Verify that a completely different valid input dataset runs without code modification."""
    # Create synthetic raw count matrix with 15 known canonical genes
    genes = [
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
    # 2 treatment vs 2 control samples
    data = {
        "gene_id": genes,
        "sample_A1": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500],
        "sample_A2": [120, 210, 310, 390, 510, 620, 690, 810, 890, 1020, 1090, 1210, 1290, 1410, 1490],
        "sample_B1": [50, 400, 150, 800, 250, 1200, 350, 1600, 450, 2000, 550, 2400, 650, 2800, 750],
        "sample_B2": [60, 390, 160, 790, 260, 1190, 360, 1590, 460, 1990, 560, 2390, 660, 2790, 760],
    }
    df = pd.DataFrame(data)

    samples_data = {
        "sample_id": ["sample_A1", "sample_A2", "sample_B1", "sample_B2"],
        "condition_type": ["control", "control", "treatment", "treatment"],
        "biological_replicate_id": ["rep_ctrl_1", "rep_ctrl_2", "rep_trt_1", "rep_trt_2"],
        "biological_source_id": ["donor_1", "donor_1", "donor_1", "donor_1"],
    }
    samples_df = pd.DataFrame(samples_data)

    chem_params = ChemicalParams(
        agent_name="Compound_X",
        agent_identifier="CID:99999",
        concentration_or_dose=50.0,
        concentration_or_dose_unit="nM",
        developmental_age="Day 45",
        exposure_duration=72.0,
        exposure_duration_unit="h",
        vehicle="DMSO",
        organoid_context_id="organoid_ctx_custom",
        organoid_type="cortical_spheroid",
        brain_region="forebrain",
        biological_source_id="donor_1",
        species="Homo sapiens",
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        counts_csv = tmp / "custom_counts.csv"
        samples_csv = tmp / "custom_samples.csv"
        out_dir = tmp / "out"

        df.to_csv(counts_csv, index=False)
        samples_df.to_csv(samples_csv, index=False)

        report = run_independent_verification(
            counts_path=counts_csv,
            samples_path=samples_csv,
            chemical_params=chem_params,
            output_dir=out_dir,
            tolerance=1e-6,
        )

        assert report["status"] == "PASS"
        assert report["genes_compared"] == report["genes_evaluated"]
        assert report["genes_compared"] > 0
        assert report["missing_genes_count"] == 0
        assert report["max_log2fc_diff"] < 1e-6
        assert report["max_wald_stat_diff"] < 1e-6
        assert report["max_pvalue_diff"] < 1e-6
        assert report["max_padj_diff"] < 1e-6
        assert report["direction_agreement_pct"] == 100.0

        assert (out_dir / "direct_r_results.csv").exists()
        assert (out_dir / "pipeline_results.csv").exists()
        assert (out_dir / "verification_report.json").exists()
