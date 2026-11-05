"""Regression test suite for evidence-based RSeQC strandedness evaluation and biological gating."""

import json
from pathlib import Path
import pytest

from src.quantification_bridge.manifest import (
    APPROVED_FEATURECOUNTS_STRANDEDNESS_FLAG,
    APPROVED_PILOT_RUNS,
    APPROVED_RUN_ACCESSIONS,
    APPROVED_STRANDEDNESS_MODE,
)
from src.quantification_bridge.pipeline import QuantificationBridgePipeline
from src.quantification_bridge.rseqc_runner import (
    StrandednessDecision,
    infer_strandedness_decision,
)
from src.quantification_bridge.downstream_bridge import execute_downstream_pipeline


def test_infer_strandedness_decision_unstranded_summary():
    """Verify that balanced forward/reverse fractions evaluate to unstranded (flag 0)."""
    summary = {
        "mean_forward_fraction": 0.4041,
        "mean_reverse_fraction": 0.3982,
        "mean_unexplained_fraction": 0.1978,
    }
    decision = infer_strandedness_decision(summary)
    assert isinstance(decision, StrandednessDecision)
    assert decision.strandedness_flag == 0
    assert decision.strand_mode == "unstranded"
    assert decision.confidence == "HIGH"
    assert "UNSTRANDED" in decision.consensus_indication
    assert "diff =" in decision.rationale


def test_infer_strandedness_decision_reverse_summary():
    """Verify that dominant antisense fractions evaluate to reverse-stranded (flag 2)."""
    summary = {
        "mean_forward_fraction": 0.052,
        "mean_reverse_fraction": 0.884,
        "mean_unexplained_fraction": 0.064,
    }
    decision = infer_strandedness_decision(summary)
    assert decision.strandedness_flag == 2
    assert decision.strand_mode == "reverse"
    assert decision.confidence == "HIGH"
    assert "REVERSE_STRANDED" in decision.consensus_indication


def test_infer_strandedness_decision_forward_summary():
    """Verify that dominant sense fractions evaluate to forward-stranded (flag 1)."""
    summary = {
        "mean_forward_fraction": 0.891,
        "mean_reverse_fraction": 0.045,
        "mean_unexplained_fraction": 0.064,
    }
    decision = infer_strandedness_decision(summary)
    assert decision.strandedness_flag == 1
    assert decision.strand_mode == "forward"
    assert decision.confidence == "HIGH"
    assert "FORWARD_STRANDED" in decision.consensus_indication


def test_infer_strandedness_decision_ambiguous():
    """Verify that ambiguous fractions halt execution with flag -1."""
    summary = {
        "mean_forward_fraction": 0.55,
        "mean_reverse_fraction": 0.20,
        "mean_unexplained_fraction": 0.25,
    }
    decision = infer_strandedness_decision(summary)
    assert decision.strandedness_flag == -1
    assert decision.strand_mode == "ambiguous"
    assert decision.confidence == "AMBIGUOUS"
    assert "AMBIGUOUS_EVIDENCE_REQUIRES_REVIEW" in decision.consensus_indication


def test_authoritative_prjna802072_evidence_artifact():
    """Verify that authoritative references/prjna802072_strandedness_evidence.json evaluates to -s 0."""
    repo_root = Path(__file__).resolve().parents[2]
    evidence_path = repo_root / "references" / "prjna802072_strandedness_evidence.json"
    assert evidence_path.exists(), f"Evidence artifact missing: {evidence_path}"

    data = json.loads(evidence_path.read_text(encoding="utf-8"))
    assert data["total_samples"] == 12
    assert len(data["per_sample_strandedness"]) == 12

    summary = data["strandedness_summary"]
    per_sample = data["per_sample_strandedness"]

    # Assert exact 12-library empirical summary statistics
    assert summary["mean_forward_fraction"] == 0.4019
    assert summary["mean_reverse_fraction"] == 0.4025
    assert summary["mean_unexplained_fraction"] == 0.1956
    assert abs(summary["mean_forward_fraction"] - summary["mean_reverse_fraction"]) < 0.001

    # Evaluate summary
    decision_summary = infer_strandedness_decision(summary)
    assert decision_summary.strandedness_flag == 0
    assert decision_summary.strand_mode == APPROVED_STRANDEDNESS_MODE
    assert decision_summary.strandedness_flag == APPROVED_FEATURECOUNTS_STRANDEDNESS_FLAG
    assert decision_summary.confidence == "HIGH"

    # Evaluate per-sample dictionary across all 12 runs
    decision_per_sample = infer_strandedness_decision(per_sample)
    assert decision_per_sample.strandedness_flag == 0
    assert decision_per_sample.strand_mode == APPROVED_STRANDEDNESS_MODE
    assert decision_per_sample.strandedness_flag == APPROVED_FEATURECOUNTS_STRANDEDNESS_FLAG
    assert decision_per_sample.confidence == "HIGH"


def test_pipeline_evaluate_strandedness_decision():
    """Verify QuantificationBridgePipeline exposes evaluate_strandedness_decision method."""
    pipeline = QuantificationBridgePipeline()
    summary = {
        "mean_forward_fraction": 0.4041,
        "mean_reverse_fraction": 0.3982,
        "mean_unexplained_fraction": 0.1978,
    }
    decision = pipeline.evaluate_strandedness_decision(summary)
    assert decision.strandedness_flag == 0
    assert decision.strand_mode == "unstranded"


def test_downstream_biological_replication_gate_full_cohort(tmp_path: Path):
    """Verify downstream pipeline gates DESeq2 for full 12-sample cohort when biological replicates are unverified."""
    matrix_file = tmp_path / "counts_12_samples.tsv"
    runs = APPROVED_RUN_ACCESSIONS
    assert len(runs) == 12

    # Construct synthetic count matrix for all 12 samples
    header = "gene_id\t" + "\t".join(runs)
    rows = [
        header,
        "ENSG00000000003\t" + "\t".join(["100"] * 12),
        "ENSG00000000005\t" + "\t".join(["0"] * 12),
        "ENSG00000000419\t" + "\t".join(["500"] * 12),
        "ENSG00000000457\t" + "\t".join(["250"] * 12),
    ]
    matrix_file.write_text("\n".join(rows) + "\n", encoding="utf-8")

    result = execute_downstream_pipeline(
        counts_matrix_path=matrix_file,
        runs=runs,
        output_dir=tmp_path / "downstream",
    )

    # 1. Pipeline status reflects partial biological replication gate
    assert result.status == "DOWNSTREAM_PARTIAL_BIOLOGICAL_REPLICATION_GATED"
    assert result.qc_passed is True

    # 2. Advisory QC flags capture SRR17834908
    assert "advisory_flags" in result.qc_metrics
    assert any("SRR17834908" in flag for flag in result.qc_metrics["advisory_flags"])

    # 3. Steps 6B and 6C executed and passed
    assert result.harmonization_metrics["status"] == "PASS"
    assert result.normalization_metrics["status"] == "PASS"

    # 4. Inferential DESeq2 and Pathway analysis are gated to prevent unverified claims
    assert result.deseq2_status == "BLOCKED_UNVERIFIED_BIOLOGICAL_REPLICATION"
    assert result.pathway_status == "BLOCKED_UNVERIFIED_BIOLOGICAL_REPLICATION"
    assert result.biological_replicate_warning is not None
    assert "Biological replication is UNVERIFIED" in result.biological_replicate_warning


REAL_12_SAMPLE_MATRIX = Path(
    "G:/My Drive/DNT_Organoid_AI_Persistence/stage18_full_cohort/count_results/prjna802072_cohort_12_samples_count_matrix.tsv"
)


@pytest.mark.skipif(not REAL_12_SAMPLE_MATRIX.exists(), reason="Real 12-sample matrix not mounted")
def test_downstream_biological_replication_gate_real_12_sample_matrix(tmp_path: Path):
    """Verify real 63,140-gene x 12-sample matrix executes QC, Harmonization, Normalization, and gates DESeq2."""
    runs = APPROVED_RUN_ACCESSIONS
    assert len(runs) == 12

    result = execute_downstream_pipeline(
        counts_matrix_path=REAL_12_SAMPLE_MATRIX,
        runs=runs,
        output_dir=tmp_path / "downstream_real",
    )

    # 1. Pipeline status reflects partial biological replication gate
    assert result.status == "DOWNSTREAM_PARTIAL_BIOLOGICAL_REPLICATION_GATED"
    assert result.qc_passed is True
    assert result.qc_metrics["total_genes"] == 63140
    assert result.qc_metrics["total_samples"] == 12

    # 2. Advisory QC flags capture SRR17834908
    assert "advisory_flags" in result.qc_metrics
    assert any("SRR17834908" in flag for flag in result.qc_metrics["advisory_flags"])

    # 3. Steps 6B and 6C executed and passed
    assert result.harmonization_metrics["status"] == "PASS"
    assert result.harmonization_metrics["mapped_genes_count"] == 63140
    assert result.normalization_metrics["status"] == "PASS"

    # 4. Inferential DESeq2 and Pathway analysis are gated to prevent unverified claims
    assert result.deseq2_status == "BLOCKED_UNVERIFIED_BIOLOGICAL_REPLICATION"
    assert result.pathway_status == "BLOCKED_UNVERIFIED_BIOLOGICAL_REPLICATION"
    assert result.biological_replicate_warning is not None
    assert "Biological replication is UNVERIFIED" in result.biological_replicate_warning
