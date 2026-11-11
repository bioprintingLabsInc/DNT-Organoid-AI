"""Unit and integration test suite for Lexogen QuantSeq 3' mRNA-seq pipeline support.

Verifies:
1. Library type recognition and lossless metadata classification.
2. fastp argument construction (QuantSeq FWD/REV vs conventional bulk).
3. Strandedness determination, single-end RSeQC output parsing, and consensus gating.
4. featureCounts flag construction (single-end vs paired-end, strandedness).
5. UMI safeguard blocking when UMI pipeline is unvalidated.
6. Cohort library-type consistency preflight checks.
7. Full pipeline integration with QuantSeq FWD single-end data.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from unittest.mock import MagicMock, patch
import pytest

from src.quantification_bridge.library_type import LibraryType
from src.quantification_bridge.fastp_runner import (
    LEXOGEN_QUANTSEQ_FWD_ADAPTER,
    run_fastp_preprocessing,
)
from src.quantification_bridge.featurecounts_runner import run_featurecounts
from src.quantification_bridge.rseqc_runner import (
    StrandednessDecision,
    infer_strandedness_decision,
    run_infer_experiment,
)
from src.quantification_bridge.universal_pipeline import (
    ReferenceConfig,
    SampleInput,
    UniversalFastqPipeline,
    UniversalPipelineConfig,
)


def make_dummy_fastq(
    path: Path,
    read_names: list[str],
    read_len: int = 50,
    mate_suffix: str = "",
    is_gz: bool = True,
) -> Path:
    """Create a structurally valid synthetic FASTQ file for preflight checks."""
    lines = []
    for rname in read_names:
        hdr = f"@{rname}{mate_suffix}"
        seq = "A" * read_len
        qual = "I" * read_len
        lines.extend([hdr, seq, "+", qual])
    content = "\n".join(lines) + "\n"
    if is_gz:
        with gzip.open(path, "wt", encoding="utf-8") as f:
            f.write(content)
    else:
        path.write_text(content, encoding="utf-8")
    return path


# -----------------------------------------------------------------------------
# 1. Library Type Recognition Tests
# -----------------------------------------------------------------------------

def test_library_type_enum_parsing():
    """Verify robust, lossless recognition of QuantSeq and conventional library types."""
    # Conventional bulk synonyms
    assert LibraryType.from_string("CONVENTIONAL_BULK") == LibraryType.CONVENTIONAL_BULK
    assert LibraryType.from_string("conventional") == LibraryType.CONVENTIONAL_BULK
    assert LibraryType.from_string("standard_bulk") == LibraryType.CONVENTIONAL_BULK
    assert LibraryType.from_string("polyA") == LibraryType.CONVENTIONAL_BULK
    assert LibraryType.from_string("TruSeq Stranded mRNA") == LibraryType.CONVENTIONAL_BULK
    assert LibraryType.from_string("bulk_rna_seq") == LibraryType.CONVENTIONAL_BULK

    # QuantSeq FWD variations (including strings mentioning Illumina)
    assert LibraryType.from_string("LEXOGEN_QUANTSEQ_3UTR_FWD") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD
    assert LibraryType.from_string("quantseq_fwd") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD
    assert LibraryType.from_string("QuantSeq 3' mRNA-Seq FWD") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD
    assert LibraryType.from_string("3' mRNA-Seq Library Prep Kit FWD for Illumina (QuantSeq-LEXOGEN)") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD
    assert LibraryType.from_string("Lexogen QuantSeq 3' FWD") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD

    # QuantSeq REV variations
    assert LibraryType.from_string("LEXOGEN_QUANTSEQ_3UTR_REV") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV
    assert LibraryType.from_string("quantseq_rev") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV
    assert LibraryType.from_string("QuantSeq 3' mRNA-Seq REV") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV
    assert LibraryType.from_string("3' mRNA-Seq Library Prep Kit REV (QuantSeq-LEXOGEN)") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV

    # UMI variants
    assert LibraryType.from_string("LEXOGEN_QUANTSEQ_3UTR_FWD_UMI") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI
    assert LibraryType.from_string("QuantSeq FWD with UMI") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI
    assert LibraryType.from_string("QuantSeq UMI FWD") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI
    assert LibraryType.from_string("LEXOGEN_QUANTSEQ_3UTR_REV_UMI") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV_UMI
    assert LibraryType.from_string("QuantSeq REV with UMI") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV_UMI

    # QuantSeq with undocumented orientation: strictly keeps orientation UNKNOWN
    assert LibraryType.from_string("LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND
    assert LibraryType.from_string("QuantSeq") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND
    assert LibraryType.from_string("Lexogen QuantSeq") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND
    assert LibraryType.from_string("3' mRNA-Seq") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND
    assert LibraryType.from_string("QuantSeq with UMI") == LibraryType.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND_UMI

    # Non-QuantSeq Lexogen kits must NOT be classified as QuantSeq
    assert LibraryType.from_string("Lexogen CORALL Total RNA-Seq") == LibraryType.CONVENTIONAL_BULK
    assert LibraryType.from_string("CORALL") == LibraryType.CONVENTIONAL_BULK
    assert LibraryType.from_string("Lexogen") == LibraryType.UNKNOWN  # Unspecified Lexogen kit is UNKNOWN

    # Unknown fallback (does not invent false kits)
    assert LibraryType.from_string(None) == LibraryType.UNKNOWN
    assert LibraryType.from_string("") == LibraryType.UNKNOWN
    assert LibraryType.from_string("custom_unspecified_library") == LibraryType.UNKNOWN


def test_library_type_properties():
    """Verify directional and methodological properties of LibraryType."""
    fwd = LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD
    assert fwd.is_quantseq is True
    assert fwd.is_fwd is True
    assert fwd.is_rev is False
    assert fwd.has_umi is False
    assert fwd.expected_strandedness_flag == 1
    assert fwd.requires_length_normalization is False

    rev = LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV
    assert rev.is_quantseq is True
    assert rev.is_fwd is False
    assert rev.is_rev is True
    assert rev.has_umi is False
    assert rev.expected_strandedness_flag == 2
    assert rev.requires_length_normalization is False

    fwd_umi = LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI
    assert fwd_umi.has_umi is True
    assert fwd_umi.is_fwd is True
    assert fwd_umi.requires_length_normalization is False

    unk_strand = LibraryType.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND
    assert unk_strand.is_quantseq is True
    assert unk_strand.is_fwd is False
    assert unk_strand.is_rev is False
    assert unk_strand.expected_strandedness_flag is None
    assert unk_strand.has_umi is False
    assert unk_strand.requires_length_normalization is False

    bulk = LibraryType.CONVENTIONAL_BULK
    assert bulk.is_quantseq is False
    assert bulk.has_umi is False
    assert bulk.requires_length_normalization is True


# -----------------------------------------------------------------------------
# 2. fastp Command Construction Tests
# -----------------------------------------------------------------------------

def test_fastp_command_conventional_bulk_unmodified():
    """Verify that conventional bulk RNA-seq fastp flags remain 100% unchanged."""
    executed_cmds = []

    def mock_run_wsl(cmd):
        executed_cmds.append(cmd)
        m = MagicMock()
        m.returncode = 0
        m.stderr = ""
        m.stdout = ""
        return m

    with patch("src.quantification_bridge.fastp_runner.run_wsl", side_effect=mock_run_wsl), \
         patch("pathlib.Path.exists", return_value=True), \
         patch("pathlib.Path.is_file", return_value=True):
        
        # Paired-end conventional bulk
        run_fastp_preprocessing(
            raw_r1="/tmp/r1.fq.gz",
            raw_r2="/tmp/r2.fq.gz",
            out_r1="/tmp/trimmed_r1.fq.gz",
            out_r2="/tmp/trimmed_r2.fq.gz",
            sample_id="test_bulk",
            threads=4,
            library_type=LibraryType.CONVENTIONAL_BULK,
        )

        pe_cmd = executed_cmds[-1]
        assert "--detect_adapter_for_pe" in pe_cmd
        assert "--trim_poly_g" in pe_cmd
        assert "--poly_g_min_len 10" in pe_cmd
        assert "--disable_quality_filtering" in pe_cmd
        assert "--disable_length_filtering" in pe_cmd
        assert "--trim_poly_x" not in pe_cmd
        assert "--adapter_sequence" not in pe_cmd


def test_fastp_command_quantseq_fwd():
    """Verify that QuantSeq FWD receives adapter, poly(A) tail trimming, and quality filters without 5' hard clipping."""
    executed_cmds = []

    def mock_run_wsl(cmd):
        executed_cmds.append(cmd)
        m = MagicMock()
        m.returncode = 0
        m.stderr = ""
        m.stdout = ""
        return m

    with patch("src.quantification_bridge.fastp_runner.run_wsl", side_effect=mock_run_wsl), \
         patch("pathlib.Path.exists", return_value=True), \
         patch("pathlib.Path.is_file", return_value=True):
        
        run_fastp_preprocessing(
            raw_r1="/tmp/r1.fq.gz",
            raw_r2=None,
            out_r1="/tmp/trimmed_r1.fq.gz",
            sample_id="test_quantseq",
            threads=4,
            library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
        )

        cmd = executed_cmds[-1]
        # Must have default Lexogen adapter
        assert f"--adapter_sequence {LEXOGEN_QUANTSEQ_FWD_ADAPTER}" in cmd
        # Must trim poly(A) tail
        assert "--trim_poly_x" in cmd
        assert "--poly_x_min_len 10" in cmd
        # Must trim poly(G)
        assert "--trim_poly_g" in cmd
        # Must require minimum read length (default 20 bp)
        assert "--length_required 20" in cmd
        assert "--disable_length_filtering" not in cmd
        # Must apply quality trimming
        assert "--cut_tail" in cmd
        assert "--cut_mean_quality 20" in cmd
        assert "--disable_quality_filtering" not in cmd
        # Must NOT hard-clip 5' by default (STAR soft-clipping handles mismatches)
        assert "--trim_front1" not in cmd


def test_fastp_command_quantseq_with_custom_trim_front():
    """Verify that 5' hard clipping can be explicitly configured when requested."""
    executed_cmds = []

    def mock_run_wsl(cmd):
        executed_cmds.append(cmd)
        m = MagicMock()
        m.returncode = 0
        return m

    with patch("src.quantification_bridge.fastp_runner.run_wsl", side_effect=mock_run_wsl), \
         patch("pathlib.Path.exists", return_value=True):
        
        run_fastp_preprocessing(
            raw_r1="/tmp/r1.fq.gz",
            out_r1="/tmp/trimmed_r1.fq.gz",
            sample_id="test_quantseq_clip",
            library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
            trim_front1=12,
        )

        cmd = executed_cmds[-1]
        assert "--trim_front1 12" in cmd


# -----------------------------------------------------------------------------
# 3. RSeQC Single-End Parsing and Strandedness Consensus Gating Tests
# -----------------------------------------------------------------------------

def test_rseqc_parsing_single_end_output():
    """Verify that run_infer_experiment parses single-end RSeQC output ('++,--' and '+-,-+')."""
    single_end_rseqc_stdout = """
This is SingleEnd Data
Fraction of reads failed to determine: 0.0520
Fraction of reads explained by "++,--": 0.9125
Fraction of reads explained by "+-,-+": 0.0355
"""
    mock_res = MagicMock()
    mock_res.returncode = 0
    mock_res.stdout = single_end_rseqc_stdout
    mock_res.stderr = ""

    with patch("src.quantification_bridge.rseqc_runner.run_wsl", return_value=mock_res), \
         patch("pathlib.Path.exists", return_value=True), \
         patch("src.quantification_bridge.rseqc_runner.sha256_file", return_value="dummy_hash"):
        
        metrics = run_infer_experiment("/tmp/sample.bam", "/tmp/ref.bed12")
        assert metrics["forward_fraction"] == 0.9125
        assert metrics["reverse_fraction"] == 0.0355
        assert metrics["unexplained_fraction"] == 0.0520


def test_strandedness_consensus_quantseq_fwd_agreement():
    """Verify consensus verified forward when QuantSeq FWD has empirical sense reads."""
    metrics = {
        "mean_forward_fraction": 0.92,
        "mean_reverse_fraction": 0.04,
        "mean_unexplained_fraction": 0.04,
    }
    decision = infer_strandedness_decision(metrics, library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD)
    assert decision.strandedness_flag == 1
    assert decision.strand_mode == "forward"
    assert decision.confidence == "HIGH"
    assert "CONSENSUS_VERIFIED_FORWARD" in decision.consensus_indication


def test_strandedness_consensus_quantseq_fwd_conflict():
    """Verify that a conflict between declared QuantSeq FWD and empirical antisense halts downstream (-1)."""
    # Empirical reverse-stranded reads despite declared QuantSeq FWD
    metrics = {
        "mean_forward_fraction": 0.03,
        "mean_reverse_fraction": 0.93,
        "mean_unexplained_fraction": 0.04,
    }
    decision = infer_strandedness_decision(metrics, library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD)
    assert decision.strandedness_flag == -1
    assert decision.confidence == "CONFLICT"
    assert "CONFLICT_DECLARED_FWD_EMPIRICAL_REVERSE" in decision.consensus_indication
    assert "Halting downstream quantification" in decision.rationale


def test_strandedness_consensus_quantseq_rev_agreement():
    """Verify consensus verified reverse when QuantSeq REV has empirical antisense reads."""
    metrics = {
        "mean_forward_fraction": 0.05,
        "mean_reverse_fraction": 0.90,
        "mean_unexplained_fraction": 0.05,
    }
    decision = infer_strandedness_decision(metrics, library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV)
    assert decision.strandedness_flag == 2
    assert decision.strand_mode == "reverse"
    assert decision.confidence == "HIGH"
    assert "CONSENSUS_VERIFIED_REVERSE" in decision.consensus_indication


def test_strandedness_consensus_quantseq_rev_conflict():
    """Verify that a conflict between declared QuantSeq REV and empirical sense halts downstream (-1)."""
    metrics = {
        "mean_forward_fraction": 0.88,
        "mean_reverse_fraction": 0.06,
        "mean_unexplained_fraction": 0.06,
    }
    decision = infer_strandedness_decision(metrics, library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV)
    assert decision.strandedness_flag == -1
    assert decision.confidence == "CONFLICT"
    assert "CONFLICT_DECLARED_REV_EMPIRICAL_FORWARD" in decision.consensus_indication


def test_strandedness_undocumented_quantseq_resolved_by_empirical_forward():
    """Verify that undocumented QuantSeq orientation is resolved to forward when empirical evidence is sense."""
    metrics = {
        "mean_forward_fraction": 0.91,
        "mean_reverse_fraction": 0.04,
        "mean_unexplained_fraction": 0.05,
    }
    decision = infer_strandedness_decision(metrics, library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND)
    assert decision.strandedness_flag == 1
    assert decision.strand_mode == "forward"
    assert decision.confidence == "HIGH"
    assert "EMPIRICAL_QUANTSEQ_FORWARD" in decision.consensus_indication


def test_strandedness_undocumented_quantseq_resolved_by_empirical_reverse():
    """Verify that undocumented QuantSeq orientation is resolved to reverse when empirical evidence is antisense."""
    metrics = {
        "mean_forward_fraction": 0.04,
        "mean_reverse_fraction": 0.91,
        "mean_unexplained_fraction": 0.05,
    }
    decision = infer_strandedness_decision(metrics, library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND)
    assert decision.strandedness_flag == 2
    assert decision.strand_mode == "reverse"
    assert decision.confidence == "HIGH"
    assert "EMPIRICAL_QUANTSEQ_REVERSE" in decision.consensus_indication


def test_strandedness_undocumented_quantseq_non_stranded_halts():
    """Verify that undocumented QuantSeq orientation halts if empirical evidence is unstranded or ambiguous."""
    metrics = {
        "mean_forward_fraction": 0.45,
        "mean_reverse_fraction": 0.45,
        "mean_unexplained_fraction": 0.10,
    }
    decision = infer_strandedness_decision(metrics, library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND)
    assert decision.strandedness_flag == -1
    assert decision.confidence == "CONFLICT"
    assert "CONFLICT_QUANTSEQ_NON_STRANDED" in decision.consensus_indication


def test_configured_strandedness_cannot_bypass_conflicting_empirical_evidence():
    """Verify that investigator-configured strandedness does not silently bypass conflicting empirical evidence."""
    # Investigator configured -s 1 (forward), but BAMs show 92% reverse antisense
    metrics = {
        "mean_forward_fraction": 0.04,
        "mean_reverse_fraction": 0.92,
        "mean_unexplained_fraction": 0.04,
    }
    decision = infer_strandedness_decision(
        metrics,
        library_type=LibraryType.CONVENTIONAL_BULK,
        configured_strandedness=1,
    )
    assert decision.strandedness_flag == -1
    assert decision.confidence == "CONFLICT"
    assert "CONFLICT_CONFIGURED_FORWARD_EMPIRICAL_REVERSE" in decision.consensus_indication
    assert "cannot silently bypass" in decision.rationale


def test_configured_strandedness_verified_by_agreeing_empirical_evidence():
    """Verify that investigator-configured strandedness succeeds when verified by empirical evidence."""
    metrics = {
        "mean_forward_fraction": 0.92,
        "mean_reverse_fraction": 0.04,
        "mean_unexplained_fraction": 0.04,
    }
    decision = infer_strandedness_decision(
        metrics,
        library_type=LibraryType.CONVENTIONAL_BULK,
        configured_strandedness=1,
    )
    assert decision.strandedness_flag == 1
    assert decision.confidence == "HIGH"
    assert "CONSENSUS_VERIFIED_USER_CONFIGURED" in decision.consensus_indication


# -----------------------------------------------------------------------------
# 4. featureCounts Single-End Flag Construction Tests
# -----------------------------------------------------------------------------

def test_featurecounts_single_end_flags(tmp_path: Path):
    """Verify single-end featureCounts omits -p, --countReadPairs, -B, -C and uses -s 1."""
    executed_cmds = []

    def mock_run_wsl(cmd):
        executed_cmds.append(cmd)
        m = MagicMock()
        m.returncode = 0
        return m

    bam = tmp_path / "sample1.bam"
    gtf = tmp_path / "annotation.gtf"
    out_counts = tmp_path / "counts.txt"
    bam.write_bytes(b"dummy_bam")
    gtf.write_text("dummy_gtf")
    out_counts.write_text("Geneid\tLength\tCounts\nG1\t100\t10\n")

    with patch("src.quantification_bridge.featurecounts_runner.run_wsl", side_effect=mock_run_wsl):
        run_featurecounts(
            bam_paths=[bam],
            gtf_path=gtf,
            output_counts_path=out_counts,
            strandedness=1,
            paired_end=False,
        )

        cmd = executed_cmds[-1]
        assert "-s 1" in cmd
        assert "-p" not in cmd
        assert "--countReadPairs" not in cmd
        assert "-B" not in cmd
        assert "-C" not in cmd


# -----------------------------------------------------------------------------
# 5. UMI Safeguard Tests
# -----------------------------------------------------------------------------

def test_umi_safeguard_fails_closed(tmp_path: Path):
    """Verify that QuantSeq with UMIs fails closed during preflight and execution."""
    cfg = UniversalPipelineConfig(
        study_id="TEST_UMI_STUDY",
        output_dir=tmp_path / "out",
        read_layout="single_end",
        library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI,
        enable_umi_pipeline=False,
    )
    ref = ReferenceConfig(
        species="Homo sapiens",
        genome_fasta=tmp_path / "genome.fa",
        gtf_annotation=tmp_path / "genes.gtf",
    )
    ref.genome_fasta.write_text(">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n")
    ref.gtf_annotation.write_text('chr1\tHAVANA\texon\t1\t30\t.\t+\t.\tgene_id "G1";\n')

    fq = tmp_path / "s1.fq.gz"
    make_dummy_fastq(fq, ["readA", "readB"], mate_suffix="", is_gz=True)

    sample = SampleInput(
        sample_id="SAMPLE_UMI",
        r1_path=fq,
        read_layout="single_end",
        library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI,
    )

    pipeline = UniversalFastqPipeline(
        config=cfg,
        references=ref,
        samples=[sample],
    )

    preflight = pipeline.validate_inputs()
    assert preflight.has_fatal_errors is True
    assert any(f.code == "QUANTSEQ_UMI_PIPELINE_UNIMPLEMENTED" for f in preflight.findings)

    # Preflight gate blocks execution
    plan = pipeline.plan()
    assert plan.can_execute is False
    res_run = pipeline.run()
    assert res_run.status == "BLOCKED_PREFLIGHT_FAILURES"
    assert any(f.code == "QUANTSEQ_UMI_PIPELINE_UNIMPLEMENTED" for f in res_run.findings)


def test_umi_safeguard_cannot_be_bypassed_by_enable_flag(tmp_path: Path):
    """Verify that setting enable_umi_pipeline=True cannot bypass the safety gate."""
    cfg = UniversalPipelineConfig(
        study_id="TEST_UMI_BYPASS",
        output_dir=tmp_path / "out",
        read_layout="single_end",
        library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI,
        enable_umi_pipeline=True,  # Attempting to bypass
    )
    ref = ReferenceConfig(
        species="Homo sapiens",
        genome_fasta=tmp_path / "genome.fa",
        gtf_annotation=tmp_path / "genes.gtf",
    )
    ref.genome_fasta.write_text(">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n")
    ref.gtf_annotation.write_text('chr1\tHAVANA\texon\t1\t30\t.\t+\t.\tgene_id "G1";\n')

    fq = tmp_path / "s1.fq.gz"
    make_dummy_fastq(fq, ["readA", "readB"], is_gz=True)

    sample = SampleInput(
        sample_id="SAMPLE_UMI_BYPASS",
        r1_path=fq,
        read_layout="single_end",
        library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI,
    )

    pipeline = UniversalFastqPipeline(
        config=cfg,
        references=ref,
        samples=[sample],
    )

    preflight = pipeline.validate_inputs()
    assert preflight.has_fatal_errors is True
    assert any(f.code == "QUANTSEQ_UMI_PIPELINE_UNIMPLEMENTED" for f in preflight.findings)

    plan = pipeline.plan()
    assert plan.can_execute is False
    res_run = pipeline.run()
    assert res_run.status == "BLOCKED_PREFLIGHT_FAILURES"
    assert any(f.code == "QUANTSEQ_UMI_PIPELINE_UNIMPLEMENTED" for f in res_run.findings)


# -----------------------------------------------------------------------------
# 6. Cohort Library-Type Consistency Preflight Tests
# -----------------------------------------------------------------------------

def test_mixed_library_types_rejected(tmp_path: Path):
    """Verify that mixed QuantSeq and conventional bulk samples within a cohort trigger preflight failure."""
    cfg = UniversalPipelineConfig(
        study_id="TEST_MIXED",
        output_dir=tmp_path / "out",
        read_layout="single_end",
    )
    ref = ReferenceConfig(
        species="Homo sapiens",
        genome_fasta=tmp_path / "genome.fa",
        gtf_annotation=tmp_path / "genes.gtf",
    )
    ref.genome_fasta.write_text(">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n")
    ref.gtf_annotation.write_text('chr1\tHAVANA\texon\t1\t30\t.\t+\t.\tgene_id "G1";\n')

    fq1 = tmp_path / "s1.fq.gz"
    fq2 = tmp_path / "s2.fq.gz"
    make_dummy_fastq(fq1, ["r1"], is_gz=True)
    make_dummy_fastq(fq2, ["r2"], is_gz=True)

    s1 = SampleInput(
        sample_id="SAMPLE_QUANTSEQ",
        r1_path=fq1,
        read_layout="single_end",
        library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
    )
    s2 = SampleInput(
        sample_id="SAMPLE_BULK",
        r1_path=fq2,
        read_layout="single_end",
        library_type=LibraryType.CONVENTIONAL_BULK,
    )

    pipeline = UniversalFastqPipeline(
        config=cfg,
        references=ref,
        samples=[s1, s2],
    )

    preflight = pipeline.validate_inputs()
    codes = [f.code for f in preflight.findings]
    assert "MIXED_LIBRARY_TYPES_UNSUPPORTED" in codes
    assert preflight.has_fatal_errors is True


def test_study_sample_library_conflict_rejected(tmp_path: Path):
    """Verify that conflict between study-level and sample-level declarations triggers fatal preflight failure."""
    cfg = UniversalPipelineConfig(
        study_id="TEST_CONFLICT",
        output_dir=tmp_path / "out",
        read_layout="single_end",
        library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
    )
    ref = ReferenceConfig(
        species="Homo sapiens",
        genome_fasta=tmp_path / "genome.fa",
        gtf_annotation=tmp_path / "genes.gtf",
    )
    ref.genome_fasta.write_text(">chr1\nACGT\n")
    ref.gtf_annotation.write_text('chr1\tHAVANA\texon\t1\t4\t.\t+\t.\tgene_id "G1";\n')

    fq = tmp_path / "s1.fq.gz"
    make_dummy_fastq(fq, ["r1"], is_gz=True)

    # Sample explicitly declares REV while study declares FWD
    s = SampleInput(
        sample_id="SAMPLE_REV_IN_FWD_STUDY",
        r1_path=fq,
        read_layout="single_end",
        library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV,
    )

    pipeline = UniversalFastqPipeline(config=cfg, references=ref, samples=[s])
    preflight = pipeline.validate_inputs()
    codes = [f.code for f in preflight.findings]
    assert "STUDY_SAMPLE_LIBRARY_TYPE_CONFLICT" in codes
    assert preflight.has_fatal_errors is True


def test_heterogeneous_cohort_library_types_rejected(tmp_path: Path):
    """Verify that heterogeneous QuantSeq types (FWD + REV) in a cohort trigger preflight failure."""
    cfg = UniversalPipelineConfig(
        study_id="TEST_HETEROGENEOUS",
        output_dir=tmp_path / "out",
        read_layout="single_end",
    )
    ref = ReferenceConfig(
        species="Homo sapiens",
        genome_fasta=tmp_path / "genome.fa",
        gtf_annotation=tmp_path / "genes.gtf",
    )
    ref.genome_fasta.write_text(">chr1\nACGT\n")
    ref.gtf_annotation.write_text('chr1\tHAVANA\texon\t1\t4\t.\t+\t.\tgene_id "G1";\n')

    fq1 = tmp_path / "s1.fq.gz"
    fq2 = tmp_path / "s2.fq.gz"
    make_dummy_fastq(fq1, ["r1"], is_gz=True)
    make_dummy_fastq(fq2, ["r2"], is_gz=True)

    s1 = SampleInput(sample_id="S1_FWD", r1_path=fq1, read_layout="single_end", library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD)
    s2 = SampleInput(sample_id="S2_REV", r1_path=fq2, read_layout="single_end", library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV)

    pipeline = UniversalFastqPipeline(config=cfg, references=ref, samples=[s1, s2])
    preflight = pipeline.validate_inputs()
    codes = [f.code for f in preflight.findings]
    assert "HETEROGENEOUS_COHORT_LIBRARY_TYPES" in codes
    assert preflight.has_fatal_errors is True


def test_configured_strandedness_conflicts_declared_library_rejected(tmp_path: Path):
    """Verify that configured strandedness (-s 2) contradicting declared library (FWD -> -s 1) triggers preflight failure."""
    cfg = UniversalPipelineConfig(
        study_id="TEST_STRAND_CONFLICT",
        output_dir=tmp_path / "out",
        read_layout="single_end",
        library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
        strandedness=2,  # Contradicts FWD (-s 1)
    )
    ref = ReferenceConfig(
        species="Homo sapiens",
        genome_fasta=tmp_path / "genome.fa",
        gtf_annotation=tmp_path / "genes.gtf",
    )
    ref.genome_fasta.write_text(">chr1\nACGT\n")
    ref.gtf_annotation.write_text('chr1\tHAVANA\texon\t1\t4\t.\t+\t.\tgene_id "G1";\n')

    fq = tmp_path / "s1.fq.gz"
    make_dummy_fastq(fq, ["r1"], is_gz=True)
    s = SampleInput(sample_id="S1", r1_path=fq, read_layout="single_end")

    pipeline = UniversalFastqPipeline(config=cfg, references=ref, samples=[s])
    preflight = pipeline.validate_inputs()
    codes = [f.code for f in preflight.findings]
    assert "CONFIGURED_STRANDEDNESS_CONFLICT_WITH_DECLARED_LIBRARY" in codes
    assert preflight.has_fatal_errors is True


# -----------------------------------------------------------------------------
# 7. Metadata Ingestion Library Attribute Detection Test
# -----------------------------------------------------------------------------

def test_metadata_library_detection_from_attributes(tmp_path: Path):
    """Verify that SampleInput correctly resolves QuantSeq from metadata attributes."""
    s = SampleInput(
        sample_id="SRR32575599",
        read_layout="single_end",
        attributes={
            "library_kit": "3' mRNA-Seq Library Prep Kit FWD for Illumina (QuantSeq-LEXOGEN)",
            "organism": "Homo sapiens",
        },
    )
    cfg = UniversalPipelineConfig(study_id="SRA_STUDY", output_dir=tmp_path / "out", read_layout="single_end")
    ref = ReferenceConfig(species="Homo sapiens", genome_fasta=tmp_path / "fa", gtf_annotation=tmp_path / "gtf")
    ref.genome_fasta.write_text(">chr1\nACGT\n")
    ref.gtf_annotation.write_text('chr1\tHAVANA\texon\t1\t4\t.\t+\t.\tgene_id "G1";\n')

    pipeline = UniversalFastqPipeline(config=cfg, references=ref, samples=[s])
    resolved_lib = pipeline.get_sample_library_type(s)
    assert resolved_lib == LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD


# -----------------------------------------------------------------------------
# 8. End-to-End QuantSeq FWD Execution Integration Test
# -----------------------------------------------------------------------------

def test_08_end_to_end_quantseq_fwd_pipeline_execution(tmp_path: Path, monkeypatch):
    """Verify complete QuantSeq FWD pipeline execution flow through all stages."""
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n")
    gtf.write_text('chr1\tHAVANA\texon\t1\t30\t.\t+\t.\tgene_id "ENSG000001";\n')

    fq1 = tmp_path / "QS_CTRL_01.fastq.gz"
    fq2 = tmp_path / "QS_TRT_01.fastq.gz"
    make_dummy_fastq(fq1, ["read1", "read2"], is_gz=True)
    make_dummy_fastq(fq2, ["read3", "read4"], is_gz=True)

    samples = [
        SampleInput(
            sample_id="QS_CTRL_01",
            r1_path=fq1,
            read_layout="single_end",
            library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
            attributes={"condition_group": "control"},
        ),
        SampleInput(
            sample_id="QS_TRT_01",
            r1_path=fq2,
            read_layout="single_end",
            library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
            attributes={"condition_group": "treatment"},
        ),
    ]

    cfg = UniversalPipelineConfig(
        study_id="STUDY_QUANTSEQ_E2E",
        work_dir=tmp_path / "work_qs",
        read_layout="single_end",
        library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
        enforce_resources=False,
    )
    refs = ReferenceConfig(
        species="Homo sapiens",
        genome_fasta=fasta,
        gtf_annotation=gtf,
    )

    # Track tool calls
    fastp_calls = []
    def mock_fastp(raw_r1, raw_r2=None, out_r1=None, out_r2=None, json_report=None, html_report=None, sample_id="sample", threads=8, **kwargs):
        fastp_calls.append({"raw_r1": raw_r1, "kwargs": kwargs, "sample_id": sample_id})
        out1 = Path(out_r1)
        out1.parent.mkdir(parents=True, exist_ok=True)
        out1.write_bytes(Path(raw_r1).read_bytes())
        if json_report:
            Path(json_report).write_text(json.dumps({"summary": {"before_filtering": {"total_reads": 2}, "after_filtering": {"total_reads": 2}}}), encoding="utf-8")
        return {"status": "FASTP_PREPROCESSING_SUCCESS", "trimmed_r1": str(out1), "metrics_summary": {"total_reads": 2}}

    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.run_fastp_preprocessing", mock_fastp)

    def mock_build_star(fasta_path, gtf_path, index_dir, **kwargs):
        idx_p = Path(index_dir)
        idx_p.mkdir(parents=True, exist_ok=True)
        (idx_p / "Genome").write_bytes(b"mock genome")
        (idx_p / "SA").write_bytes(b"mock sa")

    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.build_star_index", mock_build_star)

    def mock_align(trimmed_r1, output_prefix, **kwargs):
        pref = Path(output_prefix)
        pref.parent.mkdir(parents=True, exist_ok=True)
        bam = pref.parent / f"{pref.name}Aligned.sortedByCoord.out.bam"
        bam.write_bytes(b"mock_bam")
        log = pref.parent / f"{pref.name}Log.final.out"
        log.write_text("Uniquely mapped reads % | 92.0%\n", encoding="utf-8")
        return {"bam_path": bam, "log_final_path": log, "metrics": {"uniquely_mapped_percent": 92.0}}

    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.align_reads_star", mock_align)
    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.run_wsl", lambda cmd: MagicMock(returncode=0, stdout="", stderr=""))

    def mock_infer_exp(bam_path, bed_path, **kwargs):
        # Single-end QuantSeq FWD produces strong sense (++,--) signal
        return {
            "status": "INFER_EXPERIMENT_SUCCESS",
            "forward_fraction": 0.94,
            "reverse_fraction": 0.03,
            "unexplained_fraction": 0.03,
        }

    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.run_infer_experiment", mock_infer_exp)

    fc_calls = []
    def mock_fc(bam_paths, gtf_path, output_counts_path, strandedness=0, paired_end=True, **kwargs):
        fc_calls.append({"strandedness": strandedness, "paired_end": paired_end, "bams": bam_paths})
        out_f = Path(output_counts_path)
        out_f.parent.mkdir(parents=True, exist_ok=True)
        header = "Geneid\tChr\tStart\tEnd\tStrand\tLength\t" + "\t".join([Path(b).name for b in bam_paths]) + "\n"
        row = "ENSG000001\tchr1\t1\t30\t+\t30\t" + "\t".join(["25" for _ in bam_paths]) + "\n"
        out_f.write_text(header + row, encoding="utf-8")
        summary_f = Path(f"{output_counts_path}.summary")
        summary_header = "Status\t" + "\t".join([Path(b).name for b in bam_paths]) + "\n"
        summary_row = "Assigned\t" + "\t".join(["25" for _ in bam_paths]) + "\n"
        summary_f.write_text(summary_header + summary_row, encoding="utf-8")

    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.run_featurecounts", mock_fc)

    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)
    plan = pipeline.plan()
    assert plan.can_execute is True

    pipeline.approve(plan.approval_token)
    res = pipeline.run(approval_token=plan.approval_token)

    # 1. Pipeline successfully proceeded past strandedness and featurecounts
    assert "fastp_QS_CTRL_01" in res.completed_stages
    assert "fastp_QS_TRT_01" in res.completed_stages
    assert "star_align_QS_CTRL_01" in res.completed_stages
    assert "star_align_QS_TRT_01" in res.completed_stages
    assert "rseqc_strandedness" in res.completed_stages
    assert "featurecounts" in res.completed_stages

    # 2. Check fastp calls: QuantSeq library type passed
    for call in fastp_calls:
        assert call["kwargs"]["library_type"] == LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD

    # 3. Check featureCounts call: -s 1 (forward) and paired_end=False
    assert len(fc_calls) == 1
    assert fc_calls[0]["strandedness"] == 1
    assert fc_calls[0]["paired_end"] is False

    # 4. Check strandedness evidence artifact
    assert res.strandedness_evidence_path.is_file()
    strand_data = json.loads(res.strandedness_evidence_path.read_text(encoding="utf-8"))
    assert strand_data["library_type"] == "LEXOGEN_QUANTSEQ_3UTR_FWD"
    assert strand_data["decision"]["strandedness_flag"] == 1
    assert strand_data["decision"]["strand_mode"] == "forward"

    # 5. Check execution record artifact
    assert res.execution_record_path.is_file()
    exec_data = json.loads(res.execution_record_path.read_text(encoding="utf-8"))
    assert exec_data["library_type"] == "LEXOGEN_QUANTSEQ_3UTR_FWD"
    assert exec_data["read_layout"] == "single_end"
    assert exec_data["strandedness_flag"] == 1


# -----------------------------------------------------------------------------
# 9. Genuine Toolchain QuantSeq FWD Mini-Genome Execution
# -----------------------------------------------------------------------------

def test_09_genuine_quantseq_fwd_mini_genome_end_to_end(tmp_path: Path):
    """Execute genuine non-mocked tools (fastp -> STAR -> RSeQC -> featureCounts) on single-end QuantSeq FWD reads."""
    import random
    random.seed(42)

    work_dir = tmp_path / "genuine_qs_work"

    # 1. 100k bp synthetic reference genome
    seq_100k = "".join(random.choices(["A", "C", "G", "T"], k=100000))
    fasta = tmp_path / "mini_genome.fa"
    fasta.write_text(f">chr1\n{seq_100k}\n", encoding="utf-8")

    # 2. GTF annotation: 2 forward-strand exons
    gtf = tmp_path / "mini_genes.gtf"
    gtf_content = (
        'chr1\tsynth\texon\t1000\t5000\t.\t+\t.\tgene_id "ENSG00000000001"; gene_name "GENE1"; transcript_id "ENST00000000001";\n'
        'chr1\tsynth\texon\t10000\t15000\t.\t+\t.\tgene_id "ENSG00000000002"; gene_name "GENE2"; transcript_id "ENST00000000002";\n'
    )
    gtf.write_text(gtf_content, encoding="utf-8")

    # 3. Construct single-end QuantSeq 3' FWD reads with 3' poly(A) tails and Illumina adapter
    # Gene 1 sense sequence + poly(A) tail + Lexogen adapter
    g1_frag = seq_100k[4900:4950]  # 50 bp 3' end of gene 1
    poly_a = "A" * 15
    adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    qs_read_seq = (g1_frag + poly_a + adapter)[:76]  # 76 bp NextSeq read

    r1 = tmp_path / "qs_sample_R1.fastq.gz"
    with gzip.open(r1, "wt", encoding="utf-8") as f:
        for i in range(50):
            qual = "I" * len(qs_read_seq)
            f.write(f"@read_qs_{i}\n{qs_read_seq}\n+\n{qual}\n")

    cfg = UniversalPipelineConfig(
        study_id="STUDY_GENUINE_QS",
        work_dir=work_dir,
        read_layout="single_end",
        library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
        strandedness=1,  # explicitly configured forward strandedness for 3' FWD
        threads=2,
        require_approval=True,
        enforce_resources=False,
    )
    refs = ReferenceConfig(
        species="Homo sapiens",
        genome_fasta=fasta,
        gtf_annotation=gtf,
        sjdb_overhang=49,
    )
    samples = [
        SampleInput(
            sample_id="GENUINE_QS_SAMPLE",
            r1_path=r1,
            read_layout="single_end",
            library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
            attributes={"condition_group": "control"},
        )
    ]

    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)
    plan = pipeline.plan()
    assert plan.can_execute is True

    # Approve and execute genuine pipeline
    pipeline.approve(plan.approval_token)
    res = pipeline.run(approval_token=plan.approval_token)

    # 1. Verify fastp trimmed the poly(A) tail and adapter
    fastp_json = work_dir / "fastp_qc" / "GENUINE_QS_SAMPLE_fastp.json"
    assert fastp_json.is_file()
    fastp_metrics = json.loads(fastp_json.read_text(encoding="utf-8"))
    assert fastp_metrics["summary"]["before_filtering"]["total_reads"] == 50
    # Trimmed reads must be <= 50 bp after removing 15 bp polyA + adapter
    assert fastp_metrics["summary"]["after_filtering"]["read1_mean_length"] < 76

    # 2. Verify STAR alignment output
    bam = work_dir / "alignments" / "GENUINE_QS_SAMPLE_Aligned.sortedByCoord.out.bam"
    assert bam.is_file()
    assert bam.stat().st_size > 0

    # 3. Verify featureCounts quantification output
    counts_file = work_dir / "counts" / "counts.tsv"
    assert counts_file.is_file()
    counts_text = counts_file.read_text(encoding="utf-8")
    assert "ENSG00000000001" in counts_text
    assert "GENUINE_QS_SAMPLE" in counts_text

    # 4. Verify featureCounts summary shows assigned reads
    fc_summary = work_dir / "counts" / "featurecounts_raw.txt.summary"
    assert fc_summary.is_file()
    summary_text = fc_summary.read_text(encoding="utf-8")
    assert "Assigned" in summary_text

    # 5. Verify execution record captures QuantSeq details
    assert res.execution_record_path.is_file()
    exec_record = json.loads(res.execution_record_path.read_text(encoding="utf-8"))
    assert exec_record["library_type"] == "LEXOGEN_QUANTSEQ_3UTR_FWD"
    assert exec_record["read_layout"] == "single_end"
    assert exec_record["strandedness_flag"] == 1


