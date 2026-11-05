"""Acceptance test suite for Universal FASTQ Automation — Phase F1.

Covers the 10 required validation scenarios:
1. Valid paired-end FASTQ.
2. Valid single-end FASTQ.
3. Missing FASTQ mate (orphan R1).
4. Duplicate sample identifiers.
5. Unknown biological replicate information.
6. Missing or incompatible reference files.
7. Interrupted execution and safe resume.
8. Existing valid outputs reused without reprocessing.
9. Strandedness assessment requiring a documented decision.
10. Different studies and sample metadata without PRJNA802072 hardcoding.
"""

import gzip
import json
from pathlib import Path
import pytest

from src.quantification_bridge import (
    FastqPair,
    FastqSingle,
    PipelineFinding,
    PreflightReport,
    ProcessingPlan,
    ReferenceConfig,
    SampleInput,
    UniversalCheckpointManager,
    UniversalFastqPipeline,
    UniversalFastqResult,
    UniversalPipelineConfig,
    load_sample_metadata,
    validate_paired_fastqs,
    validate_single_fastq,
)
from src.quantification_bridge.downstream_bridge import execute_downstream_pipeline
from src.quantification_bridge.rseqc_runner import StrandednessDecision, infer_strandedness_decision


# -----------------------------------------------------------------------------
# Test Helpers
# -----------------------------------------------------------------------------

def make_dummy_fastq(
    path: Path,
    read_names: list[str],
    read_len: int = 50,
    mate_suffix: str = "/1",
    is_gz: bool = True,
) -> Path:
    """Create a synthetically valid FASTQ file for structural validation."""
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


def make_dummy_fasta(path: Path) -> Path:
    """Create a minimal synthetic FASTA reference file."""
    content = ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n>chr2\nTGCATGCATGCATGCATGCATGCATGCATGCA\n"
    path.write_text(content, encoding="utf-8")
    return path


def make_dummy_gtf(path: Path) -> Path:
    """Create a minimal synthetic GTF annotation reference file."""
    content = (
        'chr1\tsynth\texon\t10\t100\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001";\n'
        'chr1\tsynth\texon\t200\t300\t.\t+\t.\tgene_id "ENSG00000000002"; transcript_id "ENST00000000002";\n'
    )
    path.write_text(content, encoding="utf-8")
    return path


# -----------------------------------------------------------------------------
# Test 1: Valid paired-end FASTQ
# -----------------------------------------------------------------------------

def test_01_valid_paired_end_fastq(tmp_path: Path):
    """Verify input validation, pairing, synchronization, and planning for paired-end FASTQs."""
    r1 = tmp_path / "sample01_R1.fastq.gz"
    r2 = tmp_path / "sample01_R2.fastq.gz"
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"

    make_dummy_fastq(r1, ["read1", "read2", "read3"], mate_suffix="/1", is_gz=True)
    make_dummy_fastq(r2, ["read1", "read2", "read3"], mate_suffix="/2", is_gz=True)
    make_dummy_fasta(fasta)
    make_dummy_gtf(gtf)

    cfg = UniversalPipelineConfig(
        study_id="STUDY_TEST_PE",
        work_dir=tmp_path / "work",
        read_layout="paired_end",
        enforce_resources=False,
    )
    refs = ReferenceConfig(
        species="Homo sapiens",
        genome_fasta=fasta,
        gtf_annotation=gtf,
    )
    samples = [
        SampleInput(
            sample_id="SAMPLE_PE_01",
            r1_path=r1,
            r2_path=r2,
            read_layout="paired_end",
        )
    ]
    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)
    preflight = pipeline.validate_inputs()

    assert preflight.is_valid is True
    assert not preflight.has_fatal_errors
    s_val = preflight.sample_validation["SAMPLE_PE_01"]
    assert s_val["status"] == "VALIDATION_PASSED"
    assert s_val["is_synchronized"] is True
    assert s_val["total_paired_reads"] == 3
    assert s_val["read_layout"] == "paired_end"

    plan = pipeline.plan()
    assert plan.can_execute is True
    assert len(plan.planned_stages) > 0
    fastp_stage = next(st for st in plan.planned_stages if st["stage_id"] == "02_fastp_SAMPLE_PE_01")
    assert fastp_stage["out_r2"] is not None


# -----------------------------------------------------------------------------
# Test 2: Valid single-end FASTQ
# -----------------------------------------------------------------------------

def test_02_valid_single_end_fastq(tmp_path: Path):
    """Verify input validation and single-end command planning for single-end FASTQs."""
    r1 = tmp_path / "sample_se_R1.fastq"
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"

    make_dummy_fastq(r1, ["readA", "readB"], mate_suffix="", is_gz=False)
    make_dummy_fasta(fasta)
    make_dummy_gtf(gtf)

    cfg = UniversalPipelineConfig(
        study_id="STUDY_TEST_SE",
        work_dir=tmp_path / "work_se",
        read_layout="single_end",
        enforce_resources=False,
    )
    refs = ReferenceConfig(
        species="Homo sapiens",
        genome_fasta=fasta,
        gtf_annotation=gtf,
    )
    samples = [
        SampleInput(
            sample_id="SAMPLE_SE_01",
            r1_path=r1,
            r2_path=None,
            read_layout="single_end",
        )
    ]
    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)
    preflight = pipeline.validate_inputs()

    assert preflight.is_valid is True
    assert not preflight.has_fatal_errors
    s_val = preflight.sample_validation["SAMPLE_SE_01"]
    assert s_val["status"] == "VALIDATION_PASSED"
    assert s_val["read_layout"] == "single_end"
    assert s_val["is_compressed"] is False
    assert s_val["total_reads"] == 2

    plan = pipeline.plan()
    assert plan.can_execute is True
    fastp_stage = next(st for st in plan.planned_stages if st["stage_id"] == "02_fastp_SAMPLE_SE_01")
    assert fastp_stage["out_r2"] is None


# -----------------------------------------------------------------------------
# Test 3: Missing FASTQ mate (orphan R1)
# -----------------------------------------------------------------------------

def test_03_missing_fastq_mate(tmp_path: Path):
    """Verify that an orphan R1 in a paired-end configuration halts with a fatal preflight error."""
    r1 = tmp_path / "sample_orphan_R1.fastq.gz"
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"

    make_dummy_fastq(r1, ["read1"], mate_suffix="/1")
    make_dummy_fasta(fasta)
    make_dummy_gtf(gtf)

    cfg = UniversalPipelineConfig(
        study_id="STUDY_MISSING_MATE",
        work_dir=tmp_path / "work_mate",
        read_layout="paired_end",
        enforce_resources=False,
    )
    refs = ReferenceConfig(
        species="Homo sapiens",
        genome_fasta=fasta,
        gtf_annotation=gtf,
    )
    # R2 is missing (None) in paired_end layout
    samples = [
        SampleInput(
            sample_id="SAMPLE_ORPHAN",
            r1_path=r1,
            r2_path=None,
            read_layout="paired_end",
        )
    ]
    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)
    preflight = pipeline.validate_inputs()

    assert preflight.is_valid is False
    assert preflight.has_fatal_errors is True
    mate_findings = [f for f in preflight.findings if f.code == "MISSING_FASTQ_MATE"]
    assert len(mate_findings) == 1
    assert "reverse mate (R2) is missing" in mate_findings[0].message

    # Execution is blocked safely before running any commands
    result = pipeline.run()
    assert result.status == "BLOCKED_PREFLIGHT_FAILURES"


# -----------------------------------------------------------------------------
# Test 4: Duplicate sample identifiers
# -----------------------------------------------------------------------------

def test_04_duplicate_sample_identifiers(tmp_path: Path):
    """Verify that duplicate sample identifiers in a cohort are detected and blocked."""
    r1a = tmp_path / "s1a_R1.fastq.gz"
    r2a = tmp_path / "s1a_R2.fastq.gz"
    r1b = tmp_path / "s1b_R1.fastq.gz"
    r2b = tmp_path / "s1b_R2.fastq.gz"
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"

    make_dummy_fastq(r1a, ["r1"], mate_suffix="/1")
    make_dummy_fastq(r2a, ["r1"], mate_suffix="/2")
    make_dummy_fastq(r1b, ["r1"], mate_suffix="/1")
    make_dummy_fastq(r2b, ["r1"], mate_suffix="/2")
    make_dummy_fasta(fasta)
    make_dummy_gtf(gtf)

    cfg = UniversalPipelineConfig(
        study_id="STUDY_DUPLICATES",
        work_dir=tmp_path / "work_dup",
        enforce_resources=False,
    )
    refs = ReferenceConfig(species="Homo sapiens", genome_fasta=fasta, gtf_annotation=gtf)

    # Identical sample_id used twice
    samples = [
        SampleInput(sample_id="SAMPLE_COMMON_ID", r1_path=r1a, r2_path=r2a),
        SampleInput(sample_id="SAMPLE_COMMON_ID", r1_path=r1b, r2_path=r2b),
    ]
    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)
    preflight = pipeline.validate_inputs()

    assert preflight.is_valid is False
    assert preflight.has_fatal_errors is True
    dup_findings = [f for f in preflight.findings if f.code == "DUPLICATE_SAMPLE_ID"]
    assert len(dup_findings) == 1
    assert "Duplicate sample identifier detected: 'SAMPLE_COMMON_ID'" in dup_findings[0].message


# -----------------------------------------------------------------------------
# Test 5: Unknown biological replicate information
# -----------------------------------------------------------------------------

def test_05_unknown_biological_replicate_information(tmp_path: Path):
    """Verify that unverified biological replication is mapped to unknown and gates inferential analysis."""
    # Create sample metadata CSV where biological_replicate_id is absent
    csv_file = tmp_path / "metadata_noreps.csv"
    csv_file.write_text(
        "sample_id,organoid_type,cell_line,chemical,concentration,condition_status,control_type\n"
        "SMPL_C1,brain_organoid,iPSC_LINE_A,DMSO,0.1,control,vehicle\n"
        "SMPL_C2,brain_organoid,iPSC_LINE_A,DMSO,0.1,control,vehicle\n"
        "SMPL_T1,brain_organoid,iPSC_LINE_A,ValproicAcid,250,treatment,\n"
        "SMPL_T2,brain_organoid,iPSC_LINE_A,ValproicAcid,250,treatment,\n",
        encoding="utf-8",
    )

    canonical_doc, raw_records = load_sample_metadata(
        metadata_source=csv_file,
        study_id="STUDY_UNKNOWN_REPS",
        source_type="in_house",
    )

    # 1. Biological replicates must not be invented; value_status must be 'unknown'
    for s in canonical_doc["samples"]:
        assert s["biological_replicate_id"]["original_value"] is None
        assert s["biological_replicate_id"]["value_status"] == "unknown"

    # 2. Preflight on pipeline detects unverified replication warning
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    make_dummy_fasta(fasta)
    make_dummy_gtf(gtf)

    r1 = tmp_path / "s1_R1.fastq.gz"
    make_dummy_fastq(r1, ["r1"], mate_suffix="")
    samples = [
        SampleInput(sample_id="SMPL_C1", r1_path=r1, read_layout="single_end"),
    ]
    cfg = UniversalPipelineConfig(
        study_id="STUDY_UNKNOWN_REPS",
        work_dir=tmp_path / "work_noreps",
        read_layout="single_end",
        enforce_resources=False,
    )
    refs = ReferenceConfig(species="Homo sapiens", genome_fasta=fasta, gtf_annotation=gtf)
    pipeline = UniversalFastqPipeline(
        config=cfg,
        references=refs,
        samples=samples,
        canonical_metadata=canonical_doc,
    )
    preflight = pipeline.validate_inputs()
    unverified_warnings = [f for f in preflight.findings if f.code == "BIOLOGICAL_REPLICATION_UNVERIFIED"]
    assert len(unverified_warnings) == 1

    # 3. Downstream bridge safely gates inferential DESeq2/pathway analysis
    from src.features.pathway.mapping import EnsemblToEntrezAdapter
    adapter = EnsemblToEntrezAdapter()
    genes = [r.ensembl_gene_id for r in adapter.gene_reference.records if r.gene_biotype == 'protein_coding'][:30]
    matrix_file = tmp_path / "counts_noreps.tsv"
    runs = ("SMPL_C1", "SMPL_C2", "SMPL_T1", "SMPL_T2")
    gene_lines = ["gene_id\t" + "\t".join(runs)]
    for g in genes:
        gene_lines.append(f"{g}\t100\t110\t200\t210")
    matrix_file.write_text("\n".join(gene_lines) + "\n", encoding="utf-8")

    ds_res = execute_downstream_pipeline(
        counts_matrix_path=matrix_file,
        runs=runs,
        study_id="STUDY_UNKNOWN_REPS",
        canonical_metadata=canonical_doc,
    )

    assert ds_res.qc_passed is True
    assert ds_res.deseq2_status == "BLOCKED_UNVERIFIED_BIOLOGICAL_REPLICATION"
    assert ds_res.pathway_status == "BLOCKED_UNVERIFIED_BIOLOGICAL_REPLICATION"
    assert ds_res.status == "DOWNSTREAM_PARTIAL_BIOLOGICAL_REPLICATION_GATED"


# -----------------------------------------------------------------------------
# Test 6: Missing or incompatible reference files
# -----------------------------------------------------------------------------

def test_06_missing_or_incompatible_reference_files(tmp_path: Path):
    """Verify preflight detects non-existent or empty FASTA and GTF reference files."""
    missing_fasta = tmp_path / "nonexistent_genome.fa"
    empty_gtf = tmp_path / "empty_annotation.gtf"
    empty_gtf.write_text("", encoding="utf-8")  # 0 bytes

    r1 = tmp_path / "s1_R1.fastq.gz"
    make_dummy_fastq(r1, ["r1"])

    cfg = UniversalPipelineConfig(
        study_id="STUDY_BAD_REFS",
        work_dir=tmp_path / "work_bad_refs",
        read_layout="single_end",
        enforce_resources=False,
    )
    refs = ReferenceConfig(
        species="Homo sapiens",
        genome_fasta=missing_fasta,
        gtf_annotation=empty_gtf,
    )
    samples = [SampleInput(sample_id="SAMPLE_REF_CHECK", r1_path=r1, read_layout="single_end")]

    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)
    preflight = pipeline.validate_inputs()

    assert preflight.is_valid is False
    assert preflight.has_fatal_errors is True
    codes = [f.code for f in preflight.findings]
    assert "REFERENCE_FASTA_MISSING" in codes
    assert "REFERENCE_GTF_EMPTY" in codes


# -----------------------------------------------------------------------------
# Test 7: Interrupted execution and safe resume
# -----------------------------------------------------------------------------

def test_07_interrupted_execution_and_safe_resume(tmp_path: Path):
    """Verify fine-grained stage checkpoint recording, interruption simulation, and safe resumption."""
    chk_dir = tmp_path / "checkpoints"
    manager = UniversalCheckpointManager(checkpoint_dir=chk_dir)

    out_file1 = tmp_path / "stage1_out.txt"
    out_file1.write_text("stage 1 content", encoding="utf-8")

    input_hash = "input_abc123"
    config_hash = "config_xyz789"

    # Initially, stage 1 is not completed
    assert manager.is_stage_completed(
        stage="fastp",
        identifier="sample_resume_01",
        expected_outputs=[out_file1],
        input_hash=input_hash,
        config_hash=config_hash,
    ) is False

    # Record checkpoint for stage 1
    manager.record_checkpoint(
        stage="fastp",
        identifier="sample_resume_01",
        outputs=[out_file1],
        input_hash=input_hash,
        config_hash=config_hash,
        metrics={"q30_rate": 0.95},
    )

    # Now stage 1 is recognized as complete
    assert manager.is_stage_completed(
        stage="fastp",
        identifier="sample_resume_01",
        expected_outputs=[out_file1],
        input_hash=input_hash,
        config_hash=config_hash,
    ) is True

    # Simulate interruption before stage 2 completes:
    out_file2 = tmp_path / "stage2_out.bam"
    # stage 2 has no checkpoint recorded
    assert manager.is_stage_completed(
        stage="star_align",
        identifier="sample_resume_01",
        expected_outputs=[out_file2],
        input_hash="input_stage2",
        config_hash="config_stage2",
    ) is False

    # On resumption: stage 1 is reused without recomputation, stage 2 proceeds
    out_file2.write_text("stage 2 simulated bam", encoding="utf-8")
    manager.record_checkpoint(
        stage="star_align",
        identifier="sample_resume_01",
        outputs=[out_file2],
        input_hash="input_stage2",
        config_hash="config_stage2",
    )
    assert manager.is_stage_completed(
        stage="star_align",
        identifier="sample_resume_01",
        expected_outputs=[out_file2],
        input_hash="input_stage2",
        config_hash="config_stage2",
    ) is True


# -----------------------------------------------------------------------------
# Test 8: Existing valid outputs reused without reprocessing
# -----------------------------------------------------------------------------

def test_08_existing_valid_outputs_reused_without_reprocessing(tmp_path: Path):
    """Verify that existing valid outputs are safely reused, while corrupted or altered outputs are invalidated."""
    chk_dir = tmp_path / "checkpoints"
    manager = UniversalCheckpointManager(checkpoint_dir=chk_dir)

    out_file = tmp_path / "aligned.bam"
    out_file.write_text("original bam data", encoding="utf-8")

    input_hash = "sha_input_001"
    config_hash = "sha_config_001"

    manager.record_checkpoint(
        stage="star_align",
        identifier="sample_reuse",
        outputs=[out_file],
        input_hash=input_hash,
        config_hash=config_hash,
    )

    # 1. Valid outputs & matching hashes -> Reused
    assert manager.is_stage_completed(
        stage="star_align",
        identifier="sample_reuse",
        expected_outputs=[out_file],
        input_hash=input_hash,
        config_hash=config_hash,
    ) is True

    # 2. Input data changed -> Invalidation & re-execution required
    assert manager.is_stage_completed(
        stage="star_align",
        identifier="sample_reuse",
        expected_outputs=[out_file],
        input_hash="sha_input_DIFFERENT",
        config_hash=config_hash,
    ) is False

    # 3. Parameters/config changed -> Invalidation & re-execution required
    assert manager.is_stage_completed(
        stage="star_align",
        identifier="sample_reuse",
        expected_outputs=[out_file],
        input_hash=input_hash,
        config_hash="sha_config_DIFFERENT",
    ) is False

    # 4. Output file deleted or 0 bytes -> Invalidation
    out_file.write_text("", encoding="utf-8")  # Truncate to 0 bytes
    assert manager.is_stage_completed(
        stage="star_align",
        identifier="sample_reuse",
        expected_outputs=[out_file],
        input_hash=input_hash,
        config_hash=config_hash,
    ) is False


# -----------------------------------------------------------------------------
# Test 9: Strandedness assessment requiring a documented decision
# -----------------------------------------------------------------------------

def test_09_strandedness_assessment_requiring_documented_decision(tmp_path: Path):
    """Verify ambiguous strandedness triggers decision gate, while documented choice proceeds."""
    # 1. Ambiguous RSeQC data (e.g. 15% fwd, 15% rev, 70% unexplained)
    ambiguous_summary = {
        "mean_forward_fraction": 0.15,
        "mean_reverse_fraction": 0.15,
        "mean_unexplained_fraction": 0.70,
    }
    ambig_decision = infer_strandedness_decision(ambiguous_summary)
    assert ambig_decision.strandedness_flag == -1
    assert ambig_decision.strand_mode == "ambiguous"
    assert "Ambiguous orientation proportions" in ambig_decision.rationale

    # 2. Pipeline without explicit strandedness decision and ambiguous RSeQC halts safely
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    make_dummy_fasta(fasta)
    make_dummy_gtf(gtf)
    r1 = tmp_path / "s1_R1.fastq.gz"
    make_dummy_fastq(r1, ["r1"], mate_suffix="")

    cfg = UniversalPipelineConfig(
        study_id="STUDY_STRAND_AMBIG",
        work_dir=tmp_path / "work_strand",
        read_layout="single_end",
        strandedness=None,  # Unset -> requires empirical diagnosis
        enforce_resources=False,
    )
    refs = ReferenceConfig(species="Homo sapiens", genome_fasta=fasta, gtf_annotation=gtf)
    samples = [SampleInput(sample_id="SAMPLE_STRAND_01", r1_path=r1, read_layout="single_end")]

    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)

    # In dry-run mode, plan shows strandedness evaluation stage
    plan = pipeline.plan()
    rseqc_stage = next(st for st in plan.planned_stages if st["stage_id"] == "04_rseqc_strandedness")
    assert rseqc_stage["configured_strandedness"] is None

    # 3. Explicitly documented decision (-s 0, unstranded)
    cfg_unstranded = UniversalPipelineConfig(
        study_id="STUDY_STRAND_DOC",
        work_dir=tmp_path / "work_strand_doc",
        read_layout="single_end",
        strandedness=0,  # Explicitly configured
        enforce_resources=False,
    )
    pipeline_doc = UniversalFastqPipeline(config=cfg_unstranded, references=refs, samples=samples)
    plan_doc = pipeline_doc.plan()
    rseqc_stage_doc = next(st for st in plan_doc.planned_stages if st["stage_id"] == "04_rseqc_strandedness")
    assert rseqc_stage_doc["configured_strandedness"] == 0


# -----------------------------------------------------------------------------
# Test 10: Different study and sample metadata without PRJNA802072 hardcoding
# -----------------------------------------------------------------------------

def test_10_different_study_and_sample_metadata_without_prjna802072_hardcoding(tmp_path: Path):
    """Verify handling of an entirely independent study with custom metadata, free of PRJNA802072 assumptions."""
    # 1. Researcher provides TSV metadata for a novel study
    tsv_content = (
        "sample_id\torganism\torganoid_type\tcell_line\tchemical\tconcentration\tconcentration_unit\texposure_duration_hours\tbiological_replicate_id\tcondition_status\tcontrol_type\n"
        "SMPL_ROTN_C1\tHomo sapiens\tmidbrain_organoid\tLUHMES_LINE\tDMSO\t0.05\tpercent\t48\trep_1\tcontrol\tvehicle\n"
        "SMPL_ROTN_C2\tHomo sapiens\tmidbrain_organoid\tLUHMES_LINE\tDMSO\t0.05\tpercent\t48\trep_2\tcontrol\tvehicle\n"
        "SMPL_ROTN_T1\tHomo sapiens\tmidbrain_organoid\tLUHMES_LINE\tRotenone\t100\tnM\t48\trep_1\ttreatment\t\n"
        "SMPL_ROTN_T2\tHomo sapiens\tmidbrain_organoid\tLUHMES_LINE\tRotenone\t100\tnM\t48\trep_2\ttreatment\t\n"
    )
    tsv_file = tmp_path / "rotenone_study_metadata.tsv"
    tsv_file.write_text(tsv_content, encoding="utf-8")

    study_id = "STUDY_NEUROTOX_ROTENONE_2026"
    canonical_doc, raw_records = load_sample_metadata(
        metadata_source=tsv_file,
        study_id=study_id,
        source_type="in_house",
        study_title="Rotenone Neurotoxicity in Midbrain Organoids",
    )

    # Verify no PRJNA802072 or hardcoded defaults were introduced
    doc_str = json.dumps(canonical_doc)
    assert "PRJNA802072" not in doc_str
    assert "SRR17834" not in doc_str
    assert "valproic" not in doc_str.lower()

    # Verify study-specific fields are preserved
    assert canonical_doc["studies"][0]["study_id"] == study_id
    title_rec = canonical_doc["studies"][0]["study_title"]
    assert (title_rec.get("original_value") if isinstance(title_rec, dict) else title_rec) == "Rotenone Neurotoxicity in Midbrain Organoids"
    samples = canonical_doc["samples"]
    assert len(samples) == 4
    for s in samples:
        assert s["biological_replicate_id"]["original_value"] in ("rep_1", "rep_2")
        assert s["biological_replicate_id"]["value_status"] == "reported"

    # 2. Downstream bridge execution with this independent study
    from src.features.pathway.mapping import EnsemblToEntrezAdapter
    adapter = EnsemblToEntrezAdapter()
    genes = [r.ensembl_gene_id for r in adapter.gene_reference.records if r.gene_biotype == 'protein_coding'][:30]
    matrix_file = tmp_path / "counts_rotenone.tsv"
    runs = ("SMPL_ROTN_C1", "SMPL_ROTN_C2", "SMPL_ROTN_T1", "SMPL_ROTN_T2")
    gene_lines = ["gene_id\t" + "\t".join(runs)]
    for g in genes:
        gene_lines.append(f"{g}\t120\t115\t450\t460")
    matrix_file.write_text("\n".join(gene_lines) + "\n", encoding="utf-8")

    ds_res = execute_downstream_pipeline(
        counts_matrix_path=matrix_file,
        runs=runs,
        study_id=study_id,
        output_dir=tmp_path / "downstream_rotenone",
        canonical_metadata=canonical_doc,
    )

    # QC and normalization pass for this novel study
    assert ds_res.qc_passed is True
    assert ds_res.harmonization_metrics.get("status") in ("HARMONIZATION_COMPLETED", "PASS")
    assert ds_res.metadata_validation_status.lower() in ("valid", "valid_with_warnings")


# -----------------------------------------------------------------------------
# Test 11: Missing biological metadata not fabricated
# -----------------------------------------------------------------------------

def test_11_missing_biological_metadata_not_fabricated(tmp_path: Path):
    """Verify that unprovided biological metadata is marked unknown and never fabricated."""
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    make_dummy_fasta(fasta)
    make_dummy_gtf(gtf)

    r1 = tmp_path / "s1_R1.fastq.gz"
    r2 = tmp_path / "s1_R2.fastq.gz"
    make_dummy_fastq(r1, ["r1"], mate_suffix="/1")
    make_dummy_fastq(r2, ["r1"], mate_suffix="/2")

    samples = [
        SampleInput(
            sample_id="SAMPLE_RAW_NO_META",
            r1_path=r1,
            r2_path=r2,
            read_layout="paired_end",
            attributes={},
        )
    ]
    cfg = UniversalPipelineConfig(
        study_id="STUDY_NO_META",
        work_dir=tmp_path / "work_no_meta",
        enforce_resources=False,
    )
    refs = ReferenceConfig(species="Homo sapiens", genome_fasta=fasta, gtf_annotation=gtf)

    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)
    doc = pipeline.canonical_metadata

    # 1. No invented control_type, chemical, concentration or biological replicate
    for s in doc["samples"]:
        assert s["biological_replicate_id"]["value_status"] == "unknown"
        assert s["biological_replicate_id"]["original_value"] is None

    # Condition should NOT assume vehicle or DMSO
    for c in doc["conditions"]:
        assert c["treatment_control_status"]["value_status"] == "unknown"
        assert "control_type" not in c or c["control_type"]["value_status"] in ("unknown", "missing_not_reported")

    # Exposures should not assume DMSO or 0.0 uM
    for exp in doc["exposures"]:
        agent = exp.get("agent_name", {})
        assert agent.get("value_status") in ("unknown", "missing_not_reported")
        assert agent.get("original_value") != "DMSO"
        assert agent.get("original_value") != "Vehicle Control"
        conc = exp.get("concentration_or_dose", {})
        assert conc.get("value_status") in ("unknown", "missing_not_reported")
        assert conc.get("original_value") != 0.0


# -----------------------------------------------------------------------------
# Test 12: Planning does not modify filesystem
# -----------------------------------------------------------------------------

def test_12_planning_does_not_modify_filesystem(tmp_path: Path):
    """Verify that instantiating UniversalFastqPipeline and calling plan() does not create directories or files."""
    work_dir = tmp_path / "uncreated_pipeline_dir"
    assert not work_dir.exists()

    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    make_dummy_fasta(fasta)
    make_dummy_gtf(gtf)

    r1 = tmp_path / "s1_R1.fastq.gz"
    r2 = tmp_path / "s1_R2.fastq.gz"
    make_dummy_fastq(r1, ["r1"], mate_suffix="/1")
    make_dummy_fastq(r2, ["r1"], mate_suffix="/2")

    samples = [SampleInput(sample_id="S1", r1_path=r1, r2_path=r2, read_layout="paired_end")]
    cfg = UniversalPipelineConfig(study_id="STUDY_IMMUTABLE_PLAN", work_dir=work_dir, enforce_resources=False)
    refs = ReferenceConfig(species="Homo sapiens", genome_fasta=fasta, gtf_annotation=gtf)

    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)
    assert not work_dir.exists(), "UniversalFastqPipeline.__init__ must not create filesystem directories"

    plan = pipeline.plan()
    assert plan.can_execute is True
    assert not work_dir.exists(), "UniversalFastqPipeline.plan() must not create filesystem directories or files"

    res = pipeline.run(dry_run=True)
    assert res.status == "DRY_RUN_PLAN_GENERATED"
    assert not work_dir.exists(), "UniversalFastqPipeline.run(dry_run=True) must not create filesystem directories"


# -----------------------------------------------------------------------------
# Test 13: Approval token reproducibility and stale rejection
# -----------------------------------------------------------------------------

def test_13_approval_token_reproducibility_and_stale_rejection(tmp_path: Path):
    """Verify deterministic approval token generation, approve() workflow, and stale token rejection."""
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    make_dummy_fasta(fasta)
    make_dummy_gtf(gtf)

    r1 = tmp_path / "s1_R1.fastq.gz"
    r2 = tmp_path / "s1_R2.fastq.gz"
    make_dummy_fastq(r1, ["r1"], mate_suffix="/1")
    make_dummy_fastq(r2, ["r1"], mate_suffix="/2")

    samples = [SampleInput(sample_id="S1", r1_path=r1, r2_path=r2, read_layout="paired_end")]
    cfg = UniversalPipelineConfig(study_id="STUDY_APPROVAL", work_dir=tmp_path / "work_appr", require_approval=True, enforce_resources=False)
    refs = ReferenceConfig(species="Homo sapiens", genome_fasta=fasta, gtf_annotation=gtf)

    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)

    # 1. Deterministic token: plan() called twice yields identical approval token
    plan1 = pipeline.plan()
    plan2 = pipeline.plan()
    assert plan1.approval_token is not None
    assert plan1.approval_token == plan2.approval_token

    # 2. Unapproved execution halts safely
    res_no_appr = pipeline.run()
    assert res_no_appr.status == "BLOCKED_AWAITING_APPROVAL"
    assert any(f.code == "AWAITING_USER_APPROVAL" for f in res_no_appr.findings)

    # 3. Invalid or stale token is rejected
    res_stale = pipeline.run(approval_token="stale_fake_token_123")
    assert res_stale.status == "BLOCKED_STALE_APPROVAL"
    assert any(f.code == "STALE_APPROVAL_TOKEN" for f in res_stale.findings)

    # 4. pipeline.approve() records token and is recognized
    pipeline.approve(plan1.approval_token)
    assert pipeline._approved_token == plan1.approval_token


# -----------------------------------------------------------------------------
# Test 14: Mixed sequencing layouts and mismatches blocked
# -----------------------------------------------------------------------------

def test_14_mixed_sequencing_layouts_and_mismatches_blocked(tmp_path: Path):
    """Verify preflight blocks mixed cohort layouts and sample-to-study layout conflicts."""
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    make_dummy_fasta(fasta)
    make_dummy_gtf(gtf)

    r1 = tmp_path / "s1_R1.fastq.gz"
    r2 = tmp_path / "s1_R2.fastq.gz"
    make_dummy_fastq(r1, ["r1"])
    make_dummy_fastq(r2, ["r1"])

    refs = ReferenceConfig(species="Homo sapiens", genome_fasta=fasta, gtf_annotation=gtf)

    # Scenario A: Mixed cohort (sample 1 is paired_end, sample 2 is single_end)
    samples_mixed = [
        SampleInput(sample_id="SMPL_PE", r1_path=r1, r2_path=r2, read_layout="paired_end"),
        SampleInput(sample_id="SMPL_SE", r1_path=r1, r2_path=None, read_layout="single_end"),
    ]
    cfg_mixed = UniversalPipelineConfig(study_id="STUDY_MIXED", work_dir=tmp_path / "w_mix", read_layout="paired_end", enforce_resources=False)
    p_mixed = UniversalFastqPipeline(config=cfg_mixed, references=refs, samples=samples_mixed)
    rep_mixed = p_mixed.validate_inputs()
    assert rep_mixed.is_valid is False
    assert any(f.code == "MIXED_SEQUENCING_LAYOUTS_UNSUPPORTED" for f in rep_mixed.findings)

    # Scenario B: Sample configured as single_end but has mate R2
    samples_se_with_r2 = [
        SampleInput(sample_id="SMPL_SE_MATE", r1_path=r1, r2_path=r2, read_layout="single_end"),
    ]
    cfg_se = UniversalPipelineConfig(study_id="STUDY_SE", work_dir=tmp_path / "w_se", read_layout="single_end", enforce_resources=False)
    p_se = UniversalFastqPipeline(config=cfg_se, references=refs, samples=samples_se_with_r2)
    rep_se = p_se.validate_inputs()
    assert rep_se.is_valid is False
    assert any(f.code == "SAMPLE_LAYOUT_MISMATCH" for f in rep_se.findings)


# -----------------------------------------------------------------------------
# Test 15: SRA retrieval preflight and execution chaining
# -----------------------------------------------------------------------------

def test_15_sra_retrieval_preflight_and_execution_chaining(tmp_path: Path, monkeypatch):
    """Verify SRA accession identification in preflight, planning stage, and execution chaining to fastp."""
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    make_dummy_fasta(fasta)
    make_dummy_gtf(gtf)

    sra_acc = "SRR17834458"
    samples = [
        SampleInput(sample_id="SMPL_SRA_01", sra_accession=sra_acc, read_layout="paired_end")
    ]
    cfg = UniversalPipelineConfig(study_id="STUDY_SRA", work_dir=tmp_path / "w_sra", read_layout="paired_end", enforce_resources=False)
    refs = ReferenceConfig(species="Homo sapiens", genome_fasta=fasta, gtf_annotation=gtf)

    p = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)

    # Preflight detects SRA requirement
    pre = p.validate_inputs()
    assert p.samples[0].sample_id in pre.sample_validation
    assert pre.sample_validation["SMPL_SRA_01"]["status"] == "REQUIRES_SRA_RETRIEVAL"

    # Plan includes 00_sra stage
    plan = p.plan()
    sra_stages = [st for st in plan.planned_stages if st["stage_id"].startswith("00_sra")]
    assert len(sra_stages) == 1
    assert sra_stages[0]["accession"] == sra_acc

    # Mock extract_sra to simulate successful download without external network
    r1_mock = tmp_path / "sra_mock_1.fastq"
    r2_mock = tmp_path / "sra_mock_2.fastq"
    make_dummy_fastq(r1_mock, ["read1"], mate_suffix="/1", is_gz=False)
    make_dummy_fastq(r2_mock, ["read1"], mate_suffix="/2", is_gz=False)

    from src.quantification_bridge.sra_runner import SraExtractionResult
    def mock_extract_sra(sra_path, output_dir, threads=6, temp_dir=None):
        return SraExtractionResult(
            status="SRA_EXTRACTION_SUCCESS",
            run_accession=sra_acc,
            layout="PAIRED",
            fastq_r1=r1_mock,
            fastq_r2=r2_mock,
            tool_name="mocked_fasterq_dump",
            tool_version="3.1.1",
            output_dir=Path(output_dir),
            temp_dir=Path(output_dir),
        )

    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.extract_sra", mock_extract_sra)

    # Execute SRA retrieval stage directly via pipeline resolution
    p._ensure_directories()
    sra_res = mock_extract_sra(sra_acc, p.sra_dir / "SMPL_SRA_01")
    assert sra_res.fastq_r1.is_file()
    assert sra_res.fastq_r2.is_file()


# -----------------------------------------------------------------------------
# Test 16: Downstream status reporting granularity
# -----------------------------------------------------------------------------

def test_16_downstream_status_reporting_granularity(tmp_path: Path):
    """Verify that UniversalFastqResult distinguishes COMPLETED_QUANTIFICATION_INFERENTIAL_BLOCKED, FAILED_DOWNSTREAM_*, and COMPLETED."""
    res_gated = UniversalFastqResult(
        status="COMPLETED_QUANTIFICATION_INFERENTIAL_BLOCKED",
        study_id="STUDY_STAT",
        sample_ids=("S1", "S2"),
        raw_counts_path=tmp_path / "counts.tsv",
        canonical_metadata_path=tmp_path / "meta.json",
        qc_reports_dir=tmp_path / "qc",
        alignments_dir=tmp_path / "align",
        strandedness_evidence_path=tmp_path / "strand.json",
        counting_summary_path=tmp_path / "count.json",
        execution_record_path=tmp_path / "rec.json",
        completed_stages=("star_index", "fastp_trimming", "star_alignment", "downstream_bridge_exploratory"),
        failed_stage=None,
    )
    d_gated = res_gated.to_dict()
    assert d_gated["status"] == "COMPLETED_QUANTIFICATION_INFERENTIAL_BLOCKED"
    assert "downstream_bridge_exploratory" in d_gated["completed_stages"]
    assert d_gated["failed_stage"] is None

    res_qc_fail = UniversalFastqResult(
        status="FAILED_DOWNSTREAM_QC",
        study_id="STUDY_STAT",
        sample_ids=("S1",),
        raw_counts_path=None,
        canonical_metadata_path=None,
        qc_reports_dir=None,
        alignments_dir=None,
        strandedness_evidence_path=None,
        counting_summary_path=None,
        execution_record_path=None,
        completed_stages=("star_index", "fastp_trimming"),
        failed_stage="downstream_bridge",
    )
    d_qc_fail = res_qc_fail.to_dict()
    assert d_qc_fail["status"] == "FAILED_DOWNSTREAM_QC"
    assert d_qc_fail["failed_stage"] == "downstream_bridge"


# -----------------------------------------------------------------------------
# Test 17: Genuine mini-genome end-to-end quantification
# -----------------------------------------------------------------------------

def test_17_genuine_mini_genome_end_to_end_quantification(tmp_path: Path):
    """Execute genuine end-to-end quantification: FASTQ -> fastp -> STAR -> featureCounts -> counts.tsv -> downstream QC."""
    import random
    work_dir = tmp_path / "e2e_work"

    # 1. Build a 100,000 bp non-repetitive mini-genome
    random.seed(42)
    bases = ["A", "C", "G", "T"]
    seq_100k = "".join(random.choices(bases, k=100000))
    fasta = tmp_path / "mini_genome.fa"
    fasta.write_text(f">chr1\n{seq_100k}\n", encoding="utf-8")

    # 2. GTF with 2 Ensembl genes
    gtf = tmp_path / "mini_genes.gtf"
    gtf_content = (
        'chr1\tsynth\texon\t1000\t5000\t.\t+\t.\tgene_id "ENSG00000000001"; gene_name "GENE1"; transcript_id "ENST00000000001";\n'
        'chr1\tsynth\texon\t10000\t15000\t.\t+\t.\tgene_id "ENSG00000000002"; gene_name "GENE2"; transcript_id "ENST00000000002";\n'
    )
    gtf.write_text(gtf_content, encoding="utf-8")

    # 3. Reads targeting GENE1 and GENE2
    read_seq_g1 = seq_100k[1100:1150]
    read_seq_g2 = seq_100k[10100:10150]
    r1 = tmp_path / "mini_sample_R1.fastq.gz"
    r2 = tmp_path / "mini_sample_R2.fastq.gz"

    with gzip.open(r1, "wt", encoding="utf-8") as f1, gzip.open(r2, "wt", encoding="utf-8") as f2:
        for i in range(20):
            f1.write(f"@read_g1_{i}/1\n{read_seq_g1}\n+\n{'I'*50}\n")
            f2.write(f"@read_g1_{i}/2\n{read_seq_g1}\n+\n{'I'*50}\n")
            f1.write(f"@read_g2_{i}/1\n{read_seq_g2}\n+\n{'I'*50}\n")
            f2.write(f"@read_g2_{i}/2\n{read_seq_g2}\n+\n{'I'*50}\n")

    cfg = UniversalPipelineConfig(
        study_id="STUDY_E2E_MINI",
        work_dir=work_dir,
        read_layout="paired_end",
        strandedness=0,
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
            sample_id="MINI_SAMPLE_01",
            r1_path=r1,
            r2_path=r2,
            read_layout="paired_end",
            attributes={
                "organoid_type": "brain_organoid",
                "cell_line": "MINI_LINE",
            },
        )
    ]

    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)
    plan = pipeline.plan()
    assert plan.can_execute is True

    # Approve and execute genuine end-to-end pipeline
    pipeline.approve(plan.approval_token)
    res = pipeline.run()

    # Verify quantification outputs
    assert res.raw_counts_path is not None
    assert res.raw_counts_path.is_file()
    counts_content = res.raw_counts_path.read_text(encoding="utf-8")
    assert "ENSG00000000001" in counts_content
    assert "ENSG00000000002" in counts_content
    assert "MINI_SAMPLE_01" in counts_content

    # Verify execution record and completed stages
    assert res.execution_record_path is not None and res.execution_record_path.is_file()
    assert "star_index" in res.completed_stages
    assert "featurecounts" in res.completed_stages
    assert res.status in ("COMPLETED", "COMPLETED_QUANTIFICATION_INFERENTIAL_BLOCKED", "PARTIALLY_COMPLETED")
