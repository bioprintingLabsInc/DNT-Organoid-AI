"""Tests for SRA Checkpoint Verification and Resumption Integrity.

Mandatory regression test suite for Final Targeted F1.1 SRA checkpoint-resume fix.
Validates:
1. Valid unchanged R1/R2 and checkpoint -> Resume without repeating SRA extraction
2. R1 modified without changing file size -> Reuse rejected
3. R2 modified without changing file size -> Reuse rejected
4. A large FASTQ changed only in the middle -> Reuse rejected
5. One mate deleted or truncated -> Reuse rejected
6. R1/R2 exchanged or wrong files associated with sample -> Reuse rejected
7. Changed SRA accession or incompatible read layout -> Reuse rejected
8. Missing, malformed or incomplete checkpoint -> Reuse rejected safely
9. Old checkpoint without required integrity evidence -> Reuse rejected safely
10. Extra unrelated FASTQ in the directory -> Cannot be mistaken for a recorded mate
11. Valid resume followed by fastp -> Exact verified FASTQ paths reach fastp
12. Invalid resume -> No unapproved download, deletion or downstream processing
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

import pytest

from src.quantification_bridge.sra_runner import SraExtractionResult
from src.quantification_bridge.universal_pipeline import (
    CheckpointManager,
    ReferenceConfig,
    SampleInput,
    SraCheckpointVerification,
    UniversalFastqPipeline,
    UniversalPipelineConfig,
    compute_file_content_hash,
    compute_full_file_sha256,
)


# -----------------------------------------------------------------------------
# Test Helpers
# -----------------------------------------------------------------------------

def make_synthetic_fastq(
    path: Path,
    num_reads: int = 10,
    read_len: int = 50,
    is_gz: bool = True,
    mate_suffix: str = "/1",
    id_prefix: str = "READ",
) -> Path:
    """Generate syntactically valid synthetic FASTQ or FASTQ.GZ file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    seq = "A" * read_len
    qual = "I" * read_len
    for i in range(1, num_reads + 1):
        lines.append(f"@{id_prefix}_{i}{mate_suffix}")
        lines.append(seq)
        lines.append("+")
        lines.append(qual)
    text = "\n".join(lines) + "\n"
    if is_gz:
        with gzip.open(path, "wb") as gz:
            gz.write(text.encode("utf-8"))
    else:
        path.write_text(text, encoding="utf-8")
    return path


def make_dummy_references(tmp_path: Path) -> ReferenceConfig:
    """Create minimal FASTA and GTF reference files for pipeline testing."""
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    fasta.write_text(">chr1\nACGTACGTACGTACGTACGT\n", encoding="utf-8")
    gtf.write_text(
        'chr1\tHAVANA\texon\t1\t100\t.\t+\t.\tgene_id "ENSG000001"; transcript_id "ENST000001"; gene_name "GENE1";\n',
        encoding="utf-8",
    )
    return ReferenceConfig(species="Homo sapiens", genome_fasta=fasta, gtf_annotation=gtf)


# -----------------------------------------------------------------------------
# Case 1: Valid unchanged R1/R2 and checkpoint -> Resume without repeating extraction
# -----------------------------------------------------------------------------

def test_01_valid_unchanged_r1_r2_and_checkpoint_resumes(tmp_path: Path):
    """Verify that unchanged FASTQ files and valid checkpoint resume cleanly for paired-end and single-end."""
    chk_dir = tmp_path / ".checkpoints"
    sra_dir = tmp_path / "sra_extracted"
    mgr = CheckpointManager(chk_dir)

    # 1. Paired-end
    sample_id_pe = "SMPL_PE1"
    acc_pe = "SRR10000001"
    pe_dir = sra_dir / sample_id_pe
    r1_pe = pe_dir / f"{acc_pe}_1.fastq.gz"
    r2_pe = pe_dir / f"{acc_pe}_2.fastq.gz"
    make_synthetic_fastq(r1_pe, num_reads=5, mate_suffix="/1")
    make_synthetic_fastq(r2_pe, num_reads=5, mate_suffix="/2")

    input_hash_pe = "hash_in_pe"
    config_hash_pe = "hash_cfg_pe"
    mgr.record_sra_checkpoint(
        sample_id=sample_id_pe,
        sra_accession=acc_pe,
        read_layout="paired_end",
        fastq_r1=r1_pe,
        fastq_r2=r2_pe,
        input_hash=input_hash_pe,
        config_hash=config_hash_pe,
    )

    ver_pe = mgr.verify_sra_checkpoint(
        sample_id=sample_id_pe,
        sra_accession=acc_pe,
        read_layout="paired_end",
        input_hash=input_hash_pe,
        config_hash=config_hash_pe,
    )
    assert ver_pe.is_valid is True
    assert bool(ver_pe) is True
    assert ver_pe.fastq_r1 == r1_pe.resolve()
    assert ver_pe.fastq_r2 == r2_pe.resolve()

    # 2. Single-end
    sample_id_se = "SMPL_SE1"
    acc_se = "SRR10000002"
    se_dir = sra_dir / sample_id_se
    r1_se = se_dir / f"{acc_se}_1.fastq.gz"
    make_synthetic_fastq(r1_se, num_reads=5, mate_suffix="")

    input_hash_se = "hash_in_se"
    config_hash_se = "hash_cfg_se"
    mgr.record_sra_checkpoint(
        sample_id=sample_id_se,
        sra_accession=acc_se,
        read_layout="single_end",
        fastq_r1=r1_se,
        fastq_r2=None,
        input_hash=input_hash_se,
        config_hash=config_hash_se,
    )

    ver_se = mgr.verify_sra_checkpoint(
        sample_id=sample_id_se,
        sra_accession=acc_se,
        read_layout="single_end",
        input_hash=input_hash_se,
        config_hash=config_hash_se,
    )
    assert ver_se.is_valid is True
    assert ver_se.fastq_r1 == r1_se.resolve()
    assert ver_se.fastq_r2 is None


# -----------------------------------------------------------------------------
# Case 2: R1 modified without changing file size -> Reuse rejected
# -----------------------------------------------------------------------------

def test_02_r1_modified_without_changing_file_size_rejects_reuse(tmp_path: Path):
    """Verify that same-size modifications to R1 are detected and rejected."""
    chk_dir = tmp_path / ".checkpoints"
    sra_dir = tmp_path / "sra_extracted"
    mgr = CheckpointManager(chk_dir)

    sample_id = "SMPL_MOD_R1"
    acc = "SRR10000003"
    pe_dir = sra_dir / sample_id
    r1 = pe_dir / f"{acc}_1.fastq.gz"
    r2 = pe_dir / f"{acc}_2.fastq.gz"
    make_synthetic_fastq(r1, num_reads=10, mate_suffix="/1")
    make_synthetic_fastq(r2, num_reads=10, mate_suffix="/2")

    mgr.record_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        fastq_r1=r1,
        fastq_r2=r2,
        input_hash="inhash",
        config_hash="cfghash",
    )

    # Modify 1 byte in R1 in-place (keeping size strictly identical)
    raw = bytearray(r1.read_bytes())
    orig_len = len(raw)
    raw[25] = (raw[25] + 1) % 256
    r1.write_bytes(raw)
    assert len(r1.read_bytes()) == orig_len

    ver = mgr.verify_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        input_hash="inhash",
        config_hash="cfghash",
    )
    assert ver.is_valid is False
    assert "r1" in ver.reason.lower()
    assert "checksum mismatch" in ver.reason.lower()


# -----------------------------------------------------------------------------
# Case 3: R2 modified without changing file size -> Reuse rejected
# -----------------------------------------------------------------------------

def test_03_r2_modified_without_changing_file_size_rejects_reuse(tmp_path: Path):
    """Verify that same-size modifications to R2 are detected and rejected."""
    chk_dir = tmp_path / ".checkpoints"
    sra_dir = tmp_path / "sra_extracted"
    mgr = CheckpointManager(chk_dir)

    sample_id = "SMPL_MOD_R2"
    acc = "SRR10000004"
    pe_dir = sra_dir / sample_id
    r1 = pe_dir / f"{acc}_1.fastq.gz"
    r2 = pe_dir / f"{acc}_2.fastq.gz"
    make_synthetic_fastq(r1, num_reads=10, mate_suffix="/1")
    make_synthetic_fastq(r2, num_reads=10, mate_suffix="/2")

    mgr.record_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        fastq_r1=r1,
        fastq_r2=r2,
        input_hash="inhash",
        config_hash="cfghash",
    )

    # Modify 1 byte in R2 in-place (keeping size strictly identical)
    raw = bytearray(r2.read_bytes())
    orig_len = len(raw)
    raw[25] = (raw[25] + 1) % 256
    r2.write_bytes(raw)
    assert len(r2.read_bytes()) == orig_len

    ver = mgr.verify_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        input_hash="inhash",
        config_hash="cfghash",
    )
    assert ver.is_valid is False
    assert "r2" in ver.reason.lower()
    assert "checksum mismatch" in ver.reason.lower()


# -----------------------------------------------------------------------------
# Case 4: A large FASTQ changed only in the middle -> Reuse rejected
# -----------------------------------------------------------------------------

def test_04_large_fastq_changed_only_in_middle_rejects_reuse(tmp_path: Path):
    """Verify that streaming full-file SHA-256 detects middle modifications that truncated hashing misses."""
    chk_dir = tmp_path / ".checkpoints"
    sra_dir = tmp_path / "sra_extracted"
    mgr = CheckpointManager(chk_dir)

    sample_id = "SMPL_LARGE"
    acc = "SRR10000005"
    pe_dir = sra_dir / sample_id
    pe_dir.mkdir(parents=True, exist_ok=True)
    r1 = pe_dir / f"{acc}_1.fastq"

    # Create a 12 MB synthetic uncompressed FASTQ
    # Total size: ~12,000,000 bytes (> 10MB cutoff of compute_file_content_hash)
    target_size = 12 * 1024 * 1024
    block = (b"@READ_HEADER\n" + b"A" * 80 + b"\n+\n" + b"I" * 80 + b"\n") * 1000
    with r1.open("wb") as f:
        written = 0
        while written < target_size:
            f.write(block)
            written += len(block)

    orig_size = r1.stat().st_size
    assert orig_size >= 12 * 1024 * 1024

    mgr.record_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="single_end",
        fastq_r1=r1,
        fastq_r2=None,
        input_hash="inhash",
        config_hash="cfghash",
    )

    # Hash before modification
    orig_content_hash = compute_file_content_hash(r1)
    orig_full_sha = compute_full_file_sha256(r1)

    # Modify 4 bytes strictly in the middle at offset 10,500,000
    # (Beyond 10MB cutoff, and more than 64KB before EOF)
    mod_offset = 10_500_000
    with r1.open("r+b") as f:
        f.seek(mod_offset)
        f.write(b"TTTT")

    # Confirm size is completely unchanged
    assert r1.stat().st_size == orig_size

    # Prove that the general truncated compute_file_content_hash misses this middle change:
    truncated_hash_after = compute_file_content_hash(r1)
    assert truncated_hash_after == orig_content_hash, "compute_file_content_hash should miss middle-only change"

    # Prove that compute_full_file_sha256 DOES detect this middle change:
    full_sha_after = compute_full_file_sha256(r1)
    assert full_sha_after != orig_full_sha, "compute_full_file_sha256 must detect middle-only change"

    # Verify SRA checkpoint correctly rejects reuse
    ver = mgr.verify_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="single_end",
        input_hash="inhash",
        config_hash="cfghash",
    )
    assert ver.is_valid is False
    assert "checksum mismatch" in ver.reason.lower()


# -----------------------------------------------------------------------------
# Case 5: One mate deleted or truncated -> Reuse rejected
# -----------------------------------------------------------------------------

def test_05_one_mate_deleted_or_truncated_rejects_reuse(tmp_path: Path):
    """Verify that deleting or truncating either mate rejects checkpoint reuse."""
    chk_dir = tmp_path / ".checkpoints"
    sra_dir = tmp_path / "sra_extracted"
    mgr = CheckpointManager(chk_dir)

    sample_id = "SMPL_DEL_TRUNC"
    acc = "SRR10000006"
    pe_dir = sra_dir / sample_id
    r1 = pe_dir / f"{acc}_1.fastq.gz"
    r2 = pe_dir / f"{acc}_2.fastq.gz"
    make_synthetic_fastq(r1, num_reads=10, mate_suffix="/1")
    make_synthetic_fastq(r2, num_reads=10, mate_suffix="/2")

    mgr.record_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        fastq_r1=r1,
        fastq_r2=r2,
        input_hash="inhash",
        config_hash="cfghash",
    )

    # Subcase A: R2 deleted
    r2.unlink()
    ver_del = mgr.verify_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        input_hash="inhash",
        config_hash="cfghash",
    )
    assert ver_del.is_valid is False
    assert "does not exist" in ver_del.reason.lower()

    # Re-create R2 and truncate R1
    make_synthetic_fastq(r2, num_reads=10, mate_suffix="/2")
    r1.write_bytes(r1.read_bytes()[:15])  # truncated

    ver_trunc = mgr.verify_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        input_hash="inhash",
        config_hash="cfghash",
    )
    assert ver_trunc.is_valid is False
    assert ("size mismatch" in ver_trunc.reason.lower() or "structural format" in ver_trunc.reason.lower())


# -----------------------------------------------------------------------------
# Case 6: R1/R2 exchanged or wrong files associated with sample -> Reuse rejected
# -----------------------------------------------------------------------------

def test_06_r1_r2_exchanged_or_wrong_files_rejects_reuse(tmp_path: Path):
    """Verify that exchanging R1/R2 on disk or associating wrong sample files is rejected."""
    chk_dir = tmp_path / ".checkpoints"
    sra_dir = tmp_path / "sra_extracted"
    mgr = CheckpointManager(chk_dir)

    sample_id = "SMPL_EXCHANGE"
    acc = "SRR10000007"
    pe_dir = sra_dir / sample_id
    r1 = pe_dir / f"{acc}_1.fastq.gz"
    r2 = pe_dir / f"{acc}_2.fastq.gz"
    make_synthetic_fastq(r1, num_reads=10, mate_suffix="/1", id_prefix="FORWARD")
    make_synthetic_fastq(r2, num_reads=10, mate_suffix="/2", id_prefix="REVERSE")

    mgr.record_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        fastq_r1=r1,
        fastq_r2=r2,
        input_hash="inhash",
        config_hash="cfghash",
    )

    # Subcase A: Swap contents of R1 and R2 on disk
    r1_bytes = r1.read_bytes()
    r2_bytes = r2.read_bytes()
    r1.write_bytes(r2_bytes)
    r2.write_bytes(r1_bytes)

    ver_swapped = mgr.verify_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        input_hash="inhash",
        config_hash="cfghash",
    )
    assert ver_swapped.is_valid is False
    assert "checksum mismatch" in ver_swapped.reason.lower()

    # Subcase B: Wrong files associated with sample in checkpoint
    other_r1 = sra_dir / "OTHER_SAMPLE" / "SRR999999_1.fastq.gz"
    other_r2 = sra_dir / "OTHER_SAMPLE" / "SRR999999_2.fastq.gz"
    make_synthetic_fastq(other_r1, num_reads=5, mate_suffix="/1")
    make_synthetic_fastq(other_r2, num_reads=5, mate_suffix="/2")

    mgr.record_sra_checkpoint(
        sample_id="WRONG_SAMPLE",
        sra_accession="SRR10000008",
        read_layout="paired_end",
        fastq_r1=other_r1,
        fastq_r2=other_r2,
        input_hash="inhash",
        config_hash="cfghash",
    )
    ver_wrong = mgr.verify_sra_checkpoint(
        sample_id="WRONG_SAMPLE",
        sra_accession="SRR10000008",
        read_layout="paired_end",
        input_hash="inhash",
        config_hash="cfghash",
    )
    assert ver_wrong.is_valid is False
    assert "not associated" in ver_wrong.reason.lower()


# -----------------------------------------------------------------------------
# Case 7: Changed SRA accession or incompatible read layout -> Reuse rejected
# -----------------------------------------------------------------------------

def test_07_changed_sra_accession_or_incompatible_layout_rejects_reuse(tmp_path: Path):
    """Verify that querying a checkpoint with a different SRA accession or layout rejects reuse."""
    chk_dir = tmp_path / ".checkpoints"
    sra_dir = tmp_path / "sra_extracted"
    mgr = CheckpointManager(chk_dir)

    sample_id = "SMPL_MISMATCH"
    acc = "SRR10000009"
    pe_dir = sra_dir / sample_id
    r1 = pe_dir / f"{acc}_1.fastq.gz"
    r2 = pe_dir / f"{acc}_2.fastq.gz"
    make_synthetic_fastq(r1, num_reads=5, mate_suffix="/1")
    make_synthetic_fastq(r2, num_reads=5, mate_suffix="/2")

    mgr.record_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        fastq_r1=r1,
        fastq_r2=r2,
        input_hash="inhash",
        config_hash="cfghash",
    )

    # Subcase A: Different accession queried
    ver_acc = mgr.verify_sra_checkpoint(
        sample_id=sample_id,
        sra_accession="SRR99999999",
        read_layout="paired_end",
        input_hash="inhash",
        config_hash="cfghash",
    )
    assert ver_acc.is_valid is False
    assert "accession mismatch" in ver_acc.reason.lower()

    # Subcase B: Different layout queried (single_end instead of paired_end)
    ver_layout = mgr.verify_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="single_end",
        input_hash="inhash",
        config_hash="cfghash",
    )
    assert ver_layout.is_valid is False
    assert "layout mismatch" in ver_layout.reason.lower()


# -----------------------------------------------------------------------------
# Case 8: Missing, malformed or incomplete checkpoint -> Reuse rejected safely
# -----------------------------------------------------------------------------

def test_08_missing_malformed_or_incomplete_checkpoint_rejects_safely(tmp_path: Path):
    """Verify that missing, malformed, or incomplete checkpoints are safely rejected."""
    chk_dir = tmp_path / ".checkpoints"
    chk_dir.mkdir(parents=True, exist_ok=True)
    mgr = CheckpointManager(chk_dir)

    # Subcase A: Missing checkpoint
    ver_missing = mgr.verify_sra_checkpoint(
        sample_id="NON_EXISTENT",
        sra_accession="SRR10000010",
        read_layout="paired_end",
    )
    assert ver_missing.is_valid is False
    assert "not found" in ver_missing.reason.lower()

    # Subcase B: Malformed JSON
    bad_chk = chk_dir / "sra_extract_SMPL_BAD.json"
    bad_chk.write_text("{this is not valid json}", encoding="utf-8")
    ver_bad = mgr.verify_sra_checkpoint(
        sample_id="SMPL_BAD",
        sra_accession="SRR10000011",
        read_layout="paired_end",
    )
    assert ver_bad.is_valid is False
    assert "malformed" in ver_bad.reason.lower() or "invalid json" in ver_bad.reason.lower()

    # Subcase C: Incomplete checkpoint (missing files dictionary)
    inc_chk = chk_dir / "sra_extract_SMPL_INC.json"
    inc_chk.write_text(json.dumps({"stage": "sra_extract", "status": "COMPLETED"}), encoding="utf-8")
    ver_inc = mgr.verify_sra_checkpoint(
        sample_id="SMPL_INC",
        sra_accession="SRR10000012",
        read_layout="paired_end",
    )
    assert ver_inc.is_valid is False
    assert "missing" in ver_inc.reason.lower() or "mismatch" in ver_inc.reason.lower()


# -----------------------------------------------------------------------------
# Case 9: Old checkpoint without required integrity evidence -> Reuse rejected safely
# -----------------------------------------------------------------------------

def test_09_old_checkpoint_without_integrity_evidence_rejects_safely(tmp_path: Path):
    """Verify that earlier legacy checkpoints lacking exact per-file integrity evidence are rejected."""
    chk_dir = tmp_path / ".checkpoints"
    chk_dir.mkdir(parents=True, exist_ok=True)
    mgr = CheckpointManager(chk_dir)

    sample_id = "SMPL_LEGACY"
    acc = "SRR10000013"
    legacy_chk = chk_dir / f"sra_extract_{sample_id}.json"

    # Legacy format: outputs list without 'files' dict or per-file cryptographic sha256
    legacy_data = {
        "stage": "sra_extract",
        "identifier": sample_id,
        "sample_id": sample_id,
        "sra_accession": acc,
        "read_layout": "paired_end",
        "status": "COMPLETED",
        "input_hash": "inhash",
        "config_hash": "cfghash",
        "outputs": [f"/tmp/{acc}_1.fastq.gz", f"/tmp/{acc}_2.fastq.gz"],
        "output_hashes": {f"/tmp/{acc}_1.fastq.gz": "sampled_hash"},
        "metrics": {},
    }
    legacy_chk.write_text(json.dumps(legacy_data, indent=2), encoding="utf-8")

    ver = mgr.verify_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        input_hash="inhash",
        config_hash="cfghash",
    )
    assert ver.is_valid is False
    assert "legacy" in ver.reason.lower()


# -----------------------------------------------------------------------------
# Case 10: Extra unrelated FASTQ in directory -> Cannot be mistaken for recorded mate
# -----------------------------------------------------------------------------

def test_10_extra_unrelated_fastq_in_directory_ignored(tmp_path: Path):
    """Verify that exact recorded paths prevent extra files in the directory from being mistaken for mates."""
    chk_dir = tmp_path / ".checkpoints"
    sra_dir = tmp_path / "sra_extracted"
    mgr = CheckpointManager(chk_dir)

    sample_id = "SMPL_EXTRA_FILES"
    acc = "SRR10000014"
    pe_dir = sra_dir / sample_id
    r1 = pe_dir / f"{acc}_1.fastq.gz"
    r2 = pe_dir / f"{acc}_2.fastq.gz"
    make_synthetic_fastq(r1, num_reads=5, mate_suffix="/1")
    make_synthetic_fastq(r2, num_reads=5, mate_suffix="/2")

    # Add extra confusing files with matching suffixes in the same directory
    extra_1 = pe_dir / f"{acc}_extra_1.fastq.gz"
    extra_2 = pe_dir / "unrelated_sample_2.fastq.gz"
    make_synthetic_fastq(extra_1, num_reads=10, mate_suffix="/1")
    make_synthetic_fastq(extra_2, num_reads=10, mate_suffix="/2")

    mgr.record_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        fastq_r1=r1,
        fastq_r2=r2,
        input_hash="inhash",
        config_hash="cfghash",
    )

    ver = mgr.verify_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        input_hash="inhash",
        config_hash="cfghash",
    )
    assert ver.is_valid is True
    # The verified paths MUST be the exact recorded files, not the extra files
    assert ver.fastq_r1 == r1.resolve()
    assert ver.fastq_r2 == r2.resolve()
    assert ver.fastq_r1 != extra_1.resolve()
    assert ver.fastq_r2 != extra_2.resolve()


# -----------------------------------------------------------------------------
# Case 11: Valid resume followed by fastp -> Exact verified FASTQ paths reach fastp
# -----------------------------------------------------------------------------

def test_11_valid_resume_followed_by_fastp(tmp_path: Path, monkeypatch):
    """Verify that during a valid resume, exact verified SRA FASTQ paths reach fastp without re-extracting."""
    work_dir = tmp_path / "work_resume"
    refs = make_dummy_references(tmp_path)

    sample_id = "SMPL_RESUME_PE"
    acc = "SRR10000015"
    samples = [SampleInput(sample_id=sample_id, sra_accession=acc, read_layout="paired_end")]
    cfg = UniversalPipelineConfig(
        study_id="STUDY_RESUME",
        work_dir=work_dir,
        read_layout="paired_end",
        strandedness=0,
        enforce_resources=False,
    )

    # Pre-populate exact SRA outputs in designated directory
    sra_out_dir = work_dir / "sra_extracted" / sample_id
    sra_out_dir.mkdir(parents=True, exist_ok=True)
    expected_r1 = sra_out_dir / f"{acc}_1.fastq.gz"
    expected_r2 = sra_out_dir / f"{acc}_2.fastq.gz"
    make_synthetic_fastq(expected_r1, num_reads=8, mate_suffix="/1")
    make_synthetic_fastq(expected_r2, num_reads=8, mate_suffix="/2")

    # Pre-record valid SRA checkpoint
    chk_dir = work_dir / ".checkpoints"
    mgr = CheckpointManager(chk_dir)
    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)
    input_hash = pipeline._compute_hash(acc)
    config_hash = pipeline._compute_hash(sample_id, str(cfg.threads))

    mgr.record_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        fastq_r1=expected_r1,
        fastq_r2=expected_r2,
        input_hash=input_hash,
        config_hash=config_hash,
    )

    # Track invocations: SRA extraction must NOT be called; fastp MUST be called with exact paths
    sra_calls: list[dict[str, Any]] = []
    def mock_extract_sra(*args, **kwargs):
        sra_calls.append({"args": args, "kwargs": kwargs})
        raise AssertionError("extract_sra should NOT be invoked on valid resume!")

    fastp_calls: list[dict[str, Any]] = []
    def mock_fastp(raw_r1, raw_r2, out_r1, out_r2, json_report, html_report, sample_id, threads=4):
        fastp_calls.append({"raw_r1": Path(raw_r1), "raw_r2": Path(raw_r2) if raw_r2 else None})
        Path(out_r1).parent.mkdir(parents=True, exist_ok=True)
        make_synthetic_fastq(Path(out_r1), num_reads=4)
        if out_r2:
            make_synthetic_fastq(Path(out_r2), num_reads=4)
        Path(json_report).write_text(json.dumps({"summary": {"before_filtering": {"total_reads": 8}}}), encoding="utf-8")
        return {"metrics_summary": {"total_reads": 8}}

    def mock_align(trimmed_r1, trimmed_r2, index_dir, output_prefix, threads=4, profile=None, **kwargs):
        prefix_p = Path(output_prefix)
        prefix_p.parent.mkdir(parents=True, exist_ok=True)
        bam = prefix_p.parent / f"{prefix_p.name}Aligned.sortedByCoord.out.bam"
        bam.write_bytes(b"BAM\x01mock_bam_content")
        final_log = prefix_p.parent / f"{prefix_p.name}Log.final.out"
        final_log.write_text("Uniquely mapped reads % | 95.0%\nNumber of input reads | 100\n", encoding="utf-8")
        return {"bam_path": bam, "log_final_path": final_log, "metrics": {"uniquely_mapped_percent": 95.0}}

    def mock_build_star(fasta_path, gtf_path, index_dir, **kwargs):
        idx_p = Path(index_dir)
        idx_p.mkdir(parents=True, exist_ok=True)
        (idx_p / "Genome").write_bytes(b"mock genome")
        (idx_p / "SA").write_bytes(b"mock sa")

    def mock_fc(bam_paths, gtf_path, output_counts_path, **kwargs):
        out_f = Path(output_counts_path)
        out_f.parent.mkdir(parents=True, exist_ok=True)
        header = "Geneid\tChr\tStart\tEnd\tStrand\tLength\t" + "\t".join([Path(b).name for b in bam_paths]) + "\n"
        row = "ENSG1\tchr1\t1\t100\t+\t100\t" + "\t".join(["50" for _ in bam_paths]) + "\n"
        out_f.write_text(header + row, encoding="utf-8")
        summary_f = out_f.parent / f"{out_f.name}.summary"
        summary_f.write_text("Status\tSMPL\nAssigned\t50\n", encoding="utf-8")

    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.extract_sra", mock_extract_sra)
    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.run_fastp_preprocessing", mock_fastp)
    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.align_reads_star", mock_align)
    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.build_star_index", mock_build_star)
    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.run_featurecounts", mock_fc)

    # Approve and execute
    plan = pipeline.plan()
    pipeline.approve(plan.approval_token)
    res = pipeline.run()

    assert len(sra_calls) == 0, "SRA extraction must not be called during valid resume"
    assert len(fastp_calls) == 1, "fastp must be executed once"
    # Exact verified paths must reach fastp
    assert fastp_calls[0]["raw_r1"] == expected_r1.resolve()
    assert fastp_calls[0]["raw_r2"] == expected_r2.resolve()
    assert "sra_extraction" in res.completed_stages


# -----------------------------------------------------------------------------
# Case 12: Invalid resume -> No unapproved download, deletion or downstream processing
# -----------------------------------------------------------------------------

def test_12_invalid_resume_halts_without_download_deletion_or_downstream(tmp_path: Path, monkeypatch):
    """Verify that invalid SRA checkpoints halt safely without downloading, deleting files, or running downstream stages."""
    work_dir = tmp_path / "work_invalid_resume"
    refs = make_dummy_references(tmp_path)

    sample_id = "SMPL_INVALID_RESUME"
    acc = "SRR10000016"
    samples = [SampleInput(sample_id=sample_id, sra_accession=acc, read_layout="paired_end")]
    cfg = UniversalPipelineConfig(
        study_id="STUDY_INVALID_RESUME",
        work_dir=work_dir,
        read_layout="paired_end",
        strandedness=0,
        enforce_resources=False,
    )

    sra_out_dir = work_dir / "sra_extracted" / sample_id
    sra_out_dir.mkdir(parents=True, exist_ok=True)
    r1 = sra_out_dir / f"{acc}_1.fastq.gz"
    r2 = sra_out_dir / f"{acc}_2.fastq.gz"
    make_synthetic_fastq(r1, num_reads=8, mate_suffix="/1")
    make_synthetic_fastq(r2, num_reads=8, mate_suffix="/2")

    # Record checkpoint
    chk_dir = work_dir / ".checkpoints"
    mgr = CheckpointManager(chk_dir)
    pipeline = UniversalFastqPipeline(config=cfg, references=refs, samples=samples)
    input_hash = pipeline._compute_hash(acc)
    config_hash = pipeline._compute_hash(sample_id, str(cfg.threads))

    mgr.record_sra_checkpoint(
        sample_id=sample_id,
        sra_accession=acc,
        read_layout="paired_end",
        fastq_r1=r1,
        fastq_r2=r2,
        input_hash=input_hash,
        config_hash=config_hash,
    )

    # Corrupt R1 by modifying 1 byte (checksum mismatch)
    raw = bytearray(r1.read_bytes())
    raw[20] = (raw[20] + 1) % 256
    r1.write_bytes(raw)

    # Monitor downstream calls: neither extract_sra nor fastp nor STAR should be called
    def mock_extract_sra(*args, **kwargs):
        raise AssertionError("Unapproved SRA extraction was triggered!")

    def mock_fastp(*args, **kwargs):
        raise AssertionError("Unapproved fastp processing was triggered!")

    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.extract_sra", mock_extract_sra)
    monkeypatch.setattr("src.quantification_bridge.universal_pipeline.run_fastp_preprocessing", mock_fastp)

    # Record snapshot of existing files in sra_out_dir
    files_before = set(sra_out_dir.iterdir())
    assert len(files_before) == 2

    # Approve and run
    plan = pipeline.plan()
    pipeline.approve(plan.approval_token)
    res = pipeline.run()

    # Verify pipeline halted safely
    assert res.status == "BLOCKED_INVALID_SRA_CHECKPOINT"
    assert res.failed_stage == "sra_extraction"
    assert any(f.code == "INVALID_SRA_CHECKPOINT" for f in res.findings)

    # Verify existing files were NOT deleted
    files_after = set(sra_out_dir.iterdir())
    assert files_before == files_after, "Existing SRA files must not be deleted on invalid checkpoint"
