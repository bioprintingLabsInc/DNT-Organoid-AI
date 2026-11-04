"""Unit tests for CheckpointManager, capacity checks, and cryptographic file sync."""

import json
from pathlib import Path
import tempfile
import pytest

from src.quantification_bridge.checkpoints import (
    CheckpointManager,
    ChecksumVerificationError,
    InsufficientStorageCapacityError,
    ManifestVerificationError,
    UnverifiedStorageQuotaError,
    archive_and_sync_directory,
    check_storage_capacity,
    restore_and_verify_archive,
    sha256_file,
    sync_and_verify_directory,
    sync_and_verify_file,
)


def test_sha256_file():
    """Verify SHA-256 calculation against known content."""
    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "test.txt"
        test_file.write_text("DNT Organoid AI Checkpoint Verification\n", encoding="utf-8")
        h = sha256_file(test_file)
        assert len(h) == 64
        # Deterministic check
        assert h == sha256_file(test_file)


def test_check_storage_capacity_sufficient_and_shortfall():
    """Verify capacity verification passes with realistic demands and fails when capacity is exceeded."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_p = Path(tmpdir)

        # 1. Normal requirement should pass
        res = check_storage_capacity(tmp_p, required_bytes=1024 * 1024, safety_margin_gb=0.1)
        assert res["is_sufficient"] is True
        assert res["required_gb"] >= 0.0

        # 2. Absurdly large requirement (e.g. 100 PB) must raise InsufficientStorageCapacityError
        absurd_bytes = 100 * (1024**5)  # 100 PB
        with pytest.raises(InsufficientStorageCapacityError, match="Insufficient storage capacity"):
            check_storage_capacity(tmp_p, required_bytes=absurd_bytes, safety_margin_gb=1.0)


def test_sync_and_verify_file_success_and_cleanup():
    """Verify file synchronization calculates checksums and cleans source only when requested."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_p = Path(tmpdir)
        src = tmp_p / "src.dat"
        dst = tmp_p / "dst" / "synced.dat"

        payload = b"Sample genomic index chunk binary stream\x00\x01\x02\xff" * 1000
        src.write_bytes(payload)

        # Sync without source deletion
        res = sync_and_verify_file(src, dst, delete_source_after=False)
        assert res["status"] == "VERIFIED_MATCH"
        assert res["source_deleted"] is False
        assert src.exists()
        assert dst.exists()
        assert sha256_file(src) == sha256_file(dst) == res["sha256"]

        # Sync with source deletion
        dst2 = tmp_p / "dst" / "synced2.dat"
        res2 = sync_and_verify_file(src, dst2, delete_source_after=True)
        assert res2["status"] == "VERIFIED_MATCH"
        assert res2["source_deleted"] is True
        assert not src.exists()
        assert dst2.exists()


def test_sync_and_verify_file_checksum_corruption():
    """Verify that if expected checksum differs, ChecksumVerificationError is raised."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_p = Path(tmpdir)
        src = tmp_p / "src.dat"
        dst = tmp_p / "dst.dat"
        src.write_text("Hello World", encoding="utf-8")

        fake_sha = "0000000000000000000000000000000000000000000000000000000000000000"
        with pytest.raises(ChecksumVerificationError, match="failed expected checksum"):
            sync_and_verify_file(src, dst, expected_sha256=fake_sha)


def test_sync_and_verify_directory():
    """Verify recursive directory synchronization with nested checksums."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_p = Path(tmpdir)
        src_dir = tmp_p / "star_index"
        src_dir.mkdir()
        (src_dir / "chrName.txt").write_text("chr1\nchr2\n", encoding="utf-8")
        (src_dir / "chrLength.txt").write_text("1000\n2000\n", encoding="utf-8")
        sub = src_dir / "subdir"
        sub.mkdir()
        (sub / "table.tab").write_text("col1\tcol2\n", encoding="utf-8")

        dst_dir = tmp_p / "drive_backup" / "star_index"
        res = sync_and_verify_directory(src_dir, dst_dir, safety_margin_gb=0.1)

        assert res["status"] == "DIRECTORY_VERIFIED_SYNC"
        assert res["total_files"] == 3
        assert (dst_dir / "chrName.txt").is_file()
        assert (dst_dir / "chrLength.txt").is_file()
        assert (dst_dir / "subdir" / "table.tab").is_file()
        assert sha256_file(src_dir / "chrName.txt") == sha256_file(dst_dir / "chrName.txt")


def test_checkpoint_manager_lifecycle():
    """Verify CheckpointManager records stages and synchronizes between local and persistent stores."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_p = Path(tmpdir)
        local_dir = tmp_p / "local_checkpoints"
        persist_dir = tmp_p / "drive_checkpoints"

        # Initialize manager
        mgr = CheckpointManager(local_checkpoint_dir=local_dir, persistent_checkpoint_dir=persist_dir)
        assert mgr.is_stage_completed("REFERENCES_VERIFIED") is False

        # Record stage 1
        mgr.record_stage_completed(
            "REFERENCES_VERIFIED",
            artifacts=[{"name": "Homo_sapiens.GRCh38.112.gtf", "sha256": "8c87..."}],
            metrics={"duration_s": 12.3},
        )
        assert mgr.is_stage_completed("REFERENCES_VERIFIED") is True
        rec = mgr.get_stage_record("REFERENCES_VERIFIED")
        assert rec is not None
        assert rec["status"] == "COMPLETED"

        # Verify persistent store got updated
        persist_manifest = persist_dir / "checkpoint_manifest.json"
        assert persist_manifest.is_file()
        persist_data = json.loads(persist_manifest.read_text(encoding="utf-8"))
        assert "REFERENCES_VERIFIED" in persist_data["completed_stages"]

        # Simulate Colab restart with fresh local dir restoring from persistent store
        fresh_local_dir = tmp_p / "new_colab_local"
        mgr2 = CheckpointManager(local_checkpoint_dir=fresh_local_dir, persistent_checkpoint_dir=persist_dir)
        assert mgr2.is_stage_completed("REFERENCES_VERIFIED") is True
        assert (fresh_local_dir / "checkpoint_manifest.json").is_file()


def test_checkpoint_manager_graceful_without_persistent_dir():
    """Verify CheckpointManager functions normally as a purely local manager if no persistent dir is provided."""
    with tempfile.TemporaryDirectory() as tmpdir:
        local_dir = Path(tmpdir) / "local_only"
        mgr = CheckpointManager(local_checkpoint_dir=local_dir)
        mgr.record_stage_completed("LOCAL_ONLY_STAGE")
        assert mgr.is_stage_completed("LOCAL_ONLY_STAGE") is True
        assert (local_dir / "checkpoint_manifest.json").is_file()


def test_archive_and_sync_and_restore():
    """Verify single-archive synchronization and bit-for-bit file extraction verification."""
    from src.quantification_bridge.checkpoints import archive_and_sync_directory, restore_and_verify_archive

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_p = Path(tmpdir)
        src_dir = tmp_p / "star_index"
        src_dir.mkdir()
        (src_dir / "chrName.txt").write_text("chr1\nchr2\n", encoding="utf-8")
        (src_dir / "Genome").write_bytes(b"Simulated genome binary stream\x01\x02\x03" * 50)
        (src_dir / "SA").write_bytes(b"Simulated suffix array binary stream\x04\x05\x06" * 50)

        dst_dir = tmp_p / "drive_backup"

        # Archive and sync
        sync_res = archive_and_sync_directory(
            src_dir=src_dir,
            dst_dir=dst_dir,
            archive_name="test_index.tar",
            safety_margin_gb=0.1,
            account_quota_gb=5000.0,
            delete_source_after=False,
        )
        assert sync_res["status"] == "ARCHIVE_VERIFIED_SYNC"
        assert sync_res["total_files"] == 3
        assert Path(sync_res["dst_archive_path"]).is_file()
        assert Path(sync_res["dst_manifest_path"]).is_file()

        # Restore from archive into a clean directory
        restore_dir = tmp_p / "restored_index"
        rest_res = restore_and_verify_archive(
            archive_path=sync_res["dst_archive_path"],
            target_extract_dir=restore_dir,
        )
        assert rest_res["status"] == "ARCHIVE_EXTRACTED_AND_VERIFIED"
        assert rest_res["verified_files_count"] == 3
        assert (restore_dir / "chrName.txt").read_text(encoding="utf-8") == "chr1\nchr2\n"
        assert (restore_dir / "Genome").read_bytes() == (src_dir / "Genome").read_bytes()
        assert (restore_dir / "SA").read_bytes() == (src_dir / "SA").read_bytes()


def test_check_storage_capacity_gdrive_unverified_raises():
    """Verify that a Google Drive path without authorized API or manual verification raises UnverifiedStorageQuotaError."""
    fake_gdrive = Path("/content/drive/MyDrive/Test")
    with pytest.raises(UnverifiedStorageQuotaError, match="Virtual filesystem statistics cannot establish"):
        check_storage_capacity(fake_gdrive, required_bytes=1000)


def test_check_storage_capacity_gdrive_user_verified_passes():
    """Verify that a Google Drive path with user_verified_account_free_gb passes correctly."""
    fake_gdrive = Path("/content/drive/MyDrive/Test")
    # Provide explicit user confirmation (e.g. 100 GB free)
    res = check_storage_capacity(
        fake_gdrive,
        required_bytes=1024 * 1024,
        safety_margin_gb=1.0,
        user_verified_account_free_gb=100.0,
    )
    assert res["is_sufficient"] is True
    assert res["quota_source"] == "USER_VERIFIED_ACCOUNT_QUOTA"
    assert res["account_quota_info"]["user_verified_account_free_gb"] == 100.0


def test_restore_and_verify_archive_missing_manifest_raises():
    """Verify that restore_and_verify_archive rejects an archive lacking a companion manifest."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_p = Path(tmpdir)
        arc = tmp_p / "orphan_index.tar"
        import tarfile
        with tarfile.open(arc, "w") as tar:
            pass  # empty tar without manifest

        target = tmp_p / "extracted"
        with pytest.raises(ManifestVerificationError, match="Required companion manifest missing"):
            restore_and_verify_archive(arc, target, require_manifest=True)


def test_restore_and_verify_archive_malformed_manifest_raises():
    """Verify that restore_and_verify_archive rejects an archive with invalid or incomplete JSON manifest."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_p = Path(tmpdir)
        arc = tmp_p / "corrupt_manifest.tar"
        import tarfile
        with tarfile.open(arc, "w") as tar:
            pass

        manifest_p = tmp_p / "corrupt_manifest.tar.manifest.json"
        manifest_p.write_text("{invalid json", encoding="utf-8")

        target = tmp_p / "extracted"
        with pytest.raises(ManifestVerificationError, match="malformed or invalid JSON"):
            restore_and_verify_archive(arc, target, require_manifest=True)

