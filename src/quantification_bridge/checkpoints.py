"""Checkpoint and persistent storage synchronization for Quantification Bridge.

Supports Google Colab local disk execution with persistent Google Drive (or external POSIX) backup:
1. Mount and capacity preflight check (ensures target storage has sufficient free space + safety margin).
2. End-to-end cryptographic checksum verification (SHA-256) after transfer before any local cleanup.
3. Single-archive verified transfer (tar) to minimize Google Drive API calls and FUSE overhead.
4. Resilient checkpoint state machine permitting disconnect recovery.
5. Fully optional and decoupled: works seamlessly with or without Google Drive.
"""

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import shutil
import tarfile
import time
from typing import Any


class StorageError(Exception):
    """Base exception for storage and checkpoint operations."""
    pass


class InsufficientStorageCapacityError(StorageError):
    """Raised when destination storage lacks required free space plus safety margin."""
    pass


class UnverifiedStorageQuotaError(StorageError):
    """Raised when destination is a Google Drive path and account-level quota cannot be verified."""
    pass


class ManifestVerificationError(StorageError):
    """Raised when a companion archive manifest is missing, malformed, or incomplete."""
    pass


class ChecksumVerificationError(StorageError):
    """Raised when transferred file checksum does not match source checksum."""
    pass


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    """Compute SHA-256 hash of a file efficiently."""
    hasher = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def is_google_drive_path(p: Path) -> bool:
    """Detect if path is a Google Drive mount point or subfolder."""
    posix_str = p.as_posix().lower()
    return "/content/drive" in posix_str or "drive/mydrive" in posix_str or "/gdrive" in posix_str


def check_storage_capacity(
    target_path: str | Path,
    required_bytes: int,
    safety_margin_gb: float = 5.0,
    user_confirmed_account_free_gb: float | None = None,
    user_verified_account_free_gb: float | None = None,
    account_quota_gb: float | None = None,
    google_drive_service: Any = None,
) -> dict[str, Any]:
    """Verify target storage filesystem has sufficient capacity before beginning a transfer.

    For Google Drive FUSE destinations:
    FUSE filesystem statistics (statvfs) report virtual disk metrics that do not establish
    true account-level quota. Therefore, an account-level quota MUST be verified either:
    1. Via authorized Google Drive API client (google_drive_service) querying about().get(fields='storageQuota').
    2. Via explicit user confirmation of account headroom (user_confirmed_account_free_gb) verified at
       https://drive.google.com/drive/quota.

    Folder-based file scanning is strictly prohibited from being reported as verified account headroom.
    """
    p = Path(target_path).resolve()
    # Find existing parent if path does not exist yet
    probe_path = p
    while not probe_path.exists() and probe_path.parent != probe_path:
        probe_path = probe_path.parent

    if not probe_path.exists():
        raise FileNotFoundError(f"Cannot determine filesystem capacity for non-existent path: {target_path}")

    total, used, free = shutil.disk_usage(probe_path)
    safety_margin_bytes = int(safety_margin_gb * (1024**3))
    total_needed_bytes = required_bytes + safety_margin_bytes

    effective_free = free
    quota_source = "POSIX_FILESYSTEM_STATVFS"
    account_quota_info = None

    is_gdrive = is_google_drive_path(p) or is_google_drive_path(probe_path) or is_google_drive_path(Path(target_path))
    confirmed_user_free = user_confirmed_account_free_gb if user_confirmed_account_free_gb is not None else user_verified_account_free_gb

    if google_drive_service is not None:
        # 1. Query authorized Google Drive API metadata
        try:
            about = google_drive_service.about().get(fields="storageQuota").execute()
            sq = about.get("storageQuota", {})
            limit_bytes = int(sq.get("limit", 0))
            usage_bytes = int(sq.get("usage", 0))
            if limit_bytes > 0:
                drive_free_bytes = max(0, limit_bytes - usage_bytes)
                effective_free = min(free, drive_free_bytes)
                quota_source = "AUTHORIZED_GOOGLE_DRIVE_API"
                account_quota_info = {
                    "source": quota_source,
                    "limit_gb": round(limit_bytes / (1024**3), 2),
                    "usage_gb": round(usage_bytes / (1024**3), 2),
                    "account_free_gb": round(drive_free_bytes / (1024**3), 2),
                    "verified": True,
                }
        except Exception as e:
            raise StorageError(f"Failed to query authorized Google Drive API storageQuota: {e}")
    elif confirmed_user_free is not None:
        # 2. Explicit manual user verification
        account_free_bytes = int(confirmed_user_free * (1024**3))
        effective_free = min(free, account_free_bytes)
        quota_source = "USER_VERIFIED_ACCOUNT_QUOTA"
        account_quota_info = {
            "source": quota_source,
            "user_confirmed_account_free_gb": confirmed_user_free,
            "user_verified_account_free_gb": confirmed_user_free,
            "effective_free_gb": round(effective_free / (1024**3), 2),
            "verified": True,
        }
    elif is_gdrive:
        # 3. Google Drive path without authorized API or manual verification
        raise UnverifiedStorageQuotaError(
            f"Target path '{target_path}' is on Google Drive FUSE storage. "
            "Virtual filesystem statistics cannot establish total Google Drive account usage. "
            "You must either provide an authorized google_drive_service client or supply "
            "user_confirmed_account_free_gb (manually confirmed via https://drive.google.com/drive/quota). "
            "Estimating account quota from scanned folder contents or hardcoded assumptions is prohibited."
        )
    elif account_quota_gb is not None and account_quota_gb > 0:
        # Local non-Drive path with specified quota limit
        account_quota_bytes = int(account_quota_gb * (1024**3))
        effective_free = min(free, account_quota_bytes)
        quota_source = "LOCAL_QUOTA_OVERRIDE"
        account_quota_info = {
            "source": quota_source,
            "account_quota_gb": account_quota_gb,
            "effective_free_gb": round(effective_free / (1024**3), 2),
            "verified": True,
        }
    else:
        account_quota_info = {
            "source": quota_source,
            "total_gb": round(total / (1024**3), 2),
            "free_gb": round(free / (1024**3), 2),
            "verified": True,
        }

    metrics = {
        "target_path": str(p),
        "probe_path": str(probe_path),
        "quota_source": quota_source,
        "total_gb": round(total / (1024**3), 2),
        "used_gb": round(used / (1024**3), 2),
        "free_gb": round(effective_free / (1024**3), 2),
        "required_gb": round(required_bytes / (1024**3), 2),
        "safety_margin_gb": safety_margin_gb,
        "total_needed_gb": round(total_needed_bytes / (1024**3), 2),
        "is_sufficient": effective_free >= total_needed_bytes,
        "account_quota_info": account_quota_info,
    }

    if not metrics["is_sufficient"]:
        shortfall_gb = round((total_needed_bytes - effective_free) / (1024**3), 2)
        raise InsufficientStorageCapacityError(
            f"Insufficient storage capacity on {probe_path}! "
            f"Available: {metrics['free_gb']} GB, Needed (including {safety_margin_gb} GB safety margin): "
            f"{metrics['total_needed_gb']} GB. Shortfall: {shortfall_gb} GB."
        )

    return metrics


def sync_and_verify_file(
    src_file: str | Path,
    dst_file: str | Path,
    expected_sha256: str | None = None,
    delete_source_after: bool = False,
    chunk_size: int = 1024 * 1024,
    safety_margin_gb: float = 1.0,
    user_confirmed_account_free_gb: float | None = None,
    user_verified_account_free_gb: float | None = None,
    account_quota_gb: float | None = None,
    google_drive_service: Any = None,
) -> dict[str, Any]:
    """Copy a file to destination and verify cryptographic SHA-256 before optional source cleanup."""
    s = Path(src_file).resolve()
    d = Path(dst_file).resolve()

    if not s.is_file():
        raise FileNotFoundError(f"Source file not found: {s}")

    file_size = s.stat().st_size
    # Preflight destination storage capacity
    check_storage_capacity(
        d.parent,
        required_bytes=file_size,
        safety_margin_gb=safety_margin_gb,
        user_confirmed_account_free_gb=user_confirmed_account_free_gb,
        user_verified_account_free_gb=user_verified_account_free_gb,
        account_quota_gb=account_quota_gb,
        google_drive_service=google_drive_service,
    )

    d.parent.mkdir(parents=True, exist_ok=True)

    # 1. Source SHA-256
    src_sha = sha256_file(s, chunk_size=chunk_size)
    if expected_sha256 is not None and src_sha.lower() != expected_sha256.lower():
        raise ChecksumVerificationError(
            f"Source file {s.name} failed expected checksum! Expected {expected_sha256}, got {src_sha}"
        )

    # 2. Copy file
    t0 = time.time()
    shutil.copy2(s, d)
    transfer_time = time.time() - t0

    # 3. Destination SHA-256
    dst_sha = sha256_file(d, chunk_size=chunk_size)
    if dst_sha != src_sha:
        d.unlink(missing_ok=True)
        raise ChecksumVerificationError(
            f"Transfer corruption detected for {s.name}! Source SHA-256: {src_sha}, Destination SHA-256: {dst_sha}"
        )

    # 4. Optional source cleanup ONLY after destination verification
    if delete_source_after:
        s.unlink()

    return {
        "status": "VERIFIED_MATCH",
        "filename": s.name,
        "src_path": str(s),
        "dst_path": str(d),
        "size_bytes": file_size,
        "sha256": src_sha,
        "transfer_seconds": round(transfer_time, 2),
        "source_deleted": delete_source_after,
    }


def sync_and_verify_directory(
    src_dir: str | Path,
    dst_dir: str | Path,
    safety_margin_gb: float = 5.0,
    delete_source_after: bool = False,
    user_confirmed_account_free_gb: float | None = None,
    user_verified_account_free_gb: float | None = None,
    account_quota_gb: float | None = None,
    google_drive_service: Any = None,
) -> dict[str, Any]:
    """Recursively copy a directory and verify all constituent file SHA-256 checksums."""
    s_d = Path(src_dir).resolve()
    d_d = Path(dst_dir).resolve()

    if not s_d.is_dir():
        raise NotADirectoryError(f"Source directory not found: {s_d}")

    file_list: list[Path] = [p for p in s_d.rglob("*") if p.is_file()]
    total_bytes = sum(f.stat().st_size for f in file_list)

    check_storage_capacity(
        d_d,
        required_bytes=total_bytes,
        safety_margin_gb=safety_margin_gb,
        user_confirmed_account_free_gb=user_confirmed_account_free_gb,
        user_verified_account_free_gb=user_verified_account_free_gb,
        account_quota_gb=account_quota_gb,
        google_drive_service=google_drive_service,
    )

    d_d.mkdir(parents=True, exist_ok=True)

    verified_files = []
    t0 = time.time()
    for src_f in file_list:
        rel_path = src_f.relative_to(s_d)
        dst_f = d_d / rel_path
        res = sync_and_verify_file(
            src_f,
            dst_f,
            delete_source_after=False,
            user_confirmed_account_free_gb=user_confirmed_account_free_gb,
            user_verified_account_free_gb=user_verified_account_free_gb,
            account_quota_gb=account_quota_gb,
            google_drive_service=google_drive_service,
        )
        verified_files.append(res)

    total_time = time.time() - t0

    if delete_source_after:
        shutil.rmtree(s_d)

    return {
        "status": "DIRECTORY_VERIFIED_SYNC",
        "src_dir": str(s_d),
        "dst_dir": str(d_d),
        "total_files": len(verified_files),
        "total_bytes": total_bytes,
        "total_seconds": round(total_time, 2),
        "files": verified_files,
        "source_deleted": delete_source_after,
    }


def archive_and_sync_directory(
    src_dir: str | Path,
    dst_dir: str | Path,
    archive_name: str = "star_index_GRCh38_112_dense.tar",
    safety_margin_gb: float = 10.0,
    user_confirmed_account_free_gb: float | None = None,
    user_verified_account_free_gb: float | None = None,
    account_quota_gb: float | None = None,
    google_drive_service: Any = None,
    delete_source_after: bool = False,
) -> dict[str, Any]:
    """Archive a directory into a single verified tarball and sync to persistent storage with manifest."""
    s_d = Path(src_dir).resolve()
    d_d = Path(dst_dir).resolve()

    if not s_d.is_dir():
        raise NotADirectoryError(f"Source directory not found: {s_d}")

    # 1. Build file-level inventory and checksums
    file_list = sorted([p for p in s_d.rglob("*") if p.is_file()])
    manifest_files = {}
    total_uncompressed_bytes = 0
    for f in file_list:
        rel = f.relative_to(s_d).as_posix()
        sz = f.stat().st_size
        total_uncompressed_bytes += sz
        manifest_files[rel] = {
            "size_bytes": sz,
            "sha256": sha256_file(f),
        }

    # 2. Create archive locally
    local_archive = s_d.parent / archive_name
    print(f"Creating local archive: {local_archive.name} ({len(file_list)} files)...")
    t0 = time.time()
    with tarfile.open(local_archive, "w") as tar:
        for f in file_list:
            arcname = f.relative_to(s_d).as_posix()
            tar.add(f, arcname=arcname)
    archive_size_bytes = local_archive.stat().st_size
    archive_sha256 = sha256_file(local_archive)

    # 3. Create manifest json
    manifest_data = {
        "archive_name": archive_name,
        "archive_size_bytes": archive_size_bytes,
        "archive_sha256": archive_sha256,
        "total_files": len(file_list),
        "total_uncompressed_bytes": total_uncompressed_bytes,
        "created_at_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "files": manifest_files,
    }
    local_manifest = s_d.parent / f"{archive_name}.manifest.json"
    local_manifest.write_text(json.dumps(manifest_data, indent=2), encoding="utf-8")
    manifest_sha256 = sha256_file(local_manifest)

    # 4. Check destination storage capacity
    check_storage_capacity(
        d_d,
        required_bytes=archive_size_bytes + local_manifest.stat().st_size,
        safety_margin_gb=safety_margin_gb,
        user_confirmed_account_free_gb=user_confirmed_account_free_gb,
        user_verified_account_free_gb=user_verified_account_free_gb,
        account_quota_gb=account_quota_gb,
        google_drive_service=google_drive_service,
    )

    # 5. Sync archive and manifest to destination with SHA-256 validation
    d_d.mkdir(parents=True, exist_ok=True)
    dst_archive = d_d / archive_name
    dst_manifest = d_d / f"{archive_name}.manifest.json"

    res_archive = sync_and_verify_file(
        local_archive,
        dst_archive,
        expected_sha256=archive_sha256,
        delete_source_after=True,
        user_confirmed_account_free_gb=user_confirmed_account_free_gb,
        user_verified_account_free_gb=user_verified_account_free_gb,
        account_quota_gb=account_quota_gb,
        google_drive_service=google_drive_service,
    )
    res_manifest = sync_and_verify_file(
        local_manifest,
        dst_manifest,
        expected_sha256=manifest_sha256,
        delete_source_after=True,
        user_confirmed_account_free_gb=user_confirmed_account_free_gb,
        user_verified_account_free_gb=user_verified_account_free_gb,
        account_quota_gb=account_quota_gb,
        google_drive_service=google_drive_service,
    )

    if delete_source_after:
        shutil.rmtree(s_d)

    total_time = time.time() - t0
    return {
        "status": "ARCHIVE_VERIFIED_SYNC",
        "archive_name": archive_name,
        "archive_sha256": archive_sha256,
        "archive_size_gb": round(archive_size_bytes / (1024**3), 2),
        "total_files": len(file_list),
        "total_uncompressed_gb": round(total_uncompressed_bytes / (1024**3), 2),
        "total_seconds": round(total_time, 2),
        "dst_archive_path": str(dst_archive),
        "dst_manifest_path": str(dst_manifest),
    }


def restore_and_verify_archive(
    archive_path: str | Path,
    target_extract_dir: str | Path,
    expected_sha256: str | None = None,
    require_manifest: bool = True,
) -> dict[str, Any]:
    """Extract a verified tar archive into target directory and verify constituent files against companion manifest."""
    arc_p = Path(archive_path).resolve()
    tgt_p = Path(target_extract_dir).resolve()

    if not arc_p.is_file():
        raise FileNotFoundError(f"Archive not found: {arc_p}")

    # 1. Require and validate companion manifest
    manifest_p = arc_p.parent / f"{arc_p.name}.manifest.json"
    if not manifest_p.is_file():
        if require_manifest:
            raise ManifestVerificationError(
                f"Required companion manifest missing: {manifest_p}. "
                "An archive cannot be classified as verified without its companion cryptographic manifest."
            )
        manifest = None
    else:
        try:
            manifest = json.loads(manifest_p.read_text(encoding="utf-8"))
        except Exception as e:
            raise ManifestVerificationError(f"Companion manifest is malformed or invalid JSON ({manifest_p}): {e}")

        if not isinstance(manifest, dict) or "files" not in manifest or not manifest.get("files"):
            raise ManifestVerificationError(f"Companion manifest ({manifest_p}) is incomplete: missing 'files' inventory.")

        if expected_sha256 is None:
            expected_sha256 = manifest.get("archive_sha256")

    if require_manifest and not manifest:
        raise ManifestVerificationError(f"Valid companion manifest is required for archive verification: {manifest_p}")

    # 2. Verify archive checksum
    actual_sha = sha256_file(arc_p)
    if expected_sha256 and actual_sha.lower() != expected_sha256.lower():
        raise ChecksumVerificationError(
            f"Archive checksum mismatch for {arc_p.name}! Expected {expected_sha256}, got {actual_sha}"
        )

    # 3. Extract to target directory safely
    tgt_p.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with tarfile.open(arc_p, "r") as tar:
        if hasattr(tarfile, "data_filter"):
            tar.extractall(tgt_p, filter="data")
        else:
            tar.extractall(tgt_p)
    extract_time = time.time() - t0

    # 4. Strictly verify constituent files
    verified_files = 0
    if manifest and "files" in manifest:
        for rel_name, f_info in manifest["files"].items():
            extracted_f = tgt_p / rel_name
            if not extracted_f.is_file() or extracted_f.stat().st_size != f_info["size_bytes"]:
                raise ChecksumVerificationError(
                    f"Extracted file missing or size mismatch: {rel_name} (expected {f_info['size_bytes']} bytes)"
                )
            actual_f_sha = sha256_file(extracted_f)
            if actual_f_sha != f_info["sha256"]:
                raise ChecksumVerificationError(
                    f"Extracted file checksum mismatch for {rel_name}! Expected {f_info['sha256']}, got {actual_f_sha}"
                )
            verified_files += 1

    if require_manifest and verified_files == 0:
        raise ManifestVerificationError("Zero constituent files verified from companion manifest!")

    return {
        "status": "ARCHIVE_EXTRACTED_AND_VERIFIED",
        "archive_path": str(arc_p),
        "target_dir": str(tgt_p),
        "extract_seconds": round(extract_time, 2),
        "verified_files_count": verified_files,
    }


class CheckpointManager:
    """Stateful checkpoint manager supporting Colab local execution and persistent backup."""

    def __init__(
        self,
        local_checkpoint_dir: str | Path,
        persistent_checkpoint_dir: str | Path | None = None,
    ) -> None:
        self.local_dir = Path(local_checkpoint_dir).resolve()
        self.local_dir.mkdir(parents=True, exist_ok=True)
        self.local_manifest_file = self.local_dir / "checkpoint_manifest.json"

        self.persistent_dir: Path | None = None
        self.persistent_manifest_file: Path | None = None
        if persistent_checkpoint_dir:
            p_p = Path(persistent_checkpoint_dir).resolve()
            try:
                p_p.mkdir(parents=True, exist_ok=True)
                self.persistent_dir = p_p
                self.persistent_manifest_file = p_p / "checkpoint_manifest.json"
            except Exception as e:
                print(f"[WARNING] Could not initialize persistent checkpoint dir {p_p}: {e}")

        self.manifest = self._load_manifest()

    def _load_manifest(self) -> dict[str, Any]:
        """Load manifest from persistent store if available, else local store, else new."""
        if self.persistent_manifest_file and self.persistent_manifest_file.is_file():
            try:
                data = json.loads(self.persistent_manifest_file.read_text(encoding="utf-8"))
                self.local_manifest_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
                return data
            except Exception as e:
                print(f"[WARNING] Failed to read persistent manifest: {e}")

        if self.local_manifest_file.is_file():
            try:
                return json.loads(self.local_manifest_file.read_text(encoding="utf-8"))
            except Exception:
                pass

        return {
            "version": "1.0.0",
            "created_at_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
            "completed_stages": {},
        }

    def _save_manifest(self) -> None:
        """Persist manifest atomically to both local and persistent stores."""
        content = json.dumps(self.manifest, indent=2)
        tmp_local = self.local_dir / f"tmp_manifest_{os.getpid()}.json"
        tmp_local.write_text(content, encoding="utf-8")
        tmp_local.replace(self.local_manifest_file)

        if self.persistent_dir and self.persistent_manifest_file:
            try:
                tmp_persist = self.persistent_dir / f"tmp_manifest_{os.getpid()}.json"
                tmp_persist.write_text(content, encoding="utf-8")
                tmp_persist.replace(self.persistent_manifest_file)
            except Exception as e:
                print(f"[WARNING] Failed to write persistent manifest to {self.persistent_dir}: {e}")

    def is_stage_completed(self, stage_name: str) -> bool:
        """Query if a stage has been completed and marked valid."""
        return stage_name in self.manifest.get("completed_stages", {})

    def get_stage_record(self, stage_name: str) -> dict[str, Any] | None:
        """Retrieve stored metadata record for a completed stage."""
        return self.manifest.get("completed_stages", {}).get(stage_name)

    def record_stage_completed(
        self,
        stage_name: str,
        artifacts: list[dict[str, Any]] | None = None,
        metrics: dict[str, Any] | None = None,
    ) -> None:
        """Record stage completion with timestamp, validated artifacts, and execution metrics."""
        record = {
            "stage_name": stage_name,
            "status": "COMPLETED",
            "completed_at_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
            "artifacts": artifacts or [],
            "metrics": metrics or {},
        }
        self.manifest.setdefault("completed_stages", {})[stage_name] = record
        self._save_manifest()
        print(f"[CHECKPOINT] Stage '{stage_name}' recorded and persisted successfully.")
