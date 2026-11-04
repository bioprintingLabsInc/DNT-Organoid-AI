"""Checkpoint and persistent storage synchronization for Quantification Bridge.

Supports Google Colab local disk execution with persistent Google Drive (or external POSIX) backup:
1. Mount and capacity preflight check (ensures target storage has sufficient free space + safety margin).
2. End-to-end cryptographic checksum verification (SHA-256) after transfer before any local cleanup.
3. Resilient checkpoint state machine permitting disconnect recovery.
4. Fully optional and decoupled: works seamlessly with or without Google Drive.
"""

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import shutil
import time
from typing import Any


class StorageError(Exception):
    """Base exception for storage and checkpoint operations."""
    pass


class InsufficientStorageCapacityError(StorageError):
    """Raised when destination storage lacks required free space plus safety margin."""
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


def check_storage_capacity(
    target_path: str | Path,
    required_bytes: int,
    safety_margin_gb: float = 5.0,
) -> dict[str, Any]:
    """Verify target storage filesystem has sufficient capacity before beginning a transfer."""
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

    metrics = {
        "target_path": str(p),
        "probe_path": str(probe_path),
        "total_gb": round(total / (1024**3), 2),
        "used_gb": round(used / (1024**3), 2),
        "free_gb": round(free / (1024**3), 2),
        "required_gb": round(required_bytes / (1024**3), 2),
        "safety_margin_gb": safety_margin_gb,
        "total_needed_gb": round(total_needed_bytes / (1024**3), 2),
        "is_sufficient": free >= total_needed_bytes,
    }

    if not metrics["is_sufficient"]:
        shortfall_gb = round((total_needed_bytes - free) / (1024**3), 2)
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
) -> dict[str, Any]:
    """Copy a file to destination and verify cryptographic SHA-256 before optional source cleanup."""
    s = Path(src_file).resolve()
    d = Path(dst_file).resolve()

    if not s.is_file():
        raise FileNotFoundError(f"Source file not found: {s}")

    file_size = s.stat().st_size
    # Preflight destination storage capacity
    check_storage_capacity(d.parent, required_bytes=file_size, safety_margin_gb=1.0)

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
        # Transfer corruption: remove corrupted destination file
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
) -> dict[str, Any]:
    """Recursively copy a directory and verify all constituent file SHA-256 checksums."""
    s_d = Path(src_dir).resolve()
    d_d = Path(dst_dir).resolve()

    if not s_d.is_dir():
        raise NotADirectoryError(f"Source directory not found: {s_d}")

    # Calculate total size of source directory
    file_list: list[Path] = [p for p in s_d.rglob("*") if p.is_file()]
    total_bytes = sum(f.stat().st_size for f in file_list)

    # Check capacity before copying
    check_storage_capacity(d_d, required_bytes=total_bytes, safety_margin_gb=safety_margin_gb)

    d_d.mkdir(parents=True, exist_ok=True)

    verified_files = []
    t0 = time.time()
    for src_f in file_list:
        rel_path = src_f.relative_to(s_d)
        dst_f = d_d / rel_path
        res = sync_and_verify_file(src_f, dst_f, delete_source_after=False)
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
        # Try persistent store first
        if self.persistent_manifest_file and self.persistent_manifest_file.is_file():
            try:
                data = json.loads(self.persistent_manifest_file.read_text(encoding="utf-8"))
                # Mirror to local
                self.local_manifest_file.write_text(json.dumps(data, indent=2), encoding="utf-8")
                return data
            except Exception as e:
                print(f"[WARNING] Failed to read persistent manifest: {e}")

        # Try local store
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
        # Write local
        tmp_local = self.local_dir / f"tmp_manifest_{os.getpid()}.json"
        tmp_local.write_text(content, encoding="utf-8")
        tmp_local.replace(self.local_manifest_file)

        # Write persistent if configured
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
