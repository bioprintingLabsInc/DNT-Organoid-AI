"""Verified STAR archive persistence using Google Drive's resumable upload API.

This deliberately avoids copying large files through Colab's Drive FUSE mount.
Only the API-reported size and MD5 are used to verify the remote blob; the
companion manifest records SHA-256 for the archive and individual index files.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import time

from googleapiclient.http import MediaFileUpload


_CHUNK_SIZE = 16 * 1024 * 1024
_FOLDER_MIME = "application/vnd.google-apps.folder"


def _digest(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _list_named(service, parent: str, name: str, folder_only: bool = False) -> list[dict]:
    escaped = name.replace("\\", "\\\\").replace("'", "\\'")
    query = f"'{parent}' in parents and name = '{escaped}' and trashed = false"
    if folder_only:
        query += f" and mimeType = '{_FOLDER_MIME}'"
    files = []
    token = None
    while True:
        result = service.files().list(
            q=query, fields="nextPageToken,files(id,name,mimeType,size,md5Checksum)",
            pageSize=100, pageToken=token,
        ).execute()
        files.extend(result.get("files", []))
        token = result.get("nextPageToken")
        if not token:
            break
    return files


def _ensure_folder(service, parent: str, name: str) -> str:
    matches = _list_named(service, parent, name)
    if len(matches) > 1:
        raise RuntimeError(f"Ambiguous Drive folder name {name!r}: multiple matches")
    if matches:
        if matches[0]["mimeType"] != _FOLDER_MIME:
            raise RuntimeError(f"Drive item {name!r} exists but is not a folder")
        return matches[0]["id"]
    created = service.files().create(
        body={"name": name, "mimeType": _FOLDER_MIME, "parents": [parent]},
        fields="id",
    ).execute()
    return created["id"]


def _check_quota(service, needed_bytes: int, safety_margin_gb: float) -> None:
    quota = service.about().get(fields="storageQuota(limit,usage)").execute().get("storageQuota", {})
    if quota.get("limit"):
        free_bytes = int(quota["limit"]) - int(quota.get("usage", 0))
        minimum = needed_bytes + int(safety_margin_gb * 1024 ** 3)
        if free_bytes < minimum:
            raise RuntimeError(
                f"Insufficient Google Drive quota: {free_bytes / 1024**3:.2f} GiB free; "
                f"{minimum / 1024**3:.2f} GiB required including margin"
            )
    else:
        raise RuntimeError("Cannot verify Google Drive quota through the API; upload not started")


def _upload_verified(service, local_path: Path, folder_id: str) -> dict:
    size = local_path.stat().st_size
    local_md5 = _digest(local_path, "md5")
    matches = _list_named(service, folder_id, local_path.name)
    if len(matches) > 1:
        raise RuntimeError(f"Duplicate remote filenames: {local_path.name}")
    if matches:
        match = matches[0]
        if str(match.get("size", "")) != str(size) or match.get("md5Checksum") != local_md5:
            raise RuntimeError(
                f"Remote file {local_path.name} exists with a mismatched size or MD5. "
                "Do not use it as a verified backup; resolve the conflicting file first."
            )
        print(f"[VERIFIED REUSE] {local_path.name}: Drive size and MD5 match")
        return match

    media = MediaFileUpload(
        str(local_path),
        mimetype="application/x-tar" if local_path.suffix == ".tar" else "application/json",
        chunksize=_CHUNK_SIZE, resumable=True,
    )
    request = service.files().create(
        body={"name": local_path.name, "parents": [folder_id]},
        media_body=media, fields="id,name,size,md5Checksum",
    )
    response = None
    last_percent = -10
    while response is None:
        progress, response = request.next_chunk(num_retries=6)
        if progress:
            percent = int(progress.progress() * 100)
            if percent >= last_percent + 10:
                print(f"[UPLOAD] {local_path.name}: {percent}%")
                last_percent = percent
    metadata = service.files().get(
        fileId=response["id"], fields="id,name,size,md5Checksum",
    ).execute()
    if str(metadata.get("size")) != str(size) or metadata.get("md5Checksum") != local_md5:
        raise RuntimeError(
            f"Remote integrity verification failed for {local_path.name}. "
            "Local archive retained; do not mark Stage 9 completed."
        )
    print(f"[REMOTE VERIFIED] {local_path.name}: size and MD5 match")
    return metadata


def archive_and_upload_directory(
    src_dir: str | Path,
    drive_service,
    archive_name: str = "star_index_GRCh38_112_dense.tar",
    drive_folder_path: tuple[str, ...] = ("DNT_Organoid_AI_Persistence", "references"),
    safety_margin_gb: float = 10.0,
) -> dict:
    """Create/reuse a local archive, upload resumably, verify remote MD5/size.

    The original STAR index, local archive, and manifest are deliberately kept.
    No completion checkpoint is written by this function.
    """
    if drive_service is None:
        raise RuntimeError("An authorized Google Drive API client is required")
    src_dir = Path(src_dir).resolve()
    if not src_dir.is_dir():
        raise FileNotFoundError(f"Index directory missing: {src_dir}")
    if not (src_dir / "SA").is_file() or not (src_dir / "Genome").is_file():
        raise RuntimeError("Index has not passed basic artifact presence checks")

    local_archive = src_dir.parent / archive_name
    local_manifest = src_dir.parent / (archive_name + ".manifest.json")

    can_reuse = False
    if local_archive.is_file() and local_manifest.is_file():
        try:
            previous = json.loads(local_manifest.read_text(encoding="utf-8"))
            can_reuse = (
                previous["archive_size_bytes"] == local_archive.stat().st_size
                and previous["archive_sha256"] == _digest(local_archive, "sha256")
            )
        except (OSError, KeyError, ValueError, json.JSONDecodeError):
            can_reuse = False

    if not can_reuse:
        files = sorted(p for p in src_dir.rglob("*") if p.is_file())
        if not files:
            raise RuntimeError("No files found in index directory")
        inventory = {
            p.relative_to(src_dir).as_posix(): {
                "size_bytes": p.stat().st_size, "sha256": _digest(p, "sha256")
            } for p in files
        }
        estimated = sum(x["size_bytes"] for x in inventory.values())
        disk_free = shutil.disk_usage(src_dir.parent).free
        if disk_free < estimated + 10 * 1024**3:
            raise RuntimeError("Not enough local scratch storage to safely build the index archive")
        temp_archive = local_archive.with_name(local_archive.name + ".building")
        print(f"Creating local archive from {len(files)} index files...")
        with tarfile.open(temp_archive, "w") as tar:
            for p in files:
                tar.add(p, arcname=p.relative_to(src_dir).as_posix())
        temp_archive.replace(local_archive)
        manifest = {
            "archive_name": archive_name,
            "archive_size_bytes": local_archive.stat().st_size,
            "archive_sha256": _digest(local_archive, "sha256"),
            "total_files": len(files),
            "files": inventory,
            "created_at_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        }
        local_manifest.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    else:
        manifest = previous
        print("[REUSE] Existing local archive matches its SHA-256 manifest")

    _check_quota(
        drive_service,
        local_archive.stat().st_size + local_manifest.stat().st_size,
        safety_margin_gb,
    )

    folder_id = "root"
    for folder_name in drive_folder_path:
        folder_id = _ensure_folder(drive_service, folder_id, folder_name)

    print("Uploading STAR index archive through resumable Google Drive API...")
    remote_archive = _upload_verified(drive_service, local_archive, folder_id)
    remote_manifest = _upload_verified(drive_service, local_manifest, folder_id)
    return {
        "status": "ARCHIVE_VERIFIED_SYNC",
        "total_files": manifest["total_files"],
        "archive_size_gb": round(local_archive.stat().st_size / 1024**3, 2),
        "archive_sha256": manifest["archive_sha256"],
        "dst_archive_path": f"gdrive-file-id:{remote_archive['id']}",
        "dst_manifest_path": f"gdrive-file-id:{remote_manifest['id']}",
        "archive_drive_file_id": remote_archive["id"],
        "manifest_drive_file_id": remote_manifest["id"],
        "remote_validation": "Drive API size and md5Checksum matched local files",
        "local_archive_path": str(local_archive),
        "local_manifest_path": str(local_manifest),
    }
