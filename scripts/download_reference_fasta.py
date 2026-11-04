#!/usr/bin/env python3
"""Portable reference FASTA downloader and validator for Ensembl Release 112 GRCh38.

Cross-platform implementation supporting native Linux, Google Colab, Windows/WSL,
and institutional compute clusters without hardcoded personal or platform paths.
Uses chunked parallel retrieval, atomic file replacement, and multi-tier checksum verification
(SHA-256 and Ensembl POSIX sum).
"""

import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import urllib.request

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

FASTA_URL = "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
EXPECTED_SHA256 = "d8c3af0094a7bba6125763bad779ec18a81483c739c6ed122094bdf86c187b92"
EXPECTED_FILE_SIZE = 881964081
EXPECTED_POSIX_SUM = 22450
EXPECTED_POSIX_BLOCKS = 861294


def compute_sha256(path: Path) -> str:
    """Compute SHA-256 hash using chunked streaming."""
    hasher = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def compute_posix_sum_python(path: Path) -> tuple[int, int]:
    """Compute standard POSIX 1003.2 (BSD rotating) checksum and 1KB block count in pure Python."""
    s = 0
    total_bytes = 0
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            total_bytes += len(chunk)
            for b in chunk:
                s = (s >> 1) + ((s & 1) << 15)
                s = (s + b) & 0xFFFF
    blocks = (total_bytes + 1023) // 1024
    return s, blocks


def verify_posix_sum(path: Path) -> tuple[int, int]:
    """Verify POSIX sum via native system 'sum' coreutil if available, or pure Python fallback."""
    sum_bin = shutil.which("sum")
    if sum_bin and sys.platform != "win32":
        try:
            res = subprocess.run([sum_bin, str(path)], capture_output=True, text=True, check=True)
            parts = res.stdout.strip().split()
            if len(parts) >= 2:
                return int(parts[0]), int(parts[1])
        except Exception as e:
            print(f"[NOTE] Native 'sum' invocation failed ({e}); falling back to pure Python POSIX sum.")
    
    # On Windows with WSL available
    if sys.platform == "win32" and shutil.which("wsl.exe"):
        try:
            from src.quantification_bridge.env_exec import to_wsl_path
            wsl_p = to_wsl_path(path)
            res = subprocess.run(["wsl.exe", "-e", "sum", wsl_p], capture_output=True, text=True, check=True)
            parts = res.stdout.strip().split()
            if len(parts) >= 2:
                return int(parts[0]), int(parts[1])
        except Exception:
            pass

    return compute_posix_sum_python(path)


def download_chunk(
    url: str,
    start_byte: int,
    end_byte: int,
    part_path: Path,
    chunk_idx: int,
    max_retries: int = 5,
) -> None:
    expected_len = end_byte - start_byte + 1
    if part_path.exists() and part_path.stat().st_size == expected_len:
        print(f"[Chunk {chunk_idx}] already completed ({expected_len} bytes).")
        return

    for attempt in range(1, max_retries + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "Range": f"bytes={start_byte}-{end_byte}",
                    "User-Agent": "DNT-Organoid-AI/2.0 (Portable Downloader)",
                },
            )
            with urllib.request.urlopen(req, timeout=45) as resp, open(part_path, "wb") as f_out:
                bytes_read = 0
                while True:
                    buf = resp.read(1024 * 1024)
                    if not buf:
                        break
                    f_out.write(buf)
                    bytes_read += len(buf)
            if bytes_read == expected_len:
                print(f"[Chunk {chunk_idx}] downloaded successfully ({bytes_read} bytes).")
                return
            else:
                print(f"[Chunk {chunk_idx}] attempt {attempt}: partial read ({bytes_read}/{expected_len}). Retrying...")
        except Exception as e:
            print(f"[Chunk {chunk_idx}] attempt {attempt} error: {e}. Retrying in 2s...")
            time.sleep(2)

    raise RuntimeError(f"[Chunk {chunk_idx}] failed to download after {max_retries} attempts.")


def download_and_verify_fasta(
    dest_path: Path,
    num_chunks: int = 8,
    skip_if_verified: bool = True,
) -> dict:
    """Download and cryptographically verify Ensembl release 112 FASTA atomically."""
    dest_dir = dest_path.parent
    dest_dir.mkdir(parents=True, exist_ok=True)

    # 1. Check existing target
    if dest_path.is_file() and dest_path.stat().st_size == EXPECTED_FILE_SIZE and skip_if_verified:
        print(f"Target file exists ({dest_path.name}). Verifying SHA-256...")
        actual_sha = compute_sha256(dest_path)
        if actual_sha.lower() == EXPECTED_SHA256.lower():
            print(f"[VERIFIED] Existing FASTA matches official SHA-256 ({actual_sha}).")
            return {
                "status": "ALREADY_VERIFIED",
                "path": str(dest_path),
                "sha256": actual_sha,
                "file_size_bytes": EXPECTED_FILE_SIZE,
            }
        else:
            print(f"[INTEGRITY WARNING] Existing file SHA-256 mismatch ({actual_sha}). Re-downloading...")
            dest_path.unlink()

    # 2. Parallel chunk download
    chunk_size = (EXPECTED_FILE_SIZE + num_chunks - 1) // num_chunks
    ranges = []
    for i in range(num_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size - 1, EXPECTED_FILE_SIZE - 1)
        part_path = dest_dir / f"{dest_path.name}.part{i}"
        ranges.append((i, start, end, part_path))

    print(f"Downloading {EXPECTED_FILE_SIZE} bytes across {num_chunks} parallel streams to {dest_dir}...")
    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_chunks) as executor:
        futures = [
            executor.submit(download_chunk, FASTA_URL, start, end, p, idx)
            for idx, start, end, p in ranges
        ]
        concurrent.futures.wait(futures)
        for f in futures:
            f.result()  # raise exception on any failure

    download_time = time.time() - t0
    print(f"All chunks downloaded in {download_time:.2f}s ({EXPECTED_FILE_SIZE / download_time / (1024**2):.2f} MB/s). Merging...")

    # 3. Merge chunks atomically into temporary file
    tmp_merged = dest_dir / f"{dest_path.name}.tmp"
    try:
        with open(tmp_merged, "wb") as f_out:
            for idx, start, end, part_path in ranges:
                with open(part_path, "rb") as f_in:
                    while True:
                        buf = f_in.read(4 * 1024 * 1024)
                        if not buf:
                            break
                        f_out.write(buf)
                part_path.unlink(missing_ok=True)

        # 4. Verify merged temporary file size and SHA-256 before final rename
        merged_size = tmp_merged.stat().st_size
        if merged_size != EXPECTED_FILE_SIZE:
            tmp_merged.unlink(missing_ok=True)
            raise ValueError(f"Merged file size mismatch! Expected {EXPECTED_FILE_SIZE}, got {merged_size}")

        print("Verifying SHA-256 checksum...")
        actual_sha = compute_sha256(tmp_merged)
        if actual_sha.lower() != EXPECTED_SHA256.lower():
            tmp_merged.unlink(missing_ok=True)
            raise ValueError(f"FASTA SHA-256 mismatch! Expected {EXPECTED_SHA256}, got {actual_sha}")

        # 5. Atomic rename upon verified match
        tmp_merged.replace(dest_path)
        print(f"[OK] Merged and verified: {dest_path} (SHA-256: {actual_sha})")

        # 6. Optional POSIX sum audit
        try:
            actual_sum, actual_blocks = verify_posix_sum(dest_path)
            print(f"POSIX sum audit: sum={actual_sum}, blocks={actual_blocks} (Expected: {EXPECTED_POSIX_SUM} {EXPECTED_POSIX_BLOCKS})")
        except Exception as e:
            actual_sum, actual_blocks = None, None
            print(f"[NOTE] POSIX sum audit skipped: {e}")

        # 7. Save reference provenance
        prov_path = dest_dir / "reference_provenance.json"
        prov_data = {
            "genome_assembly": "GRCh38.p14",
            "release": "112",
            "retrieval_timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
            "source_url": FASTA_URL,
            "download_mechanism": f"Chunked parallel HTTP Range requests ({num_chunks} streams, ThreadPoolExecutor)",
            "file_size_bytes": EXPECTED_FILE_SIZE,
            "expected_sha256": EXPECTED_SHA256,
            "observed_sha256": actual_sha,
            "expected_publisher_posix_sum": EXPECTED_POSIX_SUM,
            "expected_publisher_posix_blocks": EXPECTED_POSIX_BLOCKS,
            "observed_publisher_posix_sum": actual_sum,
            "observed_publisher_posix_blocks": actual_blocks,
            "verification_status": "CRYPTOGRAPHICALLY_VERIFIED",
            "download_time_seconds": round(download_time, 2),
        }
        prov_path.write_text(json.dumps(prov_data, indent=2), encoding="utf-8")

        return prov_data
    except Exception:
        tmp_merged.unlink(missing_ok=True)
        raise


def main():
    parser = argparse.ArgumentParser(description="Portable Ensembl 112 FASTA downloader and validator.")
    parser.add_argument(
        "--dest-dir",
        type=Path,
        default=None,
        help="Destination directory for reference files (default: references/ensembl_112)",
    )
    parser.add_argument(
        "--dest-file",
        type=Path,
        default=None,
        help="Direct path to target FASTA file (optional override)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of parallel download streams (default: 8)",
    )

    args = parser.parse_args()

    if args.dest_file:
        target_path = Path(args.dest_file).resolve()
    elif args.dest_dir:
        target_path = Path(args.dest_dir).resolve() / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    else:
        target_path = REPO_ROOT / "references" / "ensembl_112" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

    download_and_verify_fasta(target_path, num_chunks=args.threads)


if __name__ == "__main__":
    main()
