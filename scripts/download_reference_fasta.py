"""Download and verify Ensembl Release 112 GRCh38 primary assembly FASTA with chunked parallel retrieval and POSIX sum verification."""

import concurrent.futures
import hashlib
import json
import subprocess
import time
import urllib.request
from pathlib import Path

FASTA_URL = "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
EXPECTED_POSIX_SUM = 22450
EXPECTED_POSIX_BLOCKS = 861294
EXPECTED_FILE_SIZE = 881964081

dest_dir = Path(r"C:\Users\chari\.gemini\antigravity\scratch\DNT-Organoid-AI\references\ensembl_112")
dest_dir.mkdir(parents=True, exist_ok=True)
dest_path = dest_dir / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"


def download_chunk(chunk_idx: int, start_byte: int, end_byte: int, part_path: Path, max_retries: int = 5) -> None:
    expected_len = end_byte - start_byte + 1
    if part_path.exists() and part_path.stat().st_size == expected_len:
        print(f"[Chunk {chunk_idx}] already completed ({expected_len} bytes).")
        return

    for attempt in range(1, max_retries + 1):
        try:
            req = urllib.request.Request(
                FASTA_URL,
                headers={
                    "Range": f"bytes={start_byte}-{end_byte}",
                    "User-Agent": "DNT-Organoid-AI/2.0",
                },
            )
            with urllib.request.urlopen(req, timeout=30) as resp, open(part_path, "wb") as f_out:
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


def main():
    print(f"Target FASTA: {dest_path}")
    print(f"Expected Size: {EXPECTED_FILE_SIZE} bytes")
    print(f"Expected POSIX Sum: {EXPECTED_POSIX_SUM} {EXPECTED_POSIX_BLOCKS}")

    # Check if target file already exists and is fully valid
    if dest_path.exists() and dest_path.stat().st_size == EXPECTED_FILE_SIZE:
        print("Existing file matches expected size. Testing POSIX checksum...")
        wsl_dest = f"/mnt/c/{dest_path.as_posix().replace('C:/', '')}"
        res_sum = subprocess.run(["wsl.exe", "-e", "sum", wsl_dest], capture_output=True, text=True)
        parts = res_sum.stdout.strip().split()
        if len(parts) >= 2 and int(parts[0]) == EXPECTED_POSIX_SUM and int(parts[1]) == EXPECTED_POSIX_BLOCKS:
            print("Existing file already verified with matching POSIX checksum!")
            finish_verification(dest_path, EXPECTED_FILE_SIZE, EXPECTED_POSIX_SUM, EXPECTED_POSIX_BLOCKS, 0.0)
            return

    # Remove partial single-file if existing
    if dest_path.exists():
        dest_path.unlink()

    num_chunks = 8
    chunk_size = (EXPECTED_FILE_SIZE + num_chunks - 1) // num_chunks
    ranges = []
    for i in range(num_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size - 1, EXPECTED_FILE_SIZE - 1)
        part_path = dest_dir / f"Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.part{i}"
        ranges.append((i, start, end, part_path))

    print(f"Downloading {EXPECTED_FILE_SIZE} bytes across {num_chunks} parallel connections...")
    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_chunks) as executor:
        futures = [
            executor.submit(download_chunk, idx, start, end, p)
            for idx, start, end, p in ranges
        ]
        concurrent.futures.wait(futures)
        for f in futures:
            f.result()  # raise exception if any failed

    dt = time.time() - t0
    print(f"All {num_chunks} chunks downloaded in {dt:.2f}s ({EXPECTED_FILE_SIZE / dt / 1024 / 1024:.2f} MB/s). Merging...")

    # Merge chunks
    tmp_merged = dest_dir / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.tmp"
    with open(tmp_merged, "wb") as f_out:
        for idx, start, end, part_path in ranges:
            with open(part_path, "rb") as f_in:
                while True:
                    b = f_in.read(4 * 1024 * 1024)
                    if not b:
                        break
                    f_out.write(b)
            # Remove chunk part
            part_path.unlink(missing_ok=True)

    tmp_merged.replace(dest_path)
    print(f"Merged into {dest_path} (size: {dest_path.stat().st_size} bytes).")

    finish_verification(dest_path, dest_path.stat().st_size, EXPECTED_POSIX_SUM, EXPECTED_POSIX_BLOCKS, dt)


def finish_verification(dest_path: Path, file_size: int, exp_sum: int, exp_blocks: int, download_time: float):
    # 1. Verify POSIX sum via WSL
    print("Verifying POSIX sum checksum against Ensembl CHECKSUMS...")
    wsl_dest = f"/mnt/c/{dest_path.as_posix().replace('C:/', '')}"
    res_sum = subprocess.run(
        ["wsl.exe", "-e", "sum", wsl_dest],
        capture_output=True,
        text=True,
        check=True,
    )
    parts = res_sum.stdout.strip().split()
    actual_sum = int(parts[0])
    actual_blocks = int(parts[1])
    print(f"POSIX sum output: {res_sum.stdout.strip()}")
    print(f"Expected: sum={exp_sum}, blocks={exp_blocks}")

    if actual_sum != exp_sum or actual_blocks != exp_blocks:
        print(f"ERROR: POSIX checksum mismatch! Got sum={actual_sum}, blocks={actual_blocks}")
        dest_path.unlink(missing_ok=True)
        raise AssertionError(
            f"POSIX checksum mismatch: expected {exp_sum} {exp_blocks}, got {actual_sum} {actual_blocks}"
        )

    print("PUBLISHER CHECKSUM VERIFICATION SUCCESSFUL!")

    # 2. Calculate local SHA-256
    print("Calculating local SHA-256...")
    sha = hashlib.sha256()
    with open(dest_path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            sha.update(chunk)
    local_sha = sha.hexdigest()
    print(f"Local SHA-256: {local_sha}")

    # 3. Write provenance record
    prov_path = dest_dir / "reference_provenance.json"
    prov_data = {
        "genome_assembly": "GRCh38.p14",
        "release": "112",
        "retrieval_timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "source_url": FASTA_URL,
        "download_mechanism": "Chunked parallel HTTP Range requests (8 streams, ThreadPoolExecutor)",
        "file_size_bytes": file_size,
        "expected_publisher_posix_sum": exp_sum,
        "expected_publisher_posix_blocks": exp_blocks,
        "observed_publisher_posix_sum": actual_sum,
        "observed_publisher_posix_blocks": actual_blocks,
        "publisher_verification_status": "VERIFIED_MATCH",
        "local_sha256": local_sha,
        "download_time_seconds": round(download_time, 2),
    }
    prov_path.write_text(json.dumps(prov_data, indent=2), encoding="utf-8")
    print(f"Provenance saved to {prov_path}")
    print("ALL VERIFICATIONS PASSED SUCCESSFULLY.")


if __name__ == "__main__":
    main()
