"""Download and verify Ensembl Release 112 GRCh38 primary assembly FASTA with resumable curl and POSIX sum verification."""

import hashlib
import json
import subprocess
import time
from pathlib import Path

FASTA_URL = "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
EXPECTED_POSIX_SUM = 22450
EXPECTED_POSIX_BLOCKS = 861294

dest_dir = Path(r"C:\Users\chari\.gemini\antigravity\scratch\DNT-Organoid-AI\references\ensembl_112")
dest_dir.mkdir(parents=True, exist_ok=True)
dest_path = dest_dir / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

print(f"Target FASTA: {dest_path}")
print(f"Downloading from {FASTA_URL} using resumable curl...")

wsl_dest = f"/mnt/c/{dest_path.as_posix().replace('C:/', '')}"
curl_cmd = (
    f"curl -C - --retry 5 --retry-delay 5 -L '{FASTA_URL}' -o '{wsl_dest}'"
)

t0 = time.time()
res_curl = subprocess.run(
    ["wsl.exe", "-e", "bash", "-c", curl_cmd],
    capture_output=True,
    text=True,
    check=False,
)

if res_curl.returncode != 0:
    print(f"curl failed with code {res_curl.returncode}:\n{res_curl.stderr}")
    raise RuntimeError(f"curl download failed: {res_curl.stderr}")

download_time = time.time() - t0
file_size = dest_path.stat().st_size
print(f"Download/verification complete in {download_time:.2f}s ({file_size} bytes, {file_size / max(download_time, 0.001) / 1024 / 1024:.2f} MB/s)")

# 1. Verify POSIX sum via WSL
print("Verifying POSIX sum checksum against Ensembl CHECKSUMS...")
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
print(f"Expected: sum={EXPECTED_POSIX_SUM}, blocks={EXPECTED_POSIX_BLOCKS}")

if actual_sum != EXPECTED_POSIX_SUM or actual_blocks != EXPECTED_POSIX_BLOCKS:
    print(f"ERROR: POSIX checksum mismatch! Got sum={actual_sum}, blocks={actual_blocks}")
    # Remove truncated/corrupt file
    if dest_path.exists():
        dest_path.unlink()
    raise AssertionError(
        f"POSIX checksum mismatch: expected {EXPECTED_POSIX_SUM} {EXPECTED_POSIX_BLOCKS}, got {actual_sum} {actual_blocks}"
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
    "download_mechanism": "curl -C - --retry 5 --retry-delay 5 -L (via WSL2)",
    "file_size_bytes": file_size,
    "expected_publisher_posix_sum": EXPECTED_POSIX_SUM,
    "expected_publisher_posix_blocks": EXPECTED_POSIX_BLOCKS,
    "observed_publisher_posix_sum": actual_sum,
    "observed_publisher_posix_blocks": actual_blocks,
    "publisher_verification_status": "VERIFIED_MATCH",
    "local_sha256": local_sha,
}
prov_path.write_text(json.dumps(prov_data, indent=2), encoding="utf-8")
print(f"Provenance saved to {prov_path}")
print("ALL VERIFICATIONS PASSED SUCCESSFULLY.")
