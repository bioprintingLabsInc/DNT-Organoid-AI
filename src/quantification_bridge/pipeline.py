"""Master pipeline orchestrator for PRJNA802072 Raw-Read Quantification Bridge.

Executes:
1. Manifest verification
2. Reference verification & BED12 derivation
3. FASTQ downloads & transport MD5 validation
4. FASTQ structural pairing validation
5. fastp treatment-blind preprocessing
6. STAR genome indexing
7. STAR alignment & coordinate-sorted BAM generation
8. RSeQC strandedness diagnostic & evidence JSON generation
9. Provenance serialization
"""

import json
import os
import shutil
import subprocess
import time
import urllib.request
from pathlib import Path
from typing import Any

from .adapter import convert_featurecounts_to_matrix, sha256_file
from .bed_generator import generate_bed12_from_gtf
from .fastp_runner import run_fastp_preprocessing
from .manifest import APPROVED_PILOT_RUNS, APPROVED_RUN_ACCESSIONS, BIOLOGICAL_REPLICATE_UNIT
from .rseqc_runner import run_infer_experiment
from .star_runner import align_reads_star, build_star_index
from .validator import validate_paired_fastqs


class QuantificationBridgePipeline:
    """Orchestrates governed raw-read preprocessing and alignment for PRJNA802072."""

    def __init__(self, repo_root: str | Path | None = None) -> None:
        if repo_root is not None:
            self.repo_root = Path(repo_root).resolve()
        else:
            self.repo_root = Path(__file__).resolve().parents[2]
        self.data_dir = self.repo_root / "data"
        self.raw_fastqs_dir = self.data_dir / "raw" / "prjna802072_fastqs"
        self.trimmed_dir = self.data_dir / "interim" / "prjna802072_trimmed"
        self.alignments_dir = self.data_dir / "interim" / "prjna802072_alignments"
        self.qc_reports_dir = self.data_dir / "interim" / "prjna802072_qc_reports"
        self.ref_dir = self.repo_root / "references" / "ensembl_112"
        self.bridge_out_dir = self.data_dir / "interim" / "prjna802072_bridge"

        # Ensure directories exist
        for d in (self.raw_fastqs_dir, self.trimmed_dir, self.alignments_dir, self.qc_reports_dir, self.ref_dir, self.bridge_out_dir):
            d.mkdir(parents=True, exist_ok=True)

        self.gtf_path = self.data_dir / "raw" / "references" / "Homo_sapiens.GRCh38.112.gtf.gz"
        self.fasta_path = self.ref_dir / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        self.bed12_path = self.ref_dir / "Homo_sapiens.GRCh38.112.bed"
        self.star_index_dir = self.ref_dir / "star_index"

    def download_file_with_checksum(self, url: str, dest_path: Path, expected_md5: str) -> dict[str, Any]:
        """Download a file via HTTP/FTP and verify against expected MD5."""
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        # If file exists and matches MD5, skip download
        if dest_path.exists():
            import hashlib
            md5 = hashlib.md5()
            with open(dest_path, "rb") as f:
                for chunk in iter(lambda: f.read(1024 * 1024), b""):
                    md5.update(chunk)
            if md5.hexdigest().lower() == expected_md5.lower():
                print(f"File {dest_path.name} already exists and matches expected MD5 ({expected_md5})")
                return {
                    "path": str(dest_path),
                    "file_size": dest_path.stat().st_size,
                    "md5": expected_md5,
                    "sha256": sha256_file(dest_path),
                    "download_status": "CACHED_MATCH",
                }

        # Handle FTP via HTTPS if possible or direct
        download_url = url
        if download_url.startswith("ftp://"):
            download_url = "https://" + download_url[6:]
        elif not download_url.startswith("http"):
            download_url = "https://" + download_url

        print(f"Downloading {download_url} -> {dest_path}...")
        t0 = time.time()
        req = urllib.request.Request(download_url, headers={"User-Agent": "DNT-Organoid-AI/2.0"})
        with urllib.request.urlopen(req) as resp, open(dest_path, "wb") as f_out:
            shutil.copyfileobj(resp, f_out)
        dl_time = time.time() - t0

        # Calculate MD5 and SHA-256
        import hashlib
        md5 = hashlib.md5()
        sha = hashlib.sha256()
        with open(dest_path, "rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                md5.update(chunk)
                sha.update(chunk)

        actual_md5 = md5.hexdigest()
        actual_sha = sha.hexdigest()

        if actual_md5.lower() != expected_md5.lower():
            dest_path.unlink(missing_ok=True)
            raise ValueError(f"MD5 mismatch for {dest_path.name}: expected {expected_md5}, got {actual_md5}")

        print(f"Downloaded {dest_path.name} in {dl_time:.2f}s ({dest_path.stat().st_size} bytes, MD5: {actual_md5} verified)")

        return {
            "path": str(dest_path),
            "file_size": dest_path.stat().st_size,
            "md5": actual_md5,
            "sha256": actual_sha,
            "download_status": "DOWNLOAD_VERIFIED",
        }
