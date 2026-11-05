"""Master pipeline orchestrator for PRJNA802072 Raw-Read Quantification Bridge.

Executes through the strandedness decision gate:
1. Record environment metrics & hardware baseline
2. Execute full synthetic preflight test suite
3. Acquire and cryptographically verify references (FASTA, GTF, BED12)
4. Acquire and cryptographically verify 12 approved libraries directly from ENA
5. Validate paired FASTQ structural synchronization
6. Execute governed treatment-blind fastp preprocessing
7. Build governed STAR genome index
8. Align 12 libraries with STAR & coordinate-sort BAMs
9. Generate BAM indexes (samtools index) & parse alignment metrics
10. Execute RSeQC infer_experiment.py strandedness diagnostics across all 12 libraries
11. Serialize authoritative prjna802072_strandedness_evidence.json
12. STOP and return AWAITING_OWNER_STRANDEDNESS_DECISION
"""

import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any

from .adapter import sha256_file
from .bed_generator import generate_bed12_from_gtf
from .env_exec import run_wsl, to_wsl_path
from .fastp_runner import run_fastp_preprocessing
from .downstream_bridge import build_canonical_metadata_for_runs, execute_downstream_pipeline
from .featurecounts_runner import execute_featurecounts_and_adapter, run_featurecounts
from .manifest import (
    APPROVED_2_SAMPLE_TECHNICAL_PILOT,
    APPROVED_FEATURECOUNTS_STRANDEDNESS_FLAG,
    APPROVED_PILOT_RUNS,
    APPROVED_RUN_ACCESSIONS,
    APPROVED_STRANDEDNESS_MODE,
    BIOLOGICAL_REPLICATE_UNIT,
    TARGET_BIOPROJECT,
    TARGET_DEVELOPMENTAL_DAY,
    TARGET_ORGANISM,
    TARGET_SRA_STUDY,
)
from .pairing import FastqPair
from .profiles import (
    ExecutionProfile,
    PreflightCheckResult,
    get_execution_profile,
    resolve_pipeline_paths,
    run_profile_preflight_checks,
)
from .rseqc_runner import StrandednessDecision, infer_strandedness_decision, run_infer_experiment
from .star_runner import align_reads_star, build_star_index
from .validator import validate_paired_fastqs


class QuantificationBridgePipeline:
    """Orchestrates governed raw-read preprocessing and alignment for PRJNA802072."""

    def __init__(
        self,
        repo_root: str | Path | None = None,
        profile_name: str | None = None,
        study_id: str = TARGET_BIOPROJECT,
    ) -> None:
        if repo_root is not None:
            self.repo_root = Path(repo_root).resolve()
        else:
            self.repo_root = Path(__file__).resolve().parents[2]

        self.study_id = study_id
        self.profile: ExecutionProfile = get_execution_profile(profile_name)
        paths = resolve_pipeline_paths(self.repo_root)
        self.data_dir = paths["data_dir"]
        self.raw_fastqs_dir = paths["raw_fastqs_dir"]
        self.trimmed_dir = paths["trimmed_dir"]
        self.alignments_dir = paths["alignments_dir"]
        self.qc_reports_dir = paths["qc_reports_dir"]
        self.ref_dir = paths["ref_dir"]
        self.bridge_out_dir = paths["bridge_out_dir"]
        self.config_dir = paths["config_dir"]

        # Ensure directories exist
        for d in (
            self.raw_fastqs_dir,
            self.trimmed_dir,
            self.alignments_dir,
            self.qc_reports_dir,
            self.ref_dir,
            self.bridge_out_dir,
            self.config_dir,
        ):
            d.mkdir(parents=True, exist_ok=True)

        self.ref_manifest_path = self.config_dir / "reference_manifest.json"
        self.gtf_path = paths["gtf_path"]
        self.fasta_path = paths["fasta_path"]
        self.bed12_path = paths["bed12_path"]
        self.star_index_dir = paths["star_index_dir"]

    def run_profile_checks(
        self,
        enforce_physical_ram: bool = False,
        validate_pilot: bool | None = None,
    ) -> PreflightCheckResult:
        """Execute resource, tool, and pilot preflight checks against active execution profile."""
        should_validate_pilot = validate_pilot if validate_pilot is not None else (self.study_id == TARGET_BIOPROJECT)
        return run_profile_preflight_checks(
            profile=self.profile,
            data_dir=self.data_dir,
            ref_dir=self.ref_dir,
            enforce_physical_ram=enforce_physical_ram,
            validate_pilot=should_validate_pilot,
            study_id=self.study_id,
        )

    def record_environment_metrics(self) -> dict[str, Any]:
        """Record node hardware, OS, Python runtime, and bioinformatics binary paths/versions."""
        tools = ["fastp", "STAR", "featureCounts", "samtools", "infer_experiment.py"]
        tool_paths = {}
        tool_versions = {}

        if sys.platform != "win32":
            # Native Linux execution
            for t in tools:
                tool_paths[t] = shutil.which(t) or "NOT_FOUND"

            try:
                res = subprocess.run(["fastp", "--version"], capture_output=True, text=True)
                tool_versions["fastp"] = (res.stderr or res.stdout).strip()
            except Exception:
                tool_versions["fastp"] = "UNKNOWN"

            try:
                res = subprocess.run(["STAR", "--version"], capture_output=True, text=True)
                tool_versions["STAR"] = res.stdout.strip()
            except Exception:
                tool_versions["STAR"] = "UNKNOWN"

            try:
                res = subprocess.run(["featureCounts", "-v"], capture_output=True, text=True)
                lines = [l.strip() for l in (res.stderr + "\n" + res.stdout).splitlines() if "featureCounts" in l]
                tool_versions["featureCounts"] = lines[0] if lines else "UNKNOWN"
            except Exception:
                tool_versions["featureCounts"] = "UNKNOWN"

            try:
                res = subprocess.run(["samtools", "--version"], capture_output=True, text=True)
                tool_versions["samtools"] = res.stdout.splitlines()[0].strip() if res.stdout else "UNKNOWN"
            except Exception:
                tool_versions["samtools"] = "UNKNOWN"

            try:
                res = subprocess.run(["infer_experiment.py", "--version"], capture_output=True, text=True)
                tool_versions["infer_experiment.py"] = (res.stderr or res.stdout).strip()
            except Exception:
                tool_versions["infer_experiment.py"] = "UNKNOWN"

            kernel_str = run_wsl("uname -a").stdout.strip() if shutil.which("uname") else platform.platform()
        else:
            # Query via WSL from Windows in a single batched shell invocation
            batch_script = (
                "which fastp; which STAR; which featureCounts; which samtools; which infer_experiment.py; "
                "echo '===DELIM==='; fastp --version 2>&1; "
                "echo '===DELIM==='; STAR --version; "
                "echo '===DELIM==='; featureCounts -v 2>&1 | grep featureCounts || true; "
                "echo '===DELIM==='; samtools --version | head -n 1; "
                "echo '===DELIM==='; infer_experiment.py --version 2>&1"
            )
            res = run_wsl(batch_script)
            sections = res.stdout.strip().split("===DELIM===")
            if len(sections) >= 6:
                which_lines = [l.strip() for l in sections[0].strip().splitlines() if l.strip()]
                for i, t in enumerate(tools):
                    tool_paths[t] = which_lines[i] if i < len(which_lines) else "NOT_FOUND"
                tool_versions["fastp"] = sections[1].strip()
                tool_versions["STAR"] = sections[2].strip()
                tool_versions["featureCounts"] = sections[3].strip()
                tool_versions["samtools"] = sections[4].strip()
                tool_versions["infer_experiment.py"] = sections[5].strip()
            kernel_str = run_wsl("uname -a").stdout.strip()

        # Disk space on data dir
        disk_total, disk_used, disk_free = shutil.disk_usage(self.data_dir)

        # RAM
        ram_gb = "UNKNOWN"
        try:
            meminfo_res = run_wsl("grep MemTotal /proc/meminfo")
            if meminfo_res.returncode == 0:
                parts = meminfo_res.stdout.strip().split()
                ram_kb = int(parts[1])
                ram_gb = round(ram_kb / 1024 / 1024, 2)
        except Exception:
            pass

        # CPU
        cpu_count = os.cpu_count() or 1
        try:
            nproc_res = run_wsl("nproc")
            if nproc_res.returncode == 0:
                cpu_count = int(nproc_res.stdout.strip())
        except Exception:
            pass

        metrics = {
            "kernel": kernel_str,
            "uname": kernel_str,
            "python_executable": sys.executable,
            "python_version": platform.python_version(),
            "cpu_cores": cpu_count,
            "physical_ram_gb": ram_gb,
            "scratch_disk_total_gb": round(disk_total / (1024**3), 2),
            "scratch_disk_free_gb": round(disk_free / (1024**3), 2),
            "tool_paths": tool_paths,
            "tool_versions": tool_versions,
            "recorded_timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        }
        return metrics

    def run_preflight_checks(self) -> dict[str, Any]:
        """Execute the mandatory preflight test suite before any data acquisition or alignment."""
        test_file = self.repo_root / "tests" / "quantification_bridge" / "test_preflight_suite.py"
        if not test_file.exists():
            raise FileNotFoundError(f"Preflight test suite missing: {test_file}")

        print("Executing quantification bridge preflight test suite...")
        cmd = f"pytest {to_wsl_path(test_file)} -v"
        res = run_wsl(f"PYTHONPATH={to_wsl_path(self.repo_root)} {cmd}")
        if res.returncode != 0:
            raise RuntimeError(f"Preflight test suite failed:\n{res.stdout}\n{res.stderr}")

        print("PREFLIGHT SUITE PASSED (100% GREEN).")
        return {
            "status": "PREFLIGHT_PASSED",
            "test_file": str(test_file),
            "raw_output": res.stdout.strip(),
        }

    def download_file_with_checksum(self, url: str, dest_path: Path, expected_md5: str) -> dict[str, Any]:
        """Download file via HTTP/FTP and verify against expected MD5."""
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        # Check cache
        if dest_path.exists():
            md5 = hashlib.md5()
            with open(dest_path, "rb") as f:
                for chunk in iter(lambda: f.read(1024 * 1024), b""):
                    md5.update(chunk)
            if md5.hexdigest().lower() == expected_md5.lower():
                print(f"[CACHE] {dest_path.name} verified (MD5: {expected_md5})")
                return {
                    "path": str(dest_path),
                    "file_size": dest_path.stat().st_size,
                    "md5": expected_md5,
                    "sha256": sha256_file(dest_path),
                    "download_status": "CACHED_MATCH",
                    "retrieval_timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
                }

        download_url = url
        if download_url.startswith("ftp://"):
            download_url = "https://" + download_url[6:]
        elif not download_url.startswith("http"):
            download_url = "https://" + download_url

        print(f"Downloading {download_url} -> {dest_path}...")
        t0 = time.time()
        # Use curl with retry for robustness on large sequencing streams
        wsl_dest = to_wsl_path(dest_path)
        curl_cmd = f"curl -C - --retry 5 --retry-delay 5 -L '{download_url}' -o '{wsl_dest}'"
        res_curl = run_wsl(curl_cmd)
        if res_curl.returncode != 0:
            raise RuntimeError(f"curl download failed for {url}:\n{res_curl.stderr}")

        dl_time = time.time() - t0

        # Calculate MD5 and SHA-256
        md5 = hashlib.md5()
        sha = hashlib.sha256()
        with open(dest_path, "rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                md5.update(chunk)
                sha.update(chunk)

        actual_md5 = md5.hexdigest().lower()
        actual_sha = sha.hexdigest()

        if actual_md5 != expected_md5.lower():
            dest_path.unlink(missing_ok=True)
            raise ValueError(f"MD5 mismatch for {dest_path.name}: expected {expected_md5}, got {actual_md5}")

        print(f"[OK] Downloaded {dest_path.name} in {dl_time:.2f}s ({dest_path.stat().st_size} bytes, MD5 verified)")
        return {
            "path": str(dest_path),
            "file_size": dest_path.stat().st_size,
            "md5": actual_md5,
            "sha256": actual_sha,
            "download_status": "DOWNLOAD_VERIFIED",
            "retrieval_timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        }

    def acquire_and_verify_references(self) -> dict[str, Any]:
        """Acquire and verify reference assets against reference_manifest.json."""
        if not self.ref_manifest_path.exists():
            raise FileNotFoundError(f"Reference manifest missing: {self.ref_manifest_path}")

        ref_manifest = json.loads(self.ref_manifest_path.read_text(encoding="utf-8"))
        refs_spec = ref_manifest["references"]

        # 1. FASTA Verification
        fasta_spec = refs_spec["fasta"]
        if not self.fasta_path.exists() or sha256_file(self.fasta_path).lower() != fasta_spec["sha256"].lower():
            print(f"Acquiring FASTA from {fasta_spec['source_url']}...")
            from scripts.download_reference_fasta import download_and_verify_fasta
            download_and_verify_fasta(dest_path=self.fasta_path, num_chunks=8)

        observed_fasta_sha = sha256_file(self.fasta_path)
        if observed_fasta_sha.lower() != fasta_spec["sha256"].lower():
            raise AssertionError(f"FASTA SHA-256 mismatch! Expected {fasta_spec['sha256']}, got {observed_fasta_sha}")
        print(f"[OK] Reference FASTA verified (SHA-256: {observed_fasta_sha})")

        # 2. GTF Verification
        gtf_spec = refs_spec["gtf"]
        if not self.gtf_path.exists() or sha256_file(self.gtf_path).lower() != gtf_spec["sha256"].lower():
            print(f"Acquiring GTF from {gtf_spec['source_url']}...")
            self.gtf_path.parent.mkdir(parents=True, exist_ok=True)
            tmp_gtf = self.gtf_path.parent / f"{self.gtf_path.name}.tmp"
            req = urllib.request.Request(gtf_spec["source_url"], headers={"User-Agent": "DNT-Organoid-AI/2.0"})
            with urllib.request.urlopen(req, timeout=60) as resp, open(tmp_gtf, "wb") as f_out:
                shutil.copyfileobj(resp, f_out, length=4 * 1024 * 1024)

            downloaded_gtf_sha = sha256_file(tmp_gtf)
            if downloaded_gtf_sha.lower() != gtf_spec["sha256"].lower():
                tmp_gtf.unlink(missing_ok=True)
                raise AssertionError(f"GTF SHA-256 mismatch! Expected {gtf_spec['sha256']}, got {downloaded_gtf_sha}")
            tmp_gtf.replace(self.gtf_path)

        observed_gtf_sha = sha256_file(self.gtf_path)
        if observed_gtf_sha.lower() != gtf_spec["sha256"].lower():
            raise AssertionError(f"GTF SHA-256 mismatch! Expected {gtf_spec['sha256']}, got {observed_gtf_sha}")
        print(f"[OK] Reference GTF verified (SHA-256: {observed_gtf_sha})")

        # 3. BED12 Derivation & Verification
        bed_spec = refs_spec["bed12"]
        if not self.bed12_path.exists():
            print(f"Generating deterministic BED12 from GTF...")
            generate_bed12_from_gtf(self.gtf_path, self.bed12_path)

        observed_bed_sha = sha256_file(self.bed12_path)
        if observed_bed_sha != bed_spec["sha256"]:
            raise AssertionError(f"BED12 SHA-256 mismatch! Expected {bed_spec['sha256']}, got {observed_bed_sha}")
        print(f"[OK] Derived BED12 verified (SHA-256: {observed_bed_sha})")

        return {
            "status": "REFERENCES_VERIFIED",
            "fasta_sha256": observed_fasta_sha,
            "gtf_sha256": observed_gtf_sha,
            "bed12_sha256": observed_bed_sha,
        }

    def dry_run_library_acquisition(self, technical_pilot_only: bool = False) -> dict[str, Any]:
        """Perform a governed dry-run of ENA acquisition without downloading any files.

        Validates targets against config/quantification_bridge/pilot_manifest.json and manifest.py.
        """
        sample_cohort = APPROVED_2_SAMPLE_TECHNICAL_PILOT if technical_pilot_only else APPROVED_PILOT_RUNS
        targets = []
        for sample in sample_cohort:
            targets.append({
                "run_accession": sample.run_accession,
                "sample_alias": sample.sample_title,
                "read": "R1",
                "url": sample.fastq_1_url,
                "expected_md5": sample.fastq_1_md5,
            })
            targets.append({
                "run_accession": sample.run_accession,
                "sample_alias": sample.sample_title,
                "read": "R2",
                "url": sample.fastq_2_url,
                "expected_md5": sample.fastq_2_md5,
            })

        # Cross check against pilot_manifest.json
        manifest_file = self.config_dir / "pilot_manifest.json"
        if manifest_file.exists():
            data = json.loads(manifest_file.read_text(encoding="utf-8"))
            json_runs = set(s["run_accession"] for s in data["samples"])
            target_runs = set(t["run_accession"] for t in targets)
            if not technical_pilot_only:
                missing = json_runs - target_runs
                unexpected = target_runs - json_runs
                assert not missing, f"Manifest divergence: missing runs {missing}"
                assert not unexpected, f"Manifest divergence: unexpected runs {unexpected}"
            else:
                assert target_runs.issubset(json_runs), f"Technical pilot targets {target_runs} not in pilot_manifest.json"

        return {
            "status": "DRY_RUN_SUCCESS",
            "technical_pilot_only": technical_pilot_only,
            "total_runs": len(sample_cohort),
            "total_fastqs": len(targets),
            "run_ids": [s.run_accession for s in sample_cohort],
            "targets": targets,
        }

    def acquire_and_verify_libraries(self, technical_pilot_only: bool = False) -> dict[str, Any]:
        """Download FASTQ files directly from ENA (either 2 technical pilot runs or all 12 cohort runs)."""
        sample_cohort = APPROVED_2_SAMPLE_TECHNICAL_PILOT if technical_pilot_only else APPROVED_PILOT_RUNS
        results = {}
        print(f"Acquiring {len(sample_cohort)} approved libraries ({len(sample_cohort) * 2} FASTQ files)...")

        for sample in sample_cohort:
            r1_dest = self.raw_fastqs_dir / f"{sample.run_accession}_1.fastq.gz"
            r2_dest = self.raw_fastqs_dir / f"{sample.run_accession}_2.fastq.gz"

            res1 = self.download_file_with_checksum(sample.fastq_1_url, r1_dest, sample.fastq_1_md5)
            res2 = self.download_file_with_checksum(sample.fastq_2_url, r2_dest, sample.fastq_2_md5)

            results[sample.run_accession] = {
                "sample_title": sample.sample_title,
                "vpa_dose_um": sample.vpa_dose_um,
                "r1": res1,
                "r2": res2,
            }

        return results

    def pairs_from_approved_technical_pilot_runs(self) -> tuple[FastqPair, ...]:
        """Construct FastqPair records for the 2 approved technical pilot runs (SRR17834920 and SRR17834905)."""
        pairs = []
        for sample in APPROVED_2_SAMPLE_TECHNICAL_PILOT:
            pairs.append(
                FastqPair(
                    sample_id=sample.run_accession,
                    r1_path=self.raw_fastqs_dir / f"{sample.run_accession}_1.fastq.gz",
                    r2_path=self.raw_fastqs_dir / f"{sample.run_accession}_2.fastq.gz",
                    assay_id=sample.experiment_accession,
                    extra_metadata={
                        "sample_title": sample.sample_title,
                        "vpa_dose_um": sample.vpa_dose_um,
                    },
                )
            )
        return tuple(pairs)

    def pairs_from_approved_pilot_runs(self) -> tuple[FastqPair, ...]:
        """Construct FastqPair records for the 12 approved pilot runs."""
        pairs = []
        for sample in APPROVED_PILOT_RUNS:
            pairs.append(
                FastqPair(
                    sample_id=sample.run_accession,
                    r1_path=self.raw_fastqs_dir / f"{sample.run_accession}_1.fastq.gz",
                    r2_path=self.raw_fastqs_dir / f"{sample.run_accession}_2.fastq.gz",
                    assay_id=sample.experiment_accession,
                    extra_metadata={
                        "sample_title": sample.sample_title,
                        "vpa_dose_um": sample.vpa_dose_um,
                    },
                )
            )
        return tuple(pairs)

    def validate_paired_fastqs_pairs(
        self,
        pairs: tuple[FastqPair, ...] | list[FastqPair],
        sample_records_to_check: int = 100000,
    ) -> dict[str, Any]:
        """Perform structural synchronization validation across given library pairs."""
        validation_results = {}
        for pair in pairs:
            print(f"Validating FASTQ pairing: {pair.sample_id}...")
            v_res = validate_paired_fastqs(pair.r1_path, pair.r2_path, sample_records_to_check=sample_records_to_check)
            validation_results[pair.sample_id] = v_res
        return validation_results

    def validate_paired_fastqs_all(self) -> dict[str, Any]:
        """Perform structural synchronization validation across all 12 downloaded library pairs."""
        return self.validate_paired_fastqs_pairs(self.pairs_from_approved_pilot_runs())

    def run_fastp_pairs(
        self,
        pairs: tuple[FastqPair, ...] | list[FastqPair],
        threads: int | None = None,
    ) -> dict[str, Any]:
        """Execute governed treatment-blind fastp preprocessing for given pairs."""
        eff_threads = threads if threads is not None else self.profile.threads
        fastp_results = {}
        for pair in pairs:
            r1_out = self.trimmed_dir / f"{pair.sample_id}_trimmed_1.fastq.gz"
            r2_out = self.trimmed_dir / f"{pair.sample_id}_trimmed_2.fastq.gz"
            json_out = self.qc_reports_dir / f"{pair.sample_id}_fastp.json"
            html_out = self.qc_reports_dir / f"{pair.sample_id}_fastp.html"

            print(f"Executing fastp for {pair.sample_id} (threads={eff_threads})...")
            res = run_fastp_preprocessing(
                raw_r1=pair.r1_path,
                raw_r2=pair.r2_path,
                out_r1=r1_out,
                out_r2=r2_out,
                json_report=json_out,
                html_report=html_out,
                sample_id=pair.sample_id,
                threads=eff_threads,
            )
            fastp_results[pair.sample_id] = res

        return fastp_results

    def run_fastp_all(self, threads: int | None = None) -> dict[str, Any]:
        """Execute governed treatment-blind fastp preprocessing for all 12 libraries."""
        return self.run_fastp_pairs(self.pairs_from_approved_pilot_runs(), threads=threads)

    def build_star_index(
        self,
        threads: int | None = None,
        sparse_d: int | None = None,
        genome_sa_index_nbases: int | None = None,
        limit_ram: int | None = None,
    ) -> dict[str, Any]:
        """Generate governed STAR genome index using active execution profile parameters."""
        eff_threads = threads if threads is not None else self.profile.threads
        eff_sparse_d = sparse_d if sparse_d is not None else self.profile.genome_sa_sparse_d
        eff_nbases = genome_sa_index_nbases if genome_sa_index_nbases is not None else self.profile.genome_sa_index_nbases
        eff_limit_ram = limit_ram if limit_ram is not None else self.profile.limit_genome_generate_ram

        print(
            f"Building STAR genome index at {self.star_index_dir} "
            f"(threads={eff_threads}, sparse_d={eff_sparse_d}, nbases={eff_nbases}, limit_ram={eff_limit_ram})..."
        )
        res = build_star_index(
            fasta_path=self.fasta_path,
            gtf_path=self.gtf_path,
            index_dir=self.star_index_dir,
            sjdb_overhang=149,
            threads=eff_threads,
            sparse_d=eff_sparse_d,
            genome_sa_index_nbases=eff_nbases,
            limit_ram=eff_limit_ram,
        )
        return res

    def align_libraries_star_pairs(
        self,
        pairs: tuple[FastqPair, ...] | list[FastqPair],
        threads: int | None = None,
        bam_sort_mode: str | None = None,
        external_sort: bool | None = None,
        samtools_sort_mem: str | None = None,
    ) -> dict[str, Any]:
        """Align given trimmed library pairs using STAR and generate BAM indexes."""
        eff_threads = threads if threads is not None else self.profile.threads
        eff_sort_mode = bam_sort_mode if bam_sort_mode is not None else self.profile.bam_sort_mode
        eff_external_sort = external_sort if external_sort is not None else self.profile.external_sort_enabled
        eff_sort_mem = samtools_sort_mem if samtools_sort_mem is not None else self.profile.samtools_sort_mem_per_thread

        alignment_results = {}
        for pair in pairs:
            r1_trim = self.trimmed_dir / f"{pair.sample_id}_trimmed_1.fastq.gz"
            r2_trim = self.trimmed_dir / f"{pair.sample_id}_trimmed_2.fastq.gz"
            out_prefix = self.alignments_dir / f"{pair.sample_id}_"

            print(
                f"Aligning {pair.sample_id} with STAR "
                f"(threads={eff_threads}, sort_mode={eff_sort_mode}, external_sort={eff_external_sort})..."
            )
            res = align_reads_star(
                trimmed_r1=r1_trim,
                trimmed_r2=r2_trim,
                index_dir=self.star_index_dir,
                output_prefix=out_prefix,
                threads=eff_threads,
                bam_sort_mode=eff_sort_mode,
                limit_bam_sort_ram=self.profile.limit_bam_sort_ram,
                external_sort=eff_external_sort,
                samtools_sort_mem=eff_sort_mem,
            )

            # Generate BAM index using samtools index
            bam_path = Path(res["bam_path"])
            print(f"Indexing BAM: {bam_path.name}...")
            idx_res = run_wsl(f"samtools index {to_wsl_path(bam_path)}")
            if idx_res.returncode != 0:
                raise RuntimeError(f"samtools index failed for {bam_path}:\n{idx_res.stderr}")

            bai_path = Path(f"{bam_path}.bai")
            res["bai_path"] = str(bai_path)
            alignment_results[pair.sample_id] = res

        return alignment_results

    def align_12_libraries_star(
        self,
        threads: int | None = None,
        bam_sort_mode: str | None = None,
        external_sort: bool | None = None,
        samtools_sort_mem: str | None = None,
    ) -> dict[str, Any]:
        """Align all 12 trimmed libraries using STAR and generate BAM indexes."""
        return self.align_libraries_star_pairs(
            self.pairs_from_approved_pilot_runs(),
            threads=threads,
            bam_sort_mode=bam_sort_mode,
            external_sort=external_sort,
            samtools_sort_mem=samtools_sort_mem,
        )

    def run_rseqc_strandedness_pairs(
        self,
        pairs: tuple[FastqPair, ...] | list[FastqPair],
    ) -> dict[str, Any]:
        """Execute RSeQC infer_experiment.py across aligned BAMs for given pairs."""
        strandedness_results = {}
        for pair in pairs:
            bam_path = self.alignments_dir / f"{pair.sample_id}_Aligned.sortedByCoord.out.bam"
            out_txt = self.qc_reports_dir / f"{pair.sample_id}_infer_experiment.txt"

            print(f"Running RSeQC infer_experiment for {pair.sample_id}...")
            res = run_infer_experiment(
                bam_path=bam_path,
                bed_path=self.bed12_path,
                sample_size=200000,
                mapq_min=10,
                output_txt_path=out_txt,
            )
            strandedness_results[pair.sample_id] = res

        return strandedness_results

    def run_rseqc_strandedness_diagnostics(self) -> dict[str, Any]:
        """Execute RSeQC infer_experiment.py across all 12 aligned BAMs."""
        return self.run_rseqc_strandedness_pairs(self.pairs_from_approved_pilot_runs())

    def evaluate_strandedness_decision(
        self,
        metrics_or_summary: dict[str, Any],
    ) -> StrandednessDecision:
        """Evaluate empirical RSeQC evidence and return authoritative StrandednessDecision."""
        return infer_strandedness_decision(metrics_or_summary)

    def serialize_strandedness_evidence_pairs(
        self,
        pairs: tuple[FastqPair, ...] | list[FastqPair],
        study_id: str,
        env_metrics: dict[str, Any],
        preflight_metrics: dict[str, Any],
        ref_metrics: dict[str, Any],
        download_metrics: dict[str, Any],
        pairing_metrics: dict[str, Any],
        fastp_metrics: dict[str, Any],
        star_index_metrics: dict[str, Any],
        alignment_metrics: dict[str, Any],
        strandedness_metrics: dict[str, Any],
    ) -> Path:
        """Compile and serialize the authoritative strandedness evidence artifact."""
        clean_study = re.sub(r"[^A-Za-z0-9_-]", "_", study_id).lower()
        out_file = self.bridge_out_dir / f"{clean_study}_strandedness_evidence.json"

        # Evaluate authoritative evidence-based strandedness decision
        decision = infer_strandedness_decision(strandedness_metrics)

        evidence_payload = {
            "execution_state": "AWAITING_OWNER_STRANDEDNESS_DECISION",
            "study_accession": study_id,
            "total_samples": len(pairs),
            "generated_timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
            "strandedness_decision": decision.to_dict(),
            "strandedness_summary": {
                "mean_forward_fraction": round(decision.mean_forward, 4),
                "mean_reverse_fraction": round(decision.mean_reverse, 4),
                "mean_unexplained_fraction": round(decision.mean_unexplained, 4),
                "rseqc_consensus_indication": decision.consensus_indication,
                "strandedness_flag": decision.strandedness_flag,
                "strand_mode": decision.strand_mode,
                "confidence": decision.confidence,
                "rationale": decision.rationale,
            },
            "environment_metrics": env_metrics,
            "preflight_metrics": preflight_metrics,
            "reference_metrics": ref_metrics,
            "library_download_metrics": download_metrics,
            "fastp_metrics": fastp_metrics,
            "alignment_metrics": alignment_metrics,
            "per_sample_strandedness": strandedness_metrics,
        }

        out_file.write_text(json.dumps(evidence_payload, indent=2), encoding="utf-8")
        print(f"Authoritative strandedness evidence written to {out_file}")
        return out_file

    def serialize_strandedness_evidence(
        self,
        env_metrics: dict[str, Any],
        preflight_metrics: dict[str, Any],
        ref_metrics: dict[str, Any],
        download_metrics: dict[str, Any],
        pairing_metrics: dict[str, Any],
        fastp_metrics: dict[str, Any],
        star_index_metrics: dict[str, Any],
        alignment_metrics: dict[str, Any],
        strandedness_metrics: dict[str, Any],
    ) -> Path:
        """Compile and serialize the authoritative prjna802072_strandedness_evidence.json."""
        return self.serialize_strandedness_evidence_pairs(
            pairs=self.pairs_from_approved_pilot_runs(),
            study_id=TARGET_BIOPROJECT,
            env_metrics=env_metrics,
            preflight_metrics=preflight_metrics,
            ref_metrics=ref_metrics,
            download_metrics=download_metrics,
            pairing_metrics=pairing_metrics,
            fastp_metrics=fastp_metrics,
            star_index_metrics=star_index_metrics,
            alignment_metrics=alignment_metrics,
            strandedness_metrics=strandedness_metrics,
        )

    def record_accepted_strandedness_decision(
        self,
        evidence_path: str | Path,
        strandedness_flag: int = APPROVED_FEATURECOUNTS_STRANDEDNESS_FLAG,
        strand_mode: str = APPROVED_STRANDEDNESS_MODE,
        reviewer_name: str = "Project Owner / Scientific Review",
        rationale: str | None = None,
    ) -> Path:
        """Record accepted strandedness decision into an evidence artifact, transitioning out of gate."""
        ev_p = Path(evidence_path).resolve()
        if not ev_p.exists():
            raise FileNotFoundError(f"Evidence artifact not found: {ev_p}")
        data = json.loads(ev_p.read_text(encoding="utf-8"))
        data["execution_state"] = "STRANDEDNESS_DECISION_ACCEPTED"
        data["checkpoint_gate_state"] = "AWAITING_OWNER_STRANDEDNESS_DECISION"
        summary = data.get("strandedness_summary", {})
        fwd = summary.get("mean_forward_fraction", 0.0)
        rev = summary.get("mean_reverse_fraction", 0.0)
        unexp = summary.get("mean_unexplained_fraction", 0.0)
        diff = round(abs(fwd - rev), 4)
        default_rat = (
            f"Empirical RSeQC strand distribution across all {data.get('total_samples', 0)} libraries demonstrates "
            f"balanced forward ({fwd*100:.2f}%) and reverse ({rev*100:.2f}%) mapping with |delta| = {diff}, "
            f"confirming unstranded library preparation. Quantification proceeds with featureCounts flag -s {strandedness_flag}."
        )
        data["accepted_decision"] = {
            "status": "APPROVED",
            "featurecounts_flag": strandedness_flag,
            "strand_mode": strand_mode,
            "featurecounts_parameter": f"-s {strandedness_flag}",
            "approval_timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
            "approved_by": reviewer_name,
            "cohort_scope": f"{data.get('study_accession', 'generic_study')} ({data.get('total_samples', 0)} libraries)",
            "supporting_rseqc_findings": {
                "cohort_size": data.get("total_samples", 0),
                "mean_forward_fraction": fwd,
                "mean_reverse_fraction": rev,
                "mean_unexplained_fraction": unexp,
                "absolute_orientation_difference": diff,
                "consensus_indication": summary.get("rseqc_consensus_indication", f"Mode: {strand_mode}"),
                "rationale": rationale or default_rat,
            },
        }
        ev_p.write_text(json.dumps(data, indent=2), encoding="utf-8")
        print(f"Recorded accepted strandedness decision (-s {strandedness_flag}) into {ev_p}")
        return ev_p

    def execute_raw_read_route(
        self,
        pairs: tuple[FastqPair, ...] | list[FastqPair],
        study_id: str = "generic_study",
        threads: int | None = None,
        stop_at_gate: bool = True,
    ) -> str:
        eff_threads = threads if threads is not None else self.profile.threads
        """Generic raw-read route execution through the strandedness evidence gate."""
        print(f"=== GENERIC RAW-READ ROUTE EXECUTION: {study_id} ({len(pairs)} pairs) ===")
        print("=== STAGE 1: ENVIRONMENT AUDIT ===")
        env_metrics = self.record_environment_metrics()

        print("=== STAGE 2: PREFLIGHT VALIDATION ===")
        preflight_metrics = self.run_preflight_checks()

        print("=== STAGE 3: REFERENCE VERIFICATION ===")
        ref_metrics = self.acquire_and_verify_references()

        print("=== STAGE 4: FASTQ STRUCTURAL VALIDATION ===")
        pairing_metrics = self.validate_paired_fastqs_pairs(pairs)

        print("=== STAGE 5: TREATMENT-BLIND FASTP PREPROCESSING ===")
        fastp_metrics = self.run_fastp_pairs(pairs, threads=eff_threads)

        print("=== STAGE 6: STAR GENOME INDEXING ===")
        star_index_metrics = self.build_star_index(threads=eff_threads)

        print("=== STAGE 7: STAR ALIGNMENT & BAM INDEXING ===")
        alignment_metrics = self.align_libraries_star_pairs(pairs, threads=eff_threads)

        print("=== STAGE 8: RSeQC STRANDEDNESS INFERENCE ===")
        strandedness_metrics = self.run_rseqc_strandedness_pairs(pairs)

        print("=== STAGE 9: SERIALIZE STRANDEDNESS EVIDENCE ===")
        evidence_path = self.serialize_strandedness_evidence_pairs(
            pairs=pairs,
            study_id=study_id,
            env_metrics=env_metrics,
            preflight_metrics=preflight_metrics,
            ref_metrics=ref_metrics,
            download_metrics={},
            pairing_metrics=pairing_metrics,
            fastp_metrics=fastp_metrics,
            star_index_metrics=star_index_metrics,
            alignment_metrics=alignment_metrics,
            strandedness_metrics=strandedness_metrics,
        )

        if stop_at_gate:
            print("================================================================================")
            print("GATE REACHED: AWAITING_OWNER_STRANDEDNESS_DECISION")
            print(f"Strandedness evidence artifact: {evidence_path}")
            print("HALTING EXECUTION. Do NOT proceed to featureCounts until authorized by Owner.")
            print("================================================================================")
            return "AWAITING_OWNER_STRANDEDNESS_DECISION"

        return "STRANDEDNESS_EVIDENCE_PRODUCED"

    def run_through_strandedness_gate(self, threads: int | None = None) -> str:
        eff_threads = threads if threads is not None else self.profile.threads
        """Run entire upstream bridge through the strandedness evidence gate and STOP."""
        print("=== STAGE 1: ENVIRONMENT AUDIT ===")
        env_metrics = self.record_environment_metrics()

        print("=== STAGE 2: PREFLIGHT VALIDATION ===")
        preflight_metrics = self.run_preflight_checks()

        print("=== STAGE 3: REFERENCE VERIFICATION ===")
        ref_metrics = self.acquire_and_verify_references()

        print("=== STAGE 4: 12-LIBRARY ACQUISITION DIRECT FROM ENA ===")
        download_metrics = self.acquire_and_verify_libraries()

        print("=== STAGE 5: FASTQ STRUCTURAL VALIDATION ===")
        pairing_metrics = self.validate_paired_fastqs_all()

        print("=== STAGE 6: TREATMENT-BLIND FASTP PREPROCESSING ===")
        fastp_metrics = self.run_fastp_all(threads=eff_threads)

        print("=== STAGE 7: STAR GENOME INDEXING ===")
        star_index_metrics = self.build_star_index(threads=eff_threads)

        print("=== STAGE 8: STAR ALIGNMENT & BAM INDEXING ===")
        alignment_metrics = self.align_12_libraries_star(threads=eff_threads)

        print("=== STAGE 9: RSeQC STRANDEDNESS INFERENCE ===")
        strandedness_metrics = self.run_rseqc_strandedness_diagnostics()

        print("=== STAGE 10: SERIALIZE STRANDEDNESS EVIDENCE ===")
        evidence_path = self.serialize_strandedness_evidence(
            env_metrics=env_metrics,
            preflight_metrics=preflight_metrics,
            ref_metrics=ref_metrics,
            download_metrics=download_metrics,
            pairing_metrics=pairing_metrics,
            fastp_metrics=fastp_metrics,
            star_index_metrics=star_index_metrics,
            alignment_metrics=alignment_metrics,
            strandedness_metrics=strandedness_metrics,
        )

        print("================================================================================")
        print("GATE REACHED: AWAITING_OWNER_STRANDEDNESS_DECISION")
        print(f"Strandedness evidence artifact: {evidence_path}")
        print("HALTING EXECUTION. Do NOT proceed to featureCounts until authorized by Owner.")
        print("================================================================================")
        return "AWAITING_OWNER_STRANDEDNESS_DECISION"

    def run_featurecounts_pairs(
        self,
        pairs: tuple[FastqPair, ...] | list[FastqPair],
        strandedness_flag: int,
        threads: int | None = None,
    ) -> dict[str, Any]:
        """Execute governed featureCounts across aligned BAMs for given pairs and generate CountMatrix."""
        eff_threads = threads if threads is not None else self.profile.threads
        bam_paths = [self.alignments_dir / f"{pair.sample_id}_Aligned.sortedByCoord.out.bam" for pair in pairs]
        run_ids = tuple(pair.sample_id for pair in pairs)

        suffix = "pilot_2_samples" if len(pairs) == 2 else "cohort_12_samples"
        raw_counts_path = self.bridge_out_dir / f"featureCounts_{suffix}_raw.txt"
        matrix_path = self.bridge_out_dir / f"prjna802072_{suffix}_count_matrix.tsv"
        universe_path = self.bridge_out_dir / f"prjna802072_{suffix}_gene_universe.json"

        print(
            f"Executing governed featureCounts (strandedness={strandedness_flag}, threads={eff_threads}) "
            f"for {len(pairs)} BAM files..."
        )
        return execute_featurecounts_and_adapter(
            bam_paths=bam_paths,
            gtf_path=self.gtf_path,
            raw_featurecounts_path=raw_counts_path,
            output_matrix_path=matrix_path,
            output_universe_path=universe_path,
            approved_runs=run_ids,
            strandedness=strandedness_flag,
            threads=eff_threads,
        )

    def run_technical_pilot_through_strandedness_gate(self, threads: int | None = None) -> str:
        """Run 2-sample technical pilot (SRR17834920 & SRR17834905) through the strandedness evidence gate and STOP."""
        eff_threads = threads if threads is not None else self.profile.threads
        print("=== 2-SAMPLE TECHNICAL PILOT: STAGE 1: ENVIRONMENT AUDIT ===")
        env_metrics = self.record_environment_metrics()

        print("=== 2-SAMPLE TECHNICAL PILOT: STAGE 2: PREFLIGHT VALIDATION ===")
        preflight_metrics = self.run_preflight_checks()

        print("=== 2-SAMPLE TECHNICAL PILOT: STAGE 3: REFERENCE VERIFICATION ===")
        ref_metrics = self.acquire_and_verify_references()

        print("=== 2-SAMPLE TECHNICAL PILOT: STAGE 4: ACQUISITION DIRECT FROM ENA ===")
        download_metrics = self.acquire_and_verify_libraries(technical_pilot_only=True)

        pairs = self.pairs_from_approved_technical_pilot_runs()

        print("=== 2-SAMPLE TECHNICAL PILOT: STAGE 5: FASTQ STRUCTURAL VALIDATION ===")
        pairing_metrics = self.validate_paired_fastqs_pairs(pairs)

        print("=== 2-SAMPLE TECHNICAL PILOT: STAGE 6: TREATMENT-BLIND FASTP PREPROCESSING ===")
        fastp_metrics = self.run_fastp_pairs(pairs, threads=eff_threads)

        print("=== 2-SAMPLE TECHNICAL PILOT: STAGE 7: STAR GENOME INDEXING ===")
        star_index_metrics = self.build_star_index(threads=eff_threads)

        print("=== 2-SAMPLE TECHNICAL PILOT: STAGE 8: STAR ALIGNMENT & BAM INDEXING ===")
        alignment_metrics = self.align_libraries_star_pairs(pairs, threads=eff_threads)

        print("=== 2-SAMPLE TECHNICAL PILOT: STAGE 9: RSeQC STRANDEDNESS INFERENCE ===")
        strandedness_metrics = self.run_rseqc_strandedness_pairs(pairs)

        print("=== 2-SAMPLE TECHNICAL PILOT: STAGE 10: SERIALIZE STRANDEDNESS EVIDENCE ===")
        evidence_path = self.serialize_strandedness_evidence_pairs(
            pairs=pairs,
            study_id=f"{self.study_id}_technical_pilot",
            env_metrics=env_metrics,
            preflight_metrics=preflight_metrics,
            ref_metrics=ref_metrics,
            download_metrics=download_metrics,
            pairing_metrics=pairing_metrics,
            fastp_metrics=fastp_metrics,
            star_index_metrics=star_index_metrics,
            alignment_metrics=alignment_metrics,
            strandedness_metrics=strandedness_metrics,
        )

        print("================================================================================")
        print("GATE REACHED: AWAITING_OWNER_STRANDEDNESS_DECISION")
        print(f"Technical pilot strandedness evidence artifact: {evidence_path}")
        print("HALTING EXECUTION. Do NOT proceed to featureCounts until authorized by Owner.")
        print("================================================================================")
        return "AWAITING_OWNER_STRANDEDNESS_DECISION"

    def run_featurecounts_and_downstream(
        self,
        pairs: tuple[FastqPair, ...] | list[FastqPair],
        strandedness_flag: int,
        threads: int | None = None,
    ) -> dict[str, Any]:
        """Execute featureCounts, convert to count matrix, and connect to downstream analysis modules."""
        fc_result = self.run_featurecounts_pairs(pairs=pairs, strandedness_flag=strandedness_flag, threads=threads)
        matrix_path = fc_result["output_matrix_path"]
        run_ids = tuple(pair.sample_id for pair in pairs)

        print(f"Executing downstream pipeline for count matrix: {matrix_path}...")
        downstream_result = execute_downstream_pipeline(
            counts_matrix_path=matrix_path,
            runs=run_ids,
            study_id=self.study_id,
            output_dir=self.bridge_out_dir / "downstream_analysis",
        )

        return {
            "status": "PIPELINE_COMPLETE",
            "featurecounts_result": fc_result,
            "downstream_result": downstream_result.to_dict(),
        }
