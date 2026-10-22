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
from .manifest import (
    APPROVED_PILOT_RUNS,
    APPROVED_RUN_ACCESSIONS,
    BIOLOGICAL_REPLICATE_UNIT,
    TARGET_BIOPROJECT,
    TARGET_DEVELOPMENTAL_DAY,
    TARGET_ORGANISM,
    TARGET_SRA_STUDY,
)
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
        self.config_dir = self.repo_root / "config" / "quantification_bridge"

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
        self.gtf_path = self.data_dir / "raw" / "references" / "Homo_sapiens.GRCh38.112.gtf.gz"
        self.fasta_path = self.ref_dir / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        self.bed12_path = self.ref_dir / "Homo_sapiens.GRCh38.112.bed"
        self.star_index_dir = self.ref_dir / "star_index"

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
                "source /home/chari/miniconda3/bin/activate dnt_bridge 2>/dev/null || true; "
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
        if not self.fasta_path.exists():
            print(f"Acquiring FASTA from {fasta_spec['source_url']}...")
            dl_cmd = f"python {to_wsl_path(self.repo_root / 'scripts' / 'download_reference_fasta.py')}"
            res = run_wsl(dl_cmd)
            if res.returncode != 0:
                raise RuntimeError(f"Reference FASTA acquisition failed:\n{res.stdout}\n{res.stderr}")

        observed_fasta_sha = sha256_file(self.fasta_path)
        if observed_fasta_sha != fasta_spec["sha256"]:
            raise AssertionError(f"FASTA SHA-256 mismatch! Expected {fasta_spec['sha256']}, got {observed_fasta_sha}")
        print(f"[OK] Reference FASTA verified (SHA-256: {observed_fasta_sha})")

        # 2. GTF Verification
        gtf_spec = refs_spec["gtf"]
        if not self.gtf_path.exists():
            print(f"Acquiring GTF from {gtf_spec['source_url']}...")
            self.gtf_path.parent.mkdir(parents=True, exist_ok=True)
            wsl_gtf = to_wsl_path(self.gtf_path)
            curl_cmd = f"curl -C - --retry 5 --retry-delay 5 -L '{gtf_spec['source_url']}' -o '{wsl_gtf}'"
            res = run_wsl(curl_cmd)
            if res.returncode != 0:
                raise RuntimeError(f"GTF download failed:\n{res.stderr}")

        observed_gtf_sha = sha256_file(self.gtf_path)
        if observed_gtf_sha != gtf_spec["sha256"]:
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

    def dry_run_library_acquisition(self) -> dict[str, Any]:
        """Perform a governed dry-run of ENA acquisition without downloading any files.

        Validates all 12 runs and 24 FASTQ targets against config/quantification_bridge/pilot_manifest.json.
        """
        targets = []
        for sample in APPROVED_PILOT_RUNS:
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
            missing = json_runs - target_runs
            unexpected = target_runs - json_runs
            assert not missing, f"Manifest divergence: missing runs {missing}"
            assert not unexpected, f"Manifest divergence: unexpected runs {unexpected}"

        return {
            "status": "DRY_RUN_SUCCESS",
            "total_runs": len(APPROVED_PILOT_RUNS),
            "total_fastqs": len(targets),
            "run_ids": list(APPROVED_RUN_ACCESSIONS),
            "targets": targets,
        }

    def acquire_and_verify_libraries(self) -> dict[str, Any]:
        """Download all 24 FASTQ files for the 12 approved pilot runs directly from ENA."""
        results = {}
        print(f"Acquiring 12 approved pilot libraries ({len(APPROVED_PILOT_RUNS) * 2} FASTQ files)...")

        for sample in APPROVED_PILOT_RUNS:
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

    def validate_paired_fastqs_all(self) -> dict[str, Any]:
        """Perform structural synchronization validation across all 12 downloaded library pairs."""
        validation_results = {}
        for sample in APPROVED_PILOT_RUNS:
            r1_path = self.raw_fastqs_dir / f"{sample.run_accession}_1.fastq.gz"
            r2_path = self.raw_fastqs_dir / f"{sample.run_accession}_2.fastq.gz"
            print(f"Validating FASTQ pairing: {sample.run_accession} ({sample.sample_title})...")
            v_res = validate_paired_fastqs(r1_path, r2_path, sample_records_to_check=100000)
            validation_results[sample.run_accession] = v_res

        return validation_results

    def run_fastp_all(self, threads: int = 16) -> dict[str, Any]:
        """Execute governed treatment-blind fastp preprocessing for all 12 libraries."""
        fastp_results = {}
        for sample in APPROVED_PILOT_RUNS:
            r1_in = self.raw_fastqs_dir / f"{sample.run_accession}_1.fastq.gz"
            r2_in = self.raw_fastqs_dir / f"{sample.run_accession}_2.fastq.gz"
            r1_out = self.trimmed_dir / f"{sample.run_accession}_trimmed_1.fastq.gz"
            r2_out = self.trimmed_dir / f"{sample.run_accession}_trimmed_2.fastq.gz"
            json_out = self.qc_reports_dir / f"{sample.run_accession}_fastp.json"
            html_out = self.qc_reports_dir / f"{sample.run_accession}_fastp.html"

            print(f"Executing fastp for {sample.run_accession} ({sample.sample_title})...")
            res = run_fastp_preprocessing(
                raw_r1=r1_in,
                raw_r2=r2_in,
                out_r1=r1_out,
                out_r2=r2_out,
                json_report=json_out,
                html_report=html_out,
                sample_id=sample.run_accession,
                threads=threads,
            )
            fastp_results[sample.run_accession] = res

        return fastp_results

    def build_star_index(self, threads: int = 16) -> dict[str, Any]:
        """Generate governed STAR genome index."""
        print(f"Building STAR genome index at {self.star_index_dir} (threads={threads})...")
        res = build_star_index(
            fasta_path=self.fasta_path,
            gtf_path=self.gtf_path,
            index_dir=self.star_index_dir,
            sjdb_overhang=149,
            threads=threads,
        )
        return res

    def align_12_libraries_star(self, threads: int = 16) -> dict[str, Any]:
        """Align all 12 trimmed libraries using STAR and generate BAM indexes."""
        alignment_results = {}
        for sample in APPROVED_PILOT_RUNS:
            r1_trim = self.trimmed_dir / f"{sample.run_accession}_trimmed_1.fastq.gz"
            r2_trim = self.trimmed_dir / f"{sample.run_accession}_trimmed_2.fastq.gz"
            out_prefix = self.alignments_dir / f"{sample.run_accession}_"

            print(f"Aligning {sample.run_accession} with STAR...")
            res = align_reads_star(
                trimmed_r1=r1_trim,
                trimmed_r2=r2_trim,
                index_dir=self.star_index_dir,
                output_prefix=out_prefix,
                threads=threads,
            )

            # Generate BAM index using samtools index
            bam_path = Path(res["bam_path"])
            print(f"Indexing BAM: {bam_path.name}...")
            idx_res = run_wsl(f"samtools index {to_wsl_path(bam_path)}")
            if idx_res.returncode != 0:
                raise RuntimeError(f"samtools index failed for {bam_path}:\n{idx_res.stderr}")

            bai_path = Path(f"{bam_path}.bai")
            res["bai_path"] = str(bai_path)
            alignment_results[sample.run_accession] = res

        return alignment_results

    def run_rseqc_strandedness_diagnostics(self) -> dict[str, Any]:
        """Execute RSeQC infer_experiment.py across all 12 aligned BAMs."""
        strandedness_results = {}
        for sample in APPROVED_PILOT_RUNS:
            bam_path = self.alignments_dir / f"{sample.run_accession}_Aligned.sortedByCoord.out.bam"
            out_txt = self.qc_reports_dir / f"{sample.run_accession}_infer_experiment.txt"

            print(f"Running RSeQC infer_experiment for {sample.run_accession}...")
            res = run_infer_experiment(
                bam_path=bam_path,
                bed_path=self.bed12_path,
                sample_size=200000,
                mapq_min=10,
                output_txt_path=out_txt,
            )
            strandedness_results[sample.run_accession] = res

        return strandedness_results

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
        out_file = self.bridge_out_dir / "prjna802072_strandedness_evidence.json"

        # Determine dominant strand across all runs
        forward_fracs = [m.get("forward_fraction", 0.0) for m in strandedness_metrics.values() if m.get("forward_fraction") is not None]
        reverse_fracs = [m.get("reverse_fraction", 0.0) for m in strandedness_metrics.values() if m.get("reverse_fraction") is not None]
        unexplained_fracs = [m.get("unexplained_fraction", 0.0) for m in strandedness_metrics.values() if m.get("unexplained_fraction") is not None]

        mean_fwd = sum(forward_fracs) / len(forward_fracs) if forward_fracs else 0.0
        mean_rev = sum(reverse_fracs) / len(reverse_fracs) if reverse_fracs else 0.0
        mean_unexp = sum(unexplained_fracs) / len(unexplained_fracs) if unexplained_fracs else 0.0

        if mean_rev >= 0.70 and mean_fwd <= 0.15:
            consensus = "REVERSE_STRANDED (featureCounts -s 2)"
        elif mean_fwd >= 0.70 and mean_rev <= 0.15:
            consensus = "FORWARD_STRANDED (featureCounts -s 1)"
        elif abs(mean_fwd - mean_rev) <= 0.15:
            consensus = "UNSTRANDED (featureCounts -s 0)"
        else:
            consensus = "AMBIGUOUS_EVIDENCE_REQUIRES_REVIEW"

        evidence_payload = {
            "execution_state": "AWAITING_OWNER_STRANDEDNESS_DECISION",
            "study_accession": TARGET_BIOPROJECT,
            "sra_study": TARGET_SRA_STUDY,
            "developmental_day": TARGET_DEVELOPMENTAL_DAY,
            "organism": TARGET_ORGANISM,
            "biological_replicate_unit": BIOLOGICAL_REPLICATE_UNIT,
            "total_samples": len(APPROVED_PILOT_RUNS),
            "generated_timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
            "strandedness_summary": {
                "mean_forward_fraction": round(mean_fwd, 4),
                "mean_reverse_fraction": round(mean_rev, 4),
                "mean_unexplained_fraction": round(mean_unexp, 4),
                "rseqc_consensus_indication": consensus,
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

    def run_through_strandedness_gate(self, threads: int = 16) -> str:
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
        fastp_metrics = self.run_fastp_all(threads=threads)

        print("=== STAGE 7: STAR GENOME INDEXING ===")
        star_index_metrics = self.build_star_index(threads=threads)

        print("=== STAGE 8: STAR ALIGNMENT & BAM INDEXING ===")
        alignment_metrics = self.align_12_libraries_star(threads=threads)

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
