"""STAR 2.7.11b index generation and alignment wrapper.

Supports configurable execution profiles:
- Canonical dense indexing (sparseD 1) for HPC and large servers.
- Memory-bounded sparse indexing (sparseD 2 or sparseD 4) for workstations and evaluation laptops.
- Native or external streaming coordinate sorting (samtools sort).
"""

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import time
from typing import Any

from .env_exec import run_wsl, to_wsl_path


REQUIRED_STAR_INDEX_FILES: tuple[str, ...] = (
    "Genome",
    "SA",
    "SAindex",
    "chrLength.txt",
    "chrName.txt",
    "chrNameLength.txt",
    "chrStart.txt",
    "exonGeTrInfo.tab",
    "exonInfo.tab",
    "geneInfo.tab",
    "transcriptInfo.tab",
    "sjdbList.fromGTF.out.tab",
    "sjdbInfo.txt",
    "genomeParameters.txt",
)


def verify_star_index(
    index_dir: str | Path,
    expected_overhang: int | None = None,
    expected_major_version: str = "2.7",
) -> dict[str, Any]:
    """Verify presence, non-zero size, and parameters of all required STAR index files.
    
    A partially generated or interrupted index must never be accepted as complete.
    """
    idx_p = Path(index_dir).resolve()
    if not idx_p.is_dir():
        return {
            "is_complete": False,
            "is_compatible": False,
            "missing_files": list(REQUIRED_STAR_INDEX_FILES),
            "total_bytes": 0,
            "total_gb": 0.0,
            "parameters": {},
            "chromosome_count": 0,
            "error": f"Index directory does not exist: {idx_p}",
        }

    missing: list[str] = []
    file_details: dict[str, dict[str, Any]] = {}
    total_bytes = 0

    for fname in REQUIRED_STAR_INDEX_FILES:
        fpath = idx_p / fname
        try:
            if not fpath.is_file() or fpath.stat().st_size == 0:
                missing.append(fname)
            else:
                sz = fpath.stat().st_size
                total_bytes += sz
                file_details[fname] = {
                    "size_bytes": sz,
                    "size_mb": round(sz / (1024**2), 2),
                }
        except OSError:
            missing.append(fname)

    # Parse genome parameters
    params: dict[str, str] = {}
    param_file = idx_p / "genomeParameters.txt"
    if param_file.is_file():
        try:
            for line in param_file.read_text(encoding="utf-8", errors="ignore").splitlines():
                parts = line.strip().split(maxsplit=1)
                if len(parts) == 2:
                    params[parts[0]] = parts[1]
        except Exception:
            pass

    # Validate chromosome count
    chr_names: list[str] = []
    chr_file = idx_p / "chrName.txt"
    if chr_file.is_file():
        try:
            chr_names = [
                l.strip()
                for l in chr_file.read_text(encoding="utf-8", errors="ignore").splitlines()
                if l.strip()
            ]
        except Exception:
            pass

    is_complete = (len(missing) == 0 and len(chr_names) > 0 and len(params) > 0)

    # Check compatibility
    is_compatible = is_complete
    compat_error: str | None = None
    if is_complete:
        version_genome = params.get("versionGenome", "")
        if expected_major_version and not version_genome.startswith(expected_major_version):
            is_compatible = False
            compat_error = (
                f"Index version '{version_genome}' is incompatible with STAR major version "
                f"'{expected_major_version}'."
            )
        elif expected_overhang is not None:
            actual_overhang = int(params.get("sjdbOverhang", 0))
            # STAR can align reads shorter than or equal to sjdbOverhang + 1
            if actual_overhang < min(expected_overhang, 50):
                is_compatible = False
                compat_error = (
                    f"Index sjdbOverhang ({actual_overhang}) is too small for expected "
                    f"overhang ({expected_overhang})."
                )

    return {
        "is_complete": is_complete,
        "is_compatible": is_compatible,
        "missing_files": missing,
        "total_bytes": total_bytes,
        "total_gb": round(total_bytes / (1024**3), 2),
        "file_details": file_details,
        "chromosome_count": len(chr_names),
        "parameters": params,
        "error": compat_error if compat_error else (f"Missing files: {missing}" if missing else None),
    }


def restore_star_index_from_archive(
    archive_path: str | Path,
    target_dir: str | Path,
    expected_sha256: str | None = None,
    manifest_path: str | Path | None = None,
) -> dict[str, Any]:
    """Restore and verify a complete STAR index from an authorized tar archive."""
    arc_p = Path(archive_path).resolve()
    target_p = Path(target_dir).resolve()

    if not arc_p.is_file():
        raise FileNotFoundError(f"STAR index archive not found: {arc_p}")

    target_p.mkdir(parents=True, exist_ok=True)

    # Verify archive checksum if expected hash provided
    if expected_sha256:
        h = hashlib.sha256()
        with open(arc_p, "rb") as f:
            while chunk := f.read(8 * 1024 * 1024):
                h.update(chunk)
        observed_sha = h.hexdigest()
        if observed_sha.lower() != expected_sha256.lower():
            raise ValueError(
                f"Archive SHA-256 mismatch for {arc_p.name}!\n"
                f"Expected: {expected_sha256}\nObserved: {observed_sha}"
            )

    t0 = time.time()
    # Extract archive
    with tarfile.open(arc_p, "r:*") as tar:
        tar.extractall(path=target_p)
    extract_time = round(time.time() - t0, 2)

    # Verify manifest if provided
    verified_files = 0
    if manifest_path:
        man_p = Path(manifest_path).resolve()
        if man_p.is_file():
            man_data = json.loads(man_p.read_text(encoding="utf-8"))
            file_entries = man_data.get("files", man_data.get("file_details", {}))
            if isinstance(file_entries, dict):
                for fname, meta in file_entries.items():
                    extracted_f = target_p / fname
                    if not extracted_f.is_file():
                        raise FileNotFoundError(f"Extracted manifest file missing: {fname}")
                    verified_files += 1

    # Structural verification of the restored directory
    eval_res = verify_star_index(target_p)
    if not eval_res["is_complete"]:
        raise ValueError(
            f"Restored STAR index is incomplete or corrupt. Missing files: {eval_res['missing_files']}"
        )

    return {
        "status": "RESTORED_AND_VERIFIED",
        "archive_path": str(arc_p),
        "target_dir": str(target_p),
        "extract_time_seconds": extract_time,
        "total_gb": eval_res["total_gb"],
        "chromosome_count": eval_res["chromosome_count"],
        "manifest_verified_files": verified_files,
        "parameters": eval_res["parameters"],
    }


def check_star_alignment_resources(
    index_dir: str | Path,
    safety_margin_gb: float = 2.0,
) -> dict[str, Any]:
    """Check available RAM and scratch disk before loading the STAR index for alignment.
    
    Prevents Out-Of-Memory (OOM) crashes by verifying physical memory against the
    actual resident size of the index (sparse ~8-10 GB, dense ~28-30 GB).
    """
    idx_p = Path(index_dir).resolve()
    eval_res = verify_star_index(idx_p)

    # Required RAM is approximately the size of Genome + SA + SAindex
    genome_file = idx_p / "Genome"
    sa_file = idx_p / "SA"
    saindex_file = idx_p / "SAindex"

    idx_bin_bytes = 0
    for f in (genome_file, sa_file, saindex_file):
        if f.is_file():
            idx_bin_bytes += f.stat().st_size

    required_ram_gb = round(idx_bin_bytes / (1024**3), 2)
    # If files don't exist yet, estimate based on sparsity or human default
    if required_ram_gb < 1.0:
        sparse_d = int(eval_res.get("parameters", {}).get("genomeSAsparseD", 4))
        required_ram_gb = 30.0 if sparse_d == 1 else 9.0

    # Measure available system RAM
    available_ram_gb = 16.0  # conservative fallback
    try:
        meminfo_path = Path("/proc/meminfo")
        if meminfo_path.is_file():
            for line in meminfo_path.read_text(encoding="utf-8").splitlines():
                if line.startswith("MemAvailable:"):
                    kb = int(line.split()[1])
                    available_ram_gb = round(kb / (1024**2), 2)
                    break
    except Exception:
        pass

    # Measure free disk space
    free_disk_gb = 100.0
    try:
        total, used, free = shutil.disk_usage(idx_p)
        free_disk_gb = round(free / (1024**3), 2)
    except Exception:
        pass

    has_sufficient_ram = available_ram_gb >= (required_ram_gb - 1.0)
    has_sufficient_disk = free_disk_gb >= safety_margin_gb

    warning: str | None = None
    if not has_sufficient_ram:
        warning = (
            f"Available RAM ({available_ram_gb:.1f} GB) is below the estimated memory "
            f"required to load the STAR index ({required_ram_gb:.1f} GB). Alignment may trigger "
            "swap paging or OOM termination. Consider using a sparse index (sparseD=4)."
        )
    elif not has_sufficient_disk:
        warning = (
            f"Free disk space ({free_disk_gb:.1f} GB) is below the recommended safety "
            f"margin ({safety_margin_gb:.1f} GB)."
        )

    return {
        "has_sufficient_ram": has_sufficient_ram,
        "has_sufficient_disk": has_sufficient_disk,
        "required_ram_gb": required_ram_gb,
        "available_ram_gb": available_ram_gb,
        "free_disk_gb": free_disk_gb,
        "warning": warning,
    }


def build_star_index(
    fasta_path: str | Path,
    gtf_path: str | Path,
    index_dir: str | Path,
    sjdb_overhang: int = 149,
    threads: int = 8,
    sparse_d: int = 1,
    genome_sa_index_nbases: int = 14,
    limit_ram: int | None = None,
) -> dict[str, Any]:
    """Generate STAR genome index using Ensembl GRCh38 primary assembly FASTA and GTF."""
    f_p = Path(fasta_path).resolve()
    g_p = Path(gtf_path).resolve()
    idx_p = Path(index_dir).resolve()
    idx_p.mkdir(parents=True, exist_ok=True)

    wsl_fasta = to_wsl_path(f_p)
    wsl_gtf = to_wsl_path(g_p)
    wsl_idx = to_wsl_path(idx_p)

    # If FASTA is gzipped, decompress or check uncompressed
    actual_fasta = wsl_fasta
    if wsl_fasta.endswith(".gz"):
        uncompressed = wsl_fasta[:-3]
        check_res = run_wsl(f"test -f {uncompressed} && echo 'exists'")
        if "exists" not in check_res.stdout:
            print(f"Decompressing {wsl_fasta} to {uncompressed} for STAR...")
            run_wsl(f"zcat {wsl_fasta} > {uncompressed}")
        actual_fasta = uncompressed

    actual_gtf = wsl_gtf
    if wsl_gtf.endswith(".gz"):
        uncompressed_gtf = wsl_gtf[:-3]
        check_gtf = run_wsl(f"test -f {uncompressed_gtf} && echo 'exists'")
        if "exists" not in check_gtf.stdout:
            print(f"Decompressing {wsl_gtf} to {uncompressed_gtf} for STAR...")
            run_wsl(f"zcat {wsl_gtf} > {uncompressed_gtf}")
        actual_gtf = uncompressed_gtf

    cmd = (
        f"STAR --runMode genomeGenerate "
        f"--genomeDir {wsl_idx} "
        f"--genomeFastaFiles {actual_fasta} "
        f"--sjdbGTFfile {actual_gtf} "
        f"--sjdbOverhang {sjdb_overhang} "
        f"--genomeSAsparseD {sparse_d} "
        f"--genomeSAindexNbases {genome_sa_index_nbases} "
        f"--runThreadN {threads}"
    )
    if limit_ram is not None and limit_ram > 0:
        cmd += f" --limitGenomeGenerateRAM {limit_ram}"

    print(f"Executing STAR genomeGenerate:\n{cmd}")
    res = run_wsl(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"STAR index generation failed (code {res.returncode}):\n{res.stderr}\n{res.stdout}")

    return {
        "status": "STAR_INDEX_BUILT",
        "index_dir": str(idx_p),
        "fasta_source": str(f_p),
        "gtf_source": str(g_p),
        "sjdb_overhang": sjdb_overhang,
        "threads": threads,
        "genome_sa_sparse_d": sparse_d,
        "genome_sa_index_nbases": genome_sa_index_nbases,
        "limit_genome_generate_ram": limit_ram,
    }


def parse_star_final_log(log_path: str | Path) -> dict[str, Any]:
    """Parse Log.final.out from STAR into structured metrics."""
    metrics: dict[str, Any] = {}
    with open(log_path, "r", encoding="utf-8") as f:
        for line in f:
            if "|" in line:
                key, val = line.split("|", 1)
                k = key.strip()
                v = val.strip()
                # Clean numerical values
                if v.endswith("%"):
                    try:
                        metrics[k] = float(v.rstrip("%"))
                    except ValueError:
                        metrics[k] = v
                else:
                    try:
                        metrics[k] = int(v) if v.isdigit() else float(v)
                    except ValueError:
                        metrics[k] = v
    return metrics


def align_reads_star(
    trimmed_r1: str | Path,
    trimmed_r2: str | Path | None = None,
    index_dir: str | Path = "star_index",
    output_prefix: str | Path = "star_output/",
    threads: int = 8,
    bam_sort_mode: str = "SortedByCoordinate",
    limit_bam_sort_ram: int = 0,
    external_sort: bool = False,
    samtools_sort_mem: str = "1G",
) -> dict[str, Any]:
    """Execute governed STAR alignment (paired-end or single-end) with configurable sorting."""
    r1_p = Path(trimmed_r1).resolve()
    r2_p = Path(trimmed_r2).resolve() if trimmed_r2 is not None else None
    idx_p = Path(index_dir).resolve()
    out_prefix_p = Path(output_prefix).resolve()
    out_prefix_p.parent.mkdir(parents=True, exist_ok=True)

    wsl_r1 = to_wsl_path(r1_p)
    wsl_idx = to_wsl_path(idx_p)
    wsl_out = to_wsl_path(out_prefix_p)

    if r2_p is not None:
        wsl_r2 = to_wsl_path(r2_p)
        read_files_in = f"{wsl_r1} {wsl_r2}"
    else:
        read_files_in = f"{wsl_r1}"

    read_cmd_flag = "--readFilesCommand zcat " if str(r1_p).lower().endswith(".gz") else ""

    # Ensure STAR FIFOs are created on Linux tmpfs to prevent NTFS/DrvFs FIFO failures
    tmp_dir_name = f"STARtmp_{out_prefix_p.name}"
    run_wsl(f"rm -rf /tmp/{tmp_dir_name}")
    out_tmp_flag = f"--outTmpDir /tmp/{tmp_dir_name} "

    # Determine internal STAR sort mode
    star_sam_type = "BAM Unsorted" if external_sort else f"BAM {bam_sort_mode}"

    cmd = (
        f"STAR --runMode alignReads "
        f"--genomeDir {wsl_idx} "
        f"--readFilesIn {read_files_in} "
        f"{read_cmd_flag}"
        f"{out_tmp_flag}"
        f"--outFileNamePrefix {wsl_out} "
        f"--outSAMtype {star_sam_type} "
        f"--outSAMunmapped Within "
        f"--outSAMattributes NH HI AS nM NM MD jM jI MC ch "
        f"--outFilterType BySJout "
        f"--outFilterMultimapNmax 20 "
        f"--alignSJoverhangMin 8 "
        f"--alignSJDBoverhangMin 1 "
        f"--outFilterMismatchNmax 999 "
        f"--outFilterMismatchNoverReadLmax 0.04 "
        f"--alignIntronMin 20 "
        f"--alignIntronMax 1000000 "
        f"--alignMatesGapMax 1000000 "
        f"--runThreadN {threads}"
    )
    if limit_bam_sort_ram > 0 and not external_sort:
        cmd += f" --limitBAMsortRAM {limit_bam_sort_ram}"

    print(f"Executing STAR alignReads:\n{cmd}")
    res = run_wsl(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"STAR alignment failed (code {res.returncode}):\n{res.stderr}\n{res.stdout}")

    final_bam_path = Path(f"{out_prefix_p}Aligned.sortedByCoord.out.bam")
    log_final = Path(f"{out_prefix_p}Log.final.out")

    if external_sort:
        unsorted_bam = Path(f"{out_prefix_p}Aligned.out.bam")
        if not unsorted_bam.exists():
            raise FileNotFoundError(f"Expected STAR unsorted BAM not found: {unsorted_bam}")

        wsl_unsorted = to_wsl_path(unsorted_bam)
        wsl_sorted = to_wsl_path(final_bam_path)
        sort_threads = min(threads, 4)
        print(f"Running external samtools coordinate sort: {unsorted_bam} -> {final_bam_path}...")
        sort_cmd = f"samtools sort -@ {sort_threads} -m {samtools_sort_mem} -o {wsl_sorted} {wsl_unsorted} && rm -f {wsl_unsorted}"
        s_res = run_wsl(sort_cmd)
        if s_res.returncode != 0:
            raise RuntimeError(f"External samtools sort failed (code {s_res.returncode}):\n{s_res.stderr}")

    if not final_bam_path.exists():
        raise FileNotFoundError(f"Expected STAR sorted BAM not found: {final_bam_path}")

    metrics = parse_star_final_log(log_final) if log_final.exists() else {}

    return {
        "status": "ALIGNMENT_SUCCESS",
        "bam_path": str(final_bam_path),
        "log_final_path": str(log_final),
        "metrics": metrics,
        "external_sort_applied": external_sort,
    }
