"""STAR 2.7.11b index generation and alignment wrapper."""

import json
import re
import subprocess
from pathlib import Path
from typing import Any


from .env_exec import run_wsl, to_wsl_path


def build_star_index(
    fasta_path: str | Path,
    gtf_path: str | Path,
    index_dir: str | Path,
    sjdb_overhang: int = 149,
    threads: int = 8,
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
        # Check if uncompressed already exists
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
        f"--runThreadN {threads}"
    )

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
    trimmed_r2: str | Path,
    index_dir: str | Path,
    output_prefix: str | Path,
    threads: int = 8,
) -> dict[str, Any]:
    """Execute governed STAR paired-end alignment."""
    r1_p = Path(trimmed_r1).resolve()
    r2_p = Path(trimmed_r2).resolve()
    idx_p = Path(index_dir).resolve()
    out_prefix_p = Path(output_prefix).resolve()
    out_prefix_p.parent.mkdir(parents=True, exist_ok=True)

    wsl_r1 = to_wsl_path(r1_p)
    wsl_r2 = to_wsl_path(r2_p)
    wsl_idx = to_wsl_path(idx_p)
    wsl_out = to_wsl_path(out_prefix_p)

    cmd = (
        f"STAR --runMode alignReads "
        f"--genomeDir {wsl_idx} "
        f"--readFilesIn {wsl_r1} {wsl_r2} "
        f"--readFilesCommand zcat "
        f"--outFileNamePrefix {wsl_out} "
        f"--outSAMtype BAM SortedByCoordinate "
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

    print(f"Executing STAR alignReads:\n{cmd}")
    res = run_wsl(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"STAR alignment failed (code {res.returncode}):\n{res.stderr}\n{res.stdout}")

    bam_path = Path(f"{out_prefix_p}Aligned.sortedByCoord.out.bam")
    log_final = Path(f"{out_prefix_p}Log.final.out")

    if not bam_path.exists():
        raise FileNotFoundError(f"Expected STAR sorted BAM not found: {bam_path}")

    metrics = parse_star_final_log(log_final) if log_final.exists() else {}

    return {
        "status": "ALIGNMENT_SUCCESS",
        "bam_path": str(bam_path),
        "log_final_path": str(log_final),
        "metrics": metrics,
    }
