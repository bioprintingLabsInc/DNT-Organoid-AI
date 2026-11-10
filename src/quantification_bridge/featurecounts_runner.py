"""Governed featureCounts (subread) execution runner and metric parser.

Enforces:
1. Paired-fragment counting (-p --countReadPairs -B -C).
2. Multimapper exclusion (no -M flag).
3. Ambiguous overlap exclusion (no -O flag).
4. Governed strandedness flag propagation (-s 0, -s 1, or -s 2).
5. Clean summary metric parsing and seamless adapter integration.
"""

from pathlib import Path
import re
import subprocess
from typing import Any

from .adapter import convert_featurecounts_to_matrix, sha256_file
from .env_exec import run_wsl, to_wsl_path
from .manifest import APPROVED_RUN_ACCESSIONS


def parse_featurecounts_summary(summary_path: str | Path) -> dict[str, dict[str, int]]:
    """Parse the .summary file emitted by featureCounts into sample-wise metric counts.

    Returns:
        Dictionary mapping BAM filename/header to metric dictionary:
        {
            "sample.bam": {
                "Assigned": int,
                "Unassigned_Unmapped": int,
                "Unassigned_MultiMapping": int,
                "Unassigned_NoFeatures": int,
                "Unassigned_Ambiguity": int,
                ...
            }
        }
    """
    p = Path(summary_path).resolve()
    if not p.is_file():
        return {}

    lines = [line.strip() for line in p.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) < 2:
        return {}

    headers = lines[0].split("\t")[1:]
    metrics: dict[str, dict[str, int]] = {h: {} for h in headers}

    for line in lines[1:]:
        parts = line.split("\t")
        category = parts[0].strip()
        for idx, val in enumerate(parts[1:]):
            header = headers[idx]
            try:
                metrics[header][category] = int(val.strip())
            except ValueError:
                pass

    return metrics


def run_featurecounts(
    bam_paths: list[str | Path] | tuple[str | Path, ...],
    gtf_path: str | Path,
    output_counts_path: str | Path,
    strandedness: int = 0,
    threads: int = 8,
    paired_end: bool = True,
    both_ends_mapped: bool = True,
    exclude_chimeric: bool = True,
    feature_type: str = "exon",
    attribute_type: str = "gene_id",
) -> dict[str, Any]:
    """Execute featureCounts under frozen scientific governance constraints.

    Flags applied:
    - -p: paired-end fragments
    - --countReadPairs: count 1 fragment instead of 2 mate reads
    - -B: both ends must be mapped
    - -C: do not count chimeric fragments (mates on different chromosomes)
    - -s <strandedness>: 0 (unstranded), 1 (forward), 2 (reverse)
    - -a <gtf_path>: reference annotation
    - -t <feature_type>: feature type in GTF (default: exon)
    - -g <attribute_type>: grouping attribute in GTF (default: gene_id)
    - -T <threads>: CPU threads
    - -o <output_counts_path>: primary counts table
    """
    if strandedness not in (0, 1, 2):
        raise ValueError(f"Invalid strandedness flag {strandedness}; must be 0 (unstranded), 1 (forward), or 2 (reverse).")

    out_p = Path(output_counts_path).resolve()
    out_p.parent.mkdir(parents=True, exist_ok=True)
    gtf_p = Path(gtf_path).resolve()

    if not gtf_p.exists():
        raise FileNotFoundError(f"GTF annotation not found: {gtf_p}")

    resolved_bams = [Path(b).resolve() for b in bam_paths]
    for b in resolved_bams:
        if not b.exists():
            raise FileNotFoundError(f"Input BAM file not found: {b}")

    wsl_gtf = to_wsl_path(gtf_p)
    wsl_out = to_wsl_path(out_p)
    wsl_bams = " ".join(to_wsl_path(b) for b in resolved_bams)

    # If GTF is gzipped, decompress or check uncompressed
    actual_gtf = wsl_gtf
    if wsl_gtf.endswith(".gz"):
        uncompressed_gtf = wsl_gtf[:-3]
        check_gtf = run_wsl(f"test -f {uncompressed_gtf} && echo 'exists'")
        if "exists" not in check_gtf.stdout:
            print(f"Decompressing {wsl_gtf} to {uncompressed_gtf} for featureCounts...")
            run_wsl(f"zcat {wsl_gtf} > {uncompressed_gtf}")
        actual_gtf = uncompressed_gtf

    cmd_parts = ["featureCounts"]
    if paired_end:
        cmd_parts.extend(["-p", "--countReadPairs"])
        if both_ends_mapped:
            cmd_parts.append("-B")
        if exclude_chimeric:
            cmd_parts.append("-C")

    cmd_parts.extend([
        "-s", str(strandedness),
        "-t", feature_type,
        "-g", attribute_type,
        "-T", str(threads),
        "-a", actual_gtf,
        "-o", wsl_out,
        wsl_bams,
    ])

    full_cmd = " ".join(cmd_parts)
    print(f"Executing featureCounts:\n{full_cmd}")
    res = run_wsl(full_cmd)
    if res.returncode != 0:
        raise RuntimeError(f"featureCounts execution failed (code {res.returncode}):\n{res.stderr}\n{res.stdout}")

    if not out_p.exists() or out_p.stat().st_size == 0:
        raise FileNotFoundError(f"featureCounts failed to emit non-empty output table: {out_p}")

    summary_file = Path(f"{out_p}.summary")
    summary_metrics = parse_featurecounts_summary(summary_file) if summary_file.exists() else {}

    return {
        "status": "FEATURECOUNTS_SUCCESS",
        "output_counts_path": str(out_p),
        "output_summary_path": str(summary_file) if summary_file.exists() else None,
        "strandedness_flag": strandedness,
        "threads": threads,
        "paired_end": paired_end,
        "input_bams_count": len(resolved_bams),
        "summary_metrics": summary_metrics,
        "output_sha256": sha256_file(out_p),
    }


def execute_featurecounts_and_adapter(
    bam_paths: list[str | Path] | tuple[str | Path, ...],
    gtf_path: str | Path,
    raw_featurecounts_path: str | Path,
    output_matrix_path: str | Path,
    output_universe_path: str | Path | None = None,
    approved_runs: tuple[str, ...] = APPROVED_RUN_ACCESSIONS,
    strandedness: int = 0,
    threads: int = 8,
    paired_end: bool = True,
) -> dict[str, Any]:
    """Execute featureCounts and convert directly to frozen pipeline CountMatrix format.

    Verifies 100% numerical conservation and zero schema discrepancies.
    """
    fc_res = run_featurecounts(
        bam_paths=bam_paths,
        gtf_path=gtf_path,
        output_counts_path=raw_featurecounts_path,
        strandedness=strandedness,
        threads=threads,
        paired_end=paired_end,
    )

    adapt_res = convert_featurecounts_to_matrix(
        raw_featurecounts_path=raw_featurecounts_path,
        output_matrix_path=output_matrix_path,
        output_universe_path=output_universe_path,
        approved_runs=approved_runs,
    )

    return {
        "status": "FEATURECOUNTS_AND_ADAPTER_SUCCESS",
        "featurecounts": fc_res,
        "adapter_conversion": adapt_res,
        "output_matrix_path": str(output_matrix_path),
        "is_numerically_conserved": adapt_res["is_numerically_conserved"],
        "total_genes_emitted": adapt_res["total_genes_emitted"],
        "approved_runs": list(approved_runs),
    }
