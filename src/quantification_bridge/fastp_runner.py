"""fastp 0.23.4 execution wrapper with governed treatment-blind flags.

Supports:
- Conventional bulk RNA-seq (governed paired-end and single-end flags).
- Lexogen QuantSeq 3' mRNA-seq (FWD and REV) with adapter trimming, poly(A) tail
  trimming (--trim_poly_x), minimum length filtering, and quality trimming,
  while preserving 5' bases by default (relying on STAR soft-clipping for random primer mismatch).
"""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
from typing import Any

from .env_exec import run_wsl, to_wsl_path
from .library_type import LibraryType

# Official Lexogen QuantSeq Read 1 adapter sequence (Illumina TruSeq Read 1 adapter)
LEXOGEN_QUANTSEQ_FWD_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"


def run_fastp_preprocessing(
    raw_r1: str | Path,
    raw_r2: str | Path | None = None,
    out_r1: str | Path | None = None,
    out_r2: str | Path | None = None,
    json_report: str | Path | None = None,
    html_report: str | Path | None = None,
    sample_id: str = "sample",
    threads: int = 8,
    library_type: str | LibraryType = LibraryType.CONVENTIONAL_BULK,
    adapter_sequence: str | None = None,
    trim_poly_x: bool | None = None,
    poly_x_min_len: int = 10,
    trim_front1: int = 0,
    min_read_length: int | None = None,
    cut_tail: bool | None = None,
    cut_mean_quality: int = 20,
) -> dict[str, Any]:
    """Execute fastp with governed flags appropriate for library preparation.

    Parameters
    ----------
    raw_r1 : str | Path
        Path to raw FASTQ file (Read 1).
    raw_r2 : str | Path | None
        Path to raw FASTQ file (Read 2) for paired-end, or None for single-end.
    out_r1 : str | Path | None
        Path to output trimmed Read 1 FASTQ file.
    out_r2 : str | Path | None
        Path to output trimmed Read 2 FASTQ file.
    json_report : str | Path | None
        Path to fastp JSON execution metrics report.
    html_report : str | Path | None
        Path to fastp HTML execution metrics report.
    sample_id : str
        Sample identifier.
    threads : int
        CPU threads.
    library_type : str | LibraryType
        Recognized library preparation type (default: CONVENTIONAL_BULK).
    adapter_sequence : str | None
        Explicit adapter sequence override. For QuantSeq FWD, defaults to
        LEXOGEN_QUANTSEQ_FWD_ADAPTER if not specified.
    trim_poly_x : bool | None
        Whether to enable poly(A)/poly(X) tail trimming (--trim_poly_x).
        Defaults to True for QuantSeq, False for conventional bulk.
    poly_x_min_len : int
        Minimum poly(X) tail length to trim (default: 10).
    trim_front1 : int
        Number of bases to hard-trim from 5' end of Read 1. Defaults to 0;
        STAR soft-clipping handles random primer mismatch without discarding
        read length.
    min_read_length : int | None
        Minimum required read length after trimming. Defaults to 20 for QuantSeq,
        or None (filtering disabled) for conventional bulk.
    cut_tail : bool | None
        Whether to perform sliding window quality trimming at read tails.
        Defaults to True for QuantSeq, or False (filtering disabled) for conventional bulk.
    cut_mean_quality : int
        Mean quality threshold for sliding window trimming (default: 20).
    """
    lib_type = LibraryType.from_string(library_type)

    r1_p = Path(raw_r1).resolve()
    r2_p = Path(raw_r2).resolve() if raw_r2 is not None else None

    # Derive default output paths if not specified
    if out_r1 is not None:
        out1_p = Path(out_r1).resolve()
    else:
        out1_p = r1_p.parent / f"{sample_id}_trimmed_R1.fastq.gz"

    if r2_p is not None:
        if out_r2 is not None:
            out2_p = Path(out_r2).resolve()
        else:
            out2_p = r1_p.parent / f"{sample_id}_trimmed_R2.fastq.gz"
    else:
        out2_p = None

    json_p = Path(json_report).resolve() if json_report else out1_p.parent / f"{sample_id}_fastp.json"
    html_p = Path(html_report).resolve() if html_report else out1_p.parent / f"{sample_id}_fastp.html"

    out1_p.parent.mkdir(parents=True, exist_ok=True)
    json_p.parent.mkdir(parents=True, exist_ok=True)

    wsl_r1 = to_wsl_path(r1_p)
    wsl_out1 = to_wsl_path(out1_p)
    wsl_json = to_wsl_path(json_p)
    wsl_html = to_wsl_path(html_p)

    cmd_parts = [
        "fastp",
        f"--in1 {wsl_r1}",
        f"--out1 {wsl_out1}",
    ]

    if r2_p is not None and out2_p is not None:
        wsl_r2 = to_wsl_path(r2_p)
        wsl_out2 = to_wsl_path(out2_p)
        cmd_parts.extend([
            f"--in2 {wsl_r2}",
            f"--out2 {wsl_out2}",
        ])

    cmd_parts.extend([
        f"--json {wsl_json}",
        f"--html {wsl_html}",
        f'--report_title "{sample_id} fastp Preprocessing Report"',
    ])

    # Library-specific flag construction
    if lib_type.is_quantseq:
        # 1. Adapter trimming
        eff_adapter = adapter_sequence
        if eff_adapter is None and lib_type in (
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI,
        ):
            eff_adapter = LEXOGEN_QUANTSEQ_FWD_ADAPTER

        if eff_adapter:
            cmd_parts.append(f"--adapter_sequence {eff_adapter}")
        elif r2_p is not None:
            cmd_parts.append("--detect_adapter_for_pe")

        # 2. Poly(X) tail trimming (polyA trimming for 3' mRNA-seq)
        eff_trim_poly_x = True if trim_poly_x is None else bool(trim_poly_x)
        if eff_trim_poly_x:
            cmd_parts.extend(["--trim_poly_x", f"--poly_x_min_len {poly_x_min_len}"])

        # 3. Poly(G) trimming (Illumina two-color chemistry)
        cmd_parts.extend(["--trim_poly_g", f"--poly_g_min_len 10"])

        # 4. 5' trimming (configurable, default 0 - let STAR soft-clipping handle mismatches)
        if trim_front1 > 0:
            cmd_parts.append(f"--trim_front1 {trim_front1}")

        # 5. Length filtering
        eff_min_len = 20 if min_read_length is None else int(min_read_length)
        if eff_min_len > 0:
            cmd_parts.append(f"--length_required {eff_min_len}")
        else:
            cmd_parts.append("--disable_length_filtering")

        # 6. Quality filtering
        eff_cut_tail = True if cut_tail is None else bool(cut_tail)
        if eff_cut_tail:
            cmd_parts.extend(["--cut_tail", f"--cut_mean_quality {cut_mean_quality}"])
        else:
            cmd_parts.append("--disable_quality_filtering")

    else:
        # Standard conventional bulk RNA-seq: preserve exact historical behavior
        if r2_p is not None:
            cmd_parts.append("--detect_adapter_for_pe")

        cmd_parts.extend([
            "--trim_poly_g",
            "--poly_g_min_len 10",
            "--disable_quality_filtering",
            "--disable_length_filtering",
        ])

    cmd_parts.append(f"--thread {threads}")
    cmd = " ".join(cmd_parts)

    print(f"Executing fastp for {sample_id} ({lib_type.value}):\n{cmd}")
    res = run_wsl(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"fastp failed for {sample_id} (code {res.returncode}):\n{res.stderr}\n{res.stdout}")

    if not out1_p.exists():
        raise FileNotFoundError(f"fastp did not emit output file: {out1_p}")
    if out2_p is not None and not out2_p.exists():
        raise FileNotFoundError(f"fastp did not emit output file: {out2_p}")

    metrics = {}
    if json_p.exists():
        try:
            metrics = json.loads(json_p.read_text(encoding="utf-8"))
        except Exception as e:
            print(f"Warning: could not parse fastp JSON metrics for {sample_id}: {e}")

    return {
        "status": "FASTP_PREPROCESSING_SUCCESS",
        "sample_id": sample_id,
        "library_type": str(lib_type.value),
        "trimmed_r1": str(out1_p),
        "trimmed_r2": str(out2_p) if out2_p else None,
        "read_layout": "paired_end" if r2_p is not None else "single_end",
        "json_report": str(json_p),
        "html_report": str(html_p),
        "metrics_summary": metrics.get("summary", {}),
        "adapter_cutting": metrics.get("adapter_cutting", {}),
        "poly_g_trimming": metrics.get("poly_g_trimming", {}),
        "command_executed": cmd,
    }
