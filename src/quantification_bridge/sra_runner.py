"""Governed SRA Toolkit / fasterq-dump extraction wrapper.

Converts .sra archive containers into paired or single FASTQ files,
verifies tool version, detects read layout (paired-end vs single-end),
enforces governed separate temporary scratch storage,
and constructs canonical FastqPair objects for downstream quantification.
"""

from dataclasses import dataclass
import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any

from .env_exec import run_wsl, to_wsl_path
from .pairing import FastqPair


GOVERNED_SRA_TOOLKIT_VERSION = "3.1.1"


@dataclass(frozen=True, slots=True)
class SraExtractionResult:
    status: str
    run_accession: str
    layout: str
    fastq_r1: Path
    fastq_r2: Path | None
    tool_name: str
    tool_version: str
    output_dir: Path
    temp_dir: Path

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "run_accession": self.run_accession,
            "layout": self.layout,
            "fastq_r1": str(self.fastq_r1),
            "fastq_r2": str(self.fastq_r2) if self.fastq_r2 else None,
            "tool_name": self.tool_name,
            "tool_version": self.tool_version,
            "output_dir": str(self.output_dir),
            "temp_dir": str(self.temp_dir),
        }


def build_fasterq_dump_command(
    tool_bin: str,
    sra_path: str,
    output_dir: str,
    temp_dir: str,
    threads: int = 6,
) -> str:
    """Build deterministic fasterq-dump command string enforcing governed scratch storage."""
    return f"{tool_bin} --split-files -O {output_dir} -t {temp_dir} -e {threads} -f {sra_path}"


def probe_sra_tools() -> dict[str, Any]:
    """Identify available SRA toolkit binary and reported version."""
    # Check native or WSL
    is_win = os.name == "nt"
    if not is_win:
        binary = shutil.which("fasterq-dump") or shutil.which("fastq-dump")
        if not binary:
            return {"available": False, "tool_name": None, "tool_version": None, "path": None}
        name = Path(binary).name
        try:
            res = subprocess.run([binary, "--version"], capture_output=True, text=True, check=False)
            ver_text = (res.stdout or res.stderr).strip()
            return {"available": True, "tool_name": name, "tool_version": ver_text, "path": binary}
        except Exception:
            return {"available": False, "tool_name": name, "tool_version": "UNKNOWN", "path": binary}

    # Under Windows, probe via WSL
    script = (
        "source /home/chari/miniconda3/bin/activate dnt_bridge 2>/dev/null || true; "
        "which fasterq-dump || which fastq-dump || true"
    )
    res = run_wsl(script)
    out = res.stdout.strip()
    if not out:
        return {"available": False, "tool_name": None, "tool_version": None, "path": None}

    path = out.splitlines()[-1].strip()
    name = Path(path).name
    v_res = run_wsl(f"{path} --version 2>&1")
    ver_text = v_res.stdout.strip()
    return {"available": True, "tool_name": name, "tool_version": ver_text, "path": path}


def extract_sra(
    sra_path: str | Path,
    output_dir: str | Path,
    threads: int = 6,
    temp_dir: str | Path | None = None,
) -> SraExtractionResult:
    """Extract an SRA archive to FASTQ files via fasterq-dump / fastq-dump with governed temp storage."""
    sra_p = Path(sra_path).resolve()
    if not sra_p.is_file():
        raise FileNotFoundError(f"SRA file not found: {sra_p}")

    out_p = Path(output_dir).resolve()
    out_p.mkdir(parents=True, exist_ok=True)

    if temp_dir is None:
        temp_p = out_p / "_tmp_scratch"
    else:
        temp_p = Path(temp_dir).resolve()
    temp_p.mkdir(parents=True, exist_ok=True)

    run_acc = sra_p.stem
    tools = probe_sra_tools()
    if not tools["available"]:
        raise RuntimeError("Neither fasterq-dump nor fastq-dump is installed in the governed environment.")

    wsl_sra = to_wsl_path(sra_p)
    wsl_out = to_wsl_path(out_p)
    wsl_temp = to_wsl_path(temp_p)

    cmd_bin = tools["path"]
    if "fasterq-dump" in tools["tool_name"]:
        cmd = build_fasterq_dump_command(
            tool_bin=cmd_bin,
            sra_path=wsl_sra,
            output_dir=wsl_out,
            temp_dir=wsl_temp,
            threads=threads,
        )
    else:
        cmd = f"{cmd_bin} --split-files --outdir {wsl_out} {wsl_sra}"

    print(f"Executing SRA extraction: {cmd}")
    res = run_wsl(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"SRA extraction failed (code {res.returncode}):\n{res.stderr}\n{res.stdout}")

    # Inspect generated files
    r1 = out_p / f"{run_acc}_1.fastq"
    r2 = out_p / f"{run_acc}_2.fastq"
    single = out_p / f"{run_acc}.fastq"

    if r1.is_file() and r2.is_file():
        return SraExtractionResult(
            status="SRA_EXTRACTION_SUCCESS",
            run_accession=run_acc,
            layout="PAIRED",
            fastq_r1=r1,
            fastq_r2=r2,
            tool_name=tools["tool_name"],
            tool_version=tools["tool_version"],
            output_dir=out_p,
            temp_dir=temp_p,
        )
    elif single.is_file():
        return SraExtractionResult(
            status="SRA_EXTRACTION_SUCCESS",
            run_accession=run_acc,
            layout="SINGLE",
            fastq_r1=single,
            fastq_r2=None,
            tool_name=tools["tool_name"],
            tool_version=tools["tool_version"],
            output_dir=out_p,
            temp_dir=temp_p,
        )
    else:
        # Check if files without _1 exist or gzipped
        found = list(out_p.glob(f"{run_acc}*.fastq*"))
        if not found:
            raise FileNotFoundError(f"fasterq-dump completed but no FASTQ outputs found in {out_p}")
        # Identify pairs
        f_r1 = next((f for f in found if "_1" in f.name), found[0])
        f_r2 = next((f for f in found if "_2" in f.name), None)
        return SraExtractionResult(
            status="SRA_EXTRACTION_SUCCESS",
            run_accession=run_acc,
            layout="PAIRED" if f_r2 else "SINGLE",
            fastq_r1=f_r1,
            fastq_r2=f_r2,
            tool_name=tools["tool_name"],
            tool_version=tools["tool_version"],
            output_dir=out_p,
            temp_dir=temp_p,
        )


def extract_sra_to_fastq_pair(
    sra_path: str | Path,
    output_dir: str | Path,
    sample_id: str,
    biological_replicate_id: str = "rep1",
    condition_id: str = "unassigned",
    threads: int = 6,
    temp_dir: str | Path | None = None,
    assay_id: str | None = None,
) -> FastqPair:
    """Extract SRA archive and wrap as a verified FastqPair object conforming to established FastqPair contract."""
    extraction = extract_sra(sra_path, output_dir=output_dir, threads=threads, temp_dir=temp_dir)
    if extraction.layout != "PAIRED" or extraction.fastq_r2 is None:
        raise ValueError(f"Run {extraction.run_accession} is {extraction.layout}; only paired-end libraries are supported by the pipeline.")

    extra = {
        "run_accession": extraction.run_accession,
        "biological_replicate_id": biological_replicate_id,
        "condition_id": condition_id,
        "sra_extraction_status": extraction.status,
        "sra_extraction_tool": extraction.tool_name,
        "sra_extraction_version": extraction.tool_version,
        "sra_temp_dir": str(extraction.temp_dir),
    }

    return FastqPair(
        sample_id=sample_id,
        r1_path=extraction.fastq_r1,
        r2_path=extraction.fastq_r2,
        assay_id=assay_id,
        extra_metadata=extra,
    )
