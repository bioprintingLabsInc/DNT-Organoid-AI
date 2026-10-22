"""RSeQC infer_experiment.py strandedness diagnostic wrapper."""

import re
import subprocess
from pathlib import Path
from typing import Any

from .adapter import sha256_file
from .env_exec import run_wsl, to_wsl_path


def run_infer_experiment(
    bam_path: str | Path,
    bed_path: str | Path,
    sample_size: int = 200000,
    mapq_min: int = 10,
    output_txt_path: str | Path | None = None,
) -> dict[str, Any]:
    """Execute RSeQC infer_experiment.py with governed parameters."""
    b_p = Path(bam_path).resolve()
    bed_p = Path(bed_path).resolve()

    if not b_p.exists():
        raise FileNotFoundError(f"Input BAM not found: {b_p}")
    if not bed_p.exists():
        raise FileNotFoundError(f"Input BED not found: {bed_p}")

    bam_sha = sha256_file(b_p)
    bed_sha = sha256_file(bed_p)

    wsl_bam = to_wsl_path(b_p)
    wsl_bed = to_wsl_path(bed_p)

    cmd = f"infer_experiment.py -r {wsl_bed} -i {wsl_bam} -s {sample_size} -q {mapq_min}"
    res = run_wsl(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"infer_experiment.py failed (code {res.returncode}):\n{res.stderr}\n{res.stdout}")

    raw_output = res.stdout.strip()
    if output_txt_path:
        out_p = Path(output_txt_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        out_p.write_text(raw_output, encoding="utf-8")

    # Parse fractions:
    # Fraction of reads failed to determine: 0.0412
    # Fraction of reads explained by "1++,1--,2+-,2-+": 0.0125
    # Fraction of reads explained by "1+-,1-+,2++,2--": 0.9463
    failed_frac = None
    sense_frac = None
    antisense_frac = None

    for line in raw_output.splitlines():
        if "failed to determine" in line:
            m = re.search(r":\s*([0-9.]+)", line)
            if m:
                failed_frac = float(m.group(1))
        elif '1++,1--,2+-,2-+' in line:
            m = re.search(r":\s*([0-9.]+)", line)
            if m:
                sense_frac = float(m.group(1))
        elif '1+-,1-+,2++,2--' in line:
            m = re.search(r":\s*([0-9.]+)", line)
            if m:
                antisense_frac = float(m.group(1))

    return {
        "status": "INFER_EXPERIMENT_SUCCESS",
        "bam_path": str(b_p),
        "bam_sha256": bam_sha,
        "bed_path": str(bed_p),
        "bed_sha256": bed_sha,
        "sample_size": sample_size,
        "mapq_min": mapq_min,
        "raw_output": raw_output,
        "forward_fraction": sense_frac,
        "reverse_fraction": antisense_frac,
        "unexplained_fraction": failed_frac,
    }
