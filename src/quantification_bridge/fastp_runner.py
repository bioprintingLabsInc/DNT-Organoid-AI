"""fastp 0.23.4 execution wrapper with governed treatment-blind flags."""

import json
import subprocess
from pathlib import Path
from typing import Any


from .env_exec import run_wsl, to_wsl_path


def run_fastp_preprocessing(
    raw_r1: str | Path,
    raw_r2: str | Path,
    out_r1: str | Path,
    out_r2: str | Path,
    json_report: str | Path,
    html_report: str | Path,
    sample_id: str,
    threads: int = 8,
) -> dict[str, Any]:
    """Execute fastp with governed flags: adapter trimming, poly-G, quality/length filters disabled."""
    r1_p = Path(raw_r1).resolve()
    r2_p = Path(raw_r2).resolve()
    out1_p = Path(out_r1).resolve()
    out2_p = Path(out_r2).resolve()
    json_p = Path(json_report).resolve()
    html_p = Path(html_report).resolve()

    out1_p.parent.mkdir(parents=True, exist_ok=True)
    json_p.parent.mkdir(parents=True, exist_ok=True)

    wsl_r1 = to_wsl_path(r1_p)
    wsl_r2 = to_wsl_path(r2_p)
    wsl_out1 = to_wsl_path(out1_p)
    wsl_out2 = to_wsl_path(out2_p)
    wsl_json = to_wsl_path(json_p)
    wsl_html = to_wsl_path(html_p)

    cmd = (
        f"fastp "
        f"--in1 {wsl_r1} "
        f"--in2 {wsl_r2} "
        f"--out1 {wsl_out1} "
        f"--out2 {wsl_out2} "
        f"--json {wsl_json} "
        f"--html {wsl_html} "
        f'--report_title "{sample_id} fastp Preprocessing Report" '
        f"--detect_adapter_for_pe "
        f"--trim_poly_g "
        f"--poly_g_min_len 10 "
        f"--disable_quality_filtering "
        f"--disable_length_filtering "
        f"--thread {threads}"
    )

    print(f"Executing fastp for {sample_id}:\n{cmd}")
    res = run_wsl(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"fastp failed for {sample_id} (code {res.returncode}):\n{res.stderr}\n{res.stdout}")

    if not out1_p.exists() or not out2_p.exists():
        raise FileNotFoundError(f"fastp did not emit output files for {sample_id}")

    metrics = {}
    if json_p.exists():
        try:
            metrics = json.loads(json_p.read_text(encoding="utf-8"))
        except Exception as e:
            print(f"Warning: could not parse fastp JSON metrics for {sample_id}: {e}")

    return {
        "status": "FASTP_PREPROCESSING_SUCCESS",
        "sample_id": sample_id,
        "trimmed_r1": str(out1_p),
        "trimmed_r2": str(out2_p),
        "json_report": str(json_p),
        "html_report": str(html_p),
        "metrics_summary": metrics.get("summary", {}),
        "adapter_cutting": metrics.get("adapter_cutting", {}),
        "poly_g_trimming": metrics.get("poly_g_trimming", {}),
    }
