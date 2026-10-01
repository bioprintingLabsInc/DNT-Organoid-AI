"""Controlled headless R execution bridge for official Bioconductor DESeq2 differential expression."""

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile
from typing import Any

from src.normalization.r_bridge import REnvironmentInfo, probe_r_environment, resolve_r_binary
from src.response_builder.models import Finding, Severity

LOCKED_R_VERSION = "4.4.3"
LOCKED_BIOC_VERSION = "3.20"
LOCKED_DESEQ2_VERSION = "1.46.0"
R_CONTRAST_SCRIPT_PATH = Path("scripts/r/run_deseq2_contrast.R")


def get_r_script_checksum(script_path: Path | str = R_CONTRAST_SCRIPT_PATH) -> str:
    """Compute sha256 checksum of the R runner script."""
    p = Path(script_path)
    if p.exists():
        try:
            return hashlib.sha256(p.read_bytes()).hexdigest()
        except OSError:
            pass
    return "unknown"


def execute_deseq2_contrast_r(
    gene_ids: tuple[str, ...],
    sample_ids: tuple[str, ...],
    count_columns: tuple[tuple[int, ...], ...],
    coldata: list[dict[str, str]],
    size_factors: dict[str, float],
    design_formula: str,
    r_binary_path: str = "Rscript",
    strict_version_check: bool = True,
    r_timeout_seconds: int = 180,
    r_script_path: Path | str = R_CONTRAST_SCRIPT_PATH,
) -> tuple[dict[str, Any] | None, REnvironmentInfo, list[Finding], dict[str, Any]]:
    """Execute official Bioconductor DESeq2 differential expression via headless R runner.

    Returns:
        tuple of:
        - raw execution output dict or None on failure;
        - REnvironmentInfo capturing exact runtime versions;
        - list of Findings;
        - execution provenance dict.
    """
    findings: list[Finding] = []
    script_file = Path(r_script_path)
    script_checksum = get_r_script_checksum(script_file)

    provenance: dict[str, Any] = {
        "script_path": str(script_file),
        "script_sha256": script_checksum,
        "design_formula": design_formula,
        "input_gene_count": len(gene_ids),
        "input_sample_count": len(sample_ids),
        "strict_version_check": strict_version_check,
    }

    # 1. Probe R environment
    env_info = probe_r_environment(r_binary_path)
    provenance["r_environment"] = env_info.to_dict()

    if not env_info.is_available:
        findings.append(
            Finding(
                severity=Severity.ERROR,
                rule_id="r_environment_missing",
                entity_type="r_environment",
                entity_id=r_binary_path,
                path="r_environment",
                message=(
                    f"Official R/Bioconductor DESeq2 environment is unavailable: {env_info.error_message}. "
                    f"Step 7B differential expression requires official DESeq2 (R {LOCKED_R_VERSION}, "
                    f"Bioconductor {LOCKED_BIOC_VERSION}, DESeq2 {LOCKED_DESEQ2_VERSION})."
                ),
            )
        )
        return None, env_info, findings, provenance

    # 2. Strict version verification
    if env_info.r_version != LOCKED_R_VERSION:
        sev = Severity.ERROR if strict_version_check else Severity.WARNING
        findings.append(
            Finding(
                severity=sev,
                rule_id="r_version_mismatch",
                entity_type="r_environment",
                entity_id="R",
                path="r_version",
                message=f"Runtime R version '{env_info.r_version}' differs from locked baseline '{LOCKED_R_VERSION}'.",
            )
        )

    if env_info.bioc_version in ("unavailable", "unknown", "missing"):
        findings.append(
            Finding(
                severity=Severity.ERROR,
                rule_id="bioc_version_undetermined",
                entity_type="r_environment",
                entity_id="Bioconductor",
                path="bioc_version",
                message=(
                    f"Bioconductor version cannot be reliably determined (reported '{env_info.bioc_version}'). "
                    f"Step 7B differential expression requires Bioconductor {LOCKED_BIOC_VERSION}."
                ),
            )
        )
    elif env_info.bioc_version != LOCKED_BIOC_VERSION:
        sev = Severity.ERROR if strict_version_check else Severity.WARNING
        findings.append(
            Finding(
                severity=sev,
                rule_id="r_version_mismatch",
                entity_type="r_environment",
                entity_id="Bioconductor",
                path="bioc_version",
                message=f"Runtime Bioconductor version '{env_info.bioc_version}' differs from locked baseline '{LOCKED_BIOC_VERSION}'.",
            )
        )

    if env_info.deseq2_version != LOCKED_DESEQ2_VERSION:
        sev = Severity.ERROR if strict_version_check else Severity.WARNING
        findings.append(
            Finding(
                severity=sev,
                rule_id="r_version_mismatch",
                entity_type="r_environment",
                entity_id="DESeq2",
                path="deseq2_version",
                message=f"Runtime DESeq2 version '{env_info.deseq2_version}' differs from locked baseline '{LOCKED_DESEQ2_VERSION}'.",
            )
        )

    if any(f.severity == Severity.ERROR for f in findings):
        return None, env_info, findings, provenance

    # 3. Check script existence
    if not script_file.exists():
        findings.append(
            Finding(
                severity=Severity.ERROR,
                rule_id="r_script_missing",
                entity_type="script",
                entity_id=str(script_file),
                path="r_script_path",
                message=f"Headless R execution script '{script_file}' does not exist.",
            )
        )
        return None, env_info, findings, provenance

    # 4. Serialize input tables to temporary directory
    with tempfile.TemporaryDirectory(prefix="deseq2_de_") as tmpdir:
        tmp = Path(tmpdir)
        counts_tsv = tmp / "contrast_counts.tsv"
        coldata_tsv = tmp / "contrast_coldata.tsv"
        sizefactors_tsv = tmp / "contrast_sizefactors.tsv"
        output_json = tmp / "deseq2_output.json"

        # Write counts TSV
        with counts_tsv.open("w", encoding="utf-8", newline="") as f:
            f.write("gene_id\t" + "\t".join(sample_ids) + "\n")
            for i, gene in enumerate(gene_ids):
                row_vals = [str(count_columns[j][i]) for j in range(len(sample_ids))]
                f.write(f"{gene}\t" + "\t".join(row_vals) + "\n")

        # Write coldata TSV
        col_keys = list(coldata[0].keys())
        with coldata_tsv.open("w", encoding="utf-8", newline="") as f:
            f.write("\t".join(col_keys) + "\n")
            for row in coldata:
                f.write("\t".join(str(row.get(k, "")) for k in col_keys) + "\n")

        # Write size factors TSV
        with sizefactors_tsv.open("w", encoding="utf-8", newline="") as f:
            f.write("sample_id\tsize_factor\n")
            for sid in sample_ids:
                sf_val = size_factors[sid]
                f.write(f"{sid}\t{sf_val}\n")

        # Record input checksums
        counts_sha = hashlib.sha256(counts_tsv.read_bytes()).hexdigest()
        coldata_sha = hashlib.sha256(coldata_tsv.read_bytes()).hexdigest()
        sizefactors_sha = hashlib.sha256(sizefactors_tsv.read_bytes()).hexdigest()

        provenance["input_hashes"] = {
            "counts_matrix_sha256": counts_sha,
            "coldata_sha256": coldata_sha,
            "sizefactors_sha256": sizefactors_sha,
        }

        # 5. Invoke headless R runner
        r_bin = resolve_r_binary(r_binary_path) or r_binary_path
        cmd = [
            r_bin,
            "--vanilla",
            str(script_file.resolve()),
            "--counts",
            str(counts_tsv.resolve()),
            "--coldata",
            str(coldata_tsv.resolve()),
            "--sizefactors",
            str(sizefactors_tsv.resolve()),
            "--design",
            design_formula,
            "--output",
            str(output_json.resolve()),
        ]

        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=r_timeout_seconds,
                check=False,
            )
        except subprocess.TimeoutExpired:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="r_execution_timeout",
                    entity_type="r_process",
                    entity_id=None,
                    path="r_execution",
                    message=f"DESeq2 execution timed out after {r_timeout_seconds} seconds.",
                )
            )
            return None, env_info, findings, provenance
        except (subprocess.SubprocessError, OSError) as e:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="r_execution_failed",
                    entity_type="r_process",
                    entity_id=None,
                    path="r_execution",
                    message=f"Subprocess invocation failed: {e}",
                )
            )
            return None, env_info, findings, provenance

        provenance["r_returncode"] = proc.returncode
        if proc.stderr:
            provenance["r_stderr"] = proc.stderr.strip()[:2000]

        # 6. Parse output JSON
        if not output_json.exists():
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="r_output_missing",
                    entity_type="r_process",
                    entity_id=None,
                    path="r_output",
                    message=(
                        f"DESeq2 runner did not produce expected output JSON (exit code {proc.returncode}). "
                        f"Stderr: {proc.stderr.strip()[:500]}"
                    ),
                )
            )
            return None, env_info, findings, provenance

        try:
            with output_json.open("r", encoding="utf-8") as f:
                result_data = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="r_output_malformed",
                    entity_type="r_process",
                    entity_id=None,
                    path="r_output",
                    message=f"Failed to parse DESeq2 output JSON: {e}",
                )
            )
            return None, env_info, findings, provenance

        if result_data.get("status") != "SUCCESS":
            err_msg = result_data.get("error_message") or "Unknown error in DESeq2 runner."
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="deseq2_execution_error",
                    entity_type="deseq2_runner",
                    entity_id=design_formula,
                    path="r_execution",
                    message=f"DESeq2 execution failed in R: {err_msg}",
                )
            )
            return None, env_info, findings, provenance

        return result_data, env_info, findings, provenance
