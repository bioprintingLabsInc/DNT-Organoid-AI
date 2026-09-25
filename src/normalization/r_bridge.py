"""Controlled headless R execution bridge for official Bioconductor DESeq2 size-factor estimation."""

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

from .config import NormalizationConfig
from .models import Finding, Severity

R_SCRIPT_PATH = Path("scripts/r/estimate_size_factors.R")


@dataclass(frozen=True, slots=True)
class REnvironmentInfo:
    """Exact captured runtime environment for the R / DESeq2 execution engine."""

    r_version: str
    bioc_version: str
    deseq2_version: str
    operating_system: str
    execution_timestamp: str
    environment_identity: str
    is_available: bool
    error_message: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def get_environment_identity() -> str:
    """Return deterministic hash of the R environment definition (renv.lock)."""
    renv_path = Path("scripts/r/renv.lock")
    if renv_path.exists():
        try:
            content = renv_path.read_bytes()
            return f"renv:{hashlib.sha256(content).hexdigest()[:16]}"
        except OSError:
            pass
    return "unlocked_environment"


def resolve_r_binary(r_binary: str = "Rscript") -> str | None:
    """Resolve Rscript executable from system PATH, direct file path, or standard user-level install directory."""
    if Path(r_binary).is_file():
        return str(Path(r_binary).resolve())
    path = shutil.which(r_binary)
    if path:
        return path
    if r_binary == "Rscript":
        if platform.system() == "Windows":
            cand_x64 = Path.home() / ".r" / "R-4.4.3" / "bin" / "x64" / "Rscript.exe"
            if cand_x64.exists():
                return str(cand_x64)
            cand_bin = Path.home() / ".r" / "R-4.4.3" / "bin" / "Rscript.exe"
            if cand_bin.exists():
                return str(cand_bin)
        else:
            cand_unix = Path.home() / ".r" / "R-4.4.3" / "bin" / "Rscript"
            if cand_unix.exists():
                return str(cand_unix)
    return None


def probe_r_environment(r_binary: str = "Rscript") -> REnvironmentInfo:
    """Probe the system R environment to inspect R, Bioconductor, and DESeq2 versions."""
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    env_identity = get_environment_identity()
    os_name = platform.system()

    # 1. Check if Rscript binary exists
    binary_path = resolve_r_binary(r_binary)
    if not binary_path:
        return REnvironmentInfo(
            r_version="missing",
            bioc_version="missing",
            deseq2_version="missing",
            operating_system=os_name,
            execution_timestamp=timestamp,
            environment_identity=env_identity,
            is_available=False,
            error_message=f"R executable '{r_binary}' was not found on system PATH.",
        )

    probe_code = (
        "suppressPackageStartupMessages({"
        "r_v <- as.character(getRversion());"
        "d_v <- tryCatch(as.character(packageVersion('DESeq2')), error=function(e) 'missing');"
        "b_v <- tryCatch({"
        "if (requireNamespace('BiocManager', quietly=TRUE)) {"
        "as.character(BiocManager::version())"
        "} else if (requireNamespace('BiocVersion', quietly=TRUE)) {"
        "sub('^([0-9]+\\\\.[0-9]+).*', '\\\\1', as.character(packageVersion('BiocVersion')))"
        "} else {"
        "'unavailable'"
        "}"
        "}, error=function(e) 'unavailable');"
        "cat(paste(r_v, b_v, d_v, sep='|'));"
        "})"
    )

    try:
        proc = subprocess.run(
            [binary_path, "--vanilla", "-e", probe_code],
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
    except (subprocess.SubprocessError, OSError) as e:
        return REnvironmentInfo(
            r_version="error",
            bioc_version="error",
            deseq2_version="error",
            operating_system=os_name,
            execution_timestamp=timestamp,
            environment_identity=env_identity,
            is_available=False,
            error_message=f"Failed to execute R probe: {e}",
        )

    if proc.returncode != 0:
        return REnvironmentInfo(
            r_version="error",
            bioc_version="error",
            deseq2_version="error",
            operating_system=os_name,
            execution_timestamp=timestamp,
            environment_identity=env_identity,
            is_available=False,
            error_message=f"R probe failed with code {proc.returncode}: {proc.stderr.strip()}",
        )

    output = proc.stdout.strip()
    parts = output.split("|")
    if len(parts) != 3:
        return REnvironmentInfo(
            r_version="unknown",
            bioc_version="unknown",
            deseq2_version="unknown",
            operating_system=os_name,
            execution_timestamp=timestamp,
            environment_identity=env_identity,
            is_available=False,
            error_message=f"Malformed R probe output: '{output}'",
        )

    r_v, b_v, d_v = parts[0], parts[1], parts[2]
    is_available = d_v != "missing"
    err = None if is_available else "Official Bioconductor package 'DESeq2' is not installed in the R environment."

    return REnvironmentInfo(
        r_version=r_v,
        bioc_version=b_v,
        deseq2_version=d_v,
        operating_system=os_name,
        execution_timestamp=timestamp,
        environment_identity=env_identity,
        is_available=is_available,
        error_message=err,
    )


def estimate_size_factors_r(
    matrix_columns: tuple[tuple[int, ...], ...],
    gene_ids: tuple[str, ...],
    sample_ids: tuple[str, ...],
    method: str = "ratio",
    config: NormalizationConfig | None = None,
    r_script_path: Path | str = R_SCRIPT_PATH,
) -> tuple[dict[str, float] | None, REnvironmentInfo, list[Finding], dict[str, Any]]:
    """Execute official Bioconductor DESeq2::estimateSizeFactorsForMatrix via headless R runner.

    Returns:
        tuple of:
        - size_factors dict mapping sample_id to float, or None on failure;
        - REnvironmentInfo capturing exact runtime versions;
        - list of Findings;
        - execution provenance dict.
    """
    findings: list[Finding] = []
    cfg = config or NormalizationConfig()
    provenance: dict[str, Any] = {
        "execution_method": method,
        "input_gene_count": len(gene_ids),
        "input_sample_count": len(sample_ids),
    }

    # 1. Probe R environment
    env_info = probe_r_environment(cfg.r_binary_path)
    provenance["r_environment"] = env_info.to_dict()

    if not env_info.is_available:
        findings.append(
            Finding(
                severity=Severity.ERROR,
                rule_id="r_environment_missing",
                entity_type="r_environment",
                entity_id=cfg.r_binary_path,
                path="r_environment",
                message=(
                    f"Official R/Bioconductor DESeq2 environment is unavailable: {env_info.error_message}. "
                    "Production normalization requires the official DESeq2 package (R 4.4.3, Bioconductor 3.20, DESeq2 1.46.0)."
                ),
            )
        )
        return None, env_info, findings, provenance

    # 2. Check version matches explicitly
    if env_info.r_version != cfg.expected_r_version:
        sev = Severity.ERROR if cfg.strict_version_check else Severity.WARNING
        findings.append(
            Finding(
                severity=sev,
                rule_id="r_version_mismatch",
                entity_type="r_environment",
                entity_id="R",
                path="r_version",
                message=f"Runtime R version '{env_info.r_version}' differs from locked baseline '{cfg.expected_r_version}'.",
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
                    f"Bioconductor version cannot be reliably determined from the runtime environment (reported '{env_info.bioc_version}'). "
                    "Production normalization requires a verified programmatic Bioconductor 3.20 installation."
                ),
            )
        )
    elif env_info.bioc_version != cfg.expected_bioc_version:
        sev = Severity.ERROR if cfg.strict_version_check else Severity.WARNING
        findings.append(
            Finding(
                severity=sev,
                rule_id="r_version_mismatch",
                entity_type="r_environment",
                entity_id="Bioconductor",
                path="bioc_version",
                message=f"Runtime Bioconductor version '{env_info.bioc_version}' differs from locked baseline '{cfg.expected_bioc_version}'.",
            )
        )
    if env_info.deseq2_version != cfg.expected_deseq2_version:
        sev = Severity.ERROR if cfg.strict_version_check else Severity.WARNING
        findings.append(
            Finding(
                severity=sev,
                rule_id="r_version_mismatch",
                entity_type="r_environment",
                entity_id="DESeq2",
                path="deseq2_version",
                message=f"Runtime DESeq2 version '{env_info.deseq2_version}' differs from locked baseline '{cfg.expected_deseq2_version}'.",
            )
        )

    if any(f.severity == Severity.ERROR for f in findings):
        return None, env_info, findings, provenance

    # 3. Serialize input matrix to temporary TSV
    script_file = Path(r_script_path)
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

    with tempfile.TemporaryDirectory(prefix="deseq2_norm_") as tmpdir:
        input_tsv = Path(tmpdir) / "eligible_counts.tsv"
        output_json = Path(tmpdir) / "deseq2_result.json"

        # Format TSV: row 0 = header ("gene_id\tsample1\tsample2...")
        with input_tsv.open("w", encoding="utf-8", newline="") as f:
            f.write("gene_id\t" + "\t".join(sample_ids) + "\n")
            for i, gene in enumerate(gene_ids):
                row_vals = [str(matrix_columns[j][i]) for j in range(len(sample_ids))]
                f.write(f"{gene}\t" + "\t".join(row_vals) + "\n")

        # Record input matrix checksum
        input_bytes = input_tsv.read_bytes()
        matrix_checksum = hashlib.sha256(input_bytes).hexdigest()
        provenance["input_matrix_sha256"] = matrix_checksum

        # 4. Invoke headless R runner
        r_bin = resolve_r_binary(cfg.r_binary_path) or cfg.r_binary_path
        cmd = [
            r_bin,
            "--vanilla",
            str(script_file.resolve()),
            "--input",
            str(input_tsv.resolve()),
            "--method",
            method,
            "--output",
            str(output_json.resolve()),
        ]

        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=cfg.r_timeout_seconds,
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
                    message=f"DESeq2 execution timed out after {cfg.r_timeout_seconds} seconds.",
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

        # 5. Parse output JSON
        if not output_json.exists():
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="r_execution_failed",
                    entity_type="r_process",
                    entity_id=None,
                    path="r_execution",
                    message=f"R runner exited with code {proc.returncode} without writing output JSON. stderr: {proc.stderr.strip()}",
                )
            )
            return None, env_info, findings, provenance

        try:
            with output_json.open("r", encoding="utf-8") as f:
                res_data = json.load(f)
        except json.JSONDecodeError as e:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="r_execution_failed",
                    entity_type="r_process",
                    entity_id=None,
                    path="r_output_parsing",
                    message=f"Malformed JSON output from R runner: {e}",
                )
            )
            return None, env_info, findings, provenance

        status = res_data.get("status")
        err_msg = res_data.get("error_message") or ""

        if status == "ERROR":
            # Check for standard ratio zero-geometric-mean failure
            if "every gene contains at least one zero" in err_msg.lower() or "log geometric means" in err_msg.lower():
                findings.append(
                    Finding(
                        severity=Severity.REVIEW,
                        rule_id="standard_size_factors_failed",
                        entity_type="normalization_method",
                        entity_id=method,
                        path="size_factor_estimation",
                        message=(
                            f"Official DESeq2 standard ratio estimation failed: {err_msg}. "
                            "Every gene contains at least one zero across samples in the cohort. "
                            "Automatic fallback to poscounts is prohibited; requires an explicit auditable NormalizationReviewDisposition."
                        ),
                    )
                )
            else:
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="r_execution_failed",
                        entity_type="r_process",
                        entity_id=None,
                        path="size_factor_estimation",
                        message=f"DESeq2 size-factor estimation failed: {err_msg}",
                    )
                )
            return None, env_info, findings, provenance

        raw_sfs = res_data.get("size_factors")
        if not isinstance(raw_sfs, dict):
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="r_execution_failed",
                    entity_type="r_output",
                    entity_id=None,
                    path="size_factors",
                    message="R runner reported success but returned non-dictionary size factors.",
                )
            )
            return None, env_info, findings, provenance

        size_factors: dict[str, float] = {}
        for s in sample_ids:
            if s not in raw_sfs:
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="missing_sample_size_factor",
                        entity_type="sample",
                        entity_id=s,
                        path=f"size_factors.{s}",
                        message=f"Sample '{s}' was not returned in DESeq2 size factor output.",
                    )
                )
                return None, env_info, findings, provenance
            val = raw_sfs[s]
            if val is None or not isinstance(val, (int, float)) or val <= 0:
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="invalid_sample_size_factor",
                        entity_type="sample",
                        entity_id=s,
                        path=f"size_factors.{s}",
                        message=f"Sample '{s}' received invalid non-positive or null size factor: {val}.",
                    )
                )
                return None, env_info, findings, provenance
            size_factors[s] = float(val)

        provenance["r_exit_code"] = proc.returncode
        provenance["output_size_factors"] = size_factors
        return size_factors, env_info, findings, provenance
