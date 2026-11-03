"""Execution profile governance, environment probing, and preflight validation.

Supports cross-laboratory portability:
- Institutional HPC servers, laboratory workstations, and evaluation laptops.
- Dynamic hardware resource auditing (CPU, physical RAM, swap, disk capacity).
- Non-hardcoded path resolution with environment variable overrides.
- Explicit distinction between fatal technical blockers and advisory quality warnings.
- Sparse STAR benchmarking protocol and multi-tier acceptance criteria.
"""

from dataclasses import dataclass, field
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
from typing import Any
import yaml

from .env_exec import run_wsl, to_wsl_path
from .manifest import (
    APPROVED_2_SAMPLE_TECHNICAL_PILOT,
    APPROVED_PILOT_RUNS,
    RUN_TO_SAMPLE_MAP,
    TARGET_BIOPROJECT,
    TARGET_DEVELOPMENTAL_DAY,
    TARGET_SRA_STUDY,
)


DEFAULT_PROFILES_PATH = Path(__file__).resolve().parents[2] / "config" / "quantification_bridge" / "execution_profiles.yaml"


@dataclass(frozen=True, slots=True)
class HardwareBaseline:
    """Captured hardware metrics of the host and execution environment."""

    cpu_cores: int
    cpu_threads: int
    cpu_model: str
    physical_ram_gb: float
    swap_gb: float
    data_dir_free_gb: float
    ref_dir_free_gb: float
    os_kernel: str
    is_wsl: bool
    architecture: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "cpu_cores": self.cpu_cores,
            "cpu_threads": self.cpu_threads,
            "cpu_model": self.cpu_model,
            "physical_ram_gb": round(self.physical_ram_gb, 2),
            "swap_gb": round(self.swap_gb, 2),
            "data_dir_free_gb": round(self.data_dir_free_gb, 2),
            "ref_dir_free_gb": round(self.ref_dir_free_gb, 2),
            "os_kernel": self.os_kernel,
            "is_wsl": self.is_wsl,
            "architecture": self.architecture,
        }


@dataclass(frozen=True, slots=True)
class ToolInventory:
    """Verified software versions and paths of bioinformatics utilities."""

    tools: dict[str, dict[str, str | None]]

    def is_complete(self) -> bool:
        """True if all core required binaries are available."""
        required = ["fastp", "STAR", "featureCounts", "samtools", "infer_experiment.py", "fasterq-dump"]
        return all(self.tools.get(t, {}).get("path") not in (None, "NOT_FOUND") for t in required)

    def to_dict(self) -> dict[str, Any]:
        return self.tools


@dataclass(frozen=True, slots=True)
class ExecutionProfile:
    """Configured compute parameters for raw-read RNA-seq processing."""

    name: str
    display_name: str
    description: str
    physical_ram_gb_min: float
    cpu_cores_min: int
    storage_gb_min: float
    genome_sa_sparse_d: int
    genome_sa_index_nbases: int
    limit_genome_generate_ram: int
    bam_sort_mode: str
    limit_bam_sort_ram: int
    threads: int
    external_sort_enabled: bool
    samtools_sort_mem_per_thread: str
    reproducibility_tier: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "display_name": self.display_name,
            "description": self.description,
            "physical_ram_gb_min": self.physical_ram_gb_min,
            "cpu_cores_min": self.cpu_cores_min,
            "storage_gb_min": self.storage_gb_min,
            "genome_sa_sparse_d": self.genome_sa_sparse_d,
            "genome_sa_index_nbases": self.genome_sa_index_nbases,
            "limit_genome_generate_ram": self.limit_genome_generate_ram,
            "bam_sort_mode": self.bam_sort_mode,
            "limit_bam_sort_ram": self.limit_bam_sort_ram,
            "threads": self.threads,
            "external_sort_enabled": self.external_sort_enabled,
            "samtools_sort_mem_per_thread": self.samtools_sort_mem_per_thread,
            "reproducibility_tier": self.reproducibility_tier,
        }


@dataclass(frozen=True, slots=True)
class PreflightFinding:
    """Structured preflight finding distinguishing fatal errors from advisory warnings."""

    severity: str  # "FATAL" or "WARNING"
    code: str
    message: str
    remediation: str | None = None


@dataclass(frozen=True, slots=True)
class PreflightCheckResult:
    """Comprehensive result of resource, software, and configuration verification."""

    profile: ExecutionProfile
    is_compatible: bool
    hardware: HardwareBaseline
    tools: ToolInventory
    findings: tuple[PreflightFinding, ...]
    provenance: dict[str, Any]

    @property
    def has_fatal_errors(self) -> bool:
        return any(f.severity == "FATAL" for f in self.findings)

    @property
    def warnings_count(self) -> int:
        return sum(1 for f in self.findings if f.severity == "WARNING")


def load_execution_profiles(config_path: Path | None = None) -> dict[str, Any]:
    """Load execution profiles YAML configuration."""
    target_path = Path(config_path or DEFAULT_PROFILES_PATH).resolve()
    if not target_path.is_file():
        raise FileNotFoundError(f"Execution profiles configuration not found: {target_path}")

    with open(target_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    return data


def get_execution_profile(
    name: str | None = None,
    config_path: Path | None = None,
) -> ExecutionProfile:
    """Retrieve and parse a specific execution profile.

    Resolution order:
    1. Explicit `name` argument.
    2. Environment variable `DNT_EXECUTION_PROFILE`.
    3. Auto-detected fallback based on physical RAM.
    """
    config = load_execution_profiles(config_path)
    profiles_dict = config.get("profiles", {})

    target_name = name or os.environ.get("DNT_EXECUTION_PROFILE")
    if not target_name:
        # Auto-detect sensible profile
        target_name = _auto_select_profile(profiles_dict)

    if target_name not in profiles_dict:
        available = list(profiles_dict.keys())
        raise ValueError(f"Unknown execution profile '{target_name}'. Available: {available}")

    p_data = profiles_dict[target_name]
    hw = p_data.get("hardware_minimum", {})
    star = p_data.get("star", {})
    samtools = p_data.get("samtools", {})
    repro = p_data.get("reproducibility", {})

    return ExecutionProfile(
        name=target_name,
        display_name=p_data.get("display_name", target_name),
        description=p_data.get("description", ""),
        physical_ram_gb_min=float(hw.get("physical_ram_gb", 16.0)),
        cpu_cores_min=int(hw.get("cpu_cores", 4)),
        storage_gb_min=float(hw.get("storage_gb", 100.0)),
        genome_sa_sparse_d=int(star.get("genome_sa_sparse_d", 1)),
        genome_sa_index_nbases=int(star.get("genome_sa_index_nbases", 14)),
        limit_genome_generate_ram=int(star.get("limit_genome_generate_ram", 31000000000)),
        bam_sort_mode=str(star.get("bam_sort_mode", "SortedByCoordinate")),
        limit_bam_sort_ram=int(star.get("limit_bam_sort_ram", 0)),
        threads=int(star.get("threads", 4)),
        external_sort_enabled=bool(samtools.get("external_sort_enabled", False)),
        samtools_sort_mem_per_thread=str(samtools.get("sort_mem_per_thread", "1G")),
        reproducibility_tier=str(repro.get("tier", "standard")),
    )


def _auto_select_profile(profiles_dict: dict[str, Any]) -> str:
    """Heuristically select a suitable profile based on detected physical RAM."""
    try:
        # Check Linux / WSL RAM
        if sys.platform != "win32":
            mem_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
            ram_gb = mem_bytes / (1024**3)
        else:
            res = run_wsl("grep MemTotal /proc/meminfo")
            if res.returncode == 0:
                kb = int(res.stdout.split()[1])
                ram_gb = kb / (1024**2)
            else:
                ram_gb = 16.0
    except Exception:
        ram_gb = 16.0

    if ram_gb >= 60.0 and "institutional_hpc" in profiles_dict:
        return "institutional_hpc"
    elif ram_gb >= 28.0 and "workstation_standard" in profiles_dict:
        return "workstation_standard"
    else:
        return "workstation_evaluation"


def resolve_pipeline_paths(repo_root: Path | None = None) -> dict[str, Path]:
    """Resolve portable pipeline storage directories without machine-specific hardcoding.

    Honors environment variable overrides:
    - DNT_DATA_DIR: root data directory
    - DNT_REF_DIR: reference genome and annotation directory
    """
    root = Path(repo_root or Path(__file__).resolve().parents[2]).resolve()

    data_dir_env = os.environ.get("DNT_DATA_DIR")
    data_dir = Path(data_dir_env).resolve() if data_dir_env else root / "data"

    ref_dir_env = os.environ.get("DNT_REF_DIR")
    ref_dir = Path(ref_dir_env).resolve() if ref_dir_env else root / "references" / "ensembl_112"

    return {
        "repo_root": root,
        "data_dir": data_dir,
        "raw_fastqs_dir": data_dir / "raw" / "prjna802072_fastqs",
        "trimmed_dir": data_dir / "interim" / "prjna802072_trimmed",
        "alignments_dir": data_dir / "interim" / "prjna802072_alignments",
        "qc_reports_dir": data_dir / "interim" / "prjna802072_qc_reports",
        "bridge_out_dir": data_dir / "interim" / "prjna802072_bridge",
        "config_dir": root / "config" / "quantification_bridge",
        "ref_dir": ref_dir,
        "fasta_path": ref_dir / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        "gtf_path": data_dir / "raw" / "references" / "Homo_sapiens.GRCh38.112.gtf.gz",
        "bed12_path": ref_dir / "Homo_sapiens.GRCh38.112.bed",
        "star_index_dir": ref_dir / "star_index",
    }


def probe_hardware_baseline(
    data_dir: Path | None = None,
    ref_dir: Path | None = None,
) -> HardwareBaseline:
    """Empirically inspect CPU, RAM, swap, and disk capacity."""
    paths = resolve_pipeline_paths()
    d_dir = data_dir or paths["data_dir"]
    r_dir = ref_dir or paths["ref_dir"]

    # Ensure probe targets exist for disk usage check
    d_dir.mkdir(parents=True, exist_ok=True)
    r_dir.mkdir(parents=True, exist_ok=True)

    cpu_cores = os.cpu_count() or 4
    cpu_threads = cpu_cores
    cpu_model = platform.processor() or "x86_64"
    architecture = platform.machine()
    is_wsl = False

    # Default fallbacks
    ram_gb = 16.0
    swap_gb = 0.0
    os_kernel = platform.platform()

    if sys.platform != "win32":
        try:
            mem_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
            ram_gb = mem_bytes / (1024**3)
        except Exception:
            pass

        try:
            with open("/proc/meminfo", "r") as f:
                lines = f.readlines()
                for l in lines:
                    if l.startswith("SwapTotal:"):
                        swap_gb = int(l.split()[1]) / (1024**2)
                    elif l.startswith("model name"):
                        cpu_model = l.split(":", 1)[1].strip()
        except Exception:
            pass
    else:
        # Windows host with WSL
        is_wsl = True
        res = run_wsl("grep MemTotal /proc/meminfo; echo '==='; grep SwapTotal /proc/meminfo; echo '==='; uname -r; echo '==='; grep 'model name' /proc/cpuinfo | head -n 1")
        if res.returncode == 0:
            parts = res.stdout.strip().split("===")
            if len(parts) >= 4:
                m_line = parts[0].strip()
                s_line = parts[1].strip()
                os_kernel = f"WSL2 (kernel {parts[2].strip()})"
                cpu_line = parts[3].strip()

                if "MemTotal:" in m_line:
                    ram_gb = int(m_line.split()[1]) / (1024**2)
                if "SwapTotal:" in s_line:
                    swap_gb = int(s_line.split()[1]) / (1024**2)
                if ":" in cpu_line:
                    cpu_model = cpu_line.split(":", 1)[1].strip()

    # Disk usage
    d_usage = shutil.disk_usage(d_dir)
    r_usage = shutil.disk_usage(r_dir)

    return HardwareBaseline(
        cpu_cores=cpu_cores,
        cpu_threads=cpu_threads,
        cpu_model=cpu_model,
        physical_ram_gb=ram_gb,
        swap_gb=swap_gb,
        data_dir_free_gb=d_usage.free / (1024**3),
        ref_dir_free_gb=r_usage.free / (1024**3),
        os_kernel=os_kernel,
        is_wsl=is_wsl,
        architecture=architecture,
    )


def probe_tool_inventory() -> ToolInventory:
    """Audit installed versions and paths of required bioinformatics tools."""
    tools = ["fastp", "STAR", "featureCounts", "samtools", "infer_experiment.py", "fasterq-dump"]
    inventory: dict[str, dict[str, str | None]] = {}

    batch_script = (
        "which fastp; which STAR; which featureCounts; which samtools; which infer_experiment.py; which fasterq-dump; "
        "echo '===DELIM==='; fastp --version 2>&1; "
        "echo '===DELIM==='; STAR --version; "
        "echo '===DELIM==='; featureCounts -v 2>&1 | grep featureCounts || true; "
        "echo '===DELIM==='; samtools --version | head -n 1; "
        "echo '===DELIM==='; infer_experiment.py --version 2>&1; "
        "echo '===DELIM==='; fasterq-dump --version 2>&1"
    )

    res = run_wsl(batch_script)
    if res.returncode == 0:
        sections = res.stdout.strip().split("===DELIM===")
        if len(sections) >= 7:
            paths = [l.strip() for l in sections[0].strip().splitlines() if l.strip()]
            for i, t in enumerate(tools):
                t_path = paths[i] if i < len(paths) else None
                t_ver = sections[i + 1].strip() if (i + 1) < len(sections) else "UNKNOWN"
                inventory[t] = {
                    "path": t_path if t_path and "not found" not in t_path.lower() else None,
                    "version": t_ver if t_ver else None,
                }
    else:
        for t in tools:
            inventory[t] = {"path": None, "version": None}

    return ToolInventory(tools=inventory)


def run_profile_preflight_checks(
    profile: ExecutionProfile,
    data_dir: Path | None = None,
    ref_dir: Path | None = None,
    enforce_physical_ram: bool = False,
    validate_pilot: bool = False,
    study_id: str | None = None,
) -> PreflightCheckResult:
    """Validate resource compatibility, tool availability, and pilot consistency against an execution profile."""
    hardware = probe_hardware_baseline(data_dir, ref_dir)
    tools = probe_tool_inventory()
    findings: list[PreflightFinding] = []

    # 0. Pilot Configuration & Manifest Consistency Check (PRJNA802072 Benchmark Only)
    should_validate_pilot = validate_pilot or (study_id == TARGET_BIOPROJECT)
    if should_validate_pilot:
        try:
            validate_pilot_configuration(study_id=study_id)
        except Exception as e:
            findings.append(
                PreflightFinding(
                    severity="FATAL",
                    code="pilot_configuration_invalid",
                    message=f"PRJNA802072 technical pilot configuration consistency failure: {e}",
                    remediation="Reconcile execution_profiles.yaml technical_pilot_specification with approved manifest.py.",
                )
            )

    # 1. Tool availability (FATAL if missing)
    for t_name, info in tools.tools.items():
        if not info.get("path"):
            findings.append(
                PreflightFinding(
                    severity="FATAL",
                    code=f"tool_missing_{t_name}",
                    message=f"Required bioinformatics tool '{t_name}' not found on PATH or in conda environment.",
                    remediation=f"Activate dnt_bridge environment or run 'conda env create -f config/quantification_bridge/environment_bridge.yml'.",
                )
            )

    # 2. Physical RAM & Process Memory Verification
    # Suffix Array and index footprint calculation
    if profile.genome_sa_sparse_d >= 4:
        est_star_ram = 9.98 if profile.genome_sa_index_nbases <= 12 else (10.39 if profile.genome_sa_index_nbases <= 13 else 12.0)
    elif profile.genome_sa_sparse_d >= 2:
        est_star_ram = 18.0
    else:
        est_star_ram = 34.0

    # Total estimated process RSS: index in memory + thread alignment stacks + I/O buffers
    thread_overhead = 0.12 * profile.threads
    est_process_rss = est_star_ram + thread_overhead + 0.25
    os_overhead = 0.50
    total_node_ram_needed = est_process_rss + os_overhead

    effective_ram = hardware.physical_ram_gb
    headroom_gb = effective_ram - total_node_ram_needed

    if effective_ram < profile.physical_ram_gb_min:
        msg = (
            f"Physical RAM ({effective_ram:.1f} GB) is below the minimum recommended for profile "
            f"'{profile.name}' ({profile.physical_ram_gb_min:.1f} GB)."
        )
        if enforce_physical_ram:
            findings.append(
                PreflightFinding(
                    severity="FATAL",
                    code="insufficient_physical_ram",
                    message=msg,
                    remediation="Select a lighter execution profile (e.g. workstation_evaluation) or allocate more RAM.",
                )
            )
        else:
            findings.append(
                PreflightFinding(
                    severity="WARNING",
                    code="ram_below_profile_minimum",
                    message=f"{msg} Execution may rely on swap, which can reduce alignment throughput.",
                    remediation="Monitor memory usage or consider running under workstation_evaluation.",
                )
            )

    # Process Memory Headroom Check
    # Note: Swap availability does NOT make low physical memory acceptable for production alignment.
    if effective_ram < est_process_rss:
        deficit = est_process_rss - effective_ram
        findings.append(
            PreflightFinding(
                severity="FATAL" if enforce_physical_ram else "WARNING",
                code="insufficient_ram_for_process",
                message=(
                    f"CRITICAL MEMORY CONSTRAINT: Estimated process memory requirement ({est_process_rss:.2f} GB) "
                    f"exceeds physical RAM ({effective_ram:.2f} GB) by {deficit:.2f} GB. Swap availability does NOT "
                    f"make low physical memory acceptable for production alignment; severe I/O thrashing will occur."
                ),
                remediation="Allocate additional physical RAM to node or run on institutional HPC server.",
            )
        )
    elif headroom_gb < 0.5:
        findings.append(
            PreflightFinding(
                severity="WARNING",
                code="narrow_memory_headroom",
                message=(
                    f"Physical RAM headroom is narrow ({headroom_gb:.2f} GB remaining beyond estimated STAR process RSS "
                    f"of {est_process_rss:.2f} GB and OS overhead). Close unneeded host applications."
                ),
                remediation="Close unneeded host applications and restrict STAR threads to 2.",
            )
        )

    # Check Total Virtual Memory (RAM + Swap)
    total_vm = hardware.physical_ram_gb + hardware.swap_gb
    if total_vm < est_process_rss + 4.0:
        findings.append(
            PreflightFinding(
                severity="FATAL",
                code="insufficient_total_virtual_memory",
                message=(
                    f"Total virtual memory ({total_vm:.1f} GB = {hardware.physical_ram_gb:.1f} GB RAM + "
                    f"{hardware.swap_gb:.1f} GB swap) is insufficient for STAR ({est_process_rss + 4.0:.1f} GB required)."
                ),
                remediation="Increase swap allocation in /etc/wsl.conf or .wslconfig to at least 32 GB.",
            )
        )

    # 3. Peak Storage Safeguards
    # 2-sample pilot: sequential peak 53.26 GB, concurrent peak 58.49 GB, safe free-space threshold: 70.0 GB.
    # Full 12-sample study: 270.0 GB.
    pilot_min_disk = 55.0  # Absolute minimum execution barrier
    pilot_safe_disk = 70.0 # Recommended safe free space
    full_min_disk = 270.0

    if hardware.data_dir_free_gb < pilot_min_disk:
        findings.append(
            PreflightFinding(
                severity="FATAL",
                code="insufficient_disk_space_pilot",
                message=(
                    f"Free storage on data directory ({hardware.data_dir_free_gb:.1f} GB) is below the absolute minimum "
                    f"required for 2-sample pilot execution ({pilot_min_disk:.1f} GB required; references + index + FASTQs + BAMs)."
                ),
                remediation="Free space on data disk or set DNT_DATA_DIR to a disk volume with at least 70 GB free.",
            )
        )
    elif hardware.data_dir_free_gb < pilot_safe_disk:
        findings.append(
            PreflightFinding(
                severity="WARNING",
                code="low_disk_space_pilot_margin",
                message=(
                    f"Free storage ({hardware.data_dir_free_gb:.1f} GB) accommodates the 2-sample pilot peak (53.3 GB) "
                    f"but is below the recommended safe operating margin ({pilot_safe_disk:.1f} GB)."
                ),
                remediation="Consider freeing additional space to ensure at least 70 GB free during alignment.",
            )
        )
    elif hardware.data_dir_free_gb < full_min_disk:
        findings.append(
            PreflightFinding(
                severity="WARNING",
                code="low_disk_space_full_study",
                message=(
                    f"Free storage ({hardware.data_dir_free_gb:.1f} GB) satisfies the 2-sample pilot ({pilot_safe_disk:.1f} GB), "
                    f"but is below the full 12-sample study requirement ({full_min_disk:.1f} GB)."
                ),
                remediation="Ensure additional storage before launching the complete 12-sample cohort.",
            )
        )

    # 4. CPU core availability
    if hardware.cpu_cores < profile.cpu_cores_min:
        findings.append(
            PreflightFinding(
                severity="WARNING",
                code="low_cpu_cores",
                message=f"Detected CPU cores ({hardware.cpu_cores}) is below profile minimum ({profile.cpu_cores_min}).",
                remediation="Execution will proceed with available threads, but runtime will be extended.",
            )
        )

    is_compatible = not any(f.severity == "FATAL" for f in findings)

    provenance = {
        "profile": profile.to_dict(),
        "hardware": hardware.to_dict(),
        "tools": tools.to_dict(),
        "findings_count": len(findings),
        "fatal_count": sum(1 for f in findings if f.severity == "FATAL"),
        "warning_count": sum(1 for f in findings if f.severity == "WARNING"),
    }

    return PreflightCheckResult(
        profile=profile,
        is_compatible=is_compatible,
        hardware=hardware,
        tools=tools,
        findings=tuple(findings),
        provenance=provenance,
    )


def evaluate_qc_metric(
    metric_name: str,
    value: float,
    config: dict[str, Any] | None = None,
) -> tuple[str, str]:
    """Evaluate an RNA-seq metric against documented evaluation criteria.

    Returns:
        (status, message) where status is "PASS", "WARNING", or "FATAL"
    """
    cfg = config or load_execution_profiles()
    qc_rules = cfg.get("qc_thresholds", {})
    warnings = qc_rules.get("warnings", {})
    fatal = qc_rules.get("fatal", {})

    # Check fatal bounds
    if metric_name == "mapped_reads" and value <= fatal.get("zero_mapped_reads", 0):
        return "FATAL", f"Zero mapped reads detected ({value}). Assay processing failed."
    if metric_name == "fastq_file_size" and value <= fatal.get("empty_fastq_bytes", 0):
        return "FATAL", f"Empty FASTQ input detected ({value} bytes)."

    # Check advisory warnings
    if metric_name == "fastp_q30_bases_fraction":
        threshold = warnings.get("fastp_min_q30_bases_fraction", 0.85)
        if value < threshold:
            return "WARNING", f"Q30 base fraction ({value:.3f}) is below documented criterion ({threshold:.3f})."
    elif metric_name == "fastp_read_survival_fraction":
        threshold = warnings.get("fastp_min_read_survival_fraction", 0.90)
        if value < threshold:
            return "WARNING", f"Read survival fraction ({value:.3f}) is below documented criterion ({threshold:.3f})."
    elif metric_name == "star_unique_mapping_fraction":
        threshold = warnings.get("star_min_unique_mapping_fraction", 0.75)
        if value < threshold:
            return "WARNING", f"STAR unique mapping fraction ({value:.3f}) is below target organoid benchmark ({threshold:.3f}). Flagged for investigator review."
    elif metric_name == "star_multimapping_fraction":
        threshold = warnings.get("star_max_multimapping_fraction", 0.15)
        if value > threshold:
            return "WARNING", f"STAR multimapping fraction ({value:.3f}) exceeds advisory threshold ({threshold:.3f})."
    elif metric_name == "featurecounts_assigned_fraction":
        threshold = warnings.get("featurecounts_min_assigned_fraction", 0.60)
        if value < threshold:
            return "WARNING", f"featureCounts assigned exonic read fraction ({value:.3f}) is below target benchmark ({threshold:.3f})."
    elif metric_name == "rseqc_strand_fraction":
        threshold = warnings.get("rseqc_min_strand_specificity_fraction", 0.85)
        if value < threshold and abs(value - 0.5) > 0.10:
            return "WARNING", f"Strand specificity fraction ({value:.3f}) is ambiguous (expected > {threshold} or unstranded ~0.5)."

    return "PASS", f"Metric '{metric_name}' ({value}) satisfies evaluation criteria."


def validate_pilot_configuration(
    config_path: Path | None = None,
    study_id: str | None = None,
) -> dict[str, Any]:
    """Strictly validate pilot configuration against authoritative manifest.py and pilot_manifest.json.

    NOTE: Restrictions below are strictly scoped to the PRJNA802072 benchmark cohort.
    If another valid study (e.g. E-MTAB-15191 or GSE232218) is analyzed, PRJNA802072
    pilot constraints do not restrict it.
    """
    if study_id is not None and study_id != TARGET_BIOPROJECT:
        return {
            "status": "NON_PRJNA802072_STUDY_ACCEPTED",
            "study_accession": study_id,
            "message": f"Study '{study_id}' is not PRJNA802072; benchmark pilot constraints do not apply.",
        }

    cfg = load_execution_profiles(config_path)
    pilot_spec = cfg.get("technical_pilot_specification")
    if not pilot_spec:
        raise ValueError("Missing 'technical_pilot_specification' in execution profiles config.")

    # Validate BioProject and Study
    if pilot_spec.get("study_accession") != TARGET_BIOPROJECT:
        raise ValueError(f"Study accession mismatch: expected {TARGET_BIOPROJECT}, got {pilot_spec.get('study_accession')}")
    if pilot_spec.get("developmental_day") != TARGET_DEVELOPMENTAL_DAY:
        raise ValueError(f"Developmental day mismatch: expected {TARGET_DEVELOPMENTAL_DAY}, got {pilot_spec.get('developmental_day')}")

    # Validate Control
    ctrl = pilot_spec.get("control", {})
    if ctrl.get("run_accession") != "SRR17834920":
        raise ValueError(f"Control run mismatch: expected SRR17834920, got {ctrl.get('run_accession')}")
    if ctrl.get("dose_um") != 0.0:
        raise ValueError(f"Control dose mismatch: expected 0.0 uM, got {ctrl.get('dose_um')}")
    if ctrl.get("developmental_day") != 28:
        raise ValueError(f"Control day mismatch: expected 28, got {ctrl.get('developmental_day')}")

    # Validate Treatment
    trt = pilot_spec.get("treatment", {})
    if trt.get("run_accession") != "SRR17834905":
        raise ValueError(f"Treatment run mismatch: expected SRR17834905, got {trt.get('run_accession')}")
    if trt.get("dose_um") != 250.0:
        raise ValueError(f"Treatment dose mismatch: expected 250.0 uM, got {trt.get('dose_um')}")
    if trt.get("developmental_day") != 28:
        raise ValueError(f"Treatment day mismatch: expected 28, got {trt.get('developmental_day')}")

    # Prohibited accessions & doses check
    prohibited_runs = set(pilot_spec.get("prohibited_accessions", [])) | {"SRR17826941", "SRR17826944"}
    if ctrl.get("run_accession") in prohibited_runs or trt.get("run_accession") in prohibited_runs:
        raise ValueError("Prohibited run accession detected in pilot specification!")

    prohibited_doses = set(pilot_spec.get("prohibited_concentrations", [])) | {1000.0, 50000.0}
    if ctrl.get("dose_um") in prohibited_doses or trt.get("dose_um") in prohibited_doses:
        raise ValueError("Prohibited dose detected in pilot specification (e.g. 1 mM / 50 mM from other studies)!")

    # Verify against manifest.py records
    manifest_ctrl = RUN_TO_SAMPLE_MAP.get("SRR17834920")
    manifest_trt = RUN_TO_SAMPLE_MAP.get("SRR17834905")
    if not manifest_ctrl or not manifest_trt:
        raise ValueError("Approved pilot runs missing from RUN_TO_SAMPLE_MAP in manifest.py")

    return {
        "status": "PILOT_CONFIGURATION_VALID",
        "study_accession": TARGET_BIOPROJECT,
        "developmental_day": TARGET_DEVELOPMENTAL_DAY,
        "control_run": manifest_ctrl.run_accession,
        "control_dose_um": manifest_ctrl.vpa_dose_um,
        "treatment_run": manifest_trt.run_accession,
        "treatment_dose_um": manifest_trt.vpa_dose_um,
    }
