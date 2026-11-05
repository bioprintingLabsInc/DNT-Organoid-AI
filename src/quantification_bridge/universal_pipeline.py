"""Universal FASTQ Automation Pipeline — Phase F1.

Orchestrates reproducible, generalized RNA-seq preprocessing and quantification
across researcher-provided datasets, without Google Colab dependencies:

1. Flexible Inputs:
   - Local FASTQs (raw or gzip-compressed)
   - Paired-end and single-end libraries
   - SRA run accessions via SRA Toolkit / fasterq-dump
   - Sample metadata in CSV, TSV, JSON, or YAML format
2. Generalized Workflow:
   - fastp quality control and trimming
   - STAR genome index generation and coordinate-sorted alignment
   - RSeQC empirical strandedness diagnosis and evidence serialization
   - featureCounts quantification with empirical strandedness propagation
   - Downstream integration: Universal Ingestion v2 metadata validation,
     matrix QC, gene harmonization, normalization, and inferential gating.
3. Preflight & Resource Validation:
   - Input existence, non-emptiness, and read pairing
   - Duplicate sample ID detection
   - Reference genome and annotation compatibility
   - Software tools and version probing
   - Host CPU, physical RAM, and disk space auditing against execution profiles
   - Dry-run mode with planned stage scheduling and explicit approval gate
4. Resumable Checkpoints:
   - Fine-grained stage checkpoints with cryptographic input/config verification
   - Safe resumption of interrupted runs
   - Prevents duplicate reprocessing while blocking incompatible reuse
5. Scientific Safeguards:
   - Uninvented biological replicates (missing stays unknown)
   - Unverified replication strictly gates inferential GLM and pathway enrichment
   - Evidence-backed strandedness decisions
   - Independent GO and KEGG execution and QC tracking
"""

from dataclasses import asdict, dataclass, field
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import time
from typing import Any
import yaml

from src.metadata.integrity import MetadataIntegrityChecker
from src.metadata.mapper import CanonicalMetadataMapper, RawSampleRecord, StudyContext
from src.metadata.result import ValidationStatus
from src.metadata.validator import MetadataValidator

from .adapter import convert_featurecounts_to_matrix, sha256_file
from .bed_generator import generate_bed12_from_gtf
from .downstream_bridge import DownstreamBridgeResult, execute_downstream_pipeline
from .env_exec import run_wsl, to_wsl_path
from .fastp_runner import run_fastp_preprocessing
from .featurecounts_runner import parse_featurecounts_summary, run_featurecounts
from .pairing import FastqPair, FastqSingle, determine_read_mate, extract_sample_stem
from .profiles import (
    ExecutionProfile,
    HardwareBaseline,
    ToolInventory,
    get_execution_profile,
    run_profile_preflight_checks,
)
from .rseqc_runner import StrandednessDecision, infer_strandedness_decision, run_infer_experiment
from .sra_runner import extract_sra, probe_sra_tools
from .star_runner import align_reads_star, build_star_index, parse_star_final_log
from .validator import validate_paired_fastqs, validate_single_fastq


# -----------------------------------------------------------------------------
# Data Models
# -----------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class SampleInput:
    """Record describing a single sequencing library and its metadata."""

    sample_id: str
    r1_path: Path | None = None
    r2_path: Path | None = None
    sra_accession: str | None = None
    read_layout: str = "paired_end"  # "paired_end" or "single_end"
    attributes: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "r1_path": str(self.r1_path) if self.r1_path else None,
            "r2_path": str(self.r2_path) if self.r2_path else None,
            "sra_accession": self.sra_accession,
            "read_layout": self.read_layout,
            "attributes": self.attributes,
        }


@dataclass(frozen=True, slots=True)
class ReferenceConfig:
    """Configured reference genome and transcriptomic annotations."""

    species: str
    genome_fasta: Path
    gtf_annotation: Path
    star_index_dir: Path | None = None
    bed12_annotation: Path | None = None
    sjdb_overhang: int = 149

    def to_dict(self) -> dict[str, Any]:
        return {
            "species": self.species,
            "genome_fasta": str(self.genome_fasta),
            "gtf_annotation": str(self.gtf_annotation),
            "star_index_dir": str(self.star_index_dir) if self.star_index_dir else None,
            "bed12_annotation": str(self.bed12_annotation) if self.bed12_annotation else None,
            "sjdb_overhang": self.sjdb_overhang,
        }


@dataclass(frozen=True, slots=True)
class UniversalPipelineConfig:
    """Operational settings for universal FASTQ processing."""

    study_id: str
    output_dir: Path | None = None
    work_dir: Path | None = None
    study_title: str | None = None
    source_type: str = "in_house"
    experiment_id: str = "exp_01"
    profile_name: str | None = "workstation_evaluation"
    read_layout: str = "paired_end"  # "paired_end" or "single_end"
    strandedness: int | None = None  # None for auto RSeQC inference, or 0, 1, 2
    threads: int = 4
    require_approval: bool = True
    enforce_resources: bool = False
    min_biological_replicates: int = 2

    def __post_init__(self) -> None:
        target = self.output_dir or self.work_dir
        if target is None:
            raise ValueError("Either output_dir or work_dir must be specified for UniversalPipelineConfig.")
        object.__setattr__(self, "output_dir", Path(target).resolve())
        object.__setattr__(self, "work_dir", Path(target).resolve())

    def to_dict(self) -> dict[str, Any]:
        return {
            "study_id": self.study_id,
            "output_dir": str(self.output_dir),
            "work_dir": str(self.work_dir),
            "study_title": self.study_title,
            "source_type": self.source_type,
            "experiment_id": self.experiment_id,
            "profile_name": self.profile_name,
            "read_layout": self.read_layout,
            "strandedness": self.strandedness,
            "threads": self.threads,
            "require_approval": self.require_approval,
            "enforce_resources": self.enforce_resources,
            "min_biological_replicates": self.min_biological_replicates,
        }


@dataclass(frozen=True, slots=True)
class PipelineFinding:
    """Structured message recording validation findings and preflight checks."""

    severity: str  # "FATAL", "WARNING", "INFO"
    code: str
    message: str
    remediation: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class PreflightReport:
    """Comprehensive assessment of inputs, software, references, and hardware."""

    is_valid: bool
    findings: tuple[PipelineFinding, ...]
    tool_inventory: dict[str, Any]
    hardware_baseline: dict[str, Any]
    sample_validation: dict[str, Any]
    reference_validation: dict[str, Any]
    metadata_validation: dict[str, Any]
    reusable_checkpoints: tuple[str, ...]
    estimated_disk_gb: float
    estimated_ram_gb: float

    @property
    def has_fatal_errors(self) -> bool:
        return any(f.severity == "FATAL" for f in self.findings)

    def to_dict(self) -> dict[str, Any]:
        return {
            "is_valid": self.is_valid,
            "has_fatal_errors": self.has_fatal_errors,
            "findings": [f.to_dict() for f in self.findings],
            "tool_inventory": self.tool_inventory,
            "hardware_baseline": self.hardware_baseline,
            "sample_validation": self.sample_validation,
            "reference_validation": self.reference_validation,
            "metadata_validation": self.metadata_validation,
            "reusable_checkpoints": list(self.reusable_checkpoints),
            "estimated_disk_gb": round(self.estimated_disk_gb, 2),
            "estimated_ram_gb": round(self.estimated_ram_gb, 2),
        }


@dataclass(frozen=True, slots=True)
class ProcessingPlan:
    """Dry-run schedule specifying all planned operations before execution."""

    config: UniversalPipelineConfig
    samples: tuple[SampleInput, ...]
    references: ReferenceConfig
    preflight: PreflightReport
    planned_stages: tuple[dict[str, Any], ...]
    can_execute: bool
    approval_required: bool
    approval_token: str | None

    def to_dict(self) -> dict[str, Any]:
        return {
            "config": self.config.to_dict(),
            "samples": [s.to_dict() for s in self.samples],
            "references": self.references.to_dict(),
            "preflight": self.preflight.to_dict(),
            "planned_stages": list(self.planned_stages),
            "can_execute": self.can_execute,
            "approval_required": self.approval_required,
            "approval_token": self.approval_token,
        }


@dataclass(frozen=True, slots=True)
class UniversalFastqResult:
    """Outcome artifact containing counts, evidence, metrics, and downstream evaluation."""

    status: str
    study_id: str
    sample_ids: tuple[str, ...]
    raw_counts_path: Path | None
    canonical_metadata_path: Path | None
    qc_reports_dir: Path | None
    alignments_dir: Path | None
    strandedness_evidence_path: Path | None
    counting_summary_path: Path | None
    execution_record_path: Path | None
    downstream_result: DownstreamBridgeResult | None = None
    metrics: dict[str, Any] = field(default_factory=dict)
    findings: tuple[PipelineFinding, ...] = field(default_factory=tuple)
    completed_stages: tuple[str, ...] = field(default_factory=tuple)
    failed_stage: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "study_id": self.study_id,
            "sample_ids": list(self.sample_ids),
            "raw_counts_path": str(self.raw_counts_path) if self.raw_counts_path else None,
            "canonical_metadata_path": str(self.canonical_metadata_path) if self.canonical_metadata_path else None,
            "qc_reports_dir": str(self.qc_reports_dir) if self.qc_reports_dir else None,
            "alignments_dir": str(self.alignments_dir) if self.alignments_dir else None,
            "strandedness_evidence_path": str(self.strandedness_evidence_path) if self.strandedness_evidence_path else None,
            "counting_summary_path": str(self.counting_summary_path) if self.counting_summary_path else None,
            "execution_record_path": str(self.execution_record_path) if self.execution_record_path else None,
            "downstream_result": self.downstream_result.to_dict() if self.downstream_result else None,
            "metrics": self.metrics,
            "findings": [f.to_dict() for f in self.findings],
            "completed_stages": list(self.completed_stages),
            "failed_stage": self.failed_stage,
        }


# -----------------------------------------------------------------------------
# Metadata Ingestion Helper
# -----------------------------------------------------------------------------

def load_sample_metadata(
    metadata_source: Path | str | dict[str, Any] | list[dict[str, Any]],
    study_id: str = "study_01",
    source_type: str = "in_house",
    study_title: str | None = None,
) -> tuple[dict[str, Any], list[RawSampleRecord]]:
    """Parse researcher metadata from CSV, TSV, JSON, YAML or dict into canonical structure.

    Guarantees:
    - Preserves all source attributes without loss.
    - Missing attributes remain explicitly 'unknown' (never falsely invented).
    - Biological replicate identities are only recorded when provided.
    - Produces canonical metadata document validated against config/metadata_schema.yaml.
    """
    raw_records: list[RawSampleRecord] = []
    source_filename = "in_memory"

    # 1. Parse raw data into list of dictionaries
    if isinstance(metadata_source, (str, Path)):
        p = Path(metadata_source).resolve()
        if not p.is_file():
            raise FileNotFoundError(f"Metadata file not found: {p}")
        source_filename = p.name
        content = p.read_text(encoding="utf-8")
        ext = p.suffix.lower()

        if ext == ".csv":
            import csv
            import io
            reader = csv.DictReader(io.StringIO(content))
            raw_dicts = list(reader)
        elif ext in (".tsv", ".tab"):
            import csv
            import io
            reader = csv.DictReader(io.StringIO(content), delimiter="\t")
            raw_dicts = list(reader)
        elif ext == ".json":
            parsed = json.loads(content)
            if isinstance(parsed, dict) and "studies" in parsed and "samples" in parsed:
                # Pre-formed canonical metadata document
                return parsed, []
            elif isinstance(parsed, list):
                raw_dicts = parsed
            elif isinstance(parsed, dict) and "samples" in parsed and isinstance(parsed["samples"], list):
                raw_dicts = parsed["samples"]
            else:
                raw_dicts = [parsed]
        elif ext in (".yaml", ".yml"):
            parsed = yaml.safe_load(content)
            if isinstance(parsed, dict) and "studies" in parsed and "samples" in parsed:
                # Pre-formed canonical metadata document
                return parsed, []
            elif isinstance(parsed, list):
                raw_dicts = parsed
            elif isinstance(parsed, dict) and "samples" in parsed and isinstance(parsed["samples"], list):
                raw_dicts = parsed["samples"]
            else:
                raw_dicts = [parsed]
        else:
            raise ValueError(f"Unsupported metadata file format: {ext} (supported: .csv, .tsv, .json, .yaml, .yml)")
    elif isinstance(metadata_source, dict):
        if "studies" in metadata_source and "samples" in metadata_source:
            return metadata_source, []
        if "samples" in metadata_source and isinstance(metadata_source["samples"], list):
            raw_dicts = metadata_source["samples"]
        else:
            raw_dicts = [metadata_source]
    elif isinstance(metadata_source, list):
        raw_dicts = metadata_source
    else:
        raise TypeError(f"Invalid metadata_source type: {type(metadata_source)}")

    # 2. Extract RawSampleRecord instances
    id_keys = ("sample_id", "sample_name", "sample", "run_accession", "run", "accession", "sample_alias")
    for idx, row in enumerate(raw_dicts, start=1):
        if not isinstance(row, dict):
            continue
        sid = None
        for k in id_keys:
            if k in row and row[k]:
                sid = str(row[k]).strip()
                break
        if not sid:
            sid = f"sample_{idx:03d}"

        # Clean attributes
        clean_attrs = {str(k).strip(): v for k, v in row.items()}
        raw_records.append(
            RawSampleRecord(
                sample_id=sid,
                attributes=clean_attrs,
                source_file=source_filename,
                source_locator=f"row_{idx}",
                source_key_provenance={k: k for k in clean_attrs.keys()},
            )
        )

    # 3. Map into canonical metadata document using CanonicalMetadataMapper
    valid_source = "public" if str(source_type).lower() in ("public", "sra", "geo", "arrayexpress") else "in_house"
    mapper = CanonicalMetadataMapper()
    context = StudyContext(
        study_id=study_id,
        source_type=valid_source,
        study_title=study_title or f"Study {study_id}",
    )
    canonical_doc, mapping_findings = mapper.map_study(
        context=context,
        raw_records=raw_records,
        source_reference=source_filename,
    )

    return canonical_doc, raw_records


# -----------------------------------------------------------------------------
# Checkpoint Governance & Output Integrity
# -----------------------------------------------------------------------------

def compute_file_content_hash(path: Path | str) -> str:
    """Compute deterministic content-based hash using file content digest."""
    p = Path(path).resolve()
    if not p.exists():
        return f"missing:{p.name}"
    if p.is_dir():
        hasher = hashlib.sha256()
        for item in sorted(p.rglob("*")):
            if item.is_file():
                st = item.stat()
                hasher.update(f"{item.name}:{st.st_size}".encode("utf-8"))
                try:
                    with item.open("rb") as f:
                        hasher.update(f.read(65536))
                except Exception:
                    pass
        return hasher.hexdigest()
    st = p.stat()
    hasher = hashlib.sha256()
    hasher.update(f"{p.name}:{st.st_size}".encode("utf-8"))
    try:
        with p.open("rb") as f:
            bytes_read = 0
            while bytes_read < 10_000_000:
                chunk = f.read(65536)
                if not chunk:
                    break
                hasher.update(chunk)
                bytes_read += len(chunk)
            if st.st_size > 10_000_000:
                f.seek(-65536, os.SEEK_END)
                hasher.update(f.read(65536))
    except Exception:
        pass
    return hasher.hexdigest()


def _verify_output_integrity(out: Path | str) -> bool:
    """Verify physical output existence, non-emptiness, and format validity."""
    out_p = Path(out).resolve()
    if not out_p.exists():
        return False

    if out_p.is_dir():
        # Directory output (e.g. STAR index): check essential files exist and non-empty
        genome_f = out_p / "Genome"
        sa_f = out_p / "SA"
        if genome_f.exists() or sa_f.exists() or "star" in out_p.name.lower():
            if not (genome_f.is_file() and genome_f.stat().st_size > 0):
                return False
            if not (sa_f.is_file() and sa_f.stat().st_size > 0):
                return False
            return True
        files = [f for f in out_p.iterdir() if f.is_file() and f.stat().st_size > 0]
        return len(files) > 0

    if out_p.stat().st_size == 0:
        return False

    name_lower = out_p.name.lower()
    if name_lower.endswith(".json"):
        try:
            json.loads(out_p.read_text(encoding="utf-8"))
        except Exception:
            return False
    elif name_lower.endswith(".fastq.gz") or name_lower.endswith(".fq.gz"):
        try:
            with gzip.open(out_p, "rb") as gz:
                head = gz.read(256)
                if not head or not head.startswith(b"@"):
                    return False
        except Exception:
            return False
    elif name_lower.endswith(".gz"):
        try:
            with gzip.open(out_p, "rb") as gz:
                head = gz.read(64)
                if not head:
                    return False
        except Exception:
            return False
    elif name_lower.endswith(".fastq") or name_lower.endswith(".fq"):
        try:
            with out_p.open("rb") as f:
                head = f.read(256)
                if not head or not head.startswith(b"@"):
                    return False
        except Exception:
            return False
    elif name_lower.endswith(".bam"):
        try:
            with out_p.open("rb") as bf:
                head = bf.read(64)
            if head.startswith(b"BAM\x01"):
                return True
            if head.startswith(b"\x1f\x8b"):
                with gzip.open(out_p, "rb") as gz:
                    decomp_head = gz.read(64)
                    if decomp_head.startswith(b"BAM\x01") or b"bam" in decomp_head.lower():
                        return True
                    return False
            elif b"bam" in head.lower():
                return True
            return False
        except Exception:
            return False

    return True


class CheckpointManager:
    """Atomic file-based checkpoint store ensuring safe resumption and preventing invalid reuse."""

    def __init__(self, checkpoint_dir: Path | str) -> None:
        self.checkpoint_dir = Path(checkpoint_dir).resolve()
        # Defer directory creation to maintain filesystem immutability during dry-run / planning

    def _get_path(self, stage: str, identifier: str) -> Path:
        clean_id = re.sub(r"[^A-Za-z0-9_.-]", "_", identifier)
        return self.checkpoint_dir / f"{stage}_{clean_id}.json"

    def is_stage_completed(
        self,
        stage: str,
        identifier: str,
        expected_outputs: list[Path | str],
        input_hash: str,
        config_hash: str,
    ) -> bool:
        """Verify whether a previously executed stage can be safely reused without reprocessing."""
        if not self.checkpoint_dir.is_dir():
            return False

        chk_p = self._get_path(stage, identifier)
        if not chk_p.is_file():
            return False

        try:
            data = json.loads(chk_p.read_text(encoding="utf-8"))
        except Exception:
            return False

        if data.get("status") != "COMPLETED":
            return False
        if data.get("input_hash") != input_hash:
            return False
        if data.get("config_hash") != config_hash:
            return False

        # Verify all output files exist, are non-empty, intact, and unaltered
        recorded_hashes = data.get("output_hashes", {})
        for out in expected_outputs:
            out_p = Path(out).resolve()
            if not _verify_output_integrity(out_p):
                return False
            if str(out_p) in recorded_hashes:
                if compute_file_content_hash(out_p) != recorded_hashes[str(out_p)]:
                    return False

        return True

    def record_checkpoint(
        self,
        stage: str,
        identifier: str,
        outputs: list[Path | str],
        input_hash: str,
        config_hash: str,
        metrics: dict[str, Any] | None = None,
    ) -> None:
        """Write an atomic checkpoint record to disk."""
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        chk_p = self._get_path(stage, identifier)
        output_hashes: dict[str, str] = {}
        for o in outputs:
            p = Path(o).resolve()
            if p.exists():
                output_hashes[str(p)] = compute_file_content_hash(p)
        record = {
            "stage": stage,
            "identifier": identifier,
            "status": "COMPLETED",
            "input_hash": input_hash,
            "config_hash": config_hash,
            "outputs": [str(Path(o).resolve()) for o in outputs],
            "output_hashes": output_hashes,
            "metrics": metrics or {},
            "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        }
        chk_p.write_text(json.dumps(record, indent=2), encoding="utf-8")


PipelineCheckpointManager = CheckpointManager

_PROBED_TOOLS_CACHE: dict[str, dict[str, Any]] | None = None


# -----------------------------------------------------------------------------
# Master Universal Pipeline Orchestrator
# -----------------------------------------------------------------------------

class UniversalFastqPipeline:
    """Master orchestrator for universal FASTQ automation (Phase F1)."""

    def __init__(
        self,
        config: UniversalPipelineConfig,
        references: ReferenceConfig,
        samples: list[SampleInput],
        metadata_source: Path | str | dict[str, Any] | list[dict[str, Any]] | None = None,
        canonical_metadata: Path | str | dict[str, Any] | list[dict[str, Any]] | None = None,
    ) -> None:
        self.config = config
        self.references = references
        self.samples = samples
        prof_name = config.profile_name
        if prof_name in ("laptop", "evaluation"):
            prof_name = "workstation_evaluation"
        elif prof_name in ("standard", "workstation"):
            prof_name = "workstation_standard"
        self.profile: ExecutionProfile = get_execution_profile(prof_name)

        # Standardized output paths
        self.out_dir = Path(config.output_dir).resolve()
        self.qc_dir = self.out_dir / "fastp_qc"
        self.trimmed_dir = self.out_dir / "trimmed_fastqs"
        self.align_dir = self.out_dir / "alignments"
        self.counts_dir = self.out_dir / "counts"
        self.checkpoints_dir = self.out_dir / ".checkpoints"
        self.sra_dir = self.out_dir / "sra_extracted"
        self._approved_token: str | None = None

        # Defer directory creation to maintain filesystem immutability during dry-run / planning
        self.checkpoint_manager = CheckpointManager(self.checkpoints_dir)

        # Resolve or map metadata
        meta_src = metadata_source if metadata_source is not None else canonical_metadata
        if meta_src is not None:
            self.canonical_metadata, self.raw_sample_records = load_sample_metadata(
                metadata_source=meta_src,
                study_id=config.study_id,
                source_type=config.source_type,
                study_title=config.study_title,
            )
        else:
            # Build minimal canonical metadata from sample records without inventing defaults
            raw_dicts = []
            for s in samples:
                d = {"sample_id": s.sample_id}
                d.update(s.attributes)
                # Keep unprovided attributes explicitly missing or unknown - DO NOT invent vehicle, chemical, concentration, or unit
                raw_dicts.append(d)

            self.canonical_metadata, self.raw_sample_records = load_sample_metadata(
                metadata_source=raw_dicts,
                study_id=config.study_id,
                source_type=config.source_type,
                study_title=config.study_title,
            )

    def _ensure_directories(self) -> None:
        """Create output directories only upon approved execution initiation."""
        for d in (self.out_dir, self.qc_dir, self.trimmed_dir, self.align_dir, self.counts_dir, self.checkpoints_dir, self.sra_dir):
            d.mkdir(parents=True, exist_ok=True)

    def approve(self, token: str) -> None:
        """Explicitly record approval token for the planned execution."""
        self._approved_token = token

    # -------------------------------------------------------------------------
    # Preflight Verification & Dry-Run Planning
    # -------------------------------------------------------------------------

    def validate_inputs(self) -> PreflightReport:
        """Perform comprehensive pre-execution validation of software, files, metadata, and hardware."""
        findings: list[PipelineFinding] = []

        # 1. Hardware and software tool probing
        hw = self._probe_hardware()
        tools = self._probe_tools()

        # Check required software tools
        required_tools = ["fastp", "STAR", "featureCounts", "samtools", "infer_experiment.py"]
        has_sra = any(s.sra_accession and not s.r1_path for s in self.samples)
        if has_sra:
            required_tools.append("fasterq-dump")

        for t in required_tools:
            t_info = tools.get(t, {})
            if not t_info.get("available"):
                findings.append(
                    PipelineFinding(
                        severity="FATAL",
                        code=f"TOOL_MISSING_{t.upper()}",
                        message=f"Required bioinformatics utility '{t}' is not available in PATH or WSL.",
                        remediation=f"Install '{t}' in the active environment or Conda prefix.",
                    )
                )

        # Check hardware limits if enforced
        if self.config.enforce_resources:
            if hw.physical_ram_gb < self.profile.physical_ram_gb_min:
                findings.append(
                    PipelineFinding(
                        severity="FATAL",
                        code="INSUFFICIENT_RAM",
                        message=(
                            f"Host physical RAM ({hw.physical_ram_gb:.1f} GB) is below profile "
                            f"'{self.profile.name}' requirement ({self.profile.physical_ram_gb_min:.1f} GB)."
                        ),
                        remediation="Select a lighter execution profile (e.g. 'laptop' with sparseD index) or add RAM.",
                    )
                )

        # Estimate disk requirements: ~10 GB per sample
        est_disk_gb = max(15.0, len(self.samples) * 10.0)
        free_disk_gb = hw.data_dir_free_gb
        if free_disk_gb < est_disk_gb and self.config.enforce_resources:
            findings.append(
                PipelineFinding(
                    severity="FATAL",
                    code="INSUFFICIENT_DISK",
                    message=f"Available disk space ({free_disk_gb:.1f} GB) is below estimated requirement ({est_disk_gb:.1f} GB).",
                    remediation="Free space on the target volume or select an alternative storage location.",
                )
            )

        # 2. Reference validation
        ref_checks: dict[str, Any] = {
            "species": self.references.species,
            "fasta_exists": self.references.genome_fasta.is_file(),
            "gtf_exists": self.references.gtf_annotation.is_file(),
        }

        if not self.references.species or not self.references.species.strip():
            findings.append(
                PipelineFinding(
                    severity="FATAL",
                    code="SPECIES_UNSPECIFIED",
                    message="Species is not specified in ReferenceConfig. Species must be explicitly defined.",
                    remediation="Specify target species (e.g. 'Homo sapiens' or 'Mus musculus').",
                )
            )

        if not self.references.genome_fasta.is_file():
            findings.append(
                PipelineFinding(
                    severity="FATAL",
                    code="REFERENCE_FASTA_MISSING",
                    message=f"Genome FASTA file not found at: {self.references.genome_fasta}",
                    remediation="Provide a valid, accessible primary assembly FASTA reference file.",
                )
            )
        elif self.references.genome_fasta.stat().st_size == 0:
            findings.append(
                PipelineFinding(
                    severity="FATAL",
                    code="REFERENCE_FASTA_EMPTY",
                    message=f"Genome FASTA file is 0 bytes: {self.references.genome_fasta}",
                    remediation="Provide an intact, non-empty genome FASTA reference.",
                )
            )

        if not self.references.gtf_annotation.is_file():
            findings.append(
                PipelineFinding(
                    severity="FATAL",
                    code="REFERENCE_GTF_MISSING",
                    message=f"Gene annotation GTF file not found at: {self.references.gtf_annotation}",
                    remediation="Provide a valid, accessible GTF annotation file.",
                )
            )
        elif self.references.gtf_annotation.stat().st_size == 0:
            findings.append(
                PipelineFinding(
                    severity="FATAL",
                    code="REFERENCE_GTF_EMPTY",
                    message=f"Gene annotation GTF file is 0 bytes: {self.references.gtf_annotation}",
                    remediation="Provide an intact, non-empty GTF annotation reference.",
                )
            )

        # 3. Sample inputs & FASTQ validation
        sample_checks: dict[str, Any] = {}
        seen_sids: set[str] = set()

        if not self.samples:
            findings.append(
                PipelineFinding(
                    severity="FATAL",
                    code="NO_SAMPLES_PROVIDED",
                    message="No sequencing samples provided to UniversalFastqPipeline.",
                    remediation="Supply at least one valid SampleInput to execute or dry-run.",
                )
            )

        # Check layout consistency across samples in the cohort
        layouts = {s.read_layout for s in self.samples if s.read_layout}
        if len(layouts) > 1:
            findings.append(
                PipelineFinding(
                    severity="FATAL",
                    code="MIXED_SEQUENCING_LAYOUTS_UNSUPPORTED",
                    message=f"Cohort contains mixed sequencing layouts: {sorted(layouts)}. Mixed single-end and paired-end runs within a single study pipeline run are not supported.",
                    remediation="Split the cohort into separate single-end and paired-end pipeline runs.",
                )
            )

        for s in self.samples:
            # Duplicate sample ID check
            if s.sample_id in seen_sids:
                findings.append(
                    PipelineFinding(
                        severity="FATAL",
                        code="DUPLICATE_SAMPLE_ID",
                        message=f"Duplicate sample identifier detected: '{s.sample_id}'. Each sample must be unique.",
                        remediation="Ensure all sample identifiers in the cohort are unique.",
                    )
                )
            seen_sids.add(s.sample_id)

            # Sample layout must match pipeline configured layout
            if s.read_layout != self.config.read_layout:
                findings.append(
                    PipelineFinding(
                        severity="FATAL",
                        code="SAMPLE_LAYOUT_MISMATCH",
                        message=f"Sample '{s.sample_id}' layout '{s.read_layout}' does not match pipeline configured layout '{self.config.read_layout}'.",
                        remediation=f"Configure sample layout to match study layout '{self.config.read_layout}'.",
                    )
                )

            # Check read pairing and files
            if s.read_layout == "paired_end":
                if s.r1_path is None and s.r2_path is None and not s.sra_accession:
                    findings.append(
                        PipelineFinding(
                            severity="FATAL",
                            code="MISSING_FASTQ_INPUTS",
                            message=f"Paired-end sample '{s.sample_id}' has neither FASTQ files nor SRA accession.",
                            remediation="Provide paired FASTQ paths (r1_path, r2_path) or sra_accession.",
                        )
                    )
                elif s.r1_path is not None and s.r2_path is None:
                    findings.append(
                        PipelineFinding(
                            severity="FATAL",
                            code="MISSING_FASTQ_MATE",
                            message=f"Sample '{s.sample_id}' configured as paired_end but reverse mate (R2) is missing.",
                            remediation="Provide R2 mate path or reclassify read_layout as 'single_end'.",
                        )
                    )
                elif s.r1_path is None and s.r2_path is not None:
                    findings.append(
                        PipelineFinding(
                            severity="FATAL",
                            code="MISSING_FASTQ_MATE",
                            message=f"Sample '{s.sample_id}' configured as paired_end but forward mate (R1) is missing.",
                            remediation="Provide R1 mate path or reclassify read_layout as 'single_end'.",
                        )
                    )
                elif s.sra_accession and not s.r1_path:
                    sample_checks[s.sample_id] = {
                        "status": "REQUIRES_SRA_RETRIEVAL",
                        "accession": s.sra_accession,
                        "read_layout": "paired_end",
                    }
                elif s.r1_path and s.r2_path:
                    # Validate actual FASTQ files
                    if not s.r1_path.is_file():
                        findings.append(
                            PipelineFinding(
                                severity="FATAL",
                                code="FASTQ_R1_NOT_FOUND",
                                message=f"Read 1 FASTQ not found for '{s.sample_id}': {s.r1_path}",
                            )
                        )
                    if not s.r2_path.is_file():
                        findings.append(
                            PipelineFinding(
                                severity="FATAL",
                                code="FASTQ_R2_NOT_FOUND",
                                message=f"Read 2 FASTQ not found for '{s.sample_id}': {s.r2_path}",
                            )
                        )
                    if s.r1_path.is_file() and s.r2_path.is_file():
                        try:
                            v_res = validate_paired_fastqs(s.r1_path, s.r2_path, sample_records_to_check=5000)
                            sample_checks[s.sample_id] = v_res
                        except Exception as e:
                            findings.append(
                                PipelineFinding(
                                    severity="FATAL",
                                    code="FASTQ_VALIDATION_ERROR",
                                    message=f"Structural FASTQ validation failed for '{s.sample_id}': {e}",
                                )
                            )
            elif s.read_layout == "single_end":
                if s.r2_path is not None:
                    findings.append(
                        PipelineFinding(
                            severity="FATAL",
                            code="SAMPLE_LAYOUT_MISMATCH",
                            message=f"Single-end sample '{s.sample_id}' has mate R2 path '{s.r2_path}'. Single-end samples must not have an R2 mate.",
                            remediation="Remove R2 path or reconfigure layout to 'paired_end'.",
                        )
                    )
                if s.r1_path is None:
                    if not s.sra_accession:
                        findings.append(
                            PipelineFinding(
                                severity="FATAL",
                                code="FASTQ_R1_MISSING",
                                message=f"Single-end sample '{s.sample_id}' has neither FASTQ path nor SRA accession.",
                            )
                        )
                    else:
                        sample_checks[s.sample_id] = {
                            "status": "REQUIRES_SRA_RETRIEVAL",
                            "accession": s.sra_accession,
                            "read_layout": "single_end",
                        }
                elif not s.r1_path.is_file():
                    findings.append(
                        PipelineFinding(
                            severity="FATAL",
                            code="FASTQ_R1_NOT_FOUND",
                            message=f"Single-end FASTQ not found for '{s.sample_id}': {s.r1_path}",
                        )
                    )
                else:
                    try:
                        v_res = validate_single_fastq(s.r1_path, sample_records_to_check=5000)
                        sample_checks[s.sample_id] = v_res
                    except Exception as e:
                        findings.append(
                            PipelineFinding(
                                severity="FATAL",
                                code="FASTQ_VALIDATION_ERROR",
                                message=f"Structural single-end FASTQ validation failed for '{s.sample_id}': {e}",
                            )
                        )

        # 4. Canonical Metadata Validation
        val_res = MetadataValidator().validate(self.canonical_metadata)
        int_res = MetadataIntegrityChecker(min_biological_replicates=self.config.min_biological_replicates).check(
            self.canonical_metadata
        )

        metadata_checks = {
            "validator_status": str(val_res.status.value),
            "integrity_status": str(int_res.status.value),
            "integrity_is_valid": int_res.is_valid,
            "replicate_summary": int_res.replicate_summary,
        }

        if val_res.status == ValidationStatus.INVALID:
            findings.append(
                PipelineFinding(
                    severity="FATAL",
                    code="METADATA_INVALID",
                    message=f"Canonical metadata failed validation: {'; '.join([f.message for f in val_res.findings[:2]])}",
                    remediation="Correct required fields in sample metadata to conform with metadata schema.",
                )
            )

        if not int_res.is_valid:
            findings.append(
                PipelineFinding(
                    severity="WARNING",
                    code="BIOLOGICAL_REPLICATION_UNVERIFIED",
                    message="Biological replicate identities are unverified. Inferential GLM and pathway analysis will be safely gated.",
                )
            )

        # 5. Check reusable checkpoints
        reusable_stages = self._detect_reusable_stages()

        is_valid = not any(f.severity == "FATAL" for f in findings)
        return PreflightReport(
            is_valid=is_valid,
            findings=tuple(findings),
            tool_inventory=tools,
            hardware_baseline=hw.to_dict(),
            sample_validation=sample_checks,
            reference_validation=ref_checks,
            metadata_validation=metadata_checks,
            reusable_checkpoints=tuple(reusable_stages),
            estimated_disk_gb=est_disk_gb,
            estimated_ram_gb=self.profile.physical_ram_gb_min,
        )

    def _compute_plan_signature(self, stage_ids: list[str]) -> str:
        items = [
            f"study:{self.config.study_id}",
            f"source_type:{self.config.source_type}",
            f"experiment_id:{self.config.experiment_id}",
            f"layout:{self.config.read_layout}",
            f"strand:{self.config.strandedness}",
            f"threads:{self.config.threads}",
            f"profile:{self.profile.name}",
            f"min_bio_reps:{self.config.min_biological_replicates}",
            f"enforce_res:{self.config.enforce_resources}",
            f"species:{self.references.species}",
            f"fasta:{compute_file_content_hash(self.references.genome_fasta)}",
            f"gtf:{compute_file_content_hash(self.references.gtf_annotation)}",
            f"sjdb_overhang:{self.references.sjdb_overhang}",
            f"star_idx:{self.references.star_index_dir or ''}",
            f"stages:{','.join(stage_ids)}",
        ]
        for s in sorted(self.samples, key=lambda x: x.sample_id):
            r1_hash = compute_file_content_hash(s.r1_path) if s.r1_path and s.r1_path.is_file() else (s.sra_accession or "missing")
            r2_hash = compute_file_content_hash(s.r2_path) if s.r2_path and s.r2_path.is_file() else "none"
            items.append(f"sample:{s.sample_id}:{s.read_layout}:{s.sra_accession or 'none'}:{r1_hash}:{r2_hash}")
        return hashlib.sha256(";".join(items).encode("utf-8")).hexdigest()

    def plan(self) -> ProcessingPlan:
        """Compile a complete dry-run processing plan with planned operations and resource estimates."""
        preflight = self.validate_inputs()
        stages: list[dict[str, Any]] = []

        # 0. SRA Retrieval if applicable
        for s in self.samples:
            if s.sra_accession and not s.r1_path:
                stages.append({
                    "stage_id": f"00_sra_{s.sample_id}",
                    "name": f"SRA Extraction ({s.sample_id})",
                    "sample_id": s.sample_id,
                    "accession": s.sra_accession,
                    "target_dir": str(self.sra_dir / s.sample_id),
                })

        # 1. References setup
        star_idx_target = self.references.star_index_dir or (self.out_dir / "star_index")
        stages.append({
            "stage_id": "01_references",
            "name": "Reference Verification & STAR Index",
            "action": "Verify FASTA/GTF and build STAR index if not present",
            "target": str(star_idx_target),
            "can_reuse": star_idx_target.is_dir() and (star_idx_target / "Genome").is_file() and (star_idx_target / "SA").is_file(),
        })

        # 2. Per-sample trimming
        for s in self.samples:
            r2_target = self.trimmed_dir / f"{s.sample_id}_trimmed_R2.fastq.gz" if s.read_layout == "paired_end" else None
            stages.append({
                "stage_id": f"02_fastp_{s.sample_id}",
                "name": f"fastp QC & Trimming ({s.sample_id})",
                "sample_id": s.sample_id,
                "read_layout": s.read_layout,
                "out_r1": str(self.trimmed_dir / f"{s.sample_id}_trimmed_R1.fastq.gz"),
                "out_r2": str(r2_target) if r2_target else None,
                "report_json": str(self.qc_dir / f"{s.sample_id}_fastp.json"),
            })

        # 3. Per-sample alignment
        for s in self.samples:
            bam_out = self.align_dir / f"{s.sample_id}_Aligned.sortedByCoord.out.bam"
            stages.append({
                "stage_id": f"03_star_{s.sample_id}",
                "name": f"STAR Coordinate-Sorted Alignment ({s.sample_id})",
                "sample_id": s.sample_id,
                "target_bam": str(bam_out),
            })

        # 4. Strandedness diagnostics
        stages.append({
            "stage_id": "04_rseqc_strandedness",
            "name": "RSeQC Strandedness Evaluation",
            "target_evidence": str(self.out_dir / "strandedness_evidence.json"),
            "configured_strandedness": self.config.strandedness,
        })

        # 5. featureCounts & Count Matrix
        stages.append({
            "stage_id": "05_featurecounts",
            "name": "featureCounts Quantification & Matrix Conversion",
            "target_matrix": str(self.counts_dir / "counts.tsv"),
            "read_layout": self.config.read_layout,
        })

        # 6. Downstream Universal Ingestion v2
        stages.append({
            "stage_id": "06_downstream_bridge",
            "name": "Downstream QC, Harmonization, Normalization & Replication Gate",
            "inferential_eligible": preflight.metadata_validation.get("integrity_is_valid", False),
        })

        approval_token = self._compute_plan_signature([st["stage_id"] for st in stages])[:16]

        return ProcessingPlan(
            config=self.config,
            samples=tuple(self.samples),
            references=self.references,
            preflight=preflight,
            planned_stages=tuple(stages),
            can_execute=preflight.is_valid,
            approval_required=self.config.require_approval,
            approval_token=approval_token,
        )

    # -------------------------------------------------------------------------
    # Pipeline Execution
    # -------------------------------------------------------------------------

    def run(
        self,
        approval_token: str | None = None,
        dry_run: bool = False,
    ) -> UniversalFastqResult:
        """Execute the universal automated FASTQ pipeline through all scientific modules.

        If dry_run is True, returns planning report without executing compute steps.
        If approval_required is True and approval_token does not match, halts safely.
        """
        plan = self.plan()

        if dry_run:
            return UniversalFastqResult(
                status="DRY_RUN_PLAN_GENERATED",
                study_id=self.config.study_id,
                sample_ids=tuple(s.sample_id for s in self.samples),
                raw_counts_path=None,
                canonical_metadata_path=None,
                qc_reports_dir=self.qc_dir,
                alignments_dir=self.align_dir,
                strandedness_evidence_path=None,
                counting_summary_path=None,
                execution_record_path=None,
                metrics={"planned_stages_count": len(plan.planned_stages)},
                findings=plan.preflight.findings,
                completed_stages=(),
                failed_stage=None,
            )

        # Preflight validation gate
        if plan.preflight.has_fatal_errors:
            return UniversalFastqResult(
                status="BLOCKED_PREFLIGHT_FAILURES",
                study_id=self.config.study_id,
                sample_ids=tuple(s.sample_id for s in self.samples),
                raw_counts_path=None,
                canonical_metadata_path=None,
                qc_reports_dir=None,
                alignments_dir=None,
                strandedness_evidence_path=None,
                counting_summary_path=None,
                execution_record_path=None,
                findings=plan.preflight.findings,
                completed_stages=(),
                failed_stage="preflight_validation",
            )

        # Explicit approval gate before expensive operations
        if self.config.require_approval:
            effective_token = approval_token or self._approved_token
            if not effective_token:
                approval_finding = PipelineFinding(
                    severity="WARNING",
                    code="AWAITING_USER_APPROVAL",
                    message=(
                        f"Pipeline is configured to require explicit confirmation before execution. "
                        f"To proceed, invoke run(approval_token='{plan.approval_token}') or pipeline.approve('{plan.approval_token}')."
                    ),
                )
                return UniversalFastqResult(
                    status="BLOCKED_AWAITING_APPROVAL",
                    study_id=self.config.study_id,
                    sample_ids=tuple(s.sample_id for s in self.samples),
                    raw_counts_path=None,
                    canonical_metadata_path=None,
                    qc_reports_dir=None,
                    alignments_dir=None,
                    strandedness_evidence_path=None,
                    counting_summary_path=None,
                    execution_record_path=None,
                    findings=(approval_finding,),
                    completed_stages=(),
                    failed_stage="approval_gate",
                )
            if effective_token != plan.approval_token:
                stale_finding = PipelineFinding(
                    severity="FATAL",
                    code="STALE_APPROVAL_TOKEN",
                    message=(
                        f"Approval token '{effective_token}' is stale or invalid for the current plan signature ('{plan.approval_token}'). "
                        f"Pipeline inputs, references, or planned stages have changed since approval was obtained. "
                        f"Generate a new plan with pipeline.plan() and approve the updated token."
                    ),
                    remediation=f"Re-run plan() and obtain a fresh approval token (expected: '{plan.approval_token}').",
                )
                return UniversalFastqResult(
                    status="BLOCKED_STALE_APPROVAL",
                    study_id=self.config.study_id,
                    sample_ids=tuple(s.sample_id for s in self.samples),
                    raw_counts_path=None,
                    canonical_metadata_path=None,
                    qc_reports_dir=None,
                    alignments_dir=None,
                    strandedness_evidence_path=None,
                    counting_summary_path=None,
                    execution_record_path=None,
                    findings=(stale_finding,),
                    completed_stages=(),
                    failed_stage="approval_gate",
                )

        # Ensure output directories exist upon approved execution start
        self._ensure_directories()

        t_start = time.time()
        sample_ids = tuple(s.sample_id for s in self.samples)
        completed_stages: list[str] = []
        failed_stage: str | None = None

        # ---------------------------------------------------------------------
        # Stage 0: SRA Retrieval if applicable
        # ---------------------------------------------------------------------
        resolved_fastqs: dict[str, tuple[Path, Path | None]] = {}
        for s in self.samples:
            if s.r1_path is not None:
                resolved_fastqs[s.sample_id] = (s.r1_path, s.r2_path)
            elif s.sra_accession:
                sra_target = self.sra_dir / s.sample_id
                sra_input_hash = self._compute_hash(s.sra_accession)
                sra_config_hash = self._compute_hash(s.sample_id, str(self.config.threads))
                if not self.checkpoint_manager.is_stage_completed(
                    "sra_extract", s.sample_id, [sra_target], sra_input_hash, sra_config_hash
                ):
                    print(f"Retrieving SRA accession {s.sra_accession} for {s.sample_id}...")
                    sra_res = extract_sra(
                        sra_path=s.sra_accession,
                        output_dir=sra_target,
                        threads=self.config.threads,
                    )
                    resolved_fastqs[s.sample_id] = (sra_res.fastq_r1, sra_res.fastq_r2)
                    outputs = [sra_res.fastq_r1] + ([sra_res.fastq_r2] if sra_res.fastq_r2 else [])
                    self.checkpoint_manager.record_checkpoint(
                        stage="sra_extract",
                        identifier=s.sample_id,
                        outputs=outputs,
                        input_hash=sra_input_hash,
                        config_hash=sra_config_hash,
                        metrics=sra_res.to_dict(),
                    )
                else:
                    print(f"[REUSE] SRA extraction checkpoint valid for {s.sample_id}")
                    found = list(sra_target.glob(f"{s.sra_accession}*.fastq*"))
                    f_r1 = next((f for f in found if "_1" in f.name), found[0] if found else None)
                    f_r2 = next((f for f in found if "_2" in f.name), None)
                    if not f_r1:
                        raise FileNotFoundError(f"SRA outputs missing in {sra_target} for {s.sample_id}")
                    resolved_fastqs[s.sample_id] = (f_r1, f_r2)
        if any(s.sra_accession and not s.r1_path for s in self.samples):
            completed_stages.append("sra_extraction")

        # ---------------------------------------------------------------------
        # Stage 1: STAR Reference Index
        # ---------------------------------------------------------------------
        star_idx_dir = self.references.star_index_dir or (self.out_dir / "star_index")
        genome_file = star_idx_dir / "Genome"
        config_hash_ref = self._compute_hash(
            str(self.references.sjdb_overhang),
            str(self.profile.genome_sa_sparse_d),
            str(self.profile.genome_sa_index_nbases),
        )
        input_hash_ref = self._compute_hash(
            compute_file_content_hash(self.references.genome_fasta),
            compute_file_content_hash(self.references.gtf_annotation),
        )
        expected_star_idx = [star_idx_dir, genome_file]

        if not self.checkpoint_manager.is_stage_completed(
            "star_index", "global", expected_star_idx, input_hash_ref, config_hash_ref
        ):
            print(f"Generating STAR genome index in {star_idx_dir}...")
            build_star_index(
                fasta_path=self.references.genome_fasta,
                gtf_path=self.references.gtf_annotation,
                index_dir=star_idx_dir,
                sjdb_overhang=self.references.sjdb_overhang,
                threads=self.config.threads,
                sparse_d=self.profile.genome_sa_sparse_d,
                genome_sa_index_nbases=self.profile.genome_sa_index_nbases,
                limit_ram=self.profile.limit_genome_generate_ram,
            )
            self.checkpoint_manager.record_checkpoint(
                stage="star_index",
                identifier="global",
                outputs=expected_star_idx,
                input_hash=input_hash_ref,
                config_hash=config_hash_ref,
            )
        else:
            print(f"[REUSE] STAR index checkpoint valid in {star_idx_dir}")
        completed_stages.append("star_index")

        # ---------------------------------------------------------------------
        # Stage 2 & 3: fastp Trimming and STAR Alignment per sample
        # ---------------------------------------------------------------------
        bam_paths: list[Path] = []
        fastp_summaries: dict[str, Any] = {}
        alignment_summaries: dict[str, Any] = {}

        for sample in self.samples:
            sid = sample.sample_id
            trimmed_r1 = self.trimmed_dir / f"{sid}_trimmed_R1.fastq.gz"
            trimmed_r2 = self.trimmed_dir / f"{sid}_trimmed_R2.fastq.gz" if sample.read_layout == "paired_end" else None
            json_report = self.qc_dir / f"{sid}_fastp.json"
            html_report = self.qc_dir / f"{sid}_fastp.html"

            # 2A. fastp Preprocessing
            raw_r1, raw_r2 = resolved_fastqs[sid]
            input_hash_fastp = self._compute_hash(
                compute_file_content_hash(raw_r1),
                compute_file_content_hash(raw_r2) if raw_r2 else "",
            )
            config_hash_fastp = self._compute_hash(sid, sample.read_layout, str(self.config.threads))
            expected_fastp = [trimmed_r1, json_report] + ([trimmed_r2] if trimmed_r2 else [])

            if not self.checkpoint_manager.is_stage_completed(
                "fastp", sid, expected_fastp, input_hash_fastp, config_hash_fastp
            ):
                print(f"Executing fastp for {sid}...")
                fastp_res = run_fastp_preprocessing(
                    raw_r1=raw_r1,
                    raw_r2=raw_r2,
                    out_r1=trimmed_r1,
                    out_r2=trimmed_r2,
                    json_report=json_report,
                    html_report=html_report,
                    sample_id=sid,
                    threads=self.config.threads,
                )
                self.checkpoint_manager.record_checkpoint(
                    stage="fastp",
                    identifier=sid,
                    outputs=expected_fastp,
                    input_hash=input_hash_fastp,
                    config_hash=config_hash_fastp,
                    metrics=fastp_res.get("metrics_summary"),
                )
                fastp_summaries[sid] = fastp_res.get("metrics_summary")
            else:
                print(f"[REUSE] fastp checkpoint valid for {sid}")
                if json_report.is_file():
                    try:
                        fastp_summaries[sid] = json.loads(json_report.read_text(encoding="utf-8")).get("summary", {})
                    except Exception:
                        pass
            completed_stages.append(f"fastp_{sid}")

            # 2B. STAR Alignment
            star_prefix = self.align_dir / f"{sid}_"
            final_bam = self.align_dir / f"{sid}_Aligned.sortedByCoord.out.bam"
            final_log = self.align_dir / f"{sid}_Log.final.out"

            input_hash_align = self._compute_hash(
                compute_file_content_hash(trimmed_r1),
                compute_file_content_hash(trimmed_r2) if trimmed_r2 else "",
                input_hash_ref,
            )
            config_hash_align = self._compute_hash(
                self.profile.bam_sort_mode,
                str(self.profile.external_sort_enabled),
                str(self.config.threads),
            )
            expected_align = [final_bam, final_log]

            if not self.checkpoint_manager.is_stage_completed(
                "star_align", sid, expected_align, input_hash_align, config_hash_align
            ):
                print(f"Executing STAR alignment for {sid}...")
                align_res = align_reads_star(
                    trimmed_r1=trimmed_r1,
                    trimmed_r2=trimmed_r2,
                    index_dir=star_idx_dir,
                    output_prefix=star_prefix,
                    threads=self.config.threads,
                    bam_sort_mode=self.profile.bam_sort_mode,
                    limit_bam_sort_ram=self.profile.limit_bam_sort_ram,
                    external_sort=self.profile.external_sort_enabled,
                    samtools_sort_mem=self.profile.samtools_sort_mem_per_thread,
                )

                # Index BAM with samtools index
                wsl_bam = to_wsl_path(final_bam)
                run_wsl(f"samtools index {wsl_bam}")

                self.checkpoint_manager.record_checkpoint(
                    stage="star_align",
                    identifier=sid,
                    outputs=expected_align,
                    input_hash=input_hash_align,
                    config_hash=config_hash_align,
                    metrics=align_res.get("metrics"),
                )
                alignment_summaries[sid] = align_res.get("metrics")
            else:
                print(f"[REUSE] STAR alignment checkpoint valid for {sid}")
                if final_log.is_file():
                    alignment_summaries[sid] = parse_star_final_log(final_log)

            bam_paths.append(final_bam)
            completed_stages.append(f"star_align_{sid}")

        # ---------------------------------------------------------------------
        # Stage 4: RSeQC Strandedness Assessment
        # ---------------------------------------------------------------------
        evidence_file = self.out_dir / "strandedness_evidence.json"
        bed12_path = self.references.bed12_annotation or (self.out_dir / "annotation.bed12")
        if not bed12_path.is_file():
            print(f"Generating deterministic BED12 from GTF...")
            generate_bed12_from_gtf(self.references.gtf_annotation, bed12_path)

        per_sample_strandedness: dict[str, Any] = {}
        if self.config.strandedness is not None:
            # User or study explicit decision with justification
            decision = StrandednessDecision(
                strandedness_flag=self.config.strandedness,
                strand_mode={0: "unstranded", 1: "forward", 2: "reverse"}.get(self.config.strandedness, "unstranded"),
                confidence="HIGH",
                mean_forward=0.0,
                mean_reverse=0.0,
                mean_unexplained=0.0,
                consensus_indication=f"USER_CONFIGURED (-s {self.config.strandedness})",
                rationale="Strandedness flag explicitly configured by investigator/study protocol.",
            )
        else:
            # Empirical diagnosis via RSeQC infer_experiment.py
            print("Evaluating empirical library strandedness across BAMs...")
            for sid, bam_p in zip(sample_ids, bam_paths):
                exp_res = run_infer_experiment(bam_p, bed12_path)
                per_sample_strandedness[sid] = exp_res

            decision = infer_strandedness_decision(per_sample_strandedness)

        strandedness_evidence = {
            "study_id": self.config.study_id,
            "decision": decision.to_dict(),
            "per_sample_strandedness": per_sample_strandedness,
            "recorded_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        }
        evidence_file.write_text(json.dumps(strandedness_evidence, indent=2), encoding="utf-8")
        completed_stages.append("rseqc_strandedness")

        if decision.strandedness_flag == -1:
            return UniversalFastqResult(
                status="BLOCKED_AWAITING_STRANDEDNESS_DECISION",
                study_id=self.config.study_id,
                sample_ids=sample_ids,
                raw_counts_path=None,
                canonical_metadata_path=None,
                qc_reports_dir=self.qc_dir,
                alignments_dir=self.align_dir,
                strandedness_evidence_path=evidence_file,
                counting_summary_path=None,
                execution_record_path=None,
                metrics={"rationale": decision.rationale},
                completed_stages=tuple(completed_stages),
                failed_stage="rseqc_strandedness",
            )

        # ---------------------------------------------------------------------
        # Stage 5: featureCounts Quantification & Matrix Conversion
        # ---------------------------------------------------------------------
        raw_fc_output = self.counts_dir / "featurecounts_raw.txt"
        summary_fc_output = self.counts_dir / "featurecounts_raw.txt.summary"
        counts_tsv = self.counts_dir / "counts.tsv"

        is_pe = self.config.read_layout == "paired_end"
        print(f"Running featureCounts with strandedness -s {decision.strandedness_flag} (paired={is_pe})...")
        run_featurecounts(
            bam_paths=bam_paths,
            gtf_path=self.references.gtf_annotation,
            output_counts_path=raw_fc_output,
            strandedness=decision.strandedness_flag,
            threads=self.config.threads,
            paired_end=is_pe,
        )

        # Convert to standardized count matrix
        convert_featurecounts_to_matrix(
            raw_featurecounts_path=raw_fc_output,
            output_matrix_path=counts_tsv,
            approved_runs=sample_ids,
        )
        counting_metrics = parse_featurecounts_summary(summary_fc_output)
        summary_json = self.counts_dir / "counting_summary.json"
        summary_json.write_text(json.dumps(counting_metrics, indent=2), encoding="utf-8")
        completed_stages.append("featurecounts")

        # ---------------------------------------------------------------------
        # Stage 6: Canonical Metadata Serialization & Traceability
        # ---------------------------------------------------------------------
        meta_json = self.out_dir / "canonical_metadata.json"
        meta_json.write_text(json.dumps(self.canonical_metadata, indent=2), encoding="utf-8")
        completed_stages.append("canonical_metadata")

        # ---------------------------------------------------------------------
        # Stage 7: Downstream Bridge Integration
        # ---------------------------------------------------------------------
        downstream_res = None
        try:
            downstream_res = execute_downstream_pipeline(
                counts_matrix_path=counts_tsv,
                runs=sample_ids,
                study_id=self.config.study_id,
                output_dir=self.out_dir / "downstream_analysis",
                canonical_metadata=self.canonical_metadata,
            )
        except Exception as e:
            print(f"Downstream bridge execution error: {e}")

        # Determine granular downstream status
        failed_stage = None
        if downstream_res is None:
            final_status = "FAILED_DOWNSTREAM_EXCEPTION"
            failed_stage = "downstream_bridge"
        elif downstream_res.status in ("DOWNSTREAM_COMPLETED", "DOWNSTREAM_COMPLETED_REPLICATED"):
            final_status = "COMPLETED"
            completed_stages.append("downstream_bridge")
        elif "METADATA" in downstream_res.status or getattr(downstream_res, "metadata_validation_status", "").upper() == "INVALID":
            final_status = "FAILED_DOWNSTREAM_METADATA"
            failed_stage = "downstream_bridge"
        elif "QC" in downstream_res.status or not getattr(downstream_res, "qc_passed", True):
            final_status = "FAILED_DOWNSTREAM_QC"
            failed_stage = "downstream_bridge"
        elif "HARMONIZATION" in downstream_res.status:
            final_status = "FAILED_DOWNSTREAM_HARMONIZATION"
            failed_stage = "downstream_bridge"
        elif "NORMALIZATION" in downstream_res.status:
            final_status = "FAILED_DOWNSTREAM_NORMALIZATION"
            failed_stage = "downstream_bridge"
        elif "GATED" in downstream_res.status or "BLOCKED" in downstream_res.status or getattr(downstream_res, "deseq2_status", "") == "BLOCKED_UNVERIFIED_BIOLOGICAL_REPLICATION":
            final_status = "COMPLETED_QUANTIFICATION_INFERENTIAL_BLOCKED"
            completed_stages.append("downstream_bridge_exploratory")
        else:
            final_status = "PARTIALLY_COMPLETED"
            completed_stages.append("downstream_bridge_partial")

        # ---------------------------------------------------------------------
        # Stage 8: Execution Record Serialization
        # ---------------------------------------------------------------------
        exec_record_file = self.out_dir / "execution_record.json"
        exec_record = {
            "study_id": self.config.study_id,
            "species": self.references.species,
            "read_layout": self.config.read_layout,
            "sample_ids": list(sample_ids),
            "samples_count": len(sample_ids),
            "strandedness_flag": decision.strandedness_flag,
            "strandedness_mode": decision.strand_mode,
            "execution_duration_seconds": round(time.time() - t_start, 2),
            "profile": self.profile.to_dict(),
            "references": self.references.to_dict(),
            "downstream_status": downstream_res.status if downstream_res else "NOT_RUN",
            "completed_stages": completed_stages,
            "completed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        }
        exec_record_file.write_text(json.dumps(exec_record, indent=2), encoding="utf-8")

        return UniversalFastqResult(
            status=final_status,
            study_id=self.config.study_id,
            sample_ids=sample_ids,
            raw_counts_path=counts_tsv,
            canonical_metadata_path=meta_json,
            qc_reports_dir=self.qc_dir,
            alignments_dir=self.align_dir,
            strandedness_evidence_path=evidence_file,
            counting_summary_path=summary_json,
            execution_record_path=exec_record_file,
            downstream_result=downstream_res,
            completed_stages=tuple(completed_stages),
            failed_stage=failed_stage,
            metrics={
                "fastp": fastp_summaries,
                "star_alignment": alignment_summaries,
                "counting": counting_metrics,
                "strandedness": decision.to_dict(),
            },
        )

    # -------------------------------------------------------------------------
    # Internal Probing Helpers
    # -------------------------------------------------------------------------

    def _probe_hardware(self) -> HardwareBaseline:
        cpu_cores = os.cpu_count() or 4
        ram_gb = 16.0
        try:
            mem_res = run_wsl("grep MemTotal /proc/meminfo")
            if mem_res.returncode == 0:
                ram_kb = int(mem_res.stdout.split()[1])
                ram_gb = ram_kb / 1024 / 1024
        except Exception:
            pass

        data_free = 100.0
        try:
            _, _, free_b = shutil.disk_usage(self.out_dir)
            data_free = free_b / (1024**3)
        except Exception:
            pass

        return HardwareBaseline(
            cpu_cores=cpu_cores,
            cpu_threads=cpu_cores,
            cpu_model="Generic Host CPU",
            physical_ram_gb=ram_gb,
            swap_gb=4.0,
            data_dir_free_gb=data_free,
            ref_dir_free_gb=data_free,
            os_kernel=os.name,
            is_wsl=True,
            architecture="x86_64",
        )

    def _probe_tools(self) -> dict[str, dict[str, Any]]:
        global _PROBED_TOOLS_CACHE
        if _PROBED_TOOLS_CACHE is not None:
            return dict(_PROBED_TOOLS_CACHE)

        tools = ["fastp", "STAR", "featureCounts", "samtools", "infer_experiment.py", "fasterq-dump"]
        result: dict[str, dict[str, Any]] = {}
        for t in tools:
            w_res = run_wsl(f"which {t} || true")
            path = w_res.stdout.strip()
            if path:
                v_res = run_wsl(f"{t} --version 2>&1 || true")
                ver = v_res.stdout.strip().splitlines()[0] if v_res.stdout else "available"
                result[t] = {"available": True, "path": path, "version": ver}
            else:
                result[t] = {"available": False, "path": None, "version": None}
        _PROBED_TOOLS_CACHE = result
        return result

    def _detect_reusable_stages(self) -> list[str]:
        stages = []
        star_idx = self.references.star_index_dir or (self.out_dir / "star_index")
        if star_idx.is_dir() and (star_idx / "Genome").is_file():
            stages.append("star_index")
        return stages

    def _compute_hash(self, *items: str) -> str:
        hasher = hashlib.sha256()
        for item in items:
            hasher.update(str(item).encode("utf-8"))
        return hasher.hexdigest()
