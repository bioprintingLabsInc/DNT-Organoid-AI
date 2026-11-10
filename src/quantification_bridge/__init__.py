"""Governed Raw-Read Quantification Bridge.

Upstream preprocessing and quantification module operating strictly outside
frozen Stages 1-7, Stage 8C, and Stage 8A.1. Registered as the execution engine
for Route.RAW_READ_PIPELINE.
"""

from .checkpoints import (
    CheckpointManager,
    ChecksumVerificationError,
    InsufficientStorageCapacityError,
    ManifestVerificationError,
    UnverifiedStorageQuotaError,
    archive_and_sync_directory,
    check_storage_capacity,
    restore_and_verify_archive,
    sync_and_verify_directory,
    sync_and_verify_file,
)
from .manifest import APPROVED_2_SAMPLE_TECHNICAL_PILOT, get_approved_technical_pilot_runs
from .library_type import LibraryType
from .pairing import (
    FastqPair,
    FastqPairingError,
    FastqSingle,
    determine_read_mate,
    extract_sample_stem,
    group_fastq_pairs,
)
from .pipeline import QuantificationBridgePipeline
from .profiles import (
    ExecutionProfile,
    HardwareBaseline,
    PreflightCheckResult,
    PreflightFinding,
    ToolInventory,
    evaluate_qc_metric,
    get_execution_profile,
    load_execution_profiles,
    probe_hardware_baseline,
    probe_tool_inventory,
    resolve_pipeline_paths,
    run_profile_preflight_checks,
    validate_pilot_configuration,
)
from .universal_pipeline import (
    CheckpointManager as UniversalCheckpointManager,
    PipelineCheckpointManager,
    PipelineFinding,
    PreflightReport,
    ProcessingPlan,
    ReferenceConfig,
    SampleInput,
    SraCheckpointVerification,
    UniversalFastqPipeline,
    UniversalFastqResult,
    UniversalPipelineConfig,
    compute_full_file_sha256,
    load_sample_metadata,
)
from .validator import validate_paired_fastqs, validate_single_fastq

__version__ = "0.4.0"

__all__ = [
    "APPROVED_2_SAMPLE_TECHNICAL_PILOT",
    "CheckpointManager",
    "ChecksumVerificationError",
    "ExecutionProfile",
    "FastqPair",
    "FastqPairingError",
    "FastqSingle",
    "HardwareBaseline",
    "InsufficientStorageCapacityError",
    "LibraryType",
    "ManifestVerificationError",
    "PipelineCheckpointManager",
    "PipelineFinding",
    "PreflightCheckResult",
    "PreflightFinding",
    "PreflightReport",
    "ProcessingPlan",
    "QuantificationBridgePipeline",
    "ReferenceConfig",
    "SampleInput",
    "SraCheckpointVerification",
    "ToolInventory",
    "UniversalCheckpointManager",
    "UniversalFastqPipeline",
    "UniversalFastqResult",
    "UniversalPipelineConfig",
    "UnverifiedStorageQuotaError",
    "archive_and_sync_directory",
    "check_storage_capacity",
    "compute_full_file_sha256",
    "determine_read_mate",
    "evaluate_qc_metric",
    "extract_sample_stem",
    "get_approved_technical_pilot_runs",
    "get_execution_profile",
    "group_fastq_pairs",
    "load_execution_profiles",
    "load_sample_metadata",
    "probe_hardware_baseline",
    "probe_tool_inventory",
    "resolve_pipeline_paths",
    "restore_and_verify_archive",
    "run_profile_preflight_checks",
    "sync_and_verify_directory",
    "sync_and_verify_file",
    "validate_paired_fastqs",
    "validate_pilot_configuration",
    "validate_single_fastq",
]
