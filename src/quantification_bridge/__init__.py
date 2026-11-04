"""Governed Raw-Read Quantification Bridge.

Upstream preprocessing and quantification module operating strictly outside
frozen Stages 1-7, Stage 8C, and Stage 8A.1. Registered as the execution engine
for Route.RAW_READ_PIPELINE.
"""

from .checkpoints import (
    CheckpointManager,
    ChecksumVerificationError,
    InsufficientStorageCapacityError,
    check_storage_capacity,
    sync_and_verify_directory,
    sync_and_verify_file,
)
from .manifest import APPROVED_2_SAMPLE_TECHNICAL_PILOT, get_approved_technical_pilot_runs
from .pairing import FastqPair, FastqPairingError, determine_read_mate, extract_sample_stem, group_fastq_pairs
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

__version__ = "0.3.1"

__all__ = [
    "APPROVED_2_SAMPLE_TECHNICAL_PILOT",
    "CheckpointManager",
    "ChecksumVerificationError",
    "ExecutionProfile",
    "FastqPair",
    "FastqPairingError",
    "HardwareBaseline",
    "InsufficientStorageCapacityError",
    "PreflightCheckResult",
    "PreflightFinding",
    "QuantificationBridgePipeline",
    "ToolInventory",
    "check_storage_capacity",
    "determine_read_mate",
    "evaluate_qc_metric",
    "extract_sample_stem",
    "get_approved_technical_pilot_runs",
    "get_execution_profile",
    "group_fastq_pairs",
    "load_execution_profiles",
    "probe_hardware_baseline",
    "probe_tool_inventory",
    "resolve_pipeline_paths",
    "run_profile_preflight_checks",
    "sync_and_verify_directory",
    "sync_and_verify_file",
    "validate_pilot_configuration",
]
