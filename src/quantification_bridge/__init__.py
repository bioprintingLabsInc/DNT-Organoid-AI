"""Governed Raw-Read Quantification Bridge.

Upstream preprocessing and quantification module operating strictly outside
frozen Stages 1-7, Stage 8C, and Stage 8A.1. Registered as the execution engine
for Route.RAW_READ_PIPELINE.
"""

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
)

__version__ = "0.3.0"

__all__ = [
    "ExecutionProfile",
    "FastqPair",
    "FastqPairingError",
    "HardwareBaseline",
    "PreflightCheckResult",
    "PreflightFinding",
    "QuantificationBridgePipeline",
    "ToolInventory",
    "determine_read_mate",
    "evaluate_qc_metric",
    "extract_sample_stem",
    "get_execution_profile",
    "group_fastq_pairs",
    "load_execution_profiles",
    "probe_hardware_baseline",
    "probe_tool_inventory",
    "resolve_pipeline_paths",
    "run_profile_preflight_checks",
]
