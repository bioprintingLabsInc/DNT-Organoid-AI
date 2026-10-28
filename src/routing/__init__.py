"""Modality-aware routing and format detection for heterogeneous transcriptomic assets."""

from .detector import HDF5_MAGIC, SRA_MAGIC, detect_format
from .dispatcher import (
    RouteExecutionPlan,
    RouteExecutionResult,
    dispatch_route_execution,
    plan_route_execution,
)
from .router import route_metadata_assets

__all__ = [
    "HDF5_MAGIC",
    "SRA_MAGIC",
    "RouteExecutionPlan",
    "RouteExecutionResult",
    "detect_format",
    "dispatch_route_execution",
    "plan_route_execution",
    "route_metadata_assets",
]
