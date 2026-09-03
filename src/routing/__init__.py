"""Format detection and modality-aware routing."""

from .detector import HDF5_MAGIC, SRA_MAGIC, detect_format
from .router import route_metadata_assets

__all__ = ["HDF5_MAGIC", "SRA_MAGIC", "detect_format", "route_metadata_assets"]
