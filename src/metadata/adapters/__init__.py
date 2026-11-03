"""Source metadata adapters package."""

from .base import MappingFinding, MetadataSourceAdapter, RawSampleRecord, StudyContext
from .geo import GeoMetadataAdapter
from .json_adapter import JsonMetadataAdapter
from .magetab import MageTabAdapter
from .tabular import TabularMetadataAdapter

__all__ = [
    "MetadataSourceAdapter",
    "RawSampleRecord",
    "StudyContext",
    "MappingFinding",
    "TabularMetadataAdapter",
    "MageTabAdapter",
    "GeoMetadataAdapter",
    "JsonMetadataAdapter",
]
