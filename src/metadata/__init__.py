"""Canonical metadata validation, standardization, and import API."""

from .adapters import (
    GeoMetadataAdapter,
    JsonMetadataAdapter,
    MageTabAdapter,
    MappingFinding,
    MetadataSourceAdapter,
    RawSampleRecord,
    StudyContext,
    TabularMetadataAdapter,
)
from .integrity import IntegrityCheckResult, IntegrityFinding, MetadataIntegrityChecker
from .mapper import CanonicalMetadataMapper, FieldMappingRule
from .result import Finding, Severity, ValidationResult, ValidationStatus
from .schema import DEFAULT_SCHEMA_PATH, load_schema
from .units import (
    NormalizedDose,
    NormalizedDuration,
    UnitConversionStatus,
    normalize_concentration,
    normalize_duration,
)
from .validator import MetadataValidator

__all__ = [
    "DEFAULT_SCHEMA_PATH",
    "Finding",
    "MetadataValidator",
    "Severity",
    "ValidationResult",
    "ValidationStatus",
    "load_schema",
    "NormalizedDose",
    "NormalizedDuration",
    "UnitConversionStatus",
    "normalize_concentration",
    "normalize_duration",
    "MetadataSourceAdapter",
    "RawSampleRecord",
    "StudyContext",
    "MappingFinding",
    "TabularMetadataAdapter",
    "MageTabAdapter",
    "GeoMetadataAdapter",
    "JsonMetadataAdapter",
    "CanonicalMetadataMapper",
    "FieldMappingRule",
    "MetadataIntegrityChecker",
    "IntegrityFinding",
    "IntegrityCheckResult",
]
