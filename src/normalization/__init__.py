"""Bulk RNA-seq Normalization v1 module."""

from .cohort import NormalizationCohortBuilder
from .config import NormalizationConfig, load_config
from .models import (
    CohortNormalizationMetrics,
    CollisionDisposition,
    Finding,
    NormalizationCohort,
    NormalizationMethod,
    NormalizationReviewDisposition,
    NormalizedCohortResult,
    NormalizedCountMatrix,
    NormalizedDataset,
    ResponseEligibility,
    SampleSizeFactor,
    Severity,
    Status,
    ordered,
    status_for,
)
from .normalizer import BulkNormalizer, normalize_bulk_dataset
from .r_bridge import REnvironmentInfo, estimate_size_factors_r, probe_r_environment

__all__ = [
    "BulkNormalizer",
    "CohortNormalizationMetrics",
    "CollisionDisposition",
    "Finding",
    "NormalizationCohort",
    "NormalizationCohortBuilder",
    "NormalizationConfig",
    "NormalizationMethod",
    "NormalizationReviewDisposition",
    "NormalizedCohortResult",
    "NormalizedCountMatrix",
    "NormalizedDataset",
    "REnvironmentInfo",
    "ResponseEligibility",
    "SampleSizeFactor",
    "Severity",
    "Status",
    "estimate_size_factors_r",
    "load_config",
    "normalize_bulk_dataset",
    "ordered",
    "probe_r_environment",
    "status_for",
]
