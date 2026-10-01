"""Treatment-versus-Matched-Control Contrast Builder v1.

Constructs, validates, and gates molecular-response contrasts from explicit metadata relationships
and verified normalization cohorts for downstream negative-binomial differential expression modeling.
"""

from .builder import (
    BUILDER_VERSION,
    ResponseContrastBuilder,
    build_response_contrasts,
    extract_metadata_value,
    normalize_age_string,
    resolve_developmental_age,
)
from .models import (
    ContrastExposure,
    ContrastStatus,
    Finding,
    MolecularResponseContrast,
    ResponseContrastDataset,
    Severity,
    Status,
    ordered_findings,
    status_for,
)

__all__ = [
    "BUILDER_VERSION",
    "ContrastExposure",
    "ContrastStatus",
    "Finding",
    "MolecularResponseContrast",
    "ResponseContrastBuilder",
    "ResponseContrastDataset",
    "Severity",
    "Status",
    "build_response_contrasts",
    "extract_metadata_value",
    "normalize_age_string",
    "ordered_findings",
    "resolve_developmental_age",
    "status_for",
]
