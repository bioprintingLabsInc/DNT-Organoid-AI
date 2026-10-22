"""Stage 8A.1 — Condition-Level Gene Response Registry.

The fundamental scientific unit is:
    one gene × one treatment-control contrast

Preserves exact Stage 7 differential molecular response statistics (log2FC, SE, Wald,
p-value, adjusted p-value) alongside complete chemical, exposure, biological, and
provenance context, without performing cross-study synthesis or consuming DNT labels.
"""

from src.features.cross_study.builder import (
    build_condition_gene_response_record,
    build_records_for_contrast,
    derive_estimate_status,
    derive_response_direction,
    normalize_independent_study_id,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
    CrossStudyStatus,
    EstimateStatus,
    ResponseDirection,
)
from src.features.cross_study.qc import (
    RegistryQCSummary,
    audit_registry_qc,
)
from src.features.cross_study.registry import (
    build_condition_gene_response_registry,
)
from src.features.cross_study.serialization import (
    CANONICAL_LONG_COLUMNS,
    export_registry_to_csv,
    export_registry_to_dataframe,
)

__all__ = [
    "CANONICAL_LONG_COLUMNS",
    "ConditionGeneResponseRecord",
    "ConditionGeneResponseRegistry",
    "CrossStudyStatus",
    "EstimateStatus",
    "RegistryQCSummary",
    "ResponseDirection",
    "audit_registry_qc",
    "build_condition_gene_response_record",
    "build_condition_gene_response_registry",
    "build_records_for_contrast",
    "derive_estimate_status",
    "derive_response_direction",
    "export_registry_to_csv",
    "export_registry_to_dataframe",
    "normalize_independent_study_id",
]
