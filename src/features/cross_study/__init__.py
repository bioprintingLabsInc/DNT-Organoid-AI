"""Stage 8A.1 & 8A.2 — Condition-Level Gene Response Registry and Cross-Study Comparability Framework.

Stage 8A.1:
The fundamental scientific unit is:
    one gene x one treatment-control contrast
Preserves exact Stage 7 differential molecular response statistics (log2FC, SE, Wald,
p-value, adjusted p-value) alongside complete chemical, exposure, biological, and
provenance context, without performing cross-study synthesis or consuming DNT labels.

Stage 8A.2:
Cross-Study Comparability Evaluation and Descriptive Direction Alignment.
Enforces non-comparable context gating and strictly prohibits unweighted averaging,
pooled effect estimation, or synthetic variance manufacturing across non-equivalent regimens.
"""

from src.features.cross_study.builder import (
    build_condition_gene_response_record,
    build_records_for_contrast,
    derive_estimate_status,
    derive_response_direction,
    normalize_independent_study_id,
)
from src.features.cross_study.comparability import (
    compare_cross_study_gene_responses,
    derive_direction_concordance,
    evaluate_study_comparability,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
    CrossStudyComparisonDataset,
    CrossStudyStatus,
    DirectionConcordance,
    EstimateStatus,
    GeneCrossStudyComparisonRecord,
    ResponseDirection,
    StudyComparabilityAssessment,
)
from src.features.cross_study.qc import (
    ComparisonQCSummary,
    RegistryQCSummary,
    audit_comparison_qc,
    audit_registry_qc,
)
from src.features.cross_study.registry import (
    build_condition_gene_response_registry,
)
from src.features.cross_study.serialization import (
    CANONICAL_COMPARISON_COLUMNS,
    CANONICAL_LONG_COLUMNS,
    export_comparison_to_csv,
    export_comparison_to_dataframe,
    export_registry_to_csv,
    export_registry_to_dataframe,
)

__all__ = [
    "CANONICAL_COMPARISON_COLUMNS",
    "CANONICAL_LONG_COLUMNS",
    "ComparisonQCSummary",
    "ConditionGeneResponseRecord",
    "ConditionGeneResponseRegistry",
    "CrossStudyComparisonDataset",
    "CrossStudyStatus",
    "DirectionConcordance",
    "EstimateStatus",
    "GeneCrossStudyComparisonRecord",
    "RegistryQCSummary",
    "ResponseDirection",
    "StudyComparabilityAssessment",
    "audit_comparison_qc",
    "audit_registry_qc",
    "build_condition_gene_response_record",
    "build_condition_gene_response_registry",
    "build_records_for_contrast",
    "compare_cross_study_gene_responses",
    "derive_direction_concordance",
    "derive_estimate_status",
    "derive_response_direction",
    "evaluate_study_comparability",
    "export_comparison_to_csv",
    "export_comparison_to_dataframe",
    "export_registry_to_csv",
    "export_registry_to_dataframe",
    "normalize_independent_study_id",
]
