"""Immutable deterministic data models for Stage 8A.1 Condition-Level Gene Response Registry.

The fundamental scientific unit of Stage 8A.1 is:
    one gene × one treatment-control contrast

Preserves exact Stage 7 differential molecular response statistics (log2FC, SE, Wald,
p-value, adjusted p-value) alongside complete chemical, exposure, biological, and
provenance context, without performing cross-study synthesis or consuming DNT labels.
"""

from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from enum import StrEnum
import math
from typing import Any


class ResponseDirection(StrEnum):
    """Descriptive response direction based strictly on the sign of log2FC.

    This is non-evaluative: UP/DOWN does NOT imply biological significance or DNT outcome.
    """

    UP = "UP"
    DOWN = "DOWN"
    NO_DIRECTION = "NO_DIRECTION"
    UNDETERMINED = "UNDETERMINED"


class EstimateStatus(StrEnum):
    """Quality and estimation status of the gene-level statistical response."""

    OK = "OK"
    ZERO_COUNTS = "ZERO_COUNTS"
    COOKS_OUTLIER = "COOKS_OUTLIER"
    FIT_ERROR = "FIT_ERROR"
    INVALID_STATISTIC = "INVALID_STATISTIC"


class CrossStudyStatus(StrEnum):
    """Multi-study synthesis status for a gene/compound/regimen evidence unit.

    In Stage 8A.1 (single-dataset foundation), this is strictly SINGLE_STUDY_ONLY.
    """

    SINGLE_STUDY_ONLY = "SINGLE_STUDY_ONLY"
    MULTI_STUDY_CONTEXT_UNDER_REVIEW = "MULTI_STUDY_CONTEXT_UNDER_REVIEW"
    NONCOMPARABLE_CONTEXT = "NONCOMPARABLE_CONTEXT"
    COMPARABLE_INDEPENDENT_STUDIES = "COMPARABLE_INDEPENDENT_STUDIES"
    INSUFFICIENT_METADATA = "INSUFFICIENT_METADATA"


@dataclass(frozen=True, slots=True)
class ConditionGeneResponseRecord:
    """The fundamental condition-level gene response unit: one gene × one treatment-control contrast.

    Preserves exact Stage 7 differential statistics alongside complete experimental,
    chemical, biological, and provenance context.
    """

    # Gene Identity
    canonical_ensembl_gene_id: str
    gene_symbol: str | None

    # Study & Independence Identity
    study_id: str
    study_family_id: str
    independent_study_id: str
    experiment_id: str
    contrast_id: str

    # Treatment & Matched Control Conditions
    treatment_condition_id: str
    matched_control_condition_ids: tuple[str, ...]

    # Chemical / Compound Context
    compound_id: str
    compound_name: str

    # Exposure Context
    dose_value: float | None
    dose_unit: str | None
    dose_original: str | None
    exposure_duration: float | None
    exposure_duration_unit: str | None
    developmental_stage: str | None
    model_type: str | None

    # Biological Context
    biological_source_ids: tuple[str, ...]
    organoid_context_ids: tuple[str, ...]
    treatment_replicate_count: int
    control_replicate_count: int
    treatment_sample_ids: tuple[str, ...]
    control_sample_ids: tuple[str, ...]
    treatment_biological_replicate_ids: tuple[str, ...]
    control_biological_replicate_ids: tuple[str, ...]

    # Statistical Gene Response (Preserved exactly from Stage 7)
    base_mean: float
    log2_fold_change: float | None
    log2_fold_change_se: float | None
    wald_statistic: float | None
    p_value: float | None
    adjusted_p_value: float | None
    response_direction: ResponseDirection

    # QC, Status & Cross-Study Synthesis Boundary
    estimate_status: EstimateStatus
    missing_reason: str | None
    cross_study_status: CrossStudyStatus = CrossStudyStatus.SINGLE_STUDY_ONLY

    # Provenance
    source_stage: str = "Stage 7B"
    source_contrast_id: str = ""
    de_version: str = "1.0.0"
    registry_version: str = "1.0.0"
    input_hashes: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Convert record to dictionary with primitive types for serialization."""
        return {
            "canonical_ensembl_gene_id": self.canonical_ensembl_gene_id,
            "gene_symbol": self.gene_symbol,
            "study_id": self.study_id,
            "study_family_id": self.study_family_id,
            "independent_study_id": self.independent_study_id,
            "experiment_id": self.experiment_id,
            "contrast_id": self.contrast_id,
            "treatment_condition_id": self.treatment_condition_id,
            "matched_control_condition_ids": list(self.matched_control_condition_ids),
            "compound_id": self.compound_id,
            "compound_name": self.compound_name,
            "dose_value": self.dose_value,
            "dose_unit": self.dose_unit,
            "dose_original": self.dose_original,
            "exposure_duration": self.exposure_duration,
            "exposure_duration_unit": self.exposure_duration_unit,
            "developmental_stage": self.developmental_stage,
            "model_type": self.model_type,
            "biological_source_ids": list(self.biological_source_ids),
            "organoid_context_ids": list(self.organoid_context_ids),
            "treatment_replicate_count": self.treatment_replicate_count,
            "control_replicate_count": self.control_replicate_count,
            "treatment_sample_ids": list(self.treatment_sample_ids),
            "control_sample_ids": list(self.control_sample_ids),
            "treatment_biological_replicate_ids": list(self.treatment_biological_replicate_ids),
            "control_biological_replicate_ids": list(self.control_biological_replicate_ids),
            "base_mean": self.base_mean,
            "log2_fold_change": self.log2_fold_change,
            "log2_fold_change_se": self.log2_fold_change_se,
            "wald_statistic": self.wald_statistic,
            "p_value": self.p_value,
            "adjusted_p_value": self.adjusted_p_value,
            "response_direction": str(self.response_direction),
            "estimate_status": str(self.estimate_status),
            "missing_reason": self.missing_reason,
            "cross_study_status": str(self.cross_study_status),
            "source_stage": self.source_stage,
            "source_contrast_id": self.source_contrast_id,
            "de_version": self.de_version,
            "registry_version": self.registry_version,
            "input_hashes": dict(self.input_hashes),
        }


@dataclass(frozen=True, slots=True)
class ConditionGeneResponseRegistry:
    """Immutable collection of ConditionGeneResponseRecord objects.

    Represents the full condition-level gene-response feature space for one or more contrasts.
    """

    records: tuple[ConditionGeneResponseRecord, ...]
    registry_id: str
    registry_version: str = "1.0.0"
    created_at: str = ""
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def total_records(self) -> int:
        return len(self.records)

    @property
    def unique_canonical_genes(self) -> set[str]:
        return {r.canonical_ensembl_gene_id for r in self.records}

    @property
    def unique_contrasts(self) -> set[str]:
        return {r.contrast_id for r in self.records}

    @property
    def unique_independent_studies(self) -> set[str]:
        return {r.independent_study_id for r in self.records}

    @property
    def unique_compounds(self) -> set[str]:
        return {r.compound_id for r in self.records}

    @property
    def eligible_records(self) -> tuple[ConditionGeneResponseRecord, ...]:
        """Records with statistically valid (finite) estimates from Stage 7."""
        return tuple(r for r in self.records if r.estimate_status == EstimateStatus.OK)

    def get_records_for_gene(self, gene_id: str) -> tuple[ConditionGeneResponseRecord, ...]:
        return tuple(r for r in self.records if r.canonical_ensembl_gene_id == gene_id)

    def get_records_for_contrast(self, contrast_id: str) -> tuple[ConditionGeneResponseRecord, ...]:
        return tuple(r for r in self.records if r.contrast_id == contrast_id)

    def get_records_for_compound(self, compound_id: str) -> tuple[ConditionGeneResponseRecord, ...]:
        return tuple(r for r in self.records if r.compound_id == compound_id)
