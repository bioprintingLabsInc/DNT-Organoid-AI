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

class DirectionConcordance(StrEnum):
    """Descriptive concordance of gene response directions between two independent studies."""

    CONCORDANT_UP = "CONCORDANT_UP"
    CONCORDANT_DOWN = "CONCORDANT_DOWN"
    DISCORDANT = "DISCORDANT"
    NO_DIRECTION = "NO_DIRECTION"
    UNDETERMINED = "UNDETERMINED"


@dataclass(frozen=True, slots=True)
class StudyComparabilityAssessment:
    """Scientific comparability evaluation between two independent experimental conditions.

    Evaluates whether two independent study conditions for the same compound share
    identical/comparable regimen, developmental stage, model platform, and verified replicates.
    """

    study1_id: str
    study2_id: str
    contrast1_id: str
    contrast2_id: str
    compound_id: str
    compound_name: str
    is_independent_studies: bool
    is_compound_matched: bool
    is_dose_matched: bool
    study1_dose: float | None
    study2_dose: float | None
    dose_unit: str | None
    is_developmental_stage_matched: bool
    study1_stage: str | None
    study2_stage: str | None
    is_duration_matched: bool
    study1_duration: float | None
    study2_duration: float | None
    is_organoid_model_matched: bool
    study1_model: str | None
    study2_model: str | None
    are_replicates_verified: bool
    study1_replicates_verified: bool
    study2_replicates_verified: bool
    is_synthesis_justified: bool
    cross_study_status: CrossStudyStatus
    rationale: str
    differences: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "study1_id": self.study1_id,
            "study2_id": self.study2_id,
            "contrast1_id": self.contrast1_id,
            "contrast2_id": self.contrast2_id,
            "compound_id": self.compound_id,
            "compound_name": self.compound_name,
            "is_independent_studies": self.is_independent_studies,
            "is_compound_matched": self.is_compound_matched,
            "is_dose_matched": self.is_dose_matched,
            "study1_dose": self.study1_dose,
            "study2_dose": self.study2_dose,
            "dose_unit": self.dose_unit,
            "is_developmental_stage_matched": self.is_developmental_stage_matched,
            "study1_stage": self.study1_stage,
            "study2_stage": self.study2_stage,
            "is_duration_matched": self.is_duration_matched,
            "study1_duration": self.study1_duration,
            "study2_duration": self.study2_duration,
            "is_organoid_model_matched": self.is_organoid_model_matched,
            "study1_model": self.study1_model,
            "study2_model": self.study2_model,
            "are_replicates_verified": self.are_replicates_verified,
            "study1_replicates_verified": self.study1_replicates_verified,
            "study2_replicates_verified": self.study2_replicates_verified,
            "is_synthesis_justified": self.is_synthesis_justified,
            "cross_study_status": str(self.cross_study_status),
            "rationale": self.rationale,
            "differences": list(self.differences),
        }


@dataclass(frozen=True, slots=True)
class GeneCrossStudyComparisonRecord:
    """Gene-level cross-study comparative alignment unit.

    Preserves Study 1 and Study 2 condition-level records without averaging or collapsing.
    Evaluates descriptive direction concordance and enforces synthesis boundaries.
    """

    canonical_ensembl_gene_id: str
    gene_symbol: str | None
    study1_id: str
    study1_contrast_id: str
    study1_log2_fold_change: float | None
    study1_log2_fold_change_se: float | None
    study1_wald_statistic: float | None
    study1_p_value: float | None
    study1_adjusted_p_value: float | None
    study1_response_direction: ResponseDirection
    study1_estimate_status: EstimateStatus
    study2_id: str
    study2_contrast_id: str
    study2_log2_fold_change: float | None
    study2_log2_fold_change_se: float | None
    study2_wald_statistic: float | None
    study2_p_value: float | None
    study2_adjusted_p_value: float | None
    study2_response_direction: ResponseDirection
    study2_estimate_status: EstimateStatus
    direction_concordance: DirectionConcordance
    is_concordant: bool | None
    cross_study_status: CrossStudyStatus
    is_synthesis_justified: bool
    pooled_log2_fold_change: float | None = None
    pooled_standard_error: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "canonical_ensembl_gene_id": self.canonical_ensembl_gene_id,
            "gene_symbol": self.gene_symbol,
            "study1_id": self.study1_id,
            "study1_contrast_id": self.study1_contrast_id,
            "study1_log2_fold_change": self.study1_log2_fold_change,
            "study1_log2_fold_change_se": self.study1_log2_fold_change_se,
            "study1_wald_statistic": self.study1_wald_statistic,
            "study1_p_value": self.study1_p_value,
            "study1_adjusted_p_value": self.study1_adjusted_p_value,
            "study1_response_direction": str(self.study1_response_direction),
            "study1_estimate_status": str(self.study1_estimate_status),
            "study2_id": self.study2_id,
            "study2_contrast_id": self.study2_contrast_id,
            "study2_log2_fold_change": self.study2_log2_fold_change,
            "study2_log2_fold_change_se": self.study2_log2_fold_change_se,
            "study2_wald_statistic": self.study2_wald_statistic,
            "study2_p_value": self.study2_p_value,
            "study2_adjusted_p_value": self.study2_adjusted_p_value,
            "study2_response_direction": str(self.study2_response_direction),
            "study2_estimate_status": str(self.study2_estimate_status),
            "direction_concordance": str(self.direction_concordance),
            "is_concordant": self.is_concordant,
            "cross_study_status": str(self.cross_study_status),
            "is_synthesis_justified": self.is_synthesis_justified,
            "pooled_log2_fold_change": self.pooled_log2_fold_change,
            "pooled_standard_error": self.pooled_standard_error,
        }


@dataclass(frozen=True, slots=True)
class CrossStudyComparisonDataset:
    """Immutable collection of GeneCrossStudyComparisonRecord objects for a cross-study comparison."""

    comparison_id: str
    assessment: StudyComparabilityAssessment
    gene_comparisons: tuple[GeneCrossStudyComparisonRecord, ...]
    created_at: str = ""
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def total_genes_evaluated(self) -> int:
        return len(self.gene_comparisons)

    @property
    def shared_genes_count(self) -> int:
        return sum(
            1 for r in self.gene_comparisons
            if r.study1_estimate_status != EstimateStatus.INVALID_STATISTIC
            and r.study2_estimate_status != EstimateStatus.INVALID_STATISTIC
        )

    @property
    def study1_only_genes_count(self) -> int:
        return sum(
            1 for r in self.gene_comparisons
            if r.study1_estimate_status != EstimateStatus.INVALID_STATISTIC
            and r.study2_estimate_status == EstimateStatus.INVALID_STATISTIC
        )

    @property
    def study2_only_genes_count(self) -> int:
        return sum(
            1 for r in self.gene_comparisons
            if r.study1_estimate_status == EstimateStatus.INVALID_STATISTIC
            and r.study2_estimate_status != EstimateStatus.INVALID_STATISTIC
        )

    @property
    def concordant_up_count(self) -> int:
        return sum(1 for r in self.gene_comparisons if r.direction_concordance == DirectionConcordance.CONCORDANT_UP)

    @property
    def concordant_down_count(self) -> int:
        return sum(1 for r in self.gene_comparisons if r.direction_concordance == DirectionConcordance.CONCORDANT_DOWN)

    @property
    def discordant_count(self) -> int:
        return sum(1 for r in self.gene_comparisons if r.direction_concordance == DirectionConcordance.DISCORDANT)

    @property
    def no_direction_count(self) -> int:
        return sum(1 for r in self.gene_comparisons if r.direction_concordance == DirectionConcordance.NO_DIRECTION)

    @property
    def undetermined_count(self) -> int:
        return sum(1 for r in self.gene_comparisons if r.direction_concordance == DirectionConcordance.UNDETERMINED)

    @property
    def concordance_rate(self) -> float:
        directional_total = self.concordant_up_count + self.concordant_down_count + self.discordant_count
        if directional_total == 0:
            return 0.0
        return (self.concordant_up_count + self.concordant_down_count) / directional_total

