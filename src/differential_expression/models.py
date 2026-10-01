"""Immutable deterministic data models for Bulk Differential Expression v1.

Establishes gene-level treatment-response profiles from Step 7A contrasts
using official Bioconductor DESeq2 negative-binomial GLM with locked size factors.
"""

from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from enum import StrEnum
from typing import Any

from src.response_builder.models import ContrastExposure, ContrastStatus, Finding, Severity, Status, ordered_findings, status_for


@dataclass(frozen=True, slots=True)
class GeneDifferentialExpressionResult:
    """Differential expression statistics for a single canonical gene."""

    canonical_gene_id: str
    approved_symbol: str | None
    base_mean: float
    log2_fold_change: float | None
    lfc_standard_error: float | None
    wald_statistic: float | None
    p_value: float | None
    adjusted_p_value_bh: float | None
    result_status: str = "OK"  # "OK", "OUTLIER", "ZERO_COUNTS", "FIT_ERROR"
    result_notes: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "canonical_gene_id": self.canonical_gene_id,
            "approved_symbol": self.approved_symbol,
            "base_mean": self.base_mean,
            "log2_fold_change": self.log2_fold_change,
            "lfc_standard_error": self.lfc_standard_error,
            "wald_statistic": self.wald_statistic,
            "p_value": self.p_value,
            "adjusted_p_value_bh": self.adjusted_p_value_bh,
            "result_status": self.result_status,
            "result_notes": self.result_notes,
        }


@dataclass(frozen=True, slots=True)
class DifferentialExpressionContrastResult:
    """Complete differential-expression response profile for one Step 7A contrast."""

    contrast_id: str
    study_id: str
    experiment_id: str
    cohort_id: str
    relationship_id: str

    treatment_condition_id: str
    matched_control_condition_ids: tuple[str, ...]

    treatment_sample_ids: tuple[str, ...]
    control_sample_ids: tuple[str, ...]

    treatment_biological_replicate_ids: tuple[str, ...]
    control_biological_replicate_ids: tuple[str, ...]

    treatment_replicate_count: int
    control_replicate_count: int

    # All condition-planned exposures preserved without flattening or silent truncation
    treatment_exposures: tuple[ContrastExposure, ...] = ()
    sample_level_exposure_deviations: tuple[dict[str, Any], ...] = ()

    # Backward-compatible exposure context
    agent_name: str | None = None
    agent_identifier: str | None = None
    vehicle: str | None = None
    concentration_or_dose: float | str | None = None
    concentration_or_dose_unit: str | None = None

    exposure_start_time_or_stage: str | None = None
    developmental_age_or_stage_at_exposure: str | None = None
    exposure_duration: float | str | None = None
    exposure_duration_unit: str | None = None
    washout_or_recovery_duration: float | str | None = None
    washout_or_recovery_duration_unit: str | None = None

    # Separately preserved developmental collection ages and normalized representations
    treatment_collection_age_or_stage: str | None = None
    control_collection_age_or_stage: str | None = None
    treatment_collection_age_normalized: str | None = None
    control_collection_age_normalized: str | None = None

    biological_source_ids: tuple[str, ...] = ()
    organoid_context_ids: tuple[str, ...] = ()

    # Statistical model design
    design_formula: str = ""
    design_rationale: str = ""
    design_matrix_rank: int | None = None
    residual_degrees_of_freedom: int | None = None

    # Size factors used (locked from Step 6C)
    size_factors_used: dict[str, float] = field(default_factory=dict)

    # Gene-wise differential expression results
    gene_results: tuple[GeneDifferentialExpressionResult, ...] = ()

    # Contrast evaluation outcome
    status: ContrastStatus = ContrastStatus.BLOCKED
    findings: tuple[Finding, ...] = ()

    # R runtime & provenance
    r_environment_info: dict[str, Any] = field(default_factory=dict)
    input_hashes: dict[str, str] = field(default_factory=dict)
    configuration_version: str = "1.0.0"
    de_version: str = "1.0.0"

    @property
    def is_inferentially_eligible(self) -> bool:
        """True if the contrast DE result successfully passed inference qualification."""
        return (
            self.status == ContrastStatus.ELIGIBLE
            and self.treatment_replicate_count >= 2
            and self.control_replicate_count >= 2
            and len(self.gene_results) > 0
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "contrast_id": self.contrast_id,
            "study_id": self.study_id,
            "experiment_id": self.experiment_id,
            "cohort_id": self.cohort_id,
            "relationship_id": self.relationship_id,
            "treatment_condition_id": self.treatment_condition_id,
            "matched_control_condition_ids": list(self.matched_control_condition_ids),
            "treatment_sample_ids": list(self.treatment_sample_ids),
            "control_sample_ids": list(self.control_sample_ids),
            "treatment_biological_replicate_ids": list(self.treatment_biological_replicate_ids),
            "control_biological_replicate_ids": list(self.control_biological_replicate_ids),
            "treatment_replicate_count": self.treatment_replicate_count,
            "control_replicate_count": self.control_replicate_count,
            "treatment_exposures": [e.to_dict() for e in self.treatment_exposures],
            "sample_level_exposure_deviations": list(self.sample_level_exposure_deviations),
            "agent_name": self.agent_name,
            "agent_identifier": self.agent_identifier,
            "vehicle": self.vehicle,
            "concentration_or_dose": self.concentration_or_dose,
            "concentration_or_dose_unit": self.concentration_or_dose_unit,
            "exposure_start_time_or_stage": self.exposure_start_time_or_stage,
            "developmental_age_or_stage_at_exposure": self.developmental_age_or_stage_at_exposure,
            "exposure_duration": self.exposure_duration,
            "exposure_duration_unit": self.exposure_duration_unit,
            "washout_or_recovery_duration": self.washout_or_recovery_duration,
            "washout_or_recovery_duration_unit": self.washout_or_recovery_duration_unit,
            "treatment_collection_age_or_stage": self.treatment_collection_age_or_stage,
            "control_collection_age_or_stage": self.control_collection_age_or_stage,
            "treatment_collection_age_normalized": self.treatment_collection_age_normalized,
            "control_collection_age_normalized": self.control_collection_age_normalized,
            "biological_source_ids": list(self.biological_source_ids),
            "organoid_context_ids": list(self.organoid_context_ids),
            "design_formula": self.design_formula,
            "design_rationale": self.design_rationale,
            "design_matrix_rank": self.design_matrix_rank,
            "residual_degrees_of_freedom": self.residual_degrees_of_freedom,
            "size_factors_used": dict(self.size_factors_used),
            "gene_results": [g.to_dict() for g in self.gene_results],
            "status": str(self.status),
            "findings": [f.to_dict() for f in self.findings],
            "r_environment_info": self.r_environment_info,
            "input_hashes": self.input_hashes,
            "configuration_version": self.configuration_version,
            "de_version": self.de_version,
        }


@dataclass(frozen=True, slots=True)
class DifferentialExpressionDataset:
    """Dataset-level container aggregating all differential expression results."""

    status: Status
    findings: tuple[Finding, ...]
    contrast_results: tuple[DifferentialExpressionContrastResult, ...]
    source_asset_id: str
    reference_identity: str
    builder_version: str
    de_version: str = "1.0.0"
    execution_provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def eligible_results(self) -> tuple[DifferentialExpressionContrastResult, ...]:
        """Contrasts with successful and eligible differential expression results."""
        return tuple(r for r in self.contrast_results if r.status == ContrastStatus.ELIGIBLE)

    @property
    def needs_review_results(self) -> tuple[DifferentialExpressionContrastResult, ...]:
        """Contrasts requiring expert scientific review."""
        return tuple(r for r in self.contrast_results if r.status == ContrastStatus.NEEDS_REVIEW)

    @property
    def blocked_results(self) -> tuple[DifferentialExpressionContrastResult, ...]:
        """Contrasts blocked from differential expression inference."""
        return tuple(r for r in self.contrast_results if r.status == ContrastStatus.BLOCKED)

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": str(self.status),
            "findings": [f.to_dict() for f in self.findings],
            "contrast_results": [r.to_dict() for r in self.contrast_results],
            "number_of_contrasts": len(self.contrast_results),
            "number_of_eligible": len(self.eligible_results),
            "number_of_needs_review": len(self.needs_review_results),
            "number_of_blocked": len(self.blocked_results),
            "source_asset_id": self.source_asset_id,
            "reference_identity": self.reference_identity,
            "builder_version": self.builder_version,
            "de_version": self.de_version,
            "execution_provenance": self.execution_provenance,
        }
