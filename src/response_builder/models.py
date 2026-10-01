"""Immutable deterministic data models for Treatment-versus-Matched-Control Contrast Builder v1.

Establishes validated molecular-response contrasts from explicit metadata relationships,
preserving biological context, developmental age, replication invariants, and study/experiment boundaries.
"""

from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from enum import StrEnum
from typing import Any


class ContrastStatus(StrEnum):
    """Scientific evaluation status for an individual molecular-response contrast."""

    ELIGIBLE = "ELIGIBLE"
    NEEDS_REVIEW = "NEEDS_REVIEW"
    BLOCKED = "BLOCKED"


class Severity(StrEnum):
    """Severity levels for contrast audit findings."""

    ERROR = "ERROR"
    REVIEW = "REVIEW"
    WARNING = "WARNING"
    INFO = "INFO"


class Status(StrEnum):
    """Aggregate dataset-level status for the contrast collection."""

    FAIL = "FAIL"
    NEEDS_REVIEW = "NEEDS_REVIEW"
    PASS_WITH_WARNINGS = "PASS_WITH_WARNINGS"
    PASS = "PASS"


@dataclass(frozen=True, slots=True)
class Finding:
    """Structured audit finding capturing validation outcomes and defect details."""

    severity: Severity
    rule_id: str
    entity_type: str
    entity_id: str | None
    path: str
    message: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


_SEVERITY_ORDER = {
    Severity.ERROR: 0,
    Severity.REVIEW: 1,
    Severity.WARNING: 2,
    Severity.INFO: 3,
}


def ordered_findings(findings: Iterable[Finding]) -> tuple[Finding, ...]:
    """Sort findings deterministically by severity precedence, then rule, entity, path, message."""
    return tuple(
        sorted(
            findings,
            key=lambda item: (
                _SEVERITY_ORDER[item.severity],
                item.rule_id,
                item.entity_type,
                item.entity_id or "",
                item.path,
                item.message,
            ),
        )
    )


def status_for(findings: Iterable[Finding]) -> Status:
    """Derive dataset status based on severity precedence."""
    severities = {item.severity for item in findings}
    if Severity.ERROR in severities:
        return Status.FAIL
    if Severity.REVIEW in severities:
        return Status.NEEDS_REVIEW
    if Severity.WARNING in severities:
        return Status.PASS_WITH_WARNINGS
    return Status.PASS


@dataclass(frozen=True, slots=True)
class ContrastExposure:
    """Condition-planned exposure record preserved deterministically."""

    exposure_id: str
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

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class MolecularResponseContrast:
    """Scientifically validated contrast comparing one treatment condition to its matched controls."""

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

    # All condition-planned exposures preserved without flattening or silent truncation
    treatment_exposures: tuple[ContrastExposure, ...] = ()

    # Sample-level exposure deviations explicitly distinguished from condition-planned records
    sample_level_exposure_deviations: tuple[dict[str, Any], ...] = ()

    # Backward-compatible singular fields (only populated when len(treatment_exposures) == 1)
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

    sample_comparison_exception_ids: tuple[str, ...] = ()

    status: ContrastStatus = ContrastStatus.BLOCKED
    findings: tuple[Finding, ...] = ()

    @property
    def is_inferentially_eligible(self) -> bool:
        """True if the contrast is ELIGIBLE and has at least 2 distinct biological replicates in each arm."""
        return (
            self.status == ContrastStatus.ELIGIBLE
            and self.treatment_replicate_count >= 2
            and self.control_replicate_count >= 2
        )

    @property
    def treatment_replicate_count(self) -> int:
        """Count of distinct biological replicates in the treatment arm."""
        return len(set(self.treatment_biological_replicate_ids))

    @property
    def control_replicate_count(self) -> int:
        """Count of distinct biological replicates in the control arm."""
        return len(set(self.control_biological_replicate_ids))

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
            "is_inferentially_eligible": self.is_inferentially_eligible,
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
            "sample_comparison_exception_ids": list(self.sample_comparison_exception_ids),
            "status": str(self.status),
            "findings": [f.to_dict() for f in self.findings],
        }


@dataclass(frozen=True, slots=True)
class ResponseContrastDataset:
    """Dataset-level container aggregating all evaluated molecular-response contrasts."""

    status: Status
    findings: tuple[Finding, ...]
    contrasts: tuple[MolecularResponseContrast, ...]
    treatment_conditions_encountered: tuple[str, ...]
    source_asset_id: str
    reference_identity: str
    builder_version: str
    metadata_hash: str

    @property
    def eligible_contrasts(self) -> tuple[MolecularResponseContrast, ...]:
        """Contrasts that are ELIGIBLE and inferentially qualified."""
        return tuple(c for c in self.contrasts if c.is_inferentially_eligible)

    @property
    def needs_review_contrasts(self) -> tuple[MolecularResponseContrast, ...]:
        """Contrasts that require expert review."""
        return tuple(c for c in self.contrasts if c.status == ContrastStatus.NEEDS_REVIEW)

    @property
    def blocked_contrasts(self) -> tuple[MolecularResponseContrast, ...]:
        """Contrasts that are scientifically blocked from downstream analysis."""
        return tuple(c for c in self.contrasts if c.status == ContrastStatus.BLOCKED)

    @property
    def treatment_condition_dispositions(self) -> dict[str, ContrastStatus]:
        """Map each encountered treatment condition to its aggregated disposition status."""
        dispositions: dict[str, ContrastStatus] = {}
        for cid in self.treatment_conditions_encountered:
            c_contrasts = [c for c in self.contrasts if c.treatment_condition_id == cid]
            if not c_contrasts:
                dispositions[cid] = ContrastStatus.BLOCKED
            elif any(c.status == ContrastStatus.BLOCKED for c in c_contrasts):
                dispositions[cid] = ContrastStatus.BLOCKED
            elif any(c.status == ContrastStatus.NEEDS_REVIEW for c in c_contrasts):
                dispositions[cid] = ContrastStatus.NEEDS_REVIEW
            else:
                dispositions[cid] = ContrastStatus.ELIGIBLE
        return dispositions

    @property
    def eligible_dispositions_count(self) -> int:
        return sum(1 for s in self.treatment_condition_dispositions.values() if s == ContrastStatus.ELIGIBLE)

    @property
    def needs_review_dispositions_count(self) -> int:
        return sum(1 for s in self.treatment_condition_dispositions.values() if s == ContrastStatus.NEEDS_REVIEW)

    @property
    def blocked_dispositions_count(self) -> int:
        return sum(1 for s in self.treatment_condition_dispositions.values() if s == ContrastStatus.BLOCKED)

    def verify_treatment_condition_conservation(self) -> bool:
        """Verify that every encountered treatment condition has an explicit auditable disposition."""
        total_encountered = len(self.treatment_conditions_encountered)
        total_dispositions = (
            self.eligible_dispositions_count
            + self.needs_review_dispositions_count
            + self.blocked_dispositions_count
        )
        return total_encountered == total_dispositions

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": str(self.status),
            "findings": [f.to_dict() for f in self.findings],
            "treatment_conditions_encountered": list(self.treatment_conditions_encountered),
            "total_treatment_conditions_encountered_count": len(self.treatment_conditions_encountered),
            "total_contrasts_count": len(self.contrasts),
            "eligible_contrasts_count": len(self.eligible_contrasts),
            "needs_review_contrasts_count": len(self.needs_review_contrasts),
            "blocked_contrasts_count": len(self.blocked_contrasts),
            "eligible_dispositions_count": self.eligible_dispositions_count,
            "needs_review_dispositions_count": self.needs_review_dispositions_count,
            "blocked_dispositions_count": self.blocked_dispositions_count,
            "treatment_condition_conservation_verified": self.verify_treatment_condition_conservation(),
            "contrasts": [c.to_dict() for c in self.contrasts],
            "source_asset_id": self.source_asset_id,
            "reference_identity": self.reference_identity,
            "builder_version": self.builder_version,
            "metadata_hash": self.metadata_hash,
        }
