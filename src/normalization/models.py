"""Immutable deterministic data models for Bulk RNA-seq Normalization v1."""

from collections.abc import Iterable
from dataclasses import asdict, dataclass
from enum import StrEnum
from typing import Any

from src.qc.matrix import CountMatrix


class ResponseEligibility(StrEnum):
    """Downstream eligibility for planned treatment-versus-matched-control response construction (Step 7)."""

    ELIGIBLE_WITH_MATCHED_CONTROL = "ELIGIBLE_WITH_MATCHED_CONTROL"
    INELIGIBLE_MISSING_MATCHED_CONTROL = "INELIGIBLE_MISSING_MATCHED_CONTROL"
    CONTROL_BASELINE = "CONTROL_BASELINE"
    UNASSIGNED_OR_AMBIGUOUS = "UNASSIGNED_OR_AMBIGUOUS"


class NormalizationMethod(StrEnum):
    """Supported size-factor estimation methods."""

    RATIO = "ratio"
    POSCOUNTS = "poscounts"


class Severity(StrEnum):
    ERROR = "ERROR"
    REVIEW = "REVIEW"
    WARNING = "WARNING"
    INFO = "INFO"


class Status(StrEnum):
    FAIL = "FAIL"
    NEEDS_REVIEW = "NEEDS_REVIEW"
    PASS_WITH_WARNINGS = "PASS_WITH_WARNINGS"
    PASS = "PASS"


@dataclass(frozen=True, slots=True)
class Finding:
    severity: Severity
    rule_id: str
    entity_type: str
    entity_id: str | None
    path: str
    message: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class NormalizedCountMatrix:
    """Secondary floating-point normalized count matrix for inspection, technical auditing, and diagnostics."""

    gene_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    columns: tuple[tuple[float, ...], ...]

    def get_count(self, gene_id: str, sample_id: str) -> float | None:
        try:
            gene_idx = self.gene_ids.index(gene_id)
            sample_idx = self.sample_ids.index(sample_id)
            return self.columns[sample_idx][gene_idx]
        except ValueError:
            return None

    def to_dict(self) -> dict[str, Any]:
        return {
            "gene_ids": list(self.gene_ids),
            "sample_ids": list(self.sample_ids),
            "number_of_genes": len(self.gene_ids),
            "number_of_samples": len(self.sample_ids),
        }


@dataclass(frozen=True, slots=True)
class SampleSizeFactor:
    """Estimated size factor and downstream response eligibility for a single sample."""

    sample_id: str
    size_factor: float
    raw_library_size: int
    normalized_library_size: float
    cohort_id: str
    response_eligibility: ResponseEligibility
    matched_control_sample_ids: tuple[str, ...]
    treatment_control_status: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "size_factor": self.size_factor,
            "raw_library_size": self.raw_library_size,
            "normalized_library_size": self.normalized_library_size,
            "cohort_id": self.cohort_id,
            "response_eligibility": str(self.response_eligibility),
            "matched_control_sample_ids": list(self.matched_control_sample_ids),
            "treatment_control_status": self.treatment_control_status,
        }


@dataclass(frozen=True, slots=True)
class NormalizationCohort:
    """Scientifically validated group of samples for joint size-factor estimation."""

    cohort_id: str
    study_id: str
    experiment_id: str
    sample_ids: tuple[str, ...]
    assay_ids: tuple[str, ...]
    rna_seq_modality: str
    library_strategy: str | None
    sequencing_platform: str | None
    strandedness: str | None
    read_layout: str | None
    treatment_sample_ids: tuple[str, ...]
    control_sample_ids: tuple[str, ...] = ()
    unassigned_sample_ids: tuple[str, ...] = ()
    treatment_control_relationship_ids: tuple[str, ...] = ()
    metadata_hash: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "cohort_id": self.cohort_id,
            "study_id": self.study_id,
            "experiment_id": self.experiment_id,
            "sample_ids": list(self.sample_ids),
            "assay_ids": list(self.assay_ids),
            "rna_seq_modality": self.rna_seq_modality,
            "library_strategy": self.library_strategy,
            "sequencing_platform": self.sequencing_platform,
            "strandedness": self.strandedness,
            "read_layout": self.read_layout,
            "treatment_sample_ids": list(self.treatment_sample_ids),
            "control_sample_ids": list(self.control_sample_ids),
            "unassigned_sample_ids": list(self.unassigned_sample_ids),
            "treatment_control_relationship_ids": list(self.treatment_control_relationship_ids),
            "metadata_hash": self.metadata_hash,
        }


@dataclass(frozen=True, slots=True)
class CohortNormalizationMetrics:
    """Descriptive technical QC metrics for a normalized cohort without invented thresholds."""

    total_features_count: int
    eligible_features_count: int
    positive_geometric_mean_features_count: int
    percentage_features_evaluable: float
    size_factor_min: float
    size_factor_max: float
    size_factor_median: float
    size_factor_mean: float
    size_factor_std: float
    size_factor_ratio_max_min: float
    number_of_samples: int
    number_of_treatment_samples: int
    number_of_control_samples: int
    treatment_control_ratio: float | None
    zero_count_fraction: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class NormalizationReviewDisposition:
    """Structured, auditable reviewed disposition for NEEDS_REVIEW normalization results or poscounts authorization."""

    decision: str
    rationale: str
    reviewer: str
    review_date: str
    affected_cohort_id: str | None = None
    provenance: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        if not self.decision or not self.decision.strip():
            raise ValueError("NormalizationReviewDisposition must specify a non-empty 'decision'.")
        if not self.rationale or not self.rationale.strip():
            raise ValueError("NormalizationReviewDisposition must specify a non-empty 'rationale'.")
        if not self.reviewer or not self.reviewer.strip():
            raise ValueError("NormalizationReviewDisposition must specify a non-empty 'reviewer'.")
        if not self.review_date or not self.review_date.strip():
            raise ValueError("NormalizationReviewDisposition must specify a non-empty 'review_date'.")

    @classmethod
    def from_any(cls, disposition: Any, cohort_id: str | None = None) -> "NormalizationReviewDisposition | None":
        if disposition is None:
            return None
        if isinstance(disposition, cls):
            return disposition
        if isinstance(disposition, dict):
            decision = disposition.get("decision") or disposition.get("status") or "AUTHORIZED"
            rationale = disposition.get("rationale") or disposition.get("reason") or disposition.get("justification") or ""
            reviewer = disposition.get("reviewer") or disposition.get("review_authority") or ""
            review_date = disposition.get("review_date") or disposition.get("date") or ""
            aff_cohort = disposition.get("affected_cohort_id", cohort_id)
            prov = disposition.get("provenance")
            if not str(decision).strip() or not str(rationale).strip() or not str(reviewer).strip() or not str(review_date).strip():
                return None
            return cls(
                decision=str(decision).strip(),
                rationale=str(rationale).strip(),
                reviewer=str(reviewer).strip(),
                review_date=str(review_date).strip(),
                affected_cohort_id=str(aff_cohort).strip() if aff_cohort else None,
                provenance=prov if isinstance(prov, dict) else None,
            )
        return None

    def to_dict(self) -> dict[str, Any]:
        return {
            "decision": self.decision,
            "rationale": self.rationale,
            "reviewer": self.reviewer,
            "review_date": self.review_date,
            "affected_cohort_id": self.affected_cohort_id,
            "provenance": self.provenance,
        }


@dataclass(frozen=True, slots=True)
class CollisionDisposition:
    """Structured, auditable reviewed disposition for canonical collision gate authorization."""

    decision: str
    rationale: str
    reviewer: str
    review_date: str
    resolved_collision_ids: tuple[str, ...] = ()
    provenance: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        if not self.decision or not self.decision.strip():
            raise ValueError("CollisionDisposition must specify a non-empty 'decision'.")
        if not self.rationale or not self.rationale.strip():
            raise ValueError("CollisionDisposition must specify a non-empty 'rationale'.")
        if not self.reviewer or not self.reviewer.strip():
            raise ValueError("CollisionDisposition must specify a non-empty 'reviewer'.")
        if not self.review_date or not self.review_date.strip():
            raise ValueError("CollisionDisposition must specify a non-empty 'review_date'.")

    @classmethod
    def from_any(cls, disposition: Any) -> "CollisionDisposition | None":
        if disposition is None:
            return None
        if isinstance(disposition, cls):
            return disposition
        if isinstance(disposition, dict):
            decision = disposition.get("decision") or disposition.get("status") or "AUTHORIZED"
            rationale = disposition.get("rationale") or disposition.get("reason") or disposition.get("justification") or ""
            reviewer = disposition.get("reviewer") or disposition.get("review_authority") or ""
            review_date = disposition.get("review_date") or disposition.get("date") or ""
            resolved = disposition.get("resolved_collision_ids") or ()
            prov = disposition.get("provenance")
            if not str(decision).strip() or not str(rationale).strip() or not str(reviewer).strip() or not str(review_date).strip():
                return None
            return cls(
                decision=str(decision).strip(),
                rationale=str(rationale).strip(),
                reviewer=str(reviewer).strip(),
                review_date=str(review_date).strip(),
                resolved_collision_ids=tuple(str(x) for x in resolved),
                provenance=prov if isinstance(prov, dict) else None,
            )
        return None

    def to_dict(self) -> dict[str, Any]:
        return {
            "decision": self.decision,
            "rationale": self.rationale,
            "reviewer": self.reviewer,
            "review_date": self.review_date,
            "resolved_collision_ids": list(self.resolved_collision_ids),
            "provenance": self.provenance,
        }


@dataclass(frozen=True, slots=True)
class NormalizedCohortResult:
    """Complete, immutable normalization result for a single Normalization Cohort."""

    status: Status
    cohort: NormalizationCohort
    findings: tuple[Finding, ...]
    size_factors: tuple[SampleSizeFactor, ...]
    metrics: CohortNormalizationMetrics | None
    diagnostic_normalized_matrix: NormalizedCountMatrix | None
    raw_count_matrix: CountMatrix
    method: str
    r_environment_info: dict[str, Any]
    execution_provenance: dict[str, Any]
    review_disposition: NormalizationReviewDisposition | None = None
    collision_disposition: CollisionDisposition | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": str(self.status),
            "cohort": self.cohort.to_dict(),
            "findings": [f.to_dict() for f in self.findings],
            "size_factors": [sf.to_dict() for sf in self.size_factors],
            "metrics": self.metrics.to_dict() if self.metrics else None,
            "diagnostic_normalized_matrix": self.diagnostic_normalized_matrix.to_dict() if self.diagnostic_normalized_matrix else None,
            "raw_count_matrix_dimensions": {
                "number_of_genes": len(self.raw_count_matrix.gene_ids),
                "number_of_samples": len(self.raw_count_matrix.sample_ids),
            },
            "method": self.method,
            "r_environment_info": self.r_environment_info,
            "execution_provenance": self.execution_provenance,
            "review_disposition": self.review_disposition.to_dict() if self.review_disposition else None,
            "collision_disposition": self.collision_disposition.to_dict() if self.collision_disposition else None,
        }


@dataclass(frozen=True, slots=True)
class NormalizedDataset:
    """Dataset-level container aggregating all cohort normalization results for an input asset."""

    status: Status
    findings: tuple[Finding, ...]
    cohort_results: tuple[NormalizedCohortResult, ...]
    source_asset_id: str
    reference_identity: str
    config_version: str
    normalization_version: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": str(self.status),
            "findings": [f.to_dict() for f in self.findings],
            "cohort_results": [cr.to_dict() for cr in self.cohort_results],
            "number_of_cohorts": len(self.cohort_results),
            "source_asset_id": self.source_asset_id,
            "reference_identity": self.reference_identity,
            "config_version": self.config_version,
            "normalization_version": self.normalization_version,
        }


_SEVERITY_ORDER = {
    Severity.ERROR: 0,
    Severity.REVIEW: 1,
    Severity.WARNING: 2,
    Severity.INFO: 3,
}


def ordered(findings: Iterable[Finding]) -> tuple[Finding, ...]:
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
    """Derive normalization status based on severity precedence."""
    severities = {item.severity for item in findings}
    if Severity.ERROR in severities:
        return Status.FAIL
    if Severity.REVIEW in severities:
        return Status.NEEDS_REVIEW
    if Severity.WARNING in severities:
        return Status.PASS_WITH_WARNINGS
    return Status.PASS
