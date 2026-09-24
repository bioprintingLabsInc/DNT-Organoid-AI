"""Immutable deterministic data models for Gene Harmonization v1."""

from collections.abc import Iterable
from dataclasses import asdict, dataclass
from enum import StrEnum
from typing import Any

from src.qc.matrix import CountMatrix


class GeneMappingStatus(StrEnum):
    """Explicit mapping outcome for each source gene."""

    UNIQUELY_MAPPED = "UNIQUELY_MAPPED"
    AMBIGUOUS = "AMBIGUOUS"
    UNMAPPED = "UNMAPPED"
    INVALID = "INVALID"


class GeneIdentifierType(StrEnum):
    """Supported gene identifier types."""

    ENSEMBL_GENE_ID = "ensembl_gene_id"
    GENE_SYMBOL = "gene_symbol"
    ENTREZ_GENE_ID = "entrez_gene_id"
    HGNC_ID = "hgnc_id"
    CUSTOM = "custom"


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
class HarmonizedGene:
    """Harmonization record for a single input gene row."""

    source_index: int
    original_gene_id: str
    identifier_type: str
    canonical_gene_id: str | None
    approved_symbol: str | None
    mapping_status: GeneMappingStatus
    mapping_reason: str | None
    candidate_canonical_ids: tuple[str, ...]
    reference_id: str
    reference_version: str
    source_asset_id: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class CanonicalCollision:
    """Group of multiple source identifiers mapping to the same canonical gene."""

    canonical_gene_id: str
    original_gene_ids: tuple[str, ...]
    source_indices: tuple[int, ...]
    approved_symbol: str | None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class HarmonizationSummary:
    """Dataset-level summary metrics of the harmonization stage."""

    total_genes: int
    uniquely_mapped_genes: int
    unmapped_genes: int
    ambiguous_genes: int
    invalid_genes: int
    canonical_collision_genes: int
    colliding_source_gene_count: int
    percentage_uniquely_mapped: float
    percentage_unresolved: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class QCReviewDisposition:
    """Structured, auditable reviewed disposition for upstream NEEDS_REVIEW QC results."""

    decision: str
    reason: str
    reviewer: str | None = None
    review_date: str | None = None
    affected_qc_status: str = "NEEDS_REVIEW"
    affected_asset_id: str | None = None
    provenance: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        if not self.decision or not self.decision.strip():
            raise ValueError("QCReviewDisposition must specify a non-empty 'decision'.")
        if not self.reason or not self.reason.strip():
            raise ValueError("QCReviewDisposition must specify a non-empty 'reason' or scientific rationale.")

    @classmethod
    def from_any(
        cls,
        disposition: Any,
        asset_id: str | None = None,
        qc_status: str = "NEEDS_REVIEW",
    ) -> "QCReviewDisposition | None":
        """Coerce structured disposition from QCReviewDisposition, dict, or auditable string."""
        if disposition is None:
            return None
        if isinstance(disposition, cls):
            return disposition
        if isinstance(disposition, dict):
            decision = disposition.get("decision") or disposition.get("status") or "AUTHORIZED"
            reason = (
                disposition.get("reason")
                or disposition.get("rationale")
                or disposition.get("justification")
                or ""
            )
            reviewer = disposition.get("reviewer") or disposition.get("review_authority")
            review_date = (
                disposition.get("review_date")
                or disposition.get("review_time")
                or disposition.get("date")
            )
            affected_status = disposition.get("affected_qc_status", qc_status)
            affected_asset = disposition.get("affected_asset_id", asset_id)
            prov = disposition.get("provenance")
            if not str(decision).strip() or not str(reason).strip():
                return None
            return cls(
                decision=str(decision).strip(),
                reason=str(reason).strip(),
                reviewer=str(reviewer).strip() if reviewer else None,
                review_date=str(review_date).strip() if review_date else None,
                affected_qc_status=str(affected_status).strip(),
                affected_asset_id=str(affected_asset).strip() if affected_asset else None,
                provenance=prov if isinstance(prov, dict) else None,
            )
        if isinstance(disposition, str) and disposition.strip():
            return cls(
                decision="AUTHORIZED",
                reason=disposition.strip(),
                affected_qc_status=qc_status,
                affected_asset_id=asset_id,
            )
        return None

    def to_dict(self) -> dict[str, Any]:
        return {
            "decision": self.decision,
            "reason": self.reason,
            "reviewer": self.reviewer,
            "review_date": self.review_date,
            "affected_qc_status": self.affected_qc_status,
            "affected_asset_id": self.affected_asset_id,
            "provenance": self.provenance,
        }


@dataclass(frozen=True, slots=True)
class HarmonizedDataset:
    """Complete, immutable harmonization result for a bulk expression asset."""

    status: Status
    findings: tuple[Finding, ...]
    summary: HarmonizationSummary
    genes: tuple[HarmonizedGene, ...]
    collisions: tuple[CanonicalCollision, ...]
    matrix: CountMatrix | None
    source_asset_id: str
    assay_id: str
    sample_ids: tuple[str, ...]
    reference_identity: str
    reference_version: str
    reference_checksum: str
    config_version: str
    config_checksum: str
    harmonization_version: str
    upstream_qc_status: str | None = None
    qc_review_disposition: QCReviewDisposition | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": str(self.status),
            "upstream_qc_status": self.upstream_qc_status,
            "qc_review_disposition": self.qc_review_disposition.to_dict() if self.qc_review_disposition else None,
            "findings": [f.to_dict() for f in self.findings],
            "summary": self.summary.to_dict(),
            "genes": [g.to_dict() for g in self.genes],
            "collisions": [c.to_dict() for c in self.collisions],
            "number_of_genes": len(self.genes),
            "number_of_samples": len(self.sample_ids),
            "source_asset_id": self.source_asset_id,
            "assay_id": self.assay_id,
            "sample_ids": list(self.sample_ids),
            "reference_identity": self.reference_identity,
            "reference_version": self.reference_version,
            "reference_checksum": self.reference_checksum,
            "config_version": self.config_version,
            "config_checksum": self.config_checksum,
            "harmonization_version": self.harmonization_version,
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
    """Derive dataset harmonization status based on severity precedence."""
    severities = {item.severity for item in findings}
    if Severity.ERROR in severities:
        return Status.FAIL
    if Severity.REVIEW in severities:
        return Status.NEEDS_REVIEW
    if Severity.WARNING in severities:
        return Status.PASS_WITH_WARNINGS
    return Status.PASS
