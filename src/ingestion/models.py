"""Immutable records shared by ingestion and routing."""

from dataclasses import asdict, dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any


class AssetStatus(StrEnum):
    """READY is consistent and routable; NEEDS_REVIEW is ambiguous; INVALID has an error."""

    READY = "READY"
    NEEDS_REVIEW = "NEEDS_REVIEW"
    INVALID = "INVALID"


class Route(StrEnum):
    """Future processing destination; no route performs processing in Step 5."""

    RAW_READ_PIPELINE = "RAW_READ_PIPELINE"
    BULK_RAW_COUNT_PIPELINE = "BULK_RAW_COUNT_PIPELINE"
    SINGLE_CELL_COUNT_PIPELINE = "SINGLE_CELL_COUNT_PIPELINE"
    SINGLE_CELL_OBJECT_PIPELINE = "SINGLE_CELL_OBJECT_PIPELINE"
    PROCESSED_EXPRESSION_PIPELINE = "PROCESSED_EXPRESSION_PIPELINE"
    MANUAL_REVIEW = "MANUAL_REVIEW"


class Confidence(StrEnum):
    """Strength of the bounded evidence supporting detected_format."""

    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"
    NONE = "NONE"


@dataclass(frozen=True, slots=True)
class AssetFinding:
    severity: str
    rule: str
    path: str
    message: str


@dataclass(frozen=True, slots=True)
class RegisteredAsset:
    asset_id: str
    source_reference: str
    declared_format: str | None
    assay_ids: tuple[str, ...]
    modalities: tuple[str, ...]

    @property
    def local_path(self) -> Path:
        return Path(self.source_reference)


@dataclass(frozen=True, slots=True)
class Inspection:
    exists: bool
    readable: bool
    size: int | None
    prefix: bytes
    suffixes: tuple[str, ...]
    findings: tuple[AssetFinding, ...] = ()


@dataclass(frozen=True, slots=True)
class Detection:
    detected_format: str | None
    confidence: Confidence
    evidence: tuple[str, ...]
    findings: tuple[AssetFinding, ...] = ()


@dataclass(frozen=True, slots=True)
class AssetResult:
    asset_id: str
    source_reference: str
    declared_format: str | None
    detected_format: str | None
    confidence: Confidence
    evidence: tuple[str, ...]
    modalities: tuple[str, ...]
    route: Route
    status: AssetStatus
    findings: tuple[AssetFinding, ...]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable audit record."""
        return asdict(self)
