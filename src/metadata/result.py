"""Structured, deterministic metadata validation results."""

from dataclasses import asdict, dataclass
from enum import StrEnum
from typing import Any


class Severity(StrEnum):
    ERROR = "ERROR"
    WARNING = "WARNING"


class ValidationStatus(StrEnum):
    VALID = "VALID"
    VALID_WITH_WARNINGS = "VALID_WITH_WARNINGS"
    INVALID = "INVALID"


@dataclass(frozen=True, slots=True)
class Finding:
    severity: Severity
    entity_type: str
    entity_identifier: str | None
    path: str
    rule: str
    message: str

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable audit representation."""
        return asdict(self)


@dataclass(frozen=True, slots=True)
class ValidationResult:
    status: ValidationStatus
    findings: tuple[Finding, ...]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable result."""
        return {"status": self.status, "findings": [item.to_dict() for item in self.findings]}

