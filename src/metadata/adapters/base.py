"""Base interfaces and data structures for source metadata adapters."""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class MappingFinding:
    severity: str  # ERROR, WARNING, INFO
    rule: str
    field_path: str
    message: str


@dataclass(frozen=True, slots=True)
class RawSampleRecord:
    """Raw sample metadata extracted directly from source without loss."""

    sample_id: str
    attributes: dict[str, Any]
    source_file: str
    source_locator: str  # e.g., 'row_2', 'GSM12345', 'Sample 1'
    source_key_provenance: dict[str, str] = field(default_factory=dict)  # attr_name -> original_column_name


@dataclass(frozen=True, slots=True)
class StudyContext:
    """Study-level metadata extracted from source documents or manifests."""

    study_id: str
    source_type: str = "public"  # 'public' or 'in_house'
    study_title: str | None = None
    public_accessions: tuple[dict[str, str], ...] = ()
    experiment_id: str = "exp_01"
    experiment_title: str | None = None
    experimental_batch_id: str | None = None
    attributes: dict[str, Any] = field(default_factory=dict)
    findings: tuple[MappingFinding, ...] = ()


class MetadataSourceAdapter(ABC):
    """Abstract base class for all source metadata layout adapters."""

    @abstractmethod
    def can_read(self, source: Any) -> bool:
        """Return True if this adapter can interpret the given source file/content."""
        raise NotImplementedError

    @abstractmethod
    def parse_study(
        self,
        source: Any,
        **kwargs: Any,
    ) -> tuple[StudyContext, list[RawSampleRecord], list[MappingFinding]]:
        """Parse source metadata into StudyContext and RawSampleRecords with full provenance."""
        raise NotImplementedError
