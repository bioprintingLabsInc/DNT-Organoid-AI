"""Immutable deterministic QC result models."""
from dataclasses import dataclass
from enum import StrEnum

class Severity(StrEnum): ERROR="ERROR"; REVIEW="REVIEW"; WARNING="WARNING"; INFO="INFO"
class Status(StrEnum): FAIL="FAIL"; NEEDS_REVIEW="NEEDS_REVIEW"; PASS_WITH_WARNINGS="PASS_WITH_WARNINGS"; PASS="PASS"

@dataclass(frozen=True, slots=True)
class Finding:
    severity: Severity; rule_id: str; entity_type: str; entity_id: str | None; path: str; message: str

@dataclass(frozen=True, slots=True)
class SampleResult:
    sample_id: str; condition_id: str | None; biological_replicate_id: str | None; technical_replicate_id: str | None
    library_size: int; detected_genes: int; zero_count_genes: int; zero_fraction: float
    minimum: int; median: float; maximum: int; quantiles: tuple[float, float, float]
    dataset_library_ratio: float | None; condition_library_ratio: float | None; findings: tuple[Finding, ...]

@dataclass(frozen=True, slots=True)
class PairCorrelation:
    sample_a: str; sample_b: str; replicate_a: str; replicate_b: str; value: float | None; reason: str | None

@dataclass(frozen=True, slots=True)
class ConditionResult:
    condition_id: str; sample_count: int; biological_replicate_count: int; technical_replicates: tuple[tuple[str, tuple[str, ...]], ...]
    correlation_method: str; correlation_implementation: str; pairs: tuple[PairCorrelation, ...]
    evaluable_pairs: int; correlation_summary: tuple[float, float, float] | None; unevaluable_reasons: tuple[str, ...]; findings: tuple[Finding, ...]

@dataclass(frozen=True, slots=True)
class DatasetResult:
    status: Status; findings: tuple[Finding, ...]; number_of_genes: int; number_of_samples: int
    source_asset_id: str; assay_id: str; study_ids: tuple[str, ...]; experiment_ids: tuple[str, ...]
    ruleset_version: str; ruleset_checksum: str; calculation_version: str
    samples: tuple[SampleResult, ...]; conditions: tuple[ConditionResult, ...]

def ordered(findings):
    order={Severity.ERROR:0,Severity.REVIEW:1,Severity.WARNING:2,Severity.INFO:3}
    return tuple(sorted(findings,key=lambda x:(order[x.severity],x.rule_id,x.entity_type,x.entity_id or "",x.path,x.message)))

def status_for(findings):
    levels={f.severity for f in findings}
    return Status.FAIL if Severity.ERROR in levels else Status.NEEDS_REVIEW if Severity.REVIEW in levels else Status.PASS_WITH_WARNINGS if Severity.WARNING in levels else Status.PASS
