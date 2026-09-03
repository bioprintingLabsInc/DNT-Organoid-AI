"""Declared/detected consistency checks and modality-aware route selection."""

from typing import Any

from src.ingestion import AssetResult, AssetStatus, Confidence, Route, inspect_asset, register_assets
from src.ingestion.models import AssetFinding, Detection, RegisteredAsset
from src.metadata import MetadataValidator, ValidationStatus

from .detector import detect_format


def route_metadata_assets(metadata: dict[str, Any]) -> tuple[AssetResult, ...]:
    """Validate canonical metadata, inspect each local asset, and route safely."""
    metadata_result = MetadataValidator().validate(metadata)
    registrations = register_assets(metadata)
    return tuple(_route(asset, metadata_result.status == ValidationStatus.INVALID) for asset in registrations)


def _route(asset: RegisteredAsset, metadata_invalid: bool) -> AssetResult:
    inspection = inspect_asset(asset.local_path)
    detection = detect_format(asset, inspection)
    findings = list(detection.findings)
    if metadata_invalid:
        findings.append(AssetFinding("ERROR", "metadata_invalid", "metadata", "Canonical metadata is invalid; routing cannot be approved."))
    mismatch = bool(asset.declared_format and detection.detected_format and asset.declared_format != detection.detected_format)
    if mismatch:
        findings.append(AssetFinding("WARNING", "declared_detected_mismatch", "input_data_format", f"Declared format '{asset.declared_format}' conflicts with detected '{detection.detected_format}'."))
    route, route_findings = _select_route(asset, detection)
    if mismatch:
        route = Route.MANUAL_REVIEW
    findings.extend(route_findings)
    findings = tuple(sorted(findings, key=lambda item: (item.severity, item.rule, item.path, item.message)))
    status = AssetStatus.INVALID if any(item.severity == "ERROR" for item in findings) else AssetStatus.NEEDS_REVIEW if findings or route == Route.MANUAL_REVIEW else AssetStatus.READY
    return AssetResult(asset.asset_id, asset.source_reference, asset.declared_format, detection.detected_format, detection.confidence, detection.evidence, asset.modalities, route, status, findings)


def _select_route(asset: RegisteredAsset, detection: Detection) -> tuple[Route, tuple[AssetFinding, ...]]:
    detected = detection.detected_format
    modalities = set(asset.modalities)
    if asset.declared_format == "UMI_counts" and "bulk_rna_seq" in modalities:
        return _manual("Declared UMI_counts are incompatible with bulk_rna_seq routing", "umi_modality")
    if asset.declared_format == "raw_counts" and modalities and modalities != {"bulk_rna_seq"}:
        return _manual("Declared raw_counts cannot be routed as bulk counts for this modality", "raw_counts_modality")
    if detected in {"FASTQ", "SRA"}:
        return Route.RAW_READ_PIPELINE, ()
    if detected == "raw_counts":
        if modalities == {"bulk_rna_seq"}:
            return Route.BULK_RAW_COUNT_PIPELINE, ()
        return _manual("raw_counts require unambiguous bulk_rna_seq assay metadata", "raw_counts_modality")
    if detected == "UMI_counts":
        if modalities and modalities <= {"scRNA_seq", "snRNA_seq"}:
            return Route.SINGLE_CELL_COUNT_PIPELINE, ()
        return _manual("UMI_counts require compatible scRNA_seq or snRNA_seq metadata", "umi_modality")
    if detected == "H5AD":
        if modalities and modalities <= {"scRNA_seq", "snRNA_seq"}:
            return Route.SINGLE_CELL_OBJECT_PIPELINE, ()
        return _manual("H5AD requires compatible scRNA_seq or snRNA_seq metadata", "h5ad_modality")
    if detected == "processed_expression_matrix":
        return Route.PROCESSED_EXPRESSION_PIPELINE, ()
    if detected == "RDS":
        return _manual("RDS object contents are not interpreted automatically", "rds_manual_review")
    return Route.MANUAL_REVIEW, ()


def _manual(message: str, rule: str) -> tuple[Route, tuple[AssetFinding, ...]]:
    return Route.MANUAL_REVIEW, (AssetFinding("WARNING", rule, "route", message),)
