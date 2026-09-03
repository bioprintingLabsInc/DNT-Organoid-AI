"""Non-mutating registration of canonical input-data assets."""

from typing import Any

from .models import RegisteredAsset


def reported_value(value_record: Any) -> Any:
    """Read, but never alter, a reported canonical value_record."""
    if not isinstance(value_record, dict) or value_record.get("value_status") != "reported":
        return None
    normalized = value_record.get("normalized_value")
    return normalized if normalized is not None else value_record.get("original_value")


def register_assets(metadata: dict[str, Any]) -> tuple[RegisteredAsset, ...]:
    """Create immutable asset registrations from canonical metadata."""
    assays = {
        assay.get("assay_id"): reported_value(assay.get("rna_seq_modality"))
        for assay in metadata.get("sequencing_assays", [])
        if isinstance(assay, dict)
    }
    registrations = []
    for asset in metadata.get("input_data_assets", []):
        if not isinstance(asset, dict):
            continue
        assay_ids = tuple(asset.get("assay_ids", ())) if isinstance(asset.get("assay_ids"), list) else ()
        modalities = tuple(sorted({assays[item] for item in assay_ids if assays.get(item)}))
        registrations.append(RegisteredAsset(
            asset_id=str(asset.get("asset_id", "")),
            source_reference=str(reported_value(asset.get("source_file_identifier_or_reference")) or ""),
            declared_format=reported_value(asset.get("input_data_format")),
            assay_ids=assay_ids,
            modalities=modalities,
        ))
    return tuple(registrations)

