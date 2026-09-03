"""Loading and access helpers for the accepted metadata contract."""

from pathlib import Path
from typing import Any

import yaml


DEFAULT_SCHEMA_PATH = Path(__file__).resolve().parents[2] / "config" / "metadata_schema.yaml"


def load_schema(path: str | Path = DEFAULT_SCHEMA_PATH) -> dict[str, Any]:
    """Load the metadata contract and reject a non-mapping document."""
    with Path(path).open(encoding="utf-8") as stream:
        schema = yaml.safe_load(stream)
    if not isinstance(schema, dict):
        raise ValueError("Metadata schema must be a YAML mapping")
    return schema


def identifier_field(entity_name: str, definition: dict[str, Any]) -> str | None:
    """Return the entity's declared identifier field."""
    conventional = {
        "studies": "study_id",
        "experiments": "experiment_id",
        "conditions": "condition_id",
        "biological_sources": "biological_source_id",
        "organoid_contexts": "organoid_context_id",
        "samples": "sample_id",
        "exposures": "exposure_id",
        "treatment_control_relationships": "relationship_id",
        "sample_comparison_exceptions": "exception_id",
        "sequencing_assays": "assay_id",
        "input_data_assets": "asset_id",
        "provenance": "provenance_id",
        "dnt_reference_evidence": "evidence_id",
    }
    return conventional.get(entity_name)
