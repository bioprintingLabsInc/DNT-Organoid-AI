"""Reusable canonical metadata fixtures for validator tests."""

from copy import deepcopy
from typing import Any


def reported(value: Any) -> dict[str, Any]:
    return {
        "original_value": value,
        "normalized_value": value,
        "value_status": "reported",
        "assertion_provenance_ids": [],
    }


def missing(status: str = "missing_not_reported") -> dict[str, Any]:
    return {
        "original_value": None,
        "normalized_value": None,
        "value_status": status,
        "assertion_provenance_ids": [],
    }


def canonical_bulk_study(source_type: str = "public") -> dict[str, list[dict[str, Any]]]:
    """Return a DMSO control and compound treatment, each with three replicates."""
    document: dict[str, list[dict[str, Any]]] = {
        "studies": [{
            "study_id": "study_1",
            "source_type": source_type,
            "metadata_validation_status": "not_validated",
        }],
        "experiments": [{
            "experiment_id": "experiment_1",
            "study_id": "study_1",
            "metadata_validation_status": "not_validated",
        }],
        "conditions": [
            {
                "condition_id": "control",
                "experiment_id": "experiment_1",
                "condition_label": reported("DMSO control"),
                "treatment_control_status": reported("control"),
                "control_type": reported("vehicle"),
                "metadata_validation_status": "not_validated",
            },
            {
                "condition_id": "treatment_10um_24h_d45",
                "experiment_id": "experiment_1",
                "condition_label": reported("Compound A, 10 µM, 24 h, day 45"),
                "treatment_control_status": reported("treatment"),
                "metadata_validation_status": "not_validated",
            },
        ],
        "biological_sources": [{
            "biological_source_id": "line_1",
            "species": reported("Homo sapiens"),
            "source_material_type": reported("iPSC"),
            "cell_line_or_source_identifier": reported("line_1"),
            "metadata_validation_status": "not_validated",
        }],
        "organoid_contexts": [{
            "organoid_context_id": "cortical_1",
            "organoid_type": reported("cerebral organoid"),
            "brain_region_or_model_identity": reported("dorsal forebrain"),
            "metadata_validation_status": "not_validated",
        }],
        "samples": [],
        "exposures": [
            {
                "exposure_id": "control_exposure",
                "condition_id": "control",
                "exposure_scope": "condition_planned",
                "vehicle": reported("DMSO"),
                "developmental_age_or_stage_at_exposure": reported("day 45"),
                "metadata_validation_status": "not_validated",
            },
            {
                "exposure_id": "treatment_exposure",
                "condition_id": "treatment_10um_24h_d45",
                "exposure_scope": "condition_planned",
                "agent_name": reported("Compound A"),
                "agent_identifier": reported("example:compound-a"),
                "vehicle": reported("DMSO"),
                "concentration_or_dose": reported(10),
                "concentration_or_dose_unit": reported("µM"),
                "exposure_start_time_or_stage": reported("day 45"),
                "developmental_age_or_stage_at_exposure": reported("day 45"),
                "exposure_duration": reported(24),
                "exposure_duration_unit": reported("hour"),
                "metadata_validation_status": "not_validated",
            },
        ],
        "treatment_control_relationships": [{
            "relationship_id": "comparison_1",
            "treatment_condition_ids": ["treatment_10um_24h_d45"],
            "matched_control_condition_ids": ["control"],
            "matching_basis": reported("same line, stage, duration, and vehicle"),
            "metadata_validation_status": "not_validated",
        }],
        "sample_comparison_exceptions": [],
        "sequencing_assays": [],
        "input_data_assets": [],
        "provenance": [],
        "dnt_reference_evidence": [],
    }
    if source_type == "public":
        document["studies"][0]["public_accessions"] = [
            {"namespace": "GSE", "identifier": "GSE000001", "url": "https://example.org/GSE000001"}
        ]
    for condition_id, prefix in (("control", "c"), ("treatment_10um_24h_d45", "t")):
        for replicate in range(1, 4):
            document["samples"].append({
                "sample_id": f"sample_{prefix}{replicate}",
                "experiment_id": "experiment_1",
                "condition_id": condition_id,
                "biological_source_id": "line_1",
                "organoid_context_id": "cortical_1",
                "biological_replicate_id": reported(f"replicate_{replicate}"),
                "developmental_age_or_stage_at_collection": reported("day 46"),
                "collection_time_relative_to_exposure": reported("24 hours"),
                "metadata_validation_status": "not_validated",
            })
    sample_ids = [sample["sample_id"] for sample in document["samples"]]
    document["sequencing_assays"].append({
        "assay_id": "bulk_assay",
        "sample_ids": sample_ids,
        "rna_seq_modality": reported("bulk_rna_seq"),
        "metadata_validation_status": "not_validated",
    })
    document["input_data_assets"].append({
        "asset_id": "counts_asset",
        "assay_ids": ["bulk_assay"],
        "input_data_format": reported("raw_counts"),
        "gene_identifier_type": reported("Ensembl Gene ID"),
        "source_file_identifier_or_reference": reported("counts.tsv"),
        "availability": "available",
        "metadata_validation_status": "not_validated",
    })
    return with_provenance(document)


def with_provenance(document: dict[str, list[dict[str, Any]]]) -> dict[str, list[dict[str, Any]]]:
    """Attach deterministic field-level provenance to every reported value_record."""
    result = deepcopy(document)
    result["provenance"] = []
    identifiers = {
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
        "dnt_reference_evidence": "evidence_id",
    }

    def attach(value: Any, entity: str, entity_id: str, field: str) -> None:
        if isinstance(value, dict):
            if "value_status" in value:
                provenance_id = f"prov_{len(result['provenance']) + 1}"
                value["assertion_provenance_ids"] = [provenance_id]
                result["provenance"].append({
                    "provenance_id": provenance_id,
                    "entity_type": entity,
                    "entity_id": entity_id,
                    "field_path": field,
                    "source_kind": "source_file",
                    "source_reference": "fixture_metadata.yaml",
                })
            else:
                for nested in value.values():
                    attach(nested, entity, entity_id, field)
        elif isinstance(value, list):
            for nested in value:
                attach(nested, entity, entity_id, field)

    for entity, id_field in identifiers.items():
        for record in result[entity]:
            for field, value in record.items():
                attach(value, entity, record[id_field], field)
    return result

