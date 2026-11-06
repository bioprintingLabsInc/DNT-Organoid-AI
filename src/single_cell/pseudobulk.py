"""Strict raw-integer pseudobulk aggregation enforcing donor independence and canonical bulk count output."""

from collections import defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy import sparse

from src.qc.matrix import CountMatrix
from src.qc.result import Finding, Severity, Status
from src.single_cell.models import (
    CellMetadata,
    InvalidRawCountMatrixException,
    PseudobulkConfig,
    PseudobulkDataset,
    PseudobulkSample,
    SingleCellMatrix,
    SingleCellQCResult,
)


def aggregate_pseudobulk(
    matrix: SingleCellMatrix,
    qc_result: SingleCellQCResult,
    config: PseudobulkConfig | None = None,
    output_tsv_path: Path | str | None = None,
    study_id: str = "pseudobulk_study",
    dataset_id: str = "pseudobulk_dataset",
    metadata_template: dict[str, Any] | None = None,
) -> PseudobulkDataset:
    """Aggregate quality-passed single-cell/nucleus counts into canonical pseudobulk integer matrices.

    Strict scientific safeguards:
    1. Aggregation is by (donor_id, cell_type, condition_id).
    2. Sums raw integer UMI counts: count(g, p) = sum_{c in cells(p)} X(g, c).
       Never averages, never normalizes before aggregation.
    3. Replicate independence: all cells from a donor in a given cell_type and condition
       form exactly ONE pseudobulk sample.
    4. Minimum cell threshold: groups with fewer than config.min_cells_per_sample are excluded.
    5. Emits canonical CountMatrix ready for Route.BULK_RAW_COUNT_PIPELINE (Stage 6A/6B/6C/7B).
    """
    if config is None:
        config = PseudobulkConfig()

    findings: list[Finding] = []

    # 1. Group quality-passed cells by (donor_id, cell_type, condition_id)
    groups: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    cell_metadata_by_id = {c.cell_id: c for c in qc_result.cell_metrics}

    if qc_result.status == Status.FAIL or qc_result.passed_cells == 0:
        empty_matrix = CountMatrix(matrix.gene_ids, (), ())
        return PseudobulkDataset(
            dataset_id=dataset_id,
            cell_types=(),
            count_matrix=empty_matrix,
            per_cell_type_matrices={},
            pseudobulk_samples=(),
            metadata={},
            provenance={"aggregation_rule": "sum_raw_integers", "status": "FAILED_QC"},
            findings=qc_result.findings,
            status=Status.FAIL,
        )

    for cell in qc_result.cell_metrics:
        if cell.is_pass_qc:
            if cell.donor_id.startswith("unverified_") or cell.condition_id.startswith("unverified_"):
                continue
            key = (cell.donor_id, cell.cell_type, cell.condition_id)
            groups[key].append(cell.cell_id)

    # 2. Filter groups by minimum cell threshold
    eligible_groups: dict[tuple[str, str, str], list[str]] = {}
    for key, cell_ids in groups.items():
        donor, ct, cond = key
        if len(cell_ids) >= config.min_cells_per_sample:
            eligible_groups[key] = cell_ids
        else:
            findings.append(
                Finding(
                    Severity.REVIEW,
                    "insufficient_cells_for_pseudobulk",
                    "pseudobulk_group",
                    f"{donor}_{ct}_{cond}",
                    "cell_count",
                    f"Pseudobulk group donor='{donor}', cell_type='{ct}', condition='{cond}' has only {len(cell_ids)} cells (minimum {config.min_cells_per_sample} required); excluded from pseudobulk.",
                )
            )

    if not eligible_groups:
        findings.append(
            Finding(
                Severity.ERROR,
                "no_eligible_pseudobulk_samples",
                "pseudobulk",
                dataset_id,
                "sample_generation",
                f"No donor/cell_type/condition groups met the minimum cell count threshold ({config.min_cells_per_sample}).",
            )
        )
        empty_matrix = CountMatrix(matrix.gene_ids, (), ())
        return PseudobulkDataset(
            dataset_id=dataset_id,
            cell_types=(),
            count_matrix=empty_matrix,
            per_cell_type_matrices={},
            pseudobulk_samples=(),
            metadata={},
            provenance={"aggregation_rule": "sum_raw_integers", "status": "FAILED"},
            findings=tuple(findings),
            status=Status.FAIL,
        )

    # 3. Deterministic sample ordering and ID generation
    sorted_keys = sorted(eligible_groups.keys(), key=lambda k: (k[1], k[2], k[0]))  # cell_type, condition, donor
    sample_ids: list[str] = []
    pseudobulk_samples: list[PseudobulkSample] = []
    sample_columns: list[tuple[int, ...]] = []

    # Map cell_id to column index in single cell matrix
    cell_indices = {c: i for i, c in enumerate(matrix.cell_ids)}
    raw_mat = matrix.raw_matrix

    # Track per-cell-type sample subsets
    cell_types_seen: set[str] = set()
    per_cell_type_sample_ids: dict[str, list[str]] = defaultdict(list)
    per_cell_type_columns: dict[str, list[tuple[int, ...]]] = defaultdict(list)

    for key in sorted_keys:
        donor, ct, cond = key
        cell_types_seen.add(ct)
        # Clean sample ID: e.g. donor1_neurons_treated
        clean_donor = donor.replace(" ", "_").replace("/", "_")
        clean_ct = ct.replace(" ", "_").replace("/", "_")
        clean_cond = cond.replace(" ", "_").replace("/", "_")
        s_id = f"pb_{clean_donor}_{clean_ct}_{clean_cond}"

        # Avoid collision if multiple groups produce identical clean IDs
        suffix = 1
        base_s_id = s_id
        while s_id in sample_ids:
            suffix += 1
            s_id = f"{base_s_id}_{suffix}"

        sample_ids.append(s_id)
        per_cell_type_sample_ids[ct].append(s_id)

        # Sum raw integer counts across group cells
        cell_ids = eligible_groups[key]
        c_idxs = [cell_indices[c] for c in cell_ids]

        if sparse.issparse(raw_mat):
            sub_sparse = raw_mat[:, c_idxs]
            # Exact sum across columns (axis 1)
            summed_arr = np.asarray(sub_sparse.sum(axis=1)).flatten()
        else:
            summed_arr = np.asarray(raw_mat[:, c_idxs].sum(axis=1)).flatten()

        for v in summed_arr:
            f_val = float(v)
            if not f_val.is_integer() or f_val < 0:
                raise InvalidRawCountMatrixException(
                    f"Pseudobulk sum encountered non-integer or negative value: {v}. "
                    "Raw integer UMI counts are required; rounding or truncation is prohibited."
                )
        int_column = tuple(int(v) for v in summed_arr)
        sample_columns.append(int_column)
        per_cell_type_columns[ct].append(int_column)

        total_umi = sum(int_column)
        pseudobulk_samples.append(
            PseudobulkSample(
                sample_id=s_id,
                donor_id=donor,
                cell_type=ct,
                condition_id=cond,
                cell_count=len(cell_ids),
                total_umi_count=total_umi,
            )
        )

    # 4. Construct canonical CountMatrix (all combined)
    count_matrix = CountMatrix(
        gene_ids=matrix.gene_ids,
        sample_ids=tuple(sample_ids),
        columns=tuple(sample_columns),
    )

    # Per cell-type CountMatrices
    per_cell_type_matrices: dict[str, CountMatrix] = {}
    for ct, ct_sids in per_cell_type_sample_ids.items():
        per_cell_type_matrices[ct] = CountMatrix(
            gene_ids=matrix.gene_ids,
            sample_ids=tuple(ct_sids),
            columns=tuple(per_cell_type_columns[ct]),
        )

    # 5. Write to TSV file if requested
    emitted_path = None
    if output_tsv_path:
        out_p = Path(output_tsv_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        with out_p.open("w", encoding="utf-8") as f:
            f.write("gene\t" + "\t".join(sample_ids) + "\n")
            for g_idx, g_id in enumerate(matrix.gene_ids):
                row_vals = [str(sample_columns[s_idx][g_idx]) for s_idx in range(len(sample_ids))]
                f.write(f"{g_id}\t" + "\t".join(row_vals) + "\n")
        emitted_path = out_p

    # 6. Build canonical study metadata aligning with pseudobulk samples
    conditions_present = sorted(list({s.condition_id for s in pseudobulk_samples}))
    donors_present = sorted(list({s.donor_id for s in pseudobulk_samples}))

    def _reported(val: Any) -> dict[str, Any]:
        return {
            "original_value": val,
            "normalized_value": val,
            "value_status": "reported",
            "assertion_provenance_ids": [],
        }

    def _unknown() -> dict[str, Any]:
        return {
            "original_value": None,
            "normalized_value": None,
            "value_status": "unknown",
            "assertion_provenance_ids": [],
        }

    def _resolve_record(source: dict[str, Any] | None, field_name: str) -> dict[str, Any]:
        if not source or field_name not in source:
            return _unknown()
        val = source[field_name]
        if isinstance(val, dict) and "value_status" in val:
            return val
        if val is not None and str(val).strip() != "" and str(val).strip().lower() not in ("unknown", "none", "nan", "null"):
            return _reported(val)
        return _unknown()

    template = metadata_template or {}

    # Conditions
    template_conditions = {
        c.get("condition_id"): c
        for c in template.get("conditions", [])
        if isinstance(c, dict) and c.get("condition_id")
    }

    conditions_list = []
    for cond in conditions_present:
        t_cond = template_conditions.get(cond)
        if t_cond:
            c_label = _resolve_record(t_cond, "condition_label")
            if c_label.get("value_status") != "reported":
                c_label = _reported(cond)
            tc_status = _resolve_record(t_cond, "treatment_control_status")
            c_type = _resolve_record(t_cond, "control_type")
        else:
            c_label = _reported(cond)
            tc_status = _unknown()
            c_type = _unknown()

        c_entry: dict[str, Any] = {
            "condition_id": cond,
            "experiment_id": "exp_sc_pseudobulk",
            "condition_label": c_label,
            "treatment_control_status": tc_status,
            "metadata_validation_status": "not_validated",
        }
        if c_type.get("value_status") == "reported":
            c_entry["control_type"] = c_type
        conditions_list.append(c_entry)

    # Biological Sources
    template_sources = {
        s.get("biological_source_id"): s
        for s in template.get("biological_sources", [])
        if isinstance(s, dict) and s.get("biological_source_id")
    }
    for s in template.get("biological_sources", []):
        if isinstance(s, dict):
            d_id = s.get("donor_identifier")
            if isinstance(d_id, dict):
                d_id = d_id.get("original_value") or d_id.get("normalized_value")
            if d_id:
                template_sources[str(d_id)] = s

    sources_list = []
    for d in donors_present:
        src = template_sources.get(f"source_{d}") or template_sources.get(d)
        species_rec = _resolve_record(src, "species")
        if species_rec.get("value_status") != "reported" and template.get("species"):
            species_rec = _resolve_record(template, "species")

        mat_rec = _resolve_record(src, "source_material_type")
        if mat_rec.get("value_status") != "reported" and template.get("source_material_type"):
            mat_rec = _resolve_record(template, "source_material_type")

        cl_rec = _resolve_record(src, "cell_line_or_source_identifier")
        if cl_rec.get("value_status") != "reported":
            cl_rec = _reported(d)

        donor_rec = _resolve_record(src, "donor_identifier")
        if donor_rec.get("value_status") != "reported":
            donor_rec = _reported(d)

        sources_list.append({
            "biological_source_id": f"source_{d}",
            "species": species_rec,
            "source_material_type": mat_rec,
            "cell_line_or_source_identifier": cl_rec,
            "donor_identifier": donor_rec,
            "metadata_validation_status": "not_validated",
        })

    # Organoid Contexts
    template_contexts = [c for c in template.get("organoid_contexts", []) if isinstance(c, dict)]
    if template_contexts:
        contexts_list = template_contexts
    else:
        org_type_rec = _resolve_record(template, "organoid_type")
        region_rec = _resolve_record(template, "brain_region_or_model_identity")
        contexts_list = [{
            "organoid_context_id": "organoid_context_1",
            "organoid_type": org_type_rec,
            "brain_region_or_model_identity": region_rec,
            "metadata_validation_status": "not_validated",
        }]

    # Samples
    template_samples = {
        s.get("sample_id"): s
        for s in template.get("samples", [])
        if isinstance(s, dict) and s.get("sample_id")
    }

    samples_list = []
    for s in pseudobulk_samples:
        m_s = template_samples.get(s.sample_id)
        age_rec = _resolve_record(m_s, "developmental_age_or_stage_at_collection")
        if age_rec.get("value_status") != "reported" and template.get("developmental_age_or_stage_at_collection"):
            age_rec = _resolve_record(template, "developmental_age_or_stage_at_collection")

        dur_rec = _resolve_record(m_s, "collection_time_relative_to_exposure")
        if dur_rec.get("value_status") != "reported" and template.get("collection_time_relative_to_exposure"):
            dur_rec = _resolve_record(template, "collection_time_relative_to_exposure")

        samples_list.append({
            "sample_id": s.sample_id,
            "experiment_id": "exp_sc_pseudobulk",
            "condition_id": s.condition_id,
            "biological_source_id": f"source_{s.donor_id}",
            "organoid_context_id": contexts_list[0]["organoid_context_id"],
            "biological_replicate_id": _reported(f"{s.donor_id}_{s.cell_type}"),
            "developmental_age_or_stage_at_collection": age_rec,
            "collection_time_relative_to_exposure": dur_rec,
            "metadata_validation_status": "not_validated",
        })

    # Treatment Control Relationships
    # NEVER infer control status or fabricate relationships! Only preserve verified relationships from template.
    template_rels = [
        r for r in template.get("treatment_control_relationships", [])
        if isinstance(r, dict)
    ]
    relationships_list = template_rels

    canonical_metadata: dict[str, Any] = {
        "studies": [{
            "study_id": study_id,
            "source_type": "public",
            "metadata_validation_status": "not_validated",
        }],
        "experiments": [{
            "experiment_id": "exp_sc_pseudobulk",
            "study_id": study_id,
            "metadata_validation_status": "not_validated",
        }],
        "conditions": conditions_list,
        "biological_sources": sources_list,
        "organoid_contexts": contexts_list,
        "samples": samples_list,
        "exposures": [],
        "treatment_control_relationships": relationships_list,
        "sample_comparison_exceptions": [],
        "sequencing_assays": [{
            "assay_id": "assay_pseudobulk_bulk_counts",
            "sample_ids": sample_ids,
            "rna_seq_modality": _reported("bulk_rna_seq"),
            "metadata_validation_status": "not_validated",
        }],
        "input_data_assets": [{
            "asset_id": "asset_pseudobulk_raw_counts",
            "assay_ids": ["assay_pseudobulk_bulk_counts"],
            "input_data_format": _reported("raw_counts"),
            "gene_identifier_type": _reported("Ensembl Gene ID"),
            "source_file_identifier_or_reference": _reported(str(emitted_path) if emitted_path else "pseudobulk_counts.tsv"),
            "availability": "available",
            "metadata_validation_status": "not_validated",
        }],
        "provenance": [],
        "dnt_reference_evidence": [],
    }

    # Add field-level provenance
    for entity_name, records in canonical_metadata.items():
        if isinstance(records, list):
            for r in records:
                if isinstance(r, dict):
                    ent_id = r.get("sample_id") or r.get("condition_id") or r.get("study_id") or r.get("assay_id") or r.get("asset_id") or r.get("experiment_id") or r.get("biological_source_id") or r.get("organoid_context_id") or r.get("relationship_id")
                    for k, v in r.items():
                        if isinstance(v, dict) and "value_status" in v:
                            prov_id = f"prov_pb_{len(canonical_metadata['provenance']) + 1}"
                            v["assertion_provenance_ids"] = [prov_id]
                            canonical_metadata["provenance"].append({
                                "provenance_id": prov_id,
                                "entity_type": entity_name,
                                "entity_id": str(ent_id),
                                "field_path": k,
                                "source_kind": "pipeline_generated",
                                "source_reference": "single_cell_pseudobulk_aggregator",
                            })

    provenance_record = {
        "pipeline_phase": "Phase U4: scRNA-seq/snRNA-seq Pseudobulk Aggregation",
        "aggregation_rule": "strict_sum_raw_integer_umi_counts",
        "donor_independence_enforced": True,
        "min_cells_per_sample": config.min_cells_per_sample,
        "total_single_cells": matrix.number_of_cells,
        "qc_passed_single_cells": qc_result.passed_cells,
        "pseudobulk_sample_count": len(pseudobulk_samples),
        "cell_types": sorted(list(cell_types_seen)),
        "conditions": conditions_present,
        "donors": donors_present,
    }

    return PseudobulkDataset(
        dataset_id=dataset_id,
        cell_types=tuple(sorted(list(cell_types_seen))),
        count_matrix=count_matrix,
        per_cell_type_matrices=per_cell_type_matrices,
        pseudobulk_samples=tuple(pseudobulk_samples),
        metadata=canonical_metadata,
        provenance=provenance_record,
        findings=tuple(findings),
        status=Status.PASS,
        emitted_matrix_path=emitted_path,
    )
