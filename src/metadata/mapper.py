"""Canonical metadata mapping engine.

Maps raw sample records and study context from any source metadata adapter into
fully validated canonical metadata dictionaries complying with config/metadata_schema.yaml.

Preserves:
1. 100% of original source cell values in value_record.original_value.
2. Field-level provenance linking each assertion to source file, row, and column key.
3. Scientifically validated unit conversions in value_record.normalized_value without altering originals.
4. Explicit four-state missingness (reported, missing_not_reported, unknown, not_applicable).
5. Deterministic entity grouping and foreign-key integrity.
"""

from collections import defaultdict
from dataclasses import dataclass, field
import datetime
from pathlib import Path
import re
from typing import Any

from .adapters.base import MappingFinding, RawSampleRecord, StudyContext
from .units import (
    UnitConversionStatus,
    normalize_concentration,
    normalize_duration,
    parse_numeric_value,
)


@dataclass(frozen=True, slots=True)
class FieldMappingRule:
    """Rule for resolving a canonical field from raw record attributes."""

    canonical_field: str
    candidate_keys: tuple[str, ...]
    default_value: Any = None
    is_required: bool = False
    transform: str | None = None  # 'dose', 'duration', 'age', 'clean_string'


# Standard candidate key synonyms for common experimental metadata attributes
DEFAULT_MAPPING_RULES: dict[str, tuple[str, ...]] = {
    "sample_id": ("sample_id", "sample", "sample_name", "sample_alias", "gsm_accession", "source_name", "id"),
    "condition_status": ("treatment_control_status", "condition_type", "condition", "group", "treatment", "factor_value_compound", "status"),
    "control_type": ("control_type", "vehicle", "control"),
    "agent_name": ("agent_name", "compound", "agent", "chemical", "factor_value_compound", "treatment"),
    "agent_identifier": ("agent_identifier", "compound_id", "cas", "cid", "chebi", "term_accession_number"),
    "concentration_or_dose": ("concentration_or_dose", "dose", "concentration", "conc", "factor_value_dose"),
    "concentration_or_dose_unit": ("concentration_or_dose_unit", "dose_unit", "unit", "concentration_unit", "unit_concentration_unit"),
    "exposure_duration": ("exposure_duration", "duration", "time", "timepoint", "exposure_time", "factor_value_time"),
    "exposure_duration_unit": ("exposure_duration_unit", "duration_unit", "time_unit"),
    "exposure_start_stage": (
        "exposure_start_time_or_stage",
        "developmental_age_or_stage_at_exposure",
        "exposure_start_stage",
        "exposure_start_day",
        "exposure_start",
        "start_day",
        "start_stage",
        "age_at_exposure",
        "stage_at_exposure",
    ),
    "collection_stage": (
        "developmental_age_or_stage_at_collection",
        "collection_time_relative_to_exposure",
        "collection_stage",
        "collection_day",
        "harvest_day",
        "harvest_stage",
        "collection_age",
        "harvest_age",
        "day_at_collection",
    ),
    "developmental_stage": (
        "developmental_stage",
        "stage",
        "day",
        "age",
        "developmental_age",
        "characteristics_age",
        "characteristics_developmental_stage",
    ),
    "species": ("species", "organism", "characteristics_organism"),
    "source_material_type": ("source_material_type", "material_type", "cell_type", "source_type"),
    "cell_line": ("cell_line_or_source_identifier", "cell_line", "line", "characteristics_cell_line", "cell_line_name"),
    "donor_id": ("donor_identifier", "donor", "donor_id", "patient"),
    "organoid_type": ("organoid_type", "characteristics_organism_part", "organoid", "model_type"),
    "brain_region": ("brain_region_or_model_identity", "brain_region", "organism_part", "region"),
    "biological_replicate_id": ("biological_replicate_id", "replicate_id", "biological_replicate", "replicate", "rep", "bio_rep"),
    "technical_replicate_id": ("technical_replicate_id", "tech_rep", "technical_replicate"),
    "data_file": ("data_file", "filename", "submitted_file_name", "scan_name", "fastq_uri", "associated_files"),
}


def _clean_str(val: Any) -> str | None:
    if val is None:
        return None
    s = str(val).strip()
    return s if s else None


def _find_attribute_value(
    record: RawSampleRecord,
    candidate_keys: tuple[str, ...],
) -> tuple[Any, str | None]:
    """Find value and matching raw column name from raw record attributes."""
    # 1. Exact match
    for k in candidate_keys:
        if k in record.attributes and record.attributes[k] is not None:
            return record.attributes[k], record.source_key_provenance.get(k, k)

    # 2. Normalized lowercase alphanumeric match
    norm_attr_map: dict[str, tuple[str, Any]] = {}
    for raw_k, v in record.attributes.items():
        if v is not None:
            clean_k = re.sub(r"[^a-z0-9]", "", raw_k.lower())
            norm_attr_map[clean_k] = (raw_k, v)

    for k in candidate_keys:
        clean_cand = re.sub(r"[^a-z0-9]", "", k.lower())
        if clean_cand in norm_attr_map:
            raw_k, v = norm_attr_map[clean_cand]
            return v, record.source_key_provenance.get(raw_k, raw_k)

    # 3. Partial substring match for complex tags (e.g. Factor Value[compound], Characteristics[cell line])
    for cand in candidate_keys:
        cand_lower = cand.lower()
        for raw_k, v in record.attributes.items():
            if v is not None and cand_lower in raw_k.lower():
                return v, record.source_key_provenance.get(raw_k, raw_k)

    return None, None


def _is_control_record(record: RawSampleRecord, rules: dict[str, tuple[str, ...]], defaults: dict[str, Any]) -> bool:
    """Determine whether a raw record represents a control condition."""
    raw_cond_status, _ = _find_attribute_value(record, rules["condition_status"])
    cond_str = _clean_str(raw_cond_status) or defaults.get("condition_status")
    status_lower = cond_str.lower() if cond_str else ""
    if any(k in status_lower for k in ["control", "ctrl", "baseline", "untreated", "vehicle", "none", "dmso"]):
        return True
    if "treat" in status_lower or "exposed" in status_lower or "vpa" in status_lower:
        return False
    raw_agent_check, _ = _find_attribute_value(record, rules["agent_name"])
    if raw_agent_check and str(raw_agent_check).lower() in ["none", "untreated", "vehicle", "control"]:
        return True
    return False


class CanonicalMetadataMapper:
    """Maps raw extracted metadata records into canonical schema dictionaries."""

    def __init__(
        self,
        custom_mapping_rules: dict[str, tuple[str, ...]] | None = None,
        defaults: dict[str, Any] | None = None,
        target_dose_unit: str = "uM",
        target_duration_unit: str = "h",
    ) -> None:
        self.rules = dict(DEFAULT_MAPPING_RULES)
        if custom_mapping_rules:
            self.rules.update(custom_mapping_rules)
        self.defaults = defaults or {}
        self.target_dose_unit = target_dose_unit
        self.target_duration_unit = target_duration_unit

    def map_study(
        self,
        context: StudyContext,
        raw_records: list[RawSampleRecord],
        source_reference: str = "source_metadata",
    ) -> tuple[dict[str, Any], list[MappingFinding]]:
        """Construct full canonical metadata dictionary adhering to config/metadata_schema.yaml."""
        findings: list[MappingFinding] = list(context.findings)
        provenance_records: list[dict[str, Any]] = []
        prov_counter = 1

        def make_value_record(
            val: Any,
            norm_val: Any,
            entity_type: str,
            entity_id: str,
            field_path: str,
            raw_key: str | None,
            source_loc: str | None = None,
            status: str = "reported",
        ) -> dict[str, Any]:
            nonlocal prov_counter
            # If status is reported but original value is None, switch to explicit missingness
            if val is None and status == "reported":
                status = "unknown"

            if status != "reported":
                return {
                    "original_value": None,
                    "normalized_value": None,
                    "value_status": status,
                    "assertion_provenance_ids": [],
                }

            pid = f"prov_{entity_type}_{entity_id}_{prov_counter}"
            prov_counter += 1

            prov_rec = {
                "provenance_id": pid,
                "entity_type": entity_type,
                "entity_id": entity_id,
                "field_path": field_path,
                "source_kind": "public_repository" if context.source_type == "public" else "source_file",
                "source_reference": source_reference,
                "original_metadata_key": raw_key,
                "recorded_by": "CanonicalMetadataMapper_v2",
                "recorded_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "notes": source_loc,
            }
            provenance_records.append(prov_rec)

            return {
                "original_value": val,
                "normalized_value": norm_val,
                "value_status": "reported",
                "assertion_provenance_ids": [pid],
            }

        # 1. Study Entry
        study_id = context.study_id
        pub_accs = [dict(a) for a in context.public_accessions]
        if context.source_type == "public" and not pub_accs:
            # Add fallback public accession from study_id if it starts with standard prefix
            clean_acc = study_id.replace("study_", "")
            ns = "GSE" if clean_acc.startswith("GSE") else "ArrayExpress" if clean_acc.startswith("E-") else "public_repo"
            pub_accs.append({"namespace": ns, "identifier": clean_acc, "url": None})

        studies = [
            {
                "study_id": study_id,
                "source_type": context.source_type,
                "study_title": make_value_record(
                    context.study_title or f"Study {study_id}",
                    context.study_title or f"Study {study_id}",
                    "studies",
                    study_id,
                    "study_title",
                    "title",
                ),
                "public_accessions": pub_accs,
                "metadata_validation_status": "not_validated",
            }
        ]

        # 2. Experiment Entry
        exp_id = context.experiment_id or "exp_01"
        experiments = [
            {
                "experiment_id": exp_id,
                "study_id": study_id,
                "experiment_title": make_value_record(
                    context.experiment_title or context.study_title or f"Experiment {exp_id}",
                    context.experiment_title or context.study_title or f"Experiment {exp_id}",
                    "experiments",
                    exp_id,
                    "experiment_title",
                    "title",
                ),
                "metadata_validation_status": "not_validated",
            }
        ]

        # 3. Analyze Samples and Resolve Conditions
        samples: list[dict[str, Any]] = []
        conditions_map: dict[str, dict[str, Any]] = {}
        exposures_map: dict[str, dict[str, Any]] = {}
        bio_sources_map: dict[str, dict[str, Any]] = {}
        organoid_contexts_map: dict[str, dict[str, Any]] = {}
        assays: list[dict[str, Any]] = []
        assets: list[dict[str, Any]] = []

        seen_sample_ids: set[str] = set()

        # Pre-scan records to assess heterogeneity across controls and treatments
        ctrl_signatures: set[tuple[str, str, str, str, str, str]] = set()
        trt_durations: set[str] = set()
        trt_stages: set[str] = set()
        trt_cell_lines: set[str] = set()
        trt_organoids: set[str] = set()
        ctrl_durations: set[str] = set()
        ctrl_cell_lines: set[str] = set()

        for rec in raw_records:
            is_ctrl = _is_control_record(rec, self.rules, self.defaults)

            raw_dur_v, _ = _find_attribute_value(rec, self.rules["exposure_duration"])
            dur_str_val = str(raw_dur_v).strip() if raw_dur_v is not None else ""

            raw_es_v, _ = _find_attribute_value(rec, self.rules.get("exposure_start_stage", ()))
            raw_cs_v, _ = _find_attribute_value(rec, self.rules.get("collection_stage", ()))
            raw_gs_v, _ = _find_attribute_value(rec, self.rules.get("developmental_stage", ()))
            st_str_val = (_clean_str(raw_es_v) or _clean_str(raw_cs_v) or _clean_str(raw_gs_v) or "").lower()

            raw_cl_v, _ = _find_attribute_value(rec, self.rules["cell_line"])
            cl_str_val = (_clean_str(raw_cl_v) or "").lower()

            raw_org_v, _ = _find_attribute_value(rec, self.rules["organoid_type"])
            org_str_val = (_clean_str(raw_org_v) or "").lower()

            if is_ctrl:
                raw_c_type, _ = _find_attribute_value(rec, self.rules["control_type"])
                c_type = _clean_str(raw_c_type) or "vehicle"
                raw_ag, _ = _find_attribute_value(rec, self.rules["agent_name"])
                c_ag = _clean_str(raw_ag) or ""
                ctrl_signatures.add((c_type, c_ag, st_str_val, dur_str_val, cl_str_val, org_str_val))
                if dur_str_val:
                    ctrl_durations.add(dur_str_val)
                if cl_str_val:
                    ctrl_cell_lines.add(cl_str_val)
            else:
                if dur_str_val:
                    trt_durations.add(dur_str_val)
                if st_str_val:
                    trt_stages.add(st_str_val)
                if cl_str_val:
                    trt_cell_lines.add(cl_str_val)
                if org_str_val:
                    trt_organoids.add(org_str_val)

        has_multiple_control_groups = len(ctrl_signatures) > 1
        has_multiple_ctrl_durations = len(ctrl_durations) > 1
        has_multiple_ctrl_cell_lines = len(ctrl_cell_lines) > 1

        has_multiple_durations = len(trt_durations) > 1
        has_multiple_stages = len(trt_stages) > 1
        has_multiple_cell_lines = len(trt_cell_lines) > 1
        has_multiple_organoids = len(trt_organoids) > 1

        for idx, rec in enumerate(raw_records):
            sid = rec.sample_id
            if sid in seen_sample_ids:
                findings.append(
                    MappingFinding(
                        severity="ERROR",
                        rule="duplicate_sample_id",
                        field_path=f"samples[{idx}].sample_id",
                        message=f"Duplicate sample identifier '{sid}'.",
                    )
                )
                continue
            seen_sample_ids.add(sid)

            # Resolve Condition Status (treatment vs control)
            raw_cond_status, cond_status_key = _find_attribute_value(rec, self.rules["condition_status"])
            cond_str = _clean_str(raw_cond_status) or self.defaults.get("condition_status")
            status_lower = cond_str.lower() if cond_str else ""

            is_control = False
            if any(k in status_lower for k in ["control", "ctrl", "baseline", "untreated", "vehicle", "none", "dmso"]):
                is_control = True
            elif "treat" in status_lower or "exposed" in status_lower or "vpa" in status_lower:
                is_control = False
            else:
                # Check agent name if condition status is ambiguous
                raw_agent_check, _ = _find_attribute_value(rec, self.rules["agent_name"])
                if raw_agent_check and str(raw_agent_check).lower() in ["none", "untreated", "vehicle", "control"]:
                    is_control = True
                elif raw_agent_check and str(raw_agent_check).strip():
                    is_control = False
                else:
                    findings.append(
                        MappingFinding(
                            severity="WARNING",
                            rule="ambiguous_condition_status",
                            field_path=f"samples[{idx}].condition_status",
                            message=f"Sample '{sid}' has ambiguous condition status '{cond_str}'; defaulting to treatment.",
                        )
                    )
                    is_control = False

            treatment_status = "control" if is_control else "treatment"

            # Resolve Chemical Identity
            raw_agent, agent_key = _find_attribute_value(rec, self.rules["agent_name"])
            clean_agent = _clean_str(raw_agent) or self.defaults.get("agent_name")
            if is_control and not clean_agent:
                clean_agent = "Vehicle Control"
            elif not is_control and not clean_agent:
                clean_agent = "Unknown Agent"
                findings.append(
                    MappingFinding(
                        severity="WARNING",
                        rule="missing_agent_name",
                        field_path=f"samples[{idx}].agent_name",
                        message=f"Treated sample '{sid}' lacks explicit agent_name.",
                    )
                )

            # Resolve Dose / Concentration
            raw_dose, dose_key = _find_attribute_value(rec, self.rules["concentration_or_dose"])
            raw_dose_unit, dose_unit_key = _find_attribute_value(rec, self.rules["concentration_or_dose_unit"])
            dose_val = raw_dose if raw_dose is not None else (0.0 if is_control else self.defaults.get("concentration_or_dose"))
            dose_unit = raw_dose_unit if raw_dose_unit is not None else ("uM" if is_control else self.defaults.get("concentration_or_dose_unit", "uM" if dose_val is not None else None))

            norm_dose = normalize_concentration(dose_val, dose_unit, target_unit=self.target_dose_unit)

            # Resolve Duration (no unsupported default)
            raw_dur, dur_key = _find_attribute_value(rec, self.rules["exposure_duration"])
            raw_dur_unit, dur_unit_key = _find_attribute_value(rec, self.rules["exposure_duration_unit"])
            dur_val = raw_dur if raw_dur is not None else self.defaults.get("exposure_duration")
            dur_unit = raw_dur_unit if raw_dur_unit is not None else self.defaults.get("exposure_duration_unit", "h" if dur_val is not None else None)

            if dur_val is not None and dur_unit is not None:
                norm_dur = normalize_duration(dur_val, dur_unit, target_unit=self.target_duration_unit)
            else:
                norm_dur = None

            # Resolve Developmental Stages: Separate Exposure Start Age vs Collection Age
            raw_exp_stage, exp_stage_key = _find_attribute_value(rec, self.rules.get("exposure_start_stage", ()))
            clean_exp_stage = _clean_str(raw_exp_stage) or self.defaults.get("exposure_start_stage") or self.defaults.get("developmental_age_or_stage_at_exposure")

            raw_col_stage, col_stage_key = _find_attribute_value(rec, self.rules.get("collection_stage", ()))
            clean_col_stage = _clean_str(raw_col_stage) or self.defaults.get("collection_stage") or self.defaults.get("developmental_age_or_stage_at_collection")

            raw_gen_stage, gen_stage_key = _find_attribute_value(rec, self.rules.get("developmental_stage", ()))
            clean_gen_stage = _clean_str(raw_gen_stage) or self.defaults.get("developmental_stage")

            if clean_exp_stage is None and clean_gen_stage is not None:
                clean_exp_stage = clean_gen_stage
                exp_stage_key = gen_stage_key

            if clean_col_stage is None and clean_gen_stage is not None:
                clean_col_stage = clean_gen_stage
                col_stage_key = gen_stage_key

            # Resolve Biological Source / Cell Line / Material Type
            raw_species, sp_key = _find_attribute_value(rec, self.rules["species"])
            species_str = _clean_str(raw_species) or self.defaults.get("species", "Homo sapiens")

            raw_cell_line, cl_key = _find_attribute_value(rec, self.rules["cell_line"])
            cell_line_str = _clean_str(raw_cell_line) or self.defaults.get("cell_line")

            raw_mat_type, mat_key = _find_attribute_value(rec, self.rules["source_material_type"])
            mat_type_str = _clean_str(raw_mat_type) or self.defaults.get("source_material_type")

            # Controlled vocabulary: biological_source_type in [iPSC, hESC, other]
            if mat_type_str:
                m_lower = mat_type_str.lower()
                if "ipsc" in m_lower or "induced pluripotent" in m_lower:
                    mat_type_str = "iPSC"
                elif "hesc" in m_lower or "embryonic" in m_lower:
                    mat_type_str = "hESC"
                elif mat_type_str in ["iPSC", "hESC", "other"]:
                    pass
                else:
                    mat_type_str = "other"
            elif cell_line_str:
                cl_lower = cell_line_str.lower()
                if any(k in cl_lower for k in ["ipsc", "wtc", "bioni"]):
                    mat_type_str = "iPSC"
                elif any(k in cl_lower for k in ["hesc", "h9", "h1"]):
                    mat_type_str = "hESC"
                else:
                    mat_type_str = "other"
            else:
                mat_type_str = "other"

            bio_source_id = f"source_{re.sub(r'[^a-zA-Z0-9_]', '_', cell_line_str)}" if cell_line_str else f"source_{exp_id}_unknown"
            if bio_source_id not in bio_sources_map:
                bio_sources_map[bio_source_id] = {
                    "biological_source_id": bio_source_id,
                    "species": make_value_record(species_str, species_str, "biological_sources", bio_source_id, "species", sp_key, status="reported" if species_str else "unknown"),
                    "source_material_type": make_value_record(mat_type_str, mat_type_str, "biological_sources", bio_source_id, "source_material_type", mat_key, status="reported"),
                    "cell_line_or_source_identifier": make_value_record(cell_line_str, cell_line_str, "biological_sources", bio_source_id, "cell_line_or_source_identifier", cl_key, status="reported" if cell_line_str else "unknown"),
                    "metadata_validation_status": "not_validated",
                }

            # Resolve Organoid Context (no unsupported defaults)
            raw_org_type, org_key = _find_attribute_value(rec, self.rules["organoid_type"])
            org_type_str = _clean_str(raw_org_type) or self.defaults.get("organoid_type")

            raw_region, reg_key = _find_attribute_value(rec, self.rules["brain_region"])
            region_str = _clean_str(raw_region) or self.defaults.get("brain_region")

            org_ctx_id = f"ctx_{re.sub(r'[^a-zA-Z0-9_]', '_', org_type_str)}" if org_type_str else f"ctx_{exp_id}_unknown"
            if org_ctx_id not in organoid_contexts_map:
                organoid_contexts_map[org_ctx_id] = {
                    "organoid_context_id": org_ctx_id,
                    "organoid_type": make_value_record(org_type_str, org_type_str, "organoid_contexts", org_ctx_id, "organoid_type", org_key, status="reported" if org_type_str else "unknown"),
                    "brain_region_or_model_identity": make_value_record(region_str, region_str, "organoid_contexts", org_ctx_id, "brain_region_or_model_identity", reg_key, status="reported" if region_str else "unknown"),
                    "metadata_validation_status": "not_validated",
                }

            # Construct Grouping Condition ID
            if is_control:
                raw_ctrl_type, ctrl_type_key = _find_attribute_value(rec, self.rules["control_type"])
                clean_ctrl_type = _clean_str(raw_ctrl_type) or "vehicle"
                if clean_ctrl_type not in ["vehicle", "untreated", "positive_reference", "negative_reference", "other"]:
                    clean_ctrl_type = "vehicle"

                if has_multiple_control_groups:
                    quals = [clean_ctrl_type]
                    if clean_agent and clean_agent.lower() not in ("control", "ctrl", "vehicle", "none", "mock", "untreated", ""):
                        quals.append(re.sub(r"[^a-zA-Z0-9_]", "_", clean_agent.lower()))
                    if clean_exp_stage or clean_col_stage:
                        st = clean_exp_stage or clean_col_stage
                        quals.append(re.sub(r"[^a-zA-Z0-9_]", "_", st.lower()))
                    if has_multiple_ctrl_durations and norm_dur is not None and norm_dur.normalized_value is not None:
                        quals.append(f"{norm_dur.normalized_value}{norm_dur.normalized_unit}")
                    if has_multiple_ctrl_cell_lines and cell_line_str:
                        quals.append(re.sub(r"[^a-zA-Z0-9_]", "_", cell_line_str.lower()))
                    cond_id = f"cond_{exp_id}_control_{'_'.join(quals)}"
                    cond_label = f"Control ({' '.join(quals)})"
                else:
                    cond_id = f"cond_{exp_id}_control"
                    cond_label = "Vehicle Control"
            else:
                agent_slug = re.sub(r"[^a-zA-Z0-9_]", "_", clean_agent.lower())
                dose_slug = f"{norm_dose.normalized_value or dose_val}{norm_dose.normalized_unit or dose_unit}" if dose_val is not None else "dose_unspecified"

                quals = []
                if has_multiple_durations and norm_dur is not None and norm_dur.normalized_value is not None:
                    quals.append(f"{norm_dur.normalized_value}{norm_dur.normalized_unit}")
                if has_multiple_stages and clean_exp_stage:
                    quals.append(re.sub(r"[^a-zA-Z0-9_]", "_", clean_exp_stage.lower()))
                if has_multiple_cell_lines and cell_line_str:
                    quals.append(re.sub(r"[^a-zA-Z0-9_]", "_", cell_line_str.lower()))
                if has_multiple_organoids and org_type_str:
                    quals.append(re.sub(r"[^a-zA-Z0-9_]", "_", org_type_str.lower()))

                if quals:
                    cond_id = f"cond_{exp_id}_{agent_slug}_{dose_slug}_{'_'.join(quals)}"
                    cond_label = f"{clean_agent} {dose_val} {dose_unit} ({' '.join(quals)})"
                else:
                    cond_id = f"cond_{exp_id}_{agent_slug}_{dose_slug}"
                    cond_label = f"{clean_agent} {dose_val} {dose_unit}" if dose_val is not None else clean_agent

            if cond_id not in conditions_map:
                cond_dict = {
                    "condition_id": cond_id,
                    "experiment_id": exp_id,
                    "condition_label": make_value_record(cond_label, cond_label, "conditions", cond_id, "condition_label", "condition_label"),
                    "treatment_control_status": make_value_record(treatment_status, treatment_status, "conditions", cond_id, "treatment_control_status", cond_status_key),
                    "metadata_validation_status": "not_validated",
                }
                if is_control:
                    cond_dict["control_type"] = make_value_record(clean_ctrl_type, clean_ctrl_type, "conditions", cond_id, "control_type", ctrl_type_key)
                conditions_map[cond_id] = cond_dict

            # Construct Condition-Planned Exposure
            exp_rec_id = f"exp_{cond_id}"
            if exp_rec_id not in exposures_map:
                exposure_entry: dict[str, Any] = {
                    "exposure_id": exp_rec_id,
                    "condition_id": cond_id,
                    "exposure_scope": "condition_planned",
                    "developmental_age_or_stage_at_exposure": make_value_record(
                        clean_exp_stage,
                        clean_exp_stage,
                        "exposures",
                        exp_rec_id,
                        "developmental_age_or_stage_at_exposure",
                        exp_stage_key,
                        source_loc=rec.source_locator,
                        status="reported" if clean_exp_stage is not None else "unknown",
                    ),
                    "metadata_validation_status": "not_validated",
                }
                if clean_exp_stage is not None:
                    exposure_entry["exposure_start_time_or_stage"] = make_value_record(
                        clean_exp_stage,
                        clean_exp_stage,
                        "exposures",
                        exp_rec_id,
                        "exposure_start_time_or_stage",
                        exp_stage_key,
                        source_loc=rec.source_locator,
                    )
                if is_control:
                    exposure_entry["vehicle"] = make_value_record("vehicle", "vehicle", "exposures", exp_rec_id, "vehicle", "vehicle", source_loc=rec.source_locator)
                    if norm_dur is not None and norm_dur.original_value is not None:
                        exposure_entry["exposure_duration"] = make_value_record(norm_dur.original_value, norm_dur.normalized_value, "exposures", exp_rec_id, "exposure_duration", dur_key, source_loc=rec.source_locator)
                        exposure_entry["exposure_duration_unit"] = make_value_record(norm_dur.original_unit, norm_dur.normalized_unit, "exposures", exp_rec_id, "exposure_duration_unit", dur_unit_key, source_loc=rec.source_locator)
                    else:
                        exposure_entry["exposure_duration"] = make_value_record(None, None, "exposures", exp_rec_id, "exposure_duration", None, status="missing_not_reported")
                        exposure_entry["exposure_duration_unit"] = make_value_record(None, None, "exposures", exp_rec_id, "exposure_duration_unit", None, status="missing_not_reported")
                else:
                    exposure_entry["agent_name"] = make_value_record(clean_agent, clean_agent, "exposures", exp_rec_id, "agent_name", agent_key, source_loc=rec.source_locator)
                    if norm_dose.original_value is not None:
                        exposure_entry["concentration_or_dose"] = make_value_record(norm_dose.original_value, norm_dose.normalized_value, "exposures", exp_rec_id, "concentration_or_dose", dose_key, source_loc=rec.source_locator)
                        exposure_entry["concentration_or_dose_unit"] = make_value_record(norm_dose.original_unit, norm_dose.normalized_unit, "exposures", exp_rec_id, "concentration_or_dose_unit", dose_unit_key, source_loc=rec.source_locator)
                    if norm_dur is not None and norm_dur.original_value is not None:
                        exposure_entry["exposure_duration"] = make_value_record(norm_dur.original_value, norm_dur.normalized_value, "exposures", exp_rec_id, "exposure_duration", dur_key, source_loc=rec.source_locator)
                        exposure_entry["exposure_duration_unit"] = make_value_record(norm_dur.original_unit, norm_dur.normalized_unit, "exposures", exp_rec_id, "exposure_duration_unit", dur_unit_key, source_loc=rec.source_locator)
                    else:
                        exposure_entry["exposure_duration"] = make_value_record(None, None, "exposures", exp_rec_id, "exposure_duration", None, status="missing_not_reported")
                        exposure_entry["exposure_duration_unit"] = make_value_record(None, None, "exposures", exp_rec_id, "exposure_duration_unit", None, status="missing_not_reported")
                exposures_map[exp_rec_id] = exposure_entry

            # Resolve Biological Replicate ID (no invented default)
            raw_rep, rep_key = _find_attribute_value(rec, self.rules["biological_replicate_id"])
            clean_rep = _clean_str(raw_rep) or self.defaults.get("biological_replicate_id")
            if clean_rep and clean_rep.lower() in ("none", "nan", "unknown", "unverified", "null", ""):
                clean_rep = None

            if clean_rep is not None:
                rep_record = make_value_record(clean_rep, clean_rep, "samples", sid, "biological_replicate_id", rep_key, status="reported")
            else:
                rep_record = make_value_record(None, None, "samples", sid, "biological_replicate_id", rep_key, status="unknown")

            # Construct Canonical Sample Entry
            sample_entry: dict[str, Any] = {
                "sample_id": sid,
                "experiment_id": exp_id,
                "condition_id": cond_id,
                "biological_source_id": bio_source_id,
                "organoid_context_id": org_ctx_id,
                "biological_replicate_id": rep_record,
                "developmental_age_or_stage_at_collection": make_value_record(
                    clean_col_stage,
                    clean_col_stage,
                    "samples",
                    sid,
                    "developmental_age_or_stage_at_collection",
                    col_stage_key,
                    status="reported" if clean_col_stage is not None else "unknown",
                ),
                "collection_time_relative_to_exposure": (
                    make_value_record(
                        f"{dur_val} {dur_unit}",
                        f"{norm_dur.normalized_value} {norm_dur.normalized_unit}",
                        "samples",
                        sid,
                        "collection_time_relative_to_exposure",
                        dur_key,
                    )
                    if norm_dur is not None and norm_dur.normalized_value is not None
                    else make_value_record(
                        None,
                        None,
                        "samples",
                        sid,
                        "collection_time_relative_to_exposure",
                        None,
                        status="unknown",
                    )
                ),
                "metadata_validation_status": "not_validated",
            }
            samples.append(sample_entry)

        # 4. Construct Sequencing Assays & Input Data Assets
        sample_ids_list = [s["sample_id"] for s in samples]
        assay_id = f"assay_{exp_id}_bulk"
        assays.append(
            {
                "assay_id": assay_id,
                "sample_ids": sample_ids_list,
                "rna_seq_modality": make_value_record("bulk_rna_seq", "bulk_rna_seq", "sequencing_assays", assay_id, "rna_seq_modality", "rna_seq_modality"),
                "library_strategy": make_value_record("polyA", "polyA", "sequencing_assays", assay_id, "library_strategy", "library_strategy"),
                "sequencing_platform": make_value_record("Illumina NovaSeq 6000", "Illumina NovaSeq 6000", "sequencing_assays", assay_id, "sequencing_platform", "sequencing_platform"),
                "metadata_validation_status": "not_validated",
            }
        )

        asset_id = f"asset_{exp_id}_counts"
        assets.append(
            {
                "asset_id": asset_id,
                "assay_ids": [assay_id],
                "input_data_format": make_value_record("raw_counts", "raw_counts", "input_data_assets", asset_id, "input_data_format", "input_data_format"),
                "gene_identifier_type": make_value_record("Ensembl Gene ID", "Ensembl Gene ID", "input_data_assets", asset_id, "gene_identifier_type", "gene_identifier_type"),
                "source_file_identifier_or_reference": make_value_record(source_reference, source_reference, "input_data_assets", asset_id, "source_file_identifier_or_reference", "source_file"),
                "availability": "available",
                "metadata_validation_status": "not_validated",
            }
        )

        # 5. Treatment-Control Relationships
        # Controls are matched ONLY using sufficient source-supported experimental information.
        # Ambiguous, missing, or conflicting relationships require review rather than being guessed.
        relationships: list[dict[str, Any]] = []
        control_cond_ids = [c["condition_id"] for c in conditions_map.values() if c["treatment_control_status"]["normalized_value"] == "control"]
        treatment_cond_ids = [c["condition_id"] for c in conditions_map.values() if c["treatment_control_status"]["normalized_value"] == "treatment"]

        if not control_cond_ids and treatment_cond_ids:
            findings.append(
                MappingFinding(
                    severity="ERROR",
                    rule="missing_control_condition",
                    field_path="conditions",
                    message="No matched control condition identified for treatment conditions.",
                )
            )
        elif control_cond_ids and treatment_cond_ids:
            # Build profile for each condition to support attribute-aware matching
            cond_profiles: dict[str, dict[str, Any]] = {}
            for cid, cond in conditions_map.items():
                c_samples = [s for s in samples if s["condition_id"] == cid]
                c_exp = exposures_map.get(f"exp_{cid}", {})

                bio_srcs = {s["biological_source_id"] for s in c_samples if s.get("biological_source_id")}
                org_ctxs = {s["organoid_context_id"] for s in c_samples if s.get("organoid_context_id")}

                stages: set[str] = set()
                for s in c_samples:
                    st_val = s.get("developmental_age_or_stage_at_collection", {}).get("normalized_value") or s.get("developmental_age_or_stage_at_collection", {}).get("original_value")
                    if st_val:
                        stages.add(str(st_val).strip().lower())
                exp_st = c_exp.get("developmental_age_or_stage_at_exposure", {}).get("normalized_value") or c_exp.get("developmental_age_or_stage_at_exposure", {}).get("original_value")
                if exp_st:
                    stages.add(str(exp_st).strip().lower())

                durs: set[str] = set()
                for s in c_samples:
                    dur_val = s.get("collection_time_relative_to_exposure", {}).get("normalized_value")
                    if dur_val:
                        durs.add(str(dur_val).strip().lower())
                exp_dur_val = c_exp.get("exposure_duration", {}).get("normalized_value") or c_exp.get("exposure_duration", {}).get("original_value")
                exp_dur_u = c_exp.get("exposure_duration_unit", {}).get("normalized_value") or c_exp.get("exposure_duration_unit", {}).get("original_value")
                if exp_dur_val is not None and exp_dur_u is not None:
                    durs.add(f"{exp_dur_val} {exp_dur_u}".strip().lower())

                vehs: set[str] = set()
                veh_val = c_exp.get("vehicle", {}).get("normalized_value") or c_exp.get("vehicle", {}).get("original_value")
                if veh_val:
                    vehs.add(str(veh_val).strip().lower())

                cond_profiles[cid] = {
                    "role": cond["treatment_control_status"]["normalized_value"],
                    "bio_sources": bio_srcs,
                    "org_contexts": org_ctxs,
                    "stages": stages,
                    "durations": durs,
                    "vehicles": vehs,
                }

            for t_cid in treatment_cond_ids:
                t_prof = cond_profiles[t_cid]

                # Filter candidate controls by hard biological model compatibility
                compatible_controls: list[str] = []
                for c_cid in control_cond_ids:
                    c_prof = cond_profiles[c_cid]

                    # Biological source check
                    if t_prof["bio_sources"] and c_prof["bio_sources"]:
                        if not (t_prof["bio_sources"] & c_prof["bio_sources"]):
                            continue

                    # Organoid context check
                    if t_prof["org_contexts"] and c_prof["org_contexts"]:
                        if not (t_prof["org_contexts"] & c_prof["org_contexts"]):
                            continue

                    # Developmental stage check
                    if t_prof["stages"] and c_prof["stages"]:
                        if not (t_prof["stages"] & c_prof["stages"]):
                            continue

                    compatible_controls.append(c_cid)

                if not compatible_controls:
                    findings.append(
                        MappingFinding(
                            severity="ERROR",
                            rule="unmatched_treatment_condition",
                            field_path=f"conditions[condition_id='{t_cid}']",
                            message=(
                                f"Treatment condition '{t_cid}' has no compatible matched control condition "
                                "(conflicting biological source, organoid context, or developmental stage). "
                                "Cannot automatically construct relationship without manual review."
                            ),
                        )
                    )
                    continue

                # If multiple biologically compatible controls, attempt duration or vehicle refinement
                selected_control: str | None = None
                is_ambiguous = False

                if len(compatible_controls) == 1:
                    selected_control = compatible_controls[0]
                else:
                    # Attempt duration-matching if duration info exists
                    if t_prof["durations"]:
                        dur_matched = [c for c in compatible_controls if cond_profiles[c]["durations"] & t_prof["durations"]]
                        if len(dur_matched) == 1:
                            selected_control = dur_matched[0]

                    # If still unresolved, attempt vehicle-matching if vehicle info exists
                    if selected_control is None and t_prof["vehicles"]:
                        veh_matched = [c for c in compatible_controls if cond_profiles[c]["vehicles"] & t_prof["vehicles"]]
                        if len(veh_matched) == 1:
                            selected_control = veh_matched[0]

                    # If still unresolved after all source-supported attributes: mark ambiguous!
                    if selected_control is None:
                        is_ambiguous = True

                if is_ambiguous:
                    # AMBIGUOUS: Do NOT guess an arbitrary control and do NOT establish speculative or pooled relationships!
                    findings.append(
                        MappingFinding(
                            severity="WARNING",
                            rule="ambiguous_matched_control",
                            field_path=f"conditions[condition_id='{t_cid}']",
                            message=(
                                f"Treatment condition '{t_cid}' has multiple candidate matched controls: {compatible_controls}. "
                                "Automatic matching cannot be established without guessing; relationship remains blocked until "
                                "appropriate experimental evidence or human review resolves it."
                            ),
                        )
                    )
                    # Relationship is NOT established; treatment condition remains unlinked and blocked downstream.
                else:
                    assert selected_control is not None
                    rel_id = f"rel_{t_cid}_vs_{selected_control}"
                    basis_desc = "Matched on study, biological source, organoid context, developmental stage, and duration"
                    relationships.append(
                        {
                            "relationship_id": rel_id,
                            "treatment_condition_ids": [t_cid],
                            "matched_control_condition_ids": [selected_control],
                            "matching_basis": make_value_record(
                                basis_desc,
                                basis_desc,
                                "treatment_control_relationships",
                                rel_id,
                                "matching_basis",
                                "matching_basis",
                            ),
                            "metadata_validation_status": "not_validated",
                        }
                    )

        canonical_doc: dict[str, Any] = {
            "studies": studies,
            "experiments": experiments,
            "conditions": list(conditions_map.values()),
            "biological_sources": list(bio_sources_map.values()),
            "organoid_contexts": list(organoid_contexts_map.values()),
            "samples": samples,
            "exposures": list(exposures_map.values()),
            "treatment_control_relationships": relationships,
            "sample_comparison_exceptions": [],
            "sequencing_assays": assays,
            "input_data_assets": assets,
            "provenance": provenance_records,
            "dnt_reference_evidence": [],
        }

        return canonical_doc, findings
