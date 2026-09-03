"""Schema-driven validation plus accepted cross-record integrity rules."""

from collections.abc import Iterable
from pathlib import Path
from typing import Any

from .result import Finding, Severity, ValidationResult, ValidationStatus
from .schema import DEFAULT_SCHEMA_PATH, identifier_field, load_schema


_SEVERITY_ORDER = {Severity.ERROR: 0, Severity.WARNING: 1}
_MISSING_STATUSES = {"missing_not_reported", "unknown"}


class MetadataValidator:
    """Validate canonical metadata without mutating or completing it."""

    def __init__(self, schema_path: str | Path = DEFAULT_SCHEMA_PATH) -> None:
        self.schema = load_schema(schema_path)
        self.entities: dict[str, dict[str, Any]] = self.schema["entities"]
        self.vocabularies: dict[str, list[str]] = self.schema["controlled_vocabularies"]
        self.conventions: dict[str, dict[str, Any]] = self.schema["record_conventions"]
        self._findings: list[Finding] = []
        self._records: dict[str, list[dict[str, Any]]] = {}
        self._indexes: dict[str, dict[str, dict[str, Any]]] = {}

    def validate(self, metadata: Any) -> ValidationResult:
        """Return deterministic findings for a canonical metadata document."""
        self._findings = []
        self._records = {}
        self._indexes = {}
        if not isinstance(metadata, dict):
            self._add(Severity.ERROR, "document", None, "$", "document_type", "Metadata document must be a mapping.")
            return self._result()

        self._validate_entity_collections(metadata)
        self._build_indexes()
        self._validate_conditional_fields()
        self._validate_foreign_keys()
        self._validate_condition_sample_integrity()
        self._validate_exposures()
        self._validate_comparisons()
        self._validate_provenance_records()
        return self._result()

    def _validate_entity_collections(self, metadata: dict[str, Any]) -> None:
        for entity in sorted(set(metadata) - set(self.entities)):
            self._add(Severity.ERROR, entity, None, entity, "unknown_entity", f"Entity collection '{entity}' is not defined by the contract.")
        for entity, definition in self.entities.items():
            value = metadata.get(entity)
            if value is None:
                self._add(Severity.ERROR, entity, None, entity, "entity_presence", f"Required entity collection '{entity}' is missing.")
                self._records[entity] = []
                continue
            if not isinstance(value, list):
                self._add(Severity.ERROR, entity, None, entity, "entity_structure", f"Entity collection '{entity}' must be an array.")
                self._records[entity] = []
                continue
            records = [record for record in value if isinstance(record, dict)]
            self._records[entity] = records
            for position, record in enumerate(value):
                if not isinstance(record, dict):
                    self._add(Severity.ERROR, entity, None, f"{entity}[{position}]", "record_structure", "Entity record must be an object.")
                    continue
                self._validate_record(entity, definition, record, position)

    def _validate_record(self, entity: str, definition: dict[str, Any], record: dict[str, Any], position: int) -> None:
        id_field = identifier_field(entity, definition)
        identifier = record.get(id_field) if id_field else None
        base = f"{entity}[{position}]"
        unknown = sorted(set(record) - set(definition["fields"]))
        for field in unknown:
            self._add(Severity.ERROR, entity, self._id(identifier), f"{base}.{field}", "unknown_field", f"Field '{field}' is not defined by the contract.")
        for field in definition.get("item_required", []):
            if field not in record:
                self._add(Severity.ERROR, entity, self._id(identifier), f"{base}.{field}", "required_field", f"Required field '{field}' is missing.")
        for field, value in record.items():
            spec = definition["fields"].get(field)
            if spec:
                self._validate_field(entity, self._id(identifier), f"{base}.{field}", value, spec)

    def _validate_field(self, entity: str, identifier: str | None, path: str, value: Any, spec: dict[str, Any]) -> None:
        declared = spec.get("type")
        if declared == "value_record":
            self._validate_value_record(entity, identifier, path, value, spec)
            return
        elif declared == "array[accession]":
            self._validate_object_array(entity, identifier, path, value, "accession")
        elif declared == "array[evidence_scope_assertion]":
            self._validate_object_array(entity, identifier, path, value, "evidence_scope_assertion")
        elif declared == "array[string]":
            if not isinstance(value, list) or any(not isinstance(item, str) or not item for item in value):
                self._add(Severity.ERROR, entity, identifier, path, "field_type", "Value must be an array of non-empty strings.")
        elif not self._matches_type(value, declared):
            self._add(Severity.ERROR, entity, identifier, path, "field_type", f"Value does not match declared type {declared!r}.")
        self._validate_controlled_value(entity, identifier, path, value, spec)

    def _validate_value_record(self, entity: str, identifier: str | None, path: str, value: Any, field_spec: dict[str, Any]) -> None:
        definition = self.conventions["value_record"]
        if not isinstance(value, dict):
            self._add(Severity.ERROR, entity, identifier, path, "value_record_type", "Value must be a value_record object.")
            return
        for field in sorted(set(value) - set(definition["fields"])):
            self._add(Severity.ERROR, entity, identifier, f"{path}.{field}", "unknown_field", f"Field '{field}' is not defined for value_record.")
        for field in definition["required"]:
            if field not in value:
                self._add(Severity.ERROR, entity, identifier, f"{path}.{field}", "value_record_required", f"value_record requires '{field}'.")
        for field in ("original_value", "normalized_value"):
            if field in value and not self._matches_type(value[field], definition["fields"][field]["type"]):
                self._add(Severity.ERROR, entity, identifier, f"{path}.{field}", "field_type", f"Value does not match declared type {definition['fields'][field]['type']!r}.")
        status = value.get("value_status")
        allowed_statuses = definition["fields"]["value_status"]["allowed"]
        if status not in allowed_statuses:
            self._add(Severity.ERROR, entity, identifier, f"{path}.value_status", "missingness_vocabulary", f"value_status must be one of {allowed_statuses}.")
        elif status == "reported" and ("original_value" not in value or value.get("original_value") is None):
            self._add(Severity.ERROR, entity, identifier, f"{path}.original_value", "reported_value_required", "A reported value requires a non-null original_value.")
        elif status in _MISSING_STATUSES:
            self._add(Severity.WARNING, entity, identifier, path, "explicit_missingness", f"Biological metadata is explicitly marked '{status}'; no value was inferred.")
        elif status == "not_applicable" and self._has_non_null_value(value):
            self._add(Severity.ERROR, entity, identifier, path, "not_applicable_value", "A not_applicable value_record cannot contain a non-null original or normalized value.")
        if status != "reported" and status != "not_applicable" and self._has_non_null_value(value):
            self._add(Severity.ERROR, entity, identifier, path, "missing_value_conflict", f"A '{status}' value_record cannot contain a non-null original or normalized value.")
        provenance_ids = value.get("assertion_provenance_ids")
        if provenance_ids is not None and (not isinstance(provenance_ids, list) or any(not isinstance(item, str) or not item for item in provenance_ids)):
            self._add(Severity.ERROR, entity, identifier, f"{path}.assertion_provenance_ids", "provenance_list_type", "assertion_provenance_ids must be an array of non-empty strings.")
        self._validate_controlled_value(entity, identifier, path, value, field_spec)

    def _validate_object_array(self, entity: str, identifier: str | None, path: str, value: Any, structure_name: str) -> None:
        if not isinstance(value, list):
            self._add(Severity.ERROR, entity, identifier, path, "field_type", f"Value must be an array of {structure_name} objects.")
            return
        definition = self.conventions[structure_name]
        for position, item in enumerate(value):
            item_path = f"{path}[{position}]"
            if not isinstance(item, dict):
                self._add(Severity.ERROR, entity, identifier, item_path, "field_type", f"Value must be a {structure_name} object.")
                continue
            for field in definition["required"]:
                if field not in item:
                    self._add(Severity.ERROR, entity, identifier, f"{item_path}.{field}", "required_field", f"{structure_name} requires '{field}'.")
            for field, field_value in item.items():
                spec = definition["fields"].get(field)
                if spec is None:
                    self._add(Severity.ERROR, entity, identifier, f"{item_path}.{field}", "unknown_field", f"Field '{field}' is not defined for {structure_name}.")
                elif spec.get("type") == "value_record":
                    self._validate_value_record(entity, identifier, f"{item_path}.{field}", field_value, spec)
                elif not self._matches_type(field_value, spec.get("type")):
                    self._add(Severity.ERROR, entity, identifier, f"{item_path}.{field}", "field_type", f"Value does not match declared type {spec.get('type')!r}.")
                else:
                    self._validate_controlled_value(entity, identifier, f"{item_path}.{field}", field_value, spec)

    def _validate_controlled_value(self, entity: str, identifier: str | None, path: str, value: Any, spec: dict[str, Any]) -> None:
        reference = spec.get("allowed_ref")
        if not reference:
            return
        allowed = self.vocabularies.get(reference) or self.conventions.get(reference, {}).get("allowed")
        candidate = value
        if spec.get("type") == "value_record" and isinstance(value, dict) and value.get("value_status") == "reported":
            candidate = value.get("normalized_value") if value.get("normalized_value") is not None else value.get("original_value")
        if candidate is not None and candidate not in allowed:
            self._add(Severity.ERROR, entity, identifier, path, "controlled_vocabulary", f"Value {candidate!r} is not in controlled vocabulary '{reference}'.")

    def _build_indexes(self) -> None:
        for entity, definition in self.entities.items():
            id_field = identifier_field(entity, definition)
            index: dict[str, dict[str, Any]] = {}
            if not id_field:
                self._indexes[entity] = index
                continue
            for position, record in enumerate(self._records.get(entity, [])):
                identifier = record.get(id_field)
                if not isinstance(identifier, str) or not identifier:
                    self._add(Severity.ERROR, entity, self._id(identifier), f"{entity}[{position}].{id_field}", "identifier", "Identifier must be a non-empty string.")
                elif identifier in index:
                    self._add(Severity.ERROR, entity, identifier, f"{entity}[{position}].{id_field}", "identifier_unique", f"Duplicate identifier '{identifier}'.")
                else:
                    index[identifier] = record
            self._indexes[entity] = index

    def _validate_foreign_keys(self) -> None:
        for entity, definition in self.entities.items():
            id_field = identifier_field(entity, definition)
            for position, record in enumerate(self._records.get(entity, [])):
                identifier = self._id(record.get(id_field)) if id_field else None
                for field, spec in definition["fields"].items():
                    target = spec.get("foreign_key")
                    if not target or field not in record or record[field] is None:
                        continue
                    target_entity = target.split(".", 1)[0]
                    values: Iterable[Any] = record[field] if spec.get("type") == "array[string]" and isinstance(record[field], list) else [record[field]]
                    for value in values:
                        if isinstance(value, str) and value not in self._indexes.get(target_entity, {}):
                            self._add(Severity.ERROR, entity, identifier, f"{entity}[{position}].{field}", "foreign_key", f"Reference '{value}' does not exist in '{target_entity}'.")

    def _validate_conditional_fields(self) -> None:
        for position, study in enumerate(self._records.get("studies", [])):
            if study.get("source_type") == "public" and not study.get("public_accessions"):
                self._add(Severity.ERROR, "studies", self._id(study.get("study_id")), f"studies[{position}].public_accessions", "conditional_required", "A public study requires at least one public accession.")
        for position, condition in enumerate(self._records.get("conditions", [])):
            if self._reported_value(condition.get("treatment_control_status")) == "control" and "control_type" not in condition:
                self._add(Severity.ERROR, "conditions", self._id(condition.get("condition_id")), f"conditions[{position}].control_type", "conditional_required", "A control condition requires control_type.")
        matrix_formats = {"raw_counts", "UMI_counts", "H5AD", "RDS", "processed_expression_matrix"}
        for position, asset in enumerate(self._records.get("input_data_assets", [])):
            identifier = self._id(asset.get("asset_id"))
            base = f"input_data_assets[{position}]"
            if self._reported_value(asset.get("input_data_format")) in matrix_formats and "gene_identifier_type" not in asset:
                self._add(Severity.ERROR, "input_data_assets", identifier, f"{base}.gene_identifier_type", "conditional_required", "An expression matrix requires gene_identifier_type.")
            if ("checksum_algorithm" in asset) != ("checksum_value" in asset):
                missing = "checksum_value" if "checksum_algorithm" in asset else "checksum_algorithm"
                self._add(Severity.ERROR, "input_data_assets", identifier, f"{base}.{missing}", "conditional_required", "Checksum algorithm and value must be supplied together.")
        for position, evidence in enumerate(self._records.get("dnt_reference_evidence", [])):
            if evidence.get("scope_type") == "context_limited" and not evidence.get("scope_assertions"):
                self._add(Severity.ERROR, "dnt_reference_evidence", self._id(evidence.get("evidence_id")), f"dnt_reference_evidence[{position}].scope_assertions", "conditional_required", "Context-limited evidence requires at least one source-supported scope assertion.")

    def _validate_condition_sample_integrity(self) -> None:
        for position, sample in enumerate(self._records.get("samples", [])):
            condition = self._indexes.get("conditions", {}).get(sample.get("condition_id"))
            if condition and sample.get("experiment_id") != condition.get("experiment_id"):
                self._add(Severity.ERROR, "samples", self._id(sample.get("sample_id")), f"samples[{position}].condition_id", "condition_experiment_consistency", "Sample and condition must belong to the same experiment.")
            if condition and self._indexes.get("experiments", {}).get(sample.get("experiment_id")) is None:
                continue
            if condition and any(item.get("condition_id") == sample.get("condition_id") for item in self._records.get("exposures", [])) and "collection_time_relative_to_exposure" not in sample:
                self._add(Severity.ERROR, "samples", self._id(sample.get("sample_id")), f"samples[{position}].collection_time_relative_to_exposure", "conditional_required", "An exposed sample requires collection_time_relative_to_exposure.")

    def _validate_exposures(self) -> None:
        planned_by_condition = {
            exposure.get("condition_id")
            for exposure in self._records.get("exposures", [])
            if exposure.get("exposure_scope") == "condition_planned"
        }
        for position, exposure in enumerate(self._records.get("exposures", [])):
            identifier = self._id(exposure.get("exposure_id"))
            base = f"exposures[{position}]"
            scope = exposure.get("exposure_scope")
            condition = self._indexes.get("conditions", {}).get(exposure.get("condition_id"))
            if scope == "condition_planned":
                for field in ("sample_id", "deviation_from_exposure_id"):
                    if exposure.get(field) is not None:
                        self._add(Severity.ERROR, "exposures", identifier, f"{base}.{field}", "planned_exposure_target", f"condition_planned exposure requires null or absent {field}.")
            elif scope == "sample_deviation":
                self._require_fields("exposures", identifier, base, exposure, ("sample_id", "deviation_reason"))
                sample = self._indexes.get("samples", {}).get(exposure.get("sample_id"))
                if sample and sample.get("condition_id") != exposure.get("condition_id"):
                    self._add(Severity.ERROR, "exposures", identifier, f"{base}.sample_id", "deviation_condition", "Deviating sample must be assigned to the exposure condition.")
                if exposure.get("condition_id") not in planned_by_condition:
                    self._add(Severity.ERROR, "exposures", identifier, f"{base}.condition_id", "deviation_supplements_plan", "A sample deviation requires at least one planned exposure for its condition.")
                planned_id = exposure.get("deviation_from_exposure_id")
                if planned_id:
                    planned = self._indexes.get("exposures", {}).get(planned_id)
                    if planned and (planned.get("exposure_scope") != "condition_planned" or planned.get("condition_id") != exposure.get("condition_id")):
                        self._add(Severity.ERROR, "exposures", identifier, f"{base}.deviation_from_exposure_id", "deviation_planned_reference", "Referenced exposure must be condition_planned for the same condition.")
            if condition:
                condition_status = self._reported_value(condition.get("treatment_control_status"))
                if condition_status == "treatment":
                    self._require_fields("exposures", identifier, base, exposure, ("agent_name", "concentration_or_dose", "developmental_age_or_stage_at_exposure", "exposure_duration"))
            self._require_unit_for_reported(exposure, "concentration_or_dose", "concentration_or_dose_unit", identifier, base)
            self._require_unit_for_reported(exposure, "exposure_duration", "exposure_duration_unit", identifier, base)
            self._require_unit_for_reported(exposure, "washout_or_recovery_duration", "washout_or_recovery_duration_unit", identifier, base)

    def _validate_comparisons(self) -> None:
        conditions = self._indexes.get("conditions", {})
        for position, relationship in enumerate(self._records.get("treatment_control_relationships", [])):
            identifier = self._id(relationship.get("relationship_id"))
            experiment_ids: set[str] = set()
            for field, expected in (("treatment_condition_ids", "treatment"), ("matched_control_condition_ids", "control")):
                values = relationship.get(field)
                if not isinstance(values, list) or not values:
                    self._add(Severity.ERROR, "treatment_control_relationships", identifier, f"treatment_control_relationships[{position}].{field}", "comparison_non_empty", f"{field} must contain at least one condition.")
                    continue
                for condition_id in values:
                    condition = conditions.get(condition_id)
                    if not condition:
                        continue
                    experiment_ids.add(condition.get("experiment_id"))
                    actual = self._reported_value(condition.get("treatment_control_status"))
                    if actual != expected:
                        self._add(Severity.ERROR, "treatment_control_relationships", identifier, f"treatment_control_relationships[{position}].{field}", "comparison_condition_role", f"Condition '{condition_id}' must have status '{expected}'.")
            if len(experiment_ids) > 1:
                self._add(Severity.ERROR, "treatment_control_relationships", identifier, f"treatment_control_relationships[{position}]", "comparison_experiment", "Compared conditions must belong to the same experiment.")

        relationships = self._indexes.get("treatment_control_relationships", {})
        samples = self._indexes.get("samples", {})
        for position, exception in enumerate(self._records.get("sample_comparison_exceptions", [])):
            relationship = relationships.get(exception.get("relationship_id"))
            if not relationship:
                continue
            participating = set(relationship.get("treatment_condition_ids", [])) | set(relationship.get("matched_control_condition_ids", []))
            for field in ("affected_sample_ids", "replacement_control_sample_ids"):
                for sample_id in exception.get(field, []) if isinstance(exception.get(field, []), list) else []:
                    sample = samples.get(sample_id)
                    if sample and sample.get("condition_id") not in participating:
                        self._add(Severity.ERROR, "sample_comparison_exceptions", self._id(exception.get("exception_id")), f"sample_comparison_exceptions[{position}].{field}", "comparison_exception_scope", f"Sample '{sample_id}' is outside the referenced condition comparison.")
            for sample_id in exception.get("replacement_control_sample_ids", []) if isinstance(exception.get("replacement_control_sample_ids", []), list) else []:
                sample = samples.get(sample_id)
                if sample and sample.get("condition_id") not in set(relationship.get("matched_control_condition_ids", [])):
                    self._add(Severity.ERROR, "sample_comparison_exceptions", self._id(exception.get("exception_id")), f"sample_comparison_exceptions[{position}].replacement_control_sample_ids", "replacement_control_role", f"Replacement sample '{sample_id}' must belong to a matched control condition.")

    def _validate_provenance_records(self) -> None:
        provenance_ids = self._indexes.get("provenance", {})
        for entity, records in self._records.items():
            for position, record in enumerate(records):
                identifier = self._id(record.get(identifier_field(entity, self.entities[entity])))
                for field, value in record.items():
                    self._validate_provenance_links(entity, identifier, f"{entity}[{position}].{field}", value, provenance_ids)
        for position, record in enumerate(self._records.get("provenance", [])):
            entity = record.get("entity_type")
            identifier = self._id(record.get("provenance_id"))
            if entity not in self.entities:
                self._add(Severity.ERROR, "provenance", identifier, f"provenance[{position}].entity_type", "provenance_entity", f"Unknown provenance entity_type '{entity}'.")
                continue
            target = self._indexes.get(entity, {}).get(record.get("entity_id"))
            if target is None:
                self._add(Severity.ERROR, "provenance", identifier, f"provenance[{position}].entity_id", "provenance_entity_reference", "Provenance entity_id does not resolve.")
            elif record.get("field_path") not in self.entities[entity]["fields"]:
                self._add(Severity.ERROR, "provenance", identifier, f"provenance[{position}].field_path", "provenance_field", "Provenance field_path is not defined for the target entity.")

    def _validate_provenance_links(self, entity: str, identifier: str | None, path: str, value: Any, provenance_ids: dict[str, Any]) -> None:
        if isinstance(value, dict):
            if "assertion_provenance_ids" in value:
                links = value.get("assertion_provenance_ids")
                if isinstance(links, list):
                    if value.get("value_status") == "reported" and not links:
                        self._add(Severity.WARNING, entity, identifier, f"{path}.assertion_provenance_ids", "provenance_missing", "Reported metadata has no provenance reference.")
                    for link in links:
                        if isinstance(link, str) and link not in provenance_ids:
                            self._add(Severity.ERROR, entity, identifier, f"{path}.assertion_provenance_ids", "provenance_reference", f"Provenance reference '{link}' does not exist.")
            for key, nested in value.items():
                if key != "assertion_provenance_ids":
                    self._validate_provenance_links(entity, identifier, f"{path}.{key}", nested, provenance_ids)
        elif isinstance(value, list):
            for position, nested in enumerate(value):
                self._validate_provenance_links(entity, identifier, f"{path}[{position}]", nested, provenance_ids)

    def _require_fields(self, entity: str, identifier: str | None, base: str, record: dict[str, Any], fields: tuple[str, ...]) -> None:
        for field in fields:
            if field not in record or record[field] is None:
                self._add(Severity.ERROR, entity, identifier, f"{base}.{field}", "conditional_required", f"Conditionally required field '{field}' is missing.")

    def _require_unit_for_reported(self, record: dict[str, Any], value_field: str, unit_field: str, identifier: str | None, base: str) -> None:
        if self._is_reported(record.get(value_field)) and unit_field not in record:
            self._add(Severity.ERROR, "exposures", identifier, f"{base}.{unit_field}", "conditional_required", f"Reported '{value_field}' requires '{unit_field}'.")

    def _result(self) -> ValidationResult:
        findings = tuple(sorted(self._findings, key=lambda item: (_SEVERITY_ORDER[item.severity], item.entity_type, item.entity_identifier or "", item.path, item.rule, item.message)))
        status = ValidationStatus.INVALID if any(item.severity == Severity.ERROR for item in findings) else ValidationStatus.VALID_WITH_WARNINGS if findings else ValidationStatus.VALID
        return ValidationResult(status, findings)

    def _add(self, severity: Severity, entity: str, identifier: str | None, path: str, rule: str, message: str) -> None:
        self._findings.append(Finding(severity, entity, identifier, path, rule, message))

    @staticmethod
    def _matches_type(value: Any, declared: Any) -> bool:
        types = declared if isinstance(declared, list) else [declared]
        mapping = {"string": str, "number": (int, float), "boolean": bool, "null": type(None), "array": list, "object": dict}
        return any(name in mapping and isinstance(value, mapping[name]) and not (name == "number" and isinstance(value, bool)) for name in types)

    @staticmethod
    def _has_non_null_value(value: dict[str, Any]) -> bool:
        return value.get("original_value") is not None or value.get("normalized_value") is not None

    @staticmethod
    def _is_reported(value: Any) -> bool:
        return isinstance(value, dict) and value.get("value_status") == "reported"

    @staticmethod
    def _reported_value(value: Any) -> Any:
        if not isinstance(value, dict) or value.get("value_status") != "reported":
            return None
        return value.get("normalized_value") if value.get("normalized_value") is not None else value.get("original_value")

    @staticmethod
    def _id(value: Any) -> str | None:
        return value if isinstance(value, str) else None
