"""Focused tests for the canonical metadata validator."""

import os
import tempfile
import unittest
from copy import deepcopy

from src.metadata import MetadataValidator, ValidationStatus
from .fixtures import canonical_bulk_study, missing, reported, with_provenance


class MetadataValidatorTests(unittest.TestCase):
    def setUp(self) -> None:
        self.validator = MetadataValidator()

    def assert_valid(self, document, status=ValidationStatus.VALID):
        result = self.validator.validate(document)
        self.assertEqual(status, result.status, result.findings)
        return result

    def assert_rule(self, document, rule):
        result = self.validator.validate(document)
        self.assertEqual(ValidationStatus.INVALID, result.status)
        self.assertIn(rule, {finding.rule for finding in result.findings})

    def test_minimal_valid_public_bulk_experiment(self):
        self.assert_valid(canonical_bulk_study())

    def test_minimal_valid_in_house_experiment(self):
        self.assert_valid(canonical_bulk_study("in_house"))

    def test_treatment_condition_has_three_biological_replicates(self):
        document = canonical_bulk_study()
        samples = [item for item in document["samples"] if item["condition_id"].startswith("treatment")]
        self.assertEqual(3, len(samples))
        self.assert_valid(document)

    def test_shared_control_across_multiple_treatments(self):
        document = canonical_bulk_study()
        second = dict(document["conditions"][1], condition_id="treatment_20um_24h_d45")
        document["conditions"].append(second)
        exposure = dict(document["exposures"][1], exposure_id="treatment_exposure_2", condition_id=second["condition_id"])
        document["exposures"].append(exposure)
        document["treatment_control_relationships"].append({
            "relationship_id": "comparison_2",
            "treatment_condition_ids": [second["condition_id"]],
            "matched_control_condition_ids": ["control"],
            "matching_basis": reported("same control"),
            "metadata_validation_status": "not_validated",
        })
        self.assertNotEqual(ValidationStatus.INVALID, self.validator.validate(document).status)

    def test_multiple_concentrations(self):
        document = canonical_bulk_study()
        document["exposures"].append(dict(document["exposures"][1], exposure_id="dose_20", concentration_or_dose=reported(20)))
        self.assertNotEqual(ValidationStatus.INVALID, self.validator.validate(document).status)

    def test_multiple_durations(self):
        document = canonical_bulk_study()
        document["exposures"].append(dict(document["exposures"][1], exposure_id="duration_48", exposure_duration=reported(48)))
        self.assertNotEqual(ValidationStatus.INVALID, self.validator.validate(document).status)

    def test_multiple_developmental_stages(self):
        document = canonical_bulk_study()
        document["exposures"].append(dict(document["exposures"][1], exposure_id="stage_60", developmental_age_or_stage_at_exposure=reported("day 60")))
        self.assertNotEqual(ValidationStatus.INVALID, self.validator.validate(document).status)

    def test_multi_agent_factorial_exposure(self):
        document = canonical_bulk_study()
        document["exposures"].append(dict(document["exposures"][1], exposure_id="second_agent", agent_name=reported("Compound B"), agent_identifier=reported("example:compound-b")))
        self.assertNotEqual(ValidationStatus.INVALID, self.validator.validate(document).status)

    def test_valid_sample_specific_exposure_deviation(self):
        document = canonical_bulk_study()
        document["exposures"].append({
            "exposure_id": "sample_deviation_1",
            "condition_id": "treatment_10um_24h_d45",
            "exposure_scope": "sample_deviation",
            "sample_id": "sample_t1",
            "deviation_from_exposure_id": "treatment_exposure",
            "deviation_reason": reported("Exposure ended early"),
            "agent_name": reported("Compound A"),
            "concentration_or_dose": reported(10),
            "concentration_or_dose_unit": reported("µM"),
            "developmental_age_or_stage_at_exposure": reported("day 45"),
            "exposure_duration": reported(22),
            "exposure_duration_unit": reported("hour"),
            "metadata_validation_status": "not_validated",
        })
        self.assertNotEqual(ValidationStatus.INVALID, self.validator.validate(document).status)

    def test_invalid_sample_referencing_nonexistent_condition(self):
        document = canonical_bulk_study()
        document["samples"][0]["condition_id"] = "missing_condition"
        self.assert_rule(document, "foreign_key")

    def test_invalid_condition_experiment_mismatch(self):
        document = canonical_bulk_study()
        document["experiments"].append({"experiment_id": "experiment_2", "study_id": "study_1", "metadata_validation_status": "not_validated"})
        document["samples"][0]["experiment_id"] = "experiment_2"
        self.assert_rule(document, "condition_experiment_consistency")

    def test_invalid_controlled_vocabulary(self):
        document = canonical_bulk_study()
        document["sequencing_assays"][0]["rna_seq_modality"] = reported("spatial_rna_seq")
        self.assert_rule(document, "controlled_vocabulary")

    def test_invalid_treatment_control_foreign_key(self):
        document = canonical_bulk_study()
        document["treatment_control_relationships"][0]["matched_control_condition_ids"] = ["missing_control"]
        self.assert_rule(document, "foreign_key")

    def test_explicit_missingness_is_warning_without_inference(self):
        document = canonical_bulk_study()
        document["samples"][0]["developmental_age_or_stage_at_collection"] = missing()
        result = self.assert_valid(document, ValidationStatus.VALID_WITH_WARNINGS)
        self.assertIn("explicit_missingness", {finding.rule for finding in result.findings})

    def test_public_study_requires_accession(self):
        document = canonical_bulk_study()
        del document["studies"][0]["public_accessions"]
        self.assert_rule(document, "conditional_required")

    def test_in_house_study_does_not_require_public_accession(self):
        self.assert_valid(canonical_bulk_study("in_house"))

    def test_valid_independent_dnt_evidence(self):
        document = canonical_bulk_study()
        document["dnt_reference_evidence"].append({
            "evidence_id": "evidence_1",
            "agent_name": reported("Compound A"),
            "agent_identifier": reported("example:compound-a"),
            "reference_state": "positive",
            "evidence_source": "independent curated reference",
            "evidence_version": "1",
            "evidence_reference": "PMID:00000001",
            "scope_type": "agent_general",
            "scope_assertions": [],
            "metadata_validation_status": "not_validated",
        })
        self.assert_valid(with_provenance(document))

    def test_broken_provenance_reference(self):
        document = canonical_bulk_study()
        document["conditions"][0]["condition_label"]["assertion_provenance_ids"] = ["missing_provenance"]
        self.assert_rule(document, "provenance_reference")

    def test_reported_value_requires_original_value(self):
        document = canonical_bulk_study()
        del document["conditions"][0]["condition_label"]["original_value"]
        self.assert_rule(document, "reported_value_required")

    def test_value_record_rejects_unknown_nested_field(self):
        document = canonical_bulk_study()
        document["conditions"][0]["condition_label"]["inferred_value"] = "not allowed"
        self.assert_rule(document, "unknown_field")

    def test_valid_sample_comparison_exception(self):
        document = canonical_bulk_study()
        document["sample_comparison_exceptions"].append({
            "exception_id": "exception_1",
            "relationship_id": "comparison_1",
            "affected_sample_ids": ["sample_t1"],
            "replacement_control_sample_ids": ["sample_c1"],
            "exception_reason": reported("Documented sample-specific pairing"),
            "metadata_validation_status": "not_validated",
        })
        self.assert_valid(with_provenance(document))

    def test_invalid_sample_comparison_exception_scope(self):
        document = canonical_bulk_study()
        document["conditions"].append({
            "condition_id": "unrelated_control",
            "experiment_id": "experiment_1",
            "condition_label": reported("Unrelated control"),
            "treatment_control_status": reported("control"),
            "control_type": reported("untreated"),
            "metadata_validation_status": "not_validated",
        })
        document["samples"].append(dict(document["samples"][0], sample_id="unrelated_sample", condition_id="unrelated_control"))
        document["sample_comparison_exceptions"].append({
            "exception_id": "exception_1",
            "relationship_id": "comparison_1",
            "affected_sample_ids": ["unrelated_sample"],
            "exception_reason": reported("Invalid scope"),
            "metadata_validation_status": "not_validated",
        })
        self.assert_rule(document, "comparison_exception_scope")

    def test_invalid_dnt_reference_state(self):
        document = canonical_bulk_study()
        document["dnt_reference_evidence"].append({
            "evidence_id": "evidence_1",
            "agent_name": reported("Compound A"),
            "reference_state": "inferred_from_expression",
            "evidence_source": "invalid",
            "evidence_version": "1",
            "scope_type": "unspecified",
            "scope_assertions": [],
            "metadata_validation_status": "not_validated",
        })
        self.assert_rule(document, "controlled_vocabulary")

    def test_primary_entity_identifiers_are_unique(self):
        document = canonical_bulk_study()
        document["sample_comparison_exceptions"].append({
            "exception_id": "exception_1",
            "relationship_id": "comparison_1",
            "affected_sample_ids": ["sample_t1"],
            "exception_reason": reported("Documented exception"),
            "metadata_validation_status": "not_validated",
        })
        document["dnt_reference_evidence"].append({
            "evidence_id": "evidence_1",
            "agent_name": reported("Compound A"),
            "reference_state": "unresolved_uncertain",
            "evidence_source": "independent source",
            "evidence_version": "1",
            "scope_type": "unspecified",
            "scope_assertions": [],
            "metadata_validation_status": "not_validated",
        })
        document = with_provenance(document)
        id_fields = {
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
        for entity, id_field in id_fields.items():
            with self.subTest(entity=entity, id_field=id_field):
                duplicate = deepcopy(document)
                duplicate[entity].append(deepcopy(duplicate[entity][0]))
                self.assert_rule(duplicate, "identifier_unique")

    def test_default_schema_path_is_independent_of_working_directory(self):
        original_directory = os.getcwd()
        with tempfile.TemporaryDirectory() as other_directory:
            try:
                os.chdir(other_directory)
                validator = MetadataValidator()
                result = validator.validate(canonical_bulk_study())
            finally:
                os.chdir(original_directory)
        self.assertEqual(ValidationStatus.VALID, result.status, result.findings)

    def test_findings_are_deterministic(self):
        document = canonical_bulk_study()
        document["samples"][0]["condition_id"] = "missing_condition"
        document["samples"][1]["developmental_age_or_stage_at_collection"] = missing("unknown")
        first = self.validator.validate(document)
        second = self.validator.validate(document)
        self.assertEqual(first, second)
        self.assertEqual(first.to_dict(), second.to_dict())

    def test_validation_status_precedence(self):
        self.assertEqual(ValidationStatus.VALID, self.validator.validate(canonical_bulk_study()).status)
        warning_document = canonical_bulk_study()
        warning_document["samples"][0]["developmental_age_or_stage_at_collection"] = missing()
        self.assertEqual(ValidationStatus.VALID_WITH_WARNINGS, self.validator.validate(warning_document).status)
        invalid_document = deepcopy(warning_document)
        invalid_document["samples"][0]["condition_id"] = "missing_condition"
        self.assertEqual(ValidationStatus.INVALID, self.validator.validate(invalid_document).status)

    def test_validation_does_not_mutate_input(self):
        document = canonical_bulk_study()
        before = deepcopy(document)
        self.validator.validate(document)
        self.assertEqual(before, document)

    def test_cross_experiment_condition_comparison_is_invalid(self):
        document = canonical_bulk_study()
        document["experiments"].append({
            "experiment_id": "experiment_2",
            "study_id": "study_1",
            "metadata_validation_status": "not_validated",
        })
        document["conditions"][0]["experiment_id"] = "experiment_2"
        self.assert_rule(document, "comparison_experiment")

    def test_duplicate_provenance_id_makes_reference_ambiguous_and_invalid(self):
        document = canonical_bulk_study()
        document["provenance"].append(deepcopy(document["provenance"][0]))
        self.assert_rule(document, "identifier_unique")

    def test_missing_values_are_not_inferred(self):
        document = canonical_bulk_study()
        missing_value = missing("unknown")
        document["exposures"][1]["concentration_or_dose"] = missing_value
        before = deepcopy(document)
        result = self.validator.validate(document)
        self.assertEqual(ValidationStatus.VALID_WITH_WARNINGS, result.status)
        self.assertEqual(before, document)
        self.assertEqual("unknown", document["exposures"][1]["concentration_or_dose"]["value_status"])
        self.assertIsNone(document["exposures"][1]["concentration_or_dose"]["original_value"])


if __name__ == "__main__":
    unittest.main()
