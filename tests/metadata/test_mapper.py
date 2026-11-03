"""Unit tests for CanonicalMetadataMapper and schema compliance."""

import pytest

from src.metadata.adapters.base import RawSampleRecord, StudyContext
from src.metadata.mapper import CanonicalMetadataMapper
from src.metadata.validator import MetadataValidator, ValidationStatus


def _create_synthetic_raw_records() -> tuple[StudyContext, list[RawSampleRecord]]:
    ctx = StudyContext(
        study_id="study_TEST_001",
        study_title="Synthetic Human Cortical Organoid VPA Exposure Study",
        source_type="public",
        public_accessions=({"namespace": "GSE", "identifier": "GSE999999", "url": None},),
        experiment_id="exp_01",
    )

    recs = [
        # Control replicates (3)
        RawSampleRecord(
            sample_id="sample_c1",
            attributes={
                "condition": "control",
                "control_type": "vehicle",
                "compound": "DMSO",
                "dose": "0.0",
                "dose_unit": "uM",
                "duration": "48",
                "duration_unit": "hours",
                "stage": "Day 45",
                "species": "Homo sapiens",
                "cell_line": "WTC11",
                "organoid": "cerebral organoid",
                "region": "dorsal forebrain",
                "replicate": "rep1",
            },
            source_file="test_sheet.csv",
            source_locator="row_1",
            source_key_provenance={"dose": "dose_col", "compound": "agent_col"},
        ),
        RawSampleRecord(
            sample_id="sample_c2",
            attributes={
                "condition": "control",
                "control_type": "vehicle",
                "compound": "DMSO",
                "dose": "0.0",
                "dose_unit": "uM",
                "duration": "48",
                "duration_unit": "hours",
                "stage": "Day 45",
                "species": "Homo sapiens",
                "cell_line": "WTC11",
                "organoid": "cerebral organoid",
                "region": "dorsal forebrain",
                "replicate": "rep2",
            },
            source_file="test_sheet.csv",
            source_locator="row_2",
            source_key_provenance={"dose": "dose_col", "compound": "agent_col"},
        ),
        RawSampleRecord(
            sample_id="sample_c3",
            attributes={
                "condition": "control",
                "control_type": "vehicle",
                "compound": "DMSO",
                "dose": "0.0",
                "dose_unit": "uM",
                "duration": "48",
                "duration_unit": "hours",
                "stage": "Day 45",
                "species": "Homo sapiens",
                "cell_line": "WTC11",
                "organoid": "cerebral organoid",
                "region": "dorsal forebrain",
                "replicate": "rep3",
            },
            source_file="test_sheet.csv",
            source_locator="row_3",
            source_key_provenance={"dose": "dose_col", "compound": "agent_col"},
        ),
        # Treatment replicates (3)
        RawSampleRecord(
            sample_id="sample_t1",
            attributes={
                "condition": "treatment",
                "compound": "Valproic Acid",
                "dose": "0.2",
                "dose_unit": "mM",  # 0.2 mM = 200 uM
                "duration": "2",
                "duration_unit": "days",  # 2 days = 48 h
                "stage": "Day 45",
                "species": "Homo sapiens",
                "cell_line": "WTC11",
                "organoid": "cerebral organoid",
                "region": "dorsal forebrain",
                "replicate": "rep1",
            },
            source_file="test_sheet.csv",
            source_locator="row_4",
            source_key_provenance={"dose": "dose_col", "compound": "agent_col"},
        ),
        RawSampleRecord(
            sample_id="sample_t2",
            attributes={
                "condition": "treatment",
                "compound": "Valproic Acid",
                "dose": "0.2",
                "dose_unit": "mM",
                "duration": "2",
                "duration_unit": "days",
                "stage": "Day 45",
                "species": "Homo sapiens",
                "cell_line": "WTC11",
                "organoid": "cerebral organoid",
                "region": "dorsal forebrain",
                "replicate": "rep2",
            },
            source_file="test_sheet.csv",
            source_locator="row_5",
            source_key_provenance={"dose": "dose_col", "compound": "agent_col"},
        ),
        RawSampleRecord(
            sample_id="sample_t3",
            attributes={
                "condition": "treatment",
                "compound": "Valproic Acid",
                "dose": "0.2",
                "dose_unit": "mM",
                "duration": "2",
                "duration_unit": "days",
                "stage": "Day 45",
                "species": "Homo sapiens",
                "cell_line": "WTC11",
                "organoid": "cerebral organoid",
                "region": "dorsal forebrain",
                "replicate": "rep3",
            },
            source_file="test_sheet.csv",
            source_locator="row_6",
            source_key_provenance={"dose": "dose_col", "compound": "agent_col"},
        ),
    ]
    return ctx, recs


class TestCanonicalMetadataMapper:
    def test_canonical_mapping_and_validator_pass(self):
        ctx, recs = _create_synthetic_raw_records()
        mapper = CanonicalMetadataMapper(target_dose_unit="uM", target_duration_unit="h")
        canonical_doc, findings = mapper.map_study(ctx, recs, source_reference="test_sheet.csv")

        assert len(canonical_doc["samples"]) == 6
        assert len(canonical_doc["conditions"]) == 2
        assert len(canonical_doc["treatment_control_relationships"]) == 1

        # Check unit normalization and original value preservation
        trt_exp = next(e for e in canonical_doc["exposures"] if "treatment" in e["condition_id"] or "valproic" in e["condition_id"])
        assert trt_exp["concentration_or_dose"]["original_value"] == "0.2"
        assert trt_exp["concentration_or_dose"]["normalized_value"] == pytest.approx(200.0)
        assert trt_exp["concentration_or_dose_unit"]["original_value"] == "mM"
        assert trt_exp["concentration_or_dose_unit"]["normalized_value"] == "uM"

        assert trt_exp["exposure_duration"]["original_value"] == "2"
        assert trt_exp["exposure_duration"]["normalized_value"] == pytest.approx(48.0)
        assert trt_exp["exposure_duration_unit"]["original_value"] == "days"
        assert trt_exp["exposure_duration_unit"]["normalized_value"] == "h"

        # Check field-level provenance
        assert len(canonical_doc["provenance"]) > 0
        prov_map = {p["provenance_id"]: p for p in canonical_doc["provenance"]}
        for p_id in trt_exp["concentration_or_dose"]["assertion_provenance_ids"]:
            assert p_id in prov_map
            assert prov_map[p_id]["field_path"] == "concentration_or_dose"
            assert prov_map[p_id]["source_reference"] == "test_sheet.csv"

        # Validate with the strict MetadataValidator
        validator = MetadataValidator()
        result = validator.validate(canonical_doc)
        error_findings = [f for f in result.findings if f.severity.value == "ERROR" or f.severity == "ERROR"]
        assert not error_findings, f"Validator found errors: {error_findings}"
        assert result.status in (ValidationStatus.VALID, ValidationStatus.VALID_WITH_WARNINGS)
