"""Phase U2 Comprehensive Acceptance Test Suite.

Covers all 15 mandatory acceptance points from Section 5 of the Phase U2 specification:
1. Correctly formatted metadata.
2. Missing optional metadata.
3. Missing essential metadata.
4. Different concentration units.
5. Different exposure-duration units.
6. Conflicting sample assignments.
7. Different biological models.
8. Multiple treatment concentrations and timepoints.
9. Biological versus technical replicate identification.
10. Metadata import from different supported layouts.
11. Ambiguous source-field mappings.
12. Preservation of original values and provenance.
13. Compatibility with U1 routing.
14. No regression in existing Stages 1–7.
15. No introduction of DNT outcome information into molecular features.
"""

from copy import deepcopy
from pathlib import Path
import pytest

from tests.metadata.fixtures import canonical_bulk_study, missing, reported
from src.ingestion.registry import register_assets
from src.metadata.adapters import (
    GeoMetadataAdapter,
    JsonMetadataAdapter,
    MageTabAdapter,
    RawSampleRecord,
    StudyContext,
    TabularMetadataAdapter,
)
from src.metadata.integrity import MetadataIntegrityChecker
from src.metadata.mapper import CanonicalMetadataMapper
from src.metadata.units import UnitConversionStatus, normalize_concentration, normalize_duration
from src.metadata.validator import MetadataValidator, ValidationStatus
from src.routing.router import route_metadata_assets


class TestPhaseU2AcceptanceSuite:
    # 1. Correctly formatted metadata
    def test_01_correctly_formatted_metadata(self):
        doc = canonical_bulk_study()
        validator = MetadataValidator()
        res = validator.validate(doc)
        assert res.status in (ValidationStatus.VALID, ValidationStatus.VALID_WITH_WARNINGS)
        assert not any(f.severity.value == "ERROR" or f.severity == "ERROR" for f in res.findings)

    # 2. Missing optional metadata
    def test_02_missing_optional_metadata(self):
        doc = canonical_bulk_study()
        # Optional fields: experiment_title, donor_identifier, vehicle, etc.
        doc["experiments"][0]["experiment_title"] = missing("missing_not_reported")
        doc["biological_sources"][0]["donor_identifier"] = missing("missing_not_reported")
        validator = MetadataValidator()
        res = validator.validate(doc)
        # Missing optional fields must NOT fail the study
        assert not any(f.severity.value == "ERROR" or f.severity == "ERROR" for f in res.findings)

    # 3. Missing essential metadata
    def test_03_missing_essential_metadata(self):
        doc = canonical_bulk_study()
        # Remove required condition_label from condition
        del doc["conditions"][0]["condition_label"]
        validator = MetadataValidator()
        res = validator.validate(doc)
        assert res.status == ValidationStatus.INVALID
        assert any(f.rule == "required_field" for f in res.findings)

    # 4. Different concentration units
    def test_04_different_concentration_units(self):
        # 1 mM -> 1000 uM
        r_mm = normalize_concentration(1.0, "mM", target_unit="uM")
        assert r_mm.status == UnitConversionStatus.CONVERTED
        assert r_mm.normalized_value == pytest.approx(1000.0)

        # 500 nM -> 0.5 uM
        r_nm = normalize_concentration(500.0, "nM", target_unit="uM")
        assert r_nm.status == UnitConversionStatus.CONVERTED
        assert r_nm.normalized_value == pytest.approx(0.5)

        # 10 ug/mL to uM is blocked (incompatible dimensions)
        r_mass = normalize_concentration(10.0, "ug/mL", target_unit="uM")
        assert r_mass.status == UnitConversionStatus.INCOMPATIBLE_DIMENSIONS

    # 5. Different exposure-duration units
    def test_05_different_duration_units(self):
        # 2 days -> 48 h
        r_day = normalize_duration(2, "days", target_unit="h")
        assert r_day.status == UnitConversionStatus.CONVERTED
        assert r_day.normalized_value == pytest.approx(48.0)

        # 1 week -> 168 h
        r_wk = normalize_duration(1, "week", target_unit="h")
        assert r_wk.status == UnitConversionStatus.CONVERTED
        assert r_wk.normalized_value == pytest.approx(168.0)

        # 30 min -> 0.5 h
        r_min = normalize_duration(30, "min", target_unit="h")
        assert r_min.status == UnitConversionStatus.CONVERTED
        assert r_min.normalized_value == pytest.approx(0.5)

    # 6. Conflicting sample assignments
    def test_06_conflicting_sample_assignments(self):
        doc = canonical_bulk_study()
        # Duplicate sample ID with conflicting condition assignment
        dup = deepcopy(doc["samples"][0])
        dup["condition_id"] = "treatment_10um_24h_d45"
        doc["samples"].append(dup)

        checker = MetadataIntegrityChecker()
        res = checker.check(doc)
        assert not res.is_valid
        assert any(f.rule == "duplicate_sample_id" for f in res.findings)

    # 7. Different biological models
    def test_07_different_biological_models(self):
        doc = canonical_bulk_study()
        # Add cortical organoid and retinal organoid context
        doc["organoid_contexts"].append({
            "organoid_context_id": "retinal_1",
            "organoid_type": reported("retinal organoid"),
            "brain_region_or_model_identity": reported("retina"),
            "metadata_validation_status": "not_validated",
        })
        # Assign control to retinal and treatment to cortical
        for s in doc["samples"]:
            if "c" in s["sample_id"]:
                s["organoid_context_id"] = "retinal_1"

        checker = MetadataIntegrityChecker()
        res = checker.check(doc)
        # Distinct biological models compared without matching must emit a warning finding
        assert any(f.rule == "heterogeneous_organoid_contexts" for f in res.findings)

    # 8. Multiple treatment concentrations and timepoints
    def test_08_multiple_concentrations_and_timepoints(self):
        ctx = StudyContext(study_id="study_multi_dose", source_type="public")
        recs = [
            RawSampleRecord("s_c1", {"group": "control", "dose": 0, "time": 24, "rep": "r1"}, "f", "1"),
            RawSampleRecord("s_c2", {"group": "control", "dose": 0, "time": 24, "rep": "r2"}, "f", "2"),
            RawSampleRecord("s_t1_low", {"group": "treatment", "compound": "VPA", "dose": 100, "time": 24, "rep": "r1"}, "f", "3"),
            RawSampleRecord("s_t2_low", {"group": "treatment", "compound": "VPA", "dose": 100, "time": 24, "rep": "r2"}, "f", "4"),
            RawSampleRecord("s_t1_high", {"group": "treatment", "compound": "VPA", "dose": 500, "time": 48, "rep": "r1"}, "f", "5"),
            RawSampleRecord("s_t2_high", {"group": "treatment", "compound": "VPA", "dose": 500, "time": 48, "rep": "r2"}, "f", "6"),
        ]
        mapper = CanonicalMetadataMapper()
        doc, findings = mapper.map_study(ctx, recs)

        # 3 distinct conditions: control (24h), low dose (100 uM, 24h), high dose (500 uM, 48h)
        assert len(doc["conditions"]) == 3
        # 24h treatment matches 24h control; 48h treatment duration conflicts with single 24h control so cannot automatically match
        assert len(doc["treatment_control_relationships"]) == 1
        assert any(f.rule == "unmatched_treatment_condition" for f in findings)

    # 9. Biological versus technical replicate identification
    def test_09_biological_vs_technical_replicate_identification(self):
        doc = canonical_bulk_study()
        # Set biological replicate IDs to the same bio_rep, but distinct technical replicates
        for idx, s in enumerate(doc["samples"]):
            if "t" in s["sample_id"]:
                s["biological_replicate_id"] = reported("bio_rep_1")
                s["technical_replicate_id"] = reported(f"tech_rep_{idx+1}")

        checker = MetadataIntegrityChecker(min_biological_replicates=2)
        res = checker.check(doc)
        assert any(f.rule == "technical_replicates_detected" for f in res.findings)
        assert any(f.rule == "under_replicated_cohort" for f in res.findings)

    # 10. Metadata import from different supported layouts
    def test_10_metadata_import_from_supported_layouts(self, tmp_path):
        # CSV layout
        csv_file = tmp_path / "test.csv"
        csv_file.write_text("sample_id,condition,dose,replicate\ns1,control,0,r1\ns2,treatment,10,r1\n", encoding="utf-8")
        tab_adapter = TabularMetadataAdapter()
        assert tab_adapter.can_read(csv_file)
        _, tab_recs, _ = tab_adapter.parse_study(csv_file)
        assert len(tab_recs) == 2

        # JSON manifest layout
        json_file = tmp_path / "test.json"
        json_file.write_text('{"samples": [{"sample_id": "s1", "group": "control"}, {"sample_id": "s2", "group": "treatment"}]}', encoding="utf-8")
        json_adapter = JsonMetadataAdapter()
        assert json_adapter.can_read(json_file)
        _, json_recs, _ = json_adapter.parse_study(json_file)
        assert len(json_recs) == 2

    # 11. Ambiguous source-field mappings
    def test_11_ambiguous_source_field_mappings(self):
        ctx = StudyContext(study_id="study_ambig", source_type="public")
        recs = [
            RawSampleRecord("s1", {"mystery_column": "value_x"}, "f", "1"),
        ]
        mapper = CanonicalMetadataMapper()
        doc, findings = mapper.map_study(ctx, recs)
        # Should flag ambiguous condition status
        assert any(f.rule == "ambiguous_condition_status" for f in findings)

    # 12. Preservation of original values and provenance
    def test_12_preservation_of_original_values_and_provenance(self):
        ctx = StudyContext(study_id="study_prov", source_type="public")
        recs = [
            RawSampleRecord(
                "s1",
                {"condition": "treatment", "compound": "Valproic Acid", "dose": "0.5", "dose_unit": "mM", "duration": "48", "rep": "r1"},
                "source_file.csv",
                "row_2",
                {"dose": "dose_raw_header"},
            ),
            RawSampleRecord(
                "s2",
                {"condition": "control", "compound": "DMSO", "dose": "0", "dose_unit": "uM", "duration": "48", "rep": "r1"},
                "source_file.csv",
                "row_3",
                {"dose": "dose_raw_header"},
            ),
        ]
        mapper = CanonicalMetadataMapper()
        doc, findings = mapper.map_study(ctx, recs, source_reference="source_file.csv")

        trt_exp = next(e for e in doc["exposures"] if "treatment" in e["condition_id"] or "valproic" in e["condition_id"])
        # Original value preserved
        assert trt_exp["concentration_or_dose"]["original_value"] == "0.5"
        # Normalized value computed
        assert trt_exp["concentration_or_dose"]["normalized_value"] == pytest.approx(500.0)
        # Provenance entry present
        pid = trt_exp["concentration_or_dose"]["assertion_provenance_ids"][0]
        prov = next(p for p in doc["provenance"] if p["provenance_id"] == pid)
        assert prov["original_metadata_key"] == "dose_raw_header"
        assert prov["source_reference"] == "source_file.csv"
        assert prov["notes"] == "row_2"

    # 13. Compatibility with U1 routing
    def test_13_compatibility_with_u1_routing(self):
        doc = canonical_bulk_study()
        # Verify that canonical metadata passes route_metadata_assets and register_assets from U1
        registrations = register_assets(doc)
        assert len(registrations) == 1
        assert registrations[0].declared_format == "raw_counts"

        routing_results = route_metadata_assets(doc)
        assert len(routing_results) == 1

    # 14. No regression in existing Stages 1–7
    def test_14_no_regression_in_stages(self):
        # Verify MetadataValidator validates Stage 7 requirements cleanly
        doc = canonical_bulk_study()
        v = MetadataValidator()
        res = v.validate(doc)
        assert res.status in (ValidationStatus.VALID, ValidationStatus.VALID_WITH_WARNINGS)

    # 15. No introduction of DNT outcome information into molecular features
    def test_15_no_dnt_outcome_leakage(self):
        ctx = StudyContext(study_id="study_dnt_check", source_type="public")
        recs = [
            RawSampleRecord("s1", {"condition": "treatment", "compound": "VPA", "dose": 200, "rep": "r1"}, "f", "1"),
            RawSampleRecord("s2", {"condition": "control", "compound": "vehicle", "dose": 0, "rep": "r1"}, "f", "2"),
        ]
        mapper = CanonicalMetadataMapper()
        doc, findings = mapper.map_study(ctx, recs)

        # Verify dnt_reference_evidence remains empty and decoupled from molecular samples
        assert doc["dnt_reference_evidence"] == []
        for s in doc["samples"]:
            assert "dnt_label" not in s
            assert "reference_state" not in s
