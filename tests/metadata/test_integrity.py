"""Unit tests for MetadataIntegrityChecker (treatment-control, biological models, replicates)."""

import pytest

from tests.metadata.fixtures import canonical_bulk_study, reported
from src.metadata.integrity import MetadataIntegrityChecker
from src.metadata.result import Severity, ValidationStatus


class TestMetadataIntegrityChecker:
    def test_valid_study_passes_integrity_check(self):
        doc = canonical_bulk_study()
        checker = MetadataIntegrityChecker(min_biological_replicates=2)
        res = checker.check(doc)
        assert res.is_valid
        assert not any(f.severity == Severity.ERROR for f in res.findings)

    def test_duplicate_sample_id_detected(self):
        doc = canonical_bulk_study()
        # Duplicate sample_c1
        doc["samples"].append(dict(doc["samples"][0]))
        checker = MetadataIntegrityChecker()
        res = checker.check(doc)
        assert not res.is_valid
        assert any(f.rule == "duplicate_sample_id" for f in res.findings)

    def test_missing_biological_replicate_id_detected(self):
        doc = canonical_bulk_study()
        # Set biological_replicate_id to empty
        doc["samples"][0]["biological_replicate_id"] = reported("")
        checker = MetadataIntegrityChecker()
        res = checker.check(doc)
        assert not res.is_valid
        assert any(f.rule == "missing_biological_replicate_id" for f in res.findings)

    def test_technical_replicate_inflation_detected(self):
        doc = canonical_bulk_study()
        # Make all 3 treatment samples share biological_replicate_id 'rep1'
        for s in doc["samples"]:
            if "t" in s["sample_id"]:
                s["biological_replicate_id"] = reported("rep1")
        checker = MetadataIntegrityChecker(min_biological_replicates=2)
        res = checker.check(doc)
        # Should flag technical replicates detected and under-replicated cohort
        assert any(f.rule == "technical_replicates_detected" for f in res.findings)
        assert any(f.rule == "under_replicated_cohort" for f in res.findings)

    def test_missing_matched_control_detected(self):
        doc = canonical_bulk_study()
        doc["treatment_control_relationships"][0]["matched_control_condition_ids"] = []
        checker = MetadataIntegrityChecker()
        res = checker.check(doc)
        assert not res.is_valid
        assert any(f.rule == "missing_matched_control" for f in res.findings)

    def test_cross_species_comparison_prohibited(self):
        doc = canonical_bulk_study()
        # Add rat biological source and assign control samples to it
        doc["biological_sources"].append({
            "biological_source_id": "rat_line_1",
            "species": reported("Rattus norvegicus"),
            "source_material_type": reported("other"),
            "cell_line_or_source_identifier": reported("rat_cns"),
            "metadata_validation_status": "not_validated",
        })
        for s in doc["samples"]:
            if "c" in s["sample_id"]:
                s["biological_source_id"] = "rat_line_1"

        checker = MetadataIntegrityChecker()
        res = checker.check(doc)
        assert not res.is_valid
        assert any(f.rule == "cross_species_comparison_prohibited" for f in res.findings)

    def test_treatment_control_role_inconsistency_detected(self):
        doc = canonical_bulk_study()
        # Invert roles: assign a control condition into treatment_condition_ids
        ctrl_cid = doc["conditions"][0]["condition_id"]
        doc["treatment_control_relationships"][0]["treatment_condition_ids"] = [ctrl_cid]

        checker = MetadataIntegrityChecker()
        res = checker.check(doc)
        assert not res.is_valid
        assert any(f.rule == "treatment_control_role_inconsistency" for f in res.findings)

    def test_developmental_stage_mismatch_detected(self):
        doc = canonical_bulk_study()
        # Set control samples to Day 60 and treatment samples to Day 30
        for s in doc["samples"]:
            if "c" in s["sample_id"]:
                s["developmental_age_or_stage_at_collection"] = reported("Day 60")
            elif "t" in s["sample_id"]:
                s["developmental_age_or_stage_at_collection"] = reported("Day 30")

        checker = MetadataIntegrityChecker()
        res = checker.check(doc)
        assert not res.is_valid
        assert any(f.rule == "developmental_stage_mismatch" for f in res.findings)

    def test_ambiguous_matched_controls_detected(self):
        doc = canonical_bulk_study()
        # Add a second control condition and assign both to relationship with needs_review
        doc["conditions"].append(dict(doc["conditions"][0], condition_id="cond_ctrl_2"))
        doc["treatment_control_relationships"][0]["matched_control_condition_ids"].append("cond_ctrl_2")
        doc["treatment_control_relationships"][0]["metadata_validation_status"] = "needs_review"

        checker = MetadataIntegrityChecker()
        res = checker.check(doc)
        assert any(f.rule == "ambiguous_matched_control" for f in res.findings)
        assert any(f.rule == "relationship_requires_review" for f in res.findings)
