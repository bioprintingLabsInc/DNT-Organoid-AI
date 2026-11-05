"""Evaluation of real-world public RNA-seq study metadata (E-MTAB-15191, PRJNA802072, GSE232218)."""

from pathlib import Path
import pytest

from src.metadata.adapters import (
    GeoMetadataAdapter,
    JsonMetadataAdapter,
    MageTabAdapter,
)
from src.metadata.integrity import MetadataIntegrityChecker
from src.metadata.mapper import CanonicalMetadataMapper
from src.metadata.validator import MetadataValidator, ValidationStatus


class TestRealPublicStudiesMetadata:
    def test_real_emtab_15191_magetab_pipeline(self):
        """Evaluate real E-MTAB-15191 MAGE-TAB files from BioStudies / ArrayExpress.

        Study: Human dorsal forebrain organoids treated with valproic acid vs untreated controls.
        Total biological samples: 16 (8 treated, 8 control; 32 rows representing paired reads).
        """
        fixture_dir = Path(__file__).resolve().parents[1] / "fixtures"
        sdrf_path = fixture_dir / "E-MTAB-15191.sdrf.txt"
        idf_path = fixture_dir / "E-MTAB-15191.idf.txt"

        assert sdrf_path.exists(), "E-MTAB-15191.sdrf.txt missing from fixtures"
        assert idf_path.exists(), "E-MTAB-15191.idf.txt missing from fixtures"

        adapter = MageTabAdapter()
        ctx, recs, findings = adapter.parse_study(sdrf_path, idf_source=idf_path)

        # 16 unique biological samples consolidated from 32 rows
        assert len(recs) == 16
        assert ctx.study_id == "study_E-MTAB-15191"

        # Map to canonical schema
        mapper = CanonicalMetadataMapper(
            target_dose_unit="uM",
            target_duration_unit="h",
            custom_mapping_rules={"biological_replicate_id": ("Source Name", "source_name")},
        )
        canonical_doc, map_findings = mapper.map_study(ctx, recs, source_reference="E-MTAB-15191.sdrf.txt")

        assert len(canonical_doc["samples"]) == 16
        assert len(canonical_doc["conditions"]) == 2  # 1 treatment (VPA 50 mM), 1 control (untreated)

        # Verify unit normalization: 50 millimolar = 50,000 uM
        trt_exp = next(e for e in canonical_doc["exposures"] if "valproic" in e["condition_id"])
        assert trt_exp["concentration_or_dose"]["original_value"] == "50"
        assert trt_exp["concentration_or_dose"]["normalized_value"] == pytest.approx(50000.0)
        assert trt_exp["concentration_or_dose_unit"]["normalized_value"] == "uM"

        # Validate with strict schema validator
        validator = MetadataValidator()
        val_res = validator.validate(canonical_doc)
        error_findings = [f for f in val_res.findings if f.severity.value == "ERROR" or f.severity == "ERROR"]
        assert not error_findings, f"Validation errors found: {error_findings}"
        assert val_res.status in (ValidationStatus.VALID, ValidationStatus.VALID_WITH_WARNINGS)

        # Run integrity checks
        checker = MetadataIntegrityChecker(min_biological_replicates=2)
        integ_res = checker.check(canonical_doc)
        assert integ_res.is_valid

    def test_real_prjna802072_pilot_manifest(self):
        """Evaluate real PRJNA802072 pilot manifest from quantification bridge.

        Study: Day 28 cortical organoids exposed to VPA 250 uM vs 0 uM vehicle control.
        Total samples: 12 (6 treated, 6 control).
        """
        repo_root = Path(__file__).resolve().parents[2]
        manifest_path = repo_root / "config" / "quantification_bridge" / "pilot_manifest.json"

        assert manifest_path.exists(), "pilot_manifest.json missing"

        adapter = JsonMetadataAdapter()
        ctx, recs, findings = adapter.parse_study(manifest_path)

        assert len(recs) == 12
        assert ctx.study_id == "PRJNA802072"

        mapper = CanonicalMetadataMapper(
            target_dose_unit="uM",
            target_duration_unit="h",
            defaults={"organoid_type": "cerebral organoid", "brain_region": "cerebral cortex"},
        )
        canonical_doc, map_findings = mapper.map_study(ctx, recs, source_reference="pilot_manifest.json")

        assert len(canonical_doc["samples"]) == 12
        assert len(canonical_doc["conditions"]) == 2

        validator = MetadataValidator()
        val_res = validator.validate(canonical_doc)
        error_findings = [f for f in val_res.findings if f.severity.value == "ERROR" or f.severity == "ERROR"]
        assert not error_findings, f"Validation errors found: {error_findings}"
        assert val_res.status in (ValidationStatus.VALID, ValidationStatus.VALID_WITH_WARNINGS)
