"""Unit tests for CanonicalMetadataMapper and schema compliance."""

from copy import deepcopy
from dataclasses import replace
import random
import pytest

from src.differential_expression.pipeline import run_contrast_differential_expression
from src.metadata.adapters.base import RawSampleRecord, StudyContext
from src.metadata.mapper import CanonicalMetadataMapper
from src.metadata.validator import MetadataValidator, ValidationStatus
from src.normalization.models import NormalizedCountMatrix, NormalizedDataset
from src.qc.matrix import CountMatrix
from src.response_builder.builder import build_response_contrasts
from src.response_builder.models import ContrastStatus
from tests.differential_expression.conftest import create_synthetic_de_fixtures
from tests.response_builder.test_builder import create_normalized_dataset_and_metadata


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

    def test_unmatched_treatment_condition_due_to_stage_mismatch(self):
        """Treatment Day 30 vs Control Day 60 cannot be automatically paired; must halt without guessing."""
        ctx = StudyContext(study_id="study_mismatch", source_type="public", experiment_id="exp_01")
        recs = [
            RawSampleRecord(
                sample_id="ctrl_d60_1",
                attributes={"condition": "control", "dose": "0", "stage": "Day 60", "cell_line": "WTC11"},
                source_file="meta.csv",
                source_locator="row_1",
            ),
            RawSampleRecord(
                sample_id="trt_d30_1",
                attributes={"condition": "treatment", "compound": "VPA", "dose": "100", "dose_unit": "uM", "stage": "Day 30", "cell_line": "WTC11"},
                source_file="meta.csv",
                source_locator="row_2",
            ),
        ]
        mapper = CanonicalMetadataMapper()
        canonical_doc, findings = mapper.map_study(ctx, recs, source_reference="meta.csv")

        # Must report unmatched_treatment_condition ERROR finding and refuse to guess
        assert any(f.rule == "unmatched_treatment_condition" for f in findings)
        assert len(canonical_doc["treatment_control_relationships"]) == 0

    def test_ambiguous_matched_controls_remain_unlinked_and_schema_compliant(self):
        """Multiple candidate controls must NOT be guessed or pooled; relationship remains unlinked and schema-valid."""
        ctx = StudyContext(study_id="study_ambig", source_type="public", experiment_id="exp_01")
        recs = [
            # Two distinct controls: vehicle DMSO and untreated
            RawSampleRecord(
                sample_id="ctrl_veh_1",
                attributes={"condition": "control", "control_type": "vehicle", "compound": "DMSO", "stage": "Day 45", "cell_line": "WTC11"},
                source_file="meta.csv",
                source_locator="row_1",
            ),
            RawSampleRecord(
                sample_id="ctrl_untreat_1",
                attributes={"condition": "control", "control_type": "untreated", "compound": "untreated", "stage": "Day 45", "cell_line": "WTC11"},
                source_file="meta.csv",
                source_locator="row_2",
            ),
            RawSampleRecord(
                sample_id="trt_1",
                attributes={"condition": "treatment", "compound": "VPA", "dose": "100", "dose_unit": "uM", "stage": "Day 45", "cell_line": "WTC11"},
                source_file="meta.csv",
                source_locator="row_3",
            ),
        ]
        mapper = CanonicalMetadataMapper()
        canonical_doc, findings = mapper.map_study(ctx, recs, source_reference="meta.csv")

        # Must emit ambiguous_matched_control WARNING
        assert any(f.rule == "ambiguous_matched_control" for f in findings)
        # Relationship must NOT be formed automatically with guessed or pooled controls
        assert len(canonical_doc["treatment_control_relationships"]) == 0

        # Validate with strict MetadataValidator: must pass without vocabulary errors
        validator = MetadataValidator()
        result = validator.validate(canonical_doc)
        error_findings = [f for f in result.findings if f.severity.value == "ERROR" or f.severity == "ERROR"]
        assert not error_findings, f"Validator found errors: {error_findings}"
        assert result.status in (ValidationStatus.VALID, ValidationStatus.VALID_WITH_WARNINGS)

    def test_ambiguous_controls_blocked_downstream_until_explicit_resolution(self):
        """Ambiguous controls are isolated and blocked from DESeq2 until human review / evidence resolves them."""
        ctx = StudyContext(study_id="study_ambig_pipe", source_type="public", experiment_id="exp_01")
        recs = [
            RawSampleRecord("ctrl_veh_1", {"condition": "control", "control_type": "vehicle", "compound": "DMSO", "stage": "Day 45", "cell_line": "WTC11", "rep": "r1"}, "f", "1"),
            RawSampleRecord("ctrl_veh_2", {"condition": "control", "control_type": "vehicle", "compound": "DMSO", "stage": "Day 45", "cell_line": "WTC11", "rep": "r2"}, "f", "2"),
            RawSampleRecord("ctrl_untreat_1", {"condition": "control", "control_type": "untreated", "compound": "untreated", "stage": "Day 45", "cell_line": "WTC11", "rep": "r1"}, "f", "3"),
            RawSampleRecord("ctrl_untreat_2", {"condition": "control", "control_type": "untreated", "compound": "untreated", "stage": "Day 45", "cell_line": "WTC11", "rep": "r2"}, "f", "4"),
            RawSampleRecord("trt_1", {"condition": "treatment", "compound": "VPA", "dose": "100", "dose_unit": "uM", "stage": "Day 45", "cell_line": "WTC11", "rep": "r1"}, "f", "5"),
            RawSampleRecord("trt_2", {"condition": "treatment", "compound": "VPA", "dose": "100", "dose_unit": "uM", "stage": "Day 45", "cell_line": "WTC11", "rep": "r2"}, "f", "6"),
        ]
        mapper = CanonicalMetadataMapper()
        doc, map_findings = mapper.map_study(ctx, recs, source_reference="meta.csv")

        # 1. Ambiguous match detected by mapper
        assert any(f.rule == "ambiguous_matched_control" for f in map_findings)
        assert len(doc["treatment_control_relationships"]) == 0

        # Build normalized dataset with realistic count distributions
        _, _, harm_ds, _ = create_synthetic_de_fixtures(n_treatment=2, n_control=4, n_genes=25)
        norm_ds, _ = create_normalized_dataset_and_metadata(
            samples_spec=doc["samples"],
            conditions_spec=doc["conditions"],
            relationships_spec=doc["treatment_control_relationships"],
            study_id=ctx.study_id,
            experiment_id=ctx.experiment_id,
        )

        g_ids = tuple(g.original_gene_id for g in harm_ds.genes)
        s_ids = norm_ds.cohort_results[0].cohort.sample_ids

        mean_levels = [50, 80, 120, 200, 350, 500, 800, 1200]
        counts_by_sample = {s: [] for s in s_ids}
        for i, gid in enumerate(g_ids):
            mu = mean_levels[i % len(mean_levels)]
            gene_rng = random.Random(1000 + i)
            for s in s_ids:
                lam = gene_rng.gammavariate(10, mu / 10)
                val = max(1, int(round(gene_rng.gauss(lam, (mu / 10) ** 0.5))))
                counts_by_sample[s].append(val)

        raw_cols = tuple(tuple(counts_by_sample[s]) for s in s_ids)
        raw_matrix = CountMatrix(gene_ids=g_ids, sample_ids=s_ids, columns=raw_cols)
        norm_matrix = NormalizedCountMatrix(gene_ids=g_ids, sample_ids=s_ids, columns=tuple(tuple(float(v) for v in col) for col in raw_cols))

        cr = norm_ds.cohort_results[0]
        cr_updated = replace(cr, raw_count_matrix=raw_matrix, diagnostic_normalized_matrix=norm_matrix)
        norm_ds = replace(norm_ds, cohort_results=(cr_updated,))

        # 2. Step 7A: Without explicit resolution, contrast is unlinked and BLOCKED
        contrast_ds = build_response_contrasts(norm_ds, doc)
        assert len(contrast_ds.contrasts) == 1
        blocked_contrast = contrast_ds.contrasts[0]
        assert blocked_contrast.status == ContrastStatus.BLOCKED
        assert blocked_contrast.matched_control_condition_ids == ()
        assert blocked_contrast.control_sample_ids == ()
        assert any(f.rule_id == "missing_explicit_treatment_control_relationship" for f in blocked_contrast.findings)

        # 3. Step 7B: Contrast skips DESeq2 execution upstream
        de_res_blocked = run_contrast_differential_expression(blocked_contrast, norm_ds, harm_ds, doc)
        assert de_res_blocked.status == ContrastStatus.BLOCKED
        assert "Skipped upstream: contrast status is BLOCKED" in de_res_blocked.design_rationale
        assert len(de_res_blocked.gene_results) == 0

        # 4. Human Review / Evidence Resolution
        # Curator confirms vehicle DMSO is the intended matched reference control per protocol
        veh_cond_id = "cond_exp_01_control_vehicle_dmso_day_45"
        trt_cond_id = "cond_exp_01_vpa_100.0uM"

        doc_resolved = deepcopy(doc)
        doc_resolved["treatment_control_relationships"].append({
            "relationship_id": f"rel_{trt_cond_id}_vs_{veh_cond_id}",
            "treatment_condition_ids": [trt_cond_id],
            "matched_control_condition_ids": [veh_cond_id],
            "matching_basis": {
                "original_value": "Curator confirmed vehicle DMSO matched control per experimental design",
                "normalized_value": "Curator confirmed vehicle DMSO matched control per experimental design",
                "value_status": "reported",
                "assertion_provenance_ids": ["prov_curator_1"],
            },
            "metadata_validation_status": "valid",
        })
        doc_resolved["provenance"].append({
            "provenance_id": "prov_curator_1",
            "entity_type": "treatment_control_relationships",
            "entity_id": f"rel_{trt_cond_id}_vs_{veh_cond_id}",
            "field_path": "matched_control_condition_ids",
            "source_kind": "curator_assertion",
            "source_reference": "experimental_design_curation",
        })

        norm_ds_resolved, _ = create_normalized_dataset_and_metadata(
            samples_spec=doc_resolved["samples"],
            conditions_spec=doc_resolved["conditions"],
            relationships_spec=doc_resolved["treatment_control_relationships"],
            study_id=ctx.study_id,
            experiment_id=ctx.experiment_id,
        )
        cr_res = norm_ds_resolved.cohort_results[0]
        cr_res_updated = replace(cr_res, raw_count_matrix=raw_matrix, diagnostic_normalized_matrix=norm_matrix)
        norm_ds_resolved = replace(norm_ds_resolved, cohort_results=(cr_res_updated,))

        # 5. Re-evaluate Step 7A: now ELIGIBLE against the single confirmed matched control
        contrast_ds_resolved = build_response_contrasts(norm_ds_resolved, doc_resolved)
        assert len(contrast_ds_resolved.contrasts) == 1
        resolved_contrast = contrast_ds_resolved.contrasts[0]
        assert resolved_contrast.status == ContrastStatus.ELIGIBLE
        assert resolved_contrast.matched_control_condition_ids == (veh_cond_id,)
        assert set(resolved_contrast.control_sample_ids) == {"ctrl_veh_1", "ctrl_veh_2"}

        # 6. Re-evaluate Step 7B: DESeq2 executes successfully
        de_res_resolved = run_contrast_differential_expression(resolved_contrast, norm_ds_resolved, harm_ds, doc_resolved)
        assert de_res_resolved.status == ContrastStatus.ELIGIBLE
        assert len(de_res_resolved.gene_results) == 25
        assert de_res_resolved.design_formula == "~ condition"
