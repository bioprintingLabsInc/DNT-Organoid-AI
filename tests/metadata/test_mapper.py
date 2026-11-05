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
            RawSampleRecord("ctrl_veh_1", {"condition": "control", "control_type": "vehicle", "compound": "DMSO", "stage": "Day 45", "cell_line": "WTC11", "organoid": "cerebral organoid", "region": "dorsal forebrain", "rep": "r1"}, "f", "1"),
            RawSampleRecord("ctrl_veh_2", {"condition": "control", "control_type": "vehicle", "compound": "DMSO", "stage": "Day 45", "cell_line": "WTC11", "organoid": "cerebral organoid", "region": "dorsal forebrain", "rep": "r2"}, "f", "2"),
            RawSampleRecord("ctrl_untreat_1", {"condition": "control", "control_type": "untreated", "compound": "untreated", "stage": "Day 45", "cell_line": "WTC11", "organoid": "cerebral organoid", "region": "dorsal forebrain", "rep": "r1"}, "f", "3"),
            RawSampleRecord("ctrl_untreat_2", {"condition": "control", "control_type": "untreated", "compound": "untreated", "stage": "Day 45", "cell_line": "WTC11", "organoid": "cerebral organoid", "region": "dorsal forebrain", "rep": "r2"}, "f", "4"),
            RawSampleRecord("trt_1", {"condition": "treatment", "compound": "VPA", "dose": "100", "dose_unit": "uM", "stage": "Day 45", "cell_line": "WTC11", "organoid": "cerebral organoid", "region": "dorsal forebrain", "rep": "r1"}, "f", "5"),
            RawSampleRecord("trt_2", {"condition": "treatment", "compound": "VPA", "dose": "100", "dose_unit": "uM", "stage": "Day 45", "cell_line": "WTC11", "organoid": "cerebral organoid", "region": "dorsal forebrain", "rep": "r2"}, "f", "6"),
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

    def test_no_unsupported_defaults_when_metadata_missing(self):
        """When attributes are omitted, mapper must record explicit missingness without inventing defaults."""
        ctx = StudyContext(study_id="study_incomplete", source_type="public", experiment_id="exp_01")
        recs = [
            RawSampleRecord(
                sample_id="s_ctrl",
                attributes={"condition": "control", "dose": 0},
                source_file="bare.csv",
                source_locator="row_1",
            ),
            RawSampleRecord(
                sample_id="s_trt",
                attributes={"condition": "treatment", "compound": "Compound_X", "dose": 50, "dose_unit": "uM"},
                source_file="bare.csv",
                source_locator="row_2",
            ),
        ]
        mapper = CanonicalMetadataMapper()
        doc, findings = mapper.map_study(ctx, recs, source_reference="bare.csv")

        # 1. Duration: must NOT default to 48.0 hours
        trt_exp = next(e for e in doc["exposures"] if "compound_x" in e["condition_id"])
        assert trt_exp["exposure_duration"]["value_status"] == "missing_not_reported"
        assert trt_exp["exposure_duration"]["original_value"] is None
        assert trt_exp["exposure_duration"]["normalized_value"] is None

        # 2. Developmental stage: must NOT default to Day 45
        assert trt_exp["developmental_age_or_stage_at_exposure"]["value_status"] == "unknown"
        assert trt_exp["developmental_age_or_stage_at_exposure"]["original_value"] is None

        trt_samp = next(s for s in doc["samples"] if s["sample_id"] == "s_trt")
        assert trt_samp["developmental_age_or_stage_at_collection"]["value_status"] == "unknown"
        assert trt_samp["developmental_age_or_stage_at_collection"]["original_value"] is None

        # 3. Cell line: must NOT default to cell_line_01
        bio_src = doc["biological_sources"][0]
        assert bio_src["cell_line_or_source_identifier"]["value_status"] == "unknown"
        assert bio_src["cell_line_or_source_identifier"]["original_value"] is None

        # 4. Organoid context: must NOT default to cerebral organoid / dorsal forebrain
        org_ctx = doc["organoid_contexts"][0]
        assert org_ctx["organoid_type"]["value_status"] == "unknown"
        assert org_ctx["organoid_type"]["original_value"] is None
        assert org_ctx["brain_region_or_model_identity"]["value_status"] == "unknown"
        assert org_ctx["brain_region_or_model_identity"]["original_value"] is None

        # 5. Biological replicate: must NOT default to rep_{sid}
        assert trt_samp["biological_replicate_id"]["value_status"] == "unknown"
        assert trt_samp["biological_replicate_id"]["original_value"] is None

        # 6. Schema validator must accept explicit missingness with zero errors
        validator = MetadataValidator()
        res = validator.validate(doc)
        error_findings = [f for f in res.findings if f.severity.value == "ERROR" or f.severity == "ERROR"]
        assert not error_findings, f"Validator found errors: {error_findings}"
        assert res.status in (ValidationStatus.VALID, ValidationStatus.VALID_WITH_WARNINGS)

    def test_missing_replicates_not_converted_to_artificial_independence(self):
        """Missing replicate metadata cannot deceive the integrity checker into seeing independent replicates."""
        from src.metadata.integrity import MetadataIntegrityChecker

        ctx = StudyContext(study_id="study_unverified_reps", source_type="public", experiment_id="exp_01")
        recs = [
            RawSampleRecord("c1", {"condition": "control", "dose": 0, "stage": "Day 28", "cell_line": "H9"}, "f", "1"),
            RawSampleRecord("c2", {"condition": "control", "dose": 0, "stage": "Day 28", "cell_line": "H9"}, "f", "2"),
            RawSampleRecord("t1", {"condition": "treatment", "compound": "VPA", "dose": 250, "dose_unit": "uM", "stage": "Day 28", "cell_line": "H9"}, "f", "3"),
            RawSampleRecord("t2", {"condition": "treatment", "compound": "VPA", "dose": 250, "dose_unit": "uM", "stage": "Day 28", "cell_line": "H9"}, "f", "4"),
        ]
        mapper = CanonicalMetadataMapper()
        doc, _ = mapper.map_study(ctx, recs, source_reference="meta.csv")

        # Samples must have biological_replicate_id marked unknown
        for s in doc["samples"]:
            assert s["biological_replicate_id"]["value_status"] == "unknown"
            assert s["biological_replicate_id"]["normalized_value"] is None

        # Integrity checker MUST catch missing replicates and reject independence
        checker = MetadataIntegrityChecker(min_biological_replicates=2)
        integ_res = checker.check(doc)
        assert not integ_res.is_valid
        missing_rep_findings = [f for f in integ_res.findings if f.rule == "missing_biological_replicate_id"]
        assert len(missing_rep_findings) == 4

    def test_separate_exposure_start_and_collection_ages(self):
        """Exposure start age and organoid collection age must remain separate fields."""
        ctx = StudyContext(study_id="study_age_sep", source_type="public", experiment_id="exp_01")
        recs = [
            RawSampleRecord(
                "s1",
                {
                    "condition": "control",
                    "dose": 0,
                    "exposure_start_day": "Day 16",
                    "harvest_day": "Day 28",
                    "cell_line": "H9",
                    "replicate": "r1",
                },
                "f", "1"
            ),
            RawSampleRecord(
                "s2",
                {
                    "condition": "treatment",
                    "compound": "VPA",
                    "dose": 250,
                    "dose_unit": "uM",
                    "exposure_start_day": "Day 16",
                    "harvest_day": "Day 28",
                    "cell_line": "H9",
                    "replicate": "r1",
                },
                "f", "2"
            ),
        ]
        mapper = CanonicalMetadataMapper()
        doc, findings = mapper.map_study(ctx, recs, source_reference="meta.csv")

        trt_exp = next(e for e in doc["exposures"] if "vpa" in e["condition_id"])
        assert trt_exp["developmental_age_or_stage_at_exposure"]["original_value"] == "Day 16"
        assert trt_exp["exposure_start_time_or_stage"]["original_value"] == "Day 16"

        trt_samp = next(s for s in doc["samples"] if s["sample_id"] == "s2")
        assert trt_samp["developmental_age_or_stage_at_collection"]["original_value"] == "Day 28"

        # Validator check
        validator = MetadataValidator()
        res = validator.validate(doc)
        assert not [f for f in res.findings if f.severity.value == "ERROR" or f.severity == "ERROR"]

    def test_samples_with_different_durations_or_stages_not_combined(self):
        """Samples with different durations or stages must not be combined into the same condition."""
        ctx = StudyContext(study_id="study_multicond", source_type="public", experiment_id="exp_01")
        recs = [
            # Control: 24h at Day 16
            RawSampleRecord("c1", {"condition": "control", "dose": 0, "duration": 24, "stage": "Day 16", "rep": "r1"}, "f", "1"),
            RawSampleRecord("c2", {"condition": "control", "dose": 0, "duration": 24, "stage": "Day 16", "rep": "r2"}, "f", "2"),
            # Treatment: 24h at Day 16
            RawSampleRecord("t24_d16_1", {"condition": "treatment", "compound": "VPA", "dose": 100, "duration": 24, "stage": "Day 16", "rep": "r1"}, "f", "3"),
            RawSampleRecord("t24_d16_2", {"condition": "treatment", "compound": "VPA", "dose": 100, "duration": 24, "stage": "Day 16", "rep": "r2"}, "f", "4"),
            # Treatment: 48h at Day 16 (different duration!)
            RawSampleRecord("t48_d16_1", {"condition": "treatment", "compound": "VPA", "dose": 100, "duration": 48, "stage": "Day 16", "rep": "r1"}, "f", "5"),
            RawSampleRecord("t48_d16_2", {"condition": "treatment", "compound": "VPA", "dose": 100, "duration": 48, "stage": "Day 16", "rep": "r2"}, "f", "6"),
            # Treatment: 24h at Day 30 (different stage!)
            RawSampleRecord("t24_d30_1", {"condition": "treatment", "compound": "VPA", "dose": 100, "duration": 24, "stage": "Day 30", "rep": "r1"}, "f", "7"),
            RawSampleRecord("t24_d30_2", {"condition": "treatment", "compound": "VPA", "dose": 100, "duration": 24, "stage": "Day 30", "rep": "r2"}, "f", "8"),
        ]
        mapper = CanonicalMetadataMapper()
        doc, _ = mapper.map_study(ctx, recs, source_reference="meta.csv")

        # Must have 4 distinct conditions: 1 control + 3 distinct treatments
        assert len(doc["conditions"]) == 4
        trt_conditions = [c for c in doc["conditions"] if c["treatment_control_status"]["normalized_value"] == "treatment"]
        assert len(trt_conditions) == 3

        # Verify sample counts per condition (each must have exactly 2 samples)
        for cond in doc["conditions"]:
            cid = cond["condition_id"]
            cond_samples = [s for s in doc["samples"] if s["condition_id"] == cid]
            assert len(cond_samples) == 2, f"Condition {cid} has {len(cond_samples)} samples instead of 2"

        # Verify each treatment condition has its correct exposure duration and stage
        cond_ids = {c["condition_id"] for c in trt_conditions}
        assert any("24.0h" in cid and "day_16" in cid for cid in cond_ids)
        assert any("48.0h" in cid and "day_16" in cid for cid in cond_ids)
        assert any("24.0h" in cid and "day_30" in cid for cid in cond_ids)

    def test_mixed_controls_matched_to_corresponding_durations(self):
        """Multi-timepoint study matches treatments to controls with corresponding exposure durations."""
        ctx = StudyContext(study_id="study_timepoint_matched", source_type="public", experiment_id="exp_01")
        recs = [
            # 24h Controls
            RawSampleRecord("c24_1", {"condition": "control", "dose": 0, "duration": 24, "stage": "Day 20", "rep": "r1"}, "f", "1"),
            RawSampleRecord("c24_2", {"condition": "control", "dose": 0, "duration": 24, "stage": "Day 20", "rep": "r2"}, "f", "2"),
            # 48h Controls
            RawSampleRecord("c48_1", {"condition": "control", "dose": 0, "duration": 48, "stage": "Day 20", "rep": "r1"}, "f", "3"),
            RawSampleRecord("c48_2", {"condition": "control", "dose": 0, "duration": 48, "stage": "Day 20", "rep": "r2"}, "f", "4"),
            # 24h Treatments
            RawSampleRecord("t24_1", {"condition": "treatment", "compound": "VPA", "dose": 100, "duration": 24, "stage": "Day 20", "rep": "r1"}, "f", "5"),
            RawSampleRecord("t24_2", {"condition": "treatment", "compound": "VPA", "dose": 100, "duration": 24, "stage": "Day 20", "rep": "r2"}, "f", "6"),
            # 48h Treatments
            RawSampleRecord("t48_1", {"condition": "treatment", "compound": "VPA", "dose": 100, "duration": 48, "stage": "Day 20", "rep": "r1"}, "f", "7"),
            RawSampleRecord("t48_2", {"condition": "treatment", "compound": "VPA", "dose": 100, "duration": 48, "stage": "Day 20", "rep": "r2"}, "f", "8"),
        ]
        mapper = CanonicalMetadataMapper()
        doc, findings = mapper.map_study(ctx, recs, source_reference="meta.csv")

        # 2 control conditions + 2 treatment conditions
        assert len(doc["conditions"]) == 4
        # 2 pairwise relationships formed unambiguously
        assert len(doc["treatment_control_relationships"]) == 2
        assert not any(f.rule == "ambiguous_matched_control" for f in findings)

        # Confirm 24h treatment matches 24h control, and 48h treatment matches 48h control
        for rel in doc["treatment_control_relationships"]:
            t_id = rel["treatment_condition_ids"][0]
            c_id = rel["matched_control_condition_ids"][0]
            if "24.0h" in t_id:
                assert "24.0h" in c_id
            elif "48.0h" in t_id:
                assert "48.0h" in c_id

    def test_unsupported_species_and_source_material_not_reported_as_assumed(self):
        """Species and source material must not be defaulted to Homo sapiens or other with reported status."""
        ctx = StudyContext(study_id="study_missing_bio", source_type="public", experiment_id="exp_01")
        recs = [
            RawSampleRecord("s_ctrl", {"condition": "control", "dose": 0, "rep": "r1"}, "f", "1"),
            RawSampleRecord("s_trt", {"condition": "treatment", "compound": "VPA", "dose": 100, "rep": "r1"}, "f", "2"),
        ]
        mapper = CanonicalMetadataMapper()
        doc, findings = mapper.map_study(ctx, recs, source_reference="meta.csv")

        bio_src = doc["biological_sources"][0]
        # Species: must be unknown, not assumed 'Homo sapiens'
        assert bio_src["species"]["value_status"] == "unknown"
        assert bio_src["species"]["original_value"] is None
        assert bio_src["species"]["normalized_value"] is None

        # Source material: must be unknown, not assumed 'other' or 'iPSC'
        assert bio_src["source_material_type"]["value_status"] == "unknown"
        assert bio_src["source_material_type"]["original_value"] is None
        assert bio_src["source_material_type"]["normalized_value"] is None

        # Must cleanly pass MetadataValidator without controlled_vocabulary error
        v = MetadataValidator()
        val_res = v.validate(doc)
        assert not any(f.rule == "controlled_vocabulary" for f in val_res.findings)

    def test_ambiguous_stage_maps_to_collection_not_exposure_start(self):
        """Ambiguous stage must map to collection age, leaving exposure start age separate and unknown."""
        ctx = StudyContext(study_id="study_ambig_stage", source_type="public", experiment_id="exp_01")
        recs = [
            RawSampleRecord("s_ctrl", {"condition": "control", "dose": 0, "stage": "Day 45", "rep": "r1"}, "f", "1"),
            RawSampleRecord("s_trt", {"condition": "treatment", "compound": "VPA", "dose": 100, "stage": "Day 45", "duration": 48, "rep": "r1"}, "f", "2"),
        ]
        mapper = CanonicalMetadataMapper()
        doc, findings = mapper.map_study(ctx, recs, source_reference="meta.csv")

        # Collection stage on sample: reported as Day 45
        trt_sample = next(s for s in doc["samples"] if s["sample_id"] == "s_trt")
        assert trt_sample["developmental_age_or_stage_at_collection"]["value_status"] == "reported"
        assert trt_sample["developmental_age_or_stage_at_collection"]["original_value"] == "Day 45"

        # Exposure start age on exposure: must remain unknown (not conflated with collection age)
        trt_exp = next(e for e in doc["exposures"] if "vpa" in e["condition_id"])
        assert trt_exp["developmental_age_or_stage_at_exposure"]["value_status"] == "unknown"
        assert trt_exp["developmental_age_or_stage_at_exposure"]["original_value"] is None

    def test_conflicting_duration_prevents_matching_even_with_single_control(self):
        """A 48h treatment must not automatically match a 24h control, even when it is the only control."""
        ctx = StudyContext(study_id="study_conflicting_duration", source_type="public", experiment_id="exp_01")
        recs = [
            RawSampleRecord("c1", {"condition": "control", "dose": 0, "duration": 24, "rep": "r1"}, "f", "1"),
            RawSampleRecord("c2", {"condition": "control", "dose": 0, "duration": 24, "rep": "r2"}, "f", "2"),
            RawSampleRecord("t1", {"condition": "treatment", "compound": "VPA", "dose": 100, "duration": 48, "rep": "r1"}, "f", "3"),
            RawSampleRecord("t2", {"condition": "treatment", "compound": "VPA", "dose": 100, "duration": 48, "rep": "r2"}, "f", "4"),
        ]
        mapper = CanonicalMetadataMapper()
        doc, findings = mapper.map_study(ctx, recs, source_reference="meta.csv")

        # Exactly 2 conditions (1 control 24h, 1 treatment 48h)
        assert len(doc["conditions"]) == 2
        # No relationship formed because durations conflict (48h != 24h)
        assert len(doc["treatment_control_relationships"]) == 0
        assert any(f.rule == "unmatched_treatment_condition" for f in findings)
