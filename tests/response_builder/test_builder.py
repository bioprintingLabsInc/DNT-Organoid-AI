"""Comprehensive scientific and architectural acceptance tests for Treatment-versus-Matched-Control Contrast Builder v1.

Exercises all core acceptance criteria and scientific hardening invariants:
1. One treatment + matched control + same age => one eligible contrast.
2. Day-30 treatment + Day-60 control => blocked; no eligible DE contrast.
3. Same-age Day-30 treatment + Day-30 control and Day-60 treatment + Day-60 control => two separate contrasts.
4. Two treatment doses sharing one explicitly declared control => two contrasts, both using the declared shared control.
5. Missing matched control => blocked.
6. One biological replicate in treatment or control => not inferentially eligible.
7. Multiple true biological replicates => retained as biological replicates, not treated as separate AI observations.
8. Technical replicates must not inflate biological replicate count.
9. Multiple treatment conditions must never be pooled into one treatment arm automatically.
10. Cross-study comparison attempt => blocked.
11. Cross-experiment comparison attempt => blocked.
12. Explicit sample comparison exception is preserved and cannot be silently ignored.
13. DNT labels/evidence do not affect contrast construction.
14. Input ordering changes do not change contrast IDs or output ordering.
15. Diagnostic normalized counts are NOT the DE input; raw count matrix is preserved for downstream modeling.
16. Treatment conditions in normalized cohort are never silently dropped (unlinked condition => auditable BLOCKED disposition).
17. Conservation invariant: total treatment conditions encountered = eligible + needs_review + blocked dispositions.
18. Multiple condition-planned exposure records preserved deterministically without flattening.
19. Exposure singular fields only populated if exactly one exposure (never hide multi-agent).
20. Sample-level exposure deviations distinguished and do not overwrite condition-planned exposures.
21. Developmental age governed normalized values preferred and establish scientific equivalence.
22. Unresolved developmental age text triggers NEEDS_REVIEW (never guessed).
23. Organoid context evaluated on scientific attributes: different IDs with compatible attributes remain eligible.
24. Organoid context mismatch (different types or brain regions) blocks contrast.
25. Unresolved organoid context triggers NEEDS_REVIEW.
26. Unpaired multi-donor biological sources remain eligible.
27. Conflicting species between treatment and control blocks contrast.
"""

import copy
import random
from typing import Any
import pytest

from src.normalization.models import (
    NormalizationCohort,
    NormalizedCohortResult,
    NormalizedCountMatrix,
    NormalizedDataset,
    ResponseEligibility,
    SampleSizeFactor,
    Status as NormalizationStatus,
)
from src.qc.matrix import CountMatrix
from src.response_builder.builder import (
    ResponseContrastBuilder,
    build_response_contrasts,
    extract_metadata_value,
    normalize_age_string,
    resolve_developmental_age,
)
from src.response_builder.models import (
    ContrastStatus,
    Severity,
    Status,
)


def create_normalized_dataset_and_metadata(
    samples_spec: list[dict[str, Any]],
    conditions_spec: list[dict[str, Any]],
    relationships_spec: list[dict[str, Any]],
    exposures_spec: list[dict[str, Any]] | None = None,
    exceptions_spec: list[dict[str, Any]] | None = None,
    organoid_contexts_spec: list[dict[str, Any]] | None = None,
    biological_sources_spec: list[dict[str, Any]] | None = None,
    study_id: str = "study_001",
    experiment_id: str = "exp_001",
    cohort_id: str = "cohort_001",
    cohort_status: NormalizationStatus = NormalizationStatus.PASS,
) -> tuple[NormalizedDataset, dict[str, Any]]:
    """Helper to assemble a scientifically coherent NormalizedDataset and canonical metadata."""
    enriched_samples: list[dict[str, Any]] = []
    for s in samples_spec:
        s_dict = dict(s)
        s_dict.setdefault("organoid_context_id", "org_cerebral")
        s_dict.setdefault("biological_source_id", "source_ipsc_01")
        enriched_samples.append(s_dict)

    sample_ids = tuple(s["sample_id"] for s in enriched_samples)
    gene_ids = ("ENSG00000000001", "ENSG00000000002", "ENSG00000000003")
    raw_columns = tuple(tuple(100 + i * 10 for i in range(len(gene_ids))) for _ in sample_ids)
    raw_matrix = CountMatrix(gene_ids=gene_ids, sample_ids=sample_ids, columns=raw_columns)

    norm_columns = tuple(tuple(100.0 for _ in range(len(gene_ids))) for _ in sample_ids)
    diagnostic_matrix = NormalizedCountMatrix(gene_ids=gene_ids, sample_ids=sample_ids, columns=norm_columns)

    # Build cohort
    t_sids = tuple(s["sample_id"] for s in enriched_samples if s.get("condition_id", "").startswith("treat"))
    c_sids = tuple(s["sample_id"] for s in enriched_samples if s.get("condition_id", "").startswith("ctrl") or "control" in s.get("condition_id", ""))

    cohort = NormalizationCohort(
        cohort_id=cohort_id,
        study_id=study_id,
        experiment_id=experiment_id,
        sample_ids=sample_ids,
        assay_ids=("assay_001",),
        rna_seq_modality="bulk_rna_seq",
        library_strategy="polyA",
        sequencing_platform="Illumina NovaSeq",
        strandedness="reverse",
        read_layout="paired_end",
        treatment_sample_ids=t_sids,
        control_sample_ids=c_sids,
        unassigned_sample_ids=(),
        treatment_control_relationship_ids=tuple(r["relationship_id"] for r in relationships_spec),
        metadata_hash="dummy_hash",
    )

    size_factors = tuple(
        SampleSizeFactor(
            sample_id=sid,
            size_factor=1.0,
            raw_library_size=300,
            normalized_library_size=300.0,
            cohort_id=cohort_id,
            response_eligibility=ResponseEligibility.ELIGIBLE_WITH_MATCHED_CONTROL if sid in t_sids else ResponseEligibility.CONTROL_BASELINE,
            matched_control_sample_ids=c_sids if sid in t_sids else (),
            treatment_control_status="treatment" if sid in t_sids else "control",
        )
        for sid in sample_ids
    )

    cohort_result = NormalizedCohortResult(
        status=cohort_status,
        cohort=cohort,
        findings=(),
        size_factors=size_factors,
        metrics=None,
        diagnostic_normalized_matrix=diagnostic_matrix,
        raw_count_matrix=raw_matrix,
        method="ratio",
        r_environment_info={"r_version": "4.4.3", "bioc_version": "3.20", "deseq2_version": "1.46.0"},
        execution_provenance={},
    )

    normalized_dataset = NormalizedDataset(
        status=cohort_status,
        findings=(),
        cohort_results=(cohort_result,),
        source_asset_id="asset_001",
        reference_identity="Ensembl Release 112",
        config_version="1.0.0",
        normalization_version="bulk_normalization_v1",
    )

    metadata: dict[str, Any] = {
        "studies": [{"study_id": study_id}],
        "experiments": [{"experiment_id": experiment_id, "study_id": study_id}],
        "conditions": conditions_spec,
        "samples": enriched_samples,
        "treatment_control_relationships": relationships_spec,
        "exposures": exposures_spec or [],
        "sample_comparison_exceptions": exceptions_spec or [],
        "organoid_contexts": organoid_contexts_spec or [{
            "organoid_context_id": "org_cerebral",
            "organoid_type": "cerebral_organoid",
            "brain_region_or_model_identity": "cerebral_cortex",
            "differentiation_protocol_identifier": "pasca_2015",
            "differentiation_protocol_version": "v1",
        }],
        "biological_sources": biological_sources_spec or [{
            "biological_source_id": "source_ipsc_01",
            "organism_or_species": "Homo sapiens",
        }],
        "sequencing_assays": [{
            "assay_id": "assay_001",
            "sample_ids": list(sample_ids),
            "rna_seq_modality": {"original_value": "bulk_rna_seq", "value_status": "reported"},
        }],
    }

    return normalized_dataset, metadata


# =========================================================================
# Acceptance Test 1: One treatment + matched control + same age => 1 eligible contrast
# =========================================================================
def test_acceptance_1_one_treatment_matched_control_same_age() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {
            "relationship_id": "rel_001",
            "treatment_condition_ids": ["cond_treat"],
            "matched_control_condition_ids": ["cond_ctrl"],
            "relationship_type": "vehicle_control",
        }
    ]
    exposures = [
        {
            "exposure_id": "exp_vpa",
            "condition_id": "cond_treat",
            "agent_name": "Valproic Acid",
            "agent_identifier": "CID_107",
            "concentration_or_dose": 250,
            "concentration_or_dose_unit": "uM",
            "exposure_start_time_or_stage": "Day 25",
            "developmental_age_or_stage_at_exposure": "Day 25",
            "exposure_duration": 5,
            "exposure_duration_unit": "days",
        }
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships, exposures_spec=exposures)
    res = build_response_contrasts(norm, meta)

    assert len(res.contrasts) == 1
    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.ELIGIBLE
    assert contrast.is_inferentially_eligible is True
    assert contrast.treatment_condition_id == "cond_treat"
    assert contrast.matched_control_condition_ids == ("cond_ctrl",)
    assert contrast.treatment_collection_age_or_stage == "Day 30"
    assert contrast.control_collection_age_or_stage == "Day 30"
    assert contrast.treatment_collection_age_normalized == "day_30"
    assert contrast.control_collection_age_normalized == "day_30"
    assert contrast.agent_name == "Valproic Acid"
    assert len(contrast.treatment_exposures) == 1


# =========================================================================
# Acceptance Test 2: Day-30 treatment + Day-60 control => blocked
# =========================================================================
def test_acceptance_2_developmental_age_mismatch_blocks_contrast() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 60"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 60"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {
            "relationship_id": "rel_mismatch_age",
            "treatment_condition_ids": ["cond_treat"],
            "matched_control_condition_ids": ["cond_ctrl"],
        }
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    assert len(res.contrasts) == 1
    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.BLOCKED
    assert contrast.is_inferentially_eligible is False
    assert any(f.rule_id == "developmental_age_mismatch" for f in contrast.findings)


# =========================================================================
# Acceptance Test 3: Same-age Day 30 and Day 60 => 2 separate contrasts
# =========================================================================
def test_acceptance_3_distinct_developmental_ages_produce_two_separate_contrasts() -> None:
    samples = [
        # Day 30 pair
        {"sample_id": "s_t30_1", "experiment_id": "exp_001", "condition_id": "treat_d30", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t30_2", "experiment_id": "exp_001", "condition_id": "treat_d30", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c30_1", "experiment_id": "exp_001", "condition_id": "ctrl_d30", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c30_2", "experiment_id": "exp_001", "condition_id": "ctrl_d30", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
        # Day 60 pair
        {"sample_id": "s_t60_1", "experiment_id": "exp_001", "condition_id": "treat_d60", "biological_replicate_id": "r3", "developmental_age_or_stage_at_collection": "Day 60"},
        {"sample_id": "s_t60_2", "experiment_id": "exp_001", "condition_id": "treat_d60", "biological_replicate_id": "r4", "developmental_age_or_stage_at_collection": "Day 60"},
        {"sample_id": "s_c60_1", "experiment_id": "exp_001", "condition_id": "ctrl_d60", "biological_replicate_id": "cr3", "developmental_age_or_stage_at_collection": "Day 60"},
        {"sample_id": "s_c60_2", "experiment_id": "exp_001", "condition_id": "ctrl_d60", "biological_replicate_id": "cr4", "developmental_age_or_stage_at_collection": "Day 60"},
    ]
    conditions = [
        {"condition_id": "treat_d30", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "ctrl_d30", "experiment_id": "exp_001", "treatment_control_status": "control"},
        {"condition_id": "treat_d60", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "ctrl_d60", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_d30", "treatment_condition_ids": ["treat_d30"], "matched_control_condition_ids": ["ctrl_d30"]},
        {"relationship_id": "rel_d60", "treatment_condition_ids": ["treat_d60"], "matched_control_condition_ids": ["ctrl_d60"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    assert len(res.contrasts) == 2
    c30 = next(c for c in res.contrasts if c.treatment_condition_id == "treat_d30")
    c60 = next(c for c in res.contrasts if c.treatment_condition_id == "treat_d60")

    assert c30.treatment_collection_age_or_stage == "Day 30"
    assert c30.control_collection_age_or_stage == "Day 30"
    assert c30.status == ContrastStatus.ELIGIBLE
    assert c30.is_inferentially_eligible is True

    assert c60.treatment_collection_age_or_stage == "Day 60"
    assert c60.control_collection_age_or_stage == "Day 60"
    assert c60.status == ContrastStatus.ELIGIBLE
    assert c60.is_inferentially_eligible is True


# =========================================================================
# Acceptance Test 4: Two doses sharing one control => 2 separate contrasts
# =========================================================================
def test_acceptance_4_multiple_doses_sharing_one_explicit_control() -> None:
    samples = [
        {"sample_id": "s_t10_1", "experiment_id": "exp_001", "condition_id": "treat_10uM", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t10_2", "experiment_id": "exp_001", "condition_id": "treat_10uM", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t50_1", "experiment_id": "exp_001", "condition_id": "treat_50uM", "biological_replicate_id": "r3", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t50_2", "experiment_id": "exp_001", "condition_id": "treat_50uM", "biological_replicate_id": "r4", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_ctrl_1", "experiment_id": "exp_001", "condition_id": "ctrl_veh", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_ctrl_2", "experiment_id": "exp_001", "condition_id": "ctrl_veh", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "treat_10uM", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "treat_50uM", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "ctrl_veh", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {
            "relationship_id": "rel_dose_response",
            "treatment_condition_ids": ["treat_10uM", "treat_50uM"],
            "matched_control_condition_ids": ["ctrl_veh"],
        }
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    assert len(res.contrasts) == 2
    c10 = next(c for c in res.contrasts if c.treatment_condition_id == "treat_10uM")
    c50 = next(c for c in res.contrasts if c.treatment_condition_id == "treat_50uM")

    assert c10.matched_control_condition_ids == ("ctrl_veh",)
    assert c10.control_sample_ids == ("s_ctrl_1", "s_ctrl_2")
    assert c10.status == ContrastStatus.ELIGIBLE

    assert c50.matched_control_condition_ids == ("ctrl_veh",)
    assert c50.control_sample_ids == ("s_ctrl_1", "s_ctrl_2")
    assert c50.status == ContrastStatus.ELIGIBLE


# =========================================================================
# Acceptance Test 5: Missing matched control => blocked
# =========================================================================
def test_acceptance_5_missing_matched_control_blocks_contrast() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
    ]
    relationships = [
        {
            "relationship_id": "rel_no_ctrl",
            "treatment_condition_ids": ["cond_treat"],
            "matched_control_condition_ids": [],
        }
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    assert len(res.contrasts) == 1
    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "missing_matched_control" for f in contrast.findings)


# =========================================================================
# Acceptance Test 6: Single biological replicate => not inferentially eligible
# =========================================================================
def test_acceptance_6_single_biological_replicate_is_not_inferentially_eligible() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "rep_single", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "c_rep1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "c_rep2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_01", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    assert len(res.contrasts) == 1
    contrast = res.contrasts[0]
    assert contrast.treatment_replicate_count == 1
    assert contrast.control_replicate_count == 2
    assert contrast.is_inferentially_eligible is False
    assert contrast.status == ContrastStatus.NEEDS_REVIEW
    assert any(f.rule_id == "insufficient_biological_replication" for f in contrast.findings)


# =========================================================================
# Acceptance Test 7: Multiple true biological replicates retained properly
# =========================================================================
def test_acceptance_7_multiple_true_biological_replicates_retained_properly() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "rep_alpha", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "rep_beta", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t3", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "rep_gamma", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "c_alpha", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "c_beta", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c3", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "c_gamma", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_01", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    assert len(res.contrasts) == 1
    contrast = res.contrasts[0]
    assert contrast.treatment_replicate_count == 3
    assert contrast.control_replicate_count == 3
    assert contrast.treatment_biological_replicate_ids == ("rep_alpha", "rep_beta", "rep_gamma")
    assert contrast.control_biological_replicate_ids == ("c_alpha", "c_beta", "c_gamma")
    assert contrast.is_inferentially_eligible is True


# =========================================================================
# Acceptance Test 8: Technical replicates do not inflate biological count
# =========================================================================
def test_acceptance_8_technical_replicates_do_not_inflate_biological_count() -> None:
    # 4 sequencing libraries from only 2 true biological replicates
    samples = [
        {"sample_id": "s_t1_tech1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "bio_rep_1", "technical_replicate_id": "tech_1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t1_tech2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "bio_rep_1", "technical_replicate_id": "tech_2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2_tech1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "bio_rep_2", "technical_replicate_id": "tech_1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2_tech2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "bio_rep_2", "technical_replicate_id": "tech_2", "developmental_age_or_stage_at_collection": "Day 30"},
        # Controls: 2 biological replicates
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "ctrl_bio_1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "ctrl_bio_2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_01", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    contrast = res.contrasts[0]
    # Treatment has 4 sample IDs, but only 2 biological replicate IDs
    assert len(contrast.treatment_sample_ids) == 4
    assert contrast.treatment_replicate_count == 2
    assert contrast.control_replicate_count == 2
    assert contrast.is_inferentially_eligible is True


# =========================================================================
# Acceptance Test 9: Multiple treatment conditions never pooled automatically
# =========================================================================
def test_acceptance_9_multiple_treatment_conditions_never_pooled() -> None:
    samples = [
        {"sample_id": "s_tA_1", "experiment_id": "exp_001", "condition_id": "treat_A", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_tA_2", "experiment_id": "exp_001", "condition_id": "treat_A", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_tB_1", "experiment_id": "exp_001", "condition_id": "treat_B", "biological_replicate_id": "r3", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_tB_2", "experiment_id": "exp_001", "condition_id": "treat_B", "biological_replicate_id": "r4", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "ctrl_veh", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "ctrl_veh", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "treat_A", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "treat_B", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "ctrl_veh", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {
            "relationship_id": "rel_multi",
            "treatment_condition_ids": ["treat_A", "treat_B"],
            "matched_control_condition_ids": ["ctrl_veh"],
        }
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    assert len(res.contrasts) == 2
    treat_cids = {c.treatment_condition_id for c in res.contrasts}
    assert treat_cids == {"treat_A", "treat_B"}


# =========================================================================
# Acceptance Test 10: Cross-study comparison blocked
# =========================================================================
def test_acceptance_10_cross_study_comparison_blocked() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_002", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_002", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_002", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_cross_study", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    # Inject second study into metadata
    meta["studies"].append({"study_id": "study_002"})
    meta["experiments"].append({"experiment_id": "exp_002", "study_id": "study_002"})

    res = build_response_contrasts(norm, meta)
    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "cross_study_comparison_blocked" for f in contrast.findings)


# =========================================================================
# Acceptance Test 11: Cross-experiment comparison blocked
# =========================================================================
def test_acceptance_11_cross_experiment_comparison_blocked() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_002", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_002", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_002", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_cross_exp", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    meta["experiments"].append({"experiment_id": "exp_002", "study_id": "study_001"})

    res = build_response_contrasts(norm, meta)
    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "cross_experiment_comparison_blocked" for f in contrast.findings)


# =========================================================================
# Acceptance Test 12: Sample comparison exception handling with and without replacement
# =========================================================================
def test_acceptance_12_sample_comparison_exception_with_and_without_replacement() -> None:
    # 1. With valid replacement control
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1_bad", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2_good", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c3_rep", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr3", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_exc", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]
    exceptions = [
        {
            "exception_id": "exc_01",
            "relationship_id": "rel_exc",
            "affected_sample_ids": ["s_c1_bad"],
            "replacement_control_sample_ids": ["s_c3_rep"],
            "exception_reason": {"original_value": "s_c1_bad has high dead cell count", "value_status": "reported"},
        }
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships, exceptions_spec=exceptions)
    res = build_response_contrasts(norm, meta)

    contrast = res.contrasts[0]
    assert "exc_01" in contrast.sample_comparison_exception_ids
    assert "s_c1_bad" not in contrast.control_sample_ids
    assert "s_c3_rep" in contrast.control_sample_ids
    assert contrast.status == ContrastStatus.ELIGIBLE
    assert any(f.rule_id == "sample_comparison_exception_applied" for f in contrast.findings)

    # 2. Without replacement control => NEEDS_REVIEW
    exceptions_no_rep = [
        {
            "exception_id": "exc_02",
            "relationship_id": "rel_exc",
            "affected_sample_ids": ["s_c1_bad"],
            "replacement_control_sample_ids": None,
            "exception_reason": {"original_value": "s_c1_bad disqualified", "value_status": "reported"},
        }
    ]
    norm2, meta2 = create_normalized_dataset_and_metadata(samples, conditions, relationships, exceptions_spec=exceptions_no_rep)
    res2 = build_response_contrasts(norm2, meta2)
    c2 = res2.contrasts[0]
    assert c2.status == ContrastStatus.NEEDS_REVIEW
    assert any(f.rule_id == "unresolved_sample_comparison_exception" for f in c2.findings)


# =========================================================================
# Acceptance Test 13: DNT labels have zero effect on contrast construction
# =========================================================================
def test_acceptance_13_dnt_labels_have_zero_effect_on_contrast_construction() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_01", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm1, meta1 = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res1 = build_response_contrasts(norm1, meta1)

    # Now add toxicant annotations and DNT labels
    meta2 = copy.deepcopy(meta1)
    meta2["dnt_labels"] = [{"condition_id": "cond_treat", "dnt_positive": True, "evidence_strength": "definitive"}]
    meta2["toxicant_annotations"] = [{"agent_name": "Valproic Acid", "cas_number": "99-66-1", "is_known_dnt": True}]

    res2 = build_response_contrasts(norm1, meta2)

    assert len(res1.contrasts) == len(res2.contrasts)
    c1 = res1.contrasts[0]
    c2 = res2.contrasts[0]
    assert c1.contrast_id == c2.contrast_id
    assert c1.status == c2.status
    assert c1.is_inferentially_eligible == c2.is_inferentially_eligible
    assert len(c1.findings) == len(c2.findings)


# =========================================================================
# Acceptance Test 14: Input ordering invariance
# =========================================================================
def test_acceptance_14_input_ordering_invariance() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_01", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm1, meta1 = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res1 = build_response_contrasts(norm1, meta1)

    # Shuffled input order
    shuffled_samples = list(samples)
    shuffled_conditions = list(conditions)
    random.Random(42).shuffle(shuffled_samples)
    random.Random(42).shuffle(shuffled_conditions)

    norm2, meta2 = create_normalized_dataset_and_metadata(shuffled_samples, shuffled_conditions, relationships)
    res2 = build_response_contrasts(norm2, meta2)

    assert len(res1.contrasts) == len(res2.contrasts)
    for c1, c2 in zip(res1.contrasts, res2.contrasts, strict=True):
        assert c1.contrast_id == c2.contrast_id
        assert c1.treatment_condition_id == c2.treatment_condition_id
        assert c1.treatment_sample_ids == c2.treatment_sample_ids
        assert c1.control_sample_ids == c2.control_sample_ids


# =========================================================================
# Acceptance Test 15: Raw count matrix preserved on normalization output
# =========================================================================
def test_acceptance_15_raw_count_matrix_is_preserved_for_downstream_de() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_01", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    contrast = res.contrasts[0]
    cohort_result = norm.cohort_results[0]
    assert contrast.cohort_id == cohort_result.cohort.cohort_id
    assert cohort_result.raw_count_matrix is not None
    for col in cohort_result.raw_count_matrix.columns:
        for val in col:
            assert isinstance(val, int)
            assert val >= 0
    assert cohort_result.diagnostic_normalized_matrix is not None
    for col in cohort_result.diagnostic_normalized_matrix.columns:
        for val in col:
            assert isinstance(val, float)


# =========================================================================
# Hardening Test 16: No silently dropped treatment conditions
# =========================================================================
def test_no_silently_dropped_treatment_conditions_unlinked_is_blocked() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat_linked", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat_linked", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
        # Unlinked treatment condition with samples in cohort but NO relationship declared
        {"sample_id": "s_unlinked_1", "experiment_id": "exp_001", "condition_id": "cond_treat_unlinked", "biological_replicate_id": "u1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_unlinked_2", "experiment_id": "exp_001", "condition_id": "cond_treat_unlinked", "biological_replicate_id": "u2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat_linked", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
        {"condition_id": "cond_treat_unlinked", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
    ]
    relationships = [
        {"relationship_id": "rel_01", "treatment_condition_ids": ["cond_treat_linked"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    # Must contain both conditions
    cids = {c.treatment_condition_id for c in res.contrasts}
    assert "cond_treat_linked" in cids
    assert "cond_treat_unlinked" in cids

    unlinked_contrast = next(c for c in res.contrasts if c.treatment_condition_id == "cond_treat_unlinked")
    assert unlinked_contrast.status == ContrastStatus.BLOCKED
    assert unlinked_contrast.is_inferentially_eligible is False
    assert unlinked_contrast.study_id == "study_001"
    assert unlinked_contrast.experiment_id == "exp_001"
    assert unlinked_contrast.treatment_sample_ids == ("s_unlinked_1", "s_unlinked_2")
    assert any(f.rule_id == "missing_explicit_treatment_control_relationship" for f in unlinked_contrast.findings)


# =========================================================================
# Hardening Test 17: Conservation Invariant for treatment conditions
# =========================================================================
def test_treatment_condition_conservation_invariant() -> None:
    samples = [
        # treat_eligible (has matched control, N=2)
        {"sample_id": "s_e1", "experiment_id": "exp_001", "condition_id": "treat_eligible", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_e2", "experiment_id": "exp_001", "condition_id": "treat_eligible", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
        # treat_under_rep (has matched control, N=1 => NEEDS_REVIEW)
        {"sample_id": "s_u1", "experiment_id": "exp_001", "condition_id": "treat_under_rep", "biological_replicate_id": "single_rep", "developmental_age_or_stage_at_collection": "Day 30"},
        # treat_blocked (missing matched control => BLOCKED)
        {"sample_id": "s_b1", "experiment_id": "exp_001", "condition_id": "treat_blocked", "biological_replicate_id": "br1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_b2", "experiment_id": "exp_001", "condition_id": "treat_blocked", "biological_replicate_id": "br2", "developmental_age_or_stage_at_collection": "Day 30"},
        # treat_unlinked (no relationship => BLOCKED)
        {"sample_id": "s_no_rel_1", "experiment_id": "exp_001", "condition_id": "treat_unlinked", "biological_replicate_id": "nr1", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "treat_eligible", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
        {"condition_id": "treat_under_rep", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "treat_blocked", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "treat_unlinked", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
    ]
    relationships = [
        {"relationship_id": "rel_el", "treatment_condition_ids": ["treat_eligible"], "matched_control_condition_ids": ["cond_ctrl"]},
        {"relationship_id": "rel_ur", "treatment_condition_ids": ["treat_under_rep"], "matched_control_condition_ids": ["cond_ctrl"]},
        {"relationship_id": "rel_bl", "treatment_condition_ids": ["treat_blocked"], "matched_control_condition_ids": []},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    assert len(res.treatment_conditions_encountered) == 4
    assert res.eligible_dispositions_count == 1
    assert res.needs_review_dispositions_count == 1
    assert res.blocked_dispositions_count == 2
    assert res.verify_treatment_condition_conservation() is True
    assert len(res.treatment_conditions_encountered) == (
        res.eligible_dispositions_count + res.needs_review_dispositions_count + res.blocked_dispositions_count
    )


# =========================================================================
# Hardening Test 18-20: Multiple Exposure Records Preservation
# =========================================================================
def test_multiple_exposure_records_preserved_deterministically() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_combo", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_combo", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_combo", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_combo", "treatment_condition_ids": ["cond_combo"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]
    # Two planned exposures: Drug + Cytokine
    exposures = [
        {
            "exposure_id": "exp_drug_vpa",
            "condition_id": "cond_combo",
            "agent_name": "Valproic Acid",
            "agent_identifier": "CID_107",
            "concentration_or_dose": 300,
            "concentration_or_dose_unit": "uM",
            "exposure_start_time_or_stage": "Day 25",
            "developmental_age_or_stage_at_exposure": "Day 25",
            "exposure_duration": 5,
            "exposure_duration_unit": "days",
        },
        {
            "exposure_id": "exp_cytokine_il6",
            "condition_id": "cond_combo",
            "agent_name": "Interleukin-6",
            "agent_identifier": "IL6_HUMAN",
            "concentration_or_dose": 50,
            "concentration_or_dose_unit": "ng/mL",
            "exposure_start_time_or_stage": "Day 25",
            "developmental_age_or_stage_at_exposure": "Day 25",
            "exposure_duration": 5,
            "exposure_duration_unit": "days",
        },
        # Sample-specific exposure deviation
        {
            "exposure_id": "exp_dev_t2",
            "sample_id": "s_t2",
            "agent_name": "Valproic Acid",
            "concentration_or_dose": 280,
            "concentration_or_dose_unit": "uM",
        },
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships, exposures_spec=exposures)
    res = build_response_contrasts(norm, meta)

    contrast = res.contrasts[0]
    # Invariant: both planned exposures retained
    assert len(contrast.treatment_exposures) == 2
    exp_ids = tuple(e.exposure_id for e in contrast.treatment_exposures)
    assert exp_ids == ("exp_cytokine_il6", "exp_drug_vpa")  # sorted deterministically
    assert contrast.treatment_exposures[0].agent_name == "Interleukin-6"
    assert contrast.treatment_exposures[1].agent_name == "Valproic Acid"

    # Multi-agent invariant: singular fields must NEVER hide multi-agent condition
    assert contrast.agent_name is None
    assert contrast.agent_identifier is None
    assert contrast.concentration_or_dose is None

    # Sample-specific deviation invariant: preserved in deviations, does not overwrite planned
    assert len(contrast.sample_level_exposure_deviations) == 1
    assert contrast.sample_level_exposure_deviations[0]["sample_id"] == "s_t2"
    assert any(f.rule_id == "sample_level_exposure_deviation_present" for f in contrast.findings)

    # Determinism test under reordered metadata exposures
    shuffled_exposures = list(reversed(exposures))
    norm2, meta2 = create_normalized_dataset_and_metadata(samples, conditions, relationships, exposures_spec=shuffled_exposures)
    res2 = build_response_contrasts(norm2, meta2)
    contrast2 = res2.contrasts[0]
    assert contrast2.contrast_id == contrast.contrast_id
    assert tuple(e.exposure_id for e in contrast2.treatment_exposures) == exp_ids


# =========================================================================
# Hardening Test 21: Governed normalized age establishes equivalence
# =========================================================================
def test_developmental_age_governed_normalized_equivalence_accepted() -> None:
    samples = [
        # Treatment has "30 days", normalized_value "day_30"
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": {"original_value": "30 days", "normalized_value": "day_30", "value_status": "reported"}},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": {"original_value": "day 30", "normalized_value": "day_30", "value_status": "reported"}},
        # Control has "Day 30", normalized_value "day_30"
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": {"original_value": "Day 30", "normalized_value": "day_30", "value_status": "reported"}},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": {"original_value": "Day 30", "normalized_value": "day_30", "value_status": "reported"}},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_age_gov", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.ELIGIBLE
    assert contrast.treatment_collection_age_normalized == "day_30"
    assert contrast.control_collection_age_normalized == "day_30"


# =========================================================================
# Hardening Test 22: Unresolved textual age triggers NEEDS_REVIEW (never guessed)
# =========================================================================
def test_unresolved_developmental_age_triggers_needs_review_not_guessed() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "mid-maturation"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "mid-maturation"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "stage 2"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "stage 2"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_unresolved", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.NEEDS_REVIEW
    assert any(f.rule_id == "unresolved_developmental_age_equivalence" for f in contrast.findings)


# =========================================================================
# Hardening Test 23: Organoid Context - Different IDs with equivalent attributes remain eligible
# =========================================================================
def test_organoid_context_different_ids_equivalent_attributes_eligible() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "ctx_treat_1"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "ctx_treat_1"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "ctx_ctrl_1"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "ctx_ctrl_1"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_compat", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]
    contexts = [
        {
            "organoid_context_id": "ctx_treat_1",
            "organoid_type": "cerebral_organoid",
            "brain_region_or_model_identity": "cerebral_cortex",
            "differentiation_protocol_identifier": "pasca_2015",
            "differentiation_protocol_version": "v2",
        },
        {
            "organoid_context_id": "ctx_ctrl_1",
            "organoid_type": "cerebral_organoid",
            "brain_region_or_model_identity": "cerebral_cortex",
            "differentiation_protocol_identifier": "pasca_2015",
            "differentiation_protocol_version": "v2",
        },
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships, organoid_contexts_spec=contexts)
    res = build_response_contrasts(norm, meta)

    contrast = res.contrasts[0]
    # Scientific invariant: Different IDs do NOT block if governed scientific attributes are compatible
    assert contrast.status == ContrastStatus.ELIGIBLE
    assert contrast.organoid_context_ids == ("ctx_ctrl_1", "ctx_treat_1")


# =========================================================================
# Hardening Test 24: Organoid context mismatch (different type/region) blocks contrast
# =========================================================================
def test_organoid_context_mismatch_blocks_contrast() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_cerebral"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_cerebral"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_hepatic"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_hepatic"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_org_mismatch", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]
    contexts = [
        {"organoid_context_id": "org_cerebral", "organoid_type": "cerebral_organoid", "brain_region_or_model_identity": "cerebral_cortex"},
        {"organoid_context_id": "org_hepatic", "organoid_type": "hepatic_organoid", "brain_region_or_model_identity": "liver"},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships, organoid_contexts_spec=contexts)
    res = build_response_contrasts(norm, meta)

    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "organoid_context_mismatch" for f in contrast.findings)


# =========================================================================
# Hardening Test 25: Unresolved organoid context triggers NEEDS_REVIEW
# =========================================================================
def test_organoid_context_unresolved_metadata_needs_review() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_ctx_missing"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_ctx_missing"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_cerebral"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_cerebral"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_missing_ctx", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.NEEDS_REVIEW
    assert any(f.rule_id == "unresolved_organoid_context_compatibility" for f in contrast.findings)


# =========================================================================
# Hardening Test 26: Biological Source - Unpaired multi-donor design remains eligible
# =========================================================================
def test_biological_source_unpaired_multi_donor_eligible() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30", "biological_source_id": "donor_A"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30", "biological_source_id": "donor_A"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30", "biological_source_id": "donor_B"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30", "biological_source_id": "donor_B"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_multi_donor", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]
    sources = [
        {"biological_source_id": "donor_A", "organism_or_species": "Homo sapiens"},
        {"biological_source_id": "donor_B", "organism_or_species": "Homo sapiens"},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships, biological_sources_spec=sources)
    res = build_response_contrasts(norm, meta)

    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.ELIGIBLE
    assert contrast.biological_source_ids == ("donor_A", "donor_B")


# =========================================================================
# Hardening Test 27: Cross-species comparison blocked
# =========================================================================
def test_cross_species_comparison_blocked() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30", "biological_source_id": "human_cell"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30", "biological_source_id": "human_cell"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30", "biological_source_id": "mouse_cell"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30", "biological_source_id": "mouse_cell"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_species_mismatch", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]
    sources = [
        {"biological_source_id": "human_cell", "organism_or_species": "Homo sapiens"},
        {"biological_source_id": "mouse_cell", "organism_or_species": "Mus musculus"},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships, biological_sources_spec=sources)
    res = build_response_contrasts(norm, meta)

    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "cross_species_comparison_blocked" for f in contrast.findings)


# =========================================================================
# Additional Scientific Invariant Tests
# =========================================================================
def test_ambiguous_developmental_age_needs_review() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 35"},  # Divergent age within condition
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_age_ambig", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.NEEDS_REVIEW
    assert any(f.rule_id == "developmental_age_ambiguous" for f in contrast.findings)


def test_treatment_control_role_inconsistency_blocks_contrast() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    # Inverted: cond_treat has status 'control', cond_ctrl has status 'treatment'
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "control"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
    ]
    relationships = [
        {"relationship_id": "rel_inverted", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    contrast = next(c for c in res.contrasts if c.treatment_condition_id == "cond_treat")
    assert contrast.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "treatment_control_role_inconsistency" for f in contrast.findings)


def test_normalize_age_string_helper() -> None:
    assert normalize_age_string("Day 30") == "day_30"
    assert normalize_age_string("day-30") == "day_30"
    assert normalize_age_string("d30") == "day_30"
    assert normalize_age_string("30 days") == "day_30"
    assert normalize_age_string(30) == "day_30"
    assert normalize_age_string("GW 10") == "gw_10"
    assert normalize_age_string("gestational week 10") == "gw_10"
    assert normalize_age_string("10 gw") == "gw_10"
    assert normalize_age_string("pcw 8") == "pcw_8"
    assert normalize_age_string("week 4") == "week_4"
    assert normalize_age_string("4 weeks") == "week_4"
    # Never invent a normalized age for unstandardized arbitrary text
    assert normalize_age_string("differentiated_stage_A") is None
    assert normalize_age_string(None) is None
    assert normalize_age_string("unknown") is None
    assert normalize_age_string("missing_not_reported") is None
