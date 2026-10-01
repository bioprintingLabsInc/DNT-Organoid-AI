"""Comprehensive scientific and architectural acceptance tests for Treatment-versus-Matched-Control Contrast Builder v1.

Exercises all 15 required acceptance criteria:
1. One treatment + matched control + same age => one eligible contrast.
2. Day-30 treatment + Day-60 control => blocked; no eligible DE contrast.
3. Same-age Day-30 treatment + Day-30 control and Day-60 treatment + Day-60 control => two separate contrasts.
4. Two treatment doses sharing one explicitly declared control => two contrasts, both using the declared shared control.
5. Missing matched control => blocked.
6. One biological replicate in treatment or control => not inferentially eligible.
7. Multiple true biological replicates => retained as biological replicates, not treated as separate AI observations.
8. Technical replicates must not inflate biological replicate count.
9. Multiple treatment conditions must never be pooled into one treatment arm automatically.
10. Cross-study pooling attempt => blocked.
11. Cross-experiment comparison attempt => blocked.
12. Explicit sample comparison exception is preserved and cannot be silently ignored.
13. DNT labels/evidence do not affect contrast construction.
14. Input ordering changes do not change contrast IDs or output ordering.
15. Diagnostic normalized counts are NOT the DE input; raw count matrix is preserved for downstream modeling.
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
    normalize_age_string,
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
    sample_ids = tuple(s["sample_id"] for s in samples_spec)
    gene_ids = ("ENSG00000000001", "ENSG00000000002", "ENSG00000000003")
    raw_columns = tuple(tuple(100 + i * 10 for i in range(len(gene_ids))) for _ in sample_ids)
    raw_matrix = CountMatrix(gene_ids=gene_ids, sample_ids=sample_ids, columns=raw_columns)

    norm_columns = tuple(tuple(100.0 for _ in range(len(gene_ids))) for _ in sample_ids)
    diagnostic_matrix = NormalizedCountMatrix(gene_ids=gene_ids, sample_ids=sample_ids, columns=norm_columns)

    # Build cohort
    t_sids = tuple(s["sample_id"] for s in samples_spec if s.get("condition_id", "").startswith("treat"))
    c_sids = tuple(s["sample_id"] for s in samples_spec if s.get("condition_id", "").startswith("ctrl") or "control" in s.get("condition_id", ""))

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
        "samples": samples_spec,
        "treatment_control_relationships": relationships_spec,
        "exposures": exposures_spec or [],
        "sample_comparison_exceptions": exceptions_spec or [],
        "organoid_contexts": organoid_contexts_spec or [{"organoid_context_id": "org_cerebral", "brain_region_or_model_identity": "cerebral_cortex"}],
        "biological_sources": biological_sources_spec or [{"biological_source_id": "source_ipsc_01", "species": "Homo sapiens"}],
        "sequencing_assays": [{
            "assay_id": "assay_001",
            "sample_ids": list(sample_ids),
            "rna_seq_modality": {"original_value": "bulk_rna_seq", "value_status": "reported"},
        }],
    }

    return normalized_dataset, metadata


# =========================================================================
# Acceptance Test 1: One treatment + matched control + same age => one eligible contrast
# =========================================================================
def test_acceptance_1_one_treatment_matched_control_same_age() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "rep1", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_cerebral", "biological_source_id": "source_01"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "rep2", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_cerebral", "biological_source_id": "source_01"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "c_rep1", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_cerebral", "biological_source_id": "source_01"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "c_rep2", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_cerebral", "biological_source_id": "source_01"},
    ]
    conditions = [
        {"condition_id": "cond_treat", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_01", "treatment_condition_ids": ["cond_treat"], "matched_control_condition_ids": ["cond_ctrl"]},
    ]
    exposures = [
        {"exposure_id": "exp_01", "condition_id": "cond_treat", "agent_name": "Valproic Acid", "concentration_or_dose": 200, "concentration_or_dose_unit": "uM", "developmental_age_or_stage_at_exposure": "Day 25", "exposure_duration": 5, "exposure_duration_unit": "days"},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships, exposures)
    res = build_response_contrasts(norm, meta)

    assert len(res.contrasts) == 1
    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.ELIGIBLE
    assert contrast.is_inferentially_eligible is True
    assert contrast.treatment_condition_id == "cond_treat"
    assert contrast.matched_control_condition_ids == ("cond_ctrl",)
    assert contrast.treatment_sample_ids == ("s_t1", "s_t2")
    assert contrast.control_sample_ids == ("s_c1", "s_c2")
    assert contrast.treatment_replicate_count == 2
    assert contrast.control_replicate_count == 2
    assert contrast.agent_name == "Valproic Acid"
    assert contrast.concentration_or_dose == 200
    assert contrast.treatment_collection_age_or_stage == "Day 30"
    assert contrast.control_collection_age_or_stage == "Day 30"
    assert len(contrast.findings) == 0
    assert res.status == Status.PASS


# =========================================================================
# Acceptance Test 2: Day-30 treatment + Day-60 control => blocked; no eligible DE contrast
# =========================================================================
def test_acceptance_2_developmental_age_mismatch_blocks_contrast() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "rep1", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_cerebral", "biological_source_id": "source_01"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "rep2", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_cerebral", "biological_source_id": "source_01"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "c_rep1", "developmental_age_or_stage_at_collection": "Day 60", "organoid_context_id": "org_cerebral", "biological_source_id": "source_01"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "c_rep2", "developmental_age_or_stage_at_collection": "Day 60", "organoid_context_id": "org_cerebral", "biological_source_id": "source_01"},
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
    assert contrast.status == ContrastStatus.BLOCKED
    assert contrast.is_inferentially_eligible is False
    assert any(f.rule_id == "developmental_age_mismatch" and f.severity == Severity.ERROR for f in contrast.findings)
    assert len(res.eligible_contrasts) == 0
    assert len(res.blocked_contrasts) == 1
    assert res.status == Status.FAIL


# =========================================================================
# Acceptance Test 3: Same-age Day-30 + Day-30 and Day-60 + Day-60 => two separate contrasts
# =========================================================================
def test_acceptance_3_distinct_developmental_ages_produce_two_separate_contrasts() -> None:
    samples = [
        # Day 30 pair
        {"sample_id": "s_t30_1", "experiment_id": "exp_001", "condition_id": "treat_d30", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_01"},
        {"sample_id": "s_t30_2", "experiment_id": "exp_001", "condition_id": "treat_d30", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_01"},
        {"sample_id": "s_c30_1", "experiment_id": "exp_001", "condition_id": "ctrl_d30", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_01"},
        {"sample_id": "s_c30_2", "experiment_id": "exp_001", "condition_id": "ctrl_d30", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30", "organoid_context_id": "org_01"},
        # Day 60 pair
        {"sample_id": "s_t60_1", "experiment_id": "exp_001", "condition_id": "treat_d60", "biological_replicate_id": "r3", "developmental_age_or_stage_at_collection": "Day 60", "organoid_context_id": "org_01"},
        {"sample_id": "s_t60_2", "experiment_id": "exp_001", "condition_id": "treat_d60", "biological_replicate_id": "r4", "developmental_age_or_stage_at_collection": "Day 60", "organoid_context_id": "org_01"},
        {"sample_id": "s_c60_1", "experiment_id": "exp_001", "condition_id": "ctrl_d60", "biological_replicate_id": "cr3", "developmental_age_or_stage_at_collection": "Day 60", "organoid_context_id": "org_01"},
        {"sample_id": "s_c60_2", "experiment_id": "exp_001", "condition_id": "ctrl_d60", "biological_replicate_id": "cr4", "developmental_age_or_stage_at_collection": "Day 60", "organoid_context_id": "org_01"},
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

    assert c30.status == ContrastStatus.ELIGIBLE
    assert c30.treatment_collection_age_or_stage == "Day 30"
    assert c30.control_collection_age_or_stage == "Day 30"

    assert c60.status == ContrastStatus.ELIGIBLE
    assert c60.treatment_collection_age_or_stage == "Day 60"
    assert c60.control_collection_age_or_stage == "Day 60"

    assert len(res.eligible_contrasts) == 2


# =========================================================================
# Acceptance Test 4: Two treatment doses sharing one explicitly declared control => two contrasts
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
            "relationship_id": "rel_dose_series",
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
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "treat_orphan", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "treat_orphan", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "treat_orphan", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
    ]
    relationships = [
        {
            "relationship_id": "rel_orphan",
            "treatment_condition_ids": ["treat_orphan"],
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
# Acceptance Test 6: One biological replicate in treatment or control => not inferentially eligible
# =========================================================================
def test_acceptance_6_single_biological_replicate_is_not_inferentially_eligible() -> None:
    samples = [
        # Treatment has only 1 biological replicate
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        # Control has 2 biological replicates
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

    assert len(res.contrasts) == 1
    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.NEEDS_REVIEW
    assert contrast.is_inferentially_eligible is False
    assert contrast.treatment_replicate_count == 1
    assert contrast.control_replicate_count == 2
    assert any(f.rule_id == "insufficient_biological_replication" for f in contrast.findings)
    assert len(res.eligible_contrasts) == 0


# =========================================================================
# Acceptance Test 7: Multiple true biological replicates => retained as biological replicates
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
# Acceptance Test 8: Technical replicates must not inflate biological replicate count
# =========================================================================
def test_acceptance_8_technical_replicates_do_not_inflate_biological_count() -> None:
    # 4 samples in treatment, but all share the SAME biological_replicate_id 'bio_rep_1'
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "bio_rep_1", "technical_replicate_id": "tech_1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "bio_rep_1", "technical_replicate_id": "tech_2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t3", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "bio_rep_1", "technical_replicate_id": "tech_3", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t4", "experiment_id": "exp_001", "condition_id": "cond_treat", "biological_replicate_id": "bio_rep_1", "technical_replicate_id": "tech_4", "developmental_age_or_stage_at_collection": "Day 30"},
        # Control has 2 true biological replicates
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "c_bio_1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "c_bio_2", "developmental_age_or_stage_at_collection": "Day 30"},
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
    assert len(contrast.treatment_sample_ids) == 4
    # Unique biological replicates must be 1, not 4
    assert contrast.treatment_replicate_count == 1
    assert contrast.is_inferentially_eligible is False
    assert contrast.status == ContrastStatus.NEEDS_REVIEW
    assert any(f.rule_id == "insufficient_biological_replication" for f in contrast.findings)


# =========================================================================
# Acceptance Test 9: Multiple treatment conditions must never be pooled into one arm
# =========================================================================
def test_acceptance_9_multiple_treatment_conditions_never_pooled() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "treat_A", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "treat_A", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t3", "experiment_id": "exp_001", "condition_id": "treat_B", "biological_replicate_id": "r3", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t4", "experiment_id": "exp_001", "condition_id": "treat_B", "biological_replicate_id": "r4", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "treat_A", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "treat_B", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {
            "relationship_id": "rel_pooled_in_spec",
            "treatment_condition_ids": ["treat_A", "treat_B"],
            "matched_control_condition_ids": ["ctrl"],
        }
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    # Must produce 2 separate contrasts, NEVER a single [treat_A + treat_B] vs ctrl
    assert len(res.contrasts) == 2
    treatment_conds_produced = {c.treatment_condition_id for c in res.contrasts}
    assert treatment_conds_produced == {"treat_A", "treat_B"}
    for c in res.contrasts:
        assert len(c.treatment_sample_ids) == 2


# =========================================================================
# Acceptance Test 10: Cross-study pooling attempt => blocked
# =========================================================================
def test_acceptance_10_cross_study_comparison_blocked() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "treat_s1", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "treat_s1", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_002", "condition_id": "ctrl_s2", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_002", "condition_id": "ctrl_s2", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "treat_s1", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "ctrl_s2", "experiment_id": "exp_002", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_cross_study", "treatment_condition_ids": ["treat_s1"], "matched_control_condition_ids": ["ctrl_s2"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    # Add second experiment belonging to second study
    meta["experiments"].append({"experiment_id": "exp_002", "study_id": "study_002"})
    meta["studies"].append({"study_id": "study_002"})

    res = build_response_contrasts(norm, meta)

    assert len(res.contrasts) == 1
    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "cross_study_comparison_blocked" for f in contrast.findings)


# =========================================================================
# Acceptance Test 11: Cross-experiment comparison attempt => blocked
# =========================================================================
def test_acceptance_11_cross_experiment_comparison_blocked() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "treat_e1", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "treat_e1", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_002", "condition_id": "ctrl_e2", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_002", "condition_id": "ctrl_e2", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "treat_e1", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "ctrl_e2", "experiment_id": "exp_002", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_cross_exp", "treatment_condition_ids": ["treat_e1"], "matched_control_condition_ids": ["ctrl_e2"]},
    ]

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    # Add second experiment in SAME study
    meta["experiments"].append({"experiment_id": "exp_002", "study_id": "study_001"})

    res = build_response_contrasts(norm, meta)

    assert len(res.contrasts) == 1
    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "cross_experiment_comparison_blocked" for f in contrast.findings)


# =========================================================================
# Acceptance Test 12: Explicit sample comparison exception preserved and not silently ignored
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

    # 2. Without valid replacement control
    exceptions_no_rep = [
        {
            "exception_id": "exc_02",
            "relationship_id": "rel_exc",
            "affected_sample_ids": ["s_c1_bad"],
            "exception_reason": {"original_value": "batch anomaly without replacement", "value_status": "reported"},
        }
    ]
    norm2, meta2 = create_normalized_dataset_and_metadata(samples, conditions, relationships, exceptions_spec=exceptions_no_rep)
    res2 = build_response_contrasts(norm2, meta2)

    contrast2 = res2.contrasts[0]
    assert "exc_02" in contrast2.sample_comparison_exception_ids
    assert contrast2.status == ContrastStatus.NEEDS_REVIEW
    assert any(f.rule_id == "unresolved_sample_comparison_exception" for f in contrast2.findings)


# =========================================================================
# Acceptance Test 13: DNT labels/evidence do not affect contrast construction
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

    # In meta2, inject extraneous DNT annotations and candidate biomarker lists
    meta2 = copy.deepcopy(meta1)
    meta2["dnt_callset_evidence"] = {"positive_controls": ["cond_treat"], "dnt_label": "POSITIVE"}
    meta2["conditions"][0]["dnt_classification"] = "DNT_POSITIVE_BENCHMARK"
    meta2["conditions"][0]["known_dnt_biomarkers"] = ["TUBB3", "SOX2", "NES"]
    meta2["conditions"][1]["dnt_classification"] = "NEGATIVE_REFERENCE"

    norm2, _ = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res2 = build_response_contrasts(norm2, meta2)

    assert len(res1.contrasts) == len(res2.contrasts)
    c1 = res1.contrasts[0]
    c2 = res2.contrasts[0]
    assert c1.contrast_id == c2.contrast_id
    assert c1.status == c2.status
    assert c1.treatment_sample_ids == c2.treatment_sample_ids
    assert c1.control_sample_ids == c2.control_sample_ids
    assert c1.treatment_replicate_count == c2.treatment_replicate_count
    assert c1.control_replicate_count == c2.control_replicate_count
    assert c1.is_inferentially_eligible == c2.is_inferentially_eligible


# =========================================================================
# Acceptance Test 14: Input ordering changes do not change contrast IDs or output ordering
# =========================================================================
def test_acceptance_14_input_ordering_invariance() -> None:
    samples = [
        {"sample_id": "s_t1", "experiment_id": "exp_001", "condition_id": "cond_treat_A", "biological_replicate_id": "r1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t2", "experiment_id": "exp_001", "condition_id": "cond_treat_A", "biological_replicate_id": "r2", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t3", "experiment_id": "exp_001", "condition_id": "cond_treat_B", "biological_replicate_id": "r3", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_t4", "experiment_id": "exp_001", "condition_id": "cond_treat_B", "biological_replicate_id": "r4", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c1", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr1", "developmental_age_or_stage_at_collection": "Day 30"},
        {"sample_id": "s_c2", "experiment_id": "exp_001", "condition_id": "cond_ctrl", "biological_replicate_id": "cr2", "developmental_age_or_stage_at_collection": "Day 30"},
    ]
    conditions = [
        {"condition_id": "cond_treat_A", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_treat_B", "experiment_id": "exp_001", "treatment_control_status": "treatment"},
        {"condition_id": "cond_ctrl", "experiment_id": "exp_001", "treatment_control_status": "control"},
    ]
    relationships = [
        {"relationship_id": "rel_01", "treatment_condition_ids": ["cond_treat_A", "cond_treat_B"], "matched_control_condition_ids": ["cond_ctrl"]},
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
    # Verify that raw counts are preserved in the cohort result linked to this contrast
    cohort_result = norm.cohort_results[0]
    assert contrast.cohort_id == cohort_result.cohort.cohort_id
    assert cohort_result.raw_count_matrix is not None
    # Confirm raw count matrix contains discrete non-negative integer counts
    for col in cohort_result.raw_count_matrix.columns:
        for val in col:
            assert isinstance(val, int)
            assert val >= 0
    # Confirm diagnostic matrix is floating-point inspection data and distinct from raw count matrix
    assert cohort_result.diagnostic_normalized_matrix is not None
    for col in cohort_result.diagnostic_normalized_matrix.columns:
        for val in col:
            assert isinstance(val, float)


# =========================================================================
# Additional Scientific Invariant Tests
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

    norm, meta = create_normalized_dataset_and_metadata(samples, conditions, relationships)
    res = build_response_contrasts(norm, meta)

    contrast = res.contrasts[0]
    assert contrast.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "organoid_context_mismatch" for f in contrast.findings)


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

    contrast = res.contrasts[0]
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
    assert normalize_age_string("differentiated_stage_A") == "differentiated_stage_a"
    assert normalize_age_string(None) is None
    assert normalize_age_string("unknown") is None
    assert normalize_age_string("missing_not_reported") is None
