"""Unit tests for Stage 8A.1 Quality Control Auditor."""

import pytest

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
    EstimateStatus,
    ResponseDirection,
)
from src.features.cross_study.qc import (
    audit_registry_qc,
)
from src.features.cross_study.registry import (
    build_condition_gene_response_registry,
)


def test_audit_passes_on_clean_registry(
    synthetic_contrast_vpa_200: DifferentialExpressionContrastResult,
) -> None:
    registry = build_condition_gene_response_registry(
        contrasts=[synthetic_contrast_vpa_200],
        registry_id="reg_clean",
    )
    qc = audit_registry_qc(registry, [synthetic_contrast_vpa_200])
    assert qc.is_valid is True
    assert qc.numerical_mismatches_count == 0
    assert qc.gene_identity_mismatches_count == 0
    assert qc.direction_mismatches_count == 0
    assert qc.context_mismatches_count == 0
    assert qc.leakage_violations_count == 0
    assert len(qc.findings) == 0


def test_audit_detects_numerical_mismatch(
    synthetic_contrast_vpa_200: DifferentialExpressionContrastResult,
) -> None:
    registry = build_condition_gene_response_registry(
        contrasts=[synthetic_contrast_vpa_200],
        registry_id="reg_dirty",
    )

    # Tamper with one record's log2_fold_change
    tampered_records = list(registry.records)
    r_bad = tampered_records[0]
    tampered_records[0] = ConditionGeneResponseRecord(
        canonical_ensembl_gene_id=r_bad.canonical_ensembl_gene_id,
        gene_symbol=r_bad.gene_symbol,
        study_id=r_bad.study_id,
        study_family_id=r_bad.study_family_id,
        independent_study_id=r_bad.independent_study_id,
        experiment_id=r_bad.experiment_id,
        contrast_id=r_bad.contrast_id,
        treatment_condition_id=r_bad.treatment_condition_id,
        matched_control_condition_ids=r_bad.matched_control_condition_ids,
        compound_id=r_bad.compound_id,
        compound_name=r_bad.compound_name,
        dose_value=r_bad.dose_value,
        dose_unit=r_bad.dose_unit,
        dose_original=r_bad.dose_original,
        exposure_duration=r_bad.exposure_duration,
        exposure_duration_unit=r_bad.exposure_duration_unit,
        developmental_stage=r_bad.developmental_stage,
        model_type=r_bad.model_type,
        biological_source_ids=r_bad.biological_source_ids,
        organoid_context_ids=r_bad.organoid_context_ids,
        treatment_replicate_count=r_bad.treatment_replicate_count,
        control_replicate_count=r_bad.control_replicate_count,
        treatment_sample_ids=r_bad.treatment_sample_ids,
        control_sample_ids=r_bad.control_sample_ids,
        treatment_biological_replicate_ids=r_bad.treatment_biological_replicate_ids,
        control_biological_replicate_ids=r_bad.control_biological_replicate_ids,
        base_mean=r_bad.base_mean,
        log2_fold_change=99.999,  # TAMPERED
        log2_fold_change_se=r_bad.log2_fold_change_se,
        wald_statistic=r_bad.wald_statistic,
        p_value=r_bad.p_value,
        adjusted_p_value=r_bad.adjusted_p_value,
        response_direction=r_bad.response_direction,
        estimate_status=r_bad.estimate_status,
        missing_reason=r_bad.missing_reason,
    )

    tampered_registry = ConditionGeneResponseRegistry(
        records=tuple(tampered_records),
        registry_id="reg_tampered",
    )

    qc = audit_registry_qc(tampered_registry, [synthetic_contrast_vpa_200])
    assert qc.is_valid is False
    assert qc.numerical_mismatches_count > 0
    assert any("log2FC mismatch" in f for f in qc.findings)
