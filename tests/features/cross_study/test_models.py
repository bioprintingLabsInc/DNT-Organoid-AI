"""Unit tests for Stage 8A.1 data models and enums."""

import pytest

from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
    CrossStudyStatus,
    EstimateStatus,
    ResponseDirection,
)


def test_response_direction_values() -> None:
    assert ResponseDirection.UP == "UP"
    assert ResponseDirection.DOWN == "DOWN"
    assert ResponseDirection.NO_DIRECTION == "NO_DIRECTION"
    assert ResponseDirection.UNDETERMINED == "UNDETERMINED"


def test_estimate_status_values() -> None:
    assert EstimateStatus.OK == "OK"
    assert EstimateStatus.ZERO_COUNTS == "ZERO_COUNTS"
    assert EstimateStatus.COOKS_OUTLIER == "COOKS_OUTLIER"
    assert EstimateStatus.FIT_ERROR == "FIT_ERROR"
    assert EstimateStatus.INVALID_STATISTIC == "INVALID_STATISTIC"


def test_cross_study_status_values() -> None:
    assert CrossStudyStatus.SINGLE_STUDY_ONLY == "SINGLE_STUDY_ONLY"
    assert CrossStudyStatus.MULTI_STUDY_CONTEXT_UNDER_REVIEW == "MULTI_STUDY_CONTEXT_UNDER_REVIEW"
    assert CrossStudyStatus.NONCOMPARABLE_CONTEXT == "NONCOMPARABLE_CONTEXT"
    assert CrossStudyStatus.COMPARABLE_INDEPENDENT_STUDIES == "COMPARABLE_INDEPENDENT_STUDIES"
    assert CrossStudyStatus.INSUFFICIENT_METADATA == "INSUFFICIENT_METADATA"


def test_record_immutability() -> None:
    rec = ConditionGeneResponseRecord(
        canonical_ensembl_gene_id="ENSG00000000001",
        gene_symbol="TEST_GENE",
        study_id="STUDY_01",
        study_family_id="STUDY_01",
        independent_study_id="STUDY_01",
        experiment_id="EXP_01",
        contrast_id="CONTRAST_01",
        treatment_condition_id="TRT",
        matched_control_condition_ids=("CTL",),
        compound_id="CID:123",
        compound_name="TestCompound",
        dose_value=10.0,
        dose_unit="uM",
        dose_original="10 uM",
        exposure_duration=24.0,
        exposure_duration_unit="h",
        developmental_stage="Day 10",
        model_type="Cortical Organoid",
        biological_source_ids=("Donor1",),
        organoid_context_ids=("ctx1",),
        treatment_replicate_count=3,
        control_replicate_count=3,
        treatment_sample_ids=("S1", "S2", "S3"),
        control_sample_ids=("S4", "S5", "S6"),
        treatment_biological_replicate_ids=("R1", "R2", "R3"),
        control_biological_replicate_ids=("R4", "R5", "R6"),
        base_mean=100.0,
        log2_fold_change=1.5,
        log2_fold_change_se=0.2,
        wald_statistic=7.5,
        p_value=1e-10,
        adjusted_p_value=1e-8,
        response_direction=ResponseDirection.UP,
        estimate_status=EstimateStatus.OK,
        missing_reason=None,
        cross_study_status=CrossStudyStatus.SINGLE_STUDY_ONLY,
    )

    with pytest.raises(Exception):  # FrozenInstanceError / AttributeError
        rec.log2_fold_change = 2.0  # type: ignore

    d = rec.to_dict()
    assert d["canonical_ensembl_gene_id"] == "ENSG00000000001"
    assert d["log2_fold_change"] == 1.5
    assert d["response_direction"] == "UP"
    assert d["estimate_status"] == "OK"
    assert d["cross_study_status"] == "SINGLE_STUDY_ONLY"


def test_registry_properties() -> None:
    rec1 = ConditionGeneResponseRecord(
        canonical_ensembl_gene_id="ENSG01",
        gene_symbol="G1",
        study_id="S1",
        study_family_id="S1",
        independent_study_id="S1",
        experiment_id="E1",
        contrast_id="C1",
        treatment_condition_id="T1",
        matched_control_condition_ids=("CTL",),
        compound_id="CMP1",
        compound_name="Compound1",
        dose_value=10.0,
        dose_unit="uM",
        dose_original="10 uM",
        exposure_duration=24.0,
        exposure_duration_unit="h",
        developmental_stage="D1",
        model_type="M1",
        biological_source_ids=(),
        organoid_context_ids=(),
        treatment_replicate_count=2,
        control_replicate_count=2,
        treatment_sample_ids=("S1", "S2"),
        control_sample_ids=("S3", "S4"),
        treatment_biological_replicate_ids=(),
        control_biological_replicate_ids=(),
        base_mean=50.0,
        log2_fold_change=1.2,
        log2_fold_change_se=0.2,
        wald_statistic=6.0,
        p_value=1e-5,
        adjusted_p_value=1e-4,
        response_direction=ResponseDirection.UP,
        estimate_status=EstimateStatus.OK,
        missing_reason=None,
    )

    rec2 = ConditionGeneResponseRecord(
        canonical_ensembl_gene_id="ENSG02",
        gene_symbol="G2",
        study_id="S1",
        study_family_id="S1",
        independent_study_id="S1",
        experiment_id="E1",
        contrast_id="C1",
        treatment_condition_id="T1",
        matched_control_condition_ids=("CTL",),
        compound_id="CMP1",
        compound_name="Compound1",
        dose_value=10.0,
        dose_unit="uM",
        dose_original="10 uM",
        exposure_duration=24.0,
        exposure_duration_unit="h",
        developmental_stage="D1",
        model_type="M1",
        biological_source_ids=(),
        organoid_context_ids=(),
        treatment_replicate_count=2,
        control_replicate_count=2,
        treatment_sample_ids=("S1", "S2"),
        control_sample_ids=("S3", "S4"),
        treatment_biological_replicate_ids=(),
        control_biological_replicate_ids=(),
        base_mean=0.0,
        log2_fold_change=None,
        log2_fold_change_se=None,
        wald_statistic=None,
        p_value=None,
        adjusted_p_value=None,
        response_direction=ResponseDirection.UNDETERMINED,
        estimate_status=EstimateStatus.ZERO_COUNTS,
        missing_reason="baseMean == 0",
    )

    registry = ConditionGeneResponseRegistry(
        records=(rec1, rec2),
        registry_id="reg_test",
    )

    assert registry.total_records == 2
    assert registry.unique_canonical_genes == {"ENSG01", "ENSG02"}
    assert registry.unique_contrasts == {"C1"}
    assert registry.unique_compounds == {"CMP1"}
    assert registry.unique_independent_studies == {"S1"}
    assert len(registry.eligible_records) == 1
    assert registry.eligible_records[0].canonical_ensembl_gene_id == "ENSG01"
    assert len(registry.get_records_for_gene("ENSG02")) == 1
    assert registry.get_records_for_gene("ENSG02")[0].estimate_status == EstimateStatus.ZERO_COUNTS
