"""Unit tests for response builder domain models."""

from src.response_builder.models import (
    ContrastExposure,
    ContrastStatus,
    Finding,
    MolecularResponseContrast,
    ResponseContrastDataset,
    Severity,
    Status,
    ordered_findings,
    status_for,
)


def test_contrast_status_enum() -> None:
    assert ContrastStatus.ELIGIBLE == "ELIGIBLE"
    assert ContrastStatus.NEEDS_REVIEW == "NEEDS_REVIEW"
    assert ContrastStatus.BLOCKED == "BLOCKED"


def test_ordered_findings() -> None:
    f_info = Finding(
        severity=Severity.INFO,
        rule_id="rule_info",
        entity_type="sample",
        entity_id="s1",
        path="samples.s1",
        message="Info finding",
    )
    f_err = Finding(
        severity=Severity.ERROR,
        rule_id="rule_err",
        entity_type="sample",
        entity_id="s2",
        path="samples.s2",
        message="Error finding",
    )
    f_rev = Finding(
        severity=Severity.REVIEW,
        rule_id="rule_rev",
        entity_type="sample",
        entity_id="s3",
        path="samples.s3",
        message="Review finding",
    )

    ordered = ordered_findings([f_info, f_err, f_rev])
    assert ordered[0].severity == Severity.ERROR
    assert ordered[1].severity == Severity.REVIEW
    assert ordered[2].severity == Severity.INFO


def test_status_for() -> None:
    f_info = Finding(Severity.INFO, "r1", "s", "1", "p", "m")
    f_warn = Finding(Severity.WARNING, "r2", "s", "2", "p", "m")
    f_rev = Finding(Severity.REVIEW, "r3", "s", "3", "p", "m")
    f_err = Finding(Severity.ERROR, "r4", "s", "4", "p", "m")

    assert status_for([]) == Status.PASS
    assert status_for([f_info]) == Status.PASS
    assert status_for([f_info, f_warn]) == Status.PASS_WITH_WARNINGS
    assert status_for([f_info, f_warn, f_rev]) == Status.NEEDS_REVIEW
    assert status_for([f_info, f_warn, f_rev, f_err]) == Status.FAIL


def test_contrast_exposure_model() -> None:
    exp = ContrastExposure(
        exposure_id="exp_01",
        agent_name="Valproic Acid",
        agent_identifier="CID_107",
        vehicle="DMSO",
        concentration_or_dose=200.0,
        concentration_or_dose_unit="uM",
        exposure_start_time_or_stage="Day 25",
        developmental_age_or_stage_at_exposure="Day 25",
        exposure_duration=5.0,
        exposure_duration_unit="days",
        washout_or_recovery_duration=0.0,
        washout_or_recovery_duration_unit="days",
    )
    d = exp.to_dict()
    assert d["exposure_id"] == "exp_01"
    assert d["agent_name"] == "Valproic Acid"
    assert d["concentration_or_dose"] == 200.0


def test_molecular_response_contrast_replication_and_eligibility() -> None:
    exp = ContrastExposure(
        exposure_id="exp_01",
        agent_name="Valproic Acid",
        agent_identifier="CID_107",
        vehicle="DMSO",
        concentration_or_dose=200,
        concentration_or_dose_unit="uM",
    )
    contrast = MolecularResponseContrast(
        contrast_id="contrast_01",
        study_id="study_01",
        experiment_id="exp_01",
        cohort_id="cohort_01",
        relationship_id="rel_01",
        treatment_condition_id="cond_treat",
        matched_control_condition_ids=("cond_ctrl",),
        treatment_sample_ids=("s_t1", "s_t2", "s_t3"),
        control_sample_ids=("s_c1", "s_c2"),
        treatment_biological_replicate_ids=("rep1", "rep2", "rep3"),
        control_biological_replicate_ids=("c_rep1", "c_rep2"),
        treatment_exposures=(exp,),
        agent_name="Valproic Acid",
        agent_identifier="CID_107",
        vehicle="DMSO",
        concentration_or_dose=200,
        concentration_or_dose_unit="uM",
        exposure_start_time_or_stage="Day 25",
        developmental_age_or_stage_at_exposure="Day 25",
        exposure_duration=5,
        exposure_duration_unit="days",
        washout_or_recovery_duration=0,
        washout_or_recovery_duration_unit="days",
        treatment_collection_age_or_stage="Day 30",
        control_collection_age_or_stage="Day 30",
        treatment_collection_age_normalized="day_30",
        control_collection_age_normalized="day_30",
        biological_source_ids=("bio_01",),
        organoid_context_ids=("org_01",),
        sample_comparison_exception_ids=(),
        status=ContrastStatus.ELIGIBLE,
        findings=(),
    )

    assert contrast.treatment_replicate_count == 3
    assert contrast.control_replicate_count == 2
    assert contrast.is_inferentially_eligible is True
    d = contrast.to_dict()
    assert d["treatment_replicate_count"] == 3
    assert d["is_inferentially_eligible"] is True
    assert len(d["treatment_exposures"]) == 1
    assert d["treatment_collection_age_normalized"] == "day_30"


def test_molecular_response_contrast_under_replicated_is_not_eligible() -> None:
    contrast = MolecularResponseContrast(
        contrast_id="contrast_under_replicated",
        study_id="study_01",
        experiment_id="exp_01",
        cohort_id="cohort_01",
        relationship_id="rel_01",
        treatment_condition_id="cond_treat",
        matched_control_condition_ids=("cond_ctrl",),
        treatment_sample_ids=("s_t1",),
        control_sample_ids=("s_c1", "s_c2"),
        treatment_biological_replicate_ids=("rep1",),
        control_biological_replicate_ids=("c_rep1", "c_rep2"),
        treatment_exposures=(),
        agent_name="Valproic Acid",
        agent_identifier="CID_107",
        vehicle="DMSO",
        concentration_or_dose=200,
        concentration_or_dose_unit="uM",
        exposure_start_time_or_stage="Day 25",
        developmental_age_or_stage_at_exposure="Day 25",
        exposure_duration=5,
        exposure_duration_unit="days",
        washout_or_recovery_duration=0,
        washout_or_recovery_duration_unit="days",
        treatment_collection_age_or_stage="Day 30",
        control_collection_age_or_stage="Day 30",
        biological_source_ids=("bio_01",),
        organoid_context_ids=("org_01",),
        sample_comparison_exception_ids=(),
        status=ContrastStatus.NEEDS_REVIEW,
        findings=(),
    )

    assert contrast.treatment_replicate_count == 1
    assert contrast.control_replicate_count == 2
    assert contrast.is_inferentially_eligible is False
